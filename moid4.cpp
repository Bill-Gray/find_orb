/* moid4.cpp: computes MOID (Minimum Orbital Intersection Distance)
between two orbits.

Copyright (C) 2010, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */

#ifdef TEST_VERSION
#include <stdio.h>
#endif

#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "comets.h"
#include "afuncs.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int debug_printf( const char *format, ...);                 /* mpc_obs.c */
double find_moid( const ELEMENTS *elem1, const ELEMENTS *elem2,  /* moid4.c */
                                     double *barbee_style_delta_v);
int setup_planet_elem( ELEMENTS *elem, const int planet_idx,
                                          const double t_cen);   /* moid4.c */

static void fill_matrix( double mat[3][3], const ELEMENTS *elem)
{
   memcpy( mat[0], elem->perih_vec, 3 * sizeof( double));
   memcpy( mat[1], elem->sideways, 3 * sizeof( double));
         /* mat[2] is the cross-product of mat[0] & mat[1]: */
   vector_cross_product( mat[2], mat[0], mat[1]);
// mat[2][0] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];
// mat[2][1] = mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2];
// mat[2][2] = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

static double compute_posn_and_derivative( const ELEMENTS *elem,
            const double true_anom, const double matrix[3][3],
            double *posn, double *vel)
{
   const double cos_true_anom = cos( true_anom);
   const double sin_true_anom = sin( true_anom);
   const double denom = 1. + elem->ecc * cos_true_anom;
   const double true_r = elem->q * (1. + elem->ecc) / denom;
   const double x = true_r * cos_true_anom;
   const double y = true_r * sin_true_anom;
   const double dx_dtheta = -y / denom;
   const double dy_dtheta = (x + elem->ecc / denom) / denom;
   int i;

   for( i = 0; i < 3; i++)
      posn[i] = x * matrix[i][0] + y * matrix[i][1];
   if( vel)
      for( i = 0; i < 3; i++)
         vel[i] = dx_dtheta * matrix[i][0] + dy_dtheta * matrix[i][1];
   return( true_r);
}

#define dot_prod( a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])

static void compute_improvement( const double *delta, const double *v1,
               const double *v2, double *d1, double *d2)
{
   const double b = 2. * dot_prod( delta, v1);
   const double c = 2. * dot_prod( delta, v2);
   const double d = 2. * dot_prod( v1, v2);
   const double e = dot_prod( v1, v1);
   const double f = dot_prod( v2, v2);

   *d1 = (d * c - 2. * f * b) / (4. * e * f - d * d);
// *d2 = (d * b - 2. * e * c) / (4. * e * f - d * d);
   *d2 = -(b + 2. * e * *d1) / d;
}

#ifdef TEST_VERSION
static double true_anomaly_to_eccentric( const double true_anom,
                                         const double ecc)
{
   const double r = (1. - ecc * ecc) / (1. + ecc * cos( true_anom));
   const double x = r * cos( true_anom) + ecc;
   const double y = r * sin( true_anom) / sqrt( 1. - ecc * ecc);
   const double ecc_anom = PI + atan2( -y, -x);

   return( ecc_anom);
}
#endif

/* In computing the MOID,  our "velocity" is really the derivative
of the object's position with respect to true anomaly.  When it's time
to compute the relative velocity of the two objects at the MOID point,
we want a for-real,  true velocity:  the vector describing,  in km/s,
the velocity of each object relative to the center.  We can then
subtract the vel2 vector from the vel1 vector,  and we'll know the
velocity adjustment required to hop from one orbit to the other.  */

static void set_true_velocity( double *vel_vect, const double r,
                        const double semimajor_axis)
{
   const double escape_velocity_at_one_au = 42.1219;     /* in km/s */
   const double space_vel =
               escape_velocity_at_one_au * sqrt( 1. / r - .5 / semimajor_axis);
   const double vel_vect_length = vector3_length( vel_vect);
   size_t idx;

   for( idx = 0; idx < 3; idx++)
      vel_vect[idx] *= space_vel / vel_vect_length;
}

#define N_STEPS 1080

/* A simple,  but effective,  MOID-finder.  It starts by computing 3x3
orthonormal matrices mat1 and mat2 that correspond to the base vectors of
the two orbits (i.e.,  each has a vector pointing toward perihelion;
one pointing 90 degrees "ahead" in the orbit;  and one perpendicular to
the plane of the orbit.)  Multiplying one by the inverse of the other
gives the transformation matrix from the first orbit to the second.
Then,  we can work as if one orbit is in the xy plane with perihelion
toward the positive x-axis.

   Optionally,  if barbee_style_delta_v != NULL,  the relative speed in
km/s at the MOID point will be determined.  The idea is that if the objects
were to pass close to one another at that point,  you could push off from
one of them at that speed and match orbits with the other.  This is a
decent approximation to figuring out how "easy" it is to get from one
object (usually the earth) to the other (usually an asteroid).  The idea
came from an exchange of e-mails with Brent W. Barbee,  of NASA's
Goddard Space Flight Center (GSFC).  It should be noted that there are
other ways of defining the encounter velocity,  including one due to
Alan Harris -- see 'elem_out.cpp' -- and one described by E. M. Shoemaker
and E. F. Helin in 1978, "Earth-Approaching Asteroids as Targets for
Exploration",  NASA CP-2053, pp. 245-256.  */

/* A small point:  'true_anomaly1' _will_ be initialized when the
first loop runs.  g++ doesn't realize that,  though,  and emits a
warning.  With gcc 4.6 or better,  this can be suppressed with a
#pragma GCC diagnostic.  With earlier versions,  the only way around
the problem is to do a superfluous initialization of 'true_anomaly1'. */

double find_moid( const ELEMENTS *elem1, const ELEMENTS *elem2,
                                     double *barbee_style_delta_v)
{
   double mat1[3][3], mat2[3][3], xform_matrix[3][3];
   const double identity_matrix[3][3] = {
            { 1., 0., 0.},
            { 0., 1., 0.},
            { 0., 0., 1.} };
   double least_dist_squared = 10000.;
   int i, j;

   fill_matrix( mat1, elem1);
   fill_matrix( mat2, elem2);
   for( i = 0; i < 3; i++)
      for( j = 0; j < 3; j++)
         xform_matrix[j][i] = dot_prod( mat1[j], mat2[i]);

   for( i = 0; i < N_STEPS; i++)
      {
      double vect1[3], vect2[3], dist_squared = 0.;
      double deriv1[3], deriv2[3], r1, r2;
      double true_anomaly2 = 2. * PI * (double)i / (double)N_STEPS;
      double delta_true1, delta_true2;
      int loop_count = 0, solution_found = 0;
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
   #pragma GCC diagnostic push              /* see comments above */
   #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
      double true_anomaly1;
   #pragma GCC diagnostic pop
#else
      double true_anomaly1 = 0.;
#endif

      do
         {
         r2 = compute_posn_and_derivative( elem2, true_anomaly2, xform_matrix,
                        vect2, deriv2);
         if( !loop_count)
            true_anomaly1 = atan2( vect2[1], vect2[0]);
         r1 = compute_posn_and_derivative( elem1, true_anomaly1, identity_matrix,
                        vect1, deriv1);
         for( j = 0; j < 3; j++)
            vect1[j] -= vect2[j];
         compute_improvement( vect1, deriv1, deriv2, &delta_true1, &delta_true2);
         true_anomaly1 += delta_true1;
         true_anomaly2 -= delta_true2;
         if( fabs( delta_true1) < 5. * PI / N_STEPS)
            if( fabs( delta_true2) < 5. * PI / N_STEPS)
               {
               for( j = 0; j < 3; j++)
                  vect1[j] += delta_true1 * deriv1[j] + delta_true2 * deriv2[j];
               solution_found = 1;
               }
         loop_count++;
//       debug_printf( "  i = %3d; loop %d; %f\n",
//                i, loop_count, sqrt( dot_prod( vect1, vect1)));
         }
         while( solution_found && loop_count < 5);
      dist_squared = dot_prod( vect1, vect1);
      if( dist_squared < least_dist_squared)
         {
         least_dist_squared = dist_squared;
         if( barbee_style_delta_v)
            {
            double delta_v[3];

            set_true_velocity( deriv1, r1, elem1->major_axis);
            set_true_velocity( deriv2, r2, elem2->major_axis);
            for( j = 0; j < 3; j++)
               delta_v[j] = deriv1[j] - deriv2[j];
            *barbee_style_delta_v = vector3_length( delta_v);
            }
         }
#ifdef TEST_VERSION
      printf( "%3d %c%8.6f%8.2f%8.2f%8.2f%8.2f%15f%15f\n", i,
                     (solution_found ? '*' : ' '),
                     sqrt( dot_prod( vect1, vect1)),
                     true_anomaly1 * 180. / PI,
                     true_anomaly2 * 180. / PI,
                     true_anomaly_to_eccentric( true_anomaly1, elem1->ecc) * 180. / PI,
                     true_anomaly_to_eccentric( true_anomaly2, elem2->ecc) * 180. / PI,
                     dot_prod( vect1, deriv1),
                     dot_prod( vect1, deriv2));
//    printf( "%3d%15f%15f%15f%15f%15f\n", i, x, y,
//                      vect[0], vect[1], vect[2]);
#endif
      }
   return( sqrt( least_dist_squared));
}

#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

#define N_PLANET_ELEMS 15
#define N_PLANET_RATES 9

int setup_planet_elem( ELEMENTS *elem, const int planet_idx,
                                             const double t_cen)
{
/* Taken straight from http://ssd.jpl.nasa.gov/elem_planets.html */
/* Gives a, ecc, incl, Omega=asc node, omega=arg per, L=longit   */
/* Slightly different values are given at:                       */
/* http://ssd.jpl.nasa.gov/txt/p_elem_t1.txt    and:             */
/* http://ssd.jpl.nasa.gov/txt/p_elem_t2.txt                     */
/* Note that for MOID-finding,  the longitude is irrelevant.  I  */
/* left it in on the assumption that we might want it someday.   */
/* Asteroid elements are from BC-405,  epoch 2451535.0.  The     */
/* "longitudes" are actually mean anomalies,  and the "LonPer"s  */
/* are actually arguments of perihelion;  the code corrects for  */
/* that last,  and the mean anomaly/longitude isn't used (yet).  */

static const double planet_elem[N_PLANET_ELEMS * 6] = {
         /*     a          eccent    inclin    AscNode   LonPer      Longit */
/* Merc */  0.38709893, .20563069,  7.00487,  48.33167,  77.45645, 252.25084,
/* Venu */  0.72333199, .00677323,  3.39471,  76.68069, 131.53298, 181.97973,
/* Eart */  1.00000011, .01671022,  0.00005, -11.26064, 102.94719, 100.46435,
/* Mars */  1.52366231, .09341233,  1.85061,  49.57854, 336.04084, 355.45332,
/* Jupi */  5.20336301, .04839266,  1.30530, 100.55615,  14.75385,  34.40438,
/* Satu */  9.53707032, .05415060,  2.48446, 113.71504,  92.43194,  49.94432,
/* Uran */ 19.19126393, .04716771,  0.76986,  74.22988, 170.96424, 313.23218,
/* Nept */ 30.06896348, .00858587,  1.76917, 131.72169,  44.97135, 304.88003,
/* Plut */ 39.48168677, .24880766, 17.14175, 110.30347, 224.06676, 238.92881,
/*  (1) */  2.7664603,  0.0783638, 10.583360,  80.494464,  73.921341,   4.036019,
/*  (2) */  2.7723257,  0.2296435, 34.846130, 173.197757, 310.264059, 350.826074,
/*  (4) */  2.3615363,  0.0900245,  7.133918, 103.951631, 149.589094, 338.305822,
/* (29) */  2.5543838,  0.0722511,  6.102741, 356.567840,  62.015715,  20.150301,
/* (16) */  2.9204983,  0.1382234,  3.093382, 150.465894, 229.122381, 333.613957,
/* (15) */  2.6437135,  0.1862108, 11.747399, 293.516504,  96.956836, 104.873024 };

static const double planet_elem_rate[N_PLANET_RATES * 6] = {
/* Merc */  0.00000066,  0.00002527, -23.51,   -446.30,   573.57, 538101628.29,
/* Venu */  0.00000092, -0.00004938,  -2.86,   -996.89,  -108.80, 210664136.06,
/* Eart */ -0.00000005, -0.00003804, -46.94, -18228.25,  1198.28, 129597740.63,
/* Mars */ -0.00007221,  0.00011902, -25.47,  -1020.19,  1560.78,  68905103.78,
/* Jupi */  0.00060737, -0.00012880,  -4.15,   1217.17,   839.93,  10925078.35,
/* Satu */ -0.00301530, -0.00036762,   6.11,  -1591.05, -1948.89,   4401052.95,
/* Uran */  0.00152025, -0.00019150,  -2.09,  -1681.40,  1312.56,   1542547.79,
/* Nept */ -0.00125196,  0.0000251,   -3.64,   -151.25,  -844.43,    786449.21,
/* Plut */ -0.00076912,  0.00006465,  11.07,    -37.33,  -132.25,    522747.90};

   const double *pdata = planet_elem + (planet_idx - 1) * 6;
   const double *rate_data = planet_elem_rate + (planet_idx - 1) * 6;
   double elem_array[6];
   int i;

   if( planet_idx >= N_PLANET_ELEMS || planet_idx < 0)
      return( -1);
   for( i = 0; i < 6; i++)
      if( i < 2)
         elem_array[i] = pdata[i];
      else
         elem_array[i] = pdata[i] * PI / 180.;
   if( planet_idx <= N_PLANET_RATES)
      {
      for( i = 0; i < 6; i++)
         if( i < 2)
            elem_array[i] += rate_data[i] * t_cen;
         else
            elem_array[i] += (rate_data[i] * t_cen / 3600.) * PI / 180.;
      }
   memset( elem, 0, sizeof( ELEMENTS));
   elem->ecc = elem_array[1];
   elem->q = (1. - elem->ecc) * elem_array[0];
   elem->incl = elem_array[2];
   elem->asc_node = elem_array[3];
               /* For planets,  the longitude of perihelion is given. */
               /* For asteroids,  it's the argument of perihelion.    */
   if( planet_idx < 9)
      elem->arg_per = elem_array[4] - elem_array[3];
   else
      elem->arg_per = elem_array[4];
      /* l = (100.46435 + (129597740.63 / 3600.) * t_cen) * PI / 180.; */
   derive_quantities( elem, SOLAR_GM);
   return( 0);
}

#ifdef TEST_VERSION
#include <stdlib.h>

static void show_elements( const ELEMENTS *elem)
{
   printf( "q=%8.5f e=%8.6f i=%8.4f asc_node=%8.4f arg_per=%8.4f\n",
               elem->q, elem->ecc,
               elem->incl * 180. / PI,
               elem->asc_node * 180. / PI,
               elem->arg_per * 180. / PI);
}

int main( const int argc, const char **argv)
{
   ELEMENTS elem, earth_elem;
   double t_cen = 0.06, barbee_style_vel;

   memset( &elem, 0, sizeof( ELEMENTS));
   if( argc == 6)
      setup_planet_elem( &earth_elem, 3, t_cen);
   else
      {
      memset( &earth_elem, 0, sizeof( ELEMENTS));
      sscanf( argv[6], "%lf,%lf,%lf,%lf,%lf",
                  &earth_elem.q,
                  &earth_elem.ecc,
                  &earth_elem.incl,
                  &earth_elem.asc_node,
                  &earth_elem.arg_per);
      earth_elem.incl *= PI / 180.;
      earth_elem.asc_node *= PI / 180.;
      earth_elem.arg_per *= PI / 180.;
      }
   elem.q = atof( argv[1]);
   elem.ecc = atof( argv[2]);
   if( elem.q < 0.)             /* actually the semimajor axis was given; */
      elem.q *= elem.ecc - 1.;  /* cvt it to a perihelion distance */
   elem.incl = atof( argv[3]) * PI / 180.;
   elem.asc_node = atof( argv[4]) * PI / 180.;
   elem.arg_per = atof( argv[5]) * PI / 180.;
   derive_quantities( &elem, SOLAR_GM);
   derive_quantities( &earth_elem, SOLAR_GM);
   printf( "MOID = %f\n", find_moid( &earth_elem, &elem, &barbee_style_vel));
   printf( "Barbee-style encounter vel = %f\n", barbee_style_vel);
   show_elements( &elem);
   if( argc != 6)
      show_elements( &earth_elem);
   return( 0);
}
#endif
