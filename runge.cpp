/* runge.cpp: numerical integration code, mostly force model stuff

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#ifdef __GNUC__
#include <unistd.h>
#endif
#include "watdefs.h"
#include "lunar.h"
#include "afuncs.h"
#include "comets.h"
#include "afuncs.h"
#include "mpc_obs.h"
#include "constant.h"

#define ldouble long double

#if defined( __WATCOMC__) || defined( __FreeBSD__)
#define ceill ceil
#define expl exp
#define fabsl fabs
#define isnanl isnan
#define powl pow
#define sqrtl sqrt
#endif

/* perturbers
   excluded_perturbers = -1;
   general_relativity_factor
   n_orbit_params
   planet_mass[],  sort of
   best_fit_planet;
   best_fit_planet_dist;
   j2_multiplier
   debug_level
   approx_planet_orientation
   object_mass
   (implicitly) planet_posn cache
*/

double object_mass = 0.;
double j2_multiplier = 1.;
extern unsigned perturbers;
extern int n_orbit_params;

extern int debug_level;
int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
int planet_posn( const int planet_no, const double jd, double *vect_2000);
                                                /* mpc_obs.cpp */
int earth_lunar_posn( const double jd, double FAR *earth_loc, double FAR *lunar_loc);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
ldouble take_rk_stepl( const ldouble jd, ELEMENTS *ref_orbit,
                 const ldouble *ival, ldouble *ovals,
                 const int n_vals, const ldouble step);     /* runge.cpp */
ldouble take_pd89_step( const ldouble jd, ELEMENTS *ref_orbit,
                 const ldouble *ival, ldouble *ovals,
                 const int n_vals, const ldouble step);    /* runge.cpp */
int symplectic_6( double jd, ELEMENTS *ref_orbit, double *vect,
                                          const double dt);
int get_planet_posn_vel( const double jd, const int planet_no,
                     double *posn, double *vel);         /* runge.cpp */
int find_best_fit_planet( const double jd, const double *ivect,
                                 double *rel_vect);      /* runge.cpp */
int detect_perturbers( const double jd, const double * __restrict xyz,
                       double *accel);          /* bc405.cpp */
void find_relative_state_vect( const double jd, const double *ivect,
               double *ovect, const int ref_planet);        /* runge.cpp */
int find_relative_orbit( const double jd, const double *ivect,
               ELEMENTS *elements, const int ref_planet);     /* runge.cpp */
static void compute_ref_state( ELEMENTS *ref_orbit, double *ref_state,
                                          const double jd);
int parallax_to_lat_alt( const double rho_cos_phi, const double rho_sin_phi,
       double *lat, double *ht_in_meters, const int planet_idx); /* ephem0.c */
void calc_approx_planet_orientation( const int planet,        /* runge.cpp */
         const int system_number, const double jde, double *matrix);
void compute_effective_solar_multiplier( const char *constraints);   /* runge.c */
double geo_potential_in_au( const double x, const double y, const double z,
                 double *derivs, const int n_terms);    /* geo_pot.c */
double shadow_check( const double *planet_loc,           /* ephem0.cpp */
                            const double *obs_posn,
                            const double planet_radius_in_au);
double comet_g_func( const ldouble r);                   /* runge.cpp */

#define N_PERTURB 19
#define IDX_MERCURY    1
#define IDX_VENUS      2
#define IDX_EARTH      3
#define IDX_MARS       4
#define IDX_JUPITER    5
#define IDX_SATURN     6
#define IDX_URANUS     7
#define IDX_NEPTUNE    8
#define IDX_PLUTO      9
#define IDX_MOON      10
#define IDX_IO        11
#define IDX_EUROPA    12
#define IDX_GANYMEDE  13
#define IDX_CALLISTO  14
#define IDX_TETHYS    15
#define IDX_DIONE     16
#define IDX_RHEA      17
#define IDX_TITAN     18
#define IDX_IAPETUS   19
#define IDX_ASTEROIDS 20

static ldouble vector3_lengthl( const ldouble *vect)
{
   return( sqrtl( vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2]));
}

static void vector_cross_productl( ldouble *xprod, const ldouble *a, const ldouble *b)
{
   xprod[0] = a[1] * b[2] - a[2] * b[1];
   xprod[1] = a[2] * b[0] - a[0] * b[2];
   xprod[2] = a[0] * b[1] - a[1] * b[0];
}

      /* The following values for Earth J2, 3, 4 and solar J2 are from
      DE-441.  They're in units of planetary radii,  and must be converted
      to AU squared. */

#define J2_IN_EARTH_UNITS (1.08262539E-03)
#define J3_IN_EARTH_UNITS (-2.53241E-06)
#define J4_IN_EARTH_UNITS (-1.619898E-06)
#define EARTH_R2 (EARTH_RADIUS_IN_AU * EARTH_RADIUS_IN_AU)
#define EARTH_J2 (J2_IN_EARTH_UNITS * EARTH_R2)
#define EARTH_J3 (J3_IN_EARTH_UNITS * EARTH_R2 * EARTH_RADIUS_IN_AU)
#define EARTH_J4 (J4_IN_EARTH_UNITS * EARTH_R2 * EARTH_R2)


         /* From MGM1025 Mars gravity model */
#define J2_IN_MARS_UNITS (1.95545132e-3)
#define J3_IN_MARS_UNITS (3.14558811e-5)
#define J4_IN_MARS_UNITS (-1.5369417e-5)
#define MARS_R (3397.0 / AU_IN_KM)
#define MARS_J2 (J2_IN_MARS_UNITS*MARS_R * MARS_R)
#define MARS_J3 (J3_IN_MARS_UNITS*MARS_R * MARS_R * MARS_R)
#define MARS_J4 (J4_IN_MARS_UNITS*MARS_R * MARS_R * MARS_R * MARS_R)

#define MERCURY_R   (2439.4 / AU_IN_KM)
#define VENUS_R     (6051. / AU_IN_KM)
#define PLUTO_R     (1500. / AU_IN_KM)
#define MOON_R      (1748.2 / AU_IN_KM)

#ifdef NOT_USED_YET
   /* I'm in no rush to include oblateness of the Moon,  Mercury,  or Venus,
   They rotate so slowly that it's almost as if they have no oblateness
   anyway,  and I'm tracking no objects that orbit them.  But it would be
   easy to add them.  Full spherical harmonic expansions are available for
   all three objects (and the following were extracted from such expansions.)
   See 'gfc_xvt.c' from the 'miscell' project for links to models.   */

#define J2_IN_MERCURY_UNITS 6e-5

            /* From MGNP180U model,  based on Magellan data */
#define J2_IN_VENUS_UNITS   4.40443532e-6
#define J3_IN_VENUS_UNITS  -2.10819981e-6
#define J4_IN_VENUS_UNITS  -2.14742625e-6

            /* From RFM_Moon_2520  */
#define J2_IN_LUNAR_UNITS   3.8917294e-4
#define J3_IN_LUNAR_UNITS  -3.0394781e-5
#define J4_IN_LUNAR_UNITS  -9.2829444e-5

            /* From MESSENGER 100x100 gravity model */
#define J2_IN_MERCURY_UNITS      5.0354217e-5
#define J3_IN_MERCURY_UNITS      1.2587789e-5
#define J4_IN_MERCURY_UNITS      1.7600889e-5

#define MERCURY_J2   (J2_IN_MERCURY_UNITS * MERCURY_R * MERCURY_R)
#define VENUS_J2     (J2_IN_VENUS_UNITS * VENUS_R * VENUS_R)
#define MOON_J2      (J2_IN_MOON_UNITS * MOON_R * MOON_R)

            /* From DE-441.  Theoretically speaking,  the sun's oblateness
            could have noticeable effects for really precise observations.
            Maybe Gaia will do the trick?  Or Solar Orbiter or BepiColombo
            or Parker Solar Probe?        */
#define J2_IN_SOLAR_UNITS 2.1961391516529825E-07
#define SOLAR_J2 (J2_IN_SOLAR_UNITS * SUN_RADIUS_IN_AU * SUN_RADIUS_IN_AU)

#endif         /* #ifdef NOT_USED_YET */

/* Start considering atmospheric drag if within 500 km of earth: */

#define ATMOSPHERIC_LIMIT (EARTH_RADIUS_IN_AU + 500. / AU_IN_KM)
/* #define ATMOSPHERIC_LIMIT 0      */

/*  Jupiter field is now from doi:10.1038/nature25776, "Measurement of
Jupiter's asymmetric gravity field".   Saturn,  Uranus, and Neptune
fields are from http://ssd.jpl.nasa.gov/?gravity_fields_op .*/

#define J2_IN_JUPITER_UNITS (.0146956572)      /* +/- 1.4e-8 */
#define J3_IN_JUPITER_UNITS (-0.042e-6)        /* +/- 0.010e-6 */
#define J4_IN_JUPITER_UNITS (-5.86609e-4)      /* +/- 4e-8 */
#ifdef NOT_USING_ANYTHING_PAST_J4_YET
   #define J5_IN_JUPITER_UNITS   (-6.9e-8)     /* +/- 0.8e-8 */
   #define J6_IN_JUPITER_UNITS   34.198e-6     /* +/- 0.9e-8 */
   #define J7_IN_JUPITER_UNITS   1.24e-7       /* +/- 0.17e-7 */

   #define J6_IN_SATURN_UNITS    86.14e-6
   #define J8_IN_SATURN_UNITS   -10.e-6
#endif
#define JUPITER_R (71492. / AU_IN_KM)
#define JUPITER_R2 (JUPITER_R * JUPITER_R)
#define JUPITER_J2 (J2_IN_JUPITER_UNITS * JUPITER_R2)
#define JUPITER_J3 (J3_IN_JUPITER_UNITS * JUPITER_R2 * JUPITER_R)
#define JUPITER_J4 (J4_IN_JUPITER_UNITS * JUPITER_R2 * JUPITER_R2)

#define J2_IN_SATURN_UNITS (.01629071)
#define J4_IN_SATURN_UNITS (-.00093583)
#define SATURN_R (60330. / AU_IN_KM)
#define SATURN_R2 (SATURN_R * SATURN_R)
#define SATURN_J2 (J2_IN_SATURN_UNITS * SATURN_R2)
#define SATURN_J3 0.
#define SATURN_J4 (J4_IN_SATURN_UNITS * SATURN_R2 * SATURN_R2)

#define J2_IN_URANUS_UNITS  3510.68e-6
#define J4_IN_URANUS_UNITS   -34.17e-6
#define URANUS_R (25559. / AU_IN_KM)
#define URANUS_R2 (URANUS_R * URANUS_R)
#define URANUS_J2 (J2_IN_URANUS_UNITS * URANUS_R * URANUS_R)
#define URANUS_J3 0.
#define URANUS_J4 (J4_IN_URANUS_UNITS * URANUS_R2 * URANUS_R2)

#define J2_IN_NEPTUNE_UNITS  3408.43e-6
#define J4_IN_NEPTUNE_UNITS   -33.40e-6
#define NEPTUNE_R (25225. / AU_IN_KM)
#define NEPTUNE_R2 (NEPTUNE_R * NEPTUNE_R)
#define NEPTUNE_J2 (J2_IN_NEPTUNE_UNITS * NEPTUNE_R * NEPTUNE_R)
#define NEPTUNE_J3 0.
#define NEPTUNE_J4 (J4_IN_NEPTUNE_UNITS * NEPTUNE_R2 * NEPTUNE_R2)

   /* One can compute the acceleration due to the first (J2) oblateness
      term analytically,  and the following function does that.  However,
      the expressions for computing the accelerations due to J3 and J4
      get ugly,  and it becomes simpler to write a function for the _potential_
      and get the accelerations via numerical differentiation.  I'm keeping
      the numerical J2 acceleration version around just for reference. */

static double jn_potential( const double *loc, const double j3,
                                        const double j4)
{
   const double r2 = loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2];
   const double r = sqrt( r2);
   const double mu = loc[2] / r;
   const double mu2 = mu * mu;
   const double p3 = mu * (2.5 * mu2 - 1.5);
   const double p4 = (35. * mu2 * mu2 - 30. * mu2 + 3.) / 8.;
/* const double p2 = 1.5 * mu2 - .5;        Danby, p. 115 */

/*    return( (j2 * p2 + j3 * p3 / r + j4 * p4 / r2) / (r * r2));  */
      return( (          j3 * p3 / r + j4 * p4 / r2) / (r * r2));
}

/* For an input planetocentric location in AU,  and a GM in AU^3/day^2,
computes the planetocentric acceleration due to J2,  J3,  and J4,  in
AU/day^2.  For the Earth,  the GGM03 model can be used;  see 'geo_pot.cpp'.

   The number of needed terms appears to scale inversely (roughly) with
height above the earth.  I put a somewhat arbitrary limit at 250 km;
below this,  the number of terms can zoom upward with little real-world
effect,  especially if the atmosphere matters.  Inside the earth,  we
keep those three terms just to avoid a nasty numerical discontinuity,
but recognize that the results are not actually meaningful.  */

static void numerical_gradient( double *grad, const double *loc,
                   const double planet_gm,
                   const double j2, const double j3, const double j4)
{
   const double r2 = loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2];
   const double r = sqrt( r2), delta = r * 0.01;
   const double r7 = r2 * r2 * r2 * sqrt( r2);
   const double tval = (1.5 * r2 - 7.5 * loc[2] * loc[2]) / r7;
   int i;

   if( j3 == EARTH_J3)
      {
      static int n_terms = -9999;
                                    /* height in units of earth radii */
      if( n_terms == -9999)
         {
         n_terms = atoi( get_environment_ptr( "GEO_TERMS"));
         if( !n_terms)
            n_terms = 3;
         }
      if( n_terms > 0)     /* n_term <= 0 -> use "usual" J2 & J3 & J4 */
         {
         double ht_above_ground = (r / EARTH_RADIUS_IN_AU) - 1.;
         int n_added_terms = 0;

         if( ht_above_ground > 0.)
            {
            if( ht_above_ground < 0.04)         /* less than about 250 km */
               ht_above_ground = 0.04;
            n_added_terms =  (int)( (double)n_terms / ht_above_ground);
            }
         geo_potential_in_au( loc[0], loc[1], loc[2], grad, 3 + n_added_terms);
         return;
         }
      }
   grad[0] = j2 * loc[0] * tval;   /* compute J2 accel analytically */
   grad[1] = j2 * loc[1] * tval;
   grad[2] = j2 * loc[2] * (3. * r2 / r7 + tval);
   for( i = 0; i < 3; i++)         /* then do J3 & J4 numerically */
      {
      double tval, tloc[3];
      int j;

      for( j = 0; j < 3; j++)
         tloc[j] = loc[j];
      tloc[i] += delta;
      tval = jn_potential( tloc, j3, j4);
      tloc[i] -= delta + delta;
      tval -= jn_potential( tloc, j3, j4);

      grad[i] += tval / (delta * 2.);
      grad[i] *= planet_gm;
      }
#ifdef OLD_DEBUGGING_CODE
   if( j3 == EARTH_J3)
      {
      double grad2[3];

      debug_printf( "Loc (in earth radii): %e %e %e; radius %e\n",
               loc[0] / EARTH_RADIUS_IN_AU,
               loc[1] / EARTH_RADIUS_IN_AU,
               loc[2] / EARTH_RADIUS_IN_AU, r / EARTH_RADIUS_IN_AU);
      debug_printf( "Gradient (numerical): %e %e %e\n",
                  grad[0], grad[1], grad[2]);

      geo_potential_in_au( loc[0], loc[1], loc[2], grad2,
                  atoi( get_environment_ptr( "GEO_TERMS")));

      debug_printf( "Gradient (analytica): %e %e %e\n",
                  grad2[0], grad2[1], grad2[2]);
      debug_printf( "Ratios                %e %e %e\n",
                  grad[0] / grad2[0],
                  grad[1] / grad2[1],
                  grad[2] / grad2[2]);
      }
#endif
}

   /* For testing purposes,  one can multiply the relativistic effect by */
   /* a constant 'general_relativity_factor'.  Set to zero,  this disables */
   /* GR;  set to one,  it mirrors the actual universe.  I've used this    */
   /* just to see how much GR matters for a given object.   */

double general_relativity_factor = 1.;

static void set_relativistic_accel( ldouble *accel, const ldouble *posnvel)
{
   int i;
   const ldouble c = AU_PER_DAY;           /* speed of light in AU per day */
   const ldouble r_squared = posnvel[0] * posnvel[0] + posnvel[1] * posnvel[1]
                                                     + posnvel[2] * posnvel[2];
   const ldouble v_squared = posnvel[3] * posnvel[3] + posnvel[4] * posnvel[4]
                                                     + posnvel[5] * posnvel[5];
   const ldouble v_dot_r   = posnvel[0] * posnvel[3] + posnvel[1] * posnvel[4]
                                                     + posnvel[2] * posnvel[5];
   const ldouble r = sqrtl( r_squared), r_cubed_c_squared = r_squared * r * c * c;
#ifndef PREVIOUS_EQUATION
   const double r_component =
                  (4. * SOLAR_GM / r - v_squared) / r_cubed_c_squared;
   const double v_component = 4. * v_dot_r / r_cubed_c_squared;
#else
   const ldouble v_component = 3. * v_dot_r / r_cubed_c_squared;
   const ldouble r_component = 0.;
#endif

   for( i = 0; i < 3; i++)
      {
      accel[i] = r_component * posnvel[i] + v_component * posnvel[i + 3];
      accel[i] *= general_relativity_factor;
      }
}

/* This is explained on page 3 of

http://www.lpi.usra.edu/books/CometsII/7009.pdf

   'Cometary Orbit Determination and Nongravitational Forces',  D. K. Yeomans,
P. W. Chodas, G. Sitarski, S. Szutowicz, M. Krolikowska, "Comets II".

   Essentially,  the idea is that the following expression gives a decent
approximation to the magnitude of comet non-gravitational accelerations.
This 'g' function is then multiplied by a comet-specific,  fitted parameter
A1 to give the radial acceleration;  by A2 to give the transverse acceleration;
and by A3 to give the "out-of-orbit" component of the acceleration,  with A1
through A3 to be fitted parameters.  The resulting accelerations should
be in AU/day^2,  which fortunately matches the units used in this program.
(Note that MPC uses 10^-8 AU/day^2.  Patrick Rocher and JPL use AU/day^2 in
their comet elements.  Find_Orb should have a toggle for the units...)

   Frequently,  only A1 and A2 are used;  A3 often doesn't appear unless
one has three apparitions.  A3 has been added to the console version of
Find_Orb,  but it appears to be just about unnecessary.

   A small bit of trickery:  for the Marsden-Sekanina formula, g(1) = 1,
requiring a normalization constant alpha.  We get that on the first call
by setting alpha = 1,  determining what value we get for g(1),  and setting
alpha to the inverse of that.

   Also,  there may be situations in which one wants to alter the
COMET_CONSTANTS interactively.  A call to comet_g_func( 0.) resets alpha
to zero,  ensuring that the next time a "real" call is made with r > 0.,
the comet constants will be reloaded and alpha recalculated.

   Also note that a simple inverse-square force is used for space junk
and some small rocks.            */

double comet_g_func( const ldouble r)
{
   if( !is_inverse_square_force_model( ))
      {                      /* default, Marsden/Sekanina formula */
      static ldouble r0 = 2.808;         /* AU */
      static ldouble m = 2.15;
      static ldouble n = 5.093;
      static ldouble k = 4.6142;
      static ldouble alpha = 0.;
      ldouble r_over_r0;

      if( !r)                 /* resetting parameters */
         {
         alpha = 0.;
         return( 0.);
         }
      if( !alpha)
         {
         const char *comet_params = get_environment_ptr( "COMET_CONSTANTS");
         double d_r0, d_m, d_n, d_k;

         sscanf( comet_params, "%lf,%lf,%lf,%lf", &d_r0, &d_m, &d_n, &d_k);
         r0 = (ldouble)d_r0;
         m = (ldouble)d_m;
         n = (ldouble)d_n;
         k = (ldouble)d_k;
         alpha = 1.;
         alpha = 1. / comet_g_func( 1.);
         }
      r_over_r0 = r / r0;
      return( alpha * powl( r_over_r0, -m) * powl( 1. + powl( r_over_r0, n), -k));
      }
   else        /* just an inverse-square force,  used for non-comets */
      return( 1. / (r * r));
}

#ifdef NOT_CURRENTLY_USED_ALT_G_FUNCTION

/* I gather JPL is now using this comet g function,  but don't know
any details as yet.  Just putting it here for future use and as
an aide memoire.       */

static inline double new_comet_g_func( const double r)
{
   const double r_squared = r * r;
   const double tval = 1. + r * r_squared / 125.;

   return( 25. * 0.04084 * exp( -2.6 * log( tval)) / r_squared);
}
#endif            /* NOT_CURRENTLY_USED_ALT_G_FUNCTION */

/* I wrote a little code to dump the following constants from DE-432.
They were given in AU^3/day^2;  from that,  I got the following,  using
the first column to fill most of the 'planet_mass' array.  Note that for
EMB through Uranus,  these are "system" masses including moons.

       mass(obj)/mass(sun)   mass(sun)/mass(obj)    GM (km^3/s^2)
Merc 1.660114153054348e-07 6.023682155592479e+06 2.203177969072598e+04
Venu 2.447838287784771e-06 4.085237186582997e+05 3.248585874397545e+05
EMB  3.040432648022641e-06 3.289005598102476e+05 4.035032298380295e+05
Mars 3.227156037554996e-07 3.098703590290707e+06 4.282837461279101e+04
Jupi 9.547919101886966e-04 1.047348630972762e+03 1.267127623546989e+08
Satu 2.858856727222416e-04 3.497901767786634e+03 3.794058466740400e+07
Uran 4.366249662744965e-05 2.290295052370693e+04 5.794556384409937e+06
Nept 5.151383772628673e-05 1.941225977597307e+04 6.836527004611366e+06
Plut 7.350487833457740e-09 1.360453921776768e+08 9.755011621830380e+02
Sun  1.000000000000000e+00 1.000000000000000e+00 1.327124381789709e+11
Eart 3.003489614792921e-06 3.329460488475656e+05 3.986004298243866e+05
Moon 3.694303322972000e-08 2.706870315119437e+07 4.902800013642883e+03
Cere 4.725582914451237e-10 2.116141052021147e+09 6.271436303937109e+01
Pall 1.018229468217943e-10 9.820968958501590e+09 1.351317153528802e+01
Juno 1.371898570800420e-11 7.289168611179182e+10 1.820680042651692e+00
Vest 1.302666122601934e-10 7.676564106868835e+09 1.728799972636488e+01
*/

#define MASS_EARTH_PLUS_MOON 3.040432648022641e-06
#define EARTH_MOON_MASS_RATIO 81.300568800524701
#define MASS_MOON (MASS_EARTH_PLUS_MOON / (1. + EARTH_MOON_MASS_RATIO))
#define MASS_EARTH (MASS_EARTH_PLUS_MOON - MASS_MOON)

#define MASS_OF_SUN_IN_KILOGRAMS 1.9891E+30

#define MASS_IO               (893.3E+20 / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_EUROPA           (479.7E+20 / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_GANYMEDE         (1482E+20  / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_CALLISTO         (1076E+20  / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_GALILEANS  (MASS_IO + MASS_EUROPA + MASS_GANYMEDE + MASS_CALLISTO)
#define MASS_JUPITER_SYSTEM  9.547919101886966e-04
#define MASS_JUPITER   (MASS_JUPITER_SYSTEM - MASS_GALILEANS)

   /* 2015 Jun 4: revised masses for Saturn's moons to Cassini values, */
   /* and added masses for Mimas,  Enceladus,  Hyperion,  and Phoebe.  */
   /* We aren't using those yet,  but if we do,  we now have them.     */
#define MASS_SATURN_SYSTEM    2.858856727222416e-04
#define MASS_MIMAS             (0.37493E+20  / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_ENCELADUS         (1.08022E+20  / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_TETHYS            (6.17449E+20  / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_DIONE             (10.95452E+20 / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_RHEA              (23.06518E+20 / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_TITAN            (1345.2E+20 / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_HYPERION          (0.056199E+20 / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_IAPETUS           (18.05635E+20 / MASS_OF_SUN_IN_KILOGRAMS)
#define MASS_PHOEBE            (0.08292E+20  / MASS_OF_SUN_IN_KILOGRAMS)
   /* "MASS_SATURN_SATS" = mass of those we're treating as perturbers: */
#define MASS_SATURN_SATS  (MASS_TETHYS + MASS_DIONE + MASS_RHEA + MASS_TITAN + MASS_IAPETUS)
#define MASS_SATURN          (MASS_SATURN_SYSTEM - MASS_SATURN_SATS)

      /* 22 Oct 2001:  replaced Uranus,  Neptune masses w/DE-405 values */
#define EARTH_MOON_RATIO 81.30056

extern const double planet_mass[N_PERTURB + 1] = { 1.,             /*  0 */
        1.660114153054348e-07,                /* mercury */        /*  1 */
        2.447838287784771e-06,                /* venus */          /*  2 */
        MASS_EARTH,                           /* Earth */          /*  3 */
        3.227156037554996e-07,                /* Mars */           /*  4 */
        MASS_JUPITER,                         /* Jupiter */        /*  5 */
        MASS_SATURN,                          /* saturn */         /*  6 */
        4.366249662744965e-05,                /* Uranus */         /*  7 */
        5.151383772628673e-05,                /* Neptune */        /*  8 */
        7.350487833457740e-09,                /* Pluto */          /*  9 */
        MASS_MOON,                            /* Moon */           /* 10 */

         MASS_IO,                                                    /* 11 */
         MASS_EUROPA,                                                /* 12 */
         MASS_GANYMEDE,                                              /* 13 */
         MASS_CALLISTO,                                              /* 14 */

         MASS_TETHYS,                                                /* 15 */
         MASS_DIONE,                                                 /* 16 */
         MASS_RHEA,                                                  /* 17 */
         MASS_TITAN,                                                 /* 18 */
         MASS_IAPETUS };                                             /* 19 */

static ldouble planetary_system_mass( const int planet_no)
{
   ldouble rval;

   if( planet_no == 5)
      rval = MASS_JUPITER_SYSTEM;
   else if( planet_no == 6)
      rval = MASS_SATURN_SYSTEM;
   else
      rval = planet_mass[planet_no];
   return( rval);
}

/* "approx_planet_orientation" returns the planet orientation matrix,  as
described in 'cospar.cpp',  for a time no more than range/2 days away.
(Which,  at present,  means half a day.)  The idea is that the orientation
of the planet's pole doesn't change much over that time scale,
and we don't want to do a lot of expensive recomputes of planet orientation
that won't have any noticeable effect.  (It's particularly expensive for
the earth,  where we do a full-meal-deal precession and nutation matrix.)
We then rotate the planet from the cached time to the input time jde (again,
less than half a day).        */

void calc_approx_planet_orientation( const int planet,
         const int system_number, const double jde, double *matrix)
{
   static double cached_matrix[9];
   static double cached_jde = 0.;
   static int cached_planet = -1;
   const double range = 1.;
   const double new_jde = floor( jde / range + .5) * range;
   const double omega = planet_rotation_rate( planet, system_number) * PI / 180.;

   if( cached_planet != planet || new_jde != cached_jde)
      {
      const double ut = new_jde - td_minus_ut( new_jde) / seconds_per_day;

      calc_planet_orientation( planet, system_number, ut, cached_matrix);
      cached_planet = planet;
      cached_jde = new_jde;
      }
   memcpy( matrix, cached_matrix, 9 * sizeof( double));
   spin_matrix( matrix, matrix + 3, omega * (jde - cached_jde));
}

/* The idea of the following is as follows.  If our distance from the sun
is more than 120% of the planet's semimajor axis (given in the 'radii'
table),  we throw its entire mass into the sun.  (Assuming it's not
turned on as a perturber.)  If our distance is less than 100%
of the semimajor axis,  we don't add in anything.  In between,  we ramp
up the mass thrown into the sun linearly.  The reason for this is to avoid
a discontinuity,  and to instead have a gradual,  linear increase in the
mass we use for the sun.         */

static ldouble include_thrown_in_planets( const ldouble r)
{
   const int n_radii = 9;
   const ldouble radii[9] = { 0., .38709927, .72333566, 1.00000261,
               1.52371034, 5.20288799, 9.53667594,  19.18916464,  30.06992276};
   const ldouble fraction = .2;
   ldouble rval = 1.;
   int i;

   for( i = 1; i < n_radii && r > radii[i]; i++)
      if( ((perturbers >> i) & 1) == 0)
         {
         if( r > radii[i] * (1. + fraction))
            rval += planet_mass[i];
         else
            rval += planet_mass[i] * (r / radii[i] - 1.) / fraction;
         }
   return( rval);
}

/* Returns a passable estimate of the atmospheric density,  in kg/m^3,
as a function of height above sea level.  This just does an interpolation
within a table giving the atmospheric density,  at ten-kilometer intervals,
from 0 to 1000 km (extrapolating off the ends).  The table is from the
"Unofficial Australian Standard Atmosphere 2000",

http://www.sworld.com.au/steven/space/atmosphere/

   which is an update to the US Standard Atmosphere of 1976.  So far, only
the earth's atmophere is handled;  but note that similar tables are
available for Mars,  Venus,  Saturn,  and Titan :

https://solarsystem.nasa.gov/docs/03b_Justus_EDLatmospheres.pdf

   The interpolation is actually done within a table of the (natural)
logs of the original table.  That eliminates taking a logarithm at
runtime and makes the code slightly simpler.  Also,  ln(rho) is
_much_,  much closer to being a linear function than rho itself.  So
the linear interpolation is much more accurate.

   The line marked '(1)' ensures that,  for negative ht_in_km,  the
returned value extrapolates correctly,  contining to increase for a few
kilometers at roughly the correct rate (to handle objects hitting the Dead
Sea),  but then levels off and drops to a near-vacuum inside the earth.
That was the easiest way of avoiding problems with the density climbing to
absurd heights as one approached the center of the earth... doing this
avoids discontinuities that it would be awkward to code around.   */

static ldouble atmospheric_density( const ldouble ht_in_km)
{
#ifdef ORIGINAL_TABLE_FOR_REFERENCE_ONLY
   const ldouble rho[] = {          1.22500e+0,  4.13510e-1,  8.89099e-2,
         1.84102e-2,  3.99568e-3,  1.02688e-3,  3.09678e-4,  8.28286e-5,
         1.84580e-5,  3.40454e-6,  5.60400e-7,  9.39781e-8,  2.22485e-8,
         8.15200e-9,  3.81812e-9,  2.07600e-9,  1.23365e-9, 7.81486e-10,
        5.19431e-10, 3.58062e-10, 2.54100e-10, 1.84452e-10, 1.36522e-10,
        1.02784e-10, 7.85246e-11, 6.07300e-11, 4.74445e-11, 3.73990e-11,
        2.97195e-11, 2.37875e-11, 1.91600e-11, 1.55181e-11, 1.26324e-11,
        1.03320e-11, 8.48753e-12, 7.00038e-12, 5.79500e-12, 4.81307e-12,
        4.00937e-12, 3.34858e-12, 2.80300e-12, 2.35089e-12, 1.97535e-12,
        1.66279e-12, 1.40213e-12, 1.18435e-12, 1.00205e-12, 8.49163e-13,
        7.20722e-13, 6.12626e-13, 5.21500e-13, 4.44559e-13, 3.79525e-13,
        3.24498e-13, 2.77890e-13, 2.38368e-13, 2.04815e-13, 1.76297e-13,
        1.52027e-13, 1.31346e-13, 1.13700e-13, 9.86238e-14, 8.57305e-14,
        7.46934e-14, 6.52351e-14, 5.71206e-14, 5.01507e-14, 4.41564e-14,
        3.89945e-14, 3.45434e-14, 3.07000e-14, 2.73763e-14, 2.44951e-14,
        2.19914e-14, 1.98103e-14, 1.79058e-14, 1.62389e-14, 1.47768e-14,
        1.34915e-14, 1.23593e-14, 1.13600e-14, 1.04763e-14, 9.69241e-15,
        8.99505e-15, 8.37279e-15, 7.81590e-15, 7.31609e-15, 6.86621e-15,
        6.46015e-15, 6.09261e-15, 5.75900e-15, 5.45526e-15, 5.17752e-15,
        4.92236e-15, 4.68680e-15, 4.46826e-15, 4.26448e-15, 4.07349e-15,
        3.89356e-15, 3.72317e-15, 3.56100e-15  };
#endif
   const ldouble log_rho[] = { 0.2029, -0.8831, -2.4201, -3.9949,
         -5.5225, -6.8812, -8.0800, -9.3987, -10.9000, -12.5904,
        -14.3946, -16.1802, -17.6210, -18.6250, -19.3835, -19.9928,
        -20.5133, -20.9698, -21.3783, -21.7503, -22.0933, -22.4136,
        -22.7145, -22.9984, -23.2676, -23.5246, -23.7715, -24.0094,
        -24.2392, -24.4619, -24.6782, -24.8890, -25.0948, -25.2958,
        -25.4924, -25.6851, -25.8740, -26.0597, -26.2424, -26.4225,
        -26.6003, -26.7762, -26.9503, -27.1225, -27.2930, -27.4618,
        -27.6290, -27.7945, -27.9585, -28.1210, -28.2821, -28.4417,
        -28.5999, -28.7565, -28.9116, -29.0650, -29.2167, -29.3666,
        -29.5147, -29.6609, -29.8052, -29.9475, -30.0876, -30.2254,
        -30.3608, -30.4936, -30.6237, -30.7510, -30.8754, -30.9966,
        -31.1145, -31.2291, -31.3403, -31.4481, -31.5526, -31.6537,
        -31.7514, -31.8457, -31.9367, -32.0244, -32.1087, -32.1897,
        -32.2674, -32.3421, -32.4138, -32.4826, -32.5487, -32.6122,
        -32.6731, -32.7317, -32.7880, -32.8422, -32.8945, -32.9450,
        -32.9940, -33.0418, -33.0885, -33.1343, -33.1795, -33.2242,
        -33.2687 };
   ldouble fraction, log_rval;
   int range = (int)( ht_in_km / 10.);
   const int table_size = sizeof( log_rho) / sizeof( log_rho[0]);  /* = 101 */

   if( range < 0)
      range = 0;
   if( range > table_size - 2)
      range = table_size - 2;
   fraction = ht_in_km / 10. - (double)range;
   log_rval = log_rho[range] + fraction * (log_rho[range + 1] - log_rho[range])
                - (fraction < 0. ? fraction * fraction : 0.);    /* (1) */
   if( !range && log_rval < -4.)
      log_rval = -4.;
   return( expl( log_rval));
}


      /* If we're closer than .1 AU,  include Galileans separately: */

#define GALILEAN_LIMIT .03
#define TITAN_LIMIT .03

unsigned excluded_perturbers = (unsigned)-1;
int best_fit_planet;
double best_fit_planet_dist;

/* The Earth and Moon pose a special problem in the following function.
The way we want things to work is this:  if the earth's perturbations are
included,  but not the moon's,  then we pretend that there is one object,
with the combined mass of the two,  at the earth-moon barycenter.  We've
got functions (JPL DE ephemerides and VSOP) to compute the EMB location,
so this is actually pretty straightforward.

   If Earth and Moon are handled as separate objects (always the case for
objects really close to us),  things are a little stickier.  When it comes
time for the Earth's position to be computed,  we call earth_lunar_posn(),
which actually computes both the earth _and_ lunar position simultaneously.
We store the latter,  and use it when it's time to compute lunar perturbations.

   Another peculiarity that should be explained:  in the real universe,
there are limits as to how close you can get to an object (except for a
black hole),  so accelerations cannot climb toward infinity.  But objects
do hit planets,  and in the preliminary steps of orbit determination,  you
may have orbits passing through planets.  To get around this,  I added
a fictional "compute_accel_multiplier".  If you're outside the planet,
this multiplier is 1,  i.e.,  no effect.  If you're within a fraction r0
radius of the center,  deep inside the planet (or sun or moon),  this
multiplier is zero,  i.e.,  there is no force (instead of acceleration
growing to infinity).

   In between,  a cubic spline linking the boundary conditions is used.
This ensures that both the accelerations and their derivatives will
be continuous.
*/

static ldouble compute_accel_multiplier( double fraction)
{
   const ldouble r0 = .8;  /* acceleration drops to zero at 80% of planet radius */
   ldouble rval;

   assert( fraction >= 0. && fraction <= 1.);
   if( fraction < r0)
      rval = 0.;
   else
      {
      fraction = (fraction - r0) / (1. - r0);
      assert( fraction >= 0. && fraction <= 1.);
      rval =  fraction * fraction * (3. - 2. * fraction);
      }
   return( rval);
}

/* THIS HAS NOTHING TO DO WITH LIGHT-TIME LAG.

Sometimes,  comet non-gravs are modelled as having a 'lagged' effect :
the magnitude of the non-gravs is determined not by the _current_
distance from the sun,  but the distance some number of days ago.
The idea is that it may take a while for the comet to heat up and
for real non-gravs to kick in.  Given a 'current' state vector and
jd,  this computes a two-body approximate distance from the sun as
of the time jd - lag.        */

static double lagged_dist( const ldouble *state_vect, const ldouble jd,
                             const ldouble lag)
{
   double svect[6], outvect[9], rval;
   size_t i;

   for( i = 0; i < 6; i++)
      svect[i] = (double)state_vect[i];
   if( !lag)      /* usually the case */
      rval = vector3_length( svect);
   else
      {
      ELEMENTS elem;

      find_relative_orbit( (double)jd, svect, &elem, 0);
      compute_ref_state( &elem, outvect, (double)( jd - lag));
      rval = vector3_length( outvect);
      }
   return( rval);
}

#define FUDGE_FACTOR .9

static ldouble planet_radius( const int idx)
{
   return( planet_radius_in_meters( idx) * FUDGE_FACTOR / AU_IN_METERS);
}

static ldouble solar_multiplier = 1.;

/* If we have constrained an area/mass ratio,  and don't have that as
a free parameter in the orbit fit,  then the above value will be slightly
less than 1,  indicating that the object is being slightly pushed away
from the sun by solar radiation pressure.      */

void compute_effective_solar_multiplier( const char *constraints)
{
   extern int force_model;
   const char *tptr;

   solar_multiplier = 1.;
   if( constraints)
      {
      if( force_model == FORCE_MODEL_NO_NONGRAVS
               || force_model == FORCE_MODEL_DELTA_V)
         if( NULL != (tptr = strstr( constraints, "A=")))
            solar_multiplier = 1. - atof( tptr + 2) * SRP1AU / SOLAR_GM;
      }
}

int planet_hit = -1;

int calc_derivativesl( const ldouble jd, const ldouble *ival, ldouble *oval,
                           const int reference_planet)
{
   ldouble r, r2 = 0., solar_accel = 1. + object_mass;
   ldouble accel_multiplier = 1.;
   int i, j;
   unsigned local_perturbers = perturbers;
   double lunar_loc[3], jupiter_loc[3], saturn_loc[3];
   ldouble relativistic_accel[3];
   double fraction_illum = 1., ival_as_double[3];
   extern int force_model;
   static const double sphere_of_influence_radius[10] = {
            10000., 0.00075, 0.00412, 0.00618,   /* sun, mer, ven, ear */
            0.00386, 0.32229, 0.36466, 0.34606,  /* mar, jup, sat, ura */
            0.57928, 0.02208 };                  /* nep, plu */

   assert( fabsl( jd) < 1e+9);
#if !defined( _WIN32) && !defined( __APPLE__)
   for( i = 0; i < 6; i++)
      if( isnanl( ival[i]))
         {
         debug_printf( "Bad derivs; jd %Lf; ref %d\n", jd, reference_planet);
         debug_printf( "%Lf %Lf %Lf\n", ival[0], ival[1], ival[2]);
         debug_printf( "%Lf %Lf %Lf\n", ival[3], ival[4], ival[5]);
         assert( 0);
         }
   assert( !isnanl( ival[0]));
   assert( !isnanl( ival[1]));
   assert( !isnanl( ival[2]));
   assert( !isnanl( ival[3]));
   assert( !isnanl( ival[4]));
   assert( !isnanl( ival[5]));
#endif
   oval[0] = ival[3];
   oval[1] = ival[4];
   oval[2] = ival[5];
   for( i = 0; i < 3; i++)
      ival_as_double[i] = (double)ival[i];
   best_fit_planet = 0;
   planet_hit = -1;
   for( i = 0; i < 3; i++)
      r2 += ival[i] * ival[i];
   r = sqrtl( r2);
   if( n_orbit_params > 6) /* decrease non-gravs when in earth's shadow */
      {
      double earth_loc[3];

      earth_lunar_posn( jd, earth_loc, NULL);
      fraction_illum = shadow_check( earth_loc, ival_as_double, EARTH_RADIUS_IN_AU);
      }
   if( force_model == FORCE_MODEL_SRP)
      solar_accel -= ival[6] * fraction_illum;
   if( r < planet_radius( 0))     /* special fudge to keep acceleration from reaching */
      {                          /* infinity inside the sun;  see above notes        */
      accel_multiplier = compute_accel_multiplier( r / planet_radius( 0));
      planet_hit = 0;
      if( debug_level)
         debug_printf( "Inside the sun: %f km\n", (double)r * AU_IN_KM);
      if( !accel_multiplier)
         {
         for( i = 3; i < 6; i++)
            oval[i] = 0.;
         return( 0);
         }
      }

   solar_accel *= -SOLAR_GM / (r2 * r);
   solar_accel *= solar_multiplier;

   if( perturbers)
      set_relativistic_accel( relativistic_accel, ival);
   else                           /* shut off relativity if no perturbers */
      for( i = 0; i < 3; i++)
         relativistic_accel[i] = 0.;

   solar_accel *= include_thrown_in_planets( r);

   for( i = 0; i < 3; i++)
      oval[i + 3] = solar_accel * ival[i]
                 + SOLAR_GM * relativistic_accel[i];

   if( (local_perturbers >> IDX_ASTEROIDS) & 1)
      if( r < 11.5 && r > 1.)
         {
         double asteroid_accel[6];

         for( i = 3; i < 6; i++)
            asteroid_accel[i] = 0.;
         detect_perturbers( jd, ival_as_double, asteroid_accel);        /* bc405.cpp */
         for( i = 3; i < 6; i++)
            oval[i] += asteroid_accel[i];
         }

   if( (n_orbit_params >= 8 && n_orbit_params <= 10 && force_model != FORCE_MODEL_DELTA_V)
                              || force_model == FORCE_MODEL_YARKO_A2)
      {                  /* Marsden & Sekanina comet formula */
      const ldouble lag = (n_orbit_params == 10 ? ival[9] : 0.);
      const ldouble g = comet_g_func( lagged_dist( ival, jd, lag)) * fraction_illum;
      ldouble transverse[3], dot_prod = 0.;

#if !defined( _WIN32) && !defined( __APPLE__)
      assert( !isnanl( g));
#endif
      memcpy( transverse, ival + 3, 3 * sizeof( ldouble));
      for( i = 0; i < 3; i++)
         dot_prod += transverse[i] * ival[i];
      for( i = 0; i < 3; i++)
         transverse[i] -= ival[i] * dot_prod / r;
      dot_prod = vector3_lengthl( transverse);
      if( force_model == FORCE_MODEL_YARKO_A2)
         for( i = 0; i < 3; i++)
            oval[i + 3] += g * (ival[6] * transverse[i] / dot_prod);
      else
         for( i = 0; i < 3; i++)
            oval[i + 3] += g * (ival[6] * ival[i] / r
                     + ival[7] * transverse[i] / dot_prod);
      if( n_orbit_params >= 9)
         {
         ldouble out_of_plane[3];

         vector_cross_productl( out_of_plane, ival, transverse);
         dot_prod = vector3_lengthl( out_of_plane);
         for( i = 0; i < 3; i++)
            oval[i + 3] += g * ival[8] * out_of_plane[i] / dot_prod;
         }
      }
   for( i = 0; i < 3; i++)       /* redundant initialization */
      jupiter_loc[i] = 0.;       /* to avoid gcc-13 warning  */

   if( perturbers)
      for( i = 1; i < N_PERTURB + 1; i++)
         if( ((local_perturbers >> i) & 1)
                   && !((excluded_perturbers >> i) & 1))
            {
            double planet_loc[15], accel[3], mass_to_use = planet_mass[i];

            r = r2 = 0.;
            if( i >= IDX_IO)       /* Galileans,  Titan */
               {
               double matrix[10], sat_loc[15];
               const double t_years = (jd - J2000) / 365.25;

               if( i >= IDX_TETHYS)         /* Saturnian satell */
                  calc_ssat_loc( jd, sat_loc,
                               ((i == IDX_IAPETUS) ? 7 : i - 13), 0L);
               else
                  {
                  calc_jsat_loc( jd, sat_loc, 1 << (i - IDX_IO), 0L);
                  memmove( sat_loc, sat_loc + (i - IDX_IO) * 3,
                                                      3 * sizeof( double));
                  }
                                 /* turn ecliptic of date to equatorial: */
               rotate_vector( sat_loc, mean_obliquity( t_years / 100.), 0);
                                 /* then to equatorial J2000: */
               setup_precession( matrix, 2000. + t_years, 2000.);
               precess_vector( matrix, sat_loc, planet_loc + 12);
                                 /* then to ecliptic J2000: */
               equatorial_to_ecliptic( planet_loc + 12);
               for( j = 0; j < 3; j++)
                  {
                  double coord;

                  if( i >= IDX_TETHYS)         /* Saturnian */
                     coord = saturn_loc[j] + planet_loc[12 + j];
                  else
                     coord = jupiter_loc[j] + planet_loc[12 + j] * JUPITER_R;
                  r2 += coord * coord;
                  planet_loc[12 + j] = coord;
                  }
               planet_loc[2] = sqrt( r2);
               }
            else
               {
               if( local_perturbers & 1024)   /* if the moon is included */
                  {
                  if( i == 3)
                     earth_lunar_posn( jd, planet_loc, lunar_loc);
                  else if( i == 10)
                     memcpy( planet_loc, lunar_loc, 3 * sizeof( double));
                  else
                     planet_posn( i, jd, planet_loc);
                  }
               else
                  planet_posn( i, jd, planet_loc);

               for( j = 0; j < 3; j++)
                  r2 += planet_loc[j] * planet_loc[j];
               memcpy( planet_loc + 12, planet_loc, 3 * sizeof( double));
               planet_loc[2] = sqrt( r2);
               }

            for( j = 0; j < 3; j++)
               {
               accel[j] = ival[j] - planet_loc[12 + j];
               r += accel[j] * accel[j];
               }
            r = sqrt( r);

            if( i == IDX_JUPITER)
               {
               memcpy( jupiter_loc, planet_loc + 12, 3 * sizeof( double));
               if( r < GALILEAN_LIMIT)
                  local_perturbers |= (15 << 11);
               else     /* "throw" Galileans into Jupiter: */
                  mass_to_use = MASS_JUPITER_SYSTEM;
               }

            if( i == IDX_SATURN)
               {
               memcpy( saturn_loc, planet_loc + 12, 3 * sizeof( double));
               if( r < TITAN_LIMIT)
                  local_perturbers |= (31 << 15);
               else        /* "throw" saturn's satellites into the primary: */
                  mass_to_use = MASS_SATURN_SYSTEM;
               }

            if( r < planet_radius( i))
               {
               accel_multiplier = compute_accel_multiplier( r / planet_radius( i));
               planet_hit = i;
               if( accel_multiplier == 0.)
                  {
                  for( i = 3; i < 6; i++)
                     oval[i] = 0.;
                  return( planet_hit);
                  }
               }
            if( i >= IDX_EARTH && i <= IDX_NEPTUNE && r < .015 && j2_multiplier)
               {          /* Within .015 AU,  we take J2 into account: */
               double grad[3], delta_j2000[3], matrix[10], delta_planet[3];
               const double j2[6] = { EARTH_J2, MARS_J2, JUPITER_J2,
                        SATURN_J2, URANUS_J2, NEPTUNE_J2 };
               const double j3[6] = { EARTH_J3, MARS_J3, JUPITER_J3,
                        SATURN_J3, URANUS_J3, NEPTUNE_J3 };
               const double j4[6] = { EARTH_J4, MARS_J4, JUPITER_J4,
                        SATURN_J4, URANUS_J4, NEPTUNE_J4 };

               calc_approx_planet_orientation( i, 0, jd, matrix);
                           /* Remembering the 'accels' are 'deltas' now... */
               memcpy( delta_j2000, accel, 3 * sizeof( double));
                           /* Cvt ecliptic to equatorial 2000...: */
               ecliptic_to_equatorial( delta_j2000);
                           /* ...then to planet-centered coords: */
               precess_vector( matrix, delta_j2000, delta_planet);

               if( j2_multiplier)
                  numerical_gradient( grad, delta_planet,
                        planet_mass[i] * SOLAR_GM,
                        j2[i - 3], j3[i - 3], j4[i - 3]);
               else     /* inside the planet */
                  for( j = 0; j < 3; j++)
                     grad[j] = 0.;
                        /* Cvt gradient from planet-centric to equatorial: */
               deprecess_vector( matrix, grad, delta_j2000);
                           /* Cvt equatorial to ecliptic: */
               equatorial_to_ecliptic( delta_j2000);
                           /* And add 'em to the output acceleration: */
               for( j = 0; j < 3; j++)
                  oval[j + 3] -= j2_multiplier * delta_j2000[j];
               if( i == IDX_EARTH && r < ATMOSPHERIC_LIMIT && n_orbit_params == 7
                           && *get_environment_ptr( "DRAG_SHUTOFF") != '1')
                  {
                  const double amr_drag = ival[6] * SOLAR_GM / SRP1AU;
                  const double rho_cos_phi =
                          hypot( delta_planet[0], delta_planet[1]) / EARTH_RADIUS_IN_AU;
                  const double rho_sin_phi = delta_planet[2] / EARTH_RADIUS_IN_AU;
                  double ht_in_meters, rho;
                  const double meters_per_km = 1000.;
                  const double dt = .5 / minutes_per_day;  /* half a minute */
                  const double Cd = 0.47;  /* drag coeff for a sphere */
                  double earth_loc[3];  /* 2nd location of earth, used for getting vel */
                  double vel[3];       /* velocity of obj relative to the earth */
                  double speed;        /* magnitude of the vel[] vector */
                  double drag[3];     /* drag acceleration,  in m/s^2 */
                  double accel_coeff;

                  parallax_to_lat_alt( rho_cos_phi, rho_sin_phi, NULL,
                                    &ht_in_meters, i);
                  rho = atmospheric_density( ht_in_meters / meters_per_km);
                  earth_lunar_posn( jd + dt, earth_loc, NULL);
                  for( j = 0; j < 3; j++)
                     {
                                    /* in degrees/day,  sidereal... */
                     const double earth_rotation_rate = 360.98564736629;
                                    /* ...and converted to radians/day:    */
                     const double earth_omega = earth_rotation_rate * (PI / 180.);
                     const double topo_vel = earth_omega *
                             (delta_planet[0] * matrix[j + 3] - delta_planet[1] * matrix[j]);

                     vel[j] = -topo_vel * accel_multiplier;
                     }
                  equatorial_to_ecliptic( vel);
                  for( j = 0; j < 3; j++)
                     {
                     const double earth_vel = (earth_loc[j] - planet_loc[j + 12]) / dt;

                     vel[j] += earth_vel - oval[j];       /* in AU/day */
                     vel[j] *= AU_IN_METERS / seconds_per_day;  /* cvt to m/s */
                     }
                  speed = vector3_length( vel);    /* also in m/s */
                  accel_coeff = rho * Cd * speed * amr_drag / 2.;
                  for( j = 0; j < 3; j++)
                     drag[j] = vel[j] * accel_coeff;
                  for( j = 0; j < 3; j++)
                     oval[j + 3] += drag[j] * seconds_per_day * seconds_per_day / AU_IN_METERS;
                  }
               }

                     /* If we're including the earth,  but the moon isn't */
                     /* being included separately,  add the mass of the moon:*/
            if( i == IDX_EARTH)
               if( !((local_perturbers >> IDX_MOON) & 1) ||
                    ((excluded_perturbers >> IDX_MOON) & 1))
                  mass_to_use += planet_mass[IDX_MOON];


            if( i < 10 && r < sphere_of_influence_radius[i]
                                && reference_planet != -1)
               {
               best_fit_planet = i;
               best_fit_planet_dist = r;
               }

/*          if( accel_multiplier)  */
               {
               const ldouble accel_factor =
                               -SOLAR_GM * mass_to_use / (r * r * r);

               for( j = 0; j < 3; j++)
                  oval[j + 3] += accel_factor * accel[j];
               }
            if( i != reference_planet)
               {
               if( reference_planet >= 0)
                  {
                  double planet_posn[3];

                  get_planet_posn_vel( jd, reference_planet, planet_posn, NULL);
                  for( j = 0; j < 3; j++)
                     planet_loc[j + 12] -= planet_posn[j];
                  r = vector3_length( planet_loc + 12);
                  }
               else
                  r = planet_loc[2];
               if( r)
                   r = -SOLAR_GM * mass_to_use / (r * r * r);
               for( j = 0; j < 3; j++)
                  oval[j + 3] += r * planet_loc[j + 12];
               }
            else        /* subtract sun's attraction to reference planet */
               {
               r = planet_loc[2];
               if( r)
                   r = SOLAR_GM / (r * r * r);
               for( j = 0; j < 3; j++)
                  oval[j + 3] += r * planet_loc[j + 12];
               }
            }
   if( planet_hit != -1)
      for( j = 3; j < 6; j++)
         oval[j] *= accel_multiplier;
   return( planet_hit);
}

int calc_derivatives( const double jd, const double *ival, double *oval,
                           const int reference_planet)
{
   size_t i;
   int rval;
   ldouble ival1[MAX_N_PARAMS], oval1[MAX_N_PARAMS];

   assert( fabs( (double)jd) < 1e+9);
   assert( n_orbit_params >= 6);
   for( i = 0; i < (size_t)n_orbit_params; i++)
      ival1[i] = (ldouble)ival[i];
   rval = calc_derivativesl( (ldouble)jd, ival1, oval1, reference_planet);
   for( i = 0; i < (size_t)n_orbit_params; i++)
      oval[i] = (double)oval1[i];
   return( rval);
}

int get_planet_posn_vel( const double jd, const int planet_no,
                     double *posn, double *vel)
{
   assert( fabs( jd) < 1e+9);
   if( posn)
      {
      if( !planet_no)       /* sun doesn't move in the heliocentric frame */
         memset( posn, 0, 3 * sizeof( double));
      else if( planet_no == 3)
         earth_lunar_posn( jd, posn, NULL);
      else if( planet_no == 10)
         earth_lunar_posn( jd, NULL, posn);
      else
         planet_posn( planet_no, jd, posn);
      }
   if( vel)
      {
      if( !planet_no)       /* sun doesn't move in the heliocentric frame */
         memset( vel, 0, 3 * sizeof( double));
      else
         {
         int i;
         double loc1[3], loc2[3];
         const double delta = 1. / minutes_per_day;    /* one minute delta... */

         get_planet_posn_vel( jd + delta, planet_no, loc2, NULL);
         get_planet_posn_vel( jd - delta, planet_no, loc1, NULL);
         for( i = 0; i < 3; i++)
            vel[i] = (loc2[i] - loc1[i]) / (2. * delta);
         }
      }
   return( 0);
}

void find_relative_state_vect( const double jd, const double *ivect,
               double *ovect, const int ref_planet)
{
   assert( ref_planet >= 0);
   memcpy( ovect, ivect, n_orbit_params * sizeof( double));
   if( ref_planet)
      {
      double planet_state[6];
      size_t i;

      get_planet_posn_vel( jd, ref_planet, planet_state, planet_state + 3);
      for( i = 0; i < 6; i++)
         ovect[i] -= planet_state[i];
      }
}

int find_relative_orbit( const double jd, const double *ivect,
               ELEMENTS *elements, const int ref_planet)
{
   double local_rel_vect[MAX_N_PARAMS];

   assert( ref_planet >= 0);
   assert( elements);
   find_relative_state_vect( jd, ivect, local_rel_vect, ref_planet);
   if( elements)
      {
      elements->gm = SOLAR_GM * planetary_system_mass( ref_planet);
      calc_classical_elements( elements, local_rel_vect, jd, 1);
      elements->epoch = jd;
      elements->central_obj = ref_planet;
      }
   return( 0);
}

int check_for_perturbers( const double t_cen, const double *vect); /* sm_vsop*/

int find_best_fit_planet( const double jd, const double *ivect,
                                 double *rel_vect)
{
   int i, j, rval = 0;
   double best_fit = 0., vel[3];
   int included = 1;          /* the Sun is always included in this */
   const int possible_perturber = check_for_perturbers( (jd - J2000) / 36525.,
                                             ivect);

   assert( ivect);
   assert( rel_vect);
   included |= (1 << possible_perturber);
   if( possible_perturber == 3)
      included |= (1 << 10);
   for( i = 0; i < 11; i++)
      if( (included >> i) & 1)
         {
         double planet_loc[3], delta[3], dist2 = 0., curr_fit;

         get_planet_posn_vel( jd, i, planet_loc, NULL);
         for( j = 0; j < 3; j++)
            {
            delta[j] = ivect[j] - planet_loc[j];
            dist2 += delta[j] * delta[j];
            }

         curr_fit = planet_mass[i] / exp( log( dist2) * 1.25);
         if( !i || curr_fit > best_fit)
            {
            rval = i;
            best_fit = curr_fit;
            memcpy( rel_vect, delta, 3 * sizeof( double));
            }
         }

   get_planet_posn_vel( jd, rval, NULL, vel);
   for( j = 0; j < 3; j++)
      rel_vect[j + 3] = ivect[j + 3] - vel[j];
   return( rval);
}

static void compute_ref_state( ELEMENTS *ref_orbit, double *ref_state,
                                          const double jd)
{
   double r2 = 0., accel;
   int i;

   if( ref_orbit->central_obj < 0)      /* using the method of Cowell */
      {
      memset( ref_state, 0, 9 * sizeof( double));
      return;
      }

   comet_posn_and_vel( ref_orbit, jd, ref_state, ref_state + 3);
   for( i = 0; i < 3; i++)
      r2 += ref_state[i] * ref_state[i];
   accel = -SOLAR_GM * planetary_system_mass( ref_orbit->central_obj) / (r2 * sqrt( r2));
   for( i = 0; i < 3; i++)
      ref_state[i + 6] = accel * ref_state[i];

   if( ref_orbit->central_obj)
      {
      double planet_state[6];

      get_planet_posn_vel( jd, ref_orbit->central_obj, planet_state, planet_state + 3);
      for( i = 0; i < 6; i++)
         ref_state[i] += planet_state[i];
      }
}

#define B_2_1     .0555555555555555555555555555556
#define B_3_1     .0208333333333333333333333333333
#define B_3_2     .0625
#define B_4_1     .03125
#define B_4_2    0.
#define B_4_3     .09375
#define B_5_1     .3125
#define B_5_2    0.
#define B_5_3   -1.171875
#define B_5_4    1.171875
#define B_6_1     .0375
#define B_6_2    0.
#define B_6_3    0.
#define B_6_4     .1875
#define B_6_5     .15
#define B_7_1     .0479101371111111111111111111111
#define B_7_2    0.
#define B_7_3    0.
#define B_7_4     .112248712777777777777777777778
#define B_7_5    -.0255056737777777777777777777778
#define B_7_6     .0128468238888888888888888888889
#define B_8_1        .016917989787292281181431107136
#define B_8_2       0.
#define B_8_3       0.
#define B_8_4       .387848278486043169526545744159
#define B_8_5       .0359773698515003278967008896348
#define B_8_6       .196970214215666060156715256072
#define B_8_7      -.172713852340501838761392997002
#define B_9_1       .0690957533591923006485645489846
#define B_9_2      0.
#define B_9_3      0.
#define B_9_4      -.634247976728854151882807874972
#define B_9_5      -.161197575224604080366876923982
#define B_9_6       .138650309458825255419866950133
#define B_9_7       .94092861403575626972423968413
#define B_9_8       .211636326481943981855372117132
#define B_10_1      .183556996839045385489806023537
#define B_10_2     0.
#define B_10_3     0.
#define B_10_4    -2.46876808431559245274431575997
#define B_10_5     -.291286887816300456388002572804
#define B_10_6     -.026473020233117375688439799466
#define B_10_7     2.84783876419280044916451825422
#define B_10_8      .281387331469849792539403641827
#define B_10_9      .123744899863314657627030212664
#define B_11_1    -1.21542481739588805916051052503
#define B_11_2     0.
#define B_11_3     0.
#define B_11_4    16.6726086659457724322804132886
#define B_11_5      .915741828416817960595718650451
#define B_11_6    -6.05660580435747094755450554309
#define B_11_7   -16.0035735941561781118417064101
#define B_11_8    14.849303086297662557545391898
#define B_11_9   -13.3715757352898493182930413962
#define B_11_10    5.13418264817963793317325361166
#define B_12_1      .258860916438264283815730932232
#define B_12_2     0.
#define B_12_3     0.
#define B_12_4    -4.77448578548920511231011750971
#define B_12_5     -.43509301377703250944070041181
#define B_12_6    -3.04948333207224150956051286631
#define B_12_7     5.57792003993609911742367663447
#define B_12_8     6.15583158986104009733868912669
#define B_12_9    -5.06210458673693837007740643391
#define B_12_10    2.19392617318067906127491429047
#define B_12_11     .134627998659334941535726237887
#define B_13_1      .822427599626507477963168204773
#define B_13_2     0.
#define B_13_3     0.
#define B_13_4   -11.6586732572776642839765530355
#define B_13_5     -.757622116690936195881116154088
#define B_13_6      .713973588159581527978269282765
#define B_13_7    12.0757749868900567395661704486
#define B_13_8    -2.12765911392040265639082085897
#define B_13_9     1.99016620704895541832807169835
#define B_13_10    -.234286471544040292660294691857
#define B_13_11     .17589857770794226507310510589
#define B_13_12    0.

/*  The coefficients BHAT(*) refer to the formula used to advance the
   integration, here the one of order 8.  The coefficients B(*) refer
   to the other formula, here the one of order 7.  */

#define CHAT_1     .0417474911415302462220859284685
#define CHAT_2    0.
#define CHAT_3    0.
#define CHAT_4    0.
#define CHAT_5    0.
#define CHAT_6    -.0554523286112393089615218946547
#define CHAT_7     .239312807201180097046747354249
#define CHAT_8     .70351066940344302305804641089
#define CHAT_9    -.759759613814460929884487677085
#define CHAT_10    .660563030922286341461378594838
#define CHAT_11    .158187482510123335529614838601
#define CHAT_12   -.238109538752862804471863555306
#define CHAT_13    .25

#define C_1      .029553213676353496981964883112
#define C_2     0.
#define C_3     0.
#define C_4     0.
#define C_5     0.
#define C_6     -.828606276487797039766805612689
#define C_7      .311240900051118327929913751627
#define C_8     2.46734519059988698196468570407
#define C_9    -2.54694165184190873912738007542
#define C_10    1.44354858367677524030187495069
#define C_11     .0794155958811272872713019541622
#define C_12     .0444444444444444444444444444445
#define C_13    0.

#define A_1     0.
#define A_2      .0555555555555555555555555555556
#define A_3      .0833333333333333333333333333334
#define A_4      .125
#define A_5      .3125
#define A_6      .375
#define A_7      .1475
#define A_8      .465
#define A_9      .564865451382259575398358501426
#define A_10     .65
#define A_11     .924656277640504446745013574318
#define A_12    1.
#define A_13      A_12

#define N_EVALS 13
#define N_EVALS_PLUS_ONE 14

ldouble take_pd89_step( const ldouble jd, ELEMENTS *ref_orbit,
                 const ldouble *ival, ldouble *ovals,
                 const int n_vals, const ldouble step)
{
   ldouble *ivals[N_EVALS_PLUS_ONE], *ivals_p[N_EVALS], rval = 0.;
   int i, j, k;
   const ldouble bvals[91] = { B_2_1,
       B_3_1, B_3_2,
       B_4_1, B_4_2, B_4_3,
       B_5_1, B_5_2, B_5_3, B_5_4,
       B_6_1, B_6_2, B_6_3, B_6_4, B_6_5,
       B_7_1, B_7_2, B_7_3, B_7_4, B_7_5, B_7_6,
       B_8_1, B_8_2, B_8_3, B_8_4, B_8_5, B_8_6, B_8_7,
       B_9_1, B_9_2, B_9_3, B_9_4, B_9_5, B_9_6, B_9_7, B_9_8,
       B_10_1, B_10_2, B_10_3, B_10_4, B_10_5, B_10_6, B_10_7, B_10_8, B_10_9,
       B_11_1, B_11_2, B_11_3, B_11_4, B_11_5, B_11_6, B_11_7, B_11_8, B_11_9, B_11_10,
       B_12_1, B_12_2, B_12_3, B_12_4, B_12_5, B_12_6, B_12_7, B_12_8, B_12_9, B_12_10, B_12_11,
       B_13_1, B_13_2, B_13_3, B_13_4, B_13_5, B_13_6, B_13_7, B_13_8, B_13_9, B_13_10, B_13_11, B_13_12,
       CHAT_1, CHAT_2, CHAT_3, CHAT_4, CHAT_5, CHAT_6, CHAT_7, CHAT_8, CHAT_9, CHAT_10, CHAT_11, CHAT_12, CHAT_13 };
   const ldouble avals[N_EVALS_PLUS_ONE] = { 0, A_1, A_2, A_3, A_4, A_5,
             A_6, A_7, A_8, A_9, A_10, A_11, A_12, A_13 };
   const ldouble *bptr = bvals;

   ivals[0] = (ldouble *)calloc( (2 * N_EVALS + 1) * n_vals, sizeof( ldouble));
   assert( ivals[0]);
   if( !ivals[0])
      return( 0.);
   for( i = 0; i < N_EVALS; i++)
      {
      ivals[i + 1] = ivals[0] + (i + 1) * n_vals;
      ivals_p[i] = ivals[0] + (i + N_EVALS + 1) * n_vals;
      }

   for( j = 0; j <= N_EVALS; j++)
      {
      ldouble ref_state_j[9], state_j[MAX_N_PARAMS];
      const ldouble jd_j = jd + step * avals[j];
      double temp_array[9];

      compute_ref_state( ref_orbit, temp_array, jd_j);
      for( i = 0; i < 9; i++)
         ref_state_j[i] = (ldouble)temp_array[i];
      if( !j)
         {
         memcpy( state_j, ival, n_orbit_params * sizeof( ldouble));
               /* subtract the analytic posn/vel from the numeric: */
         for( i = 0; i < n_vals; i++)
            ivals[0][i] = ival[i] - ref_state_j[i];
         }
      else
         for( i = 0; i < n_vals; i++)
            {
            ldouble tval = 0.;

            for( k = 0; k < j; k++)
               tval += bptr[k] * ivals_p[k][i];
            ivals[j][i] = tval * step + ivals[0][i];
            state_j[i] = ivals[j][i] + ref_state_j[i];
            }
      bptr += j;
      if( j != N_EVALS)
         {
         assert( fabsl( jd_j) < 1e+9);
         calc_derivativesl( jd_j, state_j, ivals_p[j], ref_orbit->central_obj);
         for( k = 0; k < 6; k++)
            ivals_p[j][k] -= ref_state_j[k + 3];
         }
      else     /* on last iteration,  we have our answer: */
         memcpy( ovals, state_j, n_orbit_params * sizeof( ldouble));
      }

   for( i = 0; i < n_vals; i++)
      {
      ldouble tval = 0.;
      const ldouble err_coeff[N_EVALS] = { CHAT_1 - C_1, CHAT_2 - C_2,
               CHAT_3  -  C_3, CHAT_4  -  C_4, CHAT_5  -  C_5, CHAT_6 - C_6,
               CHAT_7  -  C_7, CHAT_8  -  C_8, CHAT_9  -  C_9, CHAT_10 - C_10,
               CHAT_11 - C_11, CHAT_12 - C_12, CHAT_13 - C_13 };

      for( k = 0; k < N_EVALS; k++)
         tval += err_coeff[k] * ivals_p[k][i];
      rval += tval * tval;
      }
   return( sqrtl( rval * step * step));
}

#define ORIGINAL_FEHLBERG_CONSTANTS

         /* These "original" constants can be found in Danby, p. 298.  */
         /* The ones actually used are from _Numerical Recipes_,  and  */
         /* are for the Cash-Karp variant of Runge-Kutta.  NR says      */
         /* these constants are slightly better.  (Must admit,  I've not */
         /* done a really careful comparison!)  Dormand-Prince's method */
         /* (the 5/4 order version,  not the 9/8 one given above) may also */
         /* be worth a look.   */
#ifdef ORIGINAL_FEHLBERG_CONSTANTS
#define RKF_B21 2. / 9.
#define RKF_B31 1. / 12.
#define RKF_B32 1. / 4.
#define RKF_B41 69. / 128.
#define RKF_B42 -243. / 128.
#define RKF_B43 135. / 64.
#define RKF_B51 -17. / 12.
#define RKF_B52 27. / 4.
#define RKF_B53 -27. / 5.
#define RKF_B54 16. / 15.
#define RKF_B61 65. / 432.
#define RKF_B62 -5. / 16.
#define RKF_B63 13 / 16.
#define RKF_B64 4 / 27.
#define RKF_B65 5. / 144.
#define RKF_CHAT1 47. / 450.
#define RKF_CHAT2 0.
#define RKF_CHAT3 12 / 25.
#define RKF_CHAT4 32. / 225.
#define RKF_CHAT5 1. / 30.
#define RKF_CHAT6 6. / 25.
#define RKF_C1  1. / 9.
#define RKF_C2 0.
#define RKF_C3 9. / 20.
#define RKF_C4 16. / 45.
#define RKF_C5 1. / 12.
#define RKF_C6 0.
#define RKF_A1       0.
#define RKF_A2       2. / 9.
#define RKF_A3       1. / 3.
#define RKF_A4       .75
#define RKF_A5       1.
#define RKF_A6       5. / 6.
#else
#define RKF_B21     1. / 5.
#define RKF_B31     3. / 40.
#define RKF_B32     9. / 40.
#define RKF_B41     3. / 10.
#define RKF_B42    -9. / 10.
#define RKF_B43     6. / 5.
#define RKF_B51   -11. / 54.
#define RKF_B52     5. / 2.
#define RKF_B53   -70. / 27.
#define RKF_B54    35. / 27.
#define RKF_B61  1631. / 55296
#define RKF_B62   175. / 512.
#define RKF_B63   575. / 13824.
#define RKF_B64 44275. / 110592.
#define RKF_B65   253. / 4096.
#define RKF_CHAT1  2825. / 27648.
#define RKF_CHAT2 0.
#define RKF_CHAT3 18575. / 48384.
#define RKF_CHAT4 13525. / 55296.
#define RKF_CHAT5 277. / 14336.
#define RKF_CHAT6 .25
#define RKF_C1    37. / 378.
#define RKF_C2            0.
#define RKF_C3   250. / 621.
#define RKF_C4   125. / 594.
#define RKF_C5            0.
#define RKF_C6  512. / 1771.
#define RKF_A1           0.
#define RKF_A2            .2
#define RKF_A3            .3
#define RKF_A4            .6
#define RKF_A5           1.
#define RKF_A6            .875
#endif
/*   Butcher tableau looks like this :
A1 |
A2 | B21
A3 | B31 B32
A4 | B41 B42 B43
A5 | B51 B52 B53 B54
A6 | B61 B62 B63 B64 B65
---+---------------------
   | C1  C2  C3  C4  C5  C6
   | C^1 C^2 C^3 C^4 C^5 C^6          */

ldouble take_rk_stepl( const ldouble jd, ELEMENTS *ref_orbit,
                 const ldouble *ival, ldouble *ovals,
                 const int n_vals, const ldouble step)
{
   ldouble *ivals[7], *ivals_p[6], rval = 0.;
   int i, j, k;
            /* Revised values from _Numerical Recipes_: */
   const ldouble bvals[21] = { RKF_B21,
            RKF_B31, RKF_B32,
            RKF_B41, RKF_B42, RKF_B43,
            RKF_B51, RKF_B52, RKF_B53, RKF_B54,
            RKF_B61, RKF_B62, RKF_B63, RKF_B64, RKF_B65,
            RKF_CHAT1, RKF_CHAT2, RKF_CHAT3,
            RKF_CHAT4, RKF_CHAT5, RKF_CHAT6 };

   const ldouble avals[7] = { RKF_A1, RKF_A2, RKF_A3, RKF_A4, RKF_A5, RKF_A6, 1.};
   const ldouble *bptr = bvals;
   ldouble temp_ivals[78];

   if( n_vals > 6)
      ivals[0] = (ldouble *)calloc( 13 * n_vals, sizeof( ldouble));
   else
      ivals[0] = temp_ivals;
   assert( ivals[0]);
   if( !ivals[0])
      return( 0.);
   for( i = 0; i < 6; i++)
      {
      ivals[i + 1] = ivals[0] + (i + 1) * n_vals;
      ivals_p[i] = ivals[0] + (i + 7) * n_vals;
      }

   for( j = 0; j < 7; j++)
      {
      ldouble ref_state_j[9], state_j[MAX_N_PARAMS];
      const ldouble jd_j = jd + step * avals[j];
      double temp_array[9];

      compute_ref_state( ref_orbit, temp_array, (double)jd_j);
      for( i = 0; i < 9; i++)
         ref_state_j[i] = (ldouble)temp_array[i];
      if( !j)
         {
         memcpy( state_j, ival, n_orbit_params * sizeof( ldouble));
               /* subtract the analytic posn/vel from the numeric: */
         for( i = 0; i < n_vals; i++)
            ivals[0][i] = ival[i] - ref_state_j[i];
         }
      else
         for( i = 0; i < n_vals; i++)
            {
            ldouble tval = 0.;

            for( k = 0; k < j; k++)
               tval += bptr[k] * ivals_p[k][i];
            ivals[j][i] = tval * step + ivals[0][i];
            state_j[i] = ivals[j][i] + ref_state_j[i];
            }
      bptr += j;
      if( j != 6)
         {
#ifndef __WATCOMC__
         assert( fabsl( jd_j) < 1e+9);
#endif
         calc_derivativesl( jd_j, state_j, ivals_p[j], ref_orbit->central_obj);
         for( k = 0; k < 6; k++)
            ivals_p[j][k] -= ref_state_j[k + 3];
         }
      else     /* on last iteration,  we have our answer: */
         memcpy( ovals, state_j, n_orbit_params * sizeof( ldouble));
      }

   for( i = 0; i < n_vals; i++)
      {
      ldouble tval = 0.;
      static const ldouble err_coeffs[6] = {
            RKF_CHAT1 - RKF_C1, RKF_CHAT2 - RKF_C2, RKF_CHAT3 - RKF_C3,
            RKF_CHAT4 - RKF_C4, RKF_CHAT5 - RKF_C5, RKF_CHAT6 - RKF_C6 };

      for( k = 0; k < 6; k++)
         tval += err_coeffs[k] * ivals_p[k][i];
      rval += tval * tval;
      }

   if( n_vals > 6)
      free( ivals[0]);
   return( sqrtl( rval * step * step));
}

int symplectic_6( double jd, ELEMENTS *ref_orbit, double *vect,
                                          const double dt)
{
   int i, j;
#ifdef FOR_REFERENCE_ONLY
         /* Some compilers object to mathematically defined consts,  so  */
         /* I had to replace these lines with explicit numerical consts: */
   const double w1 = -0.117767998417887E1;
   const double w2 = 0.235573213359357E0;
   const double w3 = 0.784513610477560E0;
   const double w0 = (1-2*(w1+w2+w3));
   const double d6[7] = { w3, w2, w1, w0, w1, w2, w3 };
   const double c6[8] = { w3/2, (w3+w2)/2, (w2+w1)/2, (w1+w0)/2,
                         (w1+w0)/2, (w2+w1)/2, (w3+w2)/2, w3/2 };
#endif
   static const double d6[7] = { 0.7845136104775600,  0.2355732133593570,
            -1.1776799841788700, 1.3151863206839060, -1.1776799841788700,
             0.2355732133593570, 0.7845136104775600 };
   static const double c6[8] = { 0.3922568052387800,  0.5100434119184585,
            -0.4710533854097566, 0.0687531682525180,  0.0687531682525180,
            -0.4710533854097566, 0.5100434119184585,  0.3922568052387800 };

   for( i = 0; i < 8; i++)
      {
      double deriv[6];
      const double step = dt * c6[i];

      for( j = 0; j < 3; j++)
         vect[j] += step * vect[j + 3];
      jd += step;
      if( i != 7)
         {
         assert( fabs( jd) < 1e+9);
         for( j = 3; j < 6; j++)
            deriv[j] = 0.;
         calc_derivatives( jd, vect, deriv, ref_orbit->central_obj);
         for( j = 3; j < 6; j++)
            vect[j] += dt * d6[i] * deriv[j];
         }
      }
   return( 0);
}

