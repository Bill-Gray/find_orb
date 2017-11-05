/* orb_fun2.cpp: basic orbital element/numerical integration funcs

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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "watdefs.h"
#include "comets.h"
#include "afuncs.h"
#include "mpc_obs.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

int generate_mc_variant_from_covariance( double *orbit);    /* orb_fun2.cpp */
double improve_along_lov( double *orbit, const double epoch, const double *lov,
          const unsigned n_params, const unsigned n_obs, OBSERVE *obs);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit);
double evaluate_for_simplex_method( const OBSERVE FAR *obs,
                    const int n_obs, const double *orbit,
                    const int planet_orbiting,
                    const char *limited_orbit);     /* orb_func.cpp */

#define SIMPLEX struct simplex

SIMPLEX
   {
   double r1, r2, score;
   };

extern int available_sigmas;
extern int debug_level;
int debug_printf( const char *format, ...);


static double try_simplex_reflection( OBSERVE FAR *obs, int n_obs,
               SIMPLEX simplex[3], const double reflect, const char *constraints)
{
   double r1 = (simplex[0].r1 + simplex[1].r1) * (1. - reflect) / 2.
                     + reflect * simplex[2].r1;
   double r2 = (simplex[0].r2 + simplex[1].r2) * (1. - reflect) / 2.
                     + reflect * simplex[2].r2;
   double orbit[6], new_score;
   int i;

                   /* guard against negative/artificially low values */
   for( i = 0; i < 3; i++)
      {
      if( r1 < simplex[i].r1 / 2.)
         r1 = simplex[i].r1 / 2.;
      if( r2 < simplex[i].r2 / 2.)
         r2 = simplex[i].r2 / 2.;
      }

   herget_method( obs, n_obs, r1, r2, orbit, NULL, NULL, NULL);
   adjust_herget_results( obs, n_obs, orbit);
   new_score = evaluate_for_simplex_method( obs, n_obs, orbit, 0, constraints);
   if( new_score <= simplex[2].score)         /* step is an improvement */
      {
      simplex[2].r1 = r1;
      simplex[2].r2 = r2;
      simplex[2].score = new_score;
      }
   return( new_score);
}

int simplex_method( OBSERVE FAR *obs, int n_obs, double *orbit,
               const double r1, const double r2, const char *constraints)
{
   SIMPLEX simplex[3];
   int i, j, iter;
   int max_iter = atoi( get_environment_ptr( "SIMPLEX_ITER"));
   double temp_orbit[6];

   if( !max_iter)
      max_iter = 70;
   for( i = 0; i < 3; i++)
      {
      simplex[i].r1 = r1 * ((i & 1) ? 1. : 1.1);
      simplex[i].r2 = r2 * ((i & 2) ? 1. : 1.1);
      herget_method( obs, n_obs, simplex[i].r1, simplex[i].r2,
                                temp_orbit, NULL, NULL, NULL);
      adjust_herget_results( obs, n_obs, orbit);
      simplex[i].score = evaluate_for_simplex_method( obs, n_obs, temp_orbit, 0, constraints);
      }
   for( iter = 0; iter < max_iter; iter++)
      {
      double new_score;

      for( i = 1; i < 3; i++)
         for( j = 0; j < i; j++)    /* sort so simplex[0] = lowest-score, */
            if( simplex[i].score < simplex[j].score)  /* simplex[2] = highest */
               {
               SIMPLEX temp_simp = simplex[i];

               simplex[i] = simplex[j];
               simplex[j] = temp_simp;
               }
      if( debug_level > 2)
         {
         double dot_prod = (simplex[1].r1 - simplex[0].r1) *
                           (simplex[2].r2 - simplex[0].r2) -
                           (simplex[2].r1 - simplex[0].r1) *
                           (simplex[1].r2 - simplex[0].r2);

         debug_printf( "Simplex %d: %lg\n", iter, dot_prod);
         for( i = 0; i < 3; i++)
            debug_printf( "   r1 = %f, r2 = %f: score %f (%f %f)\n",
                        simplex[i].r1, simplex[i].r2, simplex[i].score,
                        simplex[i].r1 - simplex[0].r1,
                        simplex[i].r2 - simplex[0].r2);
         }
      new_score = try_simplex_reflection( obs, n_obs, simplex, -1., constraints);
                  /* If step was a new 'best',  try doubling it: */
      if( new_score < simplex[0].score)
         {
         try_simplex_reflection( obs, n_obs, simplex, 2., constraints);
         if( debug_level > 2)
            debug_printf( "Doubled\n");
         }
      else if( new_score > simplex[1].score)
         {
         const double contracted =
               try_simplex_reflection( obs, n_obs, simplex, .5, constraints);

         if( debug_level > 2)
            debug_printf( contracted > new_score ?
                              "Contracting\n" : "Half-con\n");
         if( contracted > new_score)   /* can't get rid of it;  try  */
            for( i = 1; i < 3; i++)     /* contracting around our best point */
               {
               simplex[i].r1 += .5 * (simplex[0].r1 - simplex[i].r1);
               simplex[i].r2 += .5 * (simplex[0].r2 - simplex[i].r2);
               herget_method( obs, n_obs, simplex[i].r1, simplex[i].r2,
                                temp_orbit, NULL, NULL, NULL);
               adjust_herget_results( obs, n_obs, orbit);
               simplex[i].score = evaluate_for_simplex_method( obs, n_obs,
                                                   temp_orbit, 0, constraints);
               }
         }
      else
         if( debug_level > 2)
            debug_printf( "Simple step\n");
      }
   herget_method( obs, n_obs, simplex[0].r1, simplex[0].r2,
                                orbit, NULL, NULL, NULL);
   adjust_herget_results( obs, n_obs, orbit);
   return( iter);
}

static void adjust_orbit_to_constraints( double *orbit, const char *constraints)
{
   if( constraints && !strcmp( constraints, "e=1"))
      {
      const double r = vector3_length( orbit);
      const double v_esc = sqrt( 2. * SOLAR_GM / r);
      const double v = vector3_length( orbit + 3);
      const double rescale = v_esc / v;
      int i;

      for( i = 3; i < 6; i++)
         orbit[i] *= rescale;
      }
}

#define SUPERPLEX struct superplex

SUPERPLEX
   {
   double orbit[6];
   double score;
   };

int set_locs( const double *orbit, double t0, OBSERVE FAR *obs, int n_obs);

static double try_superplex_reflection( OBSERVE FAR *obs, int n_obs,
               SUPERPLEX superplex[7], const double reflect,
               const char *constraints)
{
   double orbit[6], new_score;
   int i, j;
   const int n_reflect = 6;

   for( i = 0; i < 6; i++)
      {
      orbit[i] = superplex[6].orbit[i] * reflect;
      for( j = 0; j < n_reflect; j++)
         orbit[i] += (1. - reflect) * superplex[j].orbit[i] / (double)n_reflect;
      }

   adjust_orbit_to_constraints( orbit, constraints);
   set_locs( orbit, obs[0].jd, obs, n_obs);
   new_score = evaluate_for_simplex_method( obs, n_obs, orbit, 0, constraints);
   if( new_score <= superplex[6].score)         /* step is an improvement */
      {
      memcpy( superplex[6].orbit, orbit, 6 * sizeof( double));
      superplex[6].score = new_score;
      }
   return( new_score);
}

extern int debug_level;
int debug_printf( const char *format, ...);

int superplex_method( OBSERVE FAR *obs, int n_obs, double *orbit, const char *constraints)
{
   SUPERPLEX superplex[7];
   int iter, i;
   int max_iter = atoi( get_environment_ptr( "SUPERPLEX_ITER"));

   if( !max_iter)
      max_iter = 70;
   while( n_obs && !obs[n_obs - 1].is_included)
      n_obs--;
   for( i = 0; i < 7; i++)
      {
      memcpy( superplex[i].orbit, orbit, 6 * sizeof( double));
      if( i < 6)
         superplex[i].orbit[i] *= .9999;
      adjust_orbit_to_constraints( superplex[i].orbit, constraints);
      set_locs( superplex[i].orbit, obs[0].jd, obs, n_obs);
      superplex[i].score = evaluate_for_simplex_method( obs, n_obs,
                                  superplex[i].orbit, 0, constraints);
      }
   for( iter = 0; iter < max_iter; iter++)
      {
      int j;
      double new_score;

      for( i = 1; i < 7; i++)       /* sort in increasing order of score */
         {
         for( j = i; j && superplex[j - 1].score > superplex[i].score;
                              j--)
            ;
         if( j != i)
            {
            SUPERPLEX temp_simp = superplex[i];

            memmove( superplex + j + 1, superplex + j,
                        (i - j) * sizeof( SUPERPLEX));
            superplex[j] = temp_simp;
            }
         }
      new_score = try_superplex_reflection( obs, n_obs, superplex, -1., constraints);
      if( debug_level > 2)
         debug_printf( "Iter %d: score %f\n", iter, superplex[6].score);
      if( debug_level > 3)
         for( i = 0; i < 7; i++)
            {
            debug_printf( "%d (%f): ", i, superplex[i].score);
            for( j = 0; j < 6; j++)
               debug_printf( "%11.6f", superplex[i].orbit[j]);
            debug_printf( "\n");
            }
                  /* If step was a new 'best',  try doubling it: */
      if( new_score < superplex[0].score)
         {
         try_superplex_reflection( obs, n_obs, superplex, 2., constraints);
         if( debug_level > 2)
            debug_printf( "Doubled\n");
         }
      else if( new_score >= superplex[5].score)
         {
         const double contracted =
               try_superplex_reflection( obs, n_obs, superplex, .5, constraints);

         if( debug_level > 2)
            debug_printf( contracted > new_score ?
                              "Contracting\n" : "Half-con\n");
         if( contracted > new_score)   /* can't get rid of it;  try  */
            for( i = 1; i < 6; i++)     /* contracting around our best point */
               {
               for( j = 0; j < 6; j++)
                  superplex[i].orbit[j] =
                        (superplex[i].orbit[j] + superplex[0].orbit[j]) * .5;
               adjust_orbit_to_constraints( superplex[i].orbit, constraints);
               set_locs( superplex[i].orbit, obs[0].jd, obs, n_obs);
               superplex[i].score = evaluate_for_simplex_method( obs, n_obs,
                                  superplex[i].orbit, 0, constraints);
               }
         }
      else
         if( debug_level > 2)
            debug_printf( "Simple step\n");
      }
   memcpy( orbit, superplex[0].orbit, 6 * sizeof( double));
   set_locs( superplex[0].orbit, obs[0].jd, obs, n_obs);
   available_sigmas = NO_SIGMAS_AVAILABLE;
   return( iter);
}

 /* For filtering to work,  you need at least three observations with */
 /* residuals inside the desired max_residual limit.  We do one pass  */
 /* just to make sure there are at least three observations that'll   */
 /* still be active when the filtering is done.  If not,  we make     */
 /* no changes and return FILTERING_FAILED.                           */

#define FILTERING_CHANGES_MADE            1
#define FILTERING_NO_CHANGES_MADE         2
#define FILTERING_FAILED                  3

int filter_obs( OBSERVE FAR *obs, const int n_obs,
                  const double max_residual_in_arcseconds)
{
   const double max_resid =            /* cvt arcseconds to radians */
                   max_residual_in_arcseconds * PI / (180. * 3600.);
   int i, pass, n_active = 0, rval = FILTERING_NO_CHANGES_MADE;

   for( pass = 0; pass < 2; pass++)
      {
      for( i = 0; i < n_obs && rval != FILTERING_FAILED; i++)
         {
         const double dy = obs[i].dec - obs[i].computed_dec;
         const double dx = (obs[i].ra - obs[i].computed_ra) * cos( obs[i].dec);
         const int is_okay = (dx * dx + dy * dy < max_resid * max_resid);

         if( !pass && is_okay)
            {
            n_active++;
            if( n_active == 3)      /* found enough;  break out of loop */
               break;
            }
         if( pass && is_okay != obs[i].is_included)
            {
            obs[i].is_included ^= 1;
            rval = FILTERING_CHANGES_MADE;
            }
         }
      if( n_active < 3)
         rval = FILTERING_FAILED;
      }
   return( rval);
}

/* For the 'undo' functions in Find_Orb,  it's convenient to have
a way to store orbits on a stack,  then retrieve them in the usual
last-in-first-out order.  Also,  if you switch to a new object,  you
presumably don't want all the stored orbits retained;  hence the
'pop_all_orbits()' function.  In reality,  the orbits are stored
in a linked list. */

#define STORED_ORBIT struct stored_orbit

static STORED_ORBIT
   {
   STORED_ORBIT *prev;
   double epoch;
   double orbit[6];
   double solar_pressure[MAX_N_NONGRAV_PARAMS];
   int n_extra_params;
   } *stored = NULL;

extern double solar_pressure[];
extern int n_extra_params;

void push_orbit( const double epoch, const double *orbit)
{
   STORED_ORBIT *head = (STORED_ORBIT *)malloc( sizeof( STORED_ORBIT));

   assert( head);
   if( head)
      {
      head->prev = stored;
      head->epoch = epoch;
      memcpy( head->orbit, orbit, 6 * sizeof( double));
      memcpy( head->solar_pressure, solar_pressure, n_extra_params * sizeof( double));
      head->n_extra_params = n_extra_params;
      stored = head;
      }
}

int pop_orbit( double *epoch, double *orbit)
{
   const int rval = (stored ? 0 : -1);

   if( stored)
      {
      STORED_ORBIT *tptr = stored->prev;

      if( epoch)
         {
         *epoch = stored->epoch;
         memcpy( orbit, stored->orbit, 6 * sizeof( double));
         n_extra_params = stored->n_extra_params;
         memcpy( solar_pressure, stored->solar_pressure, n_extra_params * sizeof( double));
         available_sigmas = NO_SIGMAS_AVAILABLE;
         }
      free( stored);
      stored = tptr;
      }
   return( rval);
}

void pop_all_orbits( void)
{
   double unused_epoch;
   double unused_orbit[6];

   while( !pop_orbit( &unused_epoch, unused_orbit))
      ;
}

FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
double dot_product( const double *v1, const double *v2);    /* sr.c */
void set_distance( OBSERVE FAR *obs, double r);             /* orb_func.c */

/* The linear regression fit here is used to determine a perihelion distance
q,  eccentricity ecc,  and longitude of perihelion omega.  The idea is that

r = z / (1. + ecc * cos( lon - omega))

   where z = q (1 + ecc);

z / r = 1 + ecc * cos( lon - omega)
      = 1 + (ecc * cos(omega)) * cos( lon) + (ecc * sin(omega)) * sin( lon)

1/r = a + b cos(lon) + c sin(lon)

   ....where a = 1/z,  b = ecc * cos(omega)/z,  c = ecc * sin(omega)/z.  So
we can solve for best-fit versions of a, b, and c;  from these,  we can
then determine z=1/a, ecc = sqrt( b^2+c^2) * z, omega = atan2(c, b),
then q = z/(1+ecc).   */

int link_arcs( OBSERVE *obs, int n_obs, const double r1, const double r2)
{
   OBSERVE *end_obs = obs + n_obs - 1;
   double a, b, c;
   double rvect[3], rlen, *theta, *r;
   double avect[3], bvect[3], ecc, omega, q;
            /* Coeffs for three-variable multiple linear regression: */
            /* See Meeus,  _Astronomical Algorithms_, p 44 */
   double fit_m = 0., fit_p = 0., fit_q = 0., fit_r = 0., fit_s = 0.;
   double fit_t = 0., fit_u = 0., fit_v = 0., fit_w = 0., fit_d;
   double sum_of_squares = 0.;
   FILE *ofile;
   int i;

   while( n_obs && !obs->is_included)
      {
      obs++;
      n_obs--;
      }
   while( n_obs && !end_obs->is_included)
      {
      end_obs--;
      n_obs--;
      }
   if( n_obs < 3)
      return( -2);
   set_distance( obs, r1);
   set_distance( end_obs, r2);
   vector_cross_product( rvect, obs->obj_posn, end_obs->obj_posn);
   rlen = vector3_length( rvect);
   if( rvect[2] < 0.)      /* assume a prograde orbit */
      rlen = -rlen;
   for( i = 0; i < 3; i++)
      rvect[i] /= rlen;
   ofile = fopen_ext( "gauss.out", "fcwb");
   fprintf( ofile, "%.10f %.10f %.10f\n",
            rvect[0], rvect[1], rvect[2]);
   rlen = sqrt( rvect[0] * rvect[0] + rvect[1] * rvect[1]);
   fprintf( ofile, "Omega = %.5f  incl = %.5f\n",
         180. - atan2( rvect[0], rvect[1]) * 180. / PI,
         atan2( rlen, rvect[2]) * 180. / PI);

   avect[0] = -rvect[1] / rlen;
   avect[1] =  rvect[0] / rlen;
   avect[2] = 0.;
   vector_cross_product( bvect, avect, rvect);

   theta = (double *)calloc( n_obs * 2, sizeof( double));
   if( !theta)
      return( -1);
   r = theta + n_obs;
   for( i = 0; i < n_obs; i++)
      {
      int j;
      const double dval1 = dot_product( obs[i].obs_posn, rvect);
      const double dval2 = dot_product( obs[i].vect, rvect);
      const double dist_to_target = -dval1 / dval2;
      double x, y, loc[3];
      double sin_theta, cos_theta;

                  /* Find where the observed vector would intersect the */
                  /* plane defined by dot_product( loc, rvect) = 0:     */
      for( j = 0; j < 3; j++)
         loc[j] = obs[i].obs_posn[j] + obs[i].vect[j] * dist_to_target;
      fprintf( ofile, "Dist to target: %f (%f)\n", dist_to_target, obs[i].r);
      x = dot_product( loc, avect);
      y = dot_product( loc, bvect);
      theta[i] = atan2( y, x);
      sin_theta = sin( theta[i]);
      cos_theta = cos( theta[i]);
//    r[i] = obs[i].solar_r;
      r[i] = vector3_length( loc);
      fit_p += sin_theta;
      fit_q += cos_theta;
      fit_r += sin_theta * sin_theta;
      fit_s += cos_theta * sin_theta;
      fit_t += cos_theta * cos_theta;
      fit_u += 1. / r[i];
      fit_v += sin_theta / r[i];
      fit_w += cos_theta / r[i];
      }
   fit_m = (double)n_obs;
   fit_d = fit_m * fit_r * fit_t + 2. * fit_p * fit_q * fit_s
            - fit_m * fit_s * fit_s - fit_r * fit_q * fit_q
            - fit_t * fit_p * fit_p;
   a = fit_u * (fit_r * fit_t - fit_s * fit_s)
     + fit_v * (fit_q * fit_s - fit_p * fit_t)
     + fit_w * (fit_p * fit_s - fit_q * fit_r);
   b = fit_u * (fit_s * fit_q - fit_p * fit_t)
     + fit_v * (fit_m * fit_t - fit_q * fit_q)
     + fit_w * (fit_p * fit_q - fit_m * fit_s);
   c = fit_u * (fit_p * fit_s - fit_r * fit_q)
     + fit_v * (fit_p * fit_q - fit_m * fit_s)
     + fit_w * (fit_m * fit_r - fit_p * fit_p);
   a /= fit_d;
   b /= fit_d;
   c /= fit_d;
   ecc = sqrt( b * b + c * c) / a;
   omega = atan2( c, b) + 3. * PI / 2.;
   if( omega > PI + PI)
      omega -= PI + PI;
   q = 1. / (a * (1. + ecc));
   fprintf( ofile, "ecc = %f  omega = %f   q = %f   a = %f\n",
         ecc, omega * 180. / PI, q, q / (1. - ecc));
   for( i = 0; i < n_obs; i++)
      {
      const double fitted_r = 1. / (a + b * sin( theta[i]) + c * cos( theta[i]));
      const double dr = r[i] - fitted_r;

      sum_of_squares += dr * dr;
      }
   free( theta);
   fprintf( ofile, "Error %f\n", sqrt( sum_of_squares / (double)n_obs));
   fclose( ofile);
   available_sigmas = NO_SIGMAS_AVAILABLE;
   return( 0);
}

double find_r_given_solar_r( const OBSERVE FAR *obs, const double solar_r);

static int set_up_circular_orbit( OBSERVE FAR *obs1, OBSERVE FAR *obs2,
                  const double solar_r, double *dt, double *t0,
                  double *orbit)
{
   double delta[3], angle;
   int i;

   set_distance( obs1, find_r_given_solar_r( obs1, solar_r));
   set_distance( obs2, find_r_given_solar_r( obs2, solar_r));
   *t0 = solar_r * sqrt( solar_r) / GAUSS_K;
   for( i = 0; i < 3; i++)
      delta[i] = obs2->obj_posn[i] - obs1->obj_posn[i];
   angle = 2. * asin( vector3_length( delta) * .5 / solar_r);
   *dt = angle * (*t0);
   *t0 *= 2. * PI;
// debug_printf( "Solar_r %f, dt = %f, t0 = %f\n",
//          solar_r, *dt, *t0);
   if( orbit)
      {
      double xprod[30], scale, vel[30];

      vector_cross_product( xprod, obs1->obj_posn, obs2->obj_posn);
      vector_cross_product( vel, xprod, obs1->obj_posn);
               /* Circular orbit speed is GAUSS_K / sqrt( solar_r),  so: */
      scale = (GAUSS_K / sqrt( solar_r)) / vector3_length( vel);
      for( i = 0; i < 3; i++)
         {
         orbit[i] = obs1->obj_posn[i];
         orbit[i + 3] = vel[i] * scale;
         }
      }
   return( obs1->r > 0. && obs2->r > 0.);
}

int find_circular_orbits( OBSERVE FAR *obs1, OBSERVE FAR *obs2,
               double *orbit, const int desired_soln)
{
   double dt = obs2->jd - obs1->jd, t0;
   double dt_1 = 0., t0_1 = 0., dt_2, t0_2;
   double r1 = 1.02, r2, radii[50];
   int n_solutions = 0, soln_type, types[50];

   debug_printf( "Delta_t = %f\n", dt);
   while( r1 < 100. && n_solutions <= desired_soln)
      {
      r2 = r1;
      dt_2 = dt_1;
      t0_2 = t0_1;
      if( r1 < 3.5)
         r1 += .05;
      else
         r1 *= 1.05;
      set_up_circular_orbit( obs1, obs2, r1, &dt_1, &t0_1, NULL);
      if( t0_2)
         {
         int bug_out = 0;

         for( soln_type = 0; !bug_out; soln_type++)
            {
            const double delta_1 = t0_1 * (double)( soln_type / 2)
                     + ((soln_type & 1) ? t0_1 - dt_1 : dt_1) - dt;
            const double delta_2 = t0_2 * (double)( soln_type / 2)
                     + ((soln_type & 1) ? t0_2 - dt_2 : dt_2) - dt;

            if( delta_1 * delta_2 < 0.)      /* we have a zero crossing */
               {
               const double r3 = (r1 * delta_2 - r2 * delta_1) / (delta_2 - delta_1);
               double delta_3, t0_3, dt_3;

               set_up_circular_orbit( obs1, obs2, r3, &dt_3, &t0_3, NULL);
               delta_3 = t0_3 * (double)( soln_type / 2)
                     + ((soln_type & 1) ? t0_3 - dt_3 : dt_3) - dt;

               types[n_solutions] = soln_type;
               radii[n_solutions++] =
                     (r1 * delta_3 - r3 * delta_1) / (delta_3 - delta_1);
//             debug_printf( "Soln %d: type %d, r=%f\n",
//                n_solutions, soln_type, radii[n_solutions - 1]);
               }
            if( delta_1 > dt + t0_1)
               bug_out = 1;
            }
         }
      }

   if( !n_solutions)
      return( -1);
   set_up_circular_orbit( obs1, obs2, radii[desired_soln % n_solutions],
                              &dt, &t0, orbit);
   if( types[desired_soln % n_solutions] & 1)
      {
      int i;

      for( i = 3; i < 6; i++)
         orbit[i] *= -1.;
      }
   available_sigmas = NO_SIGMAS_AVAILABLE;
   return( 0);
}

/* Shamelessly copied (with minor changes) from

https://siliconandlithium.blogspot.com/2014/05/msvc-c99-mathh-header.html
https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions

   As described at the second link,  this has maximum error of 1.5x10^-7.
(Which isn't a problem here,  but some caution would be appropriate.)
It's only used in early MSVCs which lack erf(),  and in OpenWATCOM.  */

#if defined( _MSC_VER) && (_MSC_VER < 1800) || defined( __WATCOMC__)
double erf( double x)
{
    const double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    const double a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
    const int sign = (x >= 0) ? 1 : -1;
    double t, y;

    x = fabs(x);
    t = 1.0 / (1.0 + p*x);
    y = 1.0 - (((((a5 * t + a4 ) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);
    return sign*y;
}
#endif

/* In computing the inverse of the error function,  it helps that the
derivative of erf( ) is easily computed.  Easy derivatives mean easy
Newton-Raphson root-finding.  An easy second derivative,  in this case,
means an easy Halley's root-finder,  with cubic convergence.

d erf(x)
-------- = (2/sqrt(pi)) exp( -x^2)
  dx

d2 erf(x)      d erf(x)
-------- = -2x --------
  d2x            dx

   Some cancelling out of terms happens with the second derivative
that make it particularly well-suited to Halley's method.

   We start out with a rough approximation for inverf( y),  and
then do a root search for y = erf( x) with Halley's method.  At
most,  four iterations are required,  near y = +/- erf_limit.   */

static double inverf( const double y)
{
   double x, diff;
   const double erf_limit = .95;

   if( y < -erf_limit)
      return( -inverf( -y));
   else if( y > erf_limit)
      x = sqrt( -log( 1. - y)) - .34;
   else         /* a passable cubic approximation for -.9 < y < .9 */
      x = y * (.8963 + y * y * (.0889 + .494 * y * y));
   do
      {
      const double SQRT_PI =
   1.7724538509055160272981674833411451827975494561223871282138077898529113;
      const double dy = erf( x) - y;
      const double slope = (2. / SQRT_PI) * exp( -x * x);

      diff = -dy / slope;   /* Just doing this would be Newton-Raphson */
      diff = -dy / (slope + x * dy);   /* This gets us Halley's method */
      x += diff;
      }
      while( fabs( diff) > 1e-14);
   return( x);
}

/* Convert the incoming J2000 vector into the system where the z-axis
is along the observed ray,  the x-axis is perpendicular to it and in
the plane of the ecliptic,  and the y-axis is perpendicular to both.
In this system,  x/z and y/z are decent approximations to the
residuals... which we'll try to minimize in the subsequent code. */

static void rotate_obs_vect( const OBSERVE *obs, double *vect)
{
   double sideways[3], xprod[3], x, y, z;
   const double r = sqrt( obs->vect[0] * obs->vect[0]
                        + obs->vect[1] * obs->vect[1]);

   sideways[0] =  obs->vect[1] / r;
   sideways[1] = -obs->vect[0] / r;
   sideways[2] = 0.;
   vector_cross_product( xprod, sideways, obs->vect);
            /* obs->vect, sideways,  xprod now form an orthonormal matrix */
   x = dot_product( sideways, vect);
   y = dot_product( xprod, vect);
   z = dot_product( obs->vect, vect);
   *vect++ = x;
   *vect++ = y;
   *vect++ = z;
}

static double search_score( const double *loc, const double *deriv, size_t n_obs,
                             const double sigmas)
{
   double rval = 0.;

   while( n_obs--)
      {
      double xyz[3], dx, dy;
      size_t i;

      for( i = 0; i < 3; i++)
         xyz[i] = loc[i] + sigmas * deriv[i];
      dx = xyz[0] / xyz[2];
      dy = xyz[1] / xyz[2];
      loc += 3;
      deriv += 3;
      rval += dx * dx + dy * dy;
      }
   return( rval);
}

double find_parabolic_minimum_point( const double x[3], const double y[3]);

/* Theoretically speaking,  this one-dimensional minimization should use
something like Brent's method.  But the minimum really is nearly parabolic
here;  convergence is quite fast and (nearly) guaranteed. */

static double find_score_minimum( const double *xyz, const double *slopes,
               const int n_obs,
               const double start_x[3], const double start_y[3], double *xmin)
{
   double x[3], y[3];
   int iter = 5;

   memcpy( x, start_x, 3 * sizeof( double));
   memcpy( y, start_y, 3 * sizeof( double));
   while( iter--)
      {
      double new_x = find_parabolic_minimum_point( x, y);

      x[0] = x[1];
      y[0] = y[1];
      x[1] = x[2];
      y[1] = y[2];
      x[2] = new_x;
      y[2] = search_score( xyz, slopes, n_obs, new_x);
      }
   *xmin = x[2];
   return( y[2]);
}

/* The logic for searching for a "best fit" along the line of improvement
is as follows.  We sample at N_DIVS points,  using the above inverse of
the cumulative error function so that most of the searching is near
the center,  with sparser searching toward the tails of the distribution.
If,  between those points,  we find a minimum "score" (sum of the squares
of the residuals,  unweighted),  we do a parabolic search for the real
minimum.  We may find multiple minima.  We take the lowest of them.  */

#define N_DIVS 100

double improve_along_lov( double *orbit, const double epoch, const double *lov,
          const unsigned n_params, const unsigned n_obs, OBSERVE *obs)
{
   unsigned i, j;
   double x[N_DIVS], score[N_DIVS];
   double *xyz = (double *)calloc( n_obs * 6, sizeof( double));
   double *slopes = xyz + n_obs * 3;
   const double delta = 0.0001;
   double rval, lowest_score;

   assert( xyz);
   for( i = 0; i < n_obs; i++)
      for( j = 0; j < 3; j++)
         xyz[i * 3 + j] = obs[i].obj_posn[j] - obs[i].obs_posn[j];
   for( i = 0; i < n_params; i++)
      orbit[i] += delta * lov[i];
   set_locs( orbit, epoch, obs, n_obs);
   for( i = 0; i < n_obs; i++)
      for( j = 0; j < 3; j++)
         slopes[i * 3 + j] = obs[i].obj_posn[j] - obs[i].obs_posn[j];
   for( i = 0; i < n_obs * 3; i++)
      slopes[i] = (slopes[i] - xyz[i]) / delta;
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         rotate_obs_vect( obs + i, xyz + i * 3);
         rotate_obs_vect( obs + i, slopes + i * 3);
         }
      else
         for( j = 0; j < 3; j++)
            slopes[i * 3 + j] = 0.;

   for( i = 0; i < N_DIVS; i++)
      {
      x[i] = inverf( 2. * ((double)i + .5) / (double)N_DIVS - 1.);
      x[i] *= 3.;
      score[i] = search_score( xyz, slopes, n_obs, x[i]);
      }
   rval = 0.;
   lowest_score = 1e+200;
   for( i = 0; i < N_DIVS - 2; i++)
      if( score[i + 1] < score[i] && score[i + 1] < score[i + 2])
         {
         double new_x;
         double new_score = find_score_minimum( xyz, slopes, n_obs,
                                       x + i, score + i, &new_x);

         if( lowest_score > new_score)
            {
            lowest_score = new_score;
            rval = new_x;
            }
         }
   for( i = 0; i < n_params; i++)
      orbit[i] += (rval - delta) * lov[i];
   set_locs( orbit, epoch, obs, n_obs);
   free( xyz);
   return( rval);
}

double gaussian_random( void);                           /* monte0.c */

int generate_mc_variant_from_covariance( double *orbit)
{
   int i, j;
   extern double **eigenvects;

   assert( eigenvects);
   for( i = 0; i < 6 + n_extra_params; i++)
      {
      const double g_rand = gaussian_random( );

      for( j = 0; j < 6; j++)
         orbit[j] += g_rand * eigenvects[i][j];
      for( j = 0; j < n_extra_params; j++)
         solar_pressure[j] += g_rand * eigenvects[i][j + 6];
      }
   return( 0);
}
