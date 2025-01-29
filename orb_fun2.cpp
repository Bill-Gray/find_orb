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

#ifndef _WIN32
   #include <unistd.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "watdefs.h"
#include "stringex.h"
#include "comets.h"
#include "afuncs.h"
#include "mpc_obs.h"
#include "date.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

void push_orbit( const double epoch, const double *orbit);  /* orb_fun2.c */
int pop_orbit( double *epoch, double *orbit);               /* orb_fun2.c */
double generate_mc_variant_from_covariance( double *var_orbit,
                                                     const double *orbit);
double improve_along_lov( double *orbit, const double epoch, const double *lov,
          const unsigned n_params, unsigned n_obs, OBSERVE *obs);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit);
double current_jd( void);                       /* elem_out.cpp */
double evaluate_for_simplex_method( const OBSERVE FAR *obs,
                    const int n_obs, const double *orbit,
                    const int planet_orbiting,
                    const char *limited_orbit);     /* orb_func.cpp */
void init_simplex( double **vects, double *fvals,
         double (*f)( void *context, const double *vect),
               void *context, const int n);        /* simplex.c */
int get_residual_data( const OBSERVE *obs, double *xresid, double *yresid);
int simplex_step( double **vects, double *fvals,
         double (*f)( void *context, const double *vect),
               void *context, const int n);        /* simplex.c */
int apply_excluded_observations_file( OBSERVE *obs, const int n_obs);
int write_excluded_observations_file( const OBSERVE *obs, int n_obs);
char **load_file_into_memory( const char *filename, size_t *n_lines,
                        const bool fail_if_not_found);      /* mpc_obs.cpp */

extern int n_orbit_params, force_model;
extern int available_sigmas;

typedef struct
   {
   OBSERVE FAR *obs;
   int n_obs, n_params;
   const char *constraints;
   double orbit[MAX_N_PARAMS];
   } simplex_context_t;

static double simplex_scoring( void *icontext, const double *ivect)
{
   simplex_context_t *context = (simplex_context_t *)icontext;
   double rval;

   if( context->n_params == 2)
      {
      if( herget_method( context->obs, context->n_obs, ivect[0], ivect[1],
                           context->orbit, NULL, NULL, NULL))
         return( 1e+30);
      if( adjust_herget_results( context->obs, context->n_obs, context->orbit))
         return( 1e+30);
      }
   else
      {
      memcpy( context->orbit, ivect, n_orbit_params * sizeof( double));
      set_locs( context->orbit, context->obs[0].jd, context->obs,
                                                    context->n_obs);
      }
   rval = evaluate_for_simplex_method( context->obs, context->n_obs,
                        context->orbit, 0, context->constraints);
   return( rval);
}

int simplex_method( OBSERVE FAR *obs, int n_obs, double *orbit,
               const double r1, const double r2, const char *constraints)
{
   int i, iter;
   int max_iter = atoi( get_environment_ptr( "SIMPLEX_ITER"));
   double rvals[MAX_N_PARAMS], *rptr[3], scores[3];
   simplex_context_t context;

   for( i = 0; i < 3; i++)
      {
      rptr[i] = rvals + i * 2;
      rptr[i][0] = r1 * (i == 1 ? 1.5 : 1);
      rptr[i][1] = r2 * (i == 2 ? 1.5 : 1);
      }
   context.obs = obs;
   context.n_obs = n_obs;
   context.n_params = 2;
   context.constraints = constraints;
   init_simplex( rptr, scores, simplex_scoring, &context, context.n_params);

   if( !max_iter)
      max_iter = 70;
   for( iter = 0; iter < max_iter; iter++)
      simplex_step( rptr, scores, simplex_scoring, &context, context.n_params);
   memcpy( orbit, context.orbit, n_orbit_params * sizeof( double));
   available_sigmas = NO_SIGMAS_AVAILABLE;
   return( iter);
}

#ifdef NOT_USED_YET
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
#endif

int superplex_method( OBSERVE FAR *obs, int n_obs, double *orbit, const char *constraints)
{
   int i, iter;
   int max_iter = atoi( get_environment_ptr( "SUPERPLEX_ITER"));
   double *rptr[MAX_N_PARAMS + 1], scores[MAX_N_PARAMS + 1];
   double *rvals = (double *)calloc( n_orbit_params * (n_orbit_params + 1),
                                       sizeof( double));
   simplex_context_t context;

   for( i = 0; i <= n_orbit_params; i++)
      {
      rptr[i] = rvals + i * n_orbit_params;
      memcpy( rptr[i], orbit, n_orbit_params * sizeof( double));
      if( i < n_orbit_params)
         rptr[i][i] += (i < 3 ? .2 : .03);
      }
   context.obs = obs;
   context.n_obs = n_obs;
   context.n_params = n_orbit_params;
   context.constraints = constraints;
   init_simplex( rptr, scores, simplex_scoring, &context, context.n_params);

   if( !max_iter)
      max_iter = 70;
   for( iter = 0; iter < max_iter; iter++)
      simplex_step( rptr, scores, simplex_scoring, &context, context.n_params);
   memcpy( orbit, context.orbit, n_orbit_params * sizeof( double));
   available_sigmas = NO_SIGMAS_AVAILABLE;
   free( rvals);
   return( iter);
}

 /* For filtering to work,  you need at least three observations with */
 /* residuals inside the desired max_residual limit.  We do one pass  */
 /* just to make sure there are at least three observations that'll   */
 /* still be active when the filtering is done.  If not,  we make     */
 /* no changes and return -1.  If we succeed,  we return the number   */
 /* of 'flipped' observations.                                        */

int filter_obs( OBSERVE FAR *obs, const int n_obs,
                  const double max_residual, const int filter_type)
{
   int i, pass, n_active = 0, rval = 0;

   for( pass = 0; pass < 2; pass++)
      {
      double dx, dy;

      for( i = 0; i < n_obs; i++)
         if( get_residual_data( obs + i, &dx, &dy)
                            && !(obs[i].flags & OBS_DONT_USE))
            {
            int is_okay;

            if( filter_type)    /* actually filtering by arcseconds */
               {
               dx = (obs[i].ra - obs[i].computed_ra) * cos( obs[i].dec);
               dy = obs[i].dec - obs[i].computed_dec;
               dx *= 3600 * (180. / PI);
               dy *= 3600 * (180. / PI);
               }
            is_okay = (dx * dx + dy * dy < max_residual * max_residual);
            if( !pass && is_okay)
               {
               n_active++;
               if( n_active == 3)      /* found enough;  break out of loop */
                  break;
               }
            if( pass && is_okay != obs[i].is_included)
               {
               obs[i].is_included ^= 1;
               rval++;
               }
            }
      if( n_active < 3)    /* failure */
         return( -1);
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
   double orbit[MAX_N_PARAMS];
   int n_orbit_params, force_model;
   unsigned perturbers;
   } *stored = NULL;

extern unsigned perturbers;

void push_orbit( const double epoch, const double *orbit)
{
   STORED_ORBIT *head = (STORED_ORBIT *)malloc( sizeof( STORED_ORBIT));

   assert( head);
   if( head)
      {
      head->prev = stored;
      head->epoch = epoch;
      memcpy( head->orbit, orbit, n_orbit_params * sizeof( double));
      head->n_orbit_params = n_orbit_params;
      head->force_model = force_model;
      head->perturbers = perturbers;
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
         n_orbit_params = stored->n_orbit_params;
         force_model = stored->force_model;
         perturbers = stored->perturbers;
         memcpy( orbit, stored->orbit, n_orbit_params * sizeof( double));
         available_sigmas = NO_SIGMAS_AVAILABLE;
         }
      free( stored);
      stored = tptr;
      }
   return( rval);
}

void pop_all_orbits( void)
{
   while( !pop_orbit( NULL, NULL))
      ;
}

FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
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
   ofile = fopen_ext( "gauss.out", "tfcwb");
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
most,  four iterations are required,  near y = +/- erf_limit.
(With a caveat : Ralph Pass found that on his Mac,  for y=0.9993,
the iterations settled down with 'diff' alternating between
+3.074775e-14 and -3.074775e-14.  So the code now breaks out after
ten iterations.)               */

static double inverf( const double y)
{
   double x, diff;
   const double erf_limit = .915;
   int n_iterations = 0;

   if( y < -erf_limit)
      return( -inverf( -y));
   else if( y > erf_limit)
      x = sqrt( -log( 1. - y)) - .34;
   else         /* a passable cubic approximation for -.915 < y < .915 */
      x = y * (.8963 + y * y * (.0889 + .494 * y * y));
   do
      {
      const double SQRT_PI =
   1.7724538509055160272981674833411451827975494561223871282138077898529113;
      const double dy = erf( x) - y;
      const double slope = (2. / SQRT_PI) * exp( -x * x);

/*    diff = -dy / slope;      Just doing this would be Newton-Raphson */
      diff = -dy / (slope + x * dy);   /* This gets us Halley's method */
      x += diff;
      }
      while( n_iterations++ < 10 && fabs( diff) > 1e-14);
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

static inline double pseudo_resid( const double y, const double x)
{
   if( y < x && y > -x)
      return( y / x);
   else
      return( 2. - x / fabs(y));
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
      dx = pseudo_resid( xyz[0], xyz[2]);
      dy = pseudo_resid( xyz[1], xyz[2]);
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

double improve_along_lov( double *orbit, const double epoch, const double *lov,
          const unsigned n_params, unsigned n_obs, OBSERVE *obs)
{
   unsigned i, j;
   unsigned n_divs = atoi( get_environment_ptr( "IMPROVE_ALONG_LOV_DIVS"));
   double *x, *score;
   double *xyz, *slopes;
   const double delta = 0.0001;
   double rval, lowest_score;

   while( !obs->is_included && n_obs)
      {
      n_obs--;
      obs++;
      }
   while( !obs[n_obs - 1].is_included && n_obs)
      n_obs--;
   if( !n_obs)
      return( 0.);

   xyz = (double *)calloc( n_obs * 6, sizeof( double));
   slopes = xyz + n_obs * 3;
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
   if( !n_divs)
      n_divs = 10000;
   x = (double *)calloc( n_divs * 2, sizeof( double));
   assert( x);
   score = x + n_divs;
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included && obs[i].note2 != 'R')
         {
         double sigma = hypot( obs[i].posn_sigma_1, obs[i].posn_sigma_2)
                                 * 180. * 3600. / PI;

         rotate_obs_vect( obs + i, xyz + i * 3);
         rotate_obs_vect( obs + i, slopes + i * 3);
         xyz[i * 3] /= sigma;
         xyz[i * 3 + 1] /= sigma;
         slopes[i * 3] /= sigma;
         slopes[i * 3 + 1] /= sigma;
         }
      else
         {
         for( j = 0; j < 3; j++)
            xyz[i * 3 + j] = slopes[i * 3 + j] = 0.;
         xyz[i * 3 + 2] = 1.;
         }

   for( i = 0; i < n_divs; i++)
      {
      x[i] = inverf( 2. * ((double)i + .5) / (double)n_divs - 1.);
      x[i] *= 3.;
      score[i] = search_score( xyz, slopes, n_obs, x[i]);
      }
   i = 0;
   while( score[0] < score[1] && score[0] < score[2] && i++ < 20)
      {
      x[0] += x[0] - x[1];
      score[0] = search_score( xyz, slopes, n_obs, x[0]);
      }
   i = 0;
   while( score[n_divs - 1] < score[n_divs - 2] && score[n_divs - 1] < score[n_divs - 3]
                     && i++ < 20)
      {
      x[n_divs - 1] += x[n_divs - 1] - x[n_divs - 2];
      score[n_divs - 1] = search_score( xyz, slopes, n_obs, x[n_divs - 1]);
      }
   rval = 0.;
   lowest_score = 1e+200;
   for( i = 0; i < n_divs - 2; i++)
      if( score[i + 1] < score[i] && score[i + 1] < score[i + 2])
         {
         double new_x;
         const double new_score = find_score_minimum( xyz, slopes, n_obs,
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
   free( x);
   return( rval);
}

double gaussian_random( void);                           /* monte0.c */

double generate_mc_variant_from_covariance( double *var_orbit,
                                                     const double *orbit)
{
   int i, j;
   extern double **eigenvects;
   double rval = 0.;

   assert( eigenvects);
   memcpy( var_orbit, orbit, 6 * sizeof( double));
   for( i = 0; i < n_orbit_params; i++)
      {
      const double g_rand = gaussian_random( );

      for( j = 0; j < n_orbit_params; j++)
         var_orbit[j] += g_rand * eigenvects[i][j];
      rval += g_rand * g_rand;
      }
   return( rval);
}

int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */

const char *excluded_filename = "excluded.txt";

int write_excluded_observations_file( const OBSERVE *obs, int n_obs)
{
   FILE *ofile;
   int n_excluded = 0;
   size_t i, n_lines;
   char **ilines = load_file_into_memory( excluded_filename, &n_lines, false);
   char reduced_desig[13];
   bool write_it_out = true;
   size_t len;

   ofile = fopen_ext( excluded_filename, "fcwb");
   strlcpy_error( reduced_desig, obs->packed_id);
   text_search_and_replace( reduced_desig, " ", "");
   len = strlen( reduced_desig);
   if( ilines)
      {
      for( i = 0; i < n_lines; i++)
         {
         if( strstr( ilines[i], "Banned obs"))
            write_it_out = memcmp( ilines[i] + 2, reduced_desig, len)
                           || ilines[i][len + 2] != ' ';
         if( write_it_out)
            fprintf( ofile, "%s\n", ilines[i]);
         }
      free( ilines);
      }

   while( n_obs--)
      {
      if( !obs->is_included)
         {
         if( !n_excluded)
            {
            char time_buff[80];

            full_ctime( time_buff, current_jd( ), FULL_CTIME_YMD);
            fprintf( ofile, "# %s Banned obs file written %s UTC\n",
                           reduced_desig, time_buff);
            }
         n_excluded++;
         fprintf( ofile, "%s %.6f %010.6f %+010.6f\n",
                     obs->mpc_code, obs->jd - 2400000.5,
                     obs->ra * 180. / PI, obs->dec * 180. / PI);
         }
      obs++;
      }
   fclose( ofile);
   return( n_excluded);
}

double automatic_outlier_rejection_limit;
double default_automatic_outlier_rejection_limit = 3.;

int apply_excluded_observations_file( OBSERVE *obs, const int n_obs)
{
   char buff[90];
   FILE *ifile;
   int n_excluded = 0, i;

   automatic_outlier_rejection_limit = default_automatic_outlier_rejection_limit;
   ifile = fopen_ext( excluded_filename, "crb");
   if( ifile)
      {
      bool read_it_in = false;
      char reduced_desig[13];
      size_t len;

      strlcpy_error( reduced_desig, obs->packed_id);
      text_search_and_replace( reduced_desig, " ", "");
      len = strlen( reduced_desig);
      while( fgets( buff, sizeof( buff), ifile))
         {
         if( strstr( buff, "Banned obs"))
            read_it_in = !memcmp( buff + 2, reduced_desig, len)
                           && buff[len + 2] == ' ';
         if( read_it_in)
            {
            if( *buff != '#')
               {
               const double jd = 2400000.5 + atof( buff + 4);
               const double ra = atof( buff + 17) * PI / 180.;
               const double dec = atof( buff + 28) * PI / 180.;
               const double date_tolerance = 1e-6;
               const double ang_tolerance = 2e-6 * (PI / 180.);

               for( i = 0; i < n_obs; i++)
                  if( fabs( obs[i].jd - jd) < date_tolerance
                        && fabs( obs[i].ra - ra) < ang_tolerance
                        && fabs( obs[i].dec - dec) < ang_tolerance)
                     {
                     obs[i].flags |= OBS_DONT_USE;
                     obs[i].is_included = 0;
                     n_excluded++;
                     }
               }
            else if( !memcmp( buff, "# Outlier rejection ", 20))
               automatic_outlier_rejection_limit = atof( buff + 20);
            }
         }
      fclose( ifile);
      }
   return( n_excluded);
}
