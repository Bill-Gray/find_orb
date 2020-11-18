/* orb_func.cpp: basic orbital element/numerical integration funcs

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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include "watdefs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "lsquare.h"
#include "date.h"
#include "afuncs.h"
#include "lunar.h"
#include "monte0.h"

#ifndef _MSC_VER
         /* All non-Microsoft builds are for the console */
   #define CONSOLE
#endif

#ifndef _WIN32
   #include <unistd.h>
#endif

#ifdef CONSOLE
      /* In the console version of Find_Orb,  the following two functions */
      /* get remapped to Curses functions.  In the non-interactive one,   */
      /* they're mapped to 'do-nothings'.  See fo.cpp & find_orb.cpp.     */
   void refresh_console( void);
   void move_add_nstr( const int col, const int row, const char *msg, const int n_bytes);
#endif

/* MS only got around to adding 'isfinite' in VS2013 : */

#if defined( _MSC_VER) && (_MSC_VER < 1800)
   #include <float.h>
   #define isfinite _finite
#endif

#if defined( _MSC_VER) && (_MSC_VER < 1900)
                      /* For older MSVCs,  we have to supply our own  */
                      /* snprintf().  See snprintf.cpp for details.  */
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

unsigned perturbers = 0;
int integration_method = 0;
extern int debug_level;

int generic_message_box( const char *message, const char *box_type);
int find_fcct_biases( const double ra, const double dec, const char catalog,
                 const double jd, double *bias_ra, double *bias_dec);

#define AUTOMATIC_PERTURBERS  1
#define J2000 2451545.
#define JD_TO_YEAR( jd)  (2000. + ((jd) - J2000) / 365.25)
#define MAX_CONSTRAINTS 5
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
      /* GAUSS_K is a fixed constant. SOLAR_GM = 0.0002959122082855911025, */
      /* in AU^3/day^2; = 1.3271243994E+11 km3/s2 */


const double SRP1AU = 2.3e-7;
   /* "Solar radiation pressure at 1 AU",  in kg*AU^3 / (m^2*d^2),       */
   /* from a private communication from Steve Chesley.  This means       */
   /* that if you had a one-kilogram object showing one square           */
   /* meter of surface area to the sun,  and it was one AU from the      */
   /* sun,  and it absorbed all the solar radiation (and re-radiated     */
   /* it isotropically;  i.e.,  the re-radiated photons didn't cause     */
   /* any net force),  then that object would accelerate away from       */
   /* the sun at 2.3e-7 AU/day^2 (a.k.a. 4.562e-6 m/s^2).                */
   /*                                                                    */
   /*   One can derive SRP1AU  from basic principles.   The 'solar       */
   /* constant' is C = 1367.6 AU^2*W/m^2; i.e.,  if you set up a         */
   /* one square meter solar panel with 100% efficiency one AU from      */
   /* the Sun, pointed straight at the sun, it would generate            */
   /* 1367.6 watts. To get the above number, one uses                    */
   /*                                                                    */
   /* SRP1AU = C * d^2 / (meters_per_AU * c)                             */
   /*                                                                    */
   /*    C = 1367.6 AU^2*W/m^2 = 1367.6 AU^2*kg/s^3                      */
   /*    d = 86400 seconds/day                                           */
   /*    meters_per_AU = 1495978707 m/AU                                 */
   /*    c = 299792458 m/s                                               */
   /*                                                                    */
   /*    The result indicates that if the solar panel in question        */
   /* had a mass of one kilogram,  it would accelerate away from the     */
   /* sun at 2.27636e-7 AU/day^2.  I think Steve rounded off with a      */
   /* one-percent error because that's a decent match to the accuracy    */
   /* you can expect with real-world objects that reflect and re-radiate */
   /* photons,  instead of politely absorbing them and then re-radiating */
   /* them isotropically.                                                */
   /*                                                                    */
   /*   Interestingly,  this also means that an object with area/mass =  */
   /* 1286 m^2/kg would have SRP balancing the sun's gravity.  Which     */
   /* would make it big and light.  Solar sails aren't easy.             */
   /*
        A final comment : non-gravs of the A1, A2, A3 form give the
    acceleration the object would have at one AU from the sun,  in
    units of AU/day.  If the non-gravs are of the 1/r^2 model (rock-like)
    rather than,  say,  the Sekanina-Marsden model for an outgassing
    comet,  then for an area/mass ratio of z m^2/kg,  A1 = SRP1AU * z.
    Or,  alternatively,  A1 = SRP1AU * AMR.  (Roughly speaking,  anyway;
    in such scenarios,  A2 and A3 are fitted parameters,  so the radial
    component is slightly different.)  Thus,  for example,  one gets an
    area/mass ratio for 1I/`Oumuamua of 1.15 m^2/kg,  but if you fit A1
    and A2,  you get A1 = 2.39e-7 AU/day^2... somewhat lower than you'd
    expect if you just multiplied 1.15 by SRP1AU,  but pretty close;  the
    tangential components are rarely large.  */

int n_extra_params = 0, setting_outside_of_arc = 1;
double solar_pressure[MAX_N_NONGRAV_PARAMS], uncertainty_parameter = 99.;
int available_sigmas = NO_SIGMAS_AVAILABLE;
int available_sigmas_hash = 0;
static bool fail_on_hitting_planet = false;

double gaussian_random( void);                           /* monte0.c */
int get_residual_data( const OBSERVE *obs, double *xresid, double *yresid);
int debug_printf( const char *format, ...)                 /* runge.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
double initial_orbit( OBSERVE FAR *obs, int n_obs, double *orbit);
int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit);
int find_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
             const double r1, const double angle_param);   /* orb_func.cpp */
int search_for_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
              const double r1, double *angle_param);  /* orb_func.cpp */
int set_locs( const double *orbit, double t0, OBSERVE FAR *obs, int n_obs);
int find_best_fit_planet( const double jd, const double *ivect,
                                 double *rel_vect);         /* runge.cpp */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
static int evaluate_limited_orbit( const double *orbit,
                    const int planet_orbiting, const double epoch,
                    const char *limited_orbit, double *constraints);
void remove_trailing_cr_lf( char *buff);      /* ephem0.cpp */
int find_relative_orbit( const double jd, const double *ivect,
               ELEMENTS *elements, const int ref_planet);     /* runge.cpp */
void format_dist_in_buff( char *buff, const double dist_in_au); /* ephem0.c */
static inline void look_for_best_subarc( const OBSERVE FAR *obs,
       const int n_obs, const double max_arc_len, int *start, int *end);
int check_for_perturbers( const double t_cen, const double *vect); /* sm_vsop*/
int get_idx1_and_idx2( const int n_obs, const OBSERVE FAR *obs,
                                int *idx1, int *idx2);      /* elem_out.c */
void set_distance( OBSERVE FAR *obs, double r);             /* orb_func.c */
double find_r_given_solar_r( const OBSERVE FAR *obs, const double solar_r);
void attempt_extensions( OBSERVE *obs, const int n_obs, double *orbit,
                  const double epoch);                  /* orb_func.cpp */
double *get_asteroid_mass( const int astnum);   /* bc405.cpp */
char *get_file_name( char *filename, const char *template_file_name);
int compute_observer_loc( const double jde, const int planet_no,
             const double rho_cos_phi,           /* mpc_obs.cpp */
             const double rho_sin_phi, const double lon, double FAR *offset);
int compute_observer_vel( const double jde, const int planet_no,
             const double rho_cos_phi,           /* mpc_obs.cpp */
             const double rho_sin_phi, const double lon, double FAR *vel);
void get_relative_vector( const double jd, const double *ivect,
          double *relative_vect, const int planet_orbiting);  /* orb_func.c */
double get_planet_mass( const int planet_idx);                /* orb_func.c */
int compute_available_sigmas_hash( const OBSERVE FAR *obs, const int n_obs,
         const double epoch, const unsigned perturbers, const int central_obj);
double vector3_dist( const double *a, const double *b);     /* orb_func.c */
double euler_function( const OBSERVE FAR *obs1, const OBSERVE FAR *obs2);
double evaluate_initial_orbit( const OBSERVE FAR *obs,      /* orb_func.c */
                              const int n_obs, const double *orbit);
static int find_transfer_orbit( double *orbit, OBSERVE FAR *obs1,
                OBSERVE FAR *obs2,
                const int already_have_approximate_orbit);
double observation_rms( const OBSERVE FAR *obs);            /* elem_out.cpp */
double compute_weighted_rms( const OBSERVE FAR *obs, const int n_obs,
                           int *n_resids);                  /* orb_func.cpp */
double find_epoch_shown( const OBSERVE *obs, const int n_obs); /* elem_out */
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
const char *get_find_orb_text( const int index);      /* elem_out.cpp */
int snprintf_append( char *string, const size_t max_len,      /* ephem0.cpp */
                                   const char *format, ...)
#ifdef __GNUC__
         __attribute__ (( format( printf, 3, 4)))
#endif
;
void set_obs_vect( OBSERVE FAR *obs);        /* mpc_obs.h */
double improve_along_lov( double *orbit, const double epoch, const double *lov,
          const unsigned n_params, unsigned n_obs, OBSERVE *obs);
void adjust_error_ellipse_for_timing_error( double *sigma_a, double *sigma_b,
         double *angle, const double vx, const double vy);   /* errors.cpp */
void compute_error_ellipse_adjusted_for_motion( double *sigma1, double *sigma2,
                  double *posn_angle, const OBSERVE *obs,
                  const MOTION_DETAILS *m);                  /* orb_func.cpp */
double n_nearby_obs( const OBSERVE FAR *obs, const unsigned n_obs,
          const unsigned idx, const double time_span);       /* orb_func.cpp */
double find_parabolic_minimum_point( const double x[3], const double y[3]);
int orbital_monte_carlo( const double *orbit, OBSERVE *obs, const int n_obs,
         const double curr_epoch, const double epoch_shown);   /* orb_func.cpp */

void set_distance( OBSERVE FAR *obs, double r)
{
   int i;

   obs->r = r;
   for( i = 0; i < 3; i++)
      obs->obj_posn[i] = obs->obs_posn[i] + r * obs->vect[i];
   obs->solar_r = vector3_length( obs->obj_posn);
}

double vector3_dist( const double *a, const double *b)
{
   const double dx = a[0] - b[0];
   const double dy = a[1] - b[1];
   const double dz = a[2] - b[2];

   return( sqrt( dx * dx + dy * dy + dz * dz));
}

/* Euler found that for a parabolic orbit starting at distance r1 from
the sun,  ending up at a distance r2,  with a straight-line distance s
between them,  the time required can be found from

(r1 + r2 + s) ^ 1.5 +/- (r1 + r2 - s) ^ 1.5 = 6kt

   (minus sign if you're taking the "short route" -- sun isn't in the arc
you travel -- or positive sign if you go around the sun instead.)

   Our interest is in the "short route".  If the time for a transfer orbit
is less than that time,  only a hyperbolic orbit will get you from A to B
quickly enough.  Otherwise,  an elliptical orbit is possible.  (Similarly,
the time for the "long route" can tell you if an elliptical orbit taking
more than half a revolution is possible.)   */

double euler_function( const OBSERVE FAR *obs1, const OBSERVE FAR *obs2)
{
   const double s = vector3_dist( obs1->obj_posn, obs2->obj_posn);
   const double temp1 = obs1->solar_r + obs2->solar_r + s;
   const double temp2 = obs1->solar_r + obs2->solar_r - s;
   const double t1 = obs1->jd - obs1->r / AU_PER_DAY;
   const double t2 = obs2->jd - obs2->r / AU_PER_DAY;
   const double rval = temp1 * sqrt( temp1) - temp2 * sqrt( temp2)
                     - 6. * GAUSS_K * (t2 - t1);

   return( rval);
}

double find_r_given_solar_r( const OBSERVE FAR *obs, const double solar_r)
{
   double r_dot_v = 0., r_dot_r = 0., b, c, discr, rval = -1.;
   int i;

   for( i = 0; i < 3; i++)
      {
      r_dot_r += obs->obs_posn[i] * obs->obs_posn[i];
      r_dot_v += obs->obs_posn[i] * obs->vect[i];
      }
   b = 2. * r_dot_v;
   c = r_dot_r - solar_r * solar_r;
   discr = b * b - 4 * c;
   if( discr > 0.)
      rval = (-b + sqrt( discr)) / 2.;
   return( rval);
}

static double get_euler_value( const OBSERVE FAR *obs1, OBSERVE FAR *obs2,
            const double r2)
{
   set_distance( obs2, r2);
   return( euler_function( obs1, obs2));
}

/* The following function finds a zero of the quadratic passing through
the three points (x[0], y[0]), (x[1], y[1]), and (x[2], y[2]).  To do
this,  we first shift by x[0] along the x-axis,  to find the quadratic
y=ax^2 + bx + c that fits (0, y[0]), (x[1] - x[0], y[1]),
and (x[2] - x[0], y[2]).  By inspection,  in this scheme,  c = y[0].
Using the notation of the following function,  then,

dy1 = a * dx1^2 + b * dx1
dy2 = a * dx2^2 + b * dx2
dy1 / dx1 = z1 = a * dx1 + b
dy2 / dx2 = z2 = a * dx2 + b

   We then solve the resulting quadratic,  going for the positive sign
if dir != 0.   */


static double find_quadratic_zero( const double *x, const double *y,
               const int dir)
{
   const double dx1 = x[1] - x[0], dy1 = y[1] - y[0];
   const double dx2 = x[2] - x[0], dy2 = y[2] - y[0];
   const double z1 = dy1 / dx1, z2 = dy2 / dx2;
   const double a = (z1 - z2) / (dx1 - dx2);
   const double b = z1 - a * dx1;
   const double c = y[0];
   const double discriminant = b * b - 4. * a * c;
   double rval = -b;

   if( discriminant > 0.)
      rval += (dir ? 1. : -1.) * sqrt( discriminant);
   return( x[0] + rval / (2. * a));
}

int find_parabolic_orbit( OBSERVE FAR *obs, const int n_obs,
            double *orbit, const int direction)
{
   double r[3], y[3];
   int i, iter, rval = 0;
   OBSERVE FAR *obs2 = obs + n_obs - 1;
   const double thresh = 1e-10;

   for( i = 0; i < 3; i++)
      {
      r[i] = obs->r * (.9 + .1 * (double)i);
      y[i] = get_euler_value( obs, obs2, r[i]);
      }
   for( iter = 0; fabs( r[1] - r[0]) > thresh && iter < 20; iter++)
      {
      const double new_r = find_quadratic_zero( r, y, direction);

      r[2] = r[1];
      r[1] = r[0];
      r[0] = new_r;
      y[2] = y[1];
      y[1] = y[0];
      y[0] = get_euler_value( obs, obs2, new_r);
      }
   if( find_transfer_orbit( orbit, obs, obs2, 0))
      rval = -3;
   else if( set_locs( orbit, obs->jd, obs, n_obs))
      rval = -4;
   return( rval);
}

int calc_derivatives( const double jd, const double *ival, double *oval,
                           const int reference_planet);
long double take_rk_stepl( const long double jd, ELEMENTS *ref_orbit,
                 const long double *ival, long double *ovals,
                 const int n_vals, const long double step);     /* runge.cpp */
long double take_pd89_step( const long double jd, ELEMENTS *ref_orbit,
                 const long double *ival, long double *ovals,
                 const int n_vals, const long double step);    /* runge.cpp */
int symplectic_6( double jd, ELEMENTS *ref_orbit, double *vect,
                                          const double dt);
static int is_unreasonable_orbit( const double *orbit);     /* orb_func.cpp */
static int is_unreasonable_orbitl( const long double *orbit);

double integration_tolerance = 1.e-12;
double minimum_jd = 77432.5;      /* 1 Jan -4500 */
double maximum_jd = 4277757.5;    /* 1 Jan +7000 */

char *runtime_message;
int show_runtime_messages = 1;

// static int reference_planet = 0;
static unsigned perturbers_automatically_found;
extern unsigned always_included_perturbers;

static int reset_auto_perturbers( const double jd, const double *orbit)
{
   extern int forced_central_body;     /* and include asteroid perts   */
   unsigned mask;
   const int perturbing_planet = check_for_perturbers(
                                (jd - J2000) / 36525., orbit);

   mask = (1 << perturbing_planet);
   if( perturbing_planet == 3)    /* add in the moon,  too: */
      mask |= (1 << 10);
   else if( perturbing_planet == 10)    /* or vice versa:        */
      mask |= (1 << 3);
   if( forced_central_body == 100)
      mask |= (1 << 20);
   if( perturbing_planet)
      perturbers_automatically_found |= mask;
   if( forced_central_body == 100)
      perturbers_automatically_found |= (1 << 20);
   if( perturbers & AUTOMATIC_PERTURBERS)
      perturbers = mask | AUTOMATIC_PERTURBERS;
   perturbers |= always_included_perturbers;
   return( perturbing_planet);
}

clock_t integration_timeout = (clock_t)0;

#define STEP_INCREMENT 2
#define INTEGRATION_TIMED_OUT       -3
#define HIT_A_PLANET                -4

static void double_to_ldouble( long double *ovals, const double *ivals,
                                          size_t n)
{
   while( n--)
      *ovals++ = (long double)*ivals++;
}

static void ldouble_to_double( double *ovals, const long double *ivals,
                                          size_t n)
{
   while( n--)
      *ovals++ = (double)*ivals++;
}

int integrate_orbitl( long double *orbit, const long double t0, const long double t1)
{
   long double stepsize = 2.;
   static long double fixed_stepsize = -1.;
   const long double chicken = .9;
   int reset_of_elements_needed = 1;
   const long double step_increase = chicken * integration_tolerance
                 / powl( STEP_INCREMENT, (integration_method ? 9. : 5.));
   static int use_encke = -1;
   long double t = t0;
#ifdef CONSOLE
   static time_t real_time = (time_t)0;
   long double prev_t = t, last_err = 0.;
#endif
   int n_rejects = 0, rval;
   unsigned saved_perturbers = perturbers;
   int n_steps = 0, prev_n_steps = 0;
   int going_backward = (t1 < t0);
   static int n_changes;
   ELEMENTS ref_orbit;

   assert( fabsl( t0) < 1e+9);
   assert( fabsl( t1) < 1e+9);
   if( use_encke == -1)
      use_encke = atoi( get_environment_ptr( "ENCKE"));
   if( t0 > maximum_jd || t1 > maximum_jd
                       || t0 < minimum_jd || t0 < minimum_jd)
      {
      char buff[300];

      snprintf( buff, sizeof( buff), get_find_orb_text( 2027),
               (double)JD_TO_YEAR( t0), (double)JD_TO_YEAR( t1));
      generic_message_box( buff, "o");
      exit( -1);
      }
   if( debug_level > 7)
      debug_printf( "Integrating %f to %f\n", (double)t0, (double)t1);
   rval = is_unreasonable_orbitl( orbit);
   if( rval)
      {
      debug_printf( "Unreasonable %d\n", rval);
      return( -1);
      }
   ref_orbit.central_obj = -1;
   if( fixed_stepsize < 0.)
      fixed_stepsize = (long double)atof( get_environment_ptr( "FIXED_STEPSIZE"));
   if( fixed_stepsize > 0.)
      stepsize = fixed_stepsize;
   if( going_backward)
      stepsize = -stepsize;
   while( t != t1 && !rval)
      {
      long double delta_t, new_t = ceill( (t - .5) / stepsize + .5) * stepsize + .5;
      double dorbit[6];
      bool step_taken = true;

      ldouble_to_double( dorbit, orbit, 6);
      reset_auto_perturbers( t, dorbit);
      if( reset_of_elements_needed || !(n_steps % 50))
         if( use_encke)
            {
            extern int best_fit_planet;

            find_relative_orbit( t, dorbit, &ref_orbit, best_fit_planet);
            reset_of_elements_needed = 0;
            }
      n_steps++;
#ifdef CONSOLE
      if( !(n_steps % 500) && show_runtime_messages && time( NULL) != real_time)
         {
         char buff[80];
         extern int best_fit_planet, n_posns_cached;
         extern int64_t planet_ns;
         extern double best_fit_planet_dist;
// #define TEST_PLANET_CACHING_HASH_FUNCTION
#ifdef TEST_PLANET_CACHING_HASH_FUNCTION
         extern long total_n_searches, total_n_probes, max_probes_required;
#endif

         if( runtime_message)
            move_add_nstr( 9, 10, runtime_message, -1);
         sprintf( buff, "t = %.5f; %.5f to %.5f; step ",
                 (double)JD_TO_YEAR( t), (double)JD_TO_YEAR( t0), (double)JD_TO_YEAR( t1));
         if( fabsl( stepsize) > .1)
            snprintf_append( buff, sizeof( buff), "%.3f   ", (double)stepsize);
         else if( fabsl( stepsize) > .91)
            snprintf_append( buff, sizeof( buff), "%.3fm   ",
                                          (double)stepsize * minutes_per_day);
         else
            snprintf_append( buff, sizeof( buff), "%.3fs   ",
                                          (double)stepsize * seconds_per_day);
         if( prev_n_steps)    /* i.e.,  not our first time through here */
            snprintf_append( buff, sizeof( buff), "%d step/sec  ",
                        n_steps - prev_n_steps);
         move_add_nstr( 10, 10, buff, -1);
         prev_n_steps = n_steps;
         real_time = time( NULL);
         sprintf( buff, " %02d:%02d:%02d; %f; %d cached   ",
                     (int)( (real_time / 3600) % 24L),
                     (int)( (real_time / 60) % 60),
                     (int)( real_time % 60), (double)( t - prev_t),
                     n_posns_cached);
         prev_t = t;
         move_add_nstr( 11, 10, buff, -1);
         sprintf( buff, "%d steps; %d rejected", n_steps, n_rejects);
         if( best_fit_planet_dist)
            {
            snprintf_append( buff, sizeof( buff), "; center %d, ",
                            best_fit_planet);
            format_dist_in_buff( buff + strlen( buff), best_fit_planet_dist);
            }
         if( planet_ns)
            sprintf( buff + strlen( buff), "  tp:%ld.%09ld",
                  (long)( planet_ns / (int64_t)1000000000),
                  (long)( planet_ns % (int64_t)1000000000));
         strcat( buff, "  ");
         move_add_nstr( 12, 10, buff, -1);
         sprintf( buff, "last err: %.3e/%.3e  n changes: %d  ",
                        (double)last_err, (double)step_increase, n_changes);
         move_add_nstr( 13, 10, buff, -1);
         if( use_encke)
            {
            sprintf( buff, "e = %.5f; q = ", ref_orbit.ecc);
            format_dist_in_buff( buff + strlen( buff), ref_orbit.q);
            strcat( buff, "     ");
            move_add_nstr( 18, 10, buff, -1);
            }
         sprintf( buff, "Pos: %11.6f %11.6f %11.6f",
                     dorbit[0], dorbit[1], dorbit[2]);
         move_add_nstr( 14, 10, buff, -1);
         sprintf( buff, "Vel: %11.6f %11.6f %11.6f",
                     dorbit[3], dorbit[4], dorbit[5]);
         move_add_nstr( 15, 10, buff, -1);
#ifdef TEST_PLANET_CACHING_HASH_FUNCTION
         if( total_n_searches)
            {
            sprintf( buff, "%ld searches; avg %.2f max %ld     ",
                            total_n_searches,
                            (double)total_n_probes / (double)total_n_searches,
                            max_probes_required);
            move_add_nstr( 16, 10, buff, -1);
            }
#endif
         refresh_console( );
         }
#endif

               /* Make sure we don't step completely past */
               /* the time t1 we want to stop at!         */
      if( (!going_backward && new_t > t1) || (going_backward && new_t < t1))
         new_t = t1;
      delta_t = new_t - t;

      switch( integration_method)
         {
#ifdef NOT_READY_FOR_LONG_DOUBLES
         case 1:
            symplectic_6( t, &ref_orbit, orbit, delta_t);
            break;
#endif
         case 0:
         default:
            {
            long double new_vals[6];
            static long double min_stepsize;
            const double err = (integration_method ?
                   take_pd89_step( t, &ref_orbit, orbit, new_vals, 6, delta_t) :
                   take_rk_stepl( t, &ref_orbit, orbit, new_vals, 6, delta_t));

            if( !min_stepsize)
               {
               min_stepsize = (long double)atof( get_environment_ptr( "MIN_STEPSIZE"))
                                             / seconds_per_day;
               if( !min_stepsize)
                  min_stepsize = 1e-5;   /* 1e-5 day = 0.864 seconds */
               }
            if( !stepsize)
               exit( -6);
//          if( err >= integration_tolerance && fabs( stepsize) < min_stepsize)
//             debug_printf( "Err %f x tolerance; stepsize %f seconds\n",
//                         err / integration_tolerance, delta_t * seconds_per_day);
            if( err < integration_tolerance || fixed_stepsize > 0.
                        || fabs( stepsize) < min_stepsize)  /* it's good! */
               {
               memcpy( orbit, new_vals, 6 * sizeof( long double));
               if( err < step_increase && !fixed_stepsize)
                  if( fabsl( delta_t - stepsize) < fabsl( stepsize * .01))
                     {
                     n_changes++;
                     stepsize *= STEP_INCREMENT;
                     }
               }
            else           /* failed:  try again with a smaller step */
               {
               n_rejects++;
               step_taken = false;
               new_t = t;
               stepsize /= STEP_INCREMENT;
               reset_of_elements_needed = 1;
               }
#ifdef CONSOLE
            last_err = err;
#endif
            }
            break;
         }
      t = new_t;
      rval = is_unreasonable_orbitl( orbit);
      if( rval)
         {
         debug_printf( "Unreasonable %d at %.5g (%.5f to %.5f)\n",
                  rval, (double)(t - t0), (double)t0, (double)t1);
         debug_printf( "Stepsize %g\n", (double)stepsize);
         }
      else if( integration_timeout && !(n_steps % 100))
         if( clock( ) > integration_timeout)
            rval = INTEGRATION_TIMED_OUT;
      if( step_taken && fail_on_hitting_planet)
         {
         extern int planet_hit;

         if( planet_hit != -1)
            rval = HIT_A_PLANET;
         }
#ifdef PROBABLY_UNNEEDED
      if( debug_level && n_steps % 10000 == 0)
         {
         extern int best_fit_planet;
         ELEMENTS elems;

         find_relative_orbit( t, orbit, &elems, best_fit_planet);
         debug_printf( "At %f, step %g near planet %d\n", t, stepsize, best_fit_planet);
         debug_printf( "q = %f km; e = %f\n", elems.q * AU_IN_KM, elems.ecc);
         debug_printf( "Posn: %f %f %f\n", orbit[0], orbit[1], orbit[2]);
         debug_printf( "Posn: %f %f %f\n", orbit[3], orbit[4], orbit[5]);
         }
#endif
      }
   if( debug_level > 7)
      debug_printf( "Integration done: %d\n", rval);
   perturbers = saved_perturbers;
   return( rval);
}

int integrate_orbit( double *orbit, const double t0, const double t1)
{
   long double tarray[6];
   int rval;

   double_to_ldouble( tarray, orbit, 6);
   rval = integrate_orbitl( tarray, (long double)t0, (long double)t1);
   ldouble_to_double( orbit, tarray, 6);
   return( rval);
}

/* At times,  the orbits generated by 'full steps' or Herget or other methods
   are completely unreasonable.  The exact definition of 'unreasonable'
   is pretty darn fuzzy.  The following function says that if at the epoch,
   an object is more than one light-year from the Sun (about 60000 AU) or
   is travelling at faster than 5% of the speed of light,  the orbit is
   "unreasonable".
      I'm sure this definition could be tightened a lot without any real
   trouble.  The fastest natural objects in the solar system are comets
   impacting the Sun,  which they do at the escape speed of about 600 km/s,
   or .2% of the speed of light.  (With the limit being a lot tighter than
   that as one gets away from the sun;  for example,  near the earth's
   orbit,  nothing goes much faster than about 70 km/s.)  And the furthest
   objects seen orbiting the sun are of the order of 100 AU away,  as
   opposed to the 60000 AU distance given below.   */

static const double max_reasonable_dist = 1.0 * 365.25 * AU_PER_DAY;
static const double max_reasonable_speed = AU_PER_DAY * .05;

bool all_reasonable = false;

#define UNREASONABLE_TOO_FAR      0x100
#define UNREASONABLE_TOO_FAST     0x200
#define UNREASONABLE_ZERO_R       0x400
#define UNREASONABLE_ZERO_V       0x800

static int is_unreasonable_orbit( const double *orbit)
{
   int rval = 0, i;
   double r, v;

   if( all_reasonable)
      return( 0);
   r = vector3_length( orbit);
   v = vector3_length( orbit + 3);
   if( r > max_reasonable_dist)
      rval = UNREASONABLE_TOO_FAR;
   else if( !r)
      rval = UNREASONABLE_ZERO_R;
   if( v > max_reasonable_speed)
      rval |= UNREASONABLE_TOO_FAST;
   else if( !v)
      rval |= UNREASONABLE_ZERO_V;
   for( i = 0; i < 6; i++)
      if( !isfinite( orbit[i]))
         rval |= (1 << i);
   return( rval);
}

static int is_unreasonable_orbitl( const long double *orbit)
{
   double tarray[6];

   ldouble_to_double( tarray, orbit, 6);
   return( is_unreasonable_orbit( tarray));
}

/* See Explanatory Supplement,  3.26, p. 135, "Gravitational Light
Bending."  For our purposes,  what matters is the difference between
how much the object's light is bent and how much the light of
background stars is bent.  So we compute phi1 = angle between
observer,  sun,  and 'result';  and phi2 = angle between observer,
sun,  and background stars = 180 minus elongation of the object as
seen by 'observer'.  */

static void light_bending( const double *observer, double *result)
{
   const double bend_factor = 2. * SOLAR_GM / (AU_PER_DAY * AU_PER_DAY);
   size_t i;
   double p[3], plen, xprod[3], dir[3], dlen;
   const double olen = vector3_length( observer);
   const double rlen = vector3_length( result);
   double phi1, phi2, bending;

   for( i = 0; i < 3; i++)
      p[i] = result[i] - observer[i];
   plen = vector3_length( p);
   if( plen < 1e-10)  /* don't do light-bending over really */
      return;         /* short/meaningless distances */
   vector_cross_product( xprod, observer, result);
   vector_cross_product( dir, xprod, p);
   dlen = vector3_length( dir);
   if( !dlen)
      return;
   for( i = 0; i < 3; i++)
      dir[i] /= dlen;
     /* "dir" is now a unit vector perpendicular to p,  aimed away */
     /* from the sun */
   phi1 = acose( dot_product( result, observer) / (rlen * olen));
   phi2 = acose( -dot_product( p, observer) / (plen * olen));
   bending = bend_factor * (tan( phi1 / 2.) - tan( phi2 / 2.));
   bending *= plen;
   for( i = 0; i < 3; i++)
      result[i] += bending * dir[i];
}

void light_time_lag( const double *orbit, const double *observer, double *result)
{
   const double solar_r = vector3_length( orbit);
   unsigned iter;

   memcpy( result, orbit, 3 * sizeof( double));
   for( iter = 0; iter < 4; iter++)
      {
      const double r = vector3_dist( result, observer);
      const double dt = -r / AU_PER_DAY;
      const double afact = -SOLAR_GM * dt / (solar_r * solar_r * solar_r);
      unsigned i;

      for( i = 0; i < 3; i++)
         {
         result[i]     = orbit[i] + (.5 * afact * orbit[i] + orbit[i + 3]) * dt;
         result[i + 3] = orbit[i + 3] + afact * orbit[i];
         }
//    debug_printf( "iter %d: r = %.15f\n", iter, r);
      }
   light_bending( observer, result);
}

static void set_solar_r( OBSERVE FAR *ob)
{
   ob->solar_r = vector3_length( ob->obj_posn);
}

/* In integrating an orbit to compute locations for each observation,  a
little bit of trickery is used to improve performance.  First (pass = 0),
we set observations made _before_ the epoch t0,  integrating _backwards_.
If we instead,  say,  integrated from the epoch t0 to the first observation,
then integrated to each observation in order,  it would simplify the code
a bit;  but we'd sometimes be integrating (say) from an epoch in 2010 back
to 2003,  then forward again to observations made up to 2013.  That would
mean integrating over a total span of 17 years when we really only needed
to do ten years.  It would also allow numerical integration error to
build up more.

   To make matters worse,  while the initial state vector may be for 2010,
we may also want a state vector for an epoch in,  say,  2018.  (Perhaps
that's the epoch we're displaying,  or the semimajor axis is constrained
to a particular value for sometime in 2018.)  'set_locs_extended' is
bright enough to compute that second-epoch state vector at the optimal
time;  in the above case,  that would mean that after integrating forward
to 2013,  it should continue to 2018.  If the second epoch were between
2010 and 2013,  it would find the second-epoch state vector while getting
locations for the last part of the observed arc.  And so on.

   All of this does a good job of avoiding unnecessary integration and
accumulated error,  but it does make the code more opaque than would
otherwise be the case.

   If we don't actually need a state vector for a second epoch,  then
we can set orbit2 = NULL.  Or use plain ol' set_locs(),  which -- as
you can see below -- basically just calls set_locs_extended() with a
NULL orbit2.                          */

#define is_between( t1, t2, t3)  ((t2 - t1) * (t3 - t2) >= 0.)

static int set_locs_extended( const double *orbit, const double epoch_jd,
                       OBSERVE FAR *obs, const int n_obs,
                       const double epoch2, double *orbit2)
{
   int i, pass, rval = 0;

   if( is_unreasonable_orbit( orbit))
      {
      if( debug_level)
         debug_printf( "Unreasonable orbit provided to set_locs_extended: %s\n",
                        obs->packed_id);
      return( -9);
      }

   for( i = 0; i < n_obs && obs[i].jd < epoch_jd; i++)
      ;

               /* set obs[0...i-1] on pass=0, obs[i...n_obs-1] on pass=1: */
   for( pass = 0; pass < 2; pass++)
      {
      int j = (pass ? i : i - 1);
      long double curr_orbit[6];
      double curr_t = epoch_jd;

      double_to_ldouble( curr_orbit, orbit, 6);
      while( j < n_obs && j >= 0)
         {
         double light_lagged_orbit[6], temp_orbit[6];
         OBSERVE FAR *optr = obs + j;

         if( orbit2 && is_between( curr_t, epoch2, optr->jd))
            {
            rval = integrate_orbitl( curr_orbit, curr_t, epoch2);
            if( rval)
               return( rval);
            ldouble_to_double( orbit2, curr_orbit, 6);
            curr_t = epoch2;
            }
         rval = integrate_orbitl( curr_orbit, curr_t, optr->jd);
         if( rval)
            return( rval);
         curr_t = optr->jd;
         ldouble_to_double( temp_orbit, curr_orbit, 6);
         light_time_lag( temp_orbit, optr->obs_posn, light_lagged_orbit);
         FMEMCPY( optr->obj_posn, light_lagged_orbit, 3 * sizeof( double));
         FMEMCPY( optr->obj_vel, light_lagged_orbit + 3, 3 * sizeof( double));
         j += (pass ? 1 : -1);
         }
      if( orbit2)
         if( (!pass && curr_t >= epoch2) || (pass && curr_t <= epoch2))
            {
            rval = integrate_orbitl( curr_orbit, curr_t, epoch2);
            if( rval)
               return( rval);
            ldouble_to_double( orbit2, curr_orbit, 6);
            }
      }

            /* We've now set the object heliocentric positions and */
            /* velocities,  in ecliptic J2000,  for each observation */
            /* time.  Now let's go back and find observer-centric */
            /* computed RA/decs and distances to the object at those */
            /* times. */
   for( i = 0; i < n_obs; i++)
      {
      double loc[3], ra, dec, temp, r = 0.;
      const double sin_obliq_2000 = 0.397777155931913701597179975942380896684;
      const double cos_obliq_2000 = 0.917482062069181825744000384639406458043;
      int j;

      for( j = 0; j < 3; j++)
         {
         loc[j] = obs[i].obj_posn[j] - obs[i].obs_posn[j];
         r += loc[j] * loc[j];
         }
      r = sqrt( r);
      obs[i].r = r;
      temp = loc[1] * cos_obliq_2000 - loc[2] * sin_obliq_2000;
      loc[2] = loc[2] * cos_obliq_2000 + loc[1] * sin_obliq_2000;
      loc[1] = temp;
      ra = atan2( loc[1], loc[0]);
      if( r > 100000. || r <= 0.)
         debug_printf( "???? bad r: %f %f %f: %s\n",
                  loc[0], loc[1], loc[2], obs->packed_id);
      if( r)
         dec = asine( loc[2] / r);
      else
         dec = 0.;
      while( ra - obs[i].ra > PI)
         ra -= 2. * PI;
      while( ra - obs[i].ra < -PI)
         ra += 2. * PI;
      obs[i].computed_ra = ra;
      obs[i].computed_dec = dec;
      set_solar_r( obs + i);
      }
   return( 0);
}

int set_locs( const double *orbit, const double t0, OBSERVE FAR *obs,
                       const int n_obs)
{
   return( set_locs_extended( orbit, t0, obs, n_obs, t0, NULL));
}

double observation_rms( const OBSERVE FAR *obs)
{
   const double d_dec = obs->computed_dec - obs->dec;
   const double d_ra  = (obs->computed_ra  - obs->ra) * cos( obs->computed_dec);

   return( sqrt( d_dec * d_dec + d_ra * d_ra) * 3600. * 180. / PI);
}

/* 2010 Nov 4:  revised so that RMS residuals are computed as
root-mean-square of RA and dec treated separately,  meaning the
previous results needed to be multipled by sqrt(.5)  (Gareth
Williams kindly steered me the right way on this).          */

double compute_rms( const OBSERVE FAR *obs, const int n_obs)
{
   double rval = 0.;
   int i, n_included = 0;

   for( i = n_obs; i; i--, obs++)
      if( obs->is_included && obs->note2 != 'R')
         {
         const double obs_rms = observation_rms( obs);

         rval += obs_rms * obs_rms;
         n_included++;
         }
   if( !n_included)           /* avoid NaN */
      n_included = 1;
   return( sqrt( rval / (double)(n_included * 2)));
}

double compute_weighted_rms( const OBSERVE FAR *obs, const int n_obs, int *n_resids)
{
   double rval = 0;
   int i, n = 0;

   for( i = n_obs; i; i--, obs++)
      if( obs->is_included)
         {
         double xresid, yresid;

         n += get_residual_data( obs, &xresid, &yresid);

         rval += xresid * xresid + yresid * yresid;
         }
   if( n_resids)
      *n_resids = n;
   if( !n)           /* avoid NaN */
      n = 1;
   return( sqrt( rval / (double)n));
}

static double eval_3x3_determinant( const double *a, const double *b,
                         const double *c)
{
   return( a[0] * (b[1] * c[2] - c[1] * b[2])
         + b[0] * (c[1] * a[2] - a[1] * c[2])
         + c[0] * (a[1] * b[2] - b[1] * a[2]));
}

/* 'find_transfer_orbit' finds the state vector that can get an object
from the location/time described at 'obs1' to that described at 'obs2'.
It does this with the logic given in 'herget.htm#sund_xplns'. */
/* 2011 May 8:  Realized this code is responsible for much of the time
   it takes to initially compute an orbit.  One simple speed-up is to
   take advantage of the fact that in the method of Herget,  we're
   often tweaking a "nearby" orbit,  i.e.,  we already have an
   approximate orbit;  using this as our starting estimate ought to
   result in faster convergence.       */

#define XFER_TOO_FAR_AWAY          -1
#define XFER_TOO_FAST              -2
#define XFER_TOO_MANY_ITERATIONS   -3
#define XFER_INTEGRATION_FAILED    -4
#define XFER_OK                     0

clock_t t_transfer;

static int find_transfer_orbit( double *orbit, OBSERVE FAR *obs1,
                OBSERVE FAR *obs2,
                const int already_have_approximate_orbit)
{
   double r = 0.;
   const double jd1 = obs1->jd - obs1->r / AU_PER_DAY;
   const double jd2 = obs2->jd - obs2->r / AU_PER_DAY;
   const double delta_t = jd2 - jd1;
   double deriv[6];
   double orbit2[6];
   double diff_squared = 999.;
            /* Iterate until the error is less than 1e-8 of the object */
            /* distance.  This corresponds to an error of about .002 arcsec */
            /* in the actual position (with the '+ .01' allowing for */
            /* some margin for very close objects,  such as artsats). */
            /*   Assume a maximum of ten iterations,  just to be safe */
   const double target_diff = 1.e-8 * (obs2->r + .01);
   unsigned i, max_iterations = 10;
   clock_t t0 = clock( );

   assert( fabs( obs1->jd) < 1e+9);
   assert( fabs( obs2->jd) < 1e+9);
   assert( fabs( jd1) < 1e+9);
   assert( fabs( jd2) < 1e+9);
   if( obs1->r > max_reasonable_dist || obs2->r > max_reasonable_dist)
      {
      debug_printf( "Bad xfer: %f %f\n", obs1->r, obs2->r);
      return( XFER_TOO_FAR_AWAY);
      }

   set_distance( obs1, obs1->r);
   set_distance( obs2, obs2->r);

   if( !already_have_approximate_orbit)
      {
      double speed_squared = 0., speed;
      const int saved_perturbers = perturbers;
      double deriv2[6];
      unsigned pass;

      for( i = 0; i < 3; i++)
         {
         orbit[i + 3] = (obs2->obj_posn[i] - obs1->obj_posn[i]) / delta_t;
         speed_squared += orbit[i + 3] * orbit[i + 3];
         }

      if( speed_squared > max_reasonable_speed * max_reasonable_speed)
         return( XFER_TOO_FAST);
                     /* speed_squared is in (AU/day)^2, speed in km/second: */
      speed = sqrt( speed_squared) * AU_IN_KM / seconds_per_day;
      if( debug_level > 3)
         debug_printf( "Speed = %f km/s; radii %f, %f\n", speed, obs1->r, obs2->r);
      for( pass = 0; pass < 2; pass++)
         {
         OBSERVE FAR *obs = (pass ? obs2 : obs1);
         const double jd = (pass ? jd2 : jd1);

         assert( fabs( jd) < 1e+9);
         for( i = 0; i < 3; i++)
            orbit[i] = obs->obj_posn[i];
         reset_auto_perturbers( jd, orbit);
         calc_derivatives( jd, orbit, (pass ? deriv2 : deriv), -1);
         perturbers = saved_perturbers;
         }
      for( i = 0; i < 3; i++)
         orbit[3 + i] -= delta_t * (deriv[i + 3] / 3 + deriv2[i + 3] / 6.);
//       orbit[3 + i] -= deriv[i + 3] * delta_t / 2.
//                 + (deriv2[i + 3] - deriv[i + 3]) * delta_t / 6.;
      }
   for( i = 0; i < 3; i++)
      orbit[i] = obs1->obj_posn[i];

   for( i = 0; i < 3 && diff_squared > target_diff * target_diff; i++)
      {
      unsigned j;

      memcpy( orbit2, orbit, 6 * sizeof( double));
      if( integrate_orbit( orbit2, jd1, jd2))
         return( XFER_INTEGRATION_FAILED);
      diff_squared = 0;
      for( j = 0; j < 3; j++)
         {
         const double delta = orbit2[j] - obs2->obj_posn[j];

         orbit[j + 3] -=  delta / delta_t;
         diff_squared += delta * delta;
         }
      if( debug_level > 3)
         debug_printf( "Initial pass %d: %.3g%s\n", i, diff_squared,
                  (i || already_have_approximate_orbit) ? "" : " Using new method");
      }

   while( diff_squared > target_diff * target_diff && max_iterations--)
      {
      double delta[4][3], discr;
      unsigned pass;

      for( pass = 0; pass < 4; pass++)
         {
         const double h = 1.e-5;

         memcpy( orbit2, orbit, 6 * sizeof( double));
         if( pass)
            orbit2[(pass - 1) + 3] += h;
         if( integrate_orbit( orbit2, jd1, jd2))
            return( XFER_INTEGRATION_FAILED);
         for( i = 0; i < 3; i++)
            orbit2[i] -= obs2->obj_posn[i];
         memcpy( delta[pass], orbit2, 3 * sizeof( double));
         if( pass)
            for( i = 0; i < 3; i++)
               delta[pass][i] = (delta[pass][i] - delta[0][i]) / h;
         }
      discr = eval_3x3_determinant( delta[1], delta[2], delta[3]);
      if( discr)
         {
         orbit[3] -= eval_3x3_determinant( delta[0], delta[2], delta[3]) / discr;
         orbit[4] -= eval_3x3_determinant( delta[1], delta[0], delta[3]) / discr;
         orbit[5] -= eval_3x3_determinant( delta[1], delta[2], delta[0]) / discr;
         }
      diff_squared = 0.;
      for( i = 0; i < 3; i++)
         diff_squared += delta[0][i] * delta[0][i];
      if( debug_level > 3)
         debug_printf( "Transfer orbit: %d, %.3g; discr=%g\n", max_iterations,
                        diff_squared, discr);
      }
   if( !max_iterations)
      return( XFER_TOO_MANY_ITERATIONS);
               /* adjust for light-time lag: */
   r = -SOLAR_GM / (obs1->solar_r * obs1->solar_r * obs1->solar_r);
   for( i = 0; i < 3; i++)
      {
      const double time_lag = obs1->r / AU_PER_DAY;

      orbit[i] += time_lag * (orbit[3 + i] + orbit[i] * r * time_lag / 2.);
      orbit[i + 3] += orbit[i] * r * time_lag;
      }
   t_transfer += clock( ) - t0;
   return( XFER_OK);
}

static int find_parameterized_orbit( double *orbit, const double *params,
                OBSERVE obs1, OBSERVE obs2, const unsigned parameter_type,
                const int already_have_approximate_orbit)
{
   const double *ra_dec_offsets = NULL;
   int rval;

   switch( parameter_type)
      {
      case FIT_CLASSIC_HERGET:
      case FIT_HERGET_FULL:
         obs1.r *= 1. + params[0] * 100.;
         obs2.r *= 1. + params[1] * 100.;
         if( parameter_type == FIT_HERGET_FULL)
            ra_dec_offsets = params + 2;
         break;
      case FIT_FIXED_DISTANCES:
         ra_dec_offsets = params;
         break;
      case FIT_VAISALA_FULL:
         ra_dec_offsets = params + 1;
         break;
               /* figure out parabolics & circulars later */
      default:
         break;
      }
   if( ra_dec_offsets)
      {
      obs1.ra += *ra_dec_offsets++;
      obs1.dec += *ra_dec_offsets++;
      set_obs_vect( &obs1);
      obs2.ra  += *ra_dec_offsets++;
      obs2.dec += *ra_dec_offsets++;
      set_obs_vect( &obs2);
      }
   switch( parameter_type)
      {
      case FIT_VAISALA:
      case FIT_VAISALA_FULL:
         obs1.solar_r *= 1. + params[0] * 100.;
         obs1.r = find_r_given_solar_r( &obs1, obs1.solar_r);
         obs2.r = find_r_given_solar_r( &obs2, obs1.solar_r);
         break;
      default:
         break;
      }
   set_distance( &obs1, obs1.r);
   set_distance( &obs2, obs2.r);
   rval = find_transfer_orbit(  orbit, &obs1, &obs2,
                          already_have_approximate_orbit);
   return( rval);
}

int find_vaisala_orbit( double *orbit, const OBSERVE *obs1,
                     const OBSERVE *obs2, const double solar_r)
{
   double ignored_params[1];
   OBSERVE tobs = *obs1;

   ignored_params[0] = 0.;
   tobs.solar_r = solar_r;
   return( find_parameterized_orbit( orbit, ignored_params,
               tobs, *obs2, FIT_VAISALA, 0));
}

int drop_excluded_obs( OBSERVE *obs, int *n_obs)
{
   int rval = 0;

   while( *n_obs && !obs->is_included)    /* skip excluded obs at */
      {                                   /* start of arc */
      rval++;
      obs++;
      (*n_obs)--;
      }
   while( *n_obs && !obs[*n_obs - 1].is_included)   /* skip excluded obs at */
      (*n_obs)--;                                   /* end of arc */
   return( rval);
}

int extended_orbit_fit( double *orbit, OBSERVE *obs, int n_obs,
                  const unsigned fit_type, double epoch)
{
   int i, j, rval = 0, n_resids;
   int n_selected;
   const int n_params = (int)( fit_type & 0xf);
   double orbit_at_epoch[6];
   void *lsquare;
   double *resids, *slopes;
   double params[10];
   const double delta_val = 1e-6;
   OBSERVE obs1, obs2;

   obs += drop_excluded_obs( obs, &n_obs);
   n_resids = 2 * n_obs + MAX_CONSTRAINTS;
   resids = (double *)calloc( n_resids * (n_params + 1), sizeof( double));
   slopes = resids + n_resids;
   for( i = n_selected = 0; i < n_obs; i++)
      if( obs[i].flags & OBS_IS_SELECTED)
         n_selected++;
   if( n_selected == 2)
      {
      for( i = n_selected = 0; i < n_obs; i++)
         if( obs[i].flags & OBS_IS_SELECTED)
            {
            if( !n_selected)
               obs1 = obs[i];
            else
               obs2 = obs[i];
            n_selected++;
            }
      }
   else
      {
      obs1 = obs[0];
      obs2 = obs[n_obs - 1];
      }
   integrate_orbit( orbit, epoch, obs1.jd);
   for( i = -1; i < n_params; i++)
      {
      for( j = 0; j < n_params; j++)
         params[j] = 0.;
      if( i >= 0)
         params[i] = -delta_val;
      rval = find_parameterized_orbit( orbit, params, obs1, obs2,
                     fit_type, 0);
//    if( rval)
//       printf( "Iter %d, rval %d\n", i, rval);
      set_locs_extended( orbit, obs1.jd, obs, n_obs, epoch, orbit_at_epoch);
      for( j = 0; j < n_obs; j++)
         if( obs[j].is_included)
            {
            double dx, dy;

            get_residual_data( obs + j, &dx, &dy);
//          printf( "Resids: obs %d, param %d, %f %f\n", j, i, dx, dy);
            if( i == -1)
               {
               resids[j + j] = dx;
               resids[j + j + 1] = dy;
               }
            else
               {
               dx -= resids[j + j];
               dy -= resids[j + j + 1];
               slopes[(j + j    ) * n_params + i] = dx / delta_val;
               slopes[(j + j + 1) * n_params + i] = dy / delta_val;
               }
            }
      }

   lsquare = lsquare_init( n_params);
   assert( lsquare);
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         double weight = 1.;

         for( j = i + i; j < i + i + 2; j++)
            lsquare_add_observation( lsquare, resids[j], weight, slopes + j * n_params);
         }

   free( resids);
   lsquare_solve( lsquare, params);
// for( i = 0; i < n_params; i++)
//    printf( "Param %d: %.9f\n", i, params[i]);
   lsquare_free( lsquare);
   if( fit_type == FIT_CLASSIC_HERGET || fit_type == FIT_HERGET_FULL)
      {
      double rescale_adjustment = 1., tval;

      for( i = 0; i < 2; i++)
         {
         tval = fabs( params[i]) * 300.;
         if( rescale_adjustment < tval)
            rescale_adjustment = tval;
         }
      for( i = 0; i < n_params; i++)
         params[i] /= rescale_adjustment;
      }
   rval = find_parameterized_orbit( orbit, params, obs1, obs2,
                     fit_type, 0);
   set_locs_extended( orbit, obs1.jd, obs, n_obs, epoch, orbit_at_epoch);
            /* Except we really want to return the orbit at epoch : */
   memcpy( orbit_at_epoch, orbit, 6 * sizeof( double));
   integrate_orbit( orbit, obs1.jd, epoch);
   return( rval);
}

/* 'find_trial_orbit' tries to find an orbit linking the first and last
(included) observations in 'obs',  with the distance to the first object
being r1 AU and a radial velocity defined by 'angle_param'.

   First,  excluded observations at the beginning and end are skipped over
by advancing the obs pointer and/or decrementing n_obs.  If there are at
least two observations remaining,  we set the first observation to be at
distance r1.  Next,  we compute the distance in space between that point
and the ray defined by the second observation,  and we also compute the
maximum distance the object could go (assuming a parabolic orbit with
the object at escape velocity) during the time 'dt' between the two
observations.

   If the object can't get to the ray from the first observation without
going faster than escape velocity,  then we return an error code of -2.
This basically says,  "There are no non-hyperbolic orbits that satisfy
these two observations with the value you gave for r1.  Sorry."

   If the object _can_ get there,  we compute an orbit with an assumed
radial velocity linking the two observations.  If angle_param == -1,
the orbit will be one in which the object is at escape speed,  going
_toward_ us.  If angle_param == 1,  it'll be at escape speed,  going
_away_ from us.  In between,  you'll get assorted elliptical orbits.
(If angle_param > 1 or less than -1,  you'll get an hyperbolic orbit.)
*/

int find_sr_ranges( double *ranges, const double *q1, const double *p1,
                                    const double *q2, const double *p2,
                                    const double gm, const double dt);

#ifndef min
   #define min( x, y) ((x) > (y) ? (y) : (x))
#endif

int find_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
                 const double r1, const double angle_param)
{
   int rval = 0;

   uncertainty_parameter = 99.;
   obs += drop_excluded_obs( obs, &n_obs);
   if( n_obs < 2)       /* need at least two valid obs */
      rval = -1;
   else
      {
      double r2 = 0., dist2 = 0., escape_dist2;
      OBSERVE FAR *endptr = obs + n_obs - 1;
      const double dt = endptr->jd - obs->jd;
      int i;

      set_distance( obs, r1);
      for( i = 0; i < 3; i++)
         r2 += (obs->obj_posn[i] - endptr->obs_posn[i]) * endptr->vect[i];
      set_distance( endptr, r2);
      for( i = 0; i < 3; i++)
         {
         const double delta = endptr->obj_posn[i] - obs->obj_posn[i];

         dist2 += delta * delta;
         }
               /* Escape velocity for an object r AU from the sun would be */
               /* sqrt( 2 * SOLAR_GM / r).  We'll use the square of this:  */
      escape_dist2 = 2. * SOLAR_GM / min( obs->solar_r, endptr->solar_r);
               /* ...and multiply by dt squared to get the square of the   */
               /* distance the object would travel if it's at the esc speed: */
      escape_dist2 *= dt * dt;
      if( dist2 > escape_dist2)     /* only hyperbolic orbits exist */
         {
         dist2 = escape_dist2;
         rval = -2;
         }
  /* else         _Used_ to be 'else' */
         {
         set_distance( endptr, r2 + angle_param * sqrt( escape_dist2 - dist2));
         if( find_transfer_orbit( orbit, obs, endptr, 0))
            rval = -3;
         else if( set_locs( orbit, obs->jd, obs, n_obs))
            {
            debug_printf( "Set_loc fail 1\n");
            rval = -4;
            }
         else if( n_obs > 2)
            rval = adjust_herget_results( obs, n_obs, orbit);
         }
      }
   return( rval);
}

double find_parabolic_minimum_point( const double x[3], const double y[3])
{
   const double x1 = x[1] - x[0], x2 = x[2] - x[0];
   const double y1 = y[1] - y[0], y2 = y[2] - y[0];
   const double bottom = y1 * x2 - y2 * x1;
   const double top    = y1 * x2 * x2 - y2 * x1 * x1;
   const double new_x = x[0] + top / (2. * bottom);

   return( new_x);
}

int search_for_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
                 const double r1, double *angle_param)
{
   int rval = 0;
   unsigned i, n_minima = 0;
   const unsigned n_divisions = 10;
   double best_found = 0., best_rms_found = 0.;
   double rms[3];

   rms[1] = rms[2] = 0.;  /* needed only to suppress a bogus g++ warning */
   for( i = 0; i <= n_divisions; i++)
      {
      double ang_param = 2. * (double)i / (double)n_divisions - 1.;

      find_trial_orbit( orbit, obs, n_obs, r1, ang_param);
      rms[0] = rms[1];
      rms[1] = rms[2];
      rms[2] = compute_weighted_rms( obs, n_obs, NULL);
//    debug_printf( "r1 %f, ang_param %f: %f\n",
//                r1, ang_param, rms[2]);
      if( !i || best_rms_found > rms[2])
         {
         best_rms_found = rms[2];
         best_found = ang_param;
         }
      if( i > 1 && rms[1] < rms[2] && rms[1] < rms[0])
         {                             /* true parabolic minimum found */
         double x[3], y[3];
         unsigned j, iteration;

         for( j = 0; j < 3; j++)
            {
            y[j] = rms[j];
            x[j] = ang_param - 2. * (double)(2 - j) / (double)n_divisions;
            }
         for( iteration = 0; iteration < 4; iteration++)
            {
            double new_rms;
            unsigned max_idx = 0;
            const double new_x = find_parabolic_minimum_point( x, y);

            debug_printf ("x: %f %f %f; y: %f %f %f\n",
                     x[0], x[1], x[2], y[0], y[1], y[2]);
            find_trial_orbit( orbit, obs, n_obs, r1, new_x);
            new_rms = compute_weighted_rms( obs, n_obs, NULL);
            if( y[1] > y[0])
               max_idx = 1;
            if( y[2] > y[max_idx])
               max_idx = 2;
            x[max_idx] = new_x;
            y[max_idx] = new_rms;
            if( new_rms < best_rms_found)
               {
               best_rms_found = new_rms;
               best_found = new_x;
               }
            debug_printf( "Iteration %u: x=%f, rms %f\n", iteration, new_x, new_rms);
            }
         n_minima++;
         }
      }
   find_trial_orbit( orbit, obs, n_obs, r1, best_found);
   *angle_param = best_found;
   return( rval);
}

static double haltonize( unsigned idx, const unsigned radix)
{
   unsigned divisor = 1, rval = 0;

   while( idx)
      {
      rval = (rval * radix) + idx % radix;
      idx /= radix;
      divisor *= radix;
      }
   return( (double)rval / (double)divisor);
}

/* Any sigmas/covariances computed for an orbit really should only
apply for a particular epoch,  observations included,  central object,
and perturbers.  That is to say,  once you've added a new perturber
or switched from Jovicentric to heliocentric elements or shifted to
a new epoch,  the sigmas/covariances you've computed should be
ignored.

   The way I've implemented this is to say that,  when you've just
computed sigmas/covariances,  a hash is computed from the epoch,
included observations,  etc.  When elements are displayed (see
'elem_out.cpp'),  we check to see if we get the same hash.  If we
do,  then we can assume the sigmas/covariances are valid.  If not,
we avoid the blunder of displaying erroneous data.

   Part of this requires the following very primitive and slow,
but dead simple and good at avoiding collisions,  hashing algorithm:
XOR the hash value with each byte in turn,  then multiply by a big
prime.  Idea borrowed from the Linux ext3 file system.  */

static int simple_hash( int ival, const char *ibuff, int nbytes)
{
   const long big_prime = 1234567891;

   while( nbytes--)
      ival = (ival ^ (int)*ibuff++) * big_prime;
   return( ival);
}

int compute_available_sigmas_hash( const OBSERVE FAR *obs, const int n_obs,
         const double epoch, const unsigned perturber_mask, const int central_obj)
{
   int temp_array[10], rval, i;

   rval = simple_hash( available_sigmas + 1, (const char *)&epoch, sizeof( double));
   temp_array[0] = (perturber_mask & AUTOMATIC_PERTURBERS ?
               perturbers_automatically_found : perturber_mask);
   temp_array[1] = central_obj;
   if( debug_level)
      debug_printf( "Hashing: %d %x %d %f\n", n_obs, temp_array[0], central_obj, epoch);
   rval = simple_hash( rval, (const char *)temp_array, 2 * sizeof( int));
   for( i = 0; i < n_obs; i++, obs++)
      if( obs->is_included)
         {
         double temp_darray[10];

         temp_darray[0] = obs->jd;
         temp_darray[1] = obs->posn_sigma_1;
         temp_darray[2] = obs->posn_sigma_2;
         temp_darray[3] = obs->posn_sigma_theta;
         temp_darray[4] = obs->mag_sigma;
         rval = simple_hash( rval, (const char *)temp_darray, 5 * sizeof( double));
         }
   if( debug_level)
      debug_printf( "Hash result: %8x\n", (unsigned)rval);
   return( rval);
}

/* When running statistical ranging (SR) orbits,  you would like to consider
only those beyond some minimum distance.  In theory,  it's possible that you
saw the object in one place a few meters away from the telescope,  as it was
speeding away almost straight up;  then saw it again a week later,  again a
few meters in front of the telescope,  as it was coming down.

   But considering such SR variant orbits is a bad idea in several respects.
First,  the only way this can happen is with an object launched off the
earth and then coming back;  for a "real" orbit,  a basic constraint is
that the object has to actually be in orbit around the earth (perigee
greater than the earth's radius,  or the object's trajectory before and
during the observed arc can't pass through the earth).

   Ideally,  we'd simply discard any orbit that would mean the object came
out of the earth.  But at least for the moment,  I'm constraining the
minimum distance for the SR orbit instead to be at least (a) the distance at
which it would make a third of an orbit around the earth over the time in
question,  _or_ (b) the distance it would cover at 1 km/s, whichever is
greater.  This is really quite ad hoc,  and probably underestimates the
"real" distance most of the time.  (The exception would probably be certain
high-orbiting geosats and objects at the Earth-Sun L2 point.  In those
cases,  SR may fail with the following. But then other techniques will
probably succeed,  so I'm not too worried about this.)

   In the following,  I take advantage of the fact that geostationary
satellites orbit at about 42164 km,  and take one day to orbit us.  */

inline double minimum_sr_distance( const double time_span)
{
   const double geosat_dist = 42164. / AU_IN_KM;
   const double geo_orbit_dist = geosat_dist * pow( time_span * 3., 2. / 3.);
   const double lowest_speed = 1.;    /* km/s */
   const double slow_limit = time_span * lowest_speed * seconds_per_day / AU_IN_KM;

   return( slow_limit < geo_orbit_dist ? geo_orbit_dist : slow_limit);
}

double sr_min_r = 0., sr_max_r = 0.;

int n_sr_ranges;
double sr_roots[10];

double find_sr_dist( const double fraction)
{
   double total_dist = 0., dist = 0., target_dist;
   int i;

   assert( fraction >= 0.);
   assert( fraction <= 1.);
   for( i = 0; i < n_sr_ranges * 2; i += 2)
      total_dist += sr_roots[i + 1] - sr_roots[i];
   target_dist = total_dist * fraction;
   for( i = 0; i < n_sr_ranges * 2; i += 2)
      {
      const double delta = sr_roots[i + 1] - sr_roots[i];
      const double new_dist = dist + delta;

      if( new_dist < target_dist)
         dist = new_dist;
      else
         return( sr_roots[i] + target_dist - dist);
      }
   assert( 1);
   return( sr_roots[n_sr_ranges * 2 - 1]);         /* maxed out */
}

int find_nth_sr_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
                            const int orbit_number)
{
   int rval = 0;

   uncertainty_parameter = 99.;
   obs += drop_excluded_obs( obs, &n_obs);
   if( n_obs < 2)       /* need at least two valid obs */
      {
      debug_printf( "SR fail: n_obs = %d\n", n_obs);
      rval = -1;
      }
   else
      {
      double rand1 = haltonize( (unsigned)orbit_number + 1, 2);
      const double rand2 = haltonize( (unsigned)orbit_number + 1, 3);
      double dist = 0.;
      int i;

      if( !orbit_number)
         {
         if( sr_max_r)        /* user override of SR ranges */
            {
            sr_roots[0] = sr_min_r;
            sr_roots[1] = sr_max_r;
            n_sr_ranges = 1;
            }
         else        /* "normal" determination of possible ranges */
            {
            const double time_range = obs[n_obs - 1].jd - obs[0].jd;
            const double min_range = minimum_sr_distance( time_range);

            n_sr_ranges = find_sr_ranges( sr_roots,
                        obs[n_obs - 1].obs_posn, obs[n_obs - 1].vect,
                        obs[    0    ].obs_posn, obs[    0    ].vect,
                        SOLAR_GM, time_range);
            for( i = 0; i < n_sr_ranges * 2; i++)
               if( sr_roots[i] < min_range)
                  sr_roots[i] = min_range;
            }
         if( debug_level > 1)
            {
            debug_printf( "%d ranges\n", n_sr_ranges);
            for( i = 0; i < n_sr_ranges * 2; i++)
               debug_printf( "Root %d: %f\n", i, sr_roots[i]);
            }
         }
            /* We give 'rand1' a slight bias toward lower values. */
            /* It will still be in the range 0 <= rand1 < 1.      */
      rand1 = rand1 * (1. + rand1) / 2.;
      dist = find_sr_dist( rand1);
      fail_on_hitting_planet = true;
      rval = find_trial_orbit( orbit, obs, n_obs, dist, 2. * rand2 - 1.);
      fail_on_hitting_planet = false;
      }
   return( rval);
}

static int sr_orbit_compare( const void *a, const void *b)
{
   const double *ta = (const double *)a;
   const double *tb = (const double *)b;

   return( (ta[6] > tb[6]) ? 1 : -1);
}

static bool writing_sr_elems = true;

int get_sr_orbits( double *orbits, OBSERVE FAR *obs,
               const unsigned n_obs, const unsigned starting_orbit,
               const unsigned max_orbits, const double max_time,
               const double noise_in_sigmas)
{
   const clock_t end_clock =
             clock( ) + (clock_t)( max_time * (double)CLOCKS_PER_SEC);
   unsigned i, rval = 0;
   double *tptr = orbits;

   INTENTIONALLY_UNUSED_PARAMETER( noise_in_sigmas);
// perturbers = AUTOMATIC_PERTURBERS;
   for( i = 0; i < max_orbits && clock( ) < end_clock; i++)
      {
//    double *stored_ra_decs_mags_times =
//                 add_gaussian_noise_to_obs( n_obs, obs, noise_in_sigmas);

      if( !find_nth_sr_orbit( tptr, obs, n_obs, i + starting_orbit)
                   && (n_obs == 2 || !adjust_herget_results( obs, n_obs, tptr)))
         {
//       tptr[6] = compute_rms( obs, n_obs);
         tptr[6] = evaluate_initial_orbit( obs, n_obs, tptr);
         rval++;
         tptr += 7;
         }
//    restore_ra_decs_mags_times( n_obs, obs, stored_ra_decs_mags_times);
//    free( stored_ra_decs_mags_times);
      }
   qsort( orbits, rval, 7 * sizeof( double), sr_orbit_compare);
   if( writing_sr_elems)
      for( i = 0; i < rval; i++)
         {
         extern const char *elements_filename;
         const char *tname = elements_filename;
         extern int append_elements_to_element_file;
         int curr_append = append_elements_to_element_file;

         elements_filename = "sr_elems.txt";
         append_elements_to_element_file = (i ? 1 : 0);
         set_locs( orbits, obs[0].jd, obs, n_obs);
         write_out_elements_to_file( orbits, obs[0].jd,
                  find_epoch_shown( obs, n_obs),
                  obs, n_obs, "", 5,
                  1, ELEM_OUT_NO_COMMENT_DATA | ELEM_OUT_PRECISE_MEAN_RESIDS);
         orbits += 7;
         append_elements_to_element_file = curr_append;
         elements_filename = tname;
         }
   return( rval);
}

static inline void compute_sr_sigmas( const double *sr_orbits,
               const unsigned n_orbits, const double epoch,
               const double epoch_shown)
{
   unsigned i;
   double monte_data[MONTE_DATA_SIZE];
   double sigmas[MONTE_N_ENTRIES];
   FILE *monte_file;
   char filename[100];
   const int planet_orbiting = 0;      /* heliocentric only,  at least for now */
   ELEMENTS elem0;

   elem0.major_axis = elem0.ecc = 0.;     /* just to avoid uninitialized  */
   for( i = 0; i < n_orbits; i++)         /* variable warnings            */
      {
      double orbit[6];
      ELEMENTS elem;

      memset( &elem, 0, sizeof( ELEMENTS));
      elem.gm = SOLAR_GM;
      memcpy( orbit, sr_orbits + 6 * i, 6 * sizeof( double));
      integrate_orbit( orbit, epoch, epoch_shown);
      calc_classical_elements( &elem, orbit, epoch_shown, 1);
      add_monte_orbit( monte_data, &elem, i);
      if( !i)
         elem0 = elem;
      }
   monte_file = fopen_ext( get_file_name( filename, "monte.txt"), "tfcwb");

   fprintf( monte_file, "Computed from %u SR orbits\n", n_orbits);
   compute_monte_sigmas( sigmas, monte_data, n_orbits);
   uncertainty_parameter = dump_monte_data_to_file( monte_file, sigmas,
                 elem0.major_axis, elem0.ecc, planet_orbiting);
   available_sigmas = SR_SIGMAS_AVAILABLE;
   fclose( monte_file);
}

static OBSERVE *get_real_arc( OBSERVE *obs, int *n_obs,
                              int *n_real_obs)
{
   int idx1, idx2;

   *n_real_obs = get_idx1_and_idx2( *n_obs, obs, &idx1, &idx2);
   *n_obs = idx2 - idx1 + 1;
   return( obs + idx1);
}

#define RAD2SEC (180. * 3600. / PI)

#ifdef NO_LONGER_USED

int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit)
{
   int i, n_found, rval = 0;

   obs = get_real_arc( obs, &n_obs, &n_found);
   if( n_found < 2)   /* must have at least two obs */
      rval = -2;
   else if( is_unreasonable_orbit( orbit))
      rval = -3;
   else
      {
      double determ;
      double avg_x = 0., avg_xt = 0., avg_t = 0.;
      double avg_y = 0., avg_yt = 0., avg_t2 = 0.;

      for( i = 0; i < n_obs; i++)
         if( obs[i].is_included)
            {
            const double dx = obs[i].computed_ra - obs[i].ra;
            const double dy = obs[i].computed_dec - obs[i].dec;
            const double dt = obs[i].jd - obs[0].jd;

            avg_xt += dx * dt;
            avg_yt += dy * dt;
            avg_x += dx;
            avg_y += dy;
            avg_t += dt;
            avg_t2 += dt * dt;
            }

      avg_xt /= (double)n_found;
      avg_yt /= (double)n_found;
      avg_x /= (double)n_found;
      avg_y /= (double)n_found;
      avg_t /= (double)n_found;
      avg_t2 /= (double)n_found;
      determ = avg_t2 - avg_t * avg_t;
      if( determ)
         {
         const double ax = (avg_xt - avg_x * avg_t) / determ;
         const double bx = (avg_t2 * avg_x - avg_t * avg_xt) / determ;
         const double ay = (avg_yt - avg_y * avg_t) / determ;
         const double by = (avg_t2 * avg_y - avg_t * avg_yt) / determ;
         OBSERVE start = obs[0];
         OBSERVE end = obs[n_obs - 1];

         start.ra  = start.computed_ra - bx;
         end.ra    = end.computed_ra - (bx + (end.jd - start.jd) * ax);
         start.dec = start.computed_dec - by;
         end.dec   = end.computed_dec - (by + (end.jd - start.jd) * ay);
         set_obs_vect( &start);
         set_obs_vect( &end);
         set_distance( &start, start.r);
         set_distance( &end, end.r);
         uncertainty_parameter = 99.;
         if( find_transfer_orbit( orbit, &start, &end, 1))
            rval = -8;
         else
            rval = set_locs( orbit, start.jd, obs, n_obs);
         }
      else
         rval = -1;
      }
   if( rval && rval != -2)
      debug_printf( "Adjust fail %d\n", rval);
   return( rval);
}
#endif      /* #ifdef NO_LONGER_USED */

int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit)
{
   int n_found, rval = 0;

   obs = get_real_arc( obs, &n_obs, &n_found);
   if( n_found < 2)   /* must have at least two obs */
      rval = -2;
   else if( is_unreasonable_orbit( orbit))
      rval = -3;
   else
      {
      rval = extended_orbit_fit( orbit, obs, n_obs,
                     FIT_FIXED_DISTANCES, obs->jd);
      if( !rval)
         rval = set_locs( orbit, obs->jd, obs, n_obs);
      }
   return( rval);
}

int herget_method( OBSERVE FAR *obs, int n_obs, double r1, double r2,
         double *orbit, double *d_r1, double *d_r2, const char *limited_orbit)
{
   double *xresid, *yresid, *dx_dr1, *dx_dr2, *dy_dr1, *dy_dr2;
   double delta, a = 0., b = 0., c = 0., e = 0., f = 0., determ;
   int i, n_real_obs;
   const int using_pseudo_vaisala = (r1 < 0.);
   double orbit2[6];
   double orbit_offset[6], *constraint = NULL;
   int planet_orbiting = 0, n_constraints = 0;
   char tstr[80];
   OBSERVE temp_obs1, temp_obs2;

   if( limited_orbit && !*limited_orbit)
      limited_orbit = NULL;

   obs = get_real_arc( obs, &n_obs, &n_real_obs);
   if( n_real_obs < 2)        /* should never happen */
      {
      debug_printf( "??? n_obs %d, n_real_obs = %d, %s\n",
                     n_obs, n_real_obs, obs->packed_id);
      return( -1);
      }

   temp_obs1 = obs[0];
            /* Look "ahead" up to 100 days: */
   for( i = n_obs - 1; i > 0 && obs[i].jd - obs[0].jd > 100.; i--)
      ;
   temp_obs2 = obs[i];
   uncertainty_parameter = 99.;
   if( using_pseudo_vaisala)
      {
      r2 = find_r_given_solar_r( &temp_obs2, -r1);
      r1 = find_r_given_solar_r( &temp_obs1, -r1);

      if( debug_level > 7)
         debug_printf( "r1 = %f; r2 = %f: %s\n",
                  r1, r2, obs->packed_id);
      if( r1 < 0. || r2 < 0.)
         return( 1);
      if( d_r1)
         *d_r1 = r1;
      if( d_r2)
         *d_r2 = r2;
      }
   set_distance( &temp_obs1, r1);
   set_distance( &temp_obs2, r2);
   runtime_message = tstr;
   if( using_pseudo_vaisala)
      sprintf( tstr, "Vaisala %f\n", obs->solar_r);
   else
      strcpy( tstr, "H/xfer orbit (1)");
               /* Compute the trial orbit in the local orbit2 array.  That */
               /* way,  if we find it's completely stupid,  we've not      */
               /* done anything to the plain old 'orbit' vector,  which    */
               /* may already contain a valid state vector:                */
   if( find_transfer_orbit( orbit2, &temp_obs1, &temp_obs2, 0)
                   || is_unreasonable_orbit( orbit2))
      {
      runtime_message = NULL;
      return( -2);
      }
               /* But now that we know it's a good result,  let's copy:     */
   memcpy( orbit, orbit2, 6 * sizeof( double));
   available_sigmas = NO_SIGMAS_AVAILABLE;
   strcpy( tstr, using_pseudo_vaisala ? "Vaisala set_locs" : "H/set_locs (1)");
   if( set_locs( orbit, temp_obs1.jd, obs, n_obs))
      {
      runtime_message = NULL;
      return( -3);
      }
   runtime_message = NULL;
   if( !d_r1 || !d_r2 || using_pseudo_vaisala)
      return( 0);

   *d_r1 = *d_r2 = 0.;
   if( n_real_obs == 2)
      return( 0);       /* good as done... */

   xresid = (double *)calloc( 6 * (n_obs + MAX_CONSTRAINTS), sizeof( double));
   if( !xresid)
      return( -4);
   yresid = xresid + n_obs + MAX_CONSTRAINTS;
   dx_dr1 = yresid + n_obs + MAX_CONSTRAINTS;
   dx_dr2 = dx_dr1 + n_obs + MAX_CONSTRAINTS;
   dy_dr1 = dx_dr2 + n_obs + MAX_CONSTRAINTS;
   dy_dr2 = dy_dr1 + n_obs + MAX_CONSTRAINTS;

   if( limited_orbit)
      {
      planet_orbiting = find_best_fit_planet( temp_obs1.jd, orbit, orbit2);
      for( i = 0; i < 6; i++)
         orbit_offset[i] = orbit[i] - orbit2[i];
      constraint = xresid + n_obs;
      n_constraints = evaluate_limited_orbit( orbit2, planet_orbiting,
                                       temp_obs1.jd, limited_orbit, constraint);
      }

   for( i = 0; i < n_obs; i++)
      get_residual_data( obs + i, xresid + i, yresid + i);

   delta = r1 / 10000.;
   if( delta > .1) delta = .1;
   set_distance( &temp_obs1, r1 + delta);
   runtime_message = tstr;
   strcpy( tstr, "H/xfer orbit (2)");
   memcpy( orbit2, orbit, 6 * sizeof( double));
   if( find_transfer_orbit( orbit2, &temp_obs1, &temp_obs2, 1))
      {
      runtime_message = NULL;
      free( xresid);
      return( -5);
      }
   strcpy( tstr, "H/set_locs (2)");
   if( set_locs( orbit2, temp_obs1.jd, obs, n_obs))
      {
      runtime_message = NULL;
      free( xresid);
      return( -6);
      }

   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         double dx, dy;

         get_residual_data( obs + i, &dx, &dy);
         dx_dr1[i] = (dx - xresid[i]) / delta;
         dy_dr1[i] = (dy - yresid[i]) / delta;
         }

   if( limited_orbit)
      {
      double constraint2[MAX_CONSTRAINTS];

      for( i = 0; i < 6; i++)
         orbit2[i] -= orbit_offset[i];
      evaluate_limited_orbit( orbit2, planet_orbiting, temp_obs1.jd,
                                limited_orbit, constraint2);
      for( i = 0; i < n_constraints; i++)
         {
         dx_dr1[n_obs + i] = (constraint2[i] - constraint[i]) / delta;
         dy_dr1[n_obs + i] = 0.;
         }
      }

   delta = r2 / 10000.;
   if( delta > .1) delta = .1;
   set_distance( &temp_obs1, r1);
   set_distance( &temp_obs2, r2 + delta);
   strcpy( tstr, "H/xfer orbit (3)");
   memcpy( orbit2, orbit, 6 * sizeof( double));
   if( find_transfer_orbit( orbit2, &temp_obs1, &temp_obs2, 1))
      {
      runtime_message = NULL;
      free( xresid);
      return( -7);
      }
   strcpy( tstr, "H/set_locs (3)");
   if( set_locs( orbit2, temp_obs1.jd, obs, n_obs))
      {
      runtime_message = NULL;
      free( xresid);
      return( -8);
      }
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         double dx, dy;

         get_residual_data( obs + i, &dx, &dy);
         dx_dr2[i] = (dx - xresid[i]) / delta;
         dy_dr2[i] = (dy - yresid[i]) / delta;
         }

   if( limited_orbit)
      {
      double constraint2[MAX_CONSTRAINTS];

      for( i = 0; i < 6; i++)
         orbit2[i] -= orbit_offset[i];
      evaluate_limited_orbit( orbit2, planet_orbiting, temp_obs1.jd,
                                limited_orbit, constraint2);
      for( i = 0; i < n_constraints; i++)
         {
         dx_dr2[n_obs + i] = (constraint2[i] - constraint[i]) / delta;
         dy_dr2[n_obs + i] = 0.;
         }
      }

                    /* OK,  we now have all values & derivatives needed... */
   for( i = 0; i < n_obs + n_constraints; i++)
      if( obs[i].is_included || i >= n_obs)
         {
         c += xresid[i] * dx_dr1[i] + yresid[i] * dy_dr1[i];
         a += dx_dr1[i] * dx_dr1[i] + dy_dr1[i] * dy_dr1[i];
         b += dx_dr1[i] * dx_dr2[i] + dy_dr1[i] * dy_dr2[i];

         f += xresid[i] * dx_dr2[i] + yresid[i] * dy_dr2[i];
         /*  d = b;  */
         e += dx_dr2[i] * dx_dr2[i] + dy_dr2[i] * dy_dr2[i];
         }

   free( xresid);
   determ = a * e - b * b;
   runtime_message = NULL;
   if( !determ)
      return( -9);
   else
      {
      *d_r1 = -(e * c - b * f) / determ;
      *d_r2 = -(a * f - c * b) / determ;
               /* Just to keep r1 and r2 from becoming negative or zipping */
               /* off to infinity,  insist on a maximum change of 1/3:     */
      if( *d_r1 > r1 / 3.) *d_r1 = r1 / 3.;
      if( *d_r1 <-r1 / 3.) *d_r1 = -r1 / 3.;
      if( *d_r2 > r2 / 3.) *d_r2 = r2 / 3.;
      if( *d_r2 <-r2 / 3.) *d_r2 = -r2 / 3.;
      }
   return( 0);
}

#ifdef OBSOLETE

   /* The following code _used_ to be useful for determining
      parabolic orbits.  But now,  those are done just by using
      an e=1 constraint,  same as with any other type of orbit. */

static int setup_parabolic( const double *iparams, double *orbit)
{
   double r = 0., escape_vel;
   int i;
   const double sqrt_2 = 1.414213562373095048801688724209698078569671875376948;

   for( i = 0; i < 3; i++)
      {
      orbit[i] = iparams[i];
      r += iparams[i] * iparams[i];
      }
   r = sqrt( r);
   escape_vel = sqrt_2 * GAUSS_K / sqrt( r);
   orbit[3] = escape_vel * cos( iparams[3]) * cos( iparams[4]);
   orbit[4] = escape_vel * sin( iparams[3]) * cos( iparams[4]);
   orbit[5] = escape_vel *                    sin( iparams[4]);
   return( 0);
}

void improve_parabolic( OBSERVE FAR *obs, int n_obs, double *orbit,
                                                              double epoch)
{
   void *lsquare = lsquare_init( 5);
   double *xresids = (double *)calloc( 2 * n_obs + 10 * n_obs, sizeof( double));
   double *yresids = xresids + n_obs;
   double *slopes = yresids + n_obs;
   double params2[5], params[5], differences[5];
   const double delta_val = .0001;
   double v2 = 0.;
   int i, j;

   available_sigmas = NO_SIGMAS_AVAILABLE;
   uncertainty_parameter = 99.;
   memcpy( params, orbit, 3 * sizeof( double));
   params[3] = atan2( orbit[4], orbit[3]);
   for( i = 3; i < 6; i++)
      v2 += orbit[i] * orbit[i];
   params[4] = asin( orbit[5] / sqrt( v2));
   setup_parabolic( params, orbit);

   if( set_locs( orbit, epoch, obs, n_obs))
      {
      free( xresids);
      return;
      }

   for( i = 0; i < n_obs; i++)
      {
      xresids[i] = obs[i].computed_ra - obs[i].ra;
      yresids[i] = obs[i].computed_dec - obs[i].dec;
      }

   for( i = 0; i < 5; i++)
      {
      memcpy( params2, params, 5 * sizeof( double));
      params2[i] = params[i] - delta_val;
      setup_parabolic( params2, orbit);
      set_locs( orbit, epoch, obs, n_obs);
      for( j = 0; j < n_obs; j++)
         {
         slopes[i + j * 10] = obs[j].computed_ra;
         slopes[i + j * 10 + 5] = obs[j].computed_dec;
         }

      params2[i] = params[i] + delta_val;
      setup_parabolic( params2, orbit);
      set_locs( orbit, epoch, obs, n_obs);
      for( j = 0; j < n_obs; j++)
         {
         slopes[i + j * 10] -= obs[j].computed_ra;
         slopes[i + j * 10 + 5] -= obs[j].computed_dec;
         slopes[i + j * 10] /= 2. * delta_val;
         slopes[i + j * 10 + 5] /= 2. * delta_val;
         }
      }

   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         lsquare_add_observation( lsquare, xresids[i], 1., slopes + i * 10);
         lsquare_add_observation( lsquare, yresids[i], 1., slopes + i * 10 + 5);
         }

   free( xresids);
   lsquare_solve( lsquare, differences);
   lsquare_free( lsquare);

   for( i = 0; i < 5; i++)
      params[i] += differences[i];
   setup_parabolic( params, orbit);
   set_locs( orbit, epoch, obs, n_obs);
}

#endif        // #ifdef OBSOLETE

void compute_variant_orbit( double *variant, const double *ref_orbit,
                     const double n_sigmas)
{
   unsigned i;
   extern double **eigenvects;

   assert( eigenvects);
   if( eigenvects)
      for( i = 0; i < 6; i++)
         *variant++ = *ref_orbit++ + (n_sigmas * eigenvects[0][i]);
}

double get_planet_mass( const int planet_idx)
{
   double rval;
   const extern double planet_mass[];

   if( planet_idx == -1)         /* barycentric:  return entire SS mass */
      rval = 1.001341838645975;     /* in units of sun=1 */
   else if( planet_idx == 100)      /* Ceres: special 63 km^3/s^2 value */
      rval = 4.747105847157599504245142421641e-10;
   else
      rval = planet_mass[planet_idx];
   return( rval * SOLAR_GM);
}

double evaluate_for_simplex_method( const OBSERVE FAR *obs,
                    const int n_obs, const double *orbit,
                    const int planet_orbiting,
                    const char *limited_orbit)
{
   double rval = evaluate_initial_orbit( obs, n_obs, orbit);

   if( limited_orbit && *limited_orbit)
      {
      const double scale = 100.;
      double constraints[MAX_CONSTRAINTS];
      const double epoch = obs->jd;
      int i;
      const int n_constraints = evaluate_limited_orbit( orbit,
                     planet_orbiting, epoch, limited_orbit, constraints);

      for( i = n_constraints - 1; i >= 0; i--)
         rval += constraints[i] * constraints[i] * scale;
      }
   return( rval);
}

static int parse_constraint( const char *tptr)
{
   int rval = 0;

   if( *tptr && tptr[1] == '=')
      rval = *tptr;
   else if( *tptr == 'A' && tptr[1] >= '1' && tptr[1] <= '3'
                        && tptr[2] == '=')
      rval = tptr[1];
   else if( !memcmp( tptr, "Tp=", 3))
      rval = 'T';
   return( rval);
}

void get_periapsis_loc( double *ecliptic_lon, double *ecliptic_lat,
             const ELEMENTS *elem);             /* elem_out.cpp */

static int evaluate_limited_orbit( const double *orbit,
                    const int planet_orbiting, const double epoch,
                    const char *limited_orbit, double *constraints)
{
   int rval = 0;

   if( limited_orbit)
      {
      ELEMENTS elem;
      int i, variable;

      elem.gm = get_planet_mass( planet_orbiting);
      calc_classical_elements( &elem, orbit, epoch, 1);
      while( (variable = parse_constraint( limited_orbit)) != 0)
         {
         const int offset = (limited_orbit[1] == '=' ? 2 : 3);
         double value = atof( limited_orbit + offset);
         int angular_constraint = 0, tbuff_loc = 0;
         char tbuff[80];

         limited_orbit += offset;
         while( *limited_orbit && *limited_orbit != ',')
            tbuff[tbuff_loc++] = *limited_orbit++;
         tbuff[tbuff_loc] = '\0';
         if( limited_orbit[-1] == 'k' ||
                  (limited_orbit[-2] == 'k' && limited_orbit[-1] == 'm'))
            value /= AU_IN_KM;
         switch( variable)
            {
            case 'q':
               constraints[rval] = (elem.major_axis * (1. - elem.ecc) / value) - 1.;
               rval++;
               break;
            case 'Q':
               constraints[rval] = (elem.major_axis * (1. + elem.ecc) / value) - 1.;
               rval++;
               break;
            case 'e':
               assert( fabs( elem.ecc) < 1e+9);
               if( value)
                  constraints[rval++] = (elem.ecc - value) * 10000.;
               else        /* handle e=0 (circular) orbits separately: */
                  {
                  constraints[rval++] = orbit[0] * orbit[3] +
                              orbit[1] * orbit[4] + orbit[2] * orbit[5];
                  constraints[rval++] = vector3_length( orbit)
                               - elem.major_axis;
                  constraints[rval - 2] *= 10000.;
                  constraints[rval - 1] *= 10000.;
                  }
               break;
            case 'i':
               constraints[rval++] = elem.incl * 180. / PI - value;
               angular_constraint = 1;
               break;
            case 'O':
               constraints[rval++] = elem.asc_node * 180. / PI - value;
               angular_constraint = 1;
               break;
            case 'o':
               constraints[rval++] = elem.arg_per * 180. / PI - value;
               angular_constraint = 1;
               break;
            case 'P':         /* convert to a major axis */
               if( limited_orbit[-1] == 'd')    /* convert from days to yrs */
                  value /= 365.25;
               if( limited_orbit[-1] == 'h')    /* convert from hrs to yrs */
                  value /= 365.25 * 24.;
               if( limited_orbit[-1] == 'm')    /* convert from mins to yrs */
                  value /= 365.25 * 1440.;
               value = pow( value * sqrt( get_planet_mass( planet_orbiting) / SOLAR_GM),
                                                     2. / 3.);
                           /* fall-thru   */
            case 'a':
               constraints[rval++] = (elem.major_axis / value) - 1.;
               break;
            case 'n':
               value = 360. / value;         /* now value = period in days */
               value /= 365.25;              /* now value = period in years */
               value = pow( value * sqrt( get_planet_mass( planet_orbiting) / SOLAR_GM),
                                                     2. / 3.);
               constraints[rval++] = 1. / elem.major_axis - 1. / value;
               break;
            case 'A':            /* area/mass ratio */
               if( n_extra_params >= 1)
                  constraints[rval++] =
                      10000. * (solar_pressure[0] * SOLAR_GM / SRP1AU - value);
               if( n_extra_params > 1)
                  constraints[rval - 1] *= 100000.;
               break;
            case '1': case '2': case '3':
               constraints[rval++] = 1e+10 * (solar_pressure[variable - '1'] - value);
               break;
            case 'K':
               {
               extern double comet_magnitude_slope_param;

               comet_magnitude_slope_param = value;
               }
               break;
            case 'G':
               {
               extern double asteroid_magnitude_slope_param;

               asteroid_magnitude_slope_param = value;
               }
               break;
            case 'T':
               {
               const double tp = get_time_from_string( 0., tbuff, 0, NULL);

               if( tp)
                  constraints[rval++] = (tp - elem.perih_time) * 1e+5;
               }
               break;
            case 'b':      /* constraints in ecliptic lat/lon are quite */
            case 'l':      /* useful with sungrazers such as Kreutz objects */
               {
               double ecliptic_lon, ecliptic_lat, ang;

               get_periapsis_loc( &ecliptic_lon, &ecliptic_lat, &elem);
               ang = ((variable == 'b') ? ecliptic_lon : ecliptic_lat);
               constraints[rval++] = (ang * 180. / PI - value) * .001;
               }
               break;
            }
         if( angular_constraint)
            {
            while( constraints[rval - 1] > 180.)
               constraints[rval - 1] -= 360.;
            while( constraints[rval - 1] < -180.)
               constraints[rval - 1] += 360.;
            }
         if( *limited_orbit == ',')
            limited_orbit++;
         }
      for( i = 0; i < rval; i++)
         constraints[i] *= 1e+5;
//       constraints[i] *= 1e+10;
      }
   return( rval);
}

#define MAX_N_PARAMS 9

double **eigenvects;
void jacobi_eigenvalues( double *a, const int size, double *eigenvals,
                        double *eigenvects);       /* eigen.cpp */

void **calloc_double_dimension_array( const size_t x, const size_t y,
                                    const size_t obj_size)
{
   void **rval = (void **)calloc( x * sizeof( void *) + x * y * obj_size, 1);
   size_t i;

   if( rval)
      {
      rval[0] = (void *)(rval + x);
      for( i = 1; i < x; i++)
         rval[i] = (void *)( (char *)rval[0] + i * y * obj_size);
      }
   return( rval);
}

static double dotted_dist( OBSERVE FAR *obs)
{
   return( dot_product( obs->vect, obs->obj_posn) - dot_product( obs->vect, obs->obs_posn));
}

double n_nearby_obs( const OBSERVE FAR *obs, const unsigned n_obs,
          const unsigned idx, const double time_span)
{
   unsigned i = idx;
   double rval = 1., dt = 0.;
   const double span_limit = 4.;

   while( i && obs[idx].jd - obs[i - 1].jd < time_span * span_limit)
      i--;
   while( i < n_obs && dt < span_limit)
      {
      const double dt = (obs[i].jd - obs[idx].jd) / time_span;

      if( obs[i].is_included && i != idx
              && !strcmp( obs[i].mpc_code, obs[idx].mpc_code))
         rval += exp( -dt * dt / 2.);
      i++;
      }
   return( rval);
}

/* "Correcting" for over-observing is a rather complicated issue.  From 2012
Nov 8 to 2017 Sep 7,  I used a method that I got after discussion with Marco
Micheli and Alan Harris:  if an object is observed more than 'k' times in a
night,  its observations should be re-weighted by sqrt( k/Nobs), where Nobs
is the number of observations over that night.  The idea was that after
getting 'k' observations, you've beaten down the random observational
errors,  and you aren't going to do anything about the systematic errors.
This scheme was later used in FCCT15 (with k=1) and VFCC17 (with k=4),  and
is shown below in the PREVIOUS_METHOD_SEE_ABOVE_COMMENTS section.

   Since then,  I've revisited the question of what Nobs should be (number
of observations that night?  Within four hours?  etc.) and have put the
reweighting factor on a less ad hoc,  more solid mathematical foundation.
See https://www.projectpluto.com/errors.htm for details.  */

double overobserving_time_span = 0.;
unsigned overobserving_ceiling = 4;

#define IDX_ASTEROIDS 20

static double reweight_for_overobserving( const OBSERVE FAR *obs,
            const unsigned n_obs, const unsigned idx)
{
   double rval;
   const double n_nearby = n_nearby_obs( obs, n_obs, idx,
                                   overobserving_time_span);

   if( !strcmp( obs->mpc_code, "258"))    /* Gaia comes already corrected */
      return( 1.);                        /* for over-observing           */
#ifdef PREVIOUS_METHOD_SEE_ABOVE_COMMENTS
   if( n_nearby > overobserving_ceiling)
      rval = sqrt( (double)overobserving_ceiling / n_nearby);
   else
      rval = 1.;
#else                  /* newer,  less ad hoc reweighting scheme */
   const double Nmax = (double)overobserving_ceiling;

   rval = sqrt( Nmax / (n_nearby + Nmax - 1.));
#endif
   return( rval);
}

void get_relative_vector( const double jd, const double *ivect,
                  double *relative_vect, const int planet_orbiting)
{
   double planet_state_vect[6];
   int i;

   compute_observer_loc( jd, planet_orbiting, 0., 0., 0., planet_state_vect);
   compute_observer_vel( jd, planet_orbiting, 0., 0., 0., planet_state_vect + 3);
   for( i = 0; i < 6; i++)
      relative_vect[i] = ivect[i] - planet_state_vect[i];
}

/* Computing numerical derivatives in 'full_improvement' has some tricky
points to it.  One is the classic issue:  if you're computing the
derivative by either of the "standard" methods,

df(x)/dx = (f(x + h) - f(x)) / h + O(h)
df(x)/dx = (f(x + h) - f(x - h)) / 2h + O(h^2)

   how big should the 'tweak' h be?  Too small,  and roundoff errors will
cause trouble.  Too big,  and that O(h) or O(h^2) error will cause trouble.

   'full_improvement' attempts to address this issue dynamically.  It has
to do six tweaks (sometimes seven or nine,  if you're solving for additional
parameters such as non-gravitational params or an asteroid mass).  It looks
at each tweak to make sure it had an effect somewhere between 'min_change' to
'max_change' (see full_improvement() for those variables.)  If it fell outside
that range,  the tweak is re-scaled and we try again.  */

bool use_symmetric_derivatives = false;
int forced_central_body = 0;
double probability_of_blunder = 0.;
int use_blunder_method = 0;

     /* If use_blunder_method == 0: filtering is done in the "traditional" */
     /* method,  excluding obs beyond a certain number of sigmas or arcsec */
     /* If use_blunder_method == 1: filtering is to be done using blunder  */
     /* re-weighting,  with the following function                         */
     /* If use_blunder_method == 2: actually filtering within the          */
     /* full_improvement( ) function.                                      */
     /* In blunder re-weighting,  we always allow _some_ weight of 1e-90.  */
     /* This basically just guards against the situation in which all      */
     /* observations are so far off as to be considered to have zero       */
     /* probability of being "real".  In that case,  it would be as if no  */
     /* observations were made,  and full_improvement( ) would fail.       */

static double reweight_for_blunders( const double resid2, const double weight)
{
   const double resid_in_arcseconds = sqrt( resid2) * (180. / PI) * 3600.;
   const double prob_of_obs = exp( -resid_in_arcseconds * weight);
   const double minimum_rval = 1e-90;
   const double rval = weight * prob_of_obs / (probability_of_blunder + prob_of_obs);

   return( (minimum_rval > rval) ? minimum_rval : rval);
}

/* Each "observation" is really zero,  one,  or two residuals.  Most optical
observations will actually be an RA observation plus a declination
observation,  treated below as an "along-track" and "cross-track" observation
to simplify inclusion of timing errors.  Radar observations will actually be
a range observation,  a Doppler (range-rate) observation,  or both.  And some
observations are not included in the least-squares fit,  and therefore
correspond to zero observations.  The following code looks at an "observation"
and computes the residuals for zero,  one,  or two observations.  The unused
"observations" are left with zero residuals.  The number of residuals determined
is returned.

   About how timing errors are included:  they are assumed to affect optical data
only,  and only in the along-track direction.  We convert the timing error,  in
seconds (obs->time_sigma * seconds_per_day),  and multiply it by the along-track
motion in arcseconds/second (m.total_motion / 60;  m.total_motion is,  perhaps
unwisely,  in arcseconds/minute = arcminutes/hour).

   We then add the resulting 'timing_err' _in quadrature_ with the 'normal'
position error,  to get the total along-track error.  An example:  let's say
our astrometry is good to 0".7,  the object is moving along at 0".6/second,
and the uncertainty on our timing is four seconds.

   That would mean that the cross-track error would be plain old 0".7;  no need
to adjust that,  under our model that timing errors don't affect cross-track
residuals.  But during four seconds,  the object would move 2.4";  that would
be our 'timing_err'.  Added in quadrature,  that would mean that the total
along-track error would be sqrt( 0.7^2 + 2.4^2) = 2".5 (values carefully chosen
to produce an exact result).

   That sort of example,  with a circular uncertainty stretched out into an
ellipse,  is straightforward.  If the existing uncertainty is already
elliptical,  stretching _that_ out is more complex;  see 'errors.cpp'. */

void compute_error_ellipse_adjusted_for_motion( double *sigma1, double *sigma2,
                  double *posn_angle, const OBSERVE *obs,
                  const MOTION_DETAILS *m)
{
   const double timing_err_in_minutes = obs->time_sigma * minutes_per_day;
   const double dy =  timing_err_in_minutes * m->ra_motion;
   const double dx =  timing_err_in_minutes * m->dec_motion;
         /* (dx, dy) = motion of object over obs->time_sigma, in arcsec */

   *sigma1 = obs->posn_sigma_1;    /* start with "non-moving" error ellipse */
   *sigma2 = obs->posn_sigma_2;
   *posn_angle = obs->posn_sigma_theta;
   adjust_error_ellipse_for_timing_error( sigma1, sigma2, posn_angle,
                  dx, dy);
}

int get_residual_data( const OBSERVE *obs, double *xresid, double *yresid)
{
   int n_residuals = 0;

   *xresid = *yresid = 0.;
/* if( obs->is_included)   */
      {
      if( obs->note2 == 'R')
         {
         RADAR_INFO rinfo;

         compute_radar_info( obs, &rinfo);
         if( rinfo.rtt_obs)
            {
            *xresid = (rinfo.rtt_obs - rinfo.rtt_comp) / rinfo.rtt_sigma;
            n_residuals++;
            }
         if( rinfo.doppler_obs)
            {
            *yresid = (rinfo.doppler_obs - rinfo.doppler_comp) / rinfo.doppler_sigma;
            n_residuals++;
            }
         }
      else
         {
         MOTION_DETAILS m;
         double sigma_1, sigma_2, tilt;
         double cos_tilt, sin_tilt;

         compute_observation_motion_details( obs, &m);
         compute_error_ellipse_adjusted_for_motion( &sigma_1, &sigma_2,
                  &tilt, obs, &m);
         cos_tilt = cos( tilt);
         sin_tilt = sin( tilt);
         *xresid = (cos_tilt * m.xresid - sin_tilt * m.yresid);
         *xresid /= sigma_2;
         *yresid = (sin_tilt * m.xresid + cos_tilt * m.yresid);
         *yresid /= sigma_1;
         n_residuals = 2;
         }
      }
   return( n_residuals);
}

static double vect_diff2( const double *a, const double *b)
{
   size_t i;
   double rval = 0, delta;

   for( i = 3; i; i--)
      {
      delta = *a++ - *b++;
      rval += delta * delta;
      }
   return( rval);
}

static void output_json_matrix( FILE *ofile, const char *title, const double *matrix,
               const size_t dim)
{
   size_t i, j;

   fprintf( ofile, "\"%s\": [\n", title);
   for( i = 0; i < dim; i++)
      {
      fprintf( ofile, "   [");
      for( j = 0; j < dim; j++)
         fprintf( ofile, "%.8g%c ", *matrix++, (j == dim - 1) ? ' ' : ',');
      fprintf( ofile, "]%s\n", (i == dim - 1) ? "" : ",");
      }
   fprintf( ofile, "]");
}

const char *monte_label[MONTE_N_ENTRIES] = {
                           "Tp", "e", "q", "Q", "1/a", "i", "M",
                           "omega", "Omega", "MOID", "H" };

/* Describing what 'full_improvement()' does requires an entire separate
file of commentary: see 'full.txt'.  Note,  though,  that this should be
given an orbit that is somewhere within the arc of observations,  for
stability reasons.  Sigmas and/or constraints may be computed for an
entirely different second epoch,  given as "epoch2".  That is to say,
maybe your observations cover a couple of weeks in October 2001,  but
you're interested in the uncertainties as of March 2012.  'epoch' would
be in 2001.  'epoch2' would be in 2012.

   Note that if you're uninterested in sigmas or constraining the orbit,
'epoch2' won't be used for anything.         */

int full_improvement( OBSERVE FAR *obs, int n_obs, double *orbit,
                 const double epoch, const char *limited_orbit,
                 int sigmas_requested, double epoch2)
{
   double *asteroid_mass = ((limited_orbit && *limited_orbit == 'm') ?
               get_asteroid_mass( atoi( limited_orbit + 2)) : NULL);
   int n_params;
   void *lsquare;
   double FAR *xresids;
   double FAR *yresids;
   double FAR *slopes;
   double constraint_slope[MAX_CONSTRAINTS][MAX_N_PARAMS];
   double element_slopes[MAX_N_PARAMS][MONTE_N_ENTRIES];
   double elements_in_array[MONTE_N_ENTRIES];
   double differences[10];
   double original_orbit[6], original_params[MAX_N_NONGRAV_PARAMS];
   double central_obj_state[6], tvect[6];
   const double default_delta_vals[9] =
//                  { 1e-12, 1e-12, 1e-12, 1e-11, 1e-11, 1e-11,
//                  .001, .001, .1 };
                   { 1e-4, 1e-4, 1e-5, 1e-5, 1e-3, 1e-3,
                    .001, .001, .1 };
   static double delta_vals[9];
   double constraint[MAX_CONSTRAINTS];
   double sigma_squared = 0.;       /* see Danby, p. 243, (7.5.20) */
   double scale_factor = 1.;
   double integration_length;
   double before_rms;
   int planet_orbiting = forced_central_body, n_constraints = 0;
   int i, j, n_skipped_obs = 0, err_code = 0;
   int n_included_observations = 0;
   bool really_use_symmetric_derivatives;
   const int n_total_obs = n_obs;
   const char *covariance_filename = "covar.txt";
   char tstr[80];
   ELEMENTS elem;
   OBSERVE *orig_obs = NULL;
   const int showing_deltas_in_debug_file =
                      atoi( get_environment_ptr( "DEBUG_DELTAS"));
   const double r_mult = 1e+2;
   double orbit2[6];
   int set_locs_rval;
   const bool saved_fail_on_hitting_planet =
                                     fail_on_hitting_planet;

   perturbers_automatically_found = always_included_perturbers;
   if( asteroid_mass)                    /* If computing an asteroid mass, */
      {                                  /* be very sure that asteroids are */
      perturbers |= (1 << IDX_ASTEROIDS); /* actually turned on! Otherwise, */
      n_params = 7;                       /* we can loop forever...         */
      }
   else if( limited_orbit && !strncmp( limited_orbit, "np=", 3))
      n_params = atoi( limited_orbit + 3);
   else
      n_params = 6 + n_extra_params;
   if( !obs)
      {
      if( eigenvects)
         {
         free( eigenvects);
         eigenvects = NULL;
         }
      *delta_vals = 0.;
      return( 0);
      }
   if( get_idx1_and_idx2( n_obs, obs, &i, &j) < 3)
      return( -1);
   if( is_unreasonable_orbit( orbit))
      return( -2);
   if( !*delta_vals)    /* force a reset to default values */
      memcpy( delta_vals, default_delta_vals, sizeof( delta_vals));
   available_sigmas = NO_SIGMAS_AVAILABLE;
               /* We save the input orbit;  if there's an error,  we can */
               /* restore it:         */
   memcpy( original_orbit, orbit, 6 * sizeof( double));
   memcpy( original_params, solar_pressure, MAX_N_NONGRAV_PARAMS * sizeof( double));
   sprintf( tstr, "full improvement: %f  ", JD_TO_YEAR( epoch));
   runtime_message = tstr;
   for( i = 0; i < n_obs; i++)
      if( obs->note2 != 'R')
         {
         obs[i].computed_ra  = obs[i].ra;
         obs[i].computed_dec = obs[i].dec;
         }

               /* Drop unincluded observations from the start of the arc: */
   while( n_obs && !obs->is_included)
      {
      obs++;
      n_obs--;
      n_skipped_obs++;
      }
               /* ... and then from the end of the arc: */
   while( n_obs && !obs[n_obs - 1].is_included)
      n_obs--;
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         n_included_observations++;
   if( n_included_observations < 3)
      {
      debug_printf( "full_improvement fail: %d obs\n", n_included_observations);
      runtime_message = NULL;
      return( -3);
      }
   if( limited_orbit && !*limited_orbit)
      limited_orbit = NULL;
   if( n_params < 6 || asteroid_mass)
      limited_orbit = NULL;

            /* If no sigmas are requested,  we compute 'em anyway,  for */
            /* the default display epoch.                               */
   if( sigmas_requested == NO_ORBIT_SIGMAS_REQUESTED)
      {
      sigmas_requested = ORBIT_SIGMAS_REQUESTED;
      if( !limited_orbit)
         epoch2 = find_epoch_shown( obs, n_obs);
//       epoch2 = (obs[0].jd + obs[n_obs - 1].jd) / 2.;
      }

   sprintf( tstr, "fi/setting locs: %f  ", JD_TO_YEAR( epoch));
   fail_on_hitting_planet = true;
   set_locs_rval = set_locs_extended( orbit, epoch, obs, n_obs, epoch2, orbit2);
   fail_on_hitting_planet = saved_fail_on_hitting_planet;
   if( set_locs_rval)
      {
      extern int planet_hit;

      debug_printf( "Hit planet %d in full_improvement : %d\n",
                      planet_hit, set_locs_rval);
      runtime_message = NULL;
      return( -4);
      }

   if( sigmas_requested != HELIOCENTRIC_SIGMAS_ONLY)
      planet_orbiting = find_best_fit_planet( epoch2, orbit2, tvect);
   else
      get_relative_vector( epoch2, orbit2, tvect, planet_orbiting);
   for( i = 0; i < 6; i++)
      {
      central_obj_state[i] = orbit2[i] - tvect[i];
      orbit2[i] = tvect[i];
      }
            /* After the above: if the object,  at the "display epoch", */
            /* is in a heliocentric orbit,  the central object is the   */
            /* "un-moving" sun.  In a planetocentric case,  that vector */
            /* will have the state of the object we're orbiting.  In    */
            /* either case,  orbit2[] will now give the state vector,   */
            /*  at epoch2,  of our object,  relative to planet_orbiting */
            /* (which may be 0, i.e.,  heliocentric.)                   */

   if( limited_orbit)
      n_constraints = evaluate_limited_orbit( orbit2, planet_orbiting,
                                    epoch2, limited_orbit, constraint);

            /* evaluate elements of 'orbit2',  then put */
            /* into an array form: */
   elem.gm = get_planet_mass( planet_orbiting);
   elem.abs_mag = calc_absolute_magnitude( obs, n_obs);
   calc_classical_elements( &elem, orbit2, epoch2, 1);

   put_orbital_elements_in_array_form( &elem, elements_in_array);

   uncertainty_parameter = 99.;
   xresids = (double FAR *)FCALLOC( (2 + 2 * n_params) * n_obs + n_params, sizeof( double));
   yresids = xresids + n_obs;
   slopes = yresids + n_obs;

   before_rms = compute_rms( obs, n_obs);
   if( limited_orbit && *limited_orbit == 'R')
      constraint[n_constraints++] =
              r_mult * (dotted_dist( obs + n_obs - 1) - atof( limited_orbit + 2));

   sprintf( tstr, "fi/locs set: %f  ", JD_TO_YEAR( epoch));
   for( i = 0; i < n_obs; i++)
      get_residual_data( obs + i, xresids + i, yresids + i);

             /* 'integration_length' = maximum time span over which we'll */
             /* be integrating,  from the working epoch to either the first */
             /* or last observation.  We use it to scale the 'delta_val'    */
             /* tweaks.                                                     */
   integration_length = fabs( epoch - obs[n_obs - 1].jd);
   if( integration_length < fabs( epoch - obs->jd))
      integration_length = fabs( epoch - obs->jd);
               /* If the arc length is less than a week,  insist on the use */
               /* of symmetric derivatives;  it seems to add stability.  If */
               /* not,  just go with the default 'use' variable.            */
   if( obs[n_obs - 1].jd - obs[0].jd < 7.)
      really_use_symmetric_derivatives = true;
   else
      really_use_symmetric_derivatives = use_symmetric_derivatives;
   orig_obs = (OBSERVE *)calloc( n_obs, sizeof( OBSERVE));
   memcpy( orig_obs, obs, n_obs * sizeof( OBSERVE));

   for( i = 0; i < n_params; i++)
      {
      const double min_change = 0.3, max_change = 3.0, optimal_change = 1.0;
      double low_delta = 0., high_delta = 0., low_change = 0., high_change = 0.;
      int n_iterations = 0;
      const int max_iterations = 100;
      bool keep_iterating = true;

      while( keep_iterating)
         {
         double tweaked_orbit[6];
         const double original_asteroid_mass = (asteroid_mass ? *asteroid_mass : 0.);
         double delta_val =
                      delta_vals[i] / (integration_length * integration_length);
         double worst_error_in_sigmas;
         double worst_error_squared = 0, rescale;
         double original_solar_pressure[MAX_N_NONGRAV_PARAMS];
         double *slope_ptr;
         double rel_orbit[6];

                  /* for asteroid mass computations,  on first pass, */
                  /* try to set a "reasonable" delta :   */
         if( i == 6 && asteroid_mass && !n_iterations)
            delta_val = 1.e-15 + original_asteroid_mass / 100.;
         memcpy( original_solar_pressure, solar_pressure,
                               n_extra_params * sizeof( double));
         do
            {
            memcpy( tweaked_orbit, orbit, 6 * sizeof( double));
            memcpy( solar_pressure, original_solar_pressure,
                                n_extra_params * sizeof( double));
            if( i < 6)              /* adjust position/velocity */
               tweaked_orbit[i] -= delta_val;
            else
               {
               if( asteroid_mass)
                  *asteroid_mass -= delta_val;
               else
                  solar_pressure[i - 6] -= delta_val;
               }
            sprintf( tstr, "Evaluating %d of %d : iter %d   ", i + 1,
                                    n_params, n_iterations);
            if( debug_level > 4)
               debug_printf( "About to set locs #2: delta_val %f\n", delta_val);
            fail_on_hitting_planet = true;
            set_locs_rval = set_locs_extended( tweaked_orbit, epoch, obs,
                       n_obs, epoch2, rel_orbit);
            fail_on_hitting_planet = saved_fail_on_hitting_planet;
            for( j = 0; !set_locs_rval && j < n_obs; j++)
               {
               double d = vect_diff2( obs[j].obj_posn, orig_obs[j].obj_posn);

               d = sqrt( d) / obs[j].r;
               if( d > 0.001)              /* tweaked orbit too different from */
                  set_locs_rval = 1;      /* the original one; try smaller tweak */
               }
            if( debug_level > 4)
               debug_printf( "Second set done: %d\n", set_locs_rval);
            if( set_locs_rval == INTEGRATION_TIMED_OUT)
               {
               free( xresids);
               memcpy( orbit, original_orbit, 6 * sizeof( double));
               memcpy( solar_pressure, original_params, n_extra_params * sizeof( double));
               runtime_message = NULL;
               debug_printf( "Integration timeout\n");
               return( -4);
               }
            if( set_locs_rval)      /* gonna have to try again, */
               {                    /* with a smaller tweak */
               delta_val /= 2.;
               delta_vals[i] /= 2.;
               }
            }
            while( set_locs_rval);
         slope_ptr = slopes + i;
         for( j = 0; j < n_obs; j++, slope_ptr += 2 * n_params)
            get_residual_data( obs + j, slope_ptr, slope_ptr + n_params);

         for( j = 0; j < 6; j++)
            rel_orbit[j] -= central_obj_state[j];
                  /* evaluate elements of 'rel_orbit',  then  put */
                  /* into an array form: */
         elem.gm = get_planet_mass( planet_orbiting);
         elem.abs_mag = calc_absolute_magnitude( obs, n_obs);
         calc_classical_elements( &elem, rel_orbit, epoch2, 1);

         put_orbital_elements_in_array_form( &elem, element_slopes[i]);
         for( j = 0; j < MONTE_N_ENTRIES; j++)
            {
            element_slopes[i][j] -= elements_in_array[j];
            element_slopes[i][j] /= delta_val;
            }
         if( limited_orbit)
            {
            double constraint2[MAX_CONSTRAINTS];

            evaluate_limited_orbit( rel_orbit, planet_orbiting, epoch2,
                                      limited_orbit, constraint2);
            for( j = 0; j < n_constraints; j++)
               constraint_slope[j][i] =
                       (constraint2[j] - constraint[j]) / delta_val;
            if( *limited_orbit == 'R')
               {
               const double tconstraint =
                     r_mult * (dotted_dist( obs + n_obs - 1) - atof( limited_orbit + 2));

               constraint_slope[0][j] = (constraint[0] - tconstraint) / delta_val;
               }
            }
         if( really_use_symmetric_derivatives)
            {
            for( j = 0; j < 6; j++)
               tweaked_orbit[j] = 2. * orbit[j] - tweaked_orbit[j];
            if( asteroid_mass)
               *asteroid_mass = 2. * original_asteroid_mass - *asteroid_mass;
            else
               for( j = 0; j < n_extra_params; j++)
                  solar_pressure[j] = 2. * original_solar_pressure[j] -
                                    solar_pressure[j];
            memcpy( tstr, "Reverse   ", 10);
            fail_on_hitting_planet = true;
            set_locs_rval = set_locs( tweaked_orbit, epoch, obs, n_obs);
            if( set_locs_rval)      /* fall back on simple, asymmetric method */
               {
               debug_printf( "Symmetric fail : %d\n", set_locs_rval);
               memcpy( obs, orig_obs, n_obs * sizeof( OBSERVE));
               }
            fail_on_hitting_planet = saved_fail_on_hitting_planet;
            }
         else
            memcpy( obs, orig_obs, n_obs * sizeof( OBSERVE));
         slope_ptr = slopes + i;
         for( j = 0; j < n_obs; j++, slope_ptr += 2 * n_params)
            if( obs[j].is_included)
               {
               double xresidual, yresidual;

               get_residual_data( obs + j, &xresidual, &yresidual);

               slope_ptr[0] -= xresidual;
               slope_ptr[n_params] -= yresidual;
//             if( obs[j].note2 != 'R')
                  {
                  const double error_squared = slope_ptr[0] * slope_ptr[0]
                           + slope_ptr[n_params] * slope_ptr[n_params];

                  if( worst_error_squared < error_squared)
                     worst_error_squared = error_squared;
                  }
               slope_ptr[0]        /= delta_val;
               slope_ptr[n_params] /= delta_val;
               if( really_use_symmetric_derivatives && !set_locs_rval)
                  {       /* delta is actually twice the 'specified' value */
                  slope_ptr[0]        /= 2.;
                  slope_ptr[n_params] /= 2.;
                  }
               }
         worst_error_in_sigmas = sqrt( worst_error_squared);
         if( showing_deltas_in_debug_file)
            debug_printf( "Iter %d, Change param %d: %f sigmas; delta %.3e (%.3e)\n",
               n_iterations,
               i, worst_error_in_sigmas, delta_val, delta_vals[i]);
         worst_error_in_sigmas += 1e-10;        /* ensure _some change; */
                           /* evades divide-by-zero/range errors below */
         if( worst_error_in_sigmas > min_change && worst_error_in_sigmas < max_change)
            keep_iterating = false;
         if( worst_error_in_sigmas <= optimal_change)
            {
            low_delta = delta_vals[i];
            low_change = worst_error_in_sigmas;
            }
         else
            {
            high_delta = delta_vals[i];
            high_change = worst_error_in_sigmas;
            }
         rescale = optimal_change / worst_error_in_sigmas;
         if( low_change && high_change)  /* we've got it bracketed */
            {
            const double slope = log( high_delta / low_delta)
                               / log( high_change / low_change);

            rescale = exp( log( rescale) * slope);
            }
         delta_vals[i] *= rescale;
         memcpy( solar_pressure, original_solar_pressure,
                                n_extra_params * sizeof( double));
         if( asteroid_mass)
            *asteroid_mass = original_asteroid_mass;
         if( n_iterations++ >= max_iterations)
            {
            debug_printf( "Ran over iteration limit! %s\n", obs->packed_id);
            debug_printf( "Worst err %f sigmas\n", worst_error_in_sigmas);
            free( xresids);
            memcpy( orbit, original_orbit, 6 * sizeof( double));
            memcpy( solar_pressure, original_params, n_extra_params * sizeof( double));
            runtime_message = NULL;
            return( -4);
            }
         }
      }
   free( orig_obs);

   lsquare = lsquare_init( n_params);
   assert( lsquare);
   if( debug_level > 1)
      debug_printf( "Adding obs to lsquare\n");
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         {
         double loc_vals[22], weight = 1.;
         const double xresid = xresids[i];
         const double yresid = yresids[i];      /* all in _radians_ */
         const double resid2 = xresid * xresid + yresid * yresid;

         if( !(obs[i].flags & OBS_ALREADY_CORRECTED_FOR_OVEROBSERVING))
            {
            if( use_blunder_method == 2 && probability_of_blunder)
               weight = reweight_for_blunders( resid2, weight);
            if( overobserving_time_span && overobserving_ceiling)
               weight *= reweight_for_overobserving( obs, n_obs, i);
            }
         FMEMCPY( loc_vals, slopes + i * 2 * n_params,
                                         2 * n_params * sizeof( double));
         lsquare_add_observation( lsquare, xresid, weight, loc_vals);
         lsquare_add_observation( lsquare, yresid, weight, loc_vals + n_params);
         sigma_squared += weight * weight * (resid2 + 1.);
         }
   i = n_included_observations * 2 - n_params;
   if( i > 0)
      sigma_squared /= (double)i;

   if( limited_orbit)
      for( j = 0; j < n_constraints; j++)
         lsquare_add_observation( lsquare, constraint[j], 1.,
                                            constraint_slope[j]);

   if( debug_level > 1)
      debug_printf( "lsquare solve\n");
   if( !err_code)
      {
      err_code = lsquare_solve( lsquare, differences);
      if( err_code)
         debug_printf( "Failure in lsquare_solve: %d\n", err_code);
      }

   if( debug_level > 1)
      debug_printf( "Making covar file\n");
   if( !err_code && *covariance_filename)
      {
      char tbuff[200];
      FILE *ofile = fopen_ext( get_file_name( tbuff, covariance_filename), "tfcwb");
      FILE *json_ofile = fopen_ext( get_file_name( tbuff, "covar.json"), "tfcwb");
      double *matrix = lsquare_covariance_matrix( lsquare);
      double *wtw = lsquare_wtw_matrix( lsquare);
      double eigenvals[10], eigenvectors[100], element_sigmas[MONTE_N_ENTRIES];
      int pass;
      const int max_obs_in_covariance_file = 2000;

      assert( matrix);
      setvbuf( ofile, NULL, _IONBF, 0);
      fprintf( ofile, "Orbit: %.7f %.7f %.7f %.7f %.7f %.7f\nepoch JD %.5f\n",
               orbit[0], orbit[1], orbit[2],
               orbit[3], orbit[4], orbit[5], epoch2);
      jacobi_eigenvalues( wtw, n_params, eigenvals, eigenvectors);

      fprintf( ofile, "Eigenvalues computed: sigma_squared = %g\n", sigma_squared);
      assert( sigma_squared);
      for( i = 0; i < n_params; i++)
         eigenvals[i] /= sigma_squared;
      for( i = 0; i < n_params * n_params; i++)
         matrix[i] *= sigma_squared;      /* Danby, p 243, (7.5.21) */
      for( i = 0; i < n_obs && i < max_obs_in_covariance_file; i++)
         if( obs[i].is_included)
            {
            fprintf( ofile, "%4d: ", i);
            for( j = 0; j < n_params * 2; j++)
               {
               put_double_in_buff( tbuff, slopes[i * 2 * n_params + j]);
               fprintf( ofile, "%s", tbuff);
               }
            fprintf( ofile, "\n");
            }
      fprintf( json_ofile, "{ ");
      output_json_matrix( json_ofile, "covar", matrix, n_params);
      fprintf( json_ofile, ", \"state_vect\": [\n");
      for( i = 0; i < 6; i++)
         fprintf( json_ofile, "    %.12g%s\n", orbit[i],
                        (i == 5 ? " ]," : ","));
      fprintf( json_ofile, "  \"epoch\": %f\n}", epoch);
      fclose( json_ofile);
      for( pass = (matrix ? 0 : 2); pass < 4; pass++)
         {
         const char *titles[] = {  "Covariance", "Correlation",
                                               "WtW", "Eigenvectors" };

         fprintf( ofile, "%s:\n", titles[pass]);
         for( i = 0; i < n_params; i++)
            {
            for( j = 0; j < n_params; j++)
               {
               double oval;

               if( pass == 0 || pass == 1)
                  {
                  oval = matrix[i + j * n_params];
                  if( pass == 1)      /* normalize to get correlation, not covar */
                     oval /= sqrt( matrix[i * (n_params + 1)] * matrix[j * (n_params + 1)]);
                  }
               else if( pass == 2)
                  oval = wtw[i + j * n_params];
               else        /* if( pass == 3) */
                  oval = eigenvectors[j + i * n_params];
               if( pass == 1 || pass == 3)      /* correlation or eigenvects */
                  sprintf( tbuff, "%10.6f", oval);  /* values are -1 to 1 */
               else                             /* covar/WtW values can be */
                  put_double_in_buff( tbuff, oval);   /* huge or tiny */
               fprintf( ofile, "%s", tbuff);
               }
            fprintf( ofile, "\n");
            }

         if( pass != 1)
            {
            fprintf( ofile, "%s: test for orthogonality\n", titles[pass]);
            for( i = 0; i < n_params; i++)
               {
               for( j = 0; j < n_params; j++)
                  {
                  double oval = 0.;
                  unsigned k;
                  double *matrix_ptr;

                  if( pass == 0)
                     matrix_ptr = matrix;
                  else if( pass == 2)
                     matrix_ptr = wtw;
                  else
                     matrix_ptr = eigenvectors;
                  for( k = 0; k < (unsigned)n_params; k++)
                     oval += matrix_ptr[i * n_params + k] * matrix_ptr[j * n_params + k];
                  if( pass == 3)      /* eigenvects are normalized; */
                     sprintf( tbuff, "%10.6f", oval);  /* values are -1 to 1 */
                  else                             /* covar/WtW values can be */
                     put_double_in_buff( tbuff, oval);   /* huge or tiny */
                  fprintf( ofile, "%s", tbuff);
                  }
               fprintf( ofile, "\n");
               }
            }
         }
      fprintf( ofile, "Eigenvalues:\n");
      for( i = 0; i < n_params; i++)
         {
         put_double_in_buff( tbuff, eigenvals[i]);
         fprintf( ofile, "%s", tbuff);
         }
      if( eigenvects)
         {
         free( eigenvects);
         eigenvects = NULL;
         }
      if( n_params >= 6)
         {
         eigenvects       = (double **)calloc_double_dimension_array(
                                n_params, n_params, sizeof( double));

         for( i = 0; i < n_params; i++)
            for( j = 0; j < n_params; j++)
               eigenvects[j][i] = eigenvectors[i + j * n_params] / sqrt( eigenvals[j]);
         fprintf( ofile, "\nOne-sigma eigenvect:\n");
         for( i = 0; i < n_params; i++)
            {
            put_double_in_buff( tbuff, eigenvects[0][i]);
            fprintf( ofile, "%s", tbuff);
            }
         }

      fprintf( ofile, "\n\nCovariance in 'traditional' elements:\n");
      *tbuff = '\0';
      for( i = 0; i < MONTE_N_ENTRIES; i++)
         snprintf_append( tbuff, sizeof( tbuff), "   %-7s", monte_label[i]);
      fprintf( ofile, "%s\n", tbuff);    /* label element covariance at top */
                        /* also show units */
      fprintf( ofile, " (days)    (unitless) (AU)      (AU)      (1/AU)"
              "   (deg)     (deg)       (deg)     (deg)     (AU)    (mag)\n");
      for( i = 0; i < MONTE_N_ENTRIES; i++)
         {
         for( j = 0; j < MONTE_N_ENTRIES; j++)
            {
            int k, l;
            double covar_elem = 0.;
            char covar_text[20];

            for( l = 0; l < n_params; l++)
               {
               double dot = 0.;

               for( k = 0; k < n_params; k++)
                  dot += element_slopes[k][i] * matrix[k + l * n_params];
               covar_elem += dot * element_slopes[l][j];
               }
            put_double_in_buff( covar_text, covar_elem);
            fprintf( ofile, "%s", covar_text);
            if( i == j)       /* on the diagonal of the covariance matrix: */
               {              /* we can compute sigmas here */
               element_sigmas[i] = (covar_elem > 0. ? sqrt( covar_elem) : 0.);
               }              /* roundoff gets us below zero sometimes */
            }
         fprintf( ofile, " %s\n", monte_label[i]);
         }
      fprintf( ofile, "%s\n", tbuff);   /* label element covariance at end */

      for( i = 0; i < MONTE_N_ENTRIES; i++)
         {
         put_double_in_buff( tbuff, element_sigmas[i]);
         fprintf( ofile, "\n   %-6s %s %.9f", monte_label[i], tbuff, elements_in_array[i]);
         fprintf( ofile, "\n  ");
         for( j = 0; j < n_params; j++)
            {
            put_double_in_buff( tbuff, element_slopes[j][i]);
            fprintf( ofile, "%s", tbuff);
            }
         }
      for( i = 6; i < n_params; i++)
         {
         double sigma;
         char title_text[20];

         sigma_squared = matrix[i + i * n_params];
         assert( sigma_squared >= 0.);
         sigma = sqrt( sigma_squared);
         if( n_params == 7 && !asteroid_mass)
            sigma *= SOLAR_GM / SRP1AU;
         put_double_in_buff( tbuff, sigma);

         if( n_params == 7)
            strcpy( title_text, (asteroid_mass ? "Sigma_mass" : "Sigma_AMR"));
         else
            {           /* comet A1, A2, maybe A3 included */
            assert( n_params == 8 || n_params == 9);
            sprintf( title_text, "Sigma_A%d", i - 5);
            }
         fprintf( ofile, "\n%s: %s", title_text, tbuff);
         }

      fprintf( ofile, "\n\n");
      if( limited_orbit && strstr( limited_orbit, "e="))
         element_sigmas[MONTE_ECC] = 0.;
      uncertainty_parameter = dump_monte_data_to_file( ofile, element_sigmas,
            elem.major_axis, elem.ecc, planet_orbiting);
      fclose( ofile);
      free( matrix);
      free( wtw);
      if( eigenvals[0] > 0.)
         available_sigmas = COVARIANCE_AVAILABLE;
      else
         available_sigmas = NO_SIGMAS_AVAILABLE;
      available_sigmas_hash = compute_available_sigmas_hash( obs, n_obs, epoch2,
                  perturbers, planet_orbiting);
      }
   FFREE( xresids);
   lsquare_free( lsquare);
   for( i = 0; !err_code && i < 6 && i < n_params; i++)
      for( j = 0; j < 6; j++)
         {
         double max_difference = obs->r * .7, ratio;

         if( j > 2)     /* velocity component */
            max_difference /= (obs[n_obs - 1].jd - obs[0].jd) * .5;
         if( i == j)
            ratio = fabs( differences[i] / max_difference);
         else
            ratio = 0.;
         if( ratio > scale_factor)
            scale_factor = ratio;
         }
   if( debug_level > 1)
      debug_printf( "lsquare computed\n");
// if( scale_factor != 1.)
//    debug_printf( "rescale: %g\n", scale_factor);
   memcpy( orbit2, orbit, 6 * sizeof( double));
   for( i = 0; i < n_params && !err_code; i++)
      {
      if( i < 6)
         orbit[i] += differences[i] / scale_factor;
      if( i == 5)    /* is our new 'orbit' state vector reasonable?  */
         err_code = is_unreasonable_orbit( orbit);

      if( i >= 6)
         {
         if( asteroid_mass)
            *asteroid_mass += differences[i] / scale_factor;
         else
            solar_pressure[i - 6] += differences[i] / scale_factor;
         }
      }
               /* If the orbit "blew up" or otherwise failed,  restore */
               /* the original version:  */
   if( err_code || is_unreasonable_orbit( orbit))
      debug_printf( "Failed full step: %d: %s\n", err_code, obs->packed_id);
   sprintf( tstr, "Final setting of orbit    ");
   i = 6;
   do
      {
      if( setting_outside_of_arc)
         set_locs( orbit, epoch, obs - n_skipped_obs, n_total_obs);
      else
         set_locs( orbit, epoch, obs, n_obs);
      if( *get_environment_ptr( "HALF_STEPS"))
         {
         const double after_rms = compute_rms( obs, n_obs);

         sprintf( tstr, "Half-stepping %d\n", 7 - i);
         if( after_rms > before_rms * 1.5 && !limited_orbit)
            {
            for( j = 0; j < 6; j++)
               orbit[j] = (orbit[j] + orbit2[j]) * .5;
            i--;
            }
         else
            i = 0;
         }
      else
         i = 0;
      }
      while( i);

   if( debug_level > 1)
      debug_printf( "full_improve done\n");
   runtime_message = NULL;
   return( err_code);
}


/* 'score_orbit_arc' basically looks at a series of observations and
assigns a 'score' indicating how good an orbit we expect it can
produce.  After trying various schemes,  I've settled on one that
essentially measures the departure of an arc from a great circle.
To do this with three vectors A, B, C,  one would compute (AxB).C,
i.e.,  cross-product and then dot-product.  Here,  A and B are
the vectors to the first and last observations,  and we check all
intervening observations to see how far out of that plane they are.
We look for both a maximum and a minimum;  the difference between
them indicates how far out of the great-circle route the observations
go.         */

static inline double score_orbit_arc( const OBSERVE FAR *obs, const unsigned n_obs)
{
   const double five_degrees = PI * 5. / 180.;
   double xprod[3], minimum = 0., maximum = 0.;
   double max_dist = vector3_dist( obs[0].vect, obs[n_obs - 1].vect);
   unsigned i;

   if( n_obs == 1)
      return( -1.);

   vector_cross_product( xprod, obs[0].vect, obs[n_obs - 1].vect);

   for( i = 1; i < n_obs - 1; i++)
      {
      const double dot = dot_product( xprod, obs[i].vect);

      if( maximum < dot)
         maximum = dot;
      else if( minimum > dot)
         minimum = dot;
      if( max_dist < five_degrees)
         {
         const double dist1 = vector3_dist( obs[i].vect, obs[0].vect);
         const double dist2 = vector3_dist( obs[i].vect, obs[n_obs - 1].vect);

         if( max_dist < dist1)
            max_dist = dist1;
         if( max_dist < dist2)
            max_dist = dist2;
         }
      }
   if( max_dist > five_degrees)
      max_dist = five_degrees;
   return( (maximum - minimum) * max_dist);
}

/* When looking for subarcs,  we want to look just for those that are less
than 45 degrees long.  Longer than that,  and the usual initial orbit
methods (Gauss,  Herget) may diverge.  (Not necessarily -- in fact,
they'll often still work with much longer arcs -- but 45 degrees means
you can be pretty confident that a decent orbit will be found.)  A quick
dot-product of the unit vectors,  compared to cos(45 degrees) =
sqrt(2) / 2.,  suffices to check the arc length. */

static inline void look_for_best_subarc( const OBSERVE FAR *obs,
       const int n_obs, const double max_arc_len, int *start, int *end)
{
   double best_score = -999., score;
   const double cos_45_deg = 1.414213 / 2.;
   int i, j;

   *start = *end = 0;
   for( i = j = 0; i < n_obs - 1; i++)
      {
      int temp_start, temp_end;

      while( j < n_obs - 1 && obs[j + 1].jd - obs[i].jd <= max_arc_len
                  && dot_product( obs[i].vect, obs[j + 1].vect) > cos_45_deg)
         j++;
      if( obs[i].jd == obs[j].jd)   /* rare case of two or more obs with */
         j = i;            /* the same JD;  this can cause trouble,  so  */
      temp_start = i;      /* cut arc down to a single observation       */
      temp_end = j;
      while( temp_start < temp_end && (obs[temp_start].flags & OBS_DONT_USE))
         temp_start++;
      while( temp_start < temp_end && (obs[temp_end].flags & OBS_DONT_USE))
         temp_end--;
      score = score_orbit_arc( obs + temp_start, temp_end - temp_start + 1);
      if( score > best_score)
         {
         best_score = score;
         *start = temp_start;
         *end = temp_end;
         }
      }
}

/* The following was created by plotting semimajor axes versus
inclinations as a scatterplot,  from the recent 'mpcorb.dat' file.
The result is a 30x60 graph,  with 30 vertical bins representing
2 degrees each in inclination (from incl=0 at the top to incl=60
at the bottom),  and 60 bins across representing .1 AU in semimajor
axis (from 0 AU to 6 AU).  The natural logarithm of the counts in
each bin are shown;  those ranged from zero to 9,  so I didn't need
to do any scaling.  This basically reflects the scatterplot provided
by MPC at https://www.minorplanetcenter.net/iau/plot/OrbEls51.gif
(except that the vertical axis is flipped).   */

/*           0 AU      1         2         3         4         5         6 */
static const char *incl_vs_a_scattergram[30] = {
/* incl=0*/ "000000012222232223333778877777776331002520000000001332000000",
            "000000112333233333334789988888886331202630100000002442000000",
            "000000112223333233334799888876776332111520000000002442000000",
            "000000022233333333333699888876776342111520010000002442000000",
            "000000011222233233333577787877886443011510000000002442000000",
/* 10 deg*/ "000000011222223323332357678868886443101410000010001442000000",
            "000000012223332232222246689866776331000410000000002441000000",
            "000000012222323222332124688766776331100310000000001340000000",
            "000000002222222223462222466655786331011200000000001341000000",
            "000000012122223223562233355544675332000100000000000431000000",
/* 20 deg*/ "000000011122232223562255456433565210010100000000001331000000",
            "000000111222222123662366555423565210100100000000000331000000",
            "000000011122223222562356444533464100000000000000001330000000",
            "000000012112212222442244355422464100000000000000000321000000",
            "000000011112122122231222156421353110000000000000000321000000",
/* 30 deg*/ "000000011112111122221111245410132000000000000000000230000000",
            "000000011110011212211111135410122000000000000000000220000000",
            "000000001111211221221120033320011100000000000000000120000000",
            "000000000112111211201110122310000000000000000000000110000000",
            "000000000011122212201110111110000000000000000000000000000000",
/* 40 deg*/ "000000001101111121111110101100001000000000000000000000000000",
            "000000000011010001101110000001000000000000000000000000000000",
            "000000000001100001101010010100000000000000000000000000000000",
            "000000000000011011100110010000000000000000000000000000000000",
            "000000000100100000000000000000000000000000000000000000000000",
/* 50 deg*/ "000000001000100100100000000000000000000000000000000000000000",
            "000000000000100011010100000000000000000000000000000000000000",
            "000000000000011000110000000000000000000000000000000000000000",
            "000000000000010000000011000000000000000000000000000000000000",
            "000000000000000000000000000000000000000000000000000000000000" };

/* evaluate_initial_orbit( ) is supposed to give a "score" of sorts,      */
/* describing just how likely this orbit seems to be,  given its rms      */
/* errors (lower is better);  eccentricity (same,  with highly hyperbolic */
/* solutions considered very unlikely and therefore bumping up the return */
/* value);  inclination (anything above .5 radians,  or about 30 degrees, */
/* is considered unlikely).  Also,  main-belt objects (based on a) are    */
/* slightly encouraged.                                                   */
/*    Note that if the designation is that of an interstellar object,  or */
/* if the line 'COM interstellar' is in the observation file,  we don't   */
/* add a penalty for highly hyperbolic solutions.  In that case,  a high  */
/* eccentricity is exactly what we expect to see.                         */

static double adjustment_for_orbit_likelihood( const double semimajor_axis,
                 const double inclination, const double q)
{
   const int n_xbins = 60, n_ybins = 30;
   double xbin = semimajor_axis * 10.;
   double ybin = inclination * (180. / PI) / 2.;
   const int ixbin = (int)xbin, iybin = (int)ybin;
   double rval;

   if( ixbin >= 0 && iybin >= 0 && ixbin < n_xbins - 1 && iybin < n_ybins - 1
               && q > 1.1)
      {
      const char *tptr1 = incl_vs_a_scattergram[iybin] + ixbin;
      const char *tptr2 = incl_vs_a_scattergram[iybin + 1] + ixbin;

      xbin -= (double)ixbin;
      ybin -= (double)iybin;
      rval = (double)( *tptr1 - '0')
               + (double)( tptr1[1] - tptr1[0]) * xbin
               + (double)( tptr2[0] - tptr1[0]) * ybin
               + (double)( tptr2[1] + tptr1[0] - tptr1[1] - tptr2[0])
                              * xbin * ybin;
//    debug_printf( "xbin %f, ybin %f: %f\n",
//             xbin + (double)ixbin,
//             ybin + (double)iybin, rval);
      if( q < 1.2)
         rval *= (q - 1.1) / 10.;
      }
   else
      rval = 0.;
   return( rval * .1);
}

int is_interstellar = 0;

double evaluate_initial_orbit( const OBSERVE FAR *obs,
                              const int n_obs, const double *orbit)
{
   const double rms_err = compute_weighted_rms( obs, n_obs, NULL);
   double rval, rel_orbit[6], planet_radius_in_au;
   ELEMENTS elem;
   int planet_orbiting = find_best_fit_planet( obs->jd,
                                  orbit, rel_orbit);

   elem.gm = get_planet_mass( planet_orbiting);
   calc_classical_elements( &elem, rel_orbit, obs[0].jd, 1);
   if( planet_orbiting && elem.ecc > 1.01)
      {        /* it's flying past a planet:  re-evaluate as an */
      planet_orbiting = 0;               /* heliocentric object */
      elem.gm = get_planet_mass( 0);
      calc_classical_elements( &elem, orbit, obs[0].jd, 1);
      }
   planet_radius_in_au =
          planet_radius_in_meters( planet_orbiting) / AU_IN_METERS;

   rval = rms_err / 2.;
               /* For arcs with very few observations,  we really shouldn't */
               /* think too long and hard about the rms error :             */
   if( n_obs < 6)
      rval *= (double)(n_obs - 2) / 4.;
   if( !planet_orbiting)                  /* for heliocentric orbits... */
      {
      rval += elem.ecc / 2.;
      if( elem.ecc > 1.01 && !is_interstellar) /* _strongly_ discourage hyperbolics */
         rval += (elem.ecc - 1.01) * 1000.;
      if( elem.incl > .5 && elem.ecc < .8)   /* high-incl,  non-cometlike orbits */
         {                                    /* are not very likely */
         double high_inclination_penalty = (elem.incl - .5) * .2;

         if( high_inclination_penalty > .2)  /* nearly 90 degrees */
            high_inclination_penalty = .2;
         if( elem.ecc > .6)
            high_inclination_penalty *= (.8 - elem.ecc) / .2;
         rval += high_inclination_penalty;
         }
      rval -= adjustment_for_orbit_likelihood( elem.major_axis, elem.incl, elem.q);
      }
          /* strongly discourage elliptical orbits going through planets : */
   if( elem.ecc < 1. && elem.q < planet_radius_in_au)
      rval += 10000.;
   if( debug_level > 4)
      debug_printf( "Orbit with a=%f, q=%f, e=%f, i=%f: %f\n",
            elem.major_axis, elem.q, elem.ecc, elem.incl * 180. / PI, rval);
   return( rval);
}

static double attempt_improvements( double *orbit, OBSERVE *obs, const int n_obs)
{
   int method;
   double curr_score;
   OBSERVE *best_obs = (OBSERVE *)calloc( n_obs, sizeof( OBSERVE));
   double temp_orbit[6];

#ifdef CONSOLE
   if( show_runtime_messages)
      move_add_nstr( 14, 10, "Improving solution...        ", -1);
#endif
   if( set_locs( orbit, obs[0].jd, obs, n_obs))
      {
      debug_printf( "Set loc fail 17\n");
      }

   memcpy( temp_orbit, orbit, 6 * sizeof( double));
   curr_score = evaluate_initial_orbit( obs, n_obs, orbit);
   memcpy( best_obs, obs, n_obs * sizeof( OBSERVE));
   for( method = 0; method < 2; method++)
      {
      int iter = 0;
      int max_iter = 5;
      bool error_occurred = false;

                        /* We're willing to try the Herget,  then full  */
                        /* step methods, five times...                  */
      while( !error_occurred && iter++ < max_iter)
         {
         double score;

#ifdef CONSOLE
         if( show_runtime_messages)
            {
            char msg_buff[80];

            sprintf( msg_buff, "%s step: radii %f, %f",
                        (method ? "full" : "Herget"),
                        obs[0].r, obs[n_obs - 1].r);
            move_add_nstr( 14, 10, msg_buff, -1);
            }
#endif
         if( !method)         /* doing an Herget step */
            {
            double r1 = obs[0].r, r2 = obs[n_obs - 1].r;
            double d_r1, d_r2;

            if( herget_method( obs, n_obs, r1, r2, temp_orbit, &d_r1, &d_r2, NULL))
               error_occurred = true;
            else
               {
               r1 += d_r1;
               r2 += d_r2;
               herget_method( obs, n_obs, r1, r2, temp_orbit, NULL, NULL, NULL);
               if( adjust_herget_results( obs, n_obs, temp_orbit))
                  error_occurred = true;
               }
            if( debug_level > 3)
               debug_printf( "Adjusting Herget results\n");
            }
         else        /* doing a full step */
            if( full_improvement( obs, n_obs, temp_orbit, obs->jd, NULL,
                           NO_ORBIT_SIGMAS_REQUESTED, obs->jd))
               {
               debug_printf( "Full improvement failure! %s\n", obs->packed_id);
               error_occurred = true;
               }

         if( !error_occurred)
            {
            const double full_step_advantage = .3;

            score = evaluate_initial_orbit( obs, n_obs, temp_orbit);
            if( method)    /* for a full-improvement step,  make the score */
               score -= full_step_advantage;         /* just a hair better */
            if( debug_level > 2)
               debug_printf( "Method %d, run %d: score %f\n",
                                    method, iter, score);
                     /* If things totally fell apart, stop:   */
            if( score > 10000.)
               error_occurred = true;
                     /* If the step is an improvement,  or shows just a little */
                     /* fluctuation (that happens),  accept it : */
            if( score < curr_score + .01)
               {
               memcpy( orbit, temp_orbit, 6 * sizeof( double));
               curr_score = score;
               memcpy( best_obs, obs, n_obs * sizeof( OBSERVE));
               }
            else
               available_sigmas = NO_SIGMAS_AVAILABLE;
            }
         else              /* avoid showing sigmas from a previous try: */
            available_sigmas = NO_SIGMAS_AVAILABLE;
         }
      }
   memcpy( obs, best_obs, n_obs * sizeof( OBSERVE));
   free( best_obs);
   return( curr_score);
}

static void exclude_unusable_observations( OBSERVE *obs, int n_obs)
{
   while( n_obs--)
      {
      if( obs->flags & OBS_DONT_USE)
         obs->is_included = 0;
      obs++;
      }
}

static inline bool is_solar_observatory( const char *mpc_code)
{
   return( !strcmp( mpc_code, "249") || !strcmp( mpc_code, "C49")
                                     || !strcmp( mpc_code, "C50"));
}

static double find_sungrazer_orbit( OBSERVE FAR *obs, int n_obs, double *orbit)
{
   double best_score = 1e+10;
   double temp_orbit[6];
   double r;
   int i, direction, start, end;

               /* Don't try to do a multiple-orbit case at first : */
   look_for_best_subarc( obs, n_obs, 10., &start, &end);
   for( i = 0; i < start; i++)
      obs[i].is_included = 0;
   for( i = end + 1; i < n_obs; i++)
      obs[i].is_included = 0;
   obs += start;
   n_obs = end - start + 1;
   perturbers = always_included_perturbers;
   for( r = .85; r < 1.15; r += (r > .95 && r < 1.05 ? .001 : .002))
      for( direction = 0; direction < 2; direction++)
         {
         set_distance( obs, r);
         if( !find_parabolic_orbit( obs, n_obs, temp_orbit, direction))
            {
            ELEMENTS elem;
            double score;

            elem.gm = SOLAR_GM;
            calc_classical_elements( &elem, temp_orbit, obs[0].jd, 1);
            score = 10. * fabs( elem.ecc - 1.) + fabs( elem.incl * 180. / PI - 144.);
            if( best_score > score)
               {
               best_score = score;
               memcpy( orbit, temp_orbit, 6 * sizeof( double));
               }
            }
         }
   set_locs( orbit, obs[0].jd, obs, n_obs);
   return( obs[0].jd);
}

/* If there's a single observation,  generate a bogus orbit that satisfies
that observation,  putting it in a near-circular orbit with a=2.2,  with
as low an inclination as possible (z-velocity is zero).  Simply done to
make sure the program doesn't crash if only one valid observation is
supplied... this was easier than generating an error message right away.
Also used if we're just looking at the observations and don't actually
need a "real" orbit (force_bogus_orbit == true).      */

static void generate_bogus_orbit_for_single_obs( OBSERVE FAR *obs, double *orbit)
{
   unsigned i;
   double z;
   const double solar_r = 2.2;
   const double scale = GAUSS_K / sqrt( solar_r);
   const double dist_to_target =
            find_r_given_solar_r( obs, solar_r);

   obs->r = dist_to_target;
   for( i = 0; i < 3; i++)
      orbit[i] = obs->obs_posn[i] + obs->vect[i] * dist_to_target;
   z = sqrt( orbit[0] * orbit[0] + orbit[1] * orbit[1]);
   orbit[3] = -orbit[1] * scale / z;
   orbit[4] =  orbit[0] * scale / z;
   orbit[5] = 0.;
}

bool force_bogus_orbit = false;

static double only_one_position_available( OBSERVE FAR *obs,
                              const unsigned n_obs, double *orbit)
{
   unsigned i, n_usable_obs;
   double epoch = 0.;

   for( i = n_usable_obs = 0; i < n_obs; i++)
      if( !(obs[i].flags & OBS_DONT_USE))
         n_usable_obs++;
   if( n_usable_obs < 2 || force_bogus_orbit)
      {
      for( i = 0; i < n_obs && (obs[i].flags & OBS_DONT_USE); i++)
         ;
      if( i == n_obs)      /* all obs marked as unusable; */
         i = 0;            /* just use the first anyway */
      generate_bogus_orbit_for_single_obs( obs + i, orbit);
               /* Adjust epoch for light-time lag: */
      epoch = obs[i].jd - obs[i].r / AU_PER_DAY;
      set_locs( orbit, epoch, obs, n_obs);
      }
   return( epoch);
}

#define INITIAL_ORBIT_NOT_YET_FOUND       -2
#define INITIAL_ORBIT_FAILED              -1
#define INITIAL_ORBIT_FOUND                0

         /* Rather arbitrarily,  we don't use SR for spans greater */
         /* than 20 days.  Probably could drop that a lot without trouble. */
#define MAX_SR_SPAN 20

double *sr_orbits;
unsigned n_sr_orbits = 0;
unsigned max_n_sr_orbits;

double initial_orbit( OBSERVE FAR *obs, int n_obs, double *orbit)
{
   int i;
   int start = 0;
   bool dawn_based_observations = false;
   double arclen;
#ifdef CONSOLE
   char msg_buff[80];
#endif
   const double acceptable_score_limit = 5.;
   double best_score = 1e+50;
   double best_orbit[6];
   const int max_time = atoi( get_environment_ptr( "IOD_TIMEOUT"));

   for( i = 0; i < 6; i++)
      best_orbit[i] = 0.;
   if( debug_level)
      debug_printf( "initial_orbit(): %d obs;", n_obs);
   assert( orbit);
   for( i = 0; i < n_obs; i++)
      {
      obs[i].computed_ra  = obs[i].ra;
      obs[i].computed_dec = obs[i].dec;
      }

   while( n_obs && (!obs[n_obs - 1].is_included
                || (obs[n_obs - 1].flags & OBS_DONT_USE)))
      n_obs--;
   while( n_obs && (!obs->is_included || (obs->flags & OBS_DONT_USE)))
      {
      obs++;
      n_obs--;
      }
               /* We may have eliminated all observations,  or all */
               /* but one... in which case we do this :            */
   if( n_obs <= 1 || force_bogus_orbit)
      return( only_one_position_available( obs, n_obs, orbit));

   arclen = obs[n_obs - 1].jd - obs[0].jd;
   if( arclen > 730.)         /* two-year maximum */
      arclen = 730.;
                             /* following resets internals of find_orb */
   full_improvement( NULL, 0, NULL, 0., NULL, 0, 0.);
   for( i = 0; i < n_obs; i++)      /* solely to ensure a non-zero r */
      obs[i].r = 1.;
   for( i = 0; i < n_obs && is_solar_observatory( obs[i].mpc_code); i++)
      ;
   if( i == n_obs)      /* all obs are from SOHO or STEREOs */
      return( find_sungrazer_orbit( obs, n_obs, orbit));

   perturbers = AUTOMATIC_PERTURBERS | always_included_perturbers;
   if( !strcmp( obs->mpc_code, "Daw"))    /* For Dawn-based observations, */
      {                                   /* show a Ceres-centric orbit   */
      extern int forced_central_body;     /* and include asteroid perts   */

      forced_central_body = 100;
      dawn_based_observations = true;
      }
   if( debug_level)
      debug_printf( "  about to search orbits\n");

   writing_sr_elems = false;
   if( max_time)
      integration_timeout = clock( ) + (clock_t)( max_time * CLOCKS_PER_SEC);
   if( obs[n_obs - 1].jd - obs[0].jd < MAX_SR_SPAN)
      n_sr_orbits = get_sr_orbits( sr_orbits, obs, n_obs, 0, max_n_sr_orbits, .5, 0.);
   else
      n_sr_orbits = 0;     /* don't bother with SR for long time spans */
   writing_sr_elems = true;
   for( i = 0; (unsigned)i < n_sr_orbits && sr_orbits[i * 7 + 6] < .7; i++)
      ;
   n_sr_orbits = i;        /* cut orbits down to the "reasonable" ones */
   if( n_sr_orbits > 10)   /* got at least ten "reasonable" SR orbits; */
      {               /* accept the SR solution */
      const double epoch_shown = find_epoch_shown( obs, n_obs);

              /* SR orbits are stored in seven doubles,  including a score. */
      for( i = 1; i < (int)n_sr_orbits; i++)  /* Shift 'em to remove the score. */
         memmove( sr_orbits + i * 6, sr_orbits + i * 7, 6 * sizeof( double));
      memcpy( orbit, sr_orbits, 6 * sizeof( double));
      compute_sr_sigmas( sr_orbits, n_sr_orbits, obs[0].jd, epoch_shown);
      available_sigmas_hash = compute_available_sigmas_hash( obs, n_obs,
                  epoch_shown, perturbers, 0);
      set_locs( orbit, obs[0].jd, obs, n_obs);
      integration_timeout = 0;
      return( obs[0].jd);
      }

   while( best_score > acceptable_score_limit)
      {
      int end, n_subarc_obs, n_geocentric_obs = 0;
      const double max_arg_length_for_vaisala = 230.;
      double bogus_epoch;

      look_for_best_subarc( obs, n_obs, arclen, &start, &end);
      if( debug_level > 1)
         debug_printf( "  Current arc: %d to %d\n", start, end);
      bogus_epoch = only_one_position_available( obs + start, end - start + 1, orbit);
      if( bogus_epoch)
         {
         integration_timeout = 0;
         return( bogus_epoch);
         }
      for( i = 0; i < start; i++)
         obs[i].is_included = 0;
      for( i = start; i <= end; i++)
         {
         obs[i].is_included = 1;
         if( !strcmp( obs[i].mpc_code, "500"))
            n_geocentric_obs++;
         }
      exclude_unusable_observations( obs, n_obs);
      for( i = end + 1; i < n_obs; i++)
         obs[i].is_included = 0;
      arclen = obs[end].jd - obs[start].jd;
      if( debug_level)
         debug_printf( "From %f to %f (%f days)\n", obs[start].jd, obs[end].jd, arclen);
      n_subarc_obs = end - start + 1;
      fail_on_hitting_planet = true;
      if( n_subarc_obs >= 3)     /* at least three observations;  try Gauss */
         {
#ifdef CONSOLE
         if( show_runtime_messages)
            move_add_nstr( 14, 10, "In Gauss solution", -1);
#endif
         for( i = 0; i < 3; i++)
            {
            double epoch =
                convenient_gauss( obs + start, n_subarc_obs, orbit, 1., i);

            if( debug_level)
               debug_printf( "Gauss epoch: JD %f (%d)\n", epoch, i);
            if( epoch)
               {
               double score;

               if( !set_locs( orbit, epoch, obs + start, n_subarc_obs))
                  {
                  score = evaluate_initial_orbit( obs + start, n_subarc_obs, orbit);
                  if( debug_level > 2)
                     debug_printf( "Locations set; score %f (%d)\n", score, i);
                  if( score < 1000. &&
                           !integrate_orbit( orbit, epoch, obs[start].jd))
                     {
                     int pass;

                     for( pass = 0; pass < 2; pass++)
                        {
                        if( best_score > score)
                           {
                           best_score = score;
                           memcpy( best_orbit, orbit, 6 * sizeof( double));
                           if( debug_level > 2)
                              debug_printf( "A new winner from Gauss %d: %f\n", i, best_score);
                           }
                        if( !pass)
                           score = attempt_improvements( orbit, obs + start, n_subarc_obs);
                        }
                     }
                  }
               }
            else        /* break out of Gauss loop */
               i = 3;
            }
         }           /* end of trying Gauss */
#ifdef CONSOLE
      if( show_runtime_messages)
         move_add_nstr( 14, 10, "Gauss done", -1);
#endif
      if( arclen < max_arg_length_for_vaisala)
         for( i = 0; i < 2; i++)
            {
            int orbit_looks_reasonable = 1;
            double pseudo_r;

            if( i)          /* dist from observer (second) pass:  some ad hoc */
               {            /* code that says,  "for long arcs,  start farther */
                            /* from the observer".                             */
               pseudo_r = 0.004 * pow( arclen, .6666);
               if( dawn_based_observations)     /* for Dawn-based,  assume it */
                  pseudo_r = 1000. / AU_IN_KM;  /* may be a mere 1000 km away */
               if( n_geocentric_obs)            /* make sure we start outside the earth! */
                  pseudo_r += 6378. / AU_IN_KM;
               }
            else                  /* (first) Vaisala pass */
               pseudo_r = .1;

            while( pseudo_r < (i ? 5. : 200.) && orbit_looks_reasonable)
               {
               double pseudo_r_to_use;
               int herget_rval;
               double score;

               if( i)                          /* 2nd pass, dist from observer */
                  pseudo_r_to_use = pseudo_r;
               else                            /* 1st pass, dist from sun */
                  pseudo_r_to_use = -(1. + pseudo_r);
               herget_rval = herget_method( obs + start, n_subarc_obs,
                                    pseudo_r_to_use,
                                    pseudo_r_to_use,
                                    orbit, NULL, NULL, NULL);
               if( herget_rval < 0)    /* herget method failed */
                  score = 1.e+7;
               else if( herget_rval > 0)        /* vaisala method failed, */
                  score = 9e+5;                 /* but we should keep trying */
               else
                  {
                  adjust_herget_results( obs + start, n_subarc_obs, orbit);
                  score = evaluate_initial_orbit( obs + start, n_subarc_obs, orbit);
                  }
               if( debug_level > 2)
                  debug_printf( "%d, pseudo-r %f: score %f, herget rval %d\n",
                         i, pseudo_r, score, herget_rval);
               if( best_score > score)
                  {
                  if( debug_level > 2)
                     debug_printf( "A new winner\n");
                  best_score = score;
                  memcpy( best_orbit, orbit, 6 * sizeof( double));
                  }
#ifdef CONSOLE
               if( show_runtime_messages)
                  {
                  sprintf( msg_buff, "Method %d, r=%.4f", i, pseudo_r);
                  move_add_nstr( 14, 10, msg_buff, -1);
                  }
#endif
               if( score > 5e+4)   /* usually means eccentricity > 100! */
                  {
                  orbit_looks_reasonable = 0;      /* should stop looking */
                  if( debug_level > 2)
                     debug_printf( "%d: Flipped out at %f\n", i, pseudo_r);
                  }
               pseudo_r *= 1.2;
               }
            }
      if( best_score < 50. && n_obs > 2)
         {           /* maybe got a good orbit using Vaisala or Herget */
         double score;

         memcpy( orbit, best_orbit, 6 * sizeof( double));
         attempt_improvements( orbit, obs + start, n_subarc_obs);
         score = evaluate_initial_orbit( obs + start, n_subarc_obs, orbit);
         if( debug_level > 2)
            debug_printf( "After attempted improvements: score %f\n", score);
         if( score < best_score)       /* call it a success: */
            {
            memcpy( best_orbit, orbit, 6 * sizeof( double));
            best_score = score;
            }
         }
      arclen *= .7;     /* If we failed,  try again with a shorter arc */
      memcpy( orbit, best_orbit, 6 * sizeof( double));
      }

// set_locs( orbit, obs[start].jd, obs, n_obs);
   perturbers = perturbers_automatically_found & (~AUTOMATIC_PERTURBERS);
   perturbers |= always_included_perturbers;
   fail_on_hitting_planet = false;
   attempt_extensions( obs, n_obs, orbit, obs[start].jd);
// available_sigmas = NO_SIGMAS_AVAILABLE;
   integration_timeout = 0;
   if( *get_environment_ptr( "INCLUDE_ALL"))
      for( i = 0; i < n_obs; i++)
         obs[i].is_included = 1;
   return( obs[start].jd);    /* ...and return epoch = JD of first observation */
}

double generate_mc_variant_from_covariance( double *var_orbit,
                                                     const double *orbit);

int orbital_monte_carlo( const double *orbit, OBSERVE *obs, const int n_obs,
         const double curr_epoch, const double epoch_shown)
{
   unsigned i;
   extern int append_elements_to_element_file;
   extern const char *elements_filename;
   const char *saved_name = elements_filename;

   assert( sr_orbits);
   n_sr_orbits = max_n_sr_orbits;
   available_sigmas = NO_SIGMAS_AVAILABLE;
   elements_filename = "sr_elems.txt";
   for( i = 0; i < n_sr_orbits; i++)
      {
      double *torbit = sr_orbits + i * 6;
      const double sig_squared = generate_mc_variant_from_covariance( torbit, orbit);

      write_out_elements_to_file( torbit, curr_epoch, epoch_shown,
           obs, n_obs, "", 6, 1, ELEM_OUT_ALTERNATIVE_FORMAT);
      append_elements_to_element_file = 1;
      if( i < 1000)
         {
         double rms;
         int n_resids;

         set_locs( torbit, curr_epoch, obs, n_obs);
         rms = compute_weighted_rms( obs, n_obs, &n_resids);
         debug_printf( "Var %4d: %9.6f %.8f\n", i, sig_squared,
                        rms * rms * n_resids);
         }
      }
   set_locs( orbit, curr_epoch, obs, n_obs);
   compute_sr_sigmas( sr_orbits, n_sr_orbits, curr_epoch, epoch_shown);
   available_sigmas_hash = compute_available_sigmas_hash( obs, n_obs,
         epoch_shown, perturbers, 0);
   append_elements_to_element_file = 0;
   elements_filename = saved_name;
   return( 0);
}

static int count_observations_used( const OBSERVE *obs, int n_obs)
{
   int rval = 0;

   while( n_obs--)
      {
      if( obs->is_included)
         rval++;
      obs++;
      }
   return( rval);
}

void get_first_and_last_included_obs( const OBSERVE *obs,
              const int n_obs, int *first, int *last);      /* elem_out.c */

static double total_residual_err( const OBSERVE *obs)
{
   double xresid, yresid;

   get_residual_data( obs, &xresid, &yresid);
   return( sqrt( xresid * xresid + yresid * yresid));
}

/* extend_orbit_solution( ) is used if you have a solution covering part of
the arc of observations,  and wish to extend it to cover more observations.
To do this,  you might first look to see if the following or preceding
observations have errors of less than (say) a hundred arcseconds.  If so,
you should be able to toggle them on and do a full step or two to get them
to work correctly.  The first part of the code looks for such "simple" arc
extension.

   It'll also look for extensions that are no greater than the existing
arc length (i.e.,  if you have a 2.6-day arc,  it'll look 2.6 days before
the first and 2.6 days after the last observation) which have residuals
of less than two sigmas.  Such observations should be "safe" to add,  and
can help get past one or more bad observations.

   However,  you may not find such an easy extension,  either because there
are no more observations or because the adjacent observations have residuals
greater than the limit.  In that case,  the arc is extended by exactly one
observation.  If you have a "preceding" and "following" observation,  the
one that results in the least extension is chosen.

   You can set a 'time_limit' in days on how long the arc can be.  I added
this so that,  when an object is initially loaded,  some constraints can
be provided so that the code doesn't try to solve a hundred-year arc right
off the bat;  doing so can be slow,  persuading the user that the program
has locked up.

   Usually,  the number of added observations is returned (which can be zero
if both "ends" of the arc are already included,  or if any extension would
make the arc longer than time_limit).  -1 or -2 is returned if the arc was
extended by a single observation that was outside the limit. */

int extend_orbit_solution( OBSERVE FAR *obs, const int n_obs,
            const double limit, const double time_limit)
{
   int first_idx, last_idx, n_added = 0, initial_count;
   OBSERVE FAR *optr;
   double jd_low, jd_high;
   const double max_sigma = 2.;

   exclude_unusable_observations( obs, n_obs);
   initial_count = count_observations_used( obs, n_obs);
   get_first_and_last_included_obs( obs, n_obs, &first_idx, &last_idx);
   jd_low  = obs[first_idx].jd * 2. - obs[last_idx].jd;
   jd_high = obs[last_idx].jd * 2. - obs[first_idx].jd;
   optr = obs + first_idx - 1;
   while( first_idx > 0 && total_residual_err( optr) < limit
                && obs[last_idx].jd - optr->jd < time_limit)
      {
      optr->is_included = 1;
      optr--;
      first_idx--;
      n_added++;
      }
   while( first_idx > 0 && optr->jd > jd_low)
      {
      if( total_residual_err( optr) < max_sigma)
         {
         optr->is_included = 1;
         n_added++;
         }
      optr--;
      first_idx--;
      }
   optr = obs + last_idx + 1;
   while( last_idx < n_obs - 1 && total_residual_err( optr) < limit
                && optr->jd - obs[first_idx].jd < time_limit)
      {
      optr->is_included = 1;
      optr++;
      last_idx++;
      n_added++;
      }
   while( last_idx < n_obs - 1 && optr->jd < jd_high)
      {
      if( total_residual_err( optr) < max_sigma)
         {
         optr->is_included = 1;
         n_added++;
         }
      optr++;
      last_idx++;
      }

   if( !n_added)
      {
      int direction = 0;

      if( !first_idx)
         {
         if( last_idx < n_obs - 1)
            direction = 1;       /* already including first obs */
         }
      else if( last_idx == n_obs - 1)
         direction = -1;         /* already including last obs */
      else           /* could include either: */
         {
         const double time_before = obs[first_idx].jd - obs[first_idx - 1].jd;
         const double time_after  = obs[last_idx + 1].jd - obs[last_idx].jd;

         direction = (time_before > time_after ? 1 : -1);
         }
      if( direction == -1)
         if( obs[last_idx].jd - obs[first_idx - 1].jd < time_limit)
            {
            n_added = -1;
            obs[first_idx - 1].is_included = 1;
            }
      if( direction == 1)
         if( obs[last_idx + 1].jd - obs[first_idx].jd < time_limit)
            {
            n_added = -2;
            obs[last_idx + 1].is_included = 1;
            }
      }
   exclude_unusable_observations( obs, n_obs);
   if( initial_count == count_observations_used( obs, n_obs))
      n_added = 0;
// if( n_added)
//    available_sigmas = NO_SIGMAS_AVAILABLE;
   return( n_added);
}

void attempt_extensions( OBSERVE *obs, const int n_obs, double *orbit,
                                    const double epoch)
{
   double best_orbit[6];
   int best_start, best_end, i;
   bool done = false;
   const double residual_limit = 200.;   /* allow up to 200" in orbit extension */
   double arc_limit_in_days = atof( get_environment_ptr( "AUTO_ARC_LEN"));
   const int stored_setting_outside_of_arc = setting_outside_of_arc;
   int best_available_sigmas;
   unsigned best_perturbers = perturbers;

   if( !arc_limit_in_days)
      arc_limit_in_days = 3650;        /* Default to ten years at most */
   setting_outside_of_arc = 0;
         /* So far,  the "best" we've got is the orbit that was handed to */
         /* us by initial_orbit( ).  So let's store that :                */
   memcpy( best_orbit, orbit, 6 * sizeof( double));
   best_available_sigmas = available_sigmas;
   get_first_and_last_included_obs( obs, n_obs, &best_start, &best_end);
   do
      {
//    debug_printf( "Orbit: %f %f %f %f %f %f\n",
//          orbit[0], orbit[1], orbit[2],
//          orbit[3], orbit[4], orbit[5]);
      set_locs( orbit, epoch, obs, n_obs);
      if( !extend_orbit_solution( obs, n_obs, residual_limit, arc_limit_in_days))
         done = true;
      else           /* the arc was extended;  let's see if it works : */
         {
         double score[5];
         int result = 0, start, end;

         get_first_and_last_included_obs( obs, n_obs, &start, &end);
         if( debug_level)
            debug_printf( "   Try extend %d to %d (%f day arc)\n",
                     start, end, obs[end].jd - obs[start].jd);
//                printf( "   Try extend %d to %d (%f day arc)\n",
//                   start, end, obs[end].jd - obs[start].jd);
         if( obs[end].jd - obs[start].jd > 10.)
            perturbers |= (1 << 5);       /* add Jupiter for >10-day arc */
         if( obs[end].jd - obs[start].jd > 40.)
            perturbers |= (1 << 6);       /* add Saturn for >40-day arc */
         if( obs[end].jd - obs[start].jd > 200.)
            perturbers |= (1 << 2) | (1 << 3) | (1 << 4) | (1 << 7);
                    /* after 200 days,  venus,  earth,  mars,  */
                    /* and uranus can matter */
         if( obs[end].jd - obs[start].jd > 365. * 5.)
            perturbers |= (1 << 1) | (1 << 8);
                    /* after five years,  include Merc, Nept too */
         for( i = 0; i < 5 && !result; i++)
            {
            perturbers |= (perturbers_automatically_found & (~AUTOMATIC_PERTURBERS));
//          printf( "i = %d;  %d obs;  rms %f\n", i, n_obs,
//                       compute_rms( obs, n_obs));
            if( available_sigmas == COVARIANCE_AVAILABLE)
               {
               set_locs( orbit, epoch, obs, n_obs);
               improve_along_lov( orbit, epoch, eigenvects[0],
                                            n_extra_params + 6, n_obs, obs);
               }
            if( full_improvement( obs, n_obs, orbit, epoch, NULL,
                           NO_ORBIT_SIGMAS_REQUESTED, epoch))
               result = -1;   /* full improvement failed */
            else
               {
//             printf( "fully improved : rms %f\n",
//                       compute_rms( obs, n_obs));
               score[i] = evaluate_initial_orbit( obs, n_obs, orbit);
               if( i && score[i - 1] < score[i] + .1)      /* no real improvement... */
                  result = 1;                               /* we must have converged */
               if( score[i] > 10000.)      /* clearly blowing up */
                   result = -1;   /* counts as a a failure,  too */
               if( debug_level)
                  debug_printf( "     Iter %d: score %f; result %d\n",
                        i, score[i], result);
               }
            }
         if( result == -1)    /* failure of some sort;  stop trying to extend arc */
            done = true;
         if( result == 1)
            {
            if( score[i - 1] > 10.)        /* failed */
               done = true;
            else              /* looks like an improvement */
               {
               memcpy( best_orbit, orbit, 6 * sizeof( double));
               best_start = start;
               best_end = end;
               best_perturbers = perturbers;
               best_available_sigmas = available_sigmas;
               }
            }
         }
      }
      while( !done);
   setting_outside_of_arc = stored_setting_outside_of_arc;
   if( memcmp( orbit, best_orbit, 6 * sizeof( double)))
      {
      memcpy( orbit, best_orbit, 6 * sizeof( double));
      available_sigmas = best_available_sigmas;
      perturbers = best_perturbers;
#if 0
      full_improvement( obs, n_obs, orbit, epoch, NULL,
                           NO_ORBIT_SIGMAS_REQUESTED, epoch);
#endif
      }
   for( i = 0; i < best_start; i++)
      obs[i].is_included = 0;
   for( i = best_end + 1; i < n_obs; i++)
      obs[i].is_included = 0;
   set_locs( orbit, epoch, obs, n_obs);
}

int metropolis_search( OBSERVE *obs, const int n_obs, double *orbit,
               const double epoch, int n_iterations, double scale)
{
   int iter;
   double rms = compute_weighted_rms( obs, n_obs, NULL);
   extern double **eigenvects;
   double zorbit[9];

   if( !eigenvects)
      return( -1);

   memcpy( zorbit, orbit, 6 * sizeof( double));
   memcpy( zorbit + 6, solar_pressure, n_extra_params * sizeof( double));
   for( iter = 0; iter < n_iterations; iter++)
      {
      int i, j;
      double new_orbit[9], new_rms;

      memcpy( new_orbit, orbit, 9 * sizeof( double));
      for( i = 0; i < 6 + n_extra_params; i++)
//    for( i = 0; i < 1; i++)
         {
         const double n_sigmas = scale * gaussian_random( );

         for( j = 0; j < 6 + n_extra_params; j++)
            new_orbit[j] += n_sigmas * eigenvects[i][j];
         }
      memcpy( solar_pressure, new_orbit + 6, n_extra_params * sizeof( double));
      set_locs( new_orbit, epoch, obs, n_obs);
      new_rms = compute_weighted_rms( obs, n_obs, NULL);
      debug_printf( "Iter %d: prev rms %f; new rms %f\n", iter, rms, new_rms);
      if( new_rms < rms)
         {
         rms = new_rms;
         memcpy( zorbit, new_orbit, 9 * sizeof( double));
         scale /= .9;
         }
      else           /* step rejected */
         {
         memcpy( solar_pressure, zorbit + 6, n_extra_params * sizeof( double));
         scale *= .7;
         }
      }
   memcpy( orbit, zorbit, 6 * sizeof( double));
   memcpy( solar_pressure, zorbit + 6, n_extra_params * sizeof( double));
   set_locs( zorbit, epoch, obs, n_obs);
   return( 0);
}

#include "sigma.h"
#include "pl_cache.h"

void update_environ_dot_dat( void);     /* mpc_obs.cpp */
int galactic_confusion( const double ra, const double dec);
void pop_all_orbits( void);         /* orb_func2.cpp */
char *find_numbered_mp_info( const int number);    /* mpc_obs.cpp */

int clean_up_find_orb_memory( void)
{
   extern const char *temp_obs_filename;     /* miscell.cpp */
#ifndef _WIN32
   extern const char *lock_filename;      /* "/tmp/fo_lock" */
   extern FILE *lock_file;
#endif

   free_sigma_recs( );
   get_observer_data( NULL, NULL, NULL, NULL, NULL);
   get_object_name( NULL, NULL);
   planet_posn( -1, 0., NULL);
   add_gaussian_noise_to_obs( 0, NULL, 0.);
   full_improvement( NULL, 0, NULL, 0., NULL, 0, 0.);
   if( sr_orbits)
      {
      free( sr_orbits);
      sr_orbits = NULL;
      }
   find_objects_in_file( NULL, NULL, NULL);
   find_fcct_biases( 0., 0., 0, 0., NULL, NULL);
   get_find_orb_text( 0);
   load_cospar_file( NULL);
   update_environ_dot_dat( );
   load_earth_orientation_params( NULL, NULL);
   get_environment_ptr( NULL);
   pop_all_orbits( );
   galactic_confusion( -99., 0.);
   find_numbered_mp_info( 0);
#ifndef _WIN32
   fclose( lock_file);
   unlink( lock_filename);
   unlink( temp_obs_filename);
#else
   _unlink( temp_obs_filename);
#endif
   return( 0);
}
