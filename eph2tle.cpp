/* eph2tle.cpp: computes TLEs (Two-Line Elements) fitting numerically
integrated ephemerides of artificial satellites.  Executive summary
of what it does follows this GPL licence text.

Copyright (C) 2015-2017, Project Pluto

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
02110-1301, USA.

*****************************

   'eph2tle' assumes you've computed an orbit for your artsat
of interest using Find_Orb.  You've then computed an ephemeris
of geocentric J2000 equatorial state vectors with a small step
size (I almost always use 0.1 days).

   You can then run,  for example,

./eph2tle ephemeri.txt -o 13070b.tle

   to compute TLEs fitted to the ephemerides in 'ephemeri.txt'.

   By default,  the program reads in ten state vectors,  covering
a particular day,  and fits the TLE to those.  If the ephemeris
step size was,  say,  0.3 days,  and you used the '-f20' option
to tell eph2tle to read 20 state vectors at a time,  you'd get
TLEs that would each cover a six-day span.  Usually,  though,
TLEs covering longer time spans will have higher "worst errors".

   The program first fits a TLE to the middle state vector.  It
has a quick and dirty way of doing this (iterated_vector_to_tle)
which almost always converges nicely for objects with orbital
periods of less than two days.  Even if it doesn't converge,  it
provides a starting point for a downhill simplex search fitting
to all ten state vectors.  After that,  we try a least-squares
fit to get a "better" TLE.

   We then output the resulting TLE and move on to the next ten
state vectors.

   Convergence usually happens,  but it's not guaranteed.  For
one thing,  you can't fit a TLE to an hyperbolic orbit.  Fitting
a TLE to an object passing the Moon usually doesn't work.  There
are some other situations,  all involving high-flying objects,
that have horrible residuals.  These may reflect the inability
of the SDP4 model to fit every possible geocentric orbit (you
do get fits in such cases when forcing use of SGP4),  but I'm
not totally convinced of that;  it could be the program is not
doing a sufficiently exhaustive search in such cases. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "watdefs.h"
#include "date.h"
#include "norad.h"
#include "lsquare.h"
#include "stringex.h"
#include "afuncs.h"

#define AU_IN_KM 1.495978707e+8
#define AU_IN_METERS (AU_IN_KM * 1000.)
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078

#define FIT_DEFAULTS 0
#define FIT_BSTAR    1
#define FIT_EPOCH    2
#define FIT_BOTH     3

int vector_to_tle( tle_t *tle, const double *state_vect, const double epoch);

int verbose = 0;
int use_eight = 0, params_to_set = N_SAT_PARAMS;
int fitted = FIT_DEFAULTS;

#define EPHEM_TYPE_DEFAULT    '0'
#define EPHEM_TYPE_SGP        '1'
#define EPHEM_TYPE_SGP4       '2'
#define EPHEM_TYPE_SDP4       '3'
#define EPHEM_TYPE_SGP8       '4'
#define EPHEM_TYPE_SDP8       '5'
#define EPHEM_TYPE_HIGH       'H'

static int get_sxpx( const int ephem, const tle_t *tle, double *state,
                                 const double t_since_minutes)
{
   int rval = 0;

   if( tle->ephemeris_type != EPHEM_TYPE_HIGH &&
                ( tle->eo < 0. || tle->eo >= 1. || tle->xno < 0.))
      rval = -1;
   else
      {
      double params[N_SAT_PARAMS];
      int i, sxpx_rval;

      memset( params, 0, sizeof( params));
      if( ephem)
         {
         if( use_eight)
            {
            SDP8_init( params, tle);
            sxpx_rval = SDP8( t_since_minutes, tle, params, state, state + 3);
            }
         else
            {
            SDP4_init( params, tle);
            sxpx_rval = SDP4( t_since_minutes, tle, params, state, state + 3);
            }
         }
      else
         {
         SGP4_init( params, tle);
         sxpx_rval = SGP4( t_since_minutes, tle, params, state, state + 3);
         }
      if( sxpx_rval && verbose)
         {
         char buff[300];

         write_elements_in_tle_format( buff, tle);
         printf( "SXPX error: ephem %d, rval %d; e %f; tsince %f\n%s\n",
                         ephem, sxpx_rval, tle->eo, t_since_minutes, buff);
         }
      for( i = 0; i < 3; i++)
         state[i] /= AU_IN_KM;
      for( i = 3; i < 6; i++)
         state[i] *= minutes_per_day / AU_IN_KM;
      }
   return( rval);
}

#define MIN_DELTA_SQUARED 1e-22

/* A surprisingly decent way to get a TLE from a state vector is this:
Compute 'plain old Keplerian elements' from the state vector,  the
sort you would normally compute to model two-body motion,  as if you'd
never heard of TLEs or the SGP4/SDP4 orbital model.  Then,  use those
elements in a TLE and compute the corresponding state vector at epoch.

   The mismatch between two-body motion and the SGP4/SDP4 model means
that the result won't quite match the input.  However,  it'll (usually)
be fairly close,  and (usually) if you push the difference back into
the input state vector and iterate,  it will (usually) converge.

   Since it doesn't _always_ converge,  we keep track of the "best"
result (the one with the lowest root-mean-square difference from the
desired state vector).  That will usually be the last vector we compute,
but divergence happens.

   And,  of course,  the result is our best fit to the input state vector,
so we have something that may be a lovely fit to the position/velocity
at that particular epoch,  but which isn't at all good for any other
time.  Which is why the result is used only as the starting point for
a least-squares fit to the positions in an ephemeris covering the
time span of interest.  */

bool adjust_to_apogee = false;

static int iterated_vector_to_tle( tle_t *tle, const double *state_vect,
                           const double jd)
{
   int i, ephem = -1, iter = 0;
   double trial_state[6], delta = 1.;
   tle_t best_tle_yet;
   double best_delta_yet = 1e+20;
   double adjustment = 1.;
   const int max_iter = 70;
   int iterations_without_improvement = 0;

   memcpy( trial_state, state_vect, 6 * sizeof( double));
   while( iter++ < max_iter && iterations_without_improvement < 5)
      if( !vector_to_tle( tle, trial_state, jd))
         {
         double state_out[6];
         const double max_accepted_delta = .2;
         double scale = 1.;

         if( adjust_to_apogee)
            {
            if( tle->xmo > PI)
               tle->xmo -= PI + PI;
            if( tle->xmo > 0.)
               tle->epoch += (PI - tle->xmo) / (tle->xno * minutes_per_day);
            else
               tle->epoch -= (PI + tle->xmo) / (tle->xno * minutes_per_day);
            tle->xmo = PI;
            }
         if( iter < 4)
            ephem = 0;
         else
            ephem = select_ephemeris( tle);
         get_sxpx( ephem, tle, state_out, (jd - tle->epoch) * minutes_per_day);
#ifdef DEBUGGING_CODE
         printf( "%.5f %15.10f %15.10f %15.10f  %15.10f %15.10f %15.10f\n",
               jd, state_out[0], state_out[1], state_out[2],
               state_out[3], state_out[4], state_out[5]);
         printf( "%.4f  %15.10f %15.10f %15.10f  %15.10f %15.10f %15.10f\n",
               jd, state_out[0] - state_vect[0],
               state_out[1] - state_vect[1],
               state_out[2] - state_vect[2],
               state_out[3] - state_vect[3],
               state_out[4] - state_vect[4],
               state_out[5] - state_vect[5]);
#endif
         delta = 0.;
         for( i = 0; i < 6; i++)
            {
            state_out[i] -= state_vect[i];
            delta += state_out[i] * state_out[i];
            }
         if( delta > max_accepted_delta)
            scale = sqrt( max_accepted_delta / delta);
         for( i = 0; i < 6; i++)
            trial_state[i] -= state_out[i] * scale * adjustment;
         if( iter >= 4 && best_delta_yet > delta)
            {
            best_delta_yet = delta;
            best_tle_yet = *tle;
            iterations_without_improvement = 0;
            }
         else
            iterations_without_improvement++;
         if( verbose)
            printf( "Iteration %d worked : e = %f, t_per = %f, %g; ephem %d\n", iter,
                           tle->eo, 2. * PI / (tle->xno * minutes_per_day), delta * 1e+6, ephem);
         }
      else        /* Try slowing the object down in hopes of */
         {        /* getting a correct vector : */
         if( verbose)
            printf( "Iteration %d failed : e = %f, t_per = %f\n", iter,
                           tle->eo, 2. * PI / (tle->xno * minutes_per_day));
         memcpy( trial_state, state_vect, 6 * sizeof( double));
         adjustment *= .9;
         assert( iter > 2);
         }
   *tle = best_tle_yet;
   return( ephem);
}

static void error_exit( const int exit_value)
{
   printf( "Run as eph2tle <input filename> (options)\n\n");
   printf( "Options are:\n");
   printf( "   -i(international designator)     ex: -i97034A\n");
   printf( "   -n(NORAD designator)             ex: -n31415\n");
   printf( "   -v                               Verbose mode\n");
   printf( "   -o(output filename)\n");
   printf( "   -f(freq)                         Output freq (default = 10)\n");
   printf( "   -g                               Use SGP for all orbits,  never SDP\n");
   printf( "The input file is assumed to be an ephemeris of state vectors from Find_Orb.\n");
   exit( exit_value);
}

static char *fgets_trimmed( char *buff, const size_t buffsize, FILE *ifile)
{
   char *rval = fgets( buff, (int)buffsize, ifile);

   if( rval)
      {
      size_t i = strlen( rval);

      while( i && rval[i - 1] <= ' ')
         i--;
      rval[i] = '\0';
      }
   return( rval);
}

/* The six "parameters" to be set _can_ just be,  say,  inclination,  Omega,
omega, eccentricity,  semimajor axis,  and mean anomaly.  But there are
singularities in these for low inclinations and eccentricities.  To avoid
these,  we can make the six parameters the "equinoctial" elements

params[0] = h = e sin(lon_perihelion)
params[1] = k = e cos(lon_perihelion)
params[2]=  p = tan(incl/2) * sin(lon_asc_node)
params[3]=  p = tan(incl/2) * cos(lon_asc_node)
params[4] = mean longitude = omega + Omega + mean_anomaly
params[5] = semimajor axis

   Same orbit,  expressed in a manner that avoids singularities.  Except
that because TLEs can't handle hyperbolic orbits,  we do still have
singularities for negative a or e >= 1.  There's a risk that an iteration
step will take us to,  say,  e=1.2,  something SGP4/SDP4 won't understand.
So we really want any set of six real "params" to map to a valid TLE,
with e < 1.

   To eliminate those singularities as well,  we use the log of the mean
motion of a,  and revise h and k in a somewhat odd manner as well,  to
map eccentricities between 0 and 1 to the entire (h, k) plane :

params[0] = h = e sin(lon_perihelion) / (1-e)
params[1] = k = e cos(lon_perihelion) / (1-e)
params[5] = log( mean_motion)

   With these changes,  any valid TLE will map to a set of params,  and
any set of six real values put into params[] will map to a TLE.
*/

static void set_params_from_tle( double *params, const tle_t *tle)
{
   const double lon_perih = tle->omegao + tle->xnodeo;
   const double mean_lon = lon_perih + tle->xmo;
   const double r = tle->eo / (1. - tle->eo);
   const double tan_half_incl = tan( tle->xincl * .5);

   params[0] = r * sin( lon_perih);    /* (modified) h */
   params[1] = r * cos( lon_perih);    /* (modified) k */
   params[2] = tan_half_incl * sin( tle->xnodeo);     /* p */
   params[3] = tan_half_incl * cos( tle->xnodeo);     /* q */
   params[4] = mean_lon;
   params[5] = log( tle->xno);
   if( fitted & FIT_BSTAR)
      params[6] = tle->bstar;
   if( fitted & FIT_EPOCH)
      params[6 + (fitted & 1)] = tle->epoch;
}

static double zero_to_two_pi( double ival)
{
   ival = fmod( ival, 2. * PI);
   if( ival < 0.)
      ival += 2. * PI;
   return( ival);
}

static void set_tle_from_params( tle_t *tle, const double *params)
{
   const double r = sqrt( params[0] * params[0] + params[1] * params[1]);
   const double lon_perih = atan2( params[0], params[1]);
   const double tan_half_incl =
                    sqrt( params[2] * params[2] + params[3] * params[3]);

   tle->xincl  = 2 * atan( tan_half_incl);
   tle->xnodeo = atan2( params[2], params[3]);
   tle->eo     = r / (1. + r);
   tle->omegao = lon_perih - tle->xnodeo;
   tle->xmo    = params[4] - lon_perih;
   tle->xno    = exp( params[5]);
   tle->xmo = zero_to_two_pi( tle->xmo);
   tle->xnodeo = zero_to_two_pi( tle->xnodeo);
   tle->omegao = zero_to_two_pi( tle->omegao);
   if( fitted & FIT_BSTAR)
      tle->bstar = params[6];
   if( fitted & FIT_EPOCH)
      tle->epoch = params[6 + (fitted & 1)];
}

void init_simplex( double **vects, double *fvals,
         double (*f)( void *context, const double *vect),
               void *context, const int n);        /* simplex.c */
int simplex_step( double **vects, double *fvals,
         double (*f)( void *context, const double *vect),
               void *context, const int n);        /* simplex.c */

typedef struct
   {
   unsigned n_steps;
   int ephem;
   double step_size, ephem_start;
   const double *state_vect;
   double *tle_vect;
   tle_t *base_tle;
   } simplex_context_t;

static void evaluate_tle( const tle_t *tle, double *ivect,
            const double step_size, const double ephem_start,
            const int ephem, const unsigned n_steps,
            const double *ref)
{
   unsigned j;

   for( j = 0; j < n_steps; j++)
      get_sxpx( ephem, tle, ivect + j * 6,
               (ephem_start - tle->epoch) * minutes_per_day + (double)j * step_size);
   if( ref)
      for( j = 6 * n_steps; j; j--)
         *ivect++ -= *ref++;
}

static double simplex_scoring( void *icontext, const double *ivect)
{
   const simplex_context_t *context = (const simplex_context_t *)icontext;
   double err = 0.;
   size_t i, j;
   tle_t tle = *(context->base_tle);

   set_tle_from_params( &tle, ivect);
   evaluate_tle( &tle, context->tle_vect, context->step_size, context->ephem_start,
                  context->ephem, context->n_steps, context->state_vect);
   for( j = 0; j < context->n_steps; j++)
      {
      double *tptr = context->tle_vect + j * 6;

      for( i = 0; i < (context->n_steps > 1 ? 3 : 6); i++)
         err += tptr[i] * tptr[i];
      }
   err /= (double)context->n_steps;
   return( err);
}

#define MAX_PARAMS 10

static size_t max_simplex_iter = 3000;

int simplex_search( tle_t *tle, const double *starting_params,
                        const double *state_vect, const int ephem,
                        const unsigned n_steps, const double step_size,
                        const double ephem_start)
{
   double simp[MAX_PARAMS * MAX_PARAMS];
   double *vects[MAX_PARAMS], fvals[MAX_PARAMS];
   size_t i, iter;
   bool done = false;
   simplex_context_t context;
   size_t n_params = 6;

   if( fitted == FIT_BOTH)
      n_params = 8;
   else if( fitted)
      n_params = 7;
   for( i = 0; i < MAX_PARAMS; i++)
      vects[i] = simp + i * MAX_PARAMS;
   for( i = 0; i <= n_params; i++)
      {
      const double delta = .4;

      memcpy( vects[i], starting_params, n_params * sizeof( double));
      if( i)
         vects[i][i - 1] += delta;
      }
   context.n_steps = n_steps;
   context.ephem = ephem;
   context.step_size = step_size;
   context.state_vect = state_vect;
   context.tle_vect = (double *)calloc( 6 * n_steps, sizeof( double));
   context.base_tle = tle;
   context.ephem_start = ephem_start;
   init_simplex( vects, fvals, simplex_scoring, &context, (int)n_params);
   for( iter = 0; !done && iter < max_simplex_iter; iter++)
      {
      simplex_step( vects, fvals, simplex_scoring, &context, (int)n_params);
      if( fvals[n_params] / fvals[0] < 1.01 || fvals[0] < MIN_DELTA_SQUARED)
         done = true;
      }
   free( context.tle_vect);
   set_tle_from_params( tle, vects[0]);
   return( 0);
}

#ifndef PATH_MAX
   #define PATH_MAX 256
#endif

static FILE *fopen_from_findorb_dir( const char *filename, const char *permits)
{
   const char *path = getenv( "HOME");
   FILE *rval = NULL;

   if( path)
      {
      char fullname[PATH_MAX];

      assert( strlen( path) + strlen( filename) + 12 < PATH_MAX);
      strcpy( fullname, path);
      strcat( fullname, "/.find_orb/");
      strcat( fullname, filename);
      rval = fopen( fullname, permits);
      if( !rval)
         fprintf( stderr, "Couldn't open '%s'\n", fullname);
      }
   if( !rval)
      rval = fopen( filename, permits);
   return( rval);
}

static void auto_set_desigs( char *norad_desig, char *intl_desig)
{
   FILE *ifile = fopen( "/home/phred/tles/tle_list.txt", "rb");
   char buff[200];

   assert( ifile);
   while( fgets( buff, sizeof( buff), ifile))
      if( !memcmp( buff, "  note : next artsat will be ", 29))
         {
         fclose( ifile);
         memcpy( norad_desig, buff + 29, 5);
         memcpy( intl_desig, buff + 37, 8);
         printf( "Desigs set to '%s', '%s'\n", norad_desig, intl_desig);
         return;
         }
   assert( 0);       /* we shouldn't get here */
}


/* Certain objects have names (preliminary designations from surveys),  but
no NORAD or international designation.  The following ensures they get one,
which won't (we lightheartedly hope) conflict with other designations.  The
amateurs tracking classified objects use NORAD desigs starting at 90000,
so I'll use 89999 and count down... except for 9O0DC57;  the "classified"
tracking community already has designations for that.

   'sat_xref.txt' has further comments on this.     */

static void reset_desigs_by_name( const char *obj_name, tle_t *tle)
{
   FILE *ifile = fopen_from_findorb_dir( "sat_xref.txt", "rb");
   char buff[100];

   assert( ifile);
   while( fgets_trimmed( buff, sizeof( buff), ifile))
      if( *buff == ' ' && strlen( buff) > 21
                            && !strcmp( obj_name, buff + 21))
         {
         tle->norad_number = atoi( buff);
         memcpy( tle->intl_desig, buff + 12, 8);
         tle->intl_desig[8] = '\0';
         }
   fclose( ifile);
}

/* NOTE:  this precesses input J2000 state vectors to mean equator/ecliptic
of date.  I _think_ that's right,  but it's possible that nutation should be
included as well,  and even possible that SxPx assumes true orientation of
date:  i.e.,  the full set of earth orientation parameters,  including
polar motion and offsets from the IAU nutation theories,  ought to be used.

  The 'bulletin number' is in days since 2018 Jan 0,  or zero for TLEs
computed before then.  Hence,  BULLETIN_EPOCH = days from 1970 Jan 1
to 2018 Jan 0.  This will work until 2045 May 18.  At that point... well,
I'm not using the 'revolution number at epoch' field at present.  */

#define BULLETIN_EPOCH 17531
#define N_HIST_BINS 10

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen_from_findorb_dir( "eph2tle.txt", "rb");
   FILE *ofile = stdout;
   int i, j, output_freq = 10, line = 0;
   int tles_written = 0;
   int n_params = 6, n_iterations = 15;
   const int max_n_params = 8;
   char buff[200], obj_name[100];
   const char *default_intl_desig = "00000A";
   char intl_desig[9], norad_desig[6];
   double *slopes = (double *)calloc( max_n_params * 6, sizeof( double));
   double *vectors, worst_resid_in_run = 0., worst_mjd = 0.;
   double tdt = 0., *computed_vects;
   int ephem, progress_bar_freq = 2, ref_frame = -1;
   int epoch_index = -1;
   tle_t tle;
   const time_t curr_t = time( NULL);
   bool use_precession = true, archival = false;
   double step;
   unsigned n_steps, total_lines;
   int histo_counts[N_HIST_BINS];
   static int histo_divs[N_HIST_BINS] = { 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000 };
   double levenberg_marquardt_lambda0 = 0.;
   double sum_of_worst_resids = 0.;
   double dist_units = 1., time_units = 1.;
   const char *search_dist = NULL;
   long end_of_header_offset = 0L;

   if( argc < 2)
      error_exit( -1);
   strlcpy_error( intl_desig, default_intl_desig);
   strlcpy_error( norad_desig, "99999");
   setvbuf( stdout, NULL, _IONBF, 0);
   memset( &tle, 0, sizeof( tle_t));
   tle.classification = 'U';
   tle.ephemeris_type = EPHEM_TYPE_DEFAULT;
   tle.bulletin_number = (int)( curr_t / seconds_per_day - BULLETIN_EPOCH);
   *obj_name = '\0';
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'a': case 'A':
               adjust_to_apogee = true;
               break;
            case 'v': case 'V':
               verbose = 1 + atoi( argv[i] + 2);
               break;
            case 'b':
               fitted |= FIT_BSTAR;
               n_params++;
               break;
            case 'B':
               tle.bulletin_number = atoi( argv[i] + 2);
               break;
            case 'd':
               if( argv[i][2])
                  strcpy( obj_name, argv[i] + 2);
               else if( i < argc - 1)
                  strcpy( obj_name, argv[i + 1]);
               printf( "Object name reset to '%s'\n", obj_name);
               break;
            case 'e':
               fitted |= FIT_EPOCH;
               n_params++;
               break;
            case '8':
               use_eight = 1;
               break;
            case 'p': case 'P':
               params_to_set = atoi( argv[i] + 2);
               break;
            case 'o': case 'O':
               {
               const char *output_filename = argv[i] + 2;

               if( !*output_filename && i < argc - 1)
                  output_filename = argv[i + 1];
               ofile = fopen( output_filename, "wb");
               printf( "Output directed to %s\n", output_filename);
               if( !ofile)
                  {
                  perror( "Output not opened");
                  return( -1);
                  }
               setbuf( ofile, NULL);
               }
               break;
            case 'f': case 'F':
               output_freq = atoi( argv[i] + 2);
               break;
            case 'n': case 'N':
               strlcpy_error( norad_desig, argv[i] + 2);
               if( *norad_desig == 'a' && !norad_desig[1])
                  auto_set_desigs( norad_desig, intl_desig);
               break;
            case 'i': case 'I':
               strlcpy_error( intl_desig, argv[i] + 2);
               break;
            case 'l': case 'L':
               sscanf( argv[i] + 2, "%lf", &levenberg_marquardt_lambda0);
               break;
            case 'r':
               search_dist = argv[i] + 2;
               break;
            case 'R':
               srand( atoi( argv[i] + 2));
               break;
            case 'u': case 'U':
               archival = true;
               break;
            case 'x':
               epoch_index = atoi( argv[i] + 2);
               break;
            case 'y':
               max_simplex_iter = (size_t)atoi( argv[i] + 2);
               break;
            case 'z':
               n_iterations = atoi( argv[i] + 2);
               break;
            case 'g':
               tle.ephemeris_type = EPHEM_TYPE_SGP4;       /* force use of SGP4 */
               break;
            case 'h':
               tle.ephemeris_type = EPHEM_TYPE_HIGH;
               break;
            default:
               printf( "'%s' is not a valid command line option\n", argv[i]);
               error_exit( -2);
            }

   vectors = (double *)calloc( 12 * output_freq, sizeof( double));
   assert( vectors);
   computed_vects = vectors + 6 * output_freq;
   tle.norad_number = atoi( norad_desig);
   strcpy( tle.intl_desig, intl_desig);
   if( !ifile)
      {
      printf( "eph2tle.txt not found\n");
      error_exit( -4);
      }
   fprintf( ofile, "# Made by eph2tle, compiled " __DATE__ " " __TIME__ "\n");
   fprintf( ofile, "# Run at %s#\n", ctime( &curr_t));
   if( archival)
      fprintf( ofile, "# No updates     (archival TLEs)\n");
   if( search_dist)
      fprintf( ofile, "# Max error %s\n", search_dist);
   if( archival || search_dist)
      fprintf( ofile, "#\n");
   while( fgets_trimmed( buff, sizeof( buff), ifile))
      if( *buff != ';')
         fprintf( ofile, "%s\n", buff);
      else if( (!memcmp( buff, "; HEX", 5) && tle.ephemeris_type != EPHEM_TYPE_HIGH)
            || (!memcmp( buff, "; SGP", 5) && tle.ephemeris_type != EPHEM_TYPE_SGP4))
         {                 /* skip a section until we hit another comment */
         while( fgets_trimmed( buff, sizeof( buff), ifile))
            if( *buff == ';')
               break;
         }
   fclose( ifile);
   ifile = fopen( argv[1], "rb");
   if( !ifile)
      {
      printf( "%s not found\n", argv[1]);
      error_exit( -3);
      }
   if( fgets_trimmed( buff, sizeof( buff), ifile))
      {
      bool writing_data = false;
      double mjdt, mjdt_end;
      char *tptr = strstr( buff, "(500) Geocentric: ");

      while( *buff == ';')       /* skip leading comments,  if any */
         if( !fgets_trimmed( buff, sizeof( buff), ifile))
            {
            fprintf( stderr, "Nothing but comments in '%s'\n", argv[1]);
            return( -1);
            }
      if( sscanf( buff, "%lf %lf %u %d,%lf,%lf\n", &tdt, &step, &total_lines,
                  &ref_frame, &dist_units, &time_units) != 6)
         {
         fprintf( stderr, "Bad header in '%s'\n", argv[1]);
         fprintf( stderr, "Header was:\n%s\n", buff);
         return( -1);
         }
      end_of_header_offset = ftell( ifile);
      if( tdt > 3e+7)         /* probably in seconds after J2000, */
         tdt = 2451545. + tdt / seconds_per_day;  /* a SPK convention */
      if( ref_frame == -1)
         {                             /* input coords are already mean */
         use_precession = false;       /* of date;  don't precess 'em */
         ref_frame = 0;
         }
      mjdt = tdt - 2400000.5;
      mjdt_end = mjdt + step * (double)( total_lines - total_lines % output_freq);
      if( tptr && !*obj_name)
         strcpy( obj_name, tptr + 18);
      fprintf( ofile, "# Ephem range: %f %f %f\n",
            mjdt, mjdt_end, step * (double)output_freq);
      while( fgets_trimmed( buff, sizeof( buff), ifile))
         {
         if( !memcmp( buff, "Created ", 8))
            writing_data = true;
         if( writing_data && *buff != '#')
            fprintf( ofile, (*buff ? "# %s\n" : "#%s\n"), buff);
         if( !memcmp( buff, "Orbital elements: ", 18) && !*obj_name)
            strcpy( obj_name, buff + 19);
         if( !memcmp( buff, "Ephemeris", 9))
            printf( "%s\n", buff);
         }
      if( *obj_name)
         {
         printf( "Object: %s\n", obj_name);
         if( tle.norad_number == 99999)
            {
            tptr = strstr( obj_name, "NORAD ");
            if( tptr)
               tle.norad_number = atoi( tptr + 6);
            }
         if( !strcmp( intl_desig, default_intl_desig))
            for( tptr = obj_name; *tptr; tptr++)
               if( atoi( tptr) > 1900 && tptr[4] == '-' &&
                     atoi( tptr + 5) > 0)
                  {
                  memset( tle.intl_desig, 0, sizeof( tle.intl_desig));
                  memcpy( tle.intl_desig, tptr + 2, 2);    /* get year */
                  memcpy( tle.intl_desig + 2, tptr + 5, 4); /* launch # */
                  if( isalpha( tptr[9]))    /* two- or three-letter suffix */
                     {
                     tle.intl_desig[6] = tptr[9];
                     if( isalpha( tptr[10]))    /* three-letter suffix */
                        tle.intl_desig[7] = tptr[10];
                     }
                  }
         if( intl_desig == default_intl_desig)
            reset_desigs_by_name( obj_name, &tle);
         }
      }
   if( !strcmp( tle.intl_desig, default_intl_desig))
      fprintf( stderr, "WARNING: International designation left at default!\n");
   for( i = 0; i < N_HIST_BINS; i++)
      histo_counts[i] = 0;
   if( tle.ephemeris_type == EPHEM_TYPE_SGP4)
      {
      fprintf( ofile, "# SGP4 only: these TLEs are _not_ fitted to SDP4,  even for\n");
      fprintf( ofile, "# deep-space TLEs.  These may not work with your software.\n");
      }
   fprintf( ofile, "#\n");
   fprintf( ofile, "# 1 NoradU COSPAR   Epoch.epoch     dn/dt/2  d2n/dt2/6 BSTAR    T El# C\n");
   fprintf( ofile, "# 2 NoradU Inclina RAAscNode Eccent  ArgPeri MeanAno  MeanMotion Rev# C\n");
   fseek( ifile, end_of_header_offset, SEEK_SET);
   n_steps = total_lines / output_freq;
   while( n_steps--)
      {
      double *sptr = vectors;
      const double jan_1956 = 2435473.5, jan_2050 = 2469807.5;
      const double jd_utc = tdt - td_minus_utc( tdt) / seconds_per_day;

      for( i = 0; i < output_freq && fgets( buff, sizeof( buff), ifile);
                                                      i++, sptr += 6)
         {
         double jdt, utc;
         double precession_matrix[9], ivect[6];

         if( sscanf( buff, "%lf%lf%lf%lf%lf%lf%lf", &jdt,
                        ivect, ivect + 1, ivect + 2,
                        ivect + 3, ivect + 4, ivect + 5) != 7)
            {
            fprintf( stderr, "Error reading input ephem:\n%s\n", buff);
            return( -2);
            }
         if( jdt > 3e+7)         /* probably in seconds after J2000, */
            jdt = 2451545. + jdt / seconds_per_day;  /* a SPK convention */
         if( jdt < jan_1956 || jdt > jan_2050)
            {
            fprintf( stderr, "JDT %f is outside the valid range\n", jdt);
            return( -3);
            }
         utc = jdt - td_minus_utc( jdt) / seconds_per_day;
         for( j = 0; j < 6; j++)
            ivect[j] /= dist_units;
         for( j = 3; j < 6; j++)
            ivect[j] *= time_units;
         if( ref_frame == 1)
            {
            ecliptic_to_equatorial( ivect);
            ecliptic_to_equatorial( ivect + 3);
            }
         if( use_precession)
            setup_precession( precession_matrix, 2000.,
                                      2000. + (utc - 2451545.) / 365.25);
         else
            set_identity_matrix( precession_matrix);
         precess_vector( precession_matrix, ivect, sptr);
         precess_vector( precession_matrix, ivect + 3, sptr + 3);
         }
      assert( i == output_freq);

      tle.epoch = jd_utc;
      if( tle.ephemeris_type == EPHEM_TYPE_HIGH)
         {
         double *svect = &tle.xincl;

         ephem = 1;
         for( i = 0; i < 3; i++)
            {
            svect[i] = vectors[i] * AU_IN_METERS;
            svect[i + 3] = vectors[i + 3] * AU_IN_METERS / seconds_per_day;
            }
         }
      else
         {
         double start_params[max_n_params];

         if( epoch_index == -1)
            epoch_index = output_freq / 2;    /* use middle vector;  it improves the */
         tle.epoch += (double)epoch_index * step;      /* convergence for simplex */
         ephem = iterated_vector_to_tle( &tle, vectors + epoch_index * 6, tle.epoch);
         if( ephem != -1)
            ephem = select_ephemeris( &tle);
         if( verbose)
            printf( "   ephem selected = %d\n", ephem);
         if( ephem != -1 && tle.ephemeris_type == EPHEM_TYPE_SGP4)
            ephem = 0;
         set_params_from_tle( start_params, &tle);
         simplex_search( &tle, start_params, vectors, ephem,
                     output_freq, step * minutes_per_day, jd_utc);
         }


      int lsquare_rval, use_damping = 1, iter;
      char obuff[200];
      double worst_resid = 1e+20;
      tle_t tle_to_output = tle;

      if( verbose)
         printf( "   least-square fitting\n");
//    if( ephem == -1)     /* failed from the get-go */
//       failure = -1;
      for( iter = 0; iter < n_iterations; iter++)
         {
         void *lsquare = lsquare_init( n_params);
         double state0[6], params[max_n_params];
         double differences[max_n_params], rms_change = 0.;
         double this_worst_resid = 0.;
         extern double levenberg_marquardt_lambda;

         evaluate_tle( &tle, computed_vects, step * minutes_per_day, jd_utc,
                     ephem, output_freq, vectors);
         for( i = 0; i < output_freq * 6; i += 6)
            for( j = 0; j < 3; j++)
               if( this_worst_resid < fabs( computed_vects[i + j]))
                  this_worst_resid = fabs( computed_vects[i + j]);
         this_worst_resid *= AU_IN_KM;
         if( this_worst_resid < worst_resid)
            {        /* improvement,  or at least minor worsening */
            worst_resid = this_worst_resid;
            tle_to_output = tle;
            levenberg_marquardt_lambda = 0.;
            }
         else     /* not doing well:  let's try dampened iterations */
            levenberg_marquardt_lambda += levenberg_marquardt_lambda0;
         if( use_damping)
            levenberg_marquardt_lambda += levenberg_marquardt_lambda0;
         if( verbose)
            {
            write_elements_in_tle_format( obuff, &tle);
            printf( "Iter %d: worst resid %f\n%s\n", iter, this_worst_resid, obuff);
            }
         set_params_from_tle( params, &tle);
         for( j = 0; j < output_freq; j++)
            {
            double resid2 = 0.;
//          const double time_diff_in_minutes = (double)j
//                                     * step * minutes_per_day;
            const double time_diff_in_minutes = (double)(j - epoch_index)
                                       * step * minutes_per_day;

            for( i = 0; i < n_params; i++)
               {
               double state1[6], state2[6];
               double delta = (i == 6 ? 1.e-5 : 1.e-4);
               int k;

               if( tle.ephemeris_type == EPHEM_TYPE_HIGH)
                  delta = (i >= 3 ? 1e-4 : 1.);    /* one meter or 10^-4 m/s */
               params[i] -= delta;
               set_tle_from_params( &tle, params);
               get_sxpx( ephem, &tle, state1, time_diff_in_minutes);
               params[i] += delta + delta;
               set_tle_from_params( &tle, params);
               get_sxpx( ephem, &tle, state2, time_diff_in_minutes);
               params[i] -= delta;
               set_tle_from_params( &tle, params);
               for( k = 0; k < 6; k++)
                  slopes[k * n_params + i] = (state2[k] - state1[k]) / (2. * delta);
               if( verbose > 2)
                  {
                  for( k = 0; k < 6; k++)
                     printf( "%10.3g ", slopes[k * n_params + i]);
                  printf( "\n");
                  }
               }
            get_sxpx( ephem, &tle, state0, time_diff_in_minutes);
            if( verbose > 1)
               printf( "JD %f: ", jd_utc);
            for( i = 0; i < 3; i++)
               {
               const double residual = vectors[j * 6 + i] - state0[i];

               if( verbose == 2)
                  printf( "%f ", residual * AU_IN_KM);
               if( verbose == 3)
                  printf( "   %f (%f %f)\n", residual * AU_IN_KM,
                              vectors[j * 6 + i] * AU_IN_KM,
                              state0[i] * AU_IN_KM);
               resid2 += residual * residual;
               lsquare_add_observation( lsquare, residual,
                        1., slopes + i * n_params);
               }
            rms_change += resid2;
            if( resid2 > this_worst_resid)
               this_worst_resid = resid2;
            if( verbose > 1)
               printf( "\n");
            }

         lsquare_rval = lsquare_solve( lsquare, differences);
         lsquare_free( lsquare);
         use_damping = 0;
         if( lsquare_rval)
            {
            if( tle.ephemeris_type != EPHEM_TYPE_HIGH)
               {
               char date_string[30];

               full_ctime( date_string, tdt,
                    FULL_CTIME_YMD | FULL_CTIME_DATE_ONLY);
               fprintf( stderr, "ERROR %d in lsquare soln: MJD %f = %s\n",
                           lsquare_rval, tdt - 2400000.5, date_string);
               }
            use_damping = 1;
            }
         else if( tle.ephemeris_type == EPHEM_TYPE_HIGH)
            {
            for( i = 0; i < n_params; i++)
               params[i] += differences[i];
            set_tle_from_params( &tle, params);
            }
         else
            {
            rms_change = 0.;
            for( i = 0; i < 6; i++)
               rms_change += differences[i] * differences[i];
            rms_change = sqrt( rms_change);
            for( i = 0; i < n_params; i++)
               params[i] += differences[i];
            set_tle_from_params( &tle, params);
            if( verbose)
               printf( "  change in TLE = %f\n", rms_change);
            }
         }

      full_ctime( buff, tdt,
                    FULL_CTIME_YMD | FULL_CTIME_FORMAT_HH_MM);
//    if( !failure)
         {
         const double revs_per_day = tle_to_output.xno * minutes_per_day / (2. * PI);

//       if( tle.ephemeris_type != EPHEM_TYPE_HIGH)
            fprintf( ofile, "\n# Worst residual: %.2f km\n",
                          worst_resid);
//       else
//          fprintf( ofile, "\n");

         assert( tle.ephemeris_type == EPHEM_TYPE_HIGH
                  || revs_per_day < 20.);  /* allows some margin for suborbital TLEs */
         write_elements_in_tle_format( obuff, &tle_to_output);
         if( verbose)
            {
            double params[N_SAT_PARAMS], state[6];

            SDP4_init( params, &tle_to_output);
            SDP4( 0, &tle_to_output, params, state, state + 3);
            printf( "   Node: %f\n", params[25] * 180. / PI);
            printf( "   xinc: %f\n", params[27] * 180. / PI);
            printf( "   em:   %f\n", params[26]);
            printf( "%s", obuff);
            }
         fprintf( ofile, "# MJD %f (%s)\n", tdt - 2400000.5, buff);
         if( *obj_name)
            fprintf( ofile, "%s\n", obj_name);
         obuff[69] = obuff[140] = '\0';
         fprintf( ofile, "%s\n%s\n", obuff, obuff + 71);
         sum_of_worst_resids += worst_resid;
         if( worst_resid_in_run < worst_resid)
            {
            worst_resid_in_run = worst_resid;
            worst_mjd = tdt - 2400000.5;
            }
         i = 0;
         while( i < N_HIST_BINS - 1 && worst_resid > (double)histo_divs[i])
            i++;
         histo_counts[i]++;
         }
//    else
//       fprintf( ofile, "FAILED (%d) for JD %.2f = %s\n", failure,
//                      jd_utc, buff);
      tles_written++;
      line++;
      if( ofile != stdout && !(line % progress_bar_freq))
         {
         static clock_t t0;
         const clock_t t1 = clock( );

         printf( "Line %d of %u (%u%% done): %d written, JD %f   \r",
                  line, total_lines, line * 100 * output_freq / total_lines,
                  tles_written, tdt);
         fflush( stdout);
         if( t1 - t0 < CLOCKS_PER_SEC / 2)
            progress_bar_freq <<= 1;
         else if( progress_bar_freq > 1)
            progress_bar_freq >>= 1;
         t0 = t1;
         }
      tdt += step * (double)output_freq;
      }
   if( ifile)
      fclose( ifile);
   while( 1)
      {
      if( ofile != stdout)
         printf( "\n");
      fprintf( ofile, "Avg worst resid: %.2f km\n",
                                   sum_of_worst_resids / (double)tles_written);
      fprintf( ofile, "Worst residual in entire run: %.2f km on MJD %.1f\n",
                                   worst_resid_in_run, worst_mjd);
      fprintf( ofile, "       1     3     10    30    100   300"
                      "   1K    3K    10K   km\n");
      for( i = 0; i < N_HIST_BINS; i++)
         fprintf( ofile, "%6d", histo_counts[i]);
      fprintf( ofile, "\n");
      if( ofile != stdout)
         {
         fclose( ofile);
         ofile = stdout;
         }
      else
         break;
      }
   printf( "Freeing vectors\n");
   free( vectors);
   free( slopes);
   printf( "All done\n");
   return( 0);
}

