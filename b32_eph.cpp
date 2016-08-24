/* b32_eph.cpp: computes binary ephemerides,  mostly for use
with Guide

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
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "watdefs.h"
#include "lunar.h"
#include "afuncs.h"
#include "date.h"
#include "comets.h"

const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
int generic_message_box( const char *message, const char *box_type);
int earth_lunar_posn( const double jd, double FAR *earth_loc, double FAR *lunar_loc);
int planet_posn( const int planet_no, const double jd, double *vect_2000);
int integrate_orbit( double *orbit, double t0, double t1);
int add_ephemeris_details( FILE *ofile, const double start_jd,  /* b32_eph.c */
                                               const double end_jd);
int create_b32_ephemeris( const char *filename, const double epoch,
                const double *orbit, const int n_steps,         /* b32_eph.c */
                const double ephem_step, const double jd_start);
char *get_file_name( char *filename, const char *template_file_name);
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */

int add_ephemeris_details( FILE *ofile, const double start_jd,
                                               const double end_jd)
{
   time_t t0;
   char tbuff[128];
   FILE *ifile;
   extern const char *elements_filename;
   const char *vector_options = get_environment_ptr( "VECTOR_OPTS");

   t0 = time( NULL);
   fprintf( ofile, "\nCreated %s", ctime( &t0));

   full_ctime( tbuff, start_jd, CALENDAR_JULIAN_GREGORIAN);
   fprintf( ofile, "Ephemeris start: %s\n", tbuff);

   full_ctime( tbuff, end_jd, CALENDAR_JULIAN_GREGORIAN);
   fprintf( ofile, "Ephemeris end:   %s\n", tbuff);

   fprintf( ofile, "Times are all TDT\n");
   fprintf( ofile, "Positions/velocities are in %s J2000\n",
                      atoi( vector_options) ? "ecliptic" : "equatorial");

   ifile = fopen_ext( get_file_name( tbuff, elements_filename), "fcrb");
   while( fgets( tbuff, sizeof( tbuff), ifile))
      fwrite( tbuff, strlen( tbuff), 1, ofile);
   fclose( ifile);
   return( 0);
}

/* See 'interpol.cpp' for details on this.  Implementation is "clever",  and
therefore requires a fair bit of explanation. */

static double interpolate( const double *y, const double x, const int n_pts)
{
   double t = 1., c = 1., rval;
   int i;

   for( i = 0; i < n_pts; i++)
      {
      c *= x - (double)i;
      if( i)
         t *= -(double)i;
      }
   if( !c)        /* we're on an abscissa */
      rval = y[(int)( x + .5)];
   else
      {
      rval = y[0] / (t * x);
      for( i = 1; i < n_pts; i++)
         {
         t *= (double)i / (double)( i - n_pts);
         rval += y[i] / (t * (x - (double)i));
         }
      rval *= c;
      }
   return( rval);
}

#define MAX_N_PTS 20

static int best_interpolation_order( const double *ovals,
                  const int array_size, double *max_err)
{
   int i, axis, j, n_pts, rval = 0;

   *max_err = 1e+20;
   for( n_pts = 4; n_pts < MAX_N_PTS; n_pts += 2)
      {
      double worst = 0.;

      for( i = 0; (i + n_pts) * 2 < array_size; i++)
         for( axis = 0; axis < 3; axis++)
            {
            double tarray[MAX_N_PTS];
            double err;

            for( j = 0; j < n_pts; j++)
               tarray[j] = ovals[(i + j) * 6 + axis];

            if( n_pts & 1)    /* for odd number of points... */
               err = interpolate( tarray, (double)n_pts / 2., n_pts);
            else              /* for even number of points... */
               err = interpolate( tarray, (double)(n_pts - 1.) / 2., n_pts);
            err -= ovals[i * 6 + ((n_pts - 1) / 2) * 6 + 3 + axis];
            err = fabs( err);
            if( worst < err)
               worst = err;
            }
      if( worst < *max_err)
         {
         *max_err = worst;
         rval = n_pts;
         }
      }
   return( rval);
}

/* In the following,  we compute twice as many data points as we eventually
use,  going at half steps.  That lets us check the "intermediate" (odd)
points against the even ones,  to see how much error is involved;  the
order with the lowest error is then output.  We also can see which ordinate
has the highest absolute value,  and scale the long integers to fit. */

int create_b32_ephemeris( const char *filename, const double epoch,
                const double *orbit, const int n_steps,
                const double ephem_step, const double jd_start)
{
   double orbi[6], curr_jd, max_ordinate = 0., resolution;
   double *output_array = (double *)calloc( n_steps * 2,
                                 3 * sizeof( double));
   double prev_ephem_t = epoch, max_err;
   int i, j, jpl_id = 0, planet_center = 0;
   FILE *ofile;
   char tbuff[128];
   static const char *hdr_fmt = "%d %03d %10.1f %.2f %ld %d %g %d %d ";

                     /* hunt for the JPL ID: */
   for( i = 0; filename[i] && !jpl_id; i++)
      {
      jpl_id = atoi( filename + i);
      planet_center = (int)( filename[i] - '0');
      }
   memcpy( orbi, orbit, 6 * sizeof( double));
   curr_jd = jd_start;
   for( i = 0; i < n_steps * 2; i++)
      {
      double obs_posn[3];
      double *topo = output_array + i * 3;

//    ephemeris_t = curr_jd + td_minus_ut( curr_jd) / seconds_per_day;
      integrate_orbit( orbi, prev_ephem_t, curr_jd);
      prev_ephem_t = curr_jd;
      if( planet_center == 3)
         earth_lunar_posn( curr_jd, obs_posn, NULL);
      else
         planet_posn( planet_center, curr_jd, obs_posn);
      for( j = 0; j < 3; j++)
         topo[j] = orbi[j] - obs_posn[j];

      ecliptic_to_equatorial( topo);                           /* mpc_obs.cpp */
      for( j = 0; j < 3; j++)
         if( max_ordinate < fabs( topo[j]))
            max_ordinate = fabs( topo[j]);
      curr_jd += ephem_step / 2.;
      }

   resolution = max_ordinate / 2.e+9;
   ofile = fopen( filename, "wb");
   memset( tbuff, 0, 128);
   sprintf( tbuff, hdr_fmt, 128, jpl_id,
            jd_start, ephem_step, n_steps,
            best_interpolation_order( output_array, 2 * n_steps, &max_err),
            resolution, 32, 0);

   fwrite( tbuff, 128, 1, ofile);
   for( i = 0; i < n_steps; i++)
      for( j = 0; j < 3; j++)
         {
         const int32_t tval = (int32_t)( output_array[i * 6 + j] / resolution);
         fwrite( &tval, 1, sizeof( int32_t), ofile);
         }
   free( output_array);
   add_ephemeris_details( ofile, jd_start, curr_jd);
   fclose( ofile);
   sprintf( tbuff, "max err: %.3g\nEnd: ", max_err);
   full_ctime( tbuff + strlen( tbuff), curr_jd, CALENDAR_JULIAN_GREGORIAN);
   generic_message_box( tbuff, "o");
   return( 0);
}
