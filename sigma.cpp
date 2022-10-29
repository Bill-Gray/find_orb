/* sigma.cpp: handle setting of sigmas for observations

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
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "watdefs.h"
#include "sigma.h"
#include "afuncs.h"
#include "date.h"

#define SIGMA_RECORD struct sigma_record

SIGMA_RECORD
   {
   double jd1, jd2, posn_sigma;
   double mag_sigma, time_sigma;
   int mag1, mag2;
   char mpc_code[5];
   char program_code;
   };

static int n_sigma_recs;
static SIGMA_RECORD *sigma_recs;

int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */

static int parse_sigma_record( SIGMA_RECORD *w, const char *buff)
{
   int rval = 0;

   if( *buff == ' ')
      {
      int i;

      memset( w, 0, sizeof( SIGMA_RECORD));
      memcpy( w->mpc_code, buff + 1, 3);
      w->mpc_code[3] = '\0';
      w->program_code = buff[5];
      w->jd1 = -1e+10;
      w->jd2 = 1e+10;
      w->mag1 = -100;
      w->mag2 = 3000;
      if( buff[32] == '.')
         w->mag1 = (int)floor( atof( buff + 29) * 10. + .5);
      if( buff[37] == '.')
         w->mag2 = (int)floor( atof( buff + 34) * 10. + .5);
      for( i = 0; i < 2; i++)
         if( buff[i * 11 + 8] != ' ')
            {
            const double jd = (double)dmy_to_day( atoi( buff + i * 11 + 16),
                                                  atoi( buff + i * 11 + 13),
                                                  atoi( buff + i * 11 + 8),
                                                  CALENDAR_JULIAN_GREGORIAN);

            if( i)
               w->jd2 = jd;
            else
               w->jd1 = jd;
            }
      if( buff[41] != ' ')
         w->posn_sigma = atof( buff + 40);
      if( buff[48] != ' ')
         w->mag_sigma = atof( buff + 45);
      else                       /* indicate 'no mag sigma set' */
         w->mag_sigma = 0.;
      if( buff[53] != ' ')
         w->time_sigma = atof( buff + 51) / seconds_per_day;
      else                       /* indicate 'no time sigma set' */
         w->time_sigma = 0.;
      rval = 1;
      }
   return( rval);
}

const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */

int load_up_sigma_records( const char *filename)
{
   const char *override_filename = get_environment_ptr( "SIGMAS_FILE");
   FILE *ifile;

   if( *override_filename)
      ifile = fopen_ext( override_filename, "frb");
   else
      ifile = fopen_ext( filename, "fcrb");
   if( ifile)
      {
      int i, j = 0;

      for( i = 0; i < 2; i++)
         {
         char buff[120];
         SIGMA_RECORD w;

         fseek( ifile, 0L, SEEK_SET);
         while( fgets( buff, sizeof( buff), ifile))
            if( parse_sigma_record( &w, buff))
               {
               if( !i)
                  n_sigma_recs++;
               else
                  sigma_recs[j++] = w;
               }
         if( !i)
            sigma_recs = (SIGMA_RECORD *)calloc( n_sigma_recs,
                                             sizeof( SIGMA_RECORD));
         if( !sigma_recs)
            debug_printf( "%d sigma recs not alloced\n", n_sigma_recs);
         }
      fclose( ifile);
      }
   assert( n_sigma_recs > 0);
   return( n_sigma_recs);
}

void free_sigma_recs( void)
{
   if( sigma_recs)
      {
      free( sigma_recs);
      sigma_recs = NULL;
      }
   n_sigma_recs = 0;
}

/* In determining the sigmas for an observation,  we start by setting
them all to zero,  i.e.,  "undetermined".  As we go through the sigma
records (see 'sigma.txt' for a description),  we may find a sigma
record that matches the program code,  time span,  and mag range of
this observation.

   When that happens,  we adopt any position,  mag,  or time sigma
from that observation.  Many records will only set the position sigma,
or the time or magnitude sigma.  So we have to keep going through the
mag records until all sigmas are set.  'sigma.txt' has to have a
catch-all record at the end,  currently set to say that if sigmas
haven't been figured out by that point,  the observation has a
positional sigma of half an arcsecond,  a mag residual of zero
(meaning "figure it out from the number of digits given for the
mag"),  and a time residual of five seconds.   */

double get_observation_sigma( const double jd, const int mag_in_tenths,
                  const char *mpc_code, double *mag_sigma,
                  double *time_sigma, const char program_code)
{
   int i;
   double position_sigma = 0.;
   extern double minimum_observation_jd;  /* default 1100 Jan 1 */
   extern double maximum_observation_jd;  /* default 2300 Jan 1 */

   if( mag_sigma)
      *mag_sigma = 0.;
   if( time_sigma)
      *time_sigma = 0.;
   assert( n_sigma_recs);
   assert( jd > minimum_observation_jd);
   assert( jd < maximum_observation_jd);
   for( i = 0; i < n_sigma_recs; i++)
      {
      SIGMA_RECORD *w = sigma_recs + i;

      if( !memcmp( mpc_code, w->mpc_code, 3)
                              || !memcmp( "   ", w->mpc_code, 3))
         if( jd > w->jd1 && jd < w->jd2)
            if( mag_in_tenths > w->mag1 && mag_in_tenths <= w->mag2)
               if( w->program_code == ' ' || w->program_code == program_code)
                  {
                  if( mag_sigma && !*mag_sigma)
                     *mag_sigma = w->mag_sigma;
                  if( time_sigma && !*time_sigma)
                     *time_sigma = w->time_sigma;
                  if( !position_sigma)
                     position_sigma = w->posn_sigma;
                  }
      }
                  /* At this point,  all sigmas _should_ be set to */
                  /* non-zero values.  Let's make sure of this :   */
   assert( position_sigma);
   assert( !mag_sigma || *mag_sigma);
   assert( !time_sigma || *time_sigma);
   return( position_sigma);
}
