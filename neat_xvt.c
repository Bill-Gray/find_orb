/* neat_xvt.c:  conversion of NEAT pointing logs to the CSS format

Copyright (C) 2017, Project Pluto

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
#include <assert.h>
#include <math.h>

/* Code to convert NEAT pointing logs,  as found at

https://sbnarchive.psi.edu/pds4/surveys/gbo.ast.neat.survey/data_tricam/tricam_archive_masterlist.tab
https://sbnarchive.psi.edu/pds4/surveys/gbo.ast.neat.survey/data_geodss/geodss_archive_masterlist.tab
https://sbnarchive.psi.edu/pds3/neat/tricam2/data/tricam_archive_masterlist_delivery2.tab

into the .csv-formatted type of pointing log CSS uses and which Find_Orb
can use for precovery searches,  and which list_gps can use to determine where
you may have found navigation satellites.  This is similar to 'sw_xvt.c'
and 'sw_xvt2.c' (q.v.,  should be distributed with this source),  which do
the same conversion for pre-2016 and post-2016 Spacewatch pointing logs.
And there will probably be other pointing log converters to come.

Compile with

gcc -Wall -O3 -o neat_xvt neat_xvt.c -lm

NOTE that dark images have RA = 99 99 99.999,  and that some images have RAs
slightly past 24 hours/360 degrees.  Which is why we allow RAs of 24 hours
and do an fmod( ) on the RA to put it in the 0 to 360 range.

ALSO NOTE that some files have extension ',fit' and others '.fit'.  I think
the former is an error.

The times in NEAT logs are for the _start_ of the image,  so we have to
add half the exposure time to get the mid-point of the exposure.  Lucky
us,  Haleakala and Palomar are in daylight at 0:00:00 UTC.   */

static void adjust_time( char *time_buff, const int half_exposure)
{
   unsigned sec = atoi( time_buff + 11) * 3600
               + atoi( time_buff + 14) * 60 + atoi( time_buff + 17);

   sec += half_exposure;
   assert( sec < 86400 + 86400);   /* allow for day-long exposure! */
   snprintf( time_buff + 11, 9, "%02u:%02u:%02u",
               sec / 3600, (sec / 60) % 60, sec % 60);
   time_buff[19] = ' ';
}

int main( const int argc, const char **argv)
{
   int i;

   printf( "# Generated from NEAT pointing logs by 'neat_xvt.c' (q.v.)\n");
#ifdef __TIMESTAMP__
   printf( "# Version %s\n", __TIMESTAMP__);
#else
   printf( "# Version %s %s\n", __DATE__, __TIME__);
#endif
   for( i = 1; i < argc; i++)
      {
      FILE *ifile = fopen( argv[i], "rb");
      double exposure = 0.;

      if( ifile)
         {
         char buff[200];

         printf( "# Fields from '%s'\n", argv[i]);
         while( fgets( buff, sizeof( buff), ifile))
            if( atoi( buff) <= 24)
               {
               const double ra = fmod( atof( buff + 26), 360.);
               const double dec = atof( buff + 38);
               const double new_exposure = atof( buff + 53);
               char *filename = buff + 61, *tptr;
               const char *mpc_code = NULL;

               if( !memcmp( buff + 132, "GEODSS", 6))
                  mpc_code = "566";       /* Haleakala-NEAT/GEODSS  */
               if( !memcmp( buff + 132, "Tri-Cam", 7))
                  mpc_code = "644";       /* Tri-Cam at Palomar  */
               assert( mpc_code);
               if( exposure != new_exposure)
                  {
                  exposure = new_exposure;
                  printf( "# Exposure: %.2f\n", exposure);
                  }
               while( *filename == ' ')
                  filename++;
               tptr = strstr( filename, ".fit");
               if( !tptr)
                  tptr = strstr( filename, ",fit");
               if( !tptr)
                  fprintf( stderr, "No .fit\n%s", buff);
               assert( tptr);
               *tptr = '\0';
               tptr = strchr( filename, '/');
               if( !tptr)
                  fprintf( stderr, "No /: \n%s", buff);
               *tptr = ',';
               adjust_time( buff + 107, (int)( exposure / 2.));
               printf( "%.3f,%.3f,%.19s,%3s,%s\n", ra, dec, buff + 107, mpc_code, filename);
               }
         fclose( ifile);
         }
      else
         fprintf( stderr, "Couldn't open %s\n", argv[i]);
      }
   return( 0);
}
