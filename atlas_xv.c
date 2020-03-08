/* atlas_xv.c:  conversion of ATLAS pointing logs to the CSS format

Copyright (C) 2019, Project Pluto

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

/* Code to convert ATLAS pointing logs,  as found at

http://astroportal.ifa.hawaii.edu/atlas/pluto/T05.58484.log
http://astroportal.ifa.hawaii.edu/atlas/pluto/T08.58484.log

into the .csv-formatted type of pointing log CSS uses and which Find_Orb
can use for precovery searches,  and which list_gps can use to determine where
you may have found navigation satellites.  This is similar to 'sw_xvt.c'
and 'sw_xvt2.c' (q.v.,  should be distributed with this source),  which do
the same conversion for pre-2016 and post-2016 Spacewatch pointing logs;
'neat_xvt.c' for NEAT images;  and 'j95_xvt.c' for (J95) pointing logs.
And there will probably be other pointing log converters to come.

Compile with

gcc -Wextra -pedantic -Wall -O3 -o atlas_xv atlas_xv.c -lm

*/

int main( const int argc, const char **argv)
{
   int i;

   printf( "# Generated from ATLAS pointing logs by 'atlas_xv.c' (q.v.)\n");
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
         double mjd;

         printf( "# Fields from '%s'\n", argv[i]);
         while( fgets( buff, sizeof( buff), ifile))
            if( ( mjd = atof( buff + 15)) > 50000.)
               {
               const double ra = atof( buff + 88);
               const double dec = atof( buff + 97);

               if( ra && dec)
                  {
                  const double new_exposure = atof( buff + 37);
                  const double seconds_per_day = 24. * 60. * 60.;
                  char mpc_code[4];

                  assert( ra > 0. && ra < 360.);
                  assert( dec > -90. && dec < 90.);
                  if( !memcmp( buff, "01a", 3))
                     strcpy( mpc_code, "T08");
                  if( !memcmp( buff, "02a", 3))
                     strcpy( mpc_code, "T05");
                  else
                     assert( 1);
                  if( exposure != new_exposure)
                     {
                     exposure = new_exposure;
                     printf( "# Exposure: %.2f\n", exposure);
                     }
                  mjd += exposure * .5 / seconds_per_day;
                  printf( "%07.3f,%+07.3f,mjd%.5f,%.3s,%.14s\n", ra, dec, mjd,
                                 mpc_code, buff);
                  }
               }
         fclose( ifile);
         }
      else
         fprintf( stderr, "Couldn't open %s\n", argv[i]);
      }
   return( 0);
}
