/* sw_xvt2.c: converter for "current format" Spacewatch logs

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
#include <assert.h>
#include <ctype.h>

/* Code to convert Spacewatch pointing logs from 2018 May 17 to the
present,  from a file sent 2023 July 21 by Cassandry Joly as

cleaned_0_9_meter_2018_2023.txt

Also see 'cssfield.c' for comments on the
format of the output from this program;  and see 'sw_xvt.c' and 'sw_xvt2.c'
for conversions from previous eras.  Compile with

gcc -Wall -Wextra -pedantic -Werror -O3 -o sw_xvt3 sw_xvt3.c     */

static double get_base_sixty( const char *buff)
{
   return( atof( buff) + atof( buff + 3) / 60. + atof( buff + 6) / 3600.);
}

/* Looks for nn:nn:nn sequence */

const char *find_base_60( const char *buff)
{
   const char *rval = strchr( buff + 3, ':');
   int i;

   if( rval)
      {
      rval -= 2;
      for( i = 0; rval && i < 8; i++)
         if( i == 2 || i == 5)
            {
            if( rval[i] != ':')
               rval = NULL;
            }
         else if( !isdigit( rval[i]))
            rval = NULL;
      }
   return( rval);
}

/* The Spacewatch logs give the time for the start of the exposure,  and
we really want the midpoint time.  Fortunately,  Spacewatch never observes
near 00:00:00 UTC,  so we can just add half the exposure time and leave the
day/month/year ignored. */

static void add_half_exposure( char *obuff, const char *itime,
                              const double exposure_time)
{
   const double seconds = get_base_sixty( itime) * 3600. + exposure_time / 2.;
   const int millisec = (int)( seconds * 1000. + .5);

   assert( millisec >= 0 && millisec < 86400000);
   snprintf( obuff, 14, "%02d:%02d:%02d.%03d",
                  millisec / 3600000, (millisec / 60000) % 60,  /* HH MM */
                  (millisec / 1000) % 60, millisec % 1000);    /* SS.sss */
}

#define INTENTIONALLY_UNUSED_PARAMETER( param) (void)(param)

int main( const int argc, const char **argv)
{
   const char *filename = "../temp/cleaned_0_9_meter_2018_2023.txt";
   FILE *ifile = fopen( filename, "rb");
   char buff[200];
   double curr_exposure = 0.;
   const double min_exposure = 0.1;
             /* accept all exposures,  including very short focussing ones */

   INTENTIONALLY_UNUSED_PARAMETER( argv);
   INTENTIONALLY_UNUSED_PARAMETER( argc);
   assert( ifile);
   if( !fgets( buff, sizeof( buff), ifile))
      {
      fprintf( stderr, "Didn't read header\n");
      return( -1);
      }
   printf( "# Spacewatch 'latest formula' logs,  processed with sw_xvt3.c (q.v.)\n");
#ifdef __TIMESTAMP__
   printf( "# Source file date %s\n", __TIMESTAMP__);
#else
   printf( "# Version %s %s\n", __DATE__, __TIME__);
#endif
   printf( "# Input file '%s'\n", filename);
   while( fgets( buff, sizeof( buff), ifile))
      {
      char filename[100], date[50], time_text[50];
      char exposure[40], ra_text[40], dec_text[40];
      int len;

      assert( !memcmp( buff, "Mosaic_Recovery_0.9m/20", 23));
      if( 6 == sscanf( buff, "%s %s %s %s %s %s%n", filename, date,
                 time_text, exposure, ra_text, dec_text, &len))
         {
         assert( buff[len] == 10 || buff[len] == 13);
         assert( *dec_text == '+' || *dec_text == '-');

         if( curr_exposure != atof( exposure))
            {
            curr_exposure = atof( exposure);
            printf( "# Exposure: %.0f s\n", curr_exposure);
            assert( curr_exposure > min_exposure);
            }
         printf( "%.4f,%c%.4f,%sT", get_base_sixty( ra_text) * 15, *dec_text,
                                get_base_sixty( dec_text + 1), date);
         add_half_exposure( buff, time_text, curr_exposure);
         printf( "%s,691,%s\n", buff, filename + 21);
         }
      else
         printf( "# Malformed %s", buff);
      }
   fclose( ifile);
   return( 0);
}

/*
   An example of text from a "new-formula" Spacewatch pointing log :

0.9-m Log: 2016 Apr 14 UT

01   K16EF8K   11:33:50.3   +01:34:08   Exp: 272 (x 3)  M. T. Read  Not Found
   Ephem Info - Mag(V): 21.1   Rates: -093.8, -031.3   Unc: 65  Type: NEO
Recovery Info - Mag(R): N/A     Rates: N/A, N/A   Mode: Sid  O-C: N/G
Focus: Unk
 Pass       Time      FWHM(")  EL     AZ      Xoff*  Yoff*   Rot     T2(C)
Pass 1   03:52:25.88  1.2     51.8   138.1  -008.5  +035.0   Unk      Unk
Pass 2   03:58:40.49  1.1     52.7   140.2  -008.9  +035.0   Unk      Unk
Pass 3   04:04:55.12  1.1     53.5   142.4  -009.6  +034.3   Unk      Unk

   This is telling me that on 2016 Apr 14,  three images were taken centered
at RA=11:33:50.3,  dec=+01:34:08 (in epoch of date),  at 03:52:25.88, 03:58:40.49,
and 04:04:55.12 UTC. The exposures were all 272 seconds long,  with seeing
given as FWHM.

   For field sizes and much more,  see comments in 'sw_xvt.c'.
*/
