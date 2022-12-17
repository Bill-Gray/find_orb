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

/* Code to convert "current format" Spacewatch pointing logs into the
form expected by 'cssfield.c'.   Note comments at bottom of this file
on the format of that log. Also see 'cssfield.c' for comments on the
format of the output from this program.  Compile with

gcc -Wall -O3 -o sw_xvt2 sw_xvt2.c

   We read lines from the SW log.  If it starts with

0.9-m Log: yyyy mmm dd UT

   we reset the date (converted to ISO form).  If it's got an RA
and a dec in it,  we look after those for times of observation,  of
which there can be up to three.  For each of those times of observation,
we output a line with the RA and dec in decimal degrees,  the ISO
formatted observation time,  the MPC code 691,  and 'na' (the filename
is Not Applicable here.)

   That puts everything in the 'standard' CSV format inhaled by
'cssfield',  dumping to stdout.

   It assumes the log file is

Spacewatch_0.9m.log.txt

   You can specify a different name on the command line.       */

static void get_iso_date( char *iso_date, const char *buff)
{
   char month[6];
   const char *tptr, *months = "JanFebMarAprMayJunJulAugSepOctNovDec";
   int day, year, m;

   printf( "# %s", buff);
   if( sscanf( buff, "%d %4s %d", &year, month, &day) != 3)
      {
      printf( "Malformed date : '%s'\n", buff);
      exit( -1);
      }
   tptr = strstr( months, month);
   assert( tptr);
   m = tptr - months;
   assert( m % 3 == 0);
   m = m / 3 + 1;
   assert( m > 0 && m < 13);
   assert( year >= 2016 && year < 2050);
   snprintf( iso_date, 20, "%4d-%02d-%02d", year, m, day);
}

   /* convert,  say,  '31:41:59' to 31 + 41/60 + 59/3600 */

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

   snprintf( obuff, 16, "%02d:%02d:%02d.%03d",
                  millisec / 3600000, (millisec / 60000) % 60,  /* HH MM */
                  (millisec / 1000) % 60, millisec % 1000);    /* SS.sss */
}

int main( const int argc, const char **argv)
{
   const char *filename = (argc == 2 ? argv[1] :
                     "Spacewatch_0.9m.log.txt");
   FILE *ifile = fopen( filename, "rb");
   char buff[200];
   char iso_date[30];
   double exposure = 0.;
   double ra, dec;
   const double min_exposure = 5.;        /* skip focussing exposures */

   assert( ifile);
   *iso_date = '\0';
   printf( "# Spacewatch 'new formula' logs,  processed with sw_xvt2.c (q.v.)\n");
#ifdef __TIMESTAMP__
   printf( "# Source file date %s\n", __TIMESTAMP__);
#else
   printf( "# Version %s %s\n", __DATE__, __TIME__);
#endif
   printf( "# Input file '%s'\n", filename);
   while( fgets( buff, sizeof( buff), ifile))
      {
      const char *tptr = strstr( buff, "Exp: ");

      if( !memcmp( buff, "Observer", 8))
         printf( "# %s", buff);
      if( !memcmp( buff, "0.9-m Log: ", 11))
         get_iso_date( iso_date, buff + 11);
      if( tptr)
         {
         exposure = atof( tptr + 5);
         if( exposure > min_exposure)
            printf( "# Exposure: %.0f s\n", exposure);
         tptr = find_base_60( buff);
         assert( tptr);
         ra = get_base_sixty( tptr) * 15.;
         tptr = find_base_60( tptr + 8);
         assert( tptr);
         dec = get_base_sixty( tptr);
         if( tptr[-1] == '-')
            dec = -dec;
         }
      if( !memcmp( buff, "Img ", 4) || !memcmp( buff, "Pass ", 5))
         if( exposure > min_exposure)
            {
            char time_buff[20];

            tptr = find_base_60( buff);
            assert( tptr);
            add_half_exposure( time_buff, tptr, exposure);
            printf( "%.3f,%.3f,%sT%s,691,na\n", ra, dec,
                                iso_date, time_buff);
            }
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
