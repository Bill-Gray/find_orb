/* sw_xvt.c: converter for Spacewatch logs

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

/* Code to convert "old format" Spacewatch pointing log into the form
expected by 'cssfield.c'.   Bob McMillan tells me the last log using
this format was the 2016 December 24 one;  see 'sw_xvt2.c' for code to
handle the subsequent format.  Note comments at bottom of this file on
the format of these logs.  Also see 'cssfield.c' for comments on the
format of the output from this program.  Compile with

gcc -Wall -O3 -o sw_xvt sw_xvt.c

   We read lines from the SW log.  If it starts with

    UT DATE: mmm dd, yyyy

   we reset the date (converted to ISO form).  If it's got an RA
and a dec in it,  we look after those for times of observation,  of
which there can be up to three.  For each of those times of observation,
we output a line with the RA and dec in decimal degrees,  the ISO
formatted observation time,  the MPC code 691,  and 'na' (the filename
is Not Applicable here.)

   That puts everything in the 'standard' CSV format inhaled by
'cssfield',  dumping to stdout.

   It assumes the log file is

Spacewatch_mosaic.log.2003Mar23-2016Dec24.txt

   You can specify a different name on the command line.       */

static void get_iso_date( char *iso_date, const char *buff)
{
   char month[6];
   const char *tptr, *months = "JanFebMarAprMayJunJulAugSepOctNovDec";
   int day, year, m;

   printf( "# %s", buff);
   if( sscanf( buff, "%4s %d, %d", month, &day, &year) != 3)
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
                    "Spacewatch_mosaic.log.2003Mar23-2016Dec24.txt");
   FILE *ifile = fopen( filename, "rb");
   char buff[200];
   char iso_date[30];

   assert( ifile);
   *iso_date = '\0';
   printf( "# Spacewatch 'old formula' logs,  processed with sw_xvt.c (q.v.)\n");
#ifdef __TIMESTAMP__
   printf( "# Version %s\n", __TIMESTAMP__);
#else
   printf( "# Version %s %s\n", __DATE__, __TIME__);
#endif
   printf( "# Input file '%s'\n", filename);
   while( fgets( buff, sizeof( buff), ifile))
      if( strlen( buff) > 40 && buff[1] != ':')
         {
         const char *tptr = find_base_60( buff);

         if( !memcmp( buff + 1, " UT DATE: ", 9))
            get_iso_date( iso_date, buff + 10);
         else if( tptr)
            {
            const double ra = get_base_sixty( tptr) * 15.;

            tptr = find_base_60( tptr + 8);
            if( tptr)
               {
               double dec = get_base_sixty( tptr);
               const double exposure = atoi( tptr + 9);
               char time_buff[20];

               if( tptr[-1] == '-')
                  dec = -dec;
               while( (tptr = find_base_60( tptr + 8)) != NULL)
                  {
                  add_half_exposure( time_buff, tptr, exposure);
                  printf( "%.3f,%.3f,%sT%s,691,na\n", ra, dec,
                                                  iso_date, time_buff);
                  }
               }
            }
         memset( buff, 0, 120);
         }
   fclose( ifile);
   return( 0);
}

/*   (The following is an edited snip from an e-mail I sent to Bob
McMillan,  in which I attempted to puzzle out the format of the
Spacewatch pointing log.)

   An example of text from the Spacewatch pointing log :

    UT DATE: Mar 23, 2003   OBSERVER(s): J.A.Larsen
   Reg.Center    Dirs                RA          Dec     sec       Time1 FWHM      Time2 FWHM      Time3 FWHM
   ==========================================================================================================
   06.01         1,12,22       11:43:12    +05:11:48     120    06:11:20 1.7"   07:14:51 1.6"   07:52:47 1.8"

   This is telling me that on 2003 Mar 23,  three images were taken centered
at RA=11:43:12,  dec=+05:11:48 (in epoch of date),  at 06:11:20, 07:14:51,
and 07:52:47 UTC. The exposures were all 120 seconds long,  with the seeing
given as FWHM.  J. A. Larsen was the observer.

   In the MPC sky coverage files,  it looks as if this is translated into
a single rectangle,  about 1.85 degrees wide and 1.73 degrees high.  This
is also the approach used in 'sw_xvt.c':  the above line is translated
into three rectangles for the same piece of the sky,  imaged at three
different times,  plus a comment line,  resulting in :

#  Mar 23, 2003   OBSERVER(s): J.A.Larsen
175.800,5.197,2003-03-23T06:12:20.000,691,na
175.800,5.197,2003-03-23T07:15:51.000,691,na
175.800,5.197,2003-03-23T07:53:47.000,691,na

   So from the standpoint of 'sw_xvt.c',  the above is all you need to know.
In theory,  though,  the rectangles output by the program are really just a
first step.  Each of those 1.85 by 1.73 rectangles should be broken up into
eight image files,  corresponding to two images each from four CCDs.

http://adsabs.harvard.edu/abs/2007IAUS..236..329M

   shows a layout of four 4608x2048 CCDs,  three horizontal and one vertical,
covering 2.9 square degrees.

   Each of the four CCDs shown below is 4608x2048 pixels,  13.5 microns
square.  A sketch made by Marcus L. Perry,  then Chief Engineer,  shows that
each chip has an edge on the long sides of about 0.260 +/- 0.050 mm,  and
on the short sides of 0.120 +/- 0.050 mm.  Add that in,  and the "total,
real" small axis of each CCD is 2048 * 0.0135 + 2 * 0.260 = 28.168 mm;  the
real large axis is 4608 * 0.0135 + 2 * 0.120 = 62.448 mm.

   Between the three parallel CCDs,  there's a gap of 0.508 +/- 0.025 mm.
So the total distance from A to C is (3*28.168) + (2*0.508) = 85.520 mm.

   Between those three CCDs and the one shown vertically below,  there
is also a 0.508 +/- 0.025 mm gap.  So the horizontal distance from A to
B (or C to D) is width plus height plus one gap,  or 91.124 mm.


  A
   +---------------+        B
   |               +-------+
   |               |       |     (figure 1 of the above paper,  in ASCII...
   |               |       |     best seen with fixed-width fonts)
   +---------------+       |
   |               |       |
  E|          X    |       |F
   |               |       |
   +---------------+       |
   |               |       |
   |               |       |
   |               +-------+
   +---------------+        D
  C

   X,  the center point,  is equidistant from A, B, C,  and D.
Call the distance EX = x;  then x^2 + EA^2 = (EF - x)^2 + FB^2,
or EA^2 = EF^2 - 2EFx + FB^2,  or x = (FB^2 + EF^2 - EA^2) / 2EF,
where EA = 85.520 / 2 = 42.76 mm, FB = 62.448 / 2 = 31.224 mm,
and EF = 91.124 mm.  Hence,  x = 40.879 mm.
*/
