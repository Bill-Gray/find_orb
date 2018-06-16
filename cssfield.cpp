/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

/*
Support code for a new Find_Orb feature:  the ability to integrate an
orbit backward to find out on which field(s) you might have caught it
in the past.  This should support a SkyMorph-like precovery ability.

To do it,  we need a list of previously imaged fields.  Catalina Sky
Survey has kindly provided a list of (at present) over three million
fields.  This code reads in the field information (currently just CSS
fields) and generate a list in binary format,  sorted in descending
order by date.   The input files,  provided by Rob Seaman,  look
like this :

RA dec time obscode (ignorable YYYYMMDD) image_name
63.3921,26.7177,2003-03-23T02:42:04.406,703,20030323,01_03MAR23_N26021_0001.arch.H
60.2856,26.7153,2003-03-23T02:43:36.860,703,20030323,01_03MAR23_N26020_0001.arch.H
57.1807,26.7128,2003-03-23T02:45:10.025,703,20030323,01_03MAR23_N26019_0001.arch.H
54.0746,26.7103,2003-03-23T02:46:42.230,703,20030323,01_03MAR23_N26018_0001.arch.H
50.9695,26.7078,2003-03-23T02:48:15.354,703,20030323,01_03MAR23_N26017_0001.arch.H
47.8647,26.7051,2003-03-23T02:49:47.609,703,20030323,01_03MAR23_N26016_0001.arch.H
63.3961,26.7174,2003-03-23T02:52:25.776,703,20030323,01_03MAR23_N26021_0002.arch.H
60.2935,26.7142,2003-03-23T02:53:58.001,703,20030323,01_03MAR23_N26020_0002.arch.H
57.1873,26.7113,2003-03-23T02:55:31.406,703,20030323,01_03MAR23_N26019_0002.arch.H
54.0833,26.7085,2003-03-23T02:57:05.100,703,20030323,01_03MAR23_N26018_0002.arch.H

The file is read in and the fields in the following field_location_t struct
are filled out,  and the array sorted in descending order of JD.  (The idea
is that we'll usually be starting with a "recent" object,  and will search
backward for plates.  No need to integrate back to the start of CSS history
and then integrate forward.  With the fields sorted by date,  we don't
have to integrate backward,  forward,  backward,  etc.)

The resulting index is written to 'css.idx',  with a bit of text dumped
to persuade me that correct values are being written.

It may be useful to split the above data among two or more files.  (For
example, one may have a list of fields from first light up to a week ago,
and a smaller file that just has the last few nights.  You update the latter,
but the big file stays unchanged.)  For this reason,  the index actually
covers up to ten files,  css0.csv through css9.csv.  Adding more would not
be particularly difficult.

By default,  the field sizes are determined by get_field_size(),  which knows
about CSS history and can figure it out based on station code and image date.
Similar lines can be added for other codes,  given suitable data.  This can
be overridden by appending,  after the image file name,

,!(height),(width),(tilt angle)

The ! serves to notify this program that any default height/width will be
overridden by the subsequent values.  If there's just one such value,  the
field is assumed to be square.  If there are only two,  the (rectangular)
field is assumed to be aligned with the RA/dec axes of date.  If there are
three,  it's assumed that the image was tilted by that angle relative to the
RA/dec axes of date.

For Spacewatch,  for example,  three out of four fields will be in the
default "landscape" orientation.  The remaining quarter will probably use
the above scheme to specify the image height and width,  or possibly the
same height/width and a 90 degree rotation angle.

Also,  note that a generic time/date parser is used.  CSS gave FITS-style
times,  and that's what is shown in the above examples.  But one can
instead specify JD,  MJD,  and a variety of unlikely formats.  See

https://projectpluto.com/update8d.htm#time_entry

for a full list of time specification options.  */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "watdefs.h"
#include "date.h"

#define PI \
   3.1415926535897932384626433832795028841971693993751058209749445923

#pragma pack( 1)

typedef struct
{
   double ra, dec, jd;
   double height, width, tilt;
   uint32_t file_offset;
   char obscode[4];
   char file_number;
} field_location_t;

#pragma pack( )

int field_compare( const void *a, const void *b)
{
   const field_location_t *aptr = (const field_location_t *)a;
   const field_location_t *bptr = (const field_location_t *)b;

   return( aptr->jd > bptr->jd ? -1 : 1);
}                 /* sort in _descending_ order by date */

/*
CSS field sizes,  from Eric Christensen,  2017 Mar.  Updated by
Rob Seaman,  2017 Nov 04.

703 field size, 2003-2016: 2.85 x 2.85 deg.
G96, 2004 - May 2016: 1.1 x 1.1 deg.
E12 : 2.05 x 2.05 deg.

10K cameras:
G96 - May 2016 - present: 2.2 x 2.2 deg.
703 - Dec. 2016 - present: 4.4 x 4.4 deg.

(566) had a 4096-square CCD at 1.43 arcsec/pixel.

(644): The Palomar NEAT "Tri-Camera" had three separate 4096x4096
pixel CCDs, running north to south,  pixel scale 1.01 arcsec/pixel.
The images have similar RAs (probably almost identical in RA of
date),  and decs spaced out by about 1.3 degrees.  The pointing log
gives three lines for any given time,  corresponding to the northern,
middle,  and southern images.

(691) Spacewatch size roughly from MPC sky coverage files,  plus
some info from Bob McMillan.  The actual shape is a little more
complicated than this -- it's a mosaic of eight chips --  but
the 1.85 x 1.73 degree rectangle is close enough to do a rough
cut as to whether we got the object.  */

static void get_field_size( double *width, double *height, const double jd,
                        const char *obs_code)
{
   const double dec_02_2016 = 2457724.5;
   const double may_26_2016 = 2457534.5;
   const double jun_24_2005 = 2453545.5;     /* see J95 */
   static char bad_code[10];

   *height = 0.;
   switch( *obs_code)
      {
      case '5':         /* (566) Haleakala-NEAT/GEODSS  */
         *width = 4096. * 1.43 / 3600.;
         break;
      case '6':
         if( obs_code[1] == '9')       /* (691) Spacewatch */
            {
            *width = 1.85;
            *height = 1.73;
            }
         else            /* (644) NEAT at Palomar Sam'l Oschbin Schmidt */
            *width = 4096. * 1.01 / 3600.;
         break;
      case '7':         /* 703 */
         *width = (jd < dec_02_2016 ? 2.85 : 4.4);
         break;
      case 'G':         /* G96 */
         *width = (jd < may_26_2016 ? 1.1 : 2.25);
         break;
      case 'E':         /* E12 */
         *width = 2.05;
         break;
      case 'I':         /* I52:  33' field of view;  some loss in corners */
         *width = 33. / 60.;
         break;
      case 'J':         /* J95:  25' to 2005 jun 22, 18' for 2005 jun 27 on */
         *width = (jd < jun_24_2005 ? 25. / 60. : 18. / 60.);
         break;
      case 'V':         /* V06:  580" field of view */
         *width = 580. / 3600.;
         break;
      default:
         *width = 0.;
         if( memcmp( bad_code, obs_code, 4))
            {
            printf( "Bad code '%.3s', shouldn't be here\n", obs_code);
            memcpy( bad_code, obs_code, 4);
            }
         break;
      }
   if( !*height)           /* square field indicated */
      *height = *width;
   *width *= PI / 180.;
   *height *= PI / 180.;
}

static void get_field_size_from_input( double *width, double *height,
                   double *tilt, const char *buff)
{
   const char *bang = strchr( buff, '!');

   if( bang)
      {
      const int n_scanned = sscanf( bang + 1, "%lf %lf %lf", width, height, tilt);

      if( n_scanned < 1 || *width <= 0.)
         {
         printf( "Error reading field sizes\n%s\n", buff);
         exit( -1);
         }
      assert( n_scanned >= 1 && n_scanned <= 3);
      if( n_scanned == 1)     /* only one dimension given */
         *height = *width;
      *width *= PI / 180.;
      *height *= PI / 180.;
      *tilt *= PI / 180.;
      }
}

int main( const int argc, const char **argv)
{
   int file_number;
   FILE *ofile;
   char buff[200];
   field_location_t *rval = NULL;
   int n_alloced = 0, n = 0, i;
   int verbose = 0, n_duplicates = 0;

   setvbuf( stdout, NULL, _IONBF, 0);
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'v':
               verbose = 1;
               break;
            default:
               printf( "'%s' is unrecognized\n", argv[i]);
               return( -1);
            }

   for( file_number = 0; file_number < 10; file_number++)
      {
      FILE *ifile;

      sprintf( buff, "css%d.csv", file_number);
      ifile = fopen( buff, "rb");
      if( ifile)
         {
         double min_jd = 1e+10, max_jd = 0.;
         double tilt = 0.;

         printf( "%s opened;  reading fields\n", buff);
         while( fgets( buff, sizeof( buff), ifile))
            if( *buff != '#')       /* allow for comments */
               {
               char timestr[80];
               size_t i;
               int n_scanned;

               for( i = 0; buff[i]; i++)
                  if( buff[i] == ',')
                     buff[i] = ' ';
               if( n >= n_alloced)
                  {
                  n_alloced += 200 + n_alloced / 2;
                  rval = (field_location_t *)realloc( rval,
                                        n_alloced * sizeof( field_location_t));
                  }
               n_scanned = sscanf( buff, "%lf %lf %70s %3s", &rval[n].ra,
                           &rval[n].dec, timestr, (char *)&rval[n].obscode);
               assert( n_scanned == 4);
               rval[n].jd = get_time_from_string( 0., timestr, 0, NULL);
               if( rval[n].jd)      /* some lines have 'time=NULL' */
                  {
                  rval[n].ra  *= PI / 180.;
                  rval[n].dec *= PI / 180.;
                  rval[n].tilt = tilt;
                  rval[n].file_offset = (uint32_t)( ftell( ifile) - strlen( buff));
                  rval[n].file_number = (char)file_number;
                  get_field_size( &rval[n].width, &rval[n].height, rval[n].jd,
                                       rval[n].obscode);
                  get_field_size_from_input( &rval[n].width, &rval[n].height,
                                    &rval[n].tilt, buff);
                  if( min_jd > rval[n].jd)
                     min_jd = rval[n].jd;
                  if( max_jd < rval[n].jd)
                     max_jd = rval[n].jd;
                  n++;
                  }
               if( n < 10 || n % 100000 == 0)
                  printf( "%d fields read and parsed\r", n);
               }
            else if( !memcmp( buff, "# Tilt: ", 8))
               tilt = atof( buff + 8);
         fclose( ifile);
         printf( "\n%d fields found\n", n);
         full_ctime( buff, min_jd, 0);
         printf( "Fields run from %.21s to ", buff);
         full_ctime( buff, max_jd, 0);
         printf( "%.21s\n", buff);
         }
      }
   if( verbose)
      for( i = 0; i < 20 && i < n; i++)
         printf( "%f %f %f: %ld %s %f\n", rval[i].jd,
                     rval[i].ra, rval[i].dec, (long)rval[i].file_offset, rval[i].obscode,
                     rval[i].height);
   if( !n)
      {
      printf( "No fields were actually read in.  Check to see that at least one\n"
              "of the files css0.csv, css1.csv, ... exists and has valid field\n"
              "data in it.\n");
      return( -2);
      }
   printf( "Sorting...\n");
   qsort( rval, n, sizeof( rval[0]), field_compare);
   if( verbose)
      for( i = 0; i < 20; i++)
         printf( "%f %f %f: %ld %s\n", rval[i].jd,
                     rval[i].ra, rval[i].dec, (long)rval[i].file_offset, rval[i].obscode);
   full_ctime( buff, rval[n - 1].jd, 0);
   printf( "Full time span is %.21s", buff);
   full_ctime( buff, rval[  0  ].jd, 0);
   printf( " to %.21s\n", buff);
   ofile = fopen( "css.idx", "wb");
   if( !ofile)
      perror( "opening ofile failed");
   else if( fwrite( rval, sizeof( rval[0]), n, ofile) != (size_t)n)
      perror( "write failure");
   for( i = 0; i < n - 1; i++)
      if( rval[i].jd == rval[i + 1].jd
                     && !strcmp( rval[i].obscode, rval[i + 1].obscode))
         {
         n_duplicates++;
         if( verbose)
            {
            full_ctime( buff, rval[i].jd, 0);
            printf( "Duplicate found: %.30s from %.4s\n",
                           buff, rval[i].obscode);
            }
         }
   printf( "%d duplicates found\n", n_duplicates);
   return( 0);
}
