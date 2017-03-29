/*
Support code for a new Find_Orb feature:  the ability to integrate an
orbit backward to find out on which field(s) you might have caught it
in the past.  This should support a SkyMorph-like precovery ability.

To do it,  we need a list of previously imaged fields.  Catalina Sky
Survey has kindly provided a list of (at present) over three million
fields.  This code reads in the field information (currently just CSS
fields) and generate a list in binary format,  sorted in descending
order by date.   The input 'css_index.csv' file,  provided by Rob
Seaman,  looks like this :

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
*/

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "watdefs.h"
#include "date.h"

#define PI \
   3.1415926535897932384626433832795028841971693993751058209749445923

typedef struct
{
   double ra, dec, jd;
   double height, width, tilt;
   uint32_t file_offset;
   char obscode[4];
} field_location_t;

int field_compare( const void *a, const void *b)
{
   const field_location_t *aptr = (const field_location_t *)a;
   const field_location_t *bptr = (const field_location_t *)b;

   return( aptr->jd > bptr->jd ? -1 : 1);
}                 /* sort in _descending_ order by date */

/*
CSS field sizes,  from Eric Christensen,  2017 Mar :

703 field size, 2003-2016: 2.85 x 2.85 deg.
G96, 2004 - May 2016: 1.1 x 1.1 deg.
E12 : 2.05 x 2.05 deg.

10K cameras:
G96 - May 2016 - present: 2.2 x 2.2 deg.
703 - Dec. 2016 - present: 4.4 x 4.4 deg.
*/

static void get_field_size( double *width, double *height, const double jd,
                        const char *obs_code)
{
   const double dec_01_2016 = 2457723.5;
   const double may_01_2016 = 2457509.5;

   switch( *obs_code)
      {
      case '7':         /* 703 */
         *width = (jd < dec_01_2016 ? 2.85 : 4.4);
         break;
      case 'G':         /* G96 */
         *width = (jd < may_01_2016 ? 1.1 : 2.25);
         break;
      case 'E':         /* E12 */
         *width = 2.05;
         break;
      case 'I':         /* I52:  33' field of view;  some loss in corners */
         *width = 33. / 60.;
         break;
      default:
         *width = 0.;
         printf( "Bad code '%s', shouldn't be here\n", obs_code);
         break;
      }
   *width *= PI / 180.;
   *height = *width;    /* all square fields thus far */
}


int main( const int argc, const char **argv)
{
   const char *filename = "css_index.csv";
   FILE *ifile = fopen( filename, "rb");
   FILE *ofile;
   char buff[200];
   field_location_t *rval = NULL;
   int n_alloced = 0, n = 0, i;
   uint32_t file_loc = 0;
   int verbose = 0, n_duplicates = 0;

   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'f':
               filename = argv[i] + 2;
               if( !*filename && i < argc - 1)
                  filename = argv[i + 1];
               break;
            case 'v':
               verbose = 1;
               break;
            default:
               printf( "'%s' is unrecognized\n", argv[i]);
               return( -1);
            }

   ifile = fopen( filename, "rb");
   if( !ifile)
      {
      printf( "%s not opened\n", filename);
      return( -1);
      }
   printf( "%s opened;  reading fields\n", filename);
   setvbuf( stdout, NULL, _IONBF, 0);
   while( fgets( buff, sizeof( buff), ifile))
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
         rval[n].file_offset = file_loc;
         get_field_size( &rval[n].width, &rval[n].height, rval[n].jd,
                              rval[n].obscode);
         n++;
         }
      file_loc += strlen( buff);
      if( n < 10 || n % 100000 == 0)
         printf( "%d fields read and parsed\r", n);
      }
   fclose( ifile);
   printf( "\n%d fields found\n", n);
   if( verbose)
      for( i = 0; i < 20; i++)
         printf( "%f %f %f: %ld %s %f\n", rval[i].jd,
                     rval[i].ra, rval[i].dec, (long)rval[i].file_offset, rval[i].obscode,
                     rval[i].height);
   printf( "Sorting...\n");
   qsort( rval, n, sizeof( rval[0]), field_compare);
   if( verbose)
      for( i = 0; i < 20; i++)
         printf( "%f %f %f: %ld %s\n", rval[i].jd,
                     rval[i].ra, rval[i].dec, (long)rval[i].file_offset, rval[i].obscode);
   full_ctime( buff, rval[n - 1].jd, 0);
   printf( "Fields span %.30s", buff);
   full_ctime( buff, rval[  0  ].jd, 0);
   printf( " to %.30s\n", buff);
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
