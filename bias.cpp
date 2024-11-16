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

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

void ra_dec_to_xy( const double ra, const double dec, double *x, double *y);
unsigned xy_to_healpix( const double x, const double y, const unsigned N);
int find_fcct_biases( const double ra, const double dec, const char catalog,
                 const double jd, double *bias_ra, double *bias_dec);

#define BIAS_NO_DATA_FOR_THAT_CATALOG     -1
#define BIAS_NO_BIAS_FILE                 -2
#define BIAS_NOT_ALLOCATED                -3
#define BIAS_FILE_WRONG_SIZE              -4

/* Compute star catalog position/proper motion biases,  from tables by
Farnocchia,  Chesley,  Chamberlin,  Tholen,  _Icarus_ 245 (2015) 94-111,
http://adsabs.harvard.edu/abs/2015Icar..245...94F  ,  "FCCT14";  or
from its successor,  Eggl,  Farnocchia, Chamberlin, Chesley,  "EFCC18",
https://arxiv.org/abs/1909.04558.  You should really use the latter;
FCCT14 is now largely of historical interest.

   The actual debiasing files,  named 'bias.dat' for both,  are in

ftp://ssd.jpl.nasa.gov/pub/ssd/debias/debias_2014.tgz
ftp://ssd.jpl.nasa.gov/pub/ssd/debias/debias_2018.tgz

   The .tgz files contain two other files,  README.txt and tiles.dat.  As
far as Find_Orb is concerned,  you can ignore these.

   In both FCCT14 and EFCC18,  the sky is divided into 12 * 64 * 64 =
49152 "tiles" of equal area and passably square shape,  using the HEALPIX
tesselation (see 'healpix.cpp'). Within each tile, the bias of 19 (for
FCCT14) or 26 (for EFCC18) different catalogs in RA, dec,  and proper
motion in RA and dec are given.  The biases are relative to a "presumed
good reference".  For FCCT14,  that was subset of PPMXL.  ECFF18 was
able to use Gaia-DR2,  an improvement.  Those biases are then stored in
'bias.dat'.

   This code reads in all 49152 * n_cats biases (each of which is really
four biases,  one each in RA,  dec,  pm_RA,  pm_dec).  For a given
position,  it figures out the tile covering that point,  then extracts
the biases corresponding to the catalog in which we're interested and
computes the bias in RA and dec for a particular JD,  storing them in
*bias_ra and *bias_dec.       */

const char *fcct14_bias_file_name = NULL;
         /* override this if the bias file is elsewhere */

FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */

int find_fcct_biases( const double ra, const double dec, const char catalog,
                 const double jd, double *bias_ra, double *bias_dec)
{
   static int16_t *bias_data = NULL;
   static char catalog_codes[30];
   const char *catalog_loc;
         /* Number of HEALPIX tiles at level N=64 : */
   const int n_tiles = 12 * 64 * 64;
   static int n_cats = 0;
   int loc, i;
   double x, y, year;
   const double j2000 = 2451545.;      /* JD 2451545.0 = 2000 Jan 1.5 */
   static bool bias_file_known_to_be_missing = false;

   if( catalog == (char)-1)    /* free up internal data */
      {
      if( bias_data)
         free( bias_data);
      bias_data = NULL;
      return( 0);
      }

   if( bias_ra && bias_dec)
      *bias_ra = *bias_dec = 0.;
   if( bias_file_known_to_be_missing)
      return( BIAS_NO_BIAS_FILE);
   if( !bias_data)
      {
      FILE *ifile;
      const char *crunched_fcct14_file_name = "biases.bin";
      size_t nbytes, csize = 0;

      ifile = fopen_ext( crunched_fcct14_file_name, "crb");
      if( ifile)     /* not our first time:  we have a compressed file */
         {
         unsigned char ncats;

         if( fread( &ncats, 1, 1, ifile) != 1)
            return( BIAS_FILE_WRONG_SIZE);
         n_cats = (int)ncats;
         assert( n_cats == 19 || n_cats == 26);
         if( fread( catalog_codes, 1, n_cats, ifile) != (size_t)n_cats)
            return( BIAS_FILE_WRONG_SIZE);
         catalog_codes[n_cats] = '\0';
         csize = n_tiles * n_cats * 4 * sizeof( int16_t);
         bias_data = (int16_t *)malloc( csize + 1);
         if( !bias_data)
            return( BIAS_NOT_ALLOCATED);
         nbytes = fread( bias_data, 1, csize + 1, ifile);
         fclose( ifile);
         if( nbytes != csize)
            return( BIAS_FILE_WRONG_SIZE);
         }
      else
         {
         char buff[900];   /* first time: read ASCII file & binary-ize; */
         FILE *ofile;     /* store binary version for all future use   */

         ifile = NULL;
         if( fcct14_bias_file_name && * fcct14_bias_file_name)
            ifile = fopen( fcct14_bias_file_name, "rb");
         if( !ifile)
            ifile = fopen_ext( "bias.dat", "crb");
         if( !ifile)
            {
            bias_file_known_to_be_missing = true;
            free( bias_data);
            bias_data = NULL;
            return( BIAS_NO_BIAS_FILE);
            }
         loc = 0;
         while( loc < n_tiles && fgets( buff, sizeof( buff), ifile))
            if( *buff != '!' && bias_data)
               {
               char *tptr = buff;
               int16_t *bptr = bias_data + loc * n_cats * 4;

               for( i = 0; i < n_cats * 4; i++, bptr++)
                  {
                  bool is_negative = false;

                  while( *tptr == ' ')
                     tptr++;
                  if( *tptr == '-')
                     {
                     is_negative = true;
                     tptr++;
                     }
                  *bptr = 0;
                  while( *tptr != ' ')
                     {
                     if( *tptr >= '0' && *tptr <= '9')
                        *bptr = (int16_t)( *bptr * 10 + *tptr - '0');
                     tptr++;
                     }
                  if( is_negative)
                     *bptr = -*bptr;
                  }
               loc++;
               }
            else if( !memcmp( buff, "! |---", 6))
               {        /* header line */
               n_cats = (int)strlen( buff) / 29;
               assert( n_cats == 19 || n_cats == 26);
               for( i = 0; i < n_cats; i++)
                  {
                  char *tptr = buff + 22 + i * 29;

                  assert( *tptr == ' ' && tptr[1] == '-');
                  while( tptr > buff && *tptr == ' ')
                     tptr--;
                  catalog_codes[i] = *tptr;
                  }
               catalog_codes[n_cats] = '\0';
               csize = n_tiles * n_cats * 4 * sizeof( int16_t);
               bias_data = (int16_t *)malloc( csize + 1);
               if( !bias_data)
                  return( BIAS_NOT_ALLOCATED);
               }
         fclose( ifile);
         assert( loc == n_tiles);

         ofile = fopen_ext( crunched_fcct14_file_name, "fcwb");
         *buff = (char)n_cats;
         nbytes = fwrite( buff, 1, 1, ofile);
         assert( nbytes == 1);
         nbytes = fwrite( catalog_codes, 1, n_cats, ofile);
         assert( nbytes == (size_t)n_cats);
         nbytes = fwrite( bias_data, 1, csize, ofile);
         assert( nbytes == csize);
         fclose( ofile);
         }
      }
   if( n_cats == 26 && catalog_codes[25] == 'W')   /* fix a mistaken catalog */
      catalog_codes[25] = 'Y';                /* code in pre-2023 EFCC files */
   if( !catalog)        /* just inquiring as to which version we have */
      return( n_cats == 26 ? 2018 : 2014);
   catalog_loc = strchr( catalog_codes, catalog);
   if( !catalog_loc)    /* don't have bias data for this catalog */
      return( BIAS_NO_DATA_FOR_THAT_CATALOG);
   ra_dec_to_xy( ra, dec, &x, &y);
   loc = xy_to_healpix( x, y, 64);
// debug_printf( "catalog %c, tile %d: ", catalog, loc);
   loc = loc * n_cats * 4 + (int)(catalog_loc - catalog_codes) *4;
   year = (jd - j2000) / 365.25;
// debug_printf( "%d %d %d %d\n", bias_data[loc], bias_data[loc + 1],
//                   bias_data[loc + 2], bias_data[loc + 3]);
   *bias_ra  =  (double)bias_data[loc]     / 1000.     /* RA bias, in arcsec */
       + year * (double)bias_data[loc + 2] / 100000.;  /* RA pm bias, in mas/yr */
   *bias_dec =  (double)bias_data[loc + 1] / 1000.     /* dec bias, in arcsec */
       + year * (double)bias_data[loc + 3] / 100000.;  /* dec pm bias, in mas/yr */
   return( 0);
}
