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

/* Compute star catalog position/proper motion biases,  from tables by
Farnocchia,  Chesley,  Chamberlin,  Tholen,  _Icarus_ 245 (2015) 94-111,
http://adsabs.harvard.edu/abs/2015Icar..245...94F  ,  "FCCT14".  The
actual debiasing file,  'bias.dat',  is in

ftp://ssd.jpl.nasa.gov/pub/ssd/debias/debias_2014.tgz

   The .tgz file contains two other files,  README.txt and tiles.dat.  As
far as Find_Orb is concerned,  you can ignore these.

   FCCT14 divides the sky into 12 * 64 * 64 = 49152 "tiles" of equal area and
passably square shape,  using the HEALPIX tesselation (see 'healpix.cpp').
Within each tile, the bias of 19 different catalogs in RA, dec,  and proper
motion in RA and dec are given.  (The biases are relative to a "presumed good
reference" subset of PPMXL.)  Those biases are stored in 'bias.dat'.

   This code reads in all 49152 * 19 biases (each of which is really
four biases,  one each in RA,  dec,  pm_RA,  pm_dec).  For a given
position,  it figures out the tile covering that point,  then extracts
the biases corresponding to the catalog in which we're interested and
computes the bias in RA and dec for a particular JD,  storing them in
*bias_ra and *bias_dec.       */

int find_fcct_biases( const double ra, const double dec, const char catalog,
                 const double jd, double *bias_ra, double *bias_dec)
{
   static int16_t *bias_data = NULL;
   static const char *catalog_codes = "abcdegijlmopqruvwLN";
   const char *catalog_loc = strchr( catalog_codes, catalog);
         /* Number of HEALPIX tiles at level N=64 : */
   const int n_tiles = 12 * 64 * 64;
   const int n_cats = 19;
   int loc, i;
   double x, y, year;
   const double j2000 = 2451545.;      /* JD 2451545.0 = 2000 Jan 1.5 */
   static bool bias_file_known_to_be_missing = false;

   if( !bias_ra)         /* free up internal data */
      {
      if( bias_data)
         free( bias_data);
      bias_data = NULL;
      return( 0);
      }

   *bias_ra = *bias_dec = 0.;
   if( !catalog_loc)    /* don't have bias data for this catalog */
      return( BIAS_NO_DATA_FOR_THAT_CATALOG);
   if( bias_file_known_to_be_missing)
      return( BIAS_NO_BIAS_FILE);
   if( !bias_data)
      {
      FILE *ifile = fopen( "bias.dat", "rb");
      char buff[600];

      if( !ifile)
         {
         bias_file_known_to_be_missing = true;
         return( BIAS_NO_BIAS_FILE);
         }
      bias_data = (int16_t *)calloc( n_tiles, n_cats * 4 * sizeof( int16_t));
               /* The above corresponds to 12*64*64*19*4 = 3735552 shorts,  or
                  or 7471104 bytes.  The calloc() _could_ fail : */

      if( !bias_data)
         {
         fclose( ifile);
         return( BIAS_NOT_ALLOCATED);
         }
      loc = 0;
      while( loc < n_tiles && fgets( buff, sizeof( buff), ifile))
         if( *buff != '!')
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
      fclose( ifile);
      assert( loc == n_tiles);
      }
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
