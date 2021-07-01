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

#include <math.h>

/* HEALPix functions.   The algorithm to subdivide the sky in equal area
tiles is described by Gorsky et al. 2005
(http://adsabs.harvard.edu/abs/2005ApJ...622..759G).  All that is implemented
here is the conversion between RA/dec and HEALPix cell numbers,  which is all
that's needed to figure out which lines in the bias data should be read. */

#define ASIN_TWO_THIRDS 0.7297276562269663634547966598133206953965059140477136907089494914618189
#define pi 3.1415926535897932384626433832795028841971693993751058209749445923

void ra_dec_to_xy( const double ra, const double dec, double *x, double *y);
unsigned xy_to_healpix( const double x, const double y, const unsigned N);

/* The HEALPix projection uses an equal-area cylindrical projection for the
central 2/3 (by area) zone nearest the equator,  and puts each polar region
into four right triangles in an interrupted Collignon projection.  See
https://en.wikipedia.org/wiki/HEALPix for a picture of this worth many thousands
of comments.  Here,  y=0 at the equator;  y=2 at the north pole,  y=-2 at the
south pole,  and 0 to 360 in RA is projected to 0 to 8 in x.  RA and dec are
assumed to be in radians.        */

void ra_dec_to_xy( const double ra, const double dec, double *x, double *y)
{
   if( dec < -ASIN_TWO_THIRDS)
      {
      ra_dec_to_xy( ra, -dec, x, y);
      *y = -*y;
      return;
      }
   *x = fmod( ra * 4 / pi, 8.);
   if( *x < 0.)
      *x += 8.;
   if( dec < ASIN_TWO_THIRDS)
      *y = 1.5 * sin( dec);
   else
      {
      const double z = sqrt(3. * (1. - sin( dec)));
      const double center = floor( *x / 2.) * 2. + 1.;

      *y = 2. - z;
      *x = center + z * (*x - center);
      }
}

#ifdef NOT_NEEDED_YET
      /* So far,  we aren't doing reverse conversions.  Note that this is, */
      /* as yet,  untested code.  Though it's a straightforward inverse of */
      /* the above ra_dec_to_xy function.                                  */
static void xy_to_ra_dec( const double x, const double y, double *ra, double *dec)
{
   if( y < -1.)
      {
      xy_to_ra_dec( x, -y, ra, dec);
      *dec = -*dec;
      return;
      }
   if( y < 1.)
      {
      *ra = x * pi / 4.;
      *dec = asin( y * 2. / 3.);
      }
   else
      {
      const double z = 2. - y;
      const double center = floor( (x + 1.) / 2.) * 2.;

      *ra = (center + (x - center) / z) * pi / 4.;
      *dec = asin( 1. - z * z / 3.);
      }
}
#endif

/* HEALPix pixels are squares (when projected) tilted at a 45-degree angle.
It helps to rotate them as shown below,  with the origin shifted to the top
of the leftmost polar triangle.  (In the following, to show the 192 cells of
N=2 with two digits,  I mix hexadecimal and decimal; the cell numbers run
from ...98, 99, to a0, a1... and similarly up to the maximum j1 = cell 191.)


                                    03 11 23 39
                                    10 22 38 55 71
                                    21 37 54 70 87 a3
                                    36 53 69 86 a2 b9 d5
                        02 09 20 35 52 68 85 a1 b8 d4 f1 g7
i=0 1  2  3  4  5  6    08 19 34 51 67 84 a0 b7 d3 f0 g6 h9
                        18 33 50 66 83 99 b6 d2 e9 g5 h8 i7
                        32 49 65 82 98 b5 d1 e8 g4 h7 i6 j1
            01 07 17 31 48 64 81 97 b4 d0 e7 g3
            06 16 30 47 63 80 96 b3 c9 e6 g2 h6
            15 29 46 62 79 95 b2 c8 e5 g1 h5 i5      -2     (i, j) rotated to keep cell
            28 45 61 78 94 b1 c7 e4 g0 h4 i4 j0      -1     00 at the origin
00 05 14 27 44 60 77 93 b0 c6 e3 f9                 j=0
04 13 26 43 59 76 92 a9 c5 e2 f8 h3                   1
12 25 42 58 75 91 a8 c4 e1 f7 h2 i3                   2        N = 4 = power of 2
24 41 57 74 90 a7 c3 e0 f6 h1 i2 i9                   3
40 56 73 89 a6 c2 d9 f5                               4
   72 88 a5 c1 d8 f4 h0
      a4 c0 d7 f3 g9 i1
         d6 f2 g8 i0 i8

   In this scheme,  lines of constant latitude would run from lower left
to upper right;  meridians would run (mostly) from lower right to upper
left,  with some distortion in the polar zones.

NOTE:  in unrotated form,  tile 00 = north pole = (x=1, y=2);  tile 24 =
(x=0, y=1);  tile 96 = center = (x=2,y=0).  Note north-south symmetry;  if
(x, y) lands in tile n,  (8-x, y) lands in tile n_tiles - 1 - n,  where
n_tiles = 12N^2.  Which is why the code handles only the northern polar
triangles and the equatorial zone;  the southern triangles are handled
just by rotating the coordinates into the northern triangles. */

unsigned xy_to_healpix( const double x, const double y, const unsigned N)
{
   int i, j, line;
   unsigned rval = 0;

   if( y <= -1.)        /* use rotational symmetry for south triangles */
      return( 12 * N * N - 1 - xy_to_healpix( 8. - x, -y, N));
// x -= 1.;       /* translate to top of first polar triangle,  tile 00 */
// y -= 2.;
   i =  (int)floor(  (double)N * (x - y + 1.) / 2.);
   j =  (int)floor( -(double)N * (x + y - 3.) / 2.);
   line = i + j;
   if( (unsigned)line < N)       /* north polar triangles */
      rval = line * (line + 1) * 2 + i % N + (i / N) * (line + 1);
   else
      {                          /* equatorial region */
      const unsigned offset_along_line =
               (unsigned)( i - (line + 1 - N) / 2) % (4 * N);

      rval = N * (1 - N) * 2 + 4 * line * N + offset_along_line;
      }
   return( rval);
}

/* Above is for the 'ring' form of HEALPix,  used in FCCT14.  That's the
only version used in Find_Orb.  For Gaia-DR1,  I ran into the need to handle
the 'nested' HEALPix scheme,  which looks more like this :

                                    3f 3d 37 35
                                    3e 3c 36 34 4e
                                    3b 39 33 31 4b 49
                                    3a 38 32 30 4a 48 42
                        2f 2d 27 25 7f 7d 77 75 bf bd b7 b5
i=0 1  2  3  4  5  6    2e 2c 26 24 7e 7c 76 74 be bc b6 b4
                        2b 29 23 21 7b 79 73 71 bb b9 b3 b1
                        2a 28 22 20 7a 78 72 70 ba b8 b2 b0
            1f 1d 17 15 6f 6d 67 65 af ad a7 a5
            1e 1c 16 14 6e 6c 66 64 ae ac a6 a4
            1b 19 13 11 6b 69 63 61 ab a9 a3 a1      -2     (i, j) rotated to keep cell
            1a 18 12 10 6a 68 62 60 aa a8 a2 a0      -1     00 at the origin
0f 0d 07 05 5f 5d 57 55 9f 9d 97 95                 j=0
0e 0c 06 04 5e 5c 56 54 9e 9c 96 04                   1
0b 09 03 01 5b 59 53 51 9b 99 93 01                   2        N = 4 = power of 2
0a 08 02 00 5a 58 52 50 9a 98 92 00                   3
4f 4d 47 45 8f 8d 87 85                               4
   4c 46 44 8e 8c 86 84
      43 41 8b 89 83 81
         40 8a 88 82 80

   One gets a "tile number",  shifts it up by 2N,  and fits in the position
within the tile with bits interleaved.  That does lose us the rotational
symmetry we previously had.  (Its benefit is that tiles close together in
numbering are a little more likely to also be close together spatially.) */

unsigned xy_to_healpix_nested( const double x, const double y, const unsigned N)
{
   int i, j, line;
   unsigned rval = 0;

// x -= 1.;       /* translate to top of first polar triangle,  tile 00 */
// y -= 2.;
   i =  (int)floor(  (double)N * (x - y + 1.) / 2.);
   j =  (int)floor( -(double)N * (x + y - 3.) / 2.);
   line = i + j;
   if( (unsigned)line < N)       /* north polar triangles */
      rval = line * (line + 1) * 2 + i % N + (i / N) * (line + 1);
   else
      {                          /* equatorial region */
      const unsigned offset_along_line =
               (unsigned)( i - (line + 1 - N) / 2) % (4 * N);

      rval = N * (1 - N) * 2 + 4 * line * N + offset_along_line;
      }
   return( rval);
}

#include <stdio.h>
#include <stdlib.h>

#ifdef OLD_TEST_MAIN

int main( const int argc, const char **argv)
{
   const double ra  = atof( argv[1]) * pi / 180;
   const double dec = atof( argv[2]) * pi / 180;
   const int n_pow = (argc < 3 ? 2 : atoi( argv[3]));
   double x, y;
   unsigned hp;

   ra_dec_to_xy( ra, dec, &x, &y);
   printf( "x = %f  y = %f\n", x, y);
   hp = xy_to_healpix( x, y, 1 << n_pow);
   printf( "hp = %u (%x)\n", hp, hp);
   return( 0);
}
#endif      /* #ifdef OLD_TEST_MAIN */


#ifdef OTHER_TEST_MAIN
   /* Uses a list of HEALPix tile centers provided with the astrometric */
   /* bias data,  and checks that the above code gives the same tile    */
   /* numbers as the file does.  "Noise" can be added to jostle the     */
   /* positions around within the tile.                                 */
int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( "tiles.dat", "rb");
   char buff[100];
   const double noise = (argc > 1 ? atof( argv[1]) : 0) * pi / 180.;

   while( fgets( buff, sizeof( buff), ifile))
      if( *buff != '!')
         {
         unsigned tile_no, computed;
         double ra, dec, x, y;

         sscanf( buff, "%u %lf %lf", &tile_no, &ra, &dec);
         if( noise)
            {
            ra += ((2. * (double)rand( ) / (double)RAND_MAX) - 1.) * noise
                           / cos( dec);
            dec += ((2. * (double)rand( ) / (double)RAND_MAX) - 1.) * noise;
            }
         ra_dec_to_xy( ra, dec, &x, &y);
         computed = xy_to_healpix( x, y, 64);
         if( computed != tile_no)
            printf( "Computed %u\n%s", computed, buff);
         }
   fclose( ifile);
   return( 0);
}
#endif         /* #ifdef OTHER_TEST_MAIN */
