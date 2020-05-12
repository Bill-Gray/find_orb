/* expcalc.cpp: code to compute exposure time,  SNR,  or magnitude
from the other two values

Copyright (C) 2020, Project Pluto

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

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "expcalc.h"

/* A shameless clone of 'expcalc' by Frank Shelly,  leaning
heavily on his code.  Build with 'make expcalc',  which will run

g++ -o expcalc -Wall -Wextra -pedantic -DTEST_CODE expcalc.cpp

   to make the standalone calculator.

   References : http://www2.lowell.edu/rsch/LMI/etc_calc.php
(links to a paper describing underlying theory).  Note that the
formulae given in the paper for 'case 1' (compute exposure time
from SNR and magnitude) have the sign for B and C -- the linear
and constant coefficients of the quadratic to be solved --
backward;  both should be negated.  (Or A should be negated.)

   The Lowell paper makes certain assumptions about the use
of a near-optimal measuring aperture.  (Use too large a measuring
aperture,  and you're gathering lots of background noise and not
adding much signal from the asteroid;  too small an aperture,
and you may toss out some useful signal.)  They use a rather
primitive sky brightness function;  for Find_Orb,  we'll lean
on the Schaefer & Krisciunas model (see 'vislimit.cpp' in the
'lunar' repository).    */

typedef struct
{
   char band;
   double extinction;            /* mags per airmass */
   double zero_point;
} filter_t;

/* The 'zero point' is the counts from a mag 0 star in this band
per cm^2 per second.  Multiply by the collecting area of the
telescope -- don't forget to omit the central obstruction --
and you get the counts per second for your telescope;  divide
by 100^(mag/5) to get the counts per second for a star of a
given magnitude.   */

const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923;

/* A Gaussian with a particular full-width-half-maximum (FWHM) 'fwhm'
will drop to half the maximum value at a radius fwhm/2.  Since the
intensity runs as I = I0 * exp( -(r/sigma)^2 / 2),  we can say :

I0/2 = I0 * exp( -(.5 * fwhm / sigma)^2 / 2)
1/2 = exp( -(.5 * fwhm / sigma)^2 / 2)
-ln(2) = -(.5 * fwhm / sigma)^2 / 2
2 * ln(2) = (.5 * fwhm / sigma)^2
fwhm / sigma = 2 * sqrt( 2 * ln(2)) = 2.354820045030949 */

static double fraction_inside( const expcalc_config_t *c)
{
   const double full_widths_per_sigma = 2.354820045030949;
   const double r_scaled = full_widths_per_sigma * c->aperture / c->fwhm;

   return( 1. - exp( -r_scaled * r_scaled / 2.));
}

/* Sizes in centimeters */
static double effective_area( const expcalc_config_t *c)
{
   return( pi * (c->primary_diam * c->primary_diam - c->obstruction_diam * c->obstruction_diam) / 4.);
}

static double sky_electrons_per_second_per_pixel( const expcalc_config_t *c, const double zero_point)
{
   const double area = effective_area( c);
   const double sky_electrons_per_sec_per_square_arcsec =
                   pow( 10., -0.4 * c->sky_brightness) * zero_point * area * c->qe;

   return( sky_electrons_per_sec_per_square_arcsec * c->pixel_size * c->pixel_size);
}

typedef struct
{
   double extinction, zero_point;
   double n_pixels_in_aperture, s, noise2;
   double area;
} expcalc_internals_t;

static double star_electrons_per_second_per_pixel( const expcalc_config_t *c, const double mag,
              const expcalc_internals_t *e)
{
   const double mag_corr = mag + c->airmass * e->extinction;
   const double rval = pow( 10., -0.4 * mag_corr) * e->zero_point * effective_area( c)
                                    * c->qe * fraction_inside( c);
   return( rval);
}

static double solve_quadratic( const double a, const double b, const double c)
{
   const double discr = b * b - 4. * a * c;

   assert( discr >= 0.);
   return( .5 * (sqrt( discr) - b) / a);
}

#define EXPCALC_UNKNOWN_FILTER         -1
#define EXPCALC_FAILED        (-99.)

static int find_filter( expcalc_internals_t *e, const char filter)
{
   size_t i;
   static filter_t filters[] = {
     { 'U', 0.60, 5.50e+05 },
     { 'B', 0.40, 3.91e+05 },
     { 'V', 0.20, 8.66e+05 },
     { 'R', 0.10, 1.10e+06 },
     { 'I', 0.08, 6.75e+05 },
     { 'N', 0.20, 4.32e+06 },
     { 'W', 0.15, 2.00e+06 } };

   for( i = 0; i < sizeof( filters) / sizeof( filters[0]); i++)
      if( filter == filters[i].band)
         {
         e->zero_point = filters[i].zero_point;
         e->extinction = filters[i].extinction;
         return( 0);
         }
   return( -1);
}

/* For codes where we lack real information,  we'll return the first
entry in the table below : a one-meter scope similar (except in size) to
(703).  It won't be right,  of course,  but may allow comparison of what
you'd get at times t1 and t2,  or for objects A and B.  For some scopes,
we can get primary diameters from 'scopes.txt' or 'details.txt'. */

int find_expcalc_config_from_mpc_code( const char *mpc_code, expcalc_config_t *c)
{
   size_t i = 0, rval = 0;
   static expcalc_config_t configs[] = {
    /*   Code  Filt PrDi ObsDi Ape FWHM QE Rea  Pix   Sky  Airmass */
       { "500", 'N', 100., 25., 6., 3., .9, 8., 3.,    20., 1.5 },
       { "703", 'N',  72., 25., 6., 3., .9, 8., 3.,    20., 1.5 },
       { "E12", 'N',  50., 23., 6., 3., .9, 8., 1.8,   20., 1.5 },
       { "G96", 'N', 152., 48., 6., 3., .9, 8., 1.5,   20., 1.5 },
       { "I52", 'N', 100., 40., 6., 3., .9, 8., 1.036, 20., 1.5 },
       { "V06", 'N', 154., 40., 6., 3., .9, 8., 0.572, 20., 1.5 } };
   const size_t n_configs = sizeof( configs) / sizeof( configs[0]);

   while( i < n_configs && strcmp( configs[i].mpc_code, mpc_code))
      i++;
   if( i == n_configs)
      {
      rval = -1;
      c = 0;
      }
   if( c)
      *c = configs[i];
   return( rval);
}

static int set_internals( expcalc_internals_t *e, const expcalc_config_t *c)
{
   if( find_filter( e, c->filter))
      return( EXPCALC_UNKNOWN_FILTER);
   e->n_pixels_in_aperture = pi * c->aperture * c->aperture
                  / (c->pixel_size * c->pixel_size);
   e->noise2 = c->readnoise * c->readnoise;
   e->s = sky_electrons_per_second_per_pixel( c, e->zero_point);
   return( 0);
}

double mag_from_snr_and_exposure( const expcalc_config_t *c,
                              const double snr, const double exposure)
{
   const double tval = snr * snr * exposure;
   double n_star, mag_exp, mag;
   expcalc_internals_t e;

   if( set_internals( &e, c))
      return( EXPCALC_FAILED);
   n_star = solve_quadratic( exposure * exposure, -tval,
                    e.n_pixels_in_aperture * (e.noise2 - tval * e.s));
   mag_exp = n_star / (e.zero_point * effective_area( c) * c->qe * fraction_inside( c));
   mag = -2.5 * log10( mag_exp);
   mag -=  c->airmass * e.extinction;
   return( mag);
}

double snr_from_mag_and_exposure( const expcalc_config_t *c,
                              const double mag, const double exposure)
{
   expcalc_internals_t e;
   double n_star, signal, noise;

   if( set_internals( &e, c))
      return( EXPCALC_FAILED);
   n_star = star_electrons_per_second_per_pixel( c, mag, &e);
   signal = n_star * exposure;
   noise = sqrt( signal
             + e.n_pixels_in_aperture * (e.s * exposure + e.noise2));
   return( signal / noise);
}

double exposure_from_snr_and_mag( const expcalc_config_t *c,
                              const double snr, const double mag)
{
   expcalc_internals_t e;
   double n_star, exposure;

   if( set_internals( &e, c))
      return( EXPCALC_FAILED);
   n_star = star_electrons_per_second_per_pixel( c, mag, &e);
   exposure = solve_quadratic( n_star * n_star,
                        -snr * snr * (n_star + e.n_pixels_in_aperture * e.s),
                        -snr * snr * e.n_pixels_in_aperture * e.noise2);
   return( exposure);
}

#ifdef TEST_CODE

const char *usage_statement =
 "Usage:  /home/observer/bin/expcalc <switches>\n"
 "              SNR <exposure (sec)> <magnitude> |\n"
 "              mag <SNR> <exposure (sec)> |\n"
 "              exp <SNR> <magnitude>\n"
 " <switches>\n"
 "   -mpc <code>      # code to look up configuration values for (default I52)\n"
 "   -debug           # Output more verbose info to help with debugging\n"
 "   -filter <N|U|B|V|R|I|W>         # Change filter code N=none default N\n"
 "   -primary <diameter cm>          # Set primary aperature to this diameter in cm\n"
 "   -radius  <phot radius arcsc>    # Set pinhole photometry radius to this in arcsec\n"
 "   -qe <quantum efficiency>        # fractional quantum efficiency of the detector\n"
 "   -readnoise <electrons>          # detector readout noise in electrons\n"
 "   -pixelsize <arcsec>             # assume square pixel size in arcsec\n"
 "   -skybrightness <mag/arcsec/sec> # Magnitude of square arcsec sky in one second\n"
 "   -airmass <fraction>             # airmass of object\n"
 "   -fwhm <arcsec>                  # Full width half maximum in arcsec\n"
 "   -obstdiameter <cm>              # diameter in cm of central obstruction\n";

static void usage( void)
{
   fprintf( stderr, "%s", usage_statement);
   exit( -1);
}

/* Some test cases :

./expcalc -mpc V06 -primary 300 snr 60 22
SNR 4.54 with exposure time 60 seconds and magnitude 22.00

./expcalc snr 30 20.5
SNR 3.92 with exposure time 30 seconds and magnitude 20.50

./expcalc -skybrightness 19 -primary 200 -filter V snr 30 20.5
SNR 2.38 with exposure time 30 seconds and magnitude 20.50
*/

int main( const int argc, const char **argv)
{
   expcalc_config_t c;
   double mag = 0., exposure = 0., snr = 0.;
   int i;
   const char *mpc_code = "I52";
   bool debug = false;

   for( i = 1; i < argc; i++)
      if( !strcmp( argv[i], "-mpc"))
          mpc_code = argv[i + 1];
   if( find_expcalc_config_from_mpc_code( mpc_code, &c))
      {
      fprintf( stderr, "Unrecognized MPC code '%s'\n", mpc_code);
      fprintf( stderr, "At present,  this program only knows about the following codes:\n"
                       "(703) (E12) (I52) (G96) (V06)\n");
      usage( );
      }
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-' && argv[i][1])
         {
         switch( argv[i][1])
            {
            case 'a':
               c.airmass = atof( argv[i + 1]);
               break;
            case 'd':
               debug = true;
               break;
            case 'f':
               if( argv[i][2] == 'i')
                  c.filter = argv[i + 1][0];
               if( argv[i][2] == 'w')
                  c.fwhm = atof( argv[i + 1]);
               break;
            case 'o':
               c.obstruction_diam = atof( argv[i + 1]);
               break;
            case 'p':
               if( argv[i][2] == 'r')     /* primary */
                  c.primary_diam = atof( argv[i + 1]);
               if( argv[i][2] == 'i')     /* pixel size */
                  c.pixel_size = atof( argv[i + 1]);
               break;
            case 'q':
               c.qe = atof( argv[i + 1]);
               break;
            case 'r':
               c.aperture = atof( argv[i + 1]);
               break;
            case 's':
               c.sky_brightness = atof( argv[i + 1]);
               break;
            }
         }
   for( i = 1; i < argc - 2; i++)
      if( !strcmp( argv[i], "snr") || !strcmp( argv[i], "SNR"))
         {
         exposure = atof( argv[i + 1]);
         mag = atof( argv[i + 2]);
         }
      else if( !strcmp( argv[i], "exp"))
         {
         snr = atof( argv[i + 1]);
         mag = atof( argv[i + 2]);
         }
      else if( !strcmp( argv[i], "mag"))
         {
         snr = atof( argv[i + 1]);
         exposure = atof( argv[i + 2]);
         }

    if( debug)
      {
      printf( "%.2f-cm primary, %.2f-cm obstruction\n",
                  c.primary_diam, c.obstruction_diam);
      printf( "Filter %c, QE %.2f, read noise %.2f electrons/pixel\n",
                  c.filter, c.qe, c.readnoise);
      printf( "Pixels are %.2f arcsec;  aperture %.2f pixels, FWHM %.2f pixels\n",
                  c.pixel_size, c.aperture, c.fwhm);
      printf( "Sky brightness %.2f mag/arcsec^2; airmass %.2f\n",
                  c.sky_brightness, c.airmass);
      }

   if( mag && exposure)
      snr = snr_from_mag_and_exposure( &c, mag, exposure);
   else if( mag && snr)
      exposure = exposure_from_snr_and_mag( &c, snr, mag);
   else if( snr && exposure)
      mag = mag_from_snr_and_exposure( &c, snr, exposure);
   else
      usage( );
   printf( "SNR %.2f with exposure time %.0f seconds and magnitude %.2f\n",
                  snr, exposure, mag);
   return( 0);
}
#endif         /* #ifdef TEST_CODE */
