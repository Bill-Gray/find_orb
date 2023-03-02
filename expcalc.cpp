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
#include <ctype.h>
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

   Michael Richmond provides some useful background in his
class notes at

http://spiff.rit.edu/classes/ast613/lectures/signal/signal_illus.html
https://www.noao.edu/kpno/manuals/dim/#ccdtime

   The Lowell paper uses a somewhat primitive sky brightness function; for
Find_Orb,  we'll lean on the Schaefer & Krisciunas model (see 'vislimit.cpp'
in the 'lunar' repository).  See note below on aperture size.  */

typedef struct
{
   char band;
   double extinction;            /* mags per airmass */
   double zero_point;
} filter_t;

/* 'zero_point = photons per square cm per second for a mag 0 star
with a Vega-like spectrum,  as measured above the atmosphere.  See
http://spiff.rit.edu/richmond/wiyn/technotes/signal_wiyn09.cgi */

/* The Lowell paper cited above points out that the theoretical
optimal pinhole photometry aperture should be about 0.67 times
the full width at half maximum (fwhm).  More than that,  and
you're bringing in outlying pixels with less signal intensity
and your SNR will suffer.  Less than that,  and you're excluding
some pixels whose signals could make your SNR better.

   The Lowell paper also suggests "a minimum of nine pixels,
for a minimum radius of 2.5 pixels".  There's a contradiction
here,  since sqrt(9/pi) = 1.693.  I _think_ that's what they
really meant,  and that's the minimum aperture implemented below.
(But temporarily disabled;  we're using the user-specified aperture
for the nonce.)

   Finally,  after some discussions with CSS,  we agreed that the
'fwhm' for a site should be assumed to be that at the zenith,
and should be multiplied by the airmass (i.e.,  the seeing gets
worse near the horizon).  Must admit I'm not totally certain of
that;  a better understanding of how the PSF spreads as a function
of airmass would help.  Presumably,  for small scopes,  the FWHM
is just going to be proportional to wavelength/scope diameter;
you need to get beyond some point for the airmass to matter.
*/

#ifdef NOT_CURRENTLY_IN_USE
static double _optimal_aperture( const expcalc_config_t *c)
{
   double rval = c->fwhm * c->airmass * 0.67;
   const double min_aperture = c->pixel_size * 1.693;

   return( rval > min_aperture ? rval : min_aperture);
}
#endif

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
fwhm / sigma = 2 * sqrt( 2 * ln(2)) = 2.354820045030949

NOTE (2021 Nov 2) : After some discussion with CSS folks,  we agreed
that the 'fwhm' is proportional to airmass.  The value given for a
particular site is assumed to be for the zenith,  and should be
multiplied by the airmass.  (Or -- I think more likely -- the square
root of the airmass.  Due to some uncertainty about all of this,  we have
reverted this change and are _not_ adjusting for airmass here.)   */

static double _fraction_inside( const expcalc_config_t *c)
{
   const double full_widths_per_sigma = 2.354820045030949;
/* const double real_fwhm = c->fwhm * c->airmass;     */
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
   double n_star, area;
} expcalc_internals_t;

static double star_electrons_per_second_per_pixel( const expcalc_config_t *c, const double mag,
              const expcalc_internals_t *e)
{
   const double mag_corr = mag + c->airmass * e->extinction;
   const double rval = pow( 10., -0.4 * mag_corr) * e->zero_point * effective_area( c)
                                    * c->qe * _fraction_inside( c);
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

/* Both the exposure calculation code (see 'expcalc.cpp') and sky
brightness code (see 'vislimit.cpp' in the 'lunar' repository) can handle
UBVRI photometry.  For other bands,  we translate to the nearest band.
At some point,  a smarter solution may be used in which we interpolate or
extrapolate between bands.  But I expect that to be quite challenging.
This should be good enough for most purposes;  in practice,  neither
exposure times nor sky brightness can be computed accurately enough to
notice the difference.  Note that the return value is an _index_,
0=U, 1=B, 2=V, 3=R, 4=I, -1=unknown.  Also note that the exposure
calculator also knows about N and W;  those shouldn't be translated.
*/

int xlate_filter_to_ubvri( const char filter)
{
   int i;
   const char *groups[5] = { "Uu",           /* filters closest to U */
                             "B",            /* filters closest to B */
                             "Vg",           /* filters closest to V */
                             "RrWwcoN",      /* filters closest to R */
                             "Iizy" };       /* filters closest to I */

   for( i = 0; i < 5; i++)
      if( strchr( groups[i], filter))
         return( i);
   return( -1);
}

static int find_filter( expcalc_internals_t *e, char filter)
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

   if( filter != 'N' && filter != 'W')
      {
      const char *filts = "UBVRI";
      int idx = xlate_filter_to_ubvri( filter);

      if( idx == -1)    /* unrecognized filter */
         filter = 'R';
      else
         filter = filts[idx];
      }
   for( i = 0; i < sizeof( filters) / sizeof( filters[0]); i++)
      if( filter == filters[i].band)
         {
         e->zero_point = filters[i].zero_point;
         e->extinction = filters[i].extinction;
         return( 0);
         }
   assert( 0);       /* we _have_ to find a filter */
   return( -1);
}

static int set_internals( expcalc_internals_t *e, const expcalc_config_t *c)
{
   if( find_filter( e, c->filter))
      return( EXPCALC_UNKNOWN_FILTER);
   assert( c->aperture > 0.1 && c->aperture < 100.);
   assert( c->pixel_size > 0.01 && c->pixel_size < 100.);
   assert( c->readnoise > 0.);
   assert( c->primary_diam > 1.0 && c->primary_diam < 3100.);
   assert( c->obstruction_diam >= 0.0 && c->obstruction_diam < c->primary_diam);
   assert( c->qe > 0.01 && c->qe <= 1.0);
   e->n_pixels_in_aperture = pi * c->aperture * c->aperture
                  / (c->pixel_size * c->pixel_size);
   e->noise2 = c->readnoise * c->readnoise;
   e->s = sky_electrons_per_second_per_pixel( c, e->zero_point);
   return( 0);
}

static double internal_mag_from_snr_and_exposure(
                              expcalc_internals_t *e,
                              const expcalc_config_t *c,
                              const double snr, const double exposure)
{
   const double tval = snr * snr * exposure;
   double mag_exp, mag;

   if( set_internals( e, c))
      return( EXPCALC_FAILED);
   e->n_star = solve_quadratic( exposure * exposure, -tval,
                    e->n_pixels_in_aperture * (e->noise2 - tval * e->s));
   mag_exp = e->n_star / (e->zero_point * effective_area( c) * c->qe * _fraction_inside( c));
   mag = -2.5 * log10( mag_exp);
   mag -=  c->airmass * e->extinction;
   return( mag);
}

double mag_from_snr_and_exposure( const expcalc_config_t *c,
                              const double snr, const double exposure)
{
   expcalc_internals_t e;

   return( internal_mag_from_snr_and_exposure( &e, c, snr, exposure));
}

static double internal_snr_from_mag_and_exposure(
                              expcalc_internals_t *e,
                              const expcalc_config_t *c,
                              const double mag, const double exposure)
{

   double signal, noise;

   if( set_internals( e, c))
      return( EXPCALC_FAILED);
   e->n_star = star_electrons_per_second_per_pixel( c, mag, e);
   signal = e->n_star * exposure;
   noise = sqrt( signal
             + e->n_pixels_in_aperture * (e->s * exposure + e->noise2));
   return( signal / noise);
}

double snr_from_mag_and_exposure( const expcalc_config_t *c,
                              const double mag, const double exposure)
{
   expcalc_internals_t e;

   return( internal_snr_from_mag_and_exposure( &e, c, mag, exposure));
}

static double internal_exposure_from_snr_and_mag(
                              expcalc_internals_t *e,
                              const expcalc_config_t *c,
                              const double snr, const double mag)
{
   double exposure;

   if( set_internals( e, c))
      return( EXPCALC_FAILED);
   e->n_star = star_electrons_per_second_per_pixel( c, mag, e);
   exposure = solve_quadratic( e->n_star * e->n_star,
                        -snr * snr * (e->n_star + e->n_pixels_in_aperture * e->s),
                        -snr * snr * e->n_pixels_in_aperture * e->noise2);
   return( exposure);
}

double exposure_from_snr_and_mag( const expcalc_config_t *c,
                              const double snr, const double mag)
{
   expcalc_internals_t e;

   return( internal_exposure_from_snr_and_mag( &e, c, snr, mag));
}

/* See 'scope.json' and 'site_310.txt' for two different ways of storing
site data.  The former is preferred by CSS for automated processing;
'site_(mpc code).txt' files are a lot easier for humans to edit.  Both
let you do the same things. */

static const char *get_config( const char *buff, const char *config_var)
{
   const size_t len = strlen( config_var);

   if( !memcmp( buff, config_var, len) && buff[len] == ':'
                     && buff[len + 1] == ' ')
      buff += len + 2;
   else                                /* look for JSON version */
      {
      char search_str[40];

      snprintf( search_str, sizeof( search_str), "\"%s\":", config_var);
      buff = strstr( buff, search_str);
      if( buff)
         {
         buff += strlen( search_str);
         while( *buff == ' ')
            buff++;
         if( *buff == '"')
            buff++;
         }
      }
   return( buff);
}

static void set_config_double( const char *buff, const char *config_var,
            double *value)
{
   buff = get_config( buff, config_var);
   if( buff)
      *value = atof( buff);
}

#define IS_POWER_OF_TWO( n)    (((n) & ((n)-1)) == 0)

#define EXPCALC_NO_CONFIG_FOUND             -1
#define EXPCALC_GEOCENTRIC_CONFIG            0
#define EXPCALC_SITE_SPECIFIC_CONFIG         1

int find_expcalc_config_from_mpc_code( const char *mpc_code,
             FILE *ifile, expcalc_config_t *c)
{
   int rval = EXPCALC_NO_CONFIG_FOUND;
   char buff[200];
   const char *tptr;
   bool getting_data = false;

   c->airmass = 1.5;
   while( fgets( buff, sizeof( buff), ifile))
      if( get_config( buff, "500"))
         {
         rval = EXPCALC_GEOCENTRIC_CONFIG;
         getting_data = true;
         }
      else if( get_config( buff, mpc_code))
         {
         rval = EXPCALC_SITE_SPECIFIC_CONFIG;
         getting_data = true;
         }
      else if( !memcmp( buff, "Site: ", 6) && !memcmp( buff + 6, mpc_code, 3))
         {              /* 'site_(code).txt' format */
         rval = EXPCALC_SITE_SPECIFIC_CONFIG;
         getting_data = true;
         }
      else if( strchr( buff, '}'))
         getting_data = false;
      else if( getting_data)
         {
         if( (tptr = get_config( buff, "Filter")) != NULL)
            c->filter = *tptr;
         set_config_double( buff, "Primary", &c->primary_diam);
         set_config_double( buff, "Obstruction", &c->obstruction_diam);
         set_config_double( buff, "Aperture", &c->aperture);
         set_config_double( buff, "FWHM", &c->fwhm);
         set_config_double( buff, "QE", &c->qe);
         set_config_double( buff, "ReadNoise", &c->readnoise);
         set_config_double( buff, "PixelSize", &c->pixel_size);
         set_config_double( buff, "SkyBrightness", &c->sky_brightness_at_zenith);
         set_config_double( buff, "MinAlt", &c->min_alt);
         set_config_double( buff, "MaxAlt", &c->max_alt);
         set_config_double( buff, "MinDec", &c->min_dec);
         set_config_double( buff, "MaxDec", &c->max_dec);
         set_config_double( buff, "MinHA", &c->min_ha);
         set_config_double( buff, "MaxHA", &c->max_ha);
         set_config_double( buff, "MinElong", &c->min_elong);
         set_config_double( buff, "MaxElong", &c->max_elong);
         if( strstr( buff, "\"Horizon\""))      /* JSON-style horizon */
            {
            int i, n_found = 0;

            while( !strchr( buff, 0x5d) && fgets( buff, sizeof( buff), ifile))
               {
               char *endptr;

               for( i = 0; buff[i]; i++)
                  if( isdigit( buff[i]))
                     {
                     n_found++;
                     if( IS_POWER_OF_TWO( n_found))
                        c->horizon = (double *)realloc( c->horizon, n_found * 2 * sizeof( double));
                     c->horizon[n_found - 1] = strtod( buff + i, &endptr);
                     i = (int)( endptr - buff) - 1;
                     }
               }
            assert( n_found > 2);
            assert( !(n_found & 1));
            c->n_horizon_points = n_found / 2;
            }
         if( !memcmp( buff, "Horizon start", 13))     /* text-file horizon */
            {                                        /* see 'site_310.txt' */
            int n_read, n_found = 0;                 /* for explanation    */

            while( fgets( buff, sizeof( buff), ifile) &&
                              memcmp( buff, "Horizon", 7))
               if( !memcmp( buff, "Point: ", 7))
                  {
                  n_found++;
                  if( IS_POWER_OF_TWO( n_found))
                     c->horizon = (double *)realloc( c->horizon, n_found * 4 * sizeof( double));
                  n_read = sscanf( buff + 7, "%lf, %lf",
                        c->horizon + n_found * 2, c->horizon + n_found * 2 + 1);
                  assert( n_read == 2);
                  }
            assert( n_found > 1);
            c->n_horizon_points = n_found;
            }
         c->sky_brightness = c->sky_brightness_at_zenith;
         }
   return( rval);
}

void free_expcalc_config_t( expcalc_config_t *c)
{
   if( c->horizon)
      free( c->horizon);
   c->horizon = NULL;
}

static int under_horizon_slice( const double alt, const double az,
         const double *p1, const double *p2)
{
   const double x1 = fmod( az - p1[0] + 900., 360.) - 180.;
   const double x2 = fmod( az - p2[0] + 900., 360.) - 180.;
   int rval = 0;

   if( x1 * x2 <= 0. && fabs( x1 - x2) < 180.)
      {
      const double y = (p1[1] * x2 - p2[1] * x1) / (x2 - x1);

      if( y > alt)
         rval = 1;
      }
   return( rval);
}

int is_under_horizon( const double alt, const double az,
                              const expcalc_config_t *c)
{
   int i, rval = (alt < c->min_alt || alt > c->max_alt);

   if( c->horizon && !rval)
      {
      for( i = 0; i < c->n_horizon_points; i++)
         rval ^= under_horizon_slice( alt, az,
                     c->horizon + 2 * i,
                     c->horizon + 2 * ((i + 1) % c->n_horizon_points));
      }
   return( rval);
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
 "   -scope <filename>               # Reset 'scope.json' filename\n"
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
   expcalc_internals_t e;
   double mag = 0., exposure = 0., snr = 0.;
   int i;
   const char *mpc_code = "I52";
   bool debug = false;
   const char *scope_json_filename = "scope.json";
   FILE *ifile;

   for( i = 1; i < argc; i++)
      if( !strcmp( argv[i], "-mpc"))
          mpc_code = argv[i + 1];
      else if( !memcmp( argv[i], "-sc", 3))
          scope_json_filename = argv[i + 1];
   ifile = fopen( scope_json_filename, "rb");
   if( !ifile)       /* look for ~/.find_orb/scope.json */
      {
      char buff[255], *home_dir = getenv( "HOME");

      if( home_dir)
         {
         strcpy( buff, home_dir);
         strcat( buff, "/.find_orb/");
         strcat( buff, scope_json_filename);
         ifile = fopen( buff, "rb");
         }
      }
   if( !ifile)
      {
      fprintf( stderr, "Couldn't open 'scope.json'\n");
      usage( );
      }
   memset( &c, 0, sizeof( expcalc_config_t));
   switch( find_expcalc_config_from_mpc_code( mpc_code, ifile, &c))
      {
      case EXPCALC_NO_CONFIG_FOUND:
         fprintf( stderr, "No default details (not supposed to happen)\n");
         usage( );
         break;
      case EXPCALC_GEOCENTRIC_CONFIG:
         fprintf( stderr, "No details for MPC code '%s'\n", mpc_code);
         fprintf( stderr, "Using default 'geocentric' values\n");
         fprintf( stderr, "Check ~/.find_orb/scope.json to see if (%s) is listed.\n",
                       mpc_code);
         break;
      case EXPCALC_SITE_SPECIFIC_CONFIG:
         break;               /* got the details we wanted for the site */
      }
   fclose( ifile);
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
               if( argv[i][2] == 'i')   /* specify filter, <N|U|B|V|R|I|W> */
                  c.filter = argv[i + 1][0];
               else              /* full width half maximum, in arcsec */
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
               if( argv[i][2] == 'e')        /* read noise, in electrons */
                  c.readnoise = atof( argv[i + 1]);
               else              /* pinhole photometry radius in arcsec */
                  c.aperture = atof( argv[i + 1]);
               break;
            case 's':
               if( argv[i][2] != 'c')     /* resetting scope.json is handled above */
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
      printf( "Pixels are %.2f arcsec;  aperture %.2f arcsec, FWHM %.2f arcsec\n",
                  c.pixel_size, c.aperture, c.fwhm);
      printf( "Sky brightness %.2f mag/arcsec^2; airmass %.2f\n",
                  c.sky_brightness, c.airmass);
      }

   if( mag && exposure)
      snr = internal_snr_from_mag_and_exposure( &e, &c, mag, exposure);
   else if( mag && snr)
      exposure = internal_exposure_from_snr_and_mag( &e, &c, snr, mag);
   else if( snr && exposure)
      mag = internal_mag_from_snr_and_exposure( &e, &c, snr, exposure);
   else
      usage( );
   printf( "SNR %.2f with exposure time %.1f seconds and magnitude %.2f\n",
                  snr, exposure, mag);
   if( debug)
      {
      const double signal = e.n_star * exposure;
      const double sky_count = e.s * exposure * e.n_pixels_in_aperture;
      const double read_noise = e.noise2 * e.n_pixels_in_aperture;

      printf( "%.3f star electrons (total exposure)\n", signal);
      printf( "Noise from star %.3f (square root of the above line)\n", sqrt( signal));
      printf( "Noise from sky %.3f\n", sqrt( sky_count));
      printf( "Read noise %.3f (from %.3f pixels)\n", sqrt( read_noise), e.n_pixels_in_aperture);
      printf( "Total noise %.3f (above three lines added in quadrature)\n",
               sqrt( signal + sky_count + read_noise));
      }
   free_expcalc_config_t( &c);
   return( 0);
}
#endif         /* #ifdef TEST_CODE */
