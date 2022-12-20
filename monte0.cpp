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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "watdefs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "stringex.h"
#include "afuncs.h"
#include "monte0.h"

double gaussian_random( void);                           /* monte0.c */
int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
double get_planet_mass( const int planet_idx);                /* orb_func.c */
void remove_insignificant_digits( char *tbuff);          /* monte0.c */
void set_up_observation( OBSERVE FAR *obs);                 /* mpc_obs.c */
void set_obs_vect( OBSERVE FAR *obs);        /* mpc_obs.h */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

void put_orbital_elements_in_array_form( const ELEMENTS *elem,
                  double *output_array)
{
   output_array[MONTE_TP] = elem->perih_time;
   output_array[MONTE_ECC] = elem->ecc;
   output_array[MONTE_q] = elem->major_axis * (1. - elem->ecc);
   output_array[MONTE_Q] = elem->major_axis * (1. + elem->ecc);
   output_array[MONTE_INV_A] = 1. / elem->major_axis;
   output_array[MONTE_INCL] = elem->incl * 180. / PI;
   output_array[MONTE_MEAN_ANOM] = elem->mean_anomaly * 180. / PI;
   output_array[MONTE_ARG_PER] = elem->arg_per * 180. / PI;
   output_array[MONTE_ASC_NODE] = elem->asc_node * 180. / PI;
   output_array[MONTE_EARTH_MOID] = 0.;
   output_array[MONTE_H] = elem->abs_mag;
}

void add_monte_orbit( double *monte_data, const ELEMENTS *elem,
                  const int n_orbits)
{
   double tarr[MONTE_N_ENTRIES];
   double *offsets = monte_data + 2 * MONTE_N_ENTRIES;
   int i;

   put_orbital_elements_in_array_form( elem, tarr);
   if( !n_orbits)                /* initializing step */
      for( i = 0; i < MONTE_N_ENTRIES; i++)
         {
         offsets[i] = tarr[i];
         monte_data[i] = monte_data[i + MONTE_N_ENTRIES] = 0.;
         }
   else
      {
      for( i = 0; i < MONTE_N_ENTRIES; i++)
         {
         double delta = tarr[i] - offsets[i];

         if( i >= MONTE_INCL && i <= MONTE_ASC_NODE)
            {
            if( delta > 180.)            /* keep angular arguments in the */
               delta -= 360.;            /* proper range,  not wrapping   */
            else if( delta < -180.)      /* around at +/- 180 degrees     */
               delta += 360.;
            }
         monte_data[i] += delta;
         monte_data[i + MONTE_N_ENTRIES] += delta * delta;
         }
      }
}

void compute_monte_sigmas( double *sigmas, const double *monte_data,
                  const int n_orbits)
{
   int i;

   for( i = 0; i < MONTE_N_ENTRIES; i++)
      {
      const double avg_square = monte_data[i + MONTE_N_ENTRIES] / (double)n_orbits;
      const double avg_value = monte_data[i] / (double)n_orbits;

      sigmas[i] = sqrt( avg_square - avg_value * avg_value);
      }
}

static double *store_ra_decs_mags_times( unsigned n_obs, const OBSERVE *obs)
{
   double *stored_ra_decs = (double *)calloc( 4 * n_obs, sizeof( double));
   double *tptr = stored_ra_decs;

   assert( stored_ra_decs);
   if( !stored_ra_decs)
      return( NULL);
   while( n_obs--)
      {
      *tptr++ = obs->ra;
      *tptr++ = obs->dec;
      *tptr++ = obs->obs_mag;
      *tptr++ = obs->jd;
      obs++;
      }
   return( stored_ra_decs);
}

void restore_ra_decs_mags_times( unsigned n_obs, OBSERVE *obs,
                           const double *stored_ra_decs)
{
   const double *tptr = stored_ra_decs;

   assert( tptr);
   while( n_obs--)
      {
      obs->ra = *tptr++;
      obs->dec = *tptr++;
      obs->obs_mag = *tptr++;
      obs->jd = *tptr++;
      if( obs->note2 == 'S')
         set_obs_vect( obs);
      else
         set_up_observation( obs);
      obs++;
      }
}

#include <stdint.h>

/* Defining 64-bit constants portably and avoiding nuisance warnings
is rather difficult to arrange,  but can be done. */

#ifndef UINT64_C
   #ifdef _MSC_VER
      #define UINT64_C( a) (a##ui64)
   #else
      #ifdef _WIN32
         #define UINT64_C( a) (a##ULL)
      #else
         #define UINT64_C( a) ((uint64_t)(a##UL))
      #endif
   #endif
#endif

#ifdef __clang__
#pragma GCC diagnostic ignored "-Wc++11-long-long"
#endif     /* suppress nuisance warning */

/* Mostly cut & pasted from http://www.pcg-random.org/download.html */
/* Permuted Congruential Generator */

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

static uint32_t pcg32_random_r( pcg32_random_t* rng)
{
    const uint32_t xorshifted = (uint32_t)( ((rng->state >> 18u) ^ rng->state) >> 27u);
    const int rot = (int)( rng->state >> 59u);
    const uint64_t multiplier = UINT64_C( 6364136223846793005);

    /* Advance internal state                                         */
    rng->state = rng->state * multiplier + (rng->inc | 1);
    /* Calculate output function (XSH RR), uses old state for max ILP */
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

static uint64_t pcg_64_bits( pcg32_random_t* rng)
{
   return( (uint64_t)pcg32_random_r( rng)
          | ((uint64_t)pcg32_random_r( rng) << 32));
}

/* Returns a uniform random variable,  0 <= rval < 1. */

double uniform_random( const int free_up)
{
   static pcg32_random_t rng = { 314159265, 358979323 };
   const double two_to_the_63rd_power = 9223372036854775808.;

   if( free_up)         /* flag to free up memory */
      return( 0.);
   else
      return( (double)( pcg_64_bits( &rng) >> 1) / two_to_the_63rd_power);
}

/* The Box-Muller transform converts two uniformly-distributed random
variables into two Gaussian-distributed random variables.  We save one
for a subsequent call.  See

https://en.wikipedia.org/wiki/Box-Muller_transform

   Used below to add Gaussian noise to positions,  magnitudes,  and
times for Monte Carlo and statistical ranging. */

double gaussian_random( void)
{
   static double saved = 0.;
   static bool have_a_value_already = false;
   double rval;

   if( have_a_value_already)
      rval = saved;
   else
      {
      const double rt = log( 1. - uniform_random( 0));
      const double r = sqrt( -2. * rt);
      const double theta = 2. * PI * uniform_random( 0);

      rval = r * cos( theta);
      saved = r * sin( theta);
      }
   have_a_value_already = !have_a_value_already;
   return( rval);
}

/* Add some Gaussian noise to each RA/dec value,  magnitude,  and time.
Two Gaussian-distributed random numbers are generated using the
Box-Muller transform,  scaled according to the observation sigma;  this
is  then added to the RA and dec. the noise_in_arcseconds.  Then two
more Gaussians are computed and used to add noise to the magnitude and time
of observation.  The original RA/decs/mags/times are stored in an array;
calling 'restore_ra_decs_mags_times()' restores them, removing the noise.

   2012 Jul 22:  Marco Micheli pointed out that 'r' should be scaled by
the observation sigma.  It now is.  And yes,  it should have been that
way right from the beginning...

   2012 Feb 9:  switched to use of MT64 (Mersenne Twister,  64-bit version)
for generating pseudo-random numbers.  Default C-library PRNGs are sometimes
adequate and sometimes horrible.  MT64 is quite well thought of,  and should
enable us to assume that any problems are _not_ due to insufficiently random
numbers.  Also, if one uses the C-library PRNG,  you'll get different results
on different systems with different implementations.  */

double *add_gaussian_noise_to_obs( int n_obs, OBSERVE *obs,
                 const double noise_in_sigmas)
{
   const double noise_in_radians = noise_in_sigmas * PI / (180. * 3600.);
   double *rval;

   if( !obs)         /* flag to free up memory */
      {
      uniform_random( 1);
      return( NULL);
      }
   rval = store_ra_decs_mags_times( n_obs, obs);
   while( n_obs--)
      {
      const double x = gaussian_random( );
      const double y = gaussian_random( );

      obs->ra  += x * obs->posn_sigma_1 * noise_in_radians / cos( obs->dec);
      obs->dec += y * obs->posn_sigma_2 * noise_in_radians;
      if( obs->obs_mag != BLANK_MAG)
         obs->obs_mag += gaussian_random( ) * obs->mag_sigma;
      obs->jd  +=     gaussian_random( ) * obs->time_sigma;
      if( obs->note2 == 'S')
         set_obs_vect( obs);
      else
         set_up_observation( obs);
      obs++;
      }
   return( rval);
}

   /* For some time,  I displayed the covariance matrix and sigmas using the
      sprintf format specifier %10.3g.  That worked well,  except that values
      such as 1010 or 999999 were rendered as 1.01e+003 or 9.99e+005.  (999
      was left as 999,  and I've no problem with the use of scientific notation
      beyond six digits... purely a matter of personal preference;  I didn't
      want four,  five,  and six-digit numbers shown in SN.)  Also,  this
      code will remove leading zeroes in the exponent,  so that 'e+004',
      for example,  becomes 'e+4'.
         (2020 Jun 18) One can set SIGMA_
       */

static char format[10];
static int precision;

char *put_double_in_buff( char *buff, const double ival)
{
   int i;
   double low_end = .999999;

   if( isnan( ival))
      {
      strcpy( buff, "NaN");
      return( buff);
      }
   if( !precision)
      {
      precision = atoi( get_environment_ptr( "FULL_SIGMAS"));
      if( !precision)
         precision = 3;
      assert( precision > 0 && precision < 20);
      snprintf_err( format, sizeof( format), "%%%d.%dg", precision + 7, precision);
      }
   for( i = precision; i; i--)
      low_end *= 10.;

   if( fabs( ival) < low_end || fabs( ival) > low_end * 1000.)
      {
      char *tptr;

      snprintf_err( buff, 20, format, ival);
      while( (tptr = strchr( buff, 'e')) != NULL
                     &&  tptr[2] == '0')
         {           /* remove a leading zero from exponent */
         memmove( buff + 1, buff, tptr - buff + 2);
         *buff = ' ';
         }
      }
   else
      snprintf_err( buff, 20, "%*ld", precision + 7, (long)ival);
   while( *buff == ' ')
      buff++;
   return( buff);
}

   /* The following should allow display of sigmas without showing */
   /* "insignificant" digits.  For example, '0.003141" should be   */
   /* lopped down to '0.003'.  To do this,  we scan for a non-zero */
   /* digit in the string.  If that digit is a 1 or 2,  we'll add  */
   /* on one more digit;  e.g.,  "0.002131" would be lopped to     */
   /* '0.0021'.  Then we scan for more digits,  removing them.     */
   /*     Note that the decimal point requires some extra logic:   */
   /*   3141.59 -> 3000                */
   /*   27184.2 -> 27000               */
   /*   3.141   -> 3                   */
   /*   2.718   -> 2.7                 */


void remove_insignificant_digits( char *tbuff)
{
   int i;
   long total = 0, limit = 6;
   char *tptr;

   for( i = precision; i > 3; i--)
      limit *= 10;
   while( *tbuff && total < limit)
      {
      if( *tbuff >= '0' && *tbuff <= '9')
         total = total * 10 + (*tbuff - '0');
      tbuff++;
      }
   tptr = strchr( tbuff, '.');
   if( tptr)
      memset( tbuff, '0', tptr - tbuff);
   else
      tptr = tbuff;
   *tptr = '\0';
}


/* Just to explicate some of the error calculus below:

   We actually don't know sigma_a right off the bat,  since a is
not entirely defined for near-parabolic orbits.  Instead,  we have
sigma(1/a),  from which we can compute

sigma_a = sigma(1/a) * a^2

   Then,  since P_years = a^1.5,

sigma_P_years = sigma_a * 1.5 * sqrt(a)

   ...and,  since n = 360. / period in days = 360 / (days_per_year * P_years),

sigma_n = 360 * sigma_P_in_days / P_days^2
*/

double dump_monte_data_to_file( FILE *ofile, const double *sigmas,
            const double semimajor_axis, const double ecc,
            const int planet_orbiting)
{
   double uparam = 97.;   /* assume we won't get a "correct" uncertainty */
   extern const char *monte_label[];                     /* orb_func.cpp */
   static const char *units_text[MONTE_N_ENTRIES] = {
                           "days", "", "AU", "AU", "1/AU", "deg",
                           "deg", "deg", "deg", "AU", "mag" };
   const double sigma_a = sigmas[MONTE_INV_A] * semimajor_axis * semimajor_axis;
   int i;
   char tbuff[40], *tptr;
   extern int available_sigmas;

   fprintf( ofile, "Planet orbiting: %d\n", planet_orbiting);
   fprintf( ofile, "Sigmas:\n");
   for( i = 0; i < MONTE_N_ENTRIES; i++)
      {
      if( !strcmp( units_text[i], "deg"))
         {
         char zbuff[40];

         snprintf_err( zbuff, sizeof( zbuff), "%.8f", sigmas[i]);
         remove_insignificant_digits( zbuff);
         if( strlen( zbuff) == 10)   /* very low value */
            put_double_in_buff( zbuff, sigmas[i]);
         snprintf_err( tbuff, sizeof( tbuff), "%10s", zbuff);
         }
      else
         put_double_in_buff( tbuff, sigmas[i]);
      fprintf( ofile, "sigma_%-5s%s %s",
                  monte_label[i], tbuff, units_text[i]);
      if( !strcmp( units_text[i], "AU"))    /* show in km,  too */
         {
         tptr = put_double_in_buff( tbuff, sigmas[i] * AU_IN_KM);
         fprintf( ofile, " (%s km)", tptr);
         }
      fprintf( ofile, "\n");
      }
// if( semimajor_axis > sigma_a * .3)
   if( semimajor_axis > 0.)
      {
      const double GAUSS_K = .01720209895;
      const double SOLAR_GM = GAUSS_K * GAUSS_K;
      const double mass = get_planet_mass( planet_orbiting) / SOLAR_GM;
      const double per_yrs = semimajor_axis * sqrt( semimajor_axis / mass);
      const double days_per_year = 365.25;
      const double per_days = per_yrs * days_per_year;
//    const double sigma_P_in_days = days_per_year * 1.5 * sigma_a
//            * sqrt( semimajor_axis / mass);
      const double sigma_P_in_days = per_days * 1.5 * sigma_a / semimajor_axis;
      const double sigma_n = 360. * sigma_P_in_days / (per_days * per_days);
      const double runoff_coeff =
               3. * 3600. * (180. / PI) * GAUSS_K;
#ifdef ORIGINAL_FORMULATION
      const double runoff = (runoff_coeff / per_yrs) *
             (sigmas[0] * ecc + 10. * sigma_P_in_days / per_yrs);
#else
      const double runoff = runoff_coeff *
             (sigmas[0] * ecc + 15. * days_per_year * sigmas[MONTE_INV_A] * semimajor_axis) / per_yrs;
#endif
      const double uparam_const = 1.49;    /*    = ln( 648000) / 9 */

      uparam = log( runoff) / uparam_const + 1.;

      if( semimajor_axis > sigma_a * .3)
         {
         put_double_in_buff( tbuff, sigma_n);
         fprintf( ofile, "sigma_n:   %s\n", tbuff);
         put_double_in_buff( tbuff, sigma_a);
         fprintf( ofile, "sigma_a:   %s AU", tbuff);
         tptr = put_double_in_buff( tbuff, sigma_a * AU_IN_KM);
         fprintf( ofile, " (%s km)\n", tptr);
         put_double_in_buff( tbuff, sigma_P_in_days);
         fprintf( ofile, "sigma_P:   %s days\n", tbuff);
         put_double_in_buff( tbuff, sigma_P_in_days / 365.25);
         fprintf( ofile, "sigma_Py: %s years\n", tbuff);
         }
      if( uparam < 99.9 && uparam > -8.0)
         fprintf( ofile, "U=%.1f\n", uparam);
      else
         uparam = 97.;     /* 'bogus' U */
      }
   available_sigmas = MONTE_CARLO_SIGMAS_AVAILABLE;
   return( uparam);
}

/* From https://www.minorplanetcenter.net/iau/info/UValue.html ,  slightly
modified for clarity (there were operator-precedence issues) :
---------- (mostly) quoted begin ----------
The U value is calculated in the following manner. First, calculate:

      RUNOFF = (sigma_T * e + 10 * sigma_P / P) * GAUSS_K * (180. / PI) * 3600 * 3 / P
             = k * (sigma_T * e / P + 10 * sigma_P / P^2)

where sigma_T is the uncertainty in the perihelion time (in days)
      e is the eccentricity
      P is the orbital period (in years)
      sigma_P is the uncertainty in the orbital period (in days)
      GAUSS_K = Gaussian constant = 0.01720209895 radians/day
      k is three times the Gaussian constant in arcseconds
             = 3. * GAUSS_K * (180. / pi) * 3600.
             = 3. * 360 * 3600 / 365.256897...
      3600 converts to seconds of arc
      3 is a empirical factor to make the formal errors more
         closely model reality

and   RUNOFF is the in-orbit longitude runoff in seconds of
         arc per decade

---------- (mostly) quoted end ----------

   The existence of most of those constants has to do with a dog's
breakfast of units being in use.  sigma_T and sigma_P are in days,
but the orbital period P is in years.  The number of days in a year
is 2 * pi / GAUSS_K = 365.256897;  the main purpose of GAUSS_K appearing
in the formula is to handle the days/year conversion.  If the 'runoff'
were in terms of complete revolutions and we used sigma_Ty = sigma_T /
365.256897 and sigma_Py = sigma_P / 365.256897 (i.e.,  sigmas in units
of years rather than days),  we'd get

runoff_in_revolutions_per_decade =
         3 * (sigma_Ty * e / P + 10 * sigma_Py / P^2)

   ...with the factor of ten being in decades.  Note also that sigma_nR,
the uncertainty in the mean motion expressed in units of revolutions/year,
is sigma_Py / P^2,  making

runoff_in_revolutions_per_decade =
         3 * (sigma_Ty * e / P + 10 * sigma_nR)

   The sigma_Ty * e / P term is still rather mysterious to me.

   This obviously has issues for nearly-parabolic or hyperbolic orbits.
However,  we can get around that by switching from use of P to use of
the semimajor axis a = P^(2/3),  P = a ^ 1.5.  So (allowing for the
switch from years to days)

sigma_Py = sigma_a * da/dP = 1.5 * sigma_a * sqrt(a)

RUNOFF = 3 * 3600 * 360 * (sigma_Ty * e + 10 * sigma_Py / P) / P
       = 3 * 3600 * 360 * (sigma_Ty * e + 10 * 1.5 * sigma_a / a) / P

   To avoid singularities,  we're actually calculating
sigma(1/a) = sigma_a / a^2,  so (all of the following are equivalent)

RUNOFF = 3 * 3600 * 180 * (sigma_Ty * e     + 10 * 1.5 * sigma(1/a) * a) / P
RUNOFF = 3 * 3600 * 360 * (sigma_Ty * e / P + 10 * 1.5 * sigma(1/a) / sqrt( a))
RUNOFF = 3 * 3600 * 360 * (sigma_Ty * e / a + 10 * 1.5 * sigma(1/a)) / sqrt( a)

   Note that RUNOFF approaches zero as the semimajor axis approaches
infinity (i.e.,  a parabolic orbit),  and is a purely imaginary
(i.e.,  no real part) number for hyperbolic orbits.  It would be
nice to have an 'orbit quality' code that would produce useful
results for all orbits.
*/
