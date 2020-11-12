/* ephem0.cpp: low-level funcs for ephemerides & pseudo-MPECs

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

#ifdef _WIN32
   #include <direct.h>        /* for _mkdir() definition */
#else
   #include <sys/stat.h>
   #include <sys/types.h>
#endif
#ifdef __linux
   #include <linux/limits.h>
#endif
#include <sys/types.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
#include <stdbool.h>
#include "watdefs.h"
#include "afuncs.h"
#include "lunar.h"
#include "date.h"
#include "comets.h"
#include "mpc_obs.h"
#include "mpc_func.h"
#include "vislimit.h"
#include "brentmin.h"
#include "expcalc.h"
#include "rgb_defs.h"

#define J2000 2451545.0
#define JD_TO_YEAR(jd)  (2000. + ((jd)-J2000) / 365.25)
#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_MAJOR_AXIS_IN_AU (EARTH_MAJOR_AXIS / AU_IN_METERS)
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define LOG_10 2.3025850929940456840179914546843642076011014886287729760333279009675726
#define LIGHT_YEAR_IN_KM    (365.25 * seconds_per_day * SPEED_OF_LIGHT)

int save_ephemeris_file( const char *filename);       /* ephem0.cpp */
double centralize_ang( double ang);             /* elem_out.cpp */
double vector_to_polar( double *lon, double *lat, const double *vector);
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
int parallax_to_lat_alt( const double rho_cos_phi, const double rho_sin_phi,
       double *lat, double *ht_in_meters, const int planet_idx); /* ephem0.c */
double calc_obs_magnitude( const double obj_sun,
          const double obj_earth, const double earth_sun, double *phase_ang);
int lat_alt_to_parallax( const double lat, const double ht_in_meters,
             double *rho_cos_phi, double *rho_sin_phi, const int planet_idx);
int write_residuals_to_file( const char *filename, const char *ast_filename,
          const int n_obs, const OBSERVE FAR *obs_data, const int format);
void put_observer_data_in_text( const char FAR *mpc_code, char *buff);
int make_pseudo_mpec( const char *mpec_filename, const char *obj_name);
                                              /* ephem0.cpp */
int earth_lunar_posn( const double jd, double FAR *earth_loc,
                                       double FAR *lunar_loc);
void remove_trailing_cr_lf( char *buff);      /* ephem0.cpp */
void create_obs_file( const OBSERVE FAR *obs, int n_obs, const int append);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
void format_dist_in_buff( char *buff, const double dist_in_au); /* ephem0.c */
int debug_printf( const char *format, ...)                 /* runge.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
int calc_derivatives( const double jd, const double *ival, double *oval,
                           const int reference_planet);     /* runge.cpp */
char *iso_time( char *buff, const double jd, const int precision);   /* elem_out.c */
double mag_band_shift( const char mag_band);                /* elem_out.c */
char *get_file_name( char *filename, const char *template_file_name);
double current_jd( void);                       /* elem_out.cpp */
double diameter_from_abs_mag( const double abs_mag,      /* ephem0.cpp */
                                     const double optical_albedo);
double shadow_check( const double *planet_loc,           /* ephem0.cpp */
                            const double *obs_posn,
                            const double planet_radius_in_au);
int get_object_name( char *obuff, const char *packed_desig);   /* mpc_obs.c */
int get_residual_data( const OBSERVE *obs, double *xresid, double *yresid);
int snprintf_append( char *string, const size_t max_len,      /* ephem0.cpp */
                                   const char *format, ...)
#ifdef __GNUC__
         __attribute__ (( format( printf, 3, 4)))
#endif
;
int setup_planet_elem( ELEMENTS *elem, const int planet_idx,
                                          const double t_cen);   /* moid4.c */
char *mpc_station_name( char *station_data);       /* mpc_obs.cpp */
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
int remove_rgb_code( char *buff);                              /* ephem0.cpp */
static void output_signed_angle_to_buff( char *obuff, const double angle,
                               const int precision);         /* ephem0.cpp */
static void output_angle_to_buff( char *obuff, const double angle,
                               const int precision);         /* ephem0.cpp */
static void put_residual_into_text( char *text, const double resid,
                                 const int resid_format);    /* ephem0.cpp */
FILE *open_json_file( char *filename, const char *env_ptr, const char *default_name,
                  const char *packed_desig, const char *permits); /* ephem0.cpp */
size_t strlcat( char *dst, const char *src, size_t dsize);     /* miscell.cpp */
size_t strlcpy( char *dst, const char *src, size_t dsize);     /* miscell.cpp */
size_t strlcpy_err( char *dst, const char *src, size_t dsize); /* miscell.c */
size_t strlcat_err( char *dst, const char *src, size_t dsize); /* miscell.c */
void shellsort_r( void *base, const size_t n_elements, const size_t elem_size,
         int (*compare)(const void *, const void *, void *), void *context);

const char *observe_filename = "observe.txt";
const char *residual_filename = "residual.txt";
const char *ephemeris_filename = "ephemeri.txt";
bool is_default_ephem = true;
const char *elements_filename = "elements.txt";

static expcalc_config_t exposure_config;

/* Returns parallax constants (rho_cos_phi, rho_sin_phi) in AU. */

int lat_alt_to_parallax( const double lat, const double ht_in_meters,
            double *rho_cos_phi, double *rho_sin_phi, const int planet_idx)
{
   const double axis_ratio = planet_axis_ratio( planet_idx);
   const double major_axis_in_meters = planet_radius_in_meters( planet_idx);
   const double u = atan( sin( lat) * axis_ratio / cos( lat));

   *rho_sin_phi = axis_ratio * sin( u) +
                            (ht_in_meters / major_axis_in_meters) * sin( lat);
   *rho_cos_phi = cos( u) + (ht_in_meters / major_axis_in_meters) * cos( lat);
   *rho_sin_phi *= major_axis_in_meters / AU_IN_METERS;
   *rho_cos_phi *= major_axis_in_meters / AU_IN_METERS;
   return( 0);
}

/* Parallax constants (rho_cos_phi, rho_sin_phi) should be in units of the
planet's equatorial radius.  */

int parallax_to_lat_alt( const double rho_cos_phi, const double rho_sin_phi,
               double *lat, double *ht_in_meters, const int planet_idx)
{
   double talt = 0.;
   const double tlat = point_to_ellipse( 1., planet_axis_ratio( planet_idx),
                  rho_cos_phi, rho_sin_phi, &talt);

   if( lat)
      *lat = tlat;
   if( ht_in_meters)
      *ht_in_meters = talt * planet_radius_in_meters( planet_idx);
   return( 0);
}

#include <stdarg.h>
#if defined(_MSC_VER) && _MSC_VER < 1900
                      /* For older MSVCs,  we have to supply our own  */
                      /* snprintf().  See snprintf.cpp for details.  */
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

/* I find myself frequently snprintf()-ing at the end of a string,  with    */
/* something like snprintf( str + strlen( str), sizeof( str) - strlen(str), */
/* ...).  This should be a little more convenient.                          */

int snprintf_append( char *string, const size_t max_len,      /* ephem0.cpp */
                                   const char *format, ...)
{
   va_list argptr;
   int rval;
   const size_t ilen = strlen( string);

   assert( ilen <= max_len);
   va_start( argptr, format);
#if _MSC_VER <= 1100
   rval = vsprintf( string + ilen, format, argptr);
#else
   rval = vsnprintf( string + ilen, max_len - ilen, format, argptr);
#endif
   string[max_len - 1] = '\0';
   va_end( argptr);
   return( rval);
}

/* format_dist_in_buff() formats the input distance (in AU) into a
seven-byte buffer.  It does this by choosing suitable units: kilometers
if the distance is less than a million km,  AU out to 10000 AU,  then
light-years.  In actual use,  light-years indicates some sort of error,
but I wanted the program to handle that problem without just crashing.

   28 Dec 2009:  To help puzzle out a strange bug,  I extended things so
that distances beyond 10000 light-years could be shown in kilo-light-years
(KLY),  mega,  giga,  etc. units,  making up prefixes after 'yocto'.  Thus,
one "ALY" = 10^78 light-years.  Since the universe is 13.7 billion years
old,  and nothing can be seen beyond that point,  such distances are
physically unreasonable by a factor of roughly 10^68.  But for debugging
purposes,  display of these distances was useful.  The bug is now fixed,
so with any luck,  nobody should see such distances again.

   2012 Feb 13:  there's some interest in determining orbits of _very_
close objects,  such as tennis balls.  To address this,  short distances
are now shown in millimeters,  centimeters,  meters,  or .1 km,
as appropriate.

NOTE:  It used to be that I followed MPC practice in showing all distances
between .01 and ten AU in the form d.dddd , but now,  that's only true for
one to ten AU.  For distances less than 1 AU,  I'm using .ddddd,  thereby
getting an extra digit displayed.   */

static void show_dist_in_au( char *buff, const double dist_in_au)
{
   const char *fmt;

   if( dist_in_au > 999.999)
      fmt = "%7.1f";             /* " 1234.5" */
   else if( dist_in_au > 99.999)
      fmt = "%7.2f";             /* " 123.45" */
   else if( dist_in_au > 9.999)
      fmt = "%7.3f";             /* " 12.345" */
   else if( dist_in_au > .99)     /* used to be .01 */
      fmt = "%7.4f";             /* " 1.2345" */
   else
      fmt = "%7.5f";             /* " .12345" */
   snprintf( buff, 8, fmt, dist_in_au);
   *buff = ' ';   /* remove leading zero for small amounts */
}

static const char *si_prefixes = "kMGTPEZYXWVUSRQONLJIHFDCBA";
static bool use_au_only = false;

/* Given a non-negative value,  this gives a four-character output
such as '3.14', '.314', '3145', '314k', '3.1M',  etc.  */

static void show_packed_with_si_prefixes( char *buff, double ival)
{
   *buff = '\0';
   if( ival > 999.e+21)
      strcpy( buff, "!!!!");
   else if( ival > 9999.)
      {
      unsigned count = 0;

      *buff = '\0';
      while( !*buff)
         {
         ival /= 1000.;
         if( ival < 9.9)
            snprintf( buff, 5, "%3.1f%c", ival, si_prefixes[count]);
         else if( ival < 999.)
            snprintf( buff, 5, "%3u%c", (unsigned)ival, si_prefixes[count]);
         count++;
         }
      }
   else if( ival > 99.9)
      snprintf( buff, 5, "%4u", (unsigned)( ival + .5));
   else if( ival > 9.99)
      snprintf( buff, 5, "%4.1f", ival);
   else if( ival > .99)
      snprintf( buff, 5, "%4.2f", ival);
   else
      {
      char tbuff[7];

      snprintf( tbuff, sizeof( tbuff), "%5.2f", ival);
      strcpy( buff, tbuff + 1);  /* store value without leading 0 */
      }
}

void format_dist_in_buff( char *buff, const double dist_in_au)
{
   if( dist_in_au < 0.)
      strcpy( buff, " <NEG!>");
   else if( use_au_only == true)
      show_dist_in_au( buff, dist_in_au);
   else
      {
      const double dist_in_km = dist_in_au * AU_IN_KM;
      const char *fmt;

                  /* for objects within a million km (about 2.5 times   */
                  /* the distance to the moon),  switch to km/m/cm/mm:  */
      if( dist_in_km < .0099)                 /* 0 to 9900 millimeters: */
         snprintf( buff, 8, "%5.0fmm", dist_in_km * 1e+6);    /* " NNNNmm" */
      else if( dist_in_km < .099)             /* 990 to 9900 centimeters: */
         snprintf( buff, 8, "%5.0fcm", dist_in_km * 1e+5);    /* " NNNNcm" */
      else if( dist_in_km < 99.)          /* 99 to 99000 meters: */
         snprintf( buff, 8, "%6.0fm", dist_in_km * 1e+3);     /* " NNNNNm" */
      else if( dist_in_km < 999.)         /* 99.0 to 999.9 kilometers: */
         snprintf( buff, 8, "%6.1fk", dist_in_km);            /* " NNN.Nk" */
      else if( dist_in_km < 999999.)      /* 999.9 to 999999 km: */
         snprintf( buff, 8, "%7.0f", dist_in_km);
      else if( dist_in_au > 9999.999)
         {
         double dist_in_light_years =
                               dist_in_au * AU_IN_KM / LIGHT_YEAR_IN_KM;

         if( dist_in_light_years > 9999.9)
            {
            int idx = 0;

            dist_in_light_years /= 1000;
            for( idx = 0; si_prefixes[idx] && dist_in_light_years > 999.; idx++)
               dist_in_light_years /= 1000;
            if( !si_prefixes[idx])   /* can't represent this even in our */
               strcpy( buff, " <HUGE>");     /* largest made-up units */
            else
               {
               if( dist_in_light_years < 9.9)
                  snprintf( buff, 8, "%4.1fxLY", dist_in_light_years);
               else
                  snprintf( buff, 8, "%4.0fxLY", dist_in_light_years);
               buff[4] = si_prefixes[idx];
               }
            }
         else
            {
            if( dist_in_light_years > 99.999)  /* " 1234LY" */
               fmt = "%5.0fLY";
            else if( dist_in_light_years > 9.999)   /* " 12.3LY" */
               fmt = "%5.1fLY";
            else if( dist_in_light_years > .999)
               fmt = "%5.2fLY";           /* " 1.23LY" */
            else
               fmt = "%5.3fLY";           /* " .123LY" */
            snprintf( buff, 8, fmt, dist_in_light_years);
            }
         }
      else
         show_dist_in_au( buff, dist_in_au);
      *buff = ' ';                /* remove leading zero for small amts */
      }
}


      /* Input velocity is in km/s.  If it's greater than about     */
      /* three times the speed of light,  we show it in units of c. */
      /* Of course,  this will never happen with real objects;  I   */
      /* added it for debugging purposes.                           */
static void format_velocity_in_buff( char *buff, double vel)
{
   const char *format;

   if( vel < 9.999 && vel > -9.999)
      format = "%7.3f";
   else if( vel < 99.999 && vel > -99.999)
      format = "%7.2f";
   else if( vel < 999.9 && vel > -999.9)
      format = "%7.1f";
   else if( vel < 999999. && vel > -999999.)
      format = "%7.0f";
   else         /* show velocity in terms of speed of light. */
      {
      vel /= SPEED_OF_LIGHT;
      if( vel < 99.999 && vel > -99.999)
         format = "%6.1fc";
      else if( vel < 999999. && vel > -999999.)
         format = "%6.0fc";
      else        /* we give up;  it's too fast */
         format =  " !!!!!!";
      }
   snprintf( buff, 8, format, vel);
}

/* Rob Matson asked about having the program produce ECF (Earth-Centered
Fixed) coordinates,  in geometric lat/lon/altitude form.  This 'ground
track' option in ephemerides has since been expanded to allow geodetic
output,  with the possibility of creating such ephemerides from other
planets as well. */

double find_lat_lon_alt( const double ut, const double *ivect,
                  const int planet_no, double *lat_lon, const bool geometric);

/* 'get_step_size' parses input text to get a step size in days,  so that */
/* '4h' becomes .16667 days,  '30m' becomes 1/48 day,  and '10s' becomes  */
/* 10/(24*60*60) days.  The units (days, hours, minutes,  or seconds) are */
/* returned in 'step_units' if the input pointer is non-NULL.  The number */
/* of digits necessary is returned in 'step_digits' if that pointer is    */
/* non-NULL.  Both are used to ensure correct time format in ephemeris    */
/* output;  that is,  if the step size is (say) .05d,  the output times   */
/* ought to be in a format such as '2009 Mar 8.34', two places in days.   */

double get_step_size( const char *stepsize, char *step_units, int *step_digits)
{
   double step = 0.;
   char units = 'd';

   if( !strncmp( stepsize, "Obs", 3) || *stepsize == 't')     /* dummy value */
      step = 1e-6;
   if( *stepsize == 'a')
      {
      step = 1e-6;
      if( step_digits)
         *step_digits = 0;
      }
   if( sscanf( stepsize, "%lf%c", &step, &units) >= 1)
      if( step)
         {
         if( step_digits)
            {
            const char *tptr = strchr( stepsize, '.');

            *step_digits = 0;
            if( tptr)
               {
               tptr++;
               while( isdigit( *tptr++))
                  (*step_digits)++;
               }
#ifdef OBSOLETE_METHOD
            double tval = fabs( step);

            for( *step_digits = 0; tval < .999; (*step_digits)++)
               tval *= 10.;
#endif
            }
         units = tolower( units);
         if( step_units)
            *step_units = units;
         switch( units)
            {
            case 'd':
               break;
            case 'h':
               step /= hours_per_day;
               break;
            case 'm':
               step /= minutes_per_day;
               break;
            case 's':
               step /= seconds_per_day;
               break;
            case 'w':
               step *= 7.;
               break;
            case 'y':
               step *= 365.25;
               break;
            }
         }
   return( step);
}

static double centralize_ang_around_zero( double ang)
{
   ang = fmod( ang, PI + PI);
   if( ang > PI)
      ang -= PI + PI;
   else if( ang <= -PI)
      ang += PI + PI;
   return( ang);
}

typedef struct
{
   double ra, dec, jd, r;
   double sun_obj, sun_earth;
} obj_location_t;

static void setup_obj_loc( obj_location_t *p, double *orbit,
          const size_t n_orbits, const double epoch_jd, const char *mpc_code)
{
   size_t i, j;
   double rho_sin_phi = 0., rho_cos_phi = 0., lon = 0.;
   double obs_posn[3];
   int planet_no = 3;

   assert( p->jd > 2e+6);
   assert( p->jd < 3e+6);
   if( mpc_code)
      planet_no = get_observer_data( mpc_code, NULL,
                          &lon, &rho_cos_phi, &rho_sin_phi);
   compute_observer_loc( p->jd, planet_no, rho_cos_phi, rho_sin_phi,
                                    lon, obs_posn);
   for( i = 0; i < n_orbits; i++)
      {
      double topo[3];

      if( !mpc_code)
         p[i].r = 0.;
      assert( p[i].r >= 0.);
      assert( p[i].r < 100.);
      integrate_orbit( orbit, epoch_jd, p->jd - p[i].r / AU_PER_DAY);
      for( j = 0; j < 3; j++)
         topo[j] = orbit[j] - obs_posn[j];
      ecliptic_to_equatorial( topo);
      p[i].ra = atan2( topo[1], topo[0]);
      p[i].r = vector3_length( topo);
      p[i].dec = asin( topo[2] / p[i].r);
      p[i].sun_earth = vector3_length( obs_posn);
      p[i].sun_obj = vector3_length( orbit);
      orbit += 6;
      }
}

#pragma pack( 1)

typedef struct
{
   double ra, dec, jd;
   double height, width, tilt;
   uint32_t file_offset;
   char obscode[4];
   char file_number;
} field_location_t;

typedef struct
{
   double height, width;
   double min_jd, max_jd;
   char obscode[4], file_number;
} field_group_t;

#pragma pack( )

/* In searching for fields,  we may be interested in fields _within_ a certain
range of dates,  or those _outside_ that range.  We signal the former with
max_jd > min_jd.  We signal the latter by reversing them (min_jd > max_jd.) */

static inline bool jd_is_in_range( const double jd, const double min_jd,
                                                      const double max_jd)
{
   if( max_jd > min_jd)    /* looking for fields _within_ this range */
      return( jd > min_jd && jd < max_jd);
   else                    /* looking for fields _outside_ this range */
      return( jd < min_jd || jd > max_jd);
}

static int put_ephemeris_posn_angle_sigma( char *obuff, const double dist,
              const double posn_ang, const bool computer_friendly)
{
   int integer_posn_ang =
               (int)( floor( -posn_ang * 180. / PI + .5)) % 180;
   const double dist_in_arcsec = dist * 3600. * 180. / PI;
   char resid_buff[9];

   if( integer_posn_ang < 0)
      integer_posn_ang += 180;
   if( computer_friendly)
      snprintf( resid_buff, sizeof( resid_buff), "  %6u",
                           (unsigned)dist_in_arcsec);
   else
      {
      put_residual_into_text( resid_buff, dist_in_arcsec,
                                    RESIDUAL_FORMAT_OVERPRECISE);
      resid_buff[5] = '\0';
      }
   snprintf( obuff, 13, "%s %3d", resid_buff + 1, integer_posn_ang);
   return( integer_posn_ang);
}

/* Old MSVCs and OpenWATCOM lack erf() and many other math functions: */

#if defined( _MSC_VER) && (_MSC_VER < 1800) || defined( __WATCOMC__)
double erf( double x);     /* orb_fun2.cpp */
#endif

#define SWAP( A, B, TEMP)   { TEMP = A;  A = B;  B = TEMP; }

/* At present,  this assumes a 'nominal' position in p[0] and a one-sigma
variant in p[1].  We figure out the range,  in sigmas,  of that line within
the rectangle defined by 'field' (keeping in mind that the line may miss
the rectangle entirely,  and usually does).  But if the line of variation
_does_ go through the triangle,  running from sigma1 to sigma2,  then we
evaluate the probability that the object is on that field as the area
under the normal curve between sigma1 and sigma2, and return that.

   For cases where we don't have variant orbits (n_objs == 1),  we just
set up a really low difference and are basically just checking to see if
the nominal position is on the field. (Which,  with one orbit,  is all
we can do anyway.)

   For objects with low positional uncertainty,  that probability will
usually be essentially 1 or 0.  If the uncertainty is close to the field
size,  you'll start to see "maybe it's on the field and maybe it isn't"
cases.

   TO BE DONE:  figure out how to extend this to non-covariance (SR) cases.
I've ideas for this... but first,  the (more common) covariance case.  */

static double precovery_in_field( const field_location_t *field,
         const obj_location_t *p, const unsigned n_objs,
         const double margin)
{
   double sigma1 = -10., sigma2 = 10.;
   const double d_ra = (n_objs > 1 ? p[1].ra - p[0].ra : 1e-8);
   const double d_dec = (n_objs > 1 ? p[1].dec - p[0].dec : 1e-8);
   double sig1, sig2, temp, size;
   const double sqrt_2 =
      1.414213562373095048801688724209698078569671875376948073176679737990732;

   size = (field->width / 2. + margin) / cos( field->dec);
   sig1 = centralize_ang_around_zero( field->ra + size - p[0].ra) / d_ra;
   sig2 = sig1 - 2 * size / d_ra;
   if( sig1 > sig2)
      SWAP( sig1, sig2, temp);
   if( sigma1 < sig1)
      sigma1 = sig1;
   if( sigma2 > sig2)
      sigma2 = sig2;

   size = field->height / 2. + margin;
   sig1 = (field->dec + size - p[0].dec) / d_dec;
   sig2 = (field->dec - size - p[0].dec) / d_dec;
   if( sig1 > sig2)
      SWAP( sig1, sig2, temp);
   if( sigma1 < sig1)
      sigma1 = sig1;
   if( sigma2 > sig2)
      sigma2 = sig2;

   if( sigma1 >= sigma2)
      return( 0.);
   else
      return( (erf( sigma2 / sqrt_2) - erf( sigma1 / sqrt_2)) * .5);
}

static void show_precovery_extent( char *obuff, const obj_location_t *objs,
                               const int n_objs)
{
   *obuff = '\0';
   if( n_objs == 2)
      {
      double dist, posn_ang;

      calc_dist_and_posn_ang( &objs[0].ra, &objs[1].ra,
                              &dist, &posn_ang);
      put_ephemeris_posn_angle_sigma( obuff, dist, posn_ang, true);
      }
}

/* See 'precover.txt' for information on what's going on here. */

#define FIELD_BUFF_N 1024
#define COMPRESSED_FIELD_SIZE 18


static void extract_field( field_location_t *field, const char *buff,
               const field_group_t *groups)
{
   int32_t array[4];

   assert( buff[16] >= 0);
   assert( strlen( groups->obscode) == 3);
   groups += buff[16];
   assert( strlen( groups->obscode) == 3);
   memcpy( array, buff, 4 * sizeof( int32_t));
   field->ra = (double)array[0] * 2. * PI / 2e+9;
   field->dec = (double)array[1]      * PI / 2e+9;
   field->jd = groups->min_jd + (groups->max_jd - groups->min_jd)
                        * (double)array[2] / 2e+9;
   field->file_offset = array[3];
   field->tilt = (double)buff[17] * PI / 256.;
   field->height = groups->height;
   field->width  = groups->width;
   field->file_number = groups->file_number;
   strcpy( field->obscode, groups->obscode);
}


static int find_precovery_plates( OBSERVE *obs, const int n_obs,
                           const char *idx_filename,
                           FILE *ofile, const double *orbit,
                           const int n_orbits, double epoch_jd,
                           const double min_jd, const double max_jd,
                           const double limiting_mag)
{
   FILE *ifile, *original_file = NULL;
   int current_file_number = -1;
   double *orbi, stepsize = 1.;
   obj_location_t *p1, *p2, *p3;
   int n_fields_read, n;
   const double abs_mag = calc_absolute_magnitude( obs, n_obs);
   char *buff;
        /* Slightly easier to work with 'bit set means included' : */
   const int inclusion = atoi( get_environment_ptr( "FIELD_INCLUSION")) ^ 3;
   const bool show_base_60 = (*get_environment_ptr( "FIELD_DEBUG") != '\0');
   int n_groups;
   field_group_t *groups;

   if( !ofile)
      return( -1);
   ifile = fopen_ext( idx_filename, "crb");
   if( !ifile)
      {
      debug_printf( "Couldn't open %s\n", idx_filename);
      return( -2);
      }
   p1 = (obj_location_t *)calloc( 3 * n_orbits, sizeof( obj_location_t));
   p2 = p1 + n_orbits;
   p3 = p2 + n_orbits;
   buff = (char *)calloc( FIELD_BUFF_N, COMPRESSED_FIELD_SIZE);
   assert( buff);
   if( !fgets( buff, 100, ifile))
      return( -4);
   n_groups = atoi( buff);
   assert( n_groups);
   groups = (field_group_t *)calloc( n_groups, sizeof( field_group_t));
   assert( groups);
   if( fread( groups, sizeof( field_group_t), n_groups, ifile) != (size_t)n_groups)
      return( -3);
   orbi = (double *)malloc( 12 * n_orbits * sizeof( double));
   memcpy( orbi, orbit,     6 * n_orbits * sizeof( double));
   while( (n_fields_read = (int)fread( buff, COMPRESSED_FIELD_SIZE, FIELD_BUFF_N, ifile)) > 0)
      for( n = 0; n < n_fields_read; n++)
         {
         field_location_t field;

         extract_field( &field, buff + n * COMPRESSED_FIELD_SIZE, groups);
         if( jd_is_in_range( field.jd, min_jd, max_jd))
            {
            double fraction, mag;
            double margin = .1;
            const double jdt = field.jd + td_minus_utc( field.jd) / seconds_per_day;
            bool possibly_within_field = false;
            double prob;
            int i;

            while( jdt < p1->jd || jdt > p2->jd)
               {
               const double new_p2_jd = ceil( (jdt - .5) / stepsize) * stepsize + .5;
               const double scale_factor = 2.;

               if( new_p2_jd == p1->jd)
                  memcpy( p2, p1, n_orbits * sizeof( obj_location_t));
               else if( new_p2_jd != p2->jd)
                  {
                  p2->jd = new_p2_jd;
                  setup_obj_loc( p2, orbi, n_orbits, epoch_jd, NULL);
                  epoch_jd = p2->jd;
                  }
               while( stepsize > p2->r * scale_factor)
                  stepsize /= 2.;
               while( stepsize < p2->r * scale_factor)
                  stepsize *= 2.;
               p1->jd = new_p2_jd - stepsize;
               setup_obj_loc( p1, orbi, n_orbits, epoch_jd, NULL);
               epoch_jd = p1->jd;
               }
            fraction = (jdt - p1->jd) / stepsize;
            for( i = 0; i < n_orbits; i++)        /* compute approx RA/decs */
               {
               const double delta_ra = p2[i].ra - p1[i].ra;

               p3[i].ra = p1[i].ra + fraction * centralize_ang_around_zero( delta_ra);
               p3[i].dec = p1[i].dec + fraction * (p2[i].dec - p1[i].dec);
               }
            margin += EARTH_MAJOR_AXIS_IN_AU / p1->r;
            mag = abs_mag + calc_obs_magnitude(
                                 p2->sun_obj, p2->r, p2->sun_earth, NULL);
            if( mag < limiting_mag && precovery_in_field( &field, p3, n_orbits, margin) > .01)
               {                          /* approx posn is on plate;  compute */
               double *temp_orbit = orbi + 6 * n_orbits;

               memcpy( temp_orbit, orbi, 6 * n_orbits * sizeof( double));
               memcpy( p3, p2, n_orbits * sizeof( obj_location_t));
               p3->jd = jdt;
               setup_obj_loc( p3, temp_orbit, n_orbits, epoch_jd, field.obscode);
               possibly_within_field = true;
               }
            if( mag < limiting_mag && possibly_within_field
                  && (prob = precovery_in_field( &field, p3, n_orbits, 0.)) > .1)
               {
               char time_buff[40], buff[200];
               int i;
               bool matches_an_observation = false;
               bool show_it = true;
               double obj_ra = p3->ra, obj_dec = p3->dec;

               full_ctime( time_buff, field.jd, FULL_CTIME_YMD
                           | FULL_CTIME_LEADING_ZEROES
                           | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_TENTHS_SEC);
               obj_ra = centralize_ang( obj_ra);
               for( i = 0; i < n_obs; i++)
                  if( fabs( obs[i].jd - jdt) < 1.e-3
                           && !strcmp( field.obscode, obs[i].mpc_code))
                     matches_an_observation = true;
               if( matches_an_observation)
                  show_it = ((inclusion & 2) != 0);
               else
                  show_it = ((inclusion & 1) != 0);
               if( show_it)
                  {
                  obj_ra *= 180. / PI;
                  obj_dec *= 180. / PI;
                  if( !show_base_60)
                     snprintf( buff, sizeof( buff), "%8.4f %8.4f",
                                       obj_ra, obj_dec);
                  else
                     {
                     output_angle_to_buff( buff, obj_ra / 15., 3);
                     buff[12] = ' ';
                     output_signed_angle_to_buff( buff + 13, obj_dec, 2);
                     }
                  fprintf( ofile, "%c %s %4.1f %s %s",
                           (matches_an_observation ? '*' : ' '), buff, mag,
                           time_buff, field.obscode);
                  if( current_file_number != field.file_number)
                     {
                     char filename[20];

                     current_file_number = field.file_number;
                     snprintf( filename, sizeof( filename), "css_%d.csv",
                                    current_file_number);
                     if( original_file)
                        fclose( original_file);
                     original_file = fopen_ext( filename, "crb");
                     if( !original_file)
                        fprintf( ofile, "'%s' not opened\n", filename);
                     }
                  show_precovery_extent( buff, p1, n_orbits);
                  snprintf_append( buff, sizeof( buff), " %.3f ", prob);
                  fprintf( ofile, "%s", buff);

                  if( original_file)
                     {
                     fseek( original_file, field.file_offset, SEEK_SET);
                     if( fgets_trimmed( buff, sizeof( buff), original_file))
                        {
                        for( i = 0; buff[i]; i++)
                           if( buff[i] == ',')
                              buff[i] = ' ';
                        fprintf( ofile, " %s", buff);
                        }
                     else
                        fprintf( ofile, "File %d: seeked to %ld and failed",
                                 (int)field.file_number, (long)field.file_offset);
                     }
                  fprintf( ofile, "\n");
                  }
               }
            }
         }
   fclose( ifile);
   free( buff);
   free( groups);
   if( original_file)
      fclose( original_file);
   free( orbi);
   free( p1);
   return( 0);
}

/* In the following,  I'm assuming an object with H=0 and albedo=100% to
have a diameter of 1300 km.  Return value is in meters. */

double diameter_from_abs_mag( const double abs_mag,
                                     const double optical_albedo)
{
   return( 1300. * 1000. * pow( .1, abs_mag / 5.) / sqrt( optical_albedo));
}

#define RADAR_DATA struct radar_data

RADAR_DATA
   {
   double power_in_watts;
   double system_temp_deg_k;
   double gain, radar_constant;
   double altitude_limit;
   };

static inline int get_radar_data( const char *mpc_code, RADAR_DATA *rdata)
{
   char tbuff[20];
   const char *tptr;
   int rval = -1;

   memset( rdata, 0, sizeof( RADAR_DATA));
   snprintf( tbuff, sizeof( tbuff), "RADAR_%.3s", mpc_code);
   tptr = get_environment_ptr( tbuff);
   if( tptr)
      if( sscanf( tptr, "%lf,%lf,%lf,%lf,%lf",
            &rdata->power_in_watts, &rdata->system_temp_deg_k,
            &rdata->gain, &rdata->altitude_limit,
            &rdata->radar_constant) == 5)
         rval = 0;
   rdata->altitude_limit *= PI / 180.;
   return( rval);
}

/* In computing radar SNR,  we need to either measure the rotation period
(if we're lucky,  and somebody does a light curve) or take an educated guess
at it,  as a function of absolute magnitude.  In general,  smaller rocks
spin faster than bigger rocks.

   I asked about this,  and Mike Nolan replied : "I tend to use 4h for
larger than H=21, 1 for H=21-25, and .25 for H > 25. We don't really know
the distribution for the small stuff, as the biases are pretty big." To
which Lance Benner replied : "We've been adopting a somewhat more
conservative value of P = 2.1 h for objects with  H < 22.   This is based
on the observed distribution of NEA rotation periods (although there are a
couple of exceptions).  For smaller objects, a similar cutoff doesn't seem
to exist, but rapid rotators are relatively common, so we adopt P = 0.5h."

   In the following,  I went with a three-hour period for H<=21 and a .3
hour period for H>=25,  linearly interpolating for 21 < H < 25.

   Fortunately,  SNR/day in the formula below depends on the square root of
the period.  If we're out by a factor of two on our guess,  it'll mess up the
SNR/day by a factor of the square root of two.  I'd guess that errors due to
under- or over-estimating the albedo or H probably cause much more trouble.
*/

static double guessed_rotation_period_in_hours( const double abs_mag)
{
   double rval;
   const double big_limit = 21., big_period = 3.;  /* H=21 rocks spin in 3h */
   const double small_limit = 25., small_period = 0.3;
                           /* little H=25 rocks spin in a mere 18 minutes */

   if( abs_mag < big_limit)
      rval = big_period;
   else if( abs_mag < small_limit)
      rval = small_period +
               (big_period - small_period) * (small_limit - abs_mag) / (small_limit - big_limit);
   else        /* small,  fast-spinning rocks  */
      rval = small_period;
   return( rval);
}

/*   In a private communication,  Lance Benner wrote:

    "...The signal-to-noise ratio is proportional to:

[(Ptx)(radar albedo)(D^3/2)(P^1/2)(t^1/2)]  <-- numerator
------------------------------------------
[(r^4)*(Tsys)]                              <-- denominator

where Ptx is the transmitter power, the radar albedo is the radar cross section
divided by the projected area of the object, D is the diameter, P is the
rotation period, t is the integration time, r is the distance from the
telescope, and Tsys is the system temperature."

   Implicitly,  this is also all multiplied by a gain factor,  and a
constant known in this code as 'radar_constant'.

   He also mentioned that Arecibo operates at about 500 kW,  Tsys = 25 K,
gain = 10 K/Jy;  Goldstone at 450 kW,  Tsys = 15, gain = 0.94.  These
values are enshrined in the RADAR_251 and RADAR_253 lines in
'environ.dat', along with the 'radar_constant' values.  I got those
essentially by guessing, then scaling them until they matched values in
radar planning documents.  (Unfortunately, Arecibo is now at 350 kW and
gain=7 K/Jy due to Maria damage.)

   See also Lance's comments at

https://groups.io/g/mpml/message/28412

   and an excellent page on the subject at

http://www.naic.edu/~pradar/detect.php

   Lance also noted that there's also a more subtle dependence due to the subradar
latitude (you get a better echo if the pole is pointing straight at the radar).

   Jean-Luc Margot's page:

http://www2.ess.ucla.edu/~jlm/research/NEAs/intro.html

   has a graphic showing this dependence.  A 20-meter object with a four-hour
period and albedo=.1,  at a distance of .1 AU,  has SNR/day = 1.

   For Goldstone,  the SNR is about 1/20 of this (from JLM's page;  I think
this assumes Arecibo going at 900 kW.)  The same page shows a dropoff in SNR
due to declination:  the above formula holds at dec=+18,  but SNR drops to
zero at dec=-1 or dec=37.  In reality,  it should be a drop in _altitude_ :
i.e.,  Arecibo sees "perfectly" at the zenith,  and degrades as you drop from
there.  Mike Nolan tells me (private e-mail communication) that the width is
about 19.5 degrees,  so that Arecibo can "see"  down to about altitude = 70.5
degrees.   The drop factor isn't given (and is probably empirically derived),
but it looks a lot like

SNR( alt) = SNR( zenith) * ((90 - alt)/19.5) ^ .25,     70.5 < alt < 90

   Or something like that,  perhaps with a different exponent:  a function
equal to one at alt=90 degrees and zero at alt=70.5 degrees,  with roughly
the right behavior in between.

   See further comments in 'environ.def'.

   Other quantities of interest for radar:  the 'number of runs',  i.e.,
for the time the object is above the horizon,  you'll start by transmitting
power for the entire round-trip time to the object.  Then you shut off and
listen for the returning data for an equal time.  Then you switch back to
transmitting,  and back to receiving...

n_runs = (time the object is above horizon) / (twice the round-trip time)

   In bistatic mode,  one setup transmits continuously while the other
receives continuously,  so that factor of two should go away (it's as
if you're getting in twice as many runs.)  Radar ephemerides also show
"SNR/run" and "SNR/day",  where that latter quantity is SNR/run times
the square root of the number of runs (i.e.,  it's more like "SNR per
time we spend observing that day,  which will probably not be 24 hours.")

   Flaviane C. F. Venditti tells me that Arecibo "can detect objects as
close as ~0.008 AU."  Round-trip light time for an object 1 AU away is
998.01 seconds,  leading to the convenient relationship that a target
at 0.001 AU has an almost exactly one-second delay between pulse and
return.  At 0.008 AU in monostatic mode,  the transmitter has to be
toggled every eight seconds.  "The times I observed with short round
trip times (<10 seconds), the transmitter tripped frequently, so we try
not to push it too much."  I gather they _could_ go a bit closer,  but
it'd require a really,  really good reason.  */

double optical_albedo = .1;

static double radar_snr_per_day( const RADAR_DATA *rdata,
                     const double abs_mag,
                     const double radar_albedo,
                     const double dist)
{
   const double rotation_period = guessed_rotation_period_in_hours( abs_mag);
   const double diameter_in_meters =
                 diameter_from_abs_mag( abs_mag, optical_albedo);
   double snr_per_day =
            rdata->radar_constant * radar_albedo
            * sqrt( rotation_period * diameter_in_meters)
                        * diameter_in_meters / pow( dist, 4.);

   snr_per_day *= rdata->power_in_watts * rdata->gain;
   snr_per_day /= rdata->system_temp_deg_k;
   return( snr_per_day);
}

/* The following is a simple exercise in geometry:  if the sun and a
planet have specified angular radii,  and are separated in the sky
by 'sep' radians,  how much of the sun is blocked by the planet?
Simple answers are handled first;  if sep > r_sun + r_planet,  there's
no eclipse whatsoever,  and the planet is in full sunlight.  (Which
happens about 99+% of the time,  even for most artsats.  So a quick
check and returning 1 is all we usually need to do.)

   Next most common is the case where sep < r_planet - r_sun,  and
the object is totally eclipsed and we return 0;  followed by the case
where the earth/planet is fully on the sun's disk,  and the sunlight
is just decreased by that fraction.

   The rarest,  but most complicated,  case is where part of the planet's
disk,  but not all of it,  is blocking some sunlight.  Even then,  the
math to do the job is not rocket surgery... */

static double sunlight_visible( double r_sun, double r_planet, double sep)
{
   double rval;

   assert( sep > 0.);
   assert( r_planet > 0.);
   assert( r_sun > 0.);
   r_planet /= r_sun;      /* work in units of sun's apparent ang diam = 1 */
   sep /= r_sun;
   r_sun = 1.;
   if( sep >= r_planet + r_sun)
      return( 1.);               /* no actual eclipse,  not even partial */
   else if( sep <= r_planet - r_sun)
      return( 0.);               /* total eclipse */

   if( sep < r_sun - r_planet)         /* planet is transiting; */
      rval = 1. - r_planet * r_planet;      /* "annular eclipse" */
   else                  /* partial eclipse/transit */
      {
      double x, y;      /* solve for x^2+y^2=1, (x-sep)^2+y^2=r_planet^2 */
      double ang1, ang2;
      const double r_planet2 = r_planet * r_planet;

      x = (sep * sep + 1. - r_planet2) * .5 / sep;
      if( x >= 1.)     /* disks are barely touching; no real eclipse */
         return( 1.);
      if( x <= -1.)    /* disks are barely covering; totally eclipsed */
         return( 0.);
      y = sqrt( 1. - x * x);
      assert( y);
      ang1 = atan2( y, x);
      ang2 = atan2( y, x - sep);
      rval = 1 + (x * y - ang1) / PI;
      rval -= (1. - ang2 / PI) * r_planet2 + (x - sep) * y / PI;
      }
   return( rval);
}

/* See 'shadow.cpp' for an explanation of this function. */

static double solar_disc_intensity( const double r)
{
#ifdef PIERCE_SLAUGHTER_COEFFS
   const double I0 = 0.4067828571428571;   /* above integral for r=0, mu=1 */
#else
   const double I0 = 0.4062268095238096;   /* above integral for r=0, mu=1 */
#endif
   const double mu2 = 1. - r * r;
   double mu, rval = 0, power = mu2;
   size_t i;
#ifdef PIERCE_SLAUGHTER_COEFFS
   static const double coeffs[6] = { 0.30505, 1.13123, -0.78604,
                                    0.40560,  0.02297, -0.07880 };
#else       /* default Neckel & Labs coeffs */
   static const double coeffs[6] = { 0.28392, 1.36896, -1.75998,
                                    2.22154, -1.56074, 0.44630 };
#endif

   if( mu2 > 0.)        /* with roundoff,  we can be slightly negative */
      mu = sqrt( mu2);
   else
      mu = 0.;
   for( i = 0; i < sizeof( coeffs) / sizeof( coeffs[0]); i++, power *= mu)
      rval += coeffs[i] * power / (double)( i + 2);
// return( r * r);               /* for a uniformly intense solar disc */
   return( (I0 - rval) / I0);
}

/*  'shadow_check' determines how much sunlight is blocked by a planet
(thus far,  only the earth,  but could be extended to other planets and
their moons).  This is used in ephems to show "Sha" instead of a
magnitude if the object is in the earth's umbra;  to show a fainter
magnitude if it's in the penumbra; and to decrease solar radiation
pressure (SRP) in the force model as an object goes through the
earth's shadow;  see 'runge.cpp'.

   To determine if an object is in the planet's shadow:  first,  we
compute the elongation of the planet from the sun,  as seen from the
object.  If that's more than a right angle,  we can say that we're
on the daylit side and return 1.0 ("full sunlight").

   Otherwise,  we compute the angular separation between the sun and
planet as seen from the object,  and the angular radii of the sun
and planet.  We can then use the above 'sunlight_visible' function
to determine the fraction of the sun that is visible (not blocked),
almost always 100% (you don't spend much time in the earth's shadow).

   See http://www.minorplanet.info/MPB/issues/MPB_45-3.pdf for a discussion
of two instances where this matters.  2008 TC3 and 2018 LA,  impactors,
are other examples of objects passing through earth's shadow,  as are
numerous artsats.

   The partially eclipsed portion of the sun's disk is divided into a
series of rings,  equal in area, and the total intensity (after correcting
for limb darkening and fraction of the ring in eclipse) is added up.
Unless the earth's disk cuts exactly through the center of the sun, there
is usually a central part that is completely eclipsed or uneclipsed,
which is handled analytically.  See 'shadow.cpp' for details on how limb
darkening is computed.  */

double shadow_check( const double *planet_loc,
                            const double *obs_posn,
                            const double planet_radius_in_au)
{
   const double r2 = dot_product( planet_loc, planet_loc);
   const double dot = dot_product( planet_loc, obs_posn);
   const double SUN_RADIUS_IN_AU = 696000. / AU_IN_KM;
   double d = 0., angular_sep, r;
   double ang_size_sun, ang_size_planet;
   size_t i;

   if( dot < r2)     /* object is 'sunside' of the planet */
      return( 1.);
   for( i = 0; i < 3; i++)
      {
      const double diff = obs_posn[i] - planet_loc[i];

      d += diff * diff;
      }
   d = sqrt( d);        /* d = planet-object dist */
   if( d < planet_radius_in_au)     /* inside a planet, */
      return( 0.);                  /* it's dark */
   r = sqrt( r2);
   angular_sep = acos( (dot - r2) / (d * r));
   ang_size_sun = SUN_RADIUS_IN_AU / r;
   ang_size_planet = planet_radius_in_au / d;
   if( ang_size_sun + ang_size_planet < angular_sep)
      return( 1.);         /* no overlap */
   else if( ang_size_planet > ang_size_sun + angular_sep)
      return( 0.);         /* opposite extreme : total eclipse */
   else        /* partial eclipse */
      {
      double rval = 0., prev_fraction = 0., prev_intensity = 0.;
      const double r0 = fabs( angular_sep - ang_size_planet) / ang_size_sun;
      double prev_area = 0.;
      const size_t n_splits = 5;

      assert( r0 >= 0. && r0 <= 1.);
      for( i = 0; i <= n_splits; i++)
         {
         const double area = r0 * r0 + (1. - r0 * r0) * (double)i / (double)n_splits;
         const double r = sqrt( area);
         const double intensity = solar_disc_intensity( r);
         const double fraction =
                sunlight_visible( ang_size_sun * r, ang_size_planet, angular_sep) * r * r;

         rval += (intensity - prev_intensity) * (fraction - prev_fraction)
                        / (area - prev_area);
         prev_fraction = fraction;
         prev_intensity = intensity;
         prev_area = area;
         }
      return( rval);
      }
}

double vector_to_polar( double *lon, double *lat, const double *vector)
{
   const double r = vector3_length( vector);

   *lon = PI + atan2( -vector[1], -vector[0]);
   *lat =  asin( vector[2] / r);
   return( r);
}

static void format_motion( char *buff, const double motion)
{
   const char *motion_format;
   const double fabs_motion = fabs( motion);

   if( fabs_motion > 999999.)
      motion_format = "-------";
   else if( fabs_motion > 999.)
      motion_format = "%7.1f";
   else if( fabs_motion > 99.9)
      motion_format = "%7.2f";
   else
      motion_format = "%7.3f";
   snprintf( buff, 8, motion_format, motion);
}

#ifdef NOT_CURRENTLY_USED

   /* In displaying an SR scatterplot,  we might want to size it */
   /* such that,  say,  17 of 20 points appear,  with the three  */
   /* largest outliers omitted.  'find_cutoff_point' finds the   */
   /* range required to accomplish this.  It first finds the     */
   /* range required to accommodate _all_ the points,  then      */
   /* gradually reduces this until it's too small.  Then it goes */
   /* back to the size that _did_ accommodate at least that many */
   /* data points,  and returns that value.                      */

static inline double find_cutoff_point( const double *x, const double *y,
         const unsigned n_pts, const unsigned target_n_inside)
{
   double lim = 0.;
   unsigned i, n_found_inside = n_pts;

   for( i = 0; i < n_pts; i++)
      {
      if( lim < fabs( x[i]))
         lim = fabs( x[i]);
      if( lim < fabs( y[i]))
         lim = fabs( y[i]);
      }
   while( n_found_inside > target_n_inside)
      {
      lim /= 1.1;
      n_found_inside = 0;
      for( i = 0; i < n_pts; i++)
         if( fabs( x[i]) < lim && fabs( y[i]) < lim)
            n_found_inside++;
      }
   return( lim * 1.1);
}
#endif        /* #ifdef ENABLE_SCATTERPLOTS */

inline void calc_sr_dist_and_posn_ang( const DPT *ra_decs, const unsigned n_objects,
                     double *dist, double *posn_ang, FILE *ofile)
{
   unsigned i;
   double mean_x = 0., mean_y = 0.;
   const double ra0 = ra_decs[0].x, dec0 = ra_decs[0].y;
   double *x, *y;
   double sum_x2 = 0., sum_y2 = 0., sum_xy = 0.;
   double b, c, discrim, z1, z2;
   const double radians_to_arcsecs = 180. * 3600. / PI;

   x = (double *)calloc( n_objects * 2, sizeof( double));
   if( !x)
      return;
   y = x + n_objects;
   for( i = 1; i < n_objects; i++)
      {
      double dx = centralize_ang( ra_decs[i].x - ra0);
      double dy =                 ra_decs[i].y - dec0;

      if( dx > PI)
         dx -= PI + PI;
      dx *= cos( dec0) * radians_to_arcsecs;
      dy *=              radians_to_arcsecs;
      x[i] = dx;
      y[i] = dy;
      mean_x += dx;
      mean_y += dy;
      }
   mean_x /= (double)n_objects;
   mean_y /= (double)n_objects;
   for( i = 0; i < n_objects; i++)
      {
      const double dx = x[i] - mean_x;
      const double dy = y[i] - mean_y;

      sum_x2 += dx * dx;
      sum_xy += dx * dy;
      sum_y2 += dy * dy;
//    debug_printf( "  %u: %f %f (%f %f)\n", i,
//                   ra_decs[i].x * 180. / PI, ra_decs[i].y * 180. / PI,
//                   dx, dy);
      }                         /* We now have a covariance matrix: */
   sum_x2 /= (double)n_objects;        /*   / sum_x2 sum_xy \        */
   sum_xy /= (double)n_objects;        /*  (                 )       */
   sum_y2 /= (double)n_objects;        /*   \ sum_xy sum_y2 /        */
            /* The eigenvals of this are the zeroes z1, z2 of :      */
            /* z^2 - (sum_x2 + sum_y2)z + sum_x2*sum_y2 - (sum_xy)^2 */
   b = -(sum_x2 + sum_y2);
   c = sum_x2 * sum_y2 - sum_xy * sum_xy;
   discrim = b * b - 4. * c;
   assert( discrim >= 0.);
   z1 = (-b + sqrt( discrim)) * .5;
   z2 = c / z1;
   *dist = sqrt( z1);
   assert( z1 > z2);    /* z1 should be the _major_ axis */
   *posn_ang = atan2( sum_xy, sum_x2 - z2);
    if( *posn_ang < 0.)
      *posn_ang += PI;
#ifdef ENABLE_SCATTERPLOTS
// debug_printf( "Posn angles   : %f %f\n", *posn_ang * 180. / PI,
//          atan2( sum_xy, sum_x2 - z1) * 180. / PI);
// debug_printf( "Posn angles(2): %f %f\n",
//          atan2( sum_y2 - z1, sum_xy) * 180. / PI,
//          atan2( sum_y2 - z2, sum_xy) * 180. / PI);
   debug_printf( "Unc ellipse: %f x %f arcsec\n",
               sqrt( z1), sqrt( z2));
      {
      const unsigned ysize = 30, xsize = ysize * 2;
      char **histo = make_text_scattergram( x, y, n_objects, xsize, ysize,
                    find_cutoff_point( x, y, n_objects, n_objects * 9 / 10));

      for( i = 0; i <= ysize; i++)
         debug_printf( "%s\n", histo[ysize - i]);
      free( histo);
      }
#endif
   if( ofile)
      {
      fprintf( ofile, "# %u points; %.1f x %.1f ellipse at %.1f\n",
                     n_objects, *dist, sqrt( z2), 180. - *posn_ang * 180. / PI);
      for( i = 0; i < n_objects; i++)
         fprintf( ofile, "%.1f %.1f\n", x[i], y[i]);
      }
   *dist /= radians_to_arcsecs;
   free( x);
}

static inline void clean_up_json_number( char *out_text)
{
   size_t len = strlen( out_text);

   assert( len > 0);
               /* JSON numbers can't end with a decimal point... */
   if( out_text[len - 1] == '.')
      out_text[len - 1] = '\0';
               /* ...nor begin with a '+'... */
   if( *out_text == '+')
      memmove( out_text, out_text + 1, strlen( out_text));
   else if( *out_text == '-')    /* ...but a '-' is okay */
      out_text++;
               /* can't start with a zero,  unless before a decimal point */
   while( *out_text == '0' && (out_text[1] != '.' && out_text[1]))
      memmove( out_text, out_text + 1, strlen( out_text));
}

static int create_json_ephemeris( FILE *ofile, FILE *ifile, char *header,
                     const double jd_start, const double step)
{
   char buff[1024];
   int line_no = 0;

   text_search_and_replace( header, "-", "");
   text_search_and_replace( header, " RA ", " RA RA60 ");
   text_search_and_replace( header, " Dec ", " Dec Dec60 ");
   text_search_and_replace( header, "\"sigPA", "sigPos sigPA");
   text_search_and_replace( header, "ph_ang_bisector", "PABlon PABlat");
   text_search_and_replace( header, "topo_ecliptic", "topoLon topoLat");
   text_search_and_replace( header, "helio_ecliptic", "helioLon helioLat");
   text_search_and_replace( header, "RA'/hrdec", "RAvel decvel");
   text_search_and_replace( header, "'/hr PA", "motion_rate motionPA ");
   while( fgets_trimmed( buff, sizeof( buff), ifile))
      if( memcmp( buff, "....", 4) && *buff != '#')         /* skip irrelevant lines */
         {
         char *hptr = header, *bptr = buff;
         const char *preceder = (line_no ? ",\n" : "");
         const int step_number = (step ? (int)(( atof( buff) - jd_start) / step + .5) : line_no);

         fprintf( ofile, "%s      \"%d\": {", preceder, step_number);
         while( *hptr && *bptr)
            {
            int hlen = 0, blen = 0;
            char out_token[20], out_text[30];
            bool is_text = false;

            while( *bptr && *bptr == ' ')
               bptr++;
            while( *hptr && *hptr == ' ')
               hptr++;
            if( !*bptr && !*hptr)         /* end of line */
               break;
            if( !*bptr || !*hptr)
               {
               fprintf( ofile, "MISMATCH\n");
               return( -1);
               }
            while( hptr[hlen] && hptr[hlen] != ' ')
               hlen++;
            while( bptr[blen] && bptr[blen] != ' ')
               blen++;
            assert( hlen < 18);
            assert( blen < 25);
            if( hptr == header)     /* first token -> date/time */
               strcpy( out_token, "JD");
            else
               {
               if( hptr == header + 5)   /* second token -> ISO date */
                  {
                  strcpy( out_token, "ISO_time");
                  is_text = true;
                  }
               else
                  {
                  memcpy( out_token, hptr, hlen);
                  out_token[hlen] = '\0';
                  }
               fprintf( ofile, ", ");
               }
            memcpy( out_text, bptr, blen);
            out_text[blen] = '\0';
            if( out_token[hlen - 2] == '6' && out_token[hlen - 1] == '0')
               {        /* RA60 or Dec60 */
               is_text = true;
               text_search_and_replace( out_text, "_", " ");
               }
            if( !strcmp( out_token, "SM"))
               is_text = true;
            remove_trailing_cr_lf( out_text);
            if( is_text)      /* Text must be enclosed in quotes for JSON */
               {
               const size_t len = strlen( out_text);

               memmove( out_text + 1, out_text, len + 1);
               *out_text = out_text[len + 1] = '\"';
               out_text[len + 2] = '\0';
               }
            else            /* numeric quantities can't end in '.'... */
               clean_up_json_number( out_text);
            fprintf( ofile, "\"%s\": %s", out_token, out_text);
            hptr += hlen;
            bptr += blen;
            }
         fprintf( ofile, " }");
         line_no++;
         }
   if( line_no)
      fprintf( ofile, "\n");
   return( line_no);
}

#ifndef PATH_MAX
   #define PATH_MAX 256
#endif

/* Given a fully-specified filename,  this function will try to make sure
that a directory exists for it.  For example,   given a filename

/home/joe/z/k32/hi_there/thanx.txt

   it'll mkdir( "/home"), mkdir ("/home/joe"),  ... "/home/joe/z/k32/hi_there". */

static void make_path_available( const char *filename)
{
   char path[PATH_MAX];
   int i;

   for( i = 0; filename[i] && i < PATH_MAX; i++)
      {
      if( i && filename[i] == '/')
         {
         path[i] = '\0';
#ifdef _WIN32
         _mkdir( path);
#else
         mkdir( path, 0777);
#endif
         }
      path[i] = filename[i];
      }
}
#ifndef _WIN32
void fix_home_dir( char *filename)
{
   if( filename[0] == '~' && filename[1] == '/')
      text_search_and_replace( filename, "~", getenv( "HOME"));
}
#endif

char *real_packed_desig( char *obuff, const char *packed_id)
{
   strcpy( obuff, packed_id);
   text_search_and_replace( obuff, " ", "");
   if( strlen( obuff) == 12)      /* packed desig contains both perm */
      obuff[5] = '\0';            /* & prov IDs;  just use the perm */
   return( obuff);
}

static char ephem_mpc_code[70];

/* By default,  JSON files are kept with fixed names in the ~/.find_orb
directory.  However,  this can be overridden with parameters in
'environ.dat' (q.v.)  The following function opens JSON files as
specified therein (if you've actually done so) or in the default
locations (if you haven't.)  One can specify that the JSON filename
include the packed designation,  and/or the observatory code used for
the ephemeris (only works for ephems and 'combined'),  and/or a random
number based on the process ID (used for server Find_Orb to provide a
unique and semi-private ID for the JSON file). See
https://www.projectpluto.com/fo_usage#json_files for further
information. */

unsigned random_seed;

FILE *open_json_file( char *filename, const char *env_ptr, const char *default_name,
                  const char *packed_desig, const char *permits)
{
   char tbuff[100], full_permits[20];

#ifdef _WIN32
   strcpy( tbuff, "WIN_");
   strcat( tbuff, env_ptr);
   env_ptr = get_environment_ptr( tbuff);
#else
   env_ptr = get_environment_ptr( env_ptr);
#endif
   if( !*env_ptr)
      {
      get_file_name( filename, default_name);
      strcpy( full_permits, "fc");
      }
   else
      {
      extern const char *combine_all_observations;

      strcpy( filename, env_ptr);
      if( combine_all_observations && *combine_all_observations)
         strcpy( tbuff, combine_all_observations);
      else
         real_packed_desig( tbuff, packed_desig);
      text_search_and_replace( filename, "%p", tbuff);
      text_search_and_replace( filename, "%c", ephem_mpc_code);
      sprintf( tbuff, "%x", random_seed);
      text_search_and_replace( filename, "%r", tbuff);
#ifndef _WIN32
      fix_home_dir( filename);
#endif
      strcpy( full_permits, "f");
      make_path_available( filename);
      }
   strlcat( full_permits, permits, sizeof( full_permits));
   return( fopen_ext( filename, full_permits));
}

static int combine_json_elems_and_ephems( const char *packed_desig, FILE *ephem_file)
{
   char buff[1024];
   FILE *ofile, *elem_file;
   bool in_observations = false, obs_end_found = false;
   bool in_ephemerides = false;

   ofile = open_json_file( buff, "JSON_COMBINED_NAME", "combined.json", packed_desig, "wb");
   elem_file = open_json_file( buff, "JSON_ELEMENTS_NAME", "elements.json", packed_desig, "rb");
   while( !obs_end_found && fgets( buff, sizeof( buff), elem_file))
      {
      if( !strncmp( buff, "      \"observations\":", 17))
         in_observations = true;
      if( in_observations && !strncmp( buff, "      }", 7))
         {
         strcpy( buff, "      },\n");   /* JSON needs comma for 'continuation' */
         obs_end_found = true;
         }
      fwrite( buff, strlen( buff), 1, ofile);
      }
   fclose( elem_file);
   fseek( ephem_file, 0L, SEEK_SET);
   while( fgets( buff, sizeof( buff), ephem_file))
      {
      if( !strncmp( buff, "  \"ephemeris\":", 14))
         in_ephemerides = true;
      if( in_ephemerides)
         {
         fprintf( ofile, "    ");
         fwrite( buff, strlen( buff), 1, ofile);
         }
      }
   fprintf( ofile, "  }\n}");
   fclose( ofile);
   return( 0);
}

/* 'Close approach' ephemerides are generated by first computing a more
conventional ephemeris of distances to the object.  If we get three
positions in a row wherein the middle is less than the first and third,
then a minimum distance has been bracketed,  and we can use Brent's
minimization method (see the 'lunar' library for the code that does
this) to nail down the minimum.

   We look for the minimum of the _square_ of the distance,  because
that is (usually) just about quadratic in time,  and Brent's method
will converge a _lot_ faster that way.  */

static double find_closest_approach( const double *input_orbit, double jde,
                            const int planet_no, double *dist,
                            const double step, const double *prev_r)
{
   double orbit[6];
   int is_done = 0;
   brent_min_t b;

   memcpy( orbit, input_orbit, 6 * sizeof( double));
   brent_min_init( &b, jde, prev_r[0] * prev_r[0],
                jde - step, prev_r[1] * prev_r[1],
                jde - 2. * step, prev_r[2] * prev_r[2]);
   b.tolerance = 1e-6;
   while( !is_done)
      {
      double loc[3], r2, new_jde;
      int i;

      assert( b.n_iterations < 30);
      new_jde = brent_min_next( &b);
      integrate_orbit( orbit, jde, new_jde);
      jde = new_jde;
      compute_observer_loc( jde, planet_no, 0., 0., 0., loc);
      for( i = 0; i < 3; i++)
         loc[i] -= orbit[i];
      r2 = dot_product( loc, loc);
      is_done = brent_min_add( &b, r2);
      *dist = sqrt( r2);
      }
   return( jde);
}

static void add_lon_lat_to_ephem( char *buff, const size_t buflen,
               const double lon_in_radians, const double lat_in_radians)
{
   snprintf_append( buff, buflen, " %8.4f %8.4f", lon_in_radians * 180. / PI,
                                                 lat_in_radians * 180. / PI);
}

/* We have an image map,  currently based on Gaia-DR2,  of 'galactic
confusion'. See bright.cpp in my 'star_cats' repository for the code that
created this.  The image is in the .pgm (Portable GrayMap) format,
basically an array of byte values.  We convert RA/dec to a pixel
location,  and read that pixel value.  Some logic is added to save
the FILE* and cache a few bytes.

'bright2.pgm' has a one-degree resolution (and is therefore 360x180)
and is distributed by default with Find_Orb.  'bright.pgm' has ten times
the resolution,  but does require a correspondingly bigger download.
We look for the latter;  if it fails,  we try the former. */

int galactic_confusion( const double ra, const double dec)
{
   static FILE *image_file;
   static long buff_offset, hdr_offset;
   static int xsize, ysize;
   static unsigned char buff[32];
   int i, x, y;
   long loc;

   if( ra == -99.)       /* flag for "we're done here;  free everything up" */
      {
      if( image_file)
         fclose( image_file);
      image_file = NULL;
      xsize = 0;
      return( 0);
      }
   assert( ra >= 0. && ra <= 360.);
   assert( dec >= -90. && dec <= 90.);
   if( hdr_offset == -1)    /* we've already looked for a confusion image */
      return( 0);           /* and failed.  No need to fail again. */
   for( i = 0; !image_file && i < 2; i++)
      image_file = fopen_ext( i ? "bright2.pgm" : "bright.pgm", "crb");
   assert( image_file);
   if( !image_file)        /* tried both possible confusion images. */
      {
      hdr_offset = -1;
      return( 0);
      }
   if( !xsize)
      {
      char tbuff[40];

      if( fgets( tbuff, sizeof( buff), image_file))
         if( fgets( tbuff, sizeof( buff), image_file))
            {
            sscanf( tbuff, "%d %d", &xsize, &ysize);
            if( !fgets( tbuff, sizeof( buff), image_file))
               return( 0);   /* line w/max pixel value */
            assert( atoi( tbuff) == 255);
            hdr_offset = ftell( image_file);
            buff_offset = -99;
            }
      }
   assert( xsize);
   x = (int)( (360. - ra) * (double)xsize / 360.) % xsize;
   y = (int)( (90. - dec) * (double)ysize / 180.);
   loc = x + y * xsize;
   if( loc < buff_offset || loc >= buff_offset + (long)sizeof( buff))
      {
      size_t bytes_read;

      buff_offset = loc - sizeof( buff) / 2;
      if( buff_offset < 0)
         buff_offset = 0;
      fseek( image_file, buff_offset + hdr_offset, SEEK_SET);
      bytes_read = fread( buff, 1, sizeof( buff), image_file);
      assert( bytes_read >= sizeof( buff) / 2);
      }
   return( buff[loc - buff_offset]);
}

static double round_to( const double x, const double step)
{
   return( floor( x / step + .5) * step);
}


/* In 'auto-step' mode,  the step size text is a(number),  where (number)
sets an upper limit on how far the object is to move relative to the observer.
a.02 means don't move more than 2%;  a-.01 mean step backward such that the
object doesn't move by more than 1%.  In the latter case,  the object would
move by,  at most,  0.01 radians,  or about 34'.

We try to retain common step sizes such as ten minutes or two hours.  We
do this by figuring out how long it will take to move the desired fraction
(time = fraction * speed / distance),  and picking the 'common step size'
that gets us closest to that without going over.

As the object gets closer to us or goes away from us,  we should see the
auto-step vary.  */

double find_next_auto_step( const double target_diff, const bool going_backward,
                  const double curr_jd)
{
   static const double sizes[] = { 1. / 86400., 2. / 86400., 5. / 86400.,
            10. / 86400., 15. / 86400., 30. / 86400., 1. / 1440., 2. / 1440.,
            5. / 1440., 10. / 1440., 15. / 1440., 30. / 1440., 1. / 24.,
            2. / 24., 3. / 24., 4. / 24., 6. / 24., 12. / 24.,
            1., 2., 4., 8., 16., 32., 64. };
   size_t i;
   double max_diff = 0., rval = curr_jd;

   for( i = 0; i < sizeof( sizes) / sizeof( sizes[0]); i++)
      {
      double diff, new_t;

      if( going_backward)
         new_t = curr_jd - sizes[i] * .5 - .5;
      else
         new_t = curr_jd + sizes[i] * 1.5 - .5;
      if( sizes[i] > 0.9)
         new_t = sizes[i] * floor( new_t / sizes[i]);
      else        /* roundoff avoidance... */
         {
         const double int_t = floor( new_t);

         new_t = int_t + sizes[i] * floor( (new_t - int_t) / sizes[i]);
         }
      new_t += 0.5;
      diff = fabs( new_t - curr_jd);
      if( diff < target_diff && diff > max_diff && diff > sizes[i] * .9)
         {
         rval = new_t;
         max_diff = diff;
         }
      }
   return( rval);
}

static int compare_doubles( const void *a, const void *b, void *context)
{
   const double *a1 = (const double *)a;
   const double *b1 = (const double *)b;
   int rval = (*a1 > *b1 ? 1 : -1);
   const int *dir = (const int *)context;

   if( *dir < 0)
      rval = -rval;
   return( rval);
}

static double *list_of_ephem_times = NULL;

static int get_ephem_times_from_file( const char *filename)
{
   FILE *ifile;
   char buff[80];
   int n_times = 0, byte_offset = 0;
   double jd;
   int sort_order = 0, time_system = 0;

   while( *filename == ' ')
      filename++;
   ifile = fopen_ext( filename, "frb");
   while( fgets( buff, sizeof( buff), ifile))
      n_times++;
   list_of_ephem_times = (double *)calloc( n_times, sizeof( double));
   fseek( ifile, 0L, SEEK_SET);
   n_times = 0;
   while( fgets_trimmed( buff, sizeof( buff), ifile))
      if( (jd = get_time_from_string( 0., buff + byte_offset, FULL_CTIME_YMD, NULL)) > 1.)
         {
         if( time_system)   /* Input times are in TD;  cvt to UTC */
            jd -= td_minus_utc( jd) / seconds_per_day;
         list_of_ephem_times[n_times++] = jd;
         }
      else if( !memcmp( buff, "OPTION ", 7))
         switch( buff[7])
            {
            case 'T':         /* as in TD : convert times to UTC */
               time_system = 1;
               break;
            case 'A':         /* as in Ascending */
               sort_order = 1;
               break;
            case 'D':         /* as in Descending */
               sort_order = -1;
               break;
            case 'C':         /* as in (starting) Column */
               byte_offset = atoi( buff + 8) - 1;
               break;
            }
   if( sort_order)
      shellsort_r( list_of_ephem_times, n_times, sizeof( double),
                     compare_doubles, &sort_order);
   fclose( ifile);
   return( n_times);
}

#define RGB_OUTSIDE_POINTING_LIMITS      RGB_ORANGE
#define RGB_BELOW_HORIZON                RGB_DIRT

double ephemeris_mag_limit = 22.;

int ephemeris_in_a_file( const char *filename, const double *orbit,
         OBSERVE *obs, const int n_obs,
         const int planet_no,
         const double epoch_jd, const double jd_start, const char *stepsize,
         double lon, double rho_cos_phi, double rho_sin_phi,
         const int n_steps, const char *note_text,
         ephem_option_t options, unsigned n_objects)
{
   double *orbits_at_epoch, step;
   DPT *stored_ra_decs;
   double prev_ephem_t = epoch_jd, prev_r[3];
   int i, hh_mm, n_step_digits;
   int n_lines_shown = 0;
   unsigned date_format;
   const int ephem_type = ((int)(options & 7) == 6 ? 0 : (int)(options & 7));
   FILE *ofile, *computer_friendly_ofile = NULL;
   const bool computer_friendly = ((options & OPTION_COMPUTER_FRIENDLY) ? true : false);
   char step_units;
   const char *timescale = get_environment_ptr( "TT_EPHEMERIS");
   const char *override_date_format = get_environment_ptr( "DATE_FORMAT");
   double abs_mag = calc_absolute_magnitude( obs, n_obs);
   double ht_in_meters = 0., max_auto_step = 0.;
   DPT latlon;
   bool last_line_shown = true;
   RADAR_DATA rdata;
   bool show_radar_data = (get_radar_data( note_text + 1, &rdata) == 0);
   const double planet_radius_in_au =
          planet_radius_in_meters( planet_no) / AU_IN_METERS;
   int ra_format = 3, dec_format = 2;
   char buff[440], *header = NULL, alt_buff[500];
   const bool use_observation_times = !strncmp( stepsize, "Obs", 3);
   double curr_jd = jd_start, real_jd_start = jd_start;
   bool reset_lat_alt = true;

   if( (!rho_cos_phi && !rho_sin_phi && !use_observation_times)
                                 || ephem_type != OPTION_OBSERVABLES)
      options &= ~(OPTION_ALT_AZ_OUTPUT | OPTION_VISIBILITY | OPTION_MOON_ALT
                     | OPTION_MOON_AZ | OPTION_SUN_ALT | OPTION_SUN_AZ
                     | OPTION_SKY_BRIGHTNESS | OPTION_SUPPRESS_UNOBSERVABLE);
   if( n_objects == 1 || ephem_type != OPTION_OBSERVABLES)
      options &= ~OPTION_SHOW_SIGMAS;

   sscanf( get_environment_ptr( "RA_DEC_FORMAT"), "%d,%d", &ra_format, &dec_format);
   if( !use_observation_times && *stepsize != 't')
      {
      step = get_step_size( stepsize, &step_units, &n_step_digits);
      if( !step && *stepsize != 'a')
         return( -2);
      }
   else
      {
      step_units = 'd';
      n_step_digits = 6;
      step = 0.;
      }
   ofile = fopen_ext( filename, is_default_ephem ? "tfcw" : "fw");
   if( !ofile)
      return( -1);
   if( !memcmp( note_text, "(CSS)", 5))
      {
      double min_jd = jd_start;
      double max_jd = jd_start + step * (double)n_steps;
      const char *precovery_header_line =
               "    RA (J2000) dec  Mag  YYYY MM DD HH:MM:SS.s Code"
               " Sigma  PA Prob   Directory  Image Filename\n";
      int rval = -1;

      if( max_jd < min_jd)
         {
         min_jd = max_jd;
         max_jd = jd_start;
         }
      setvbuf( ofile, NULL, _IONBF, 0);
      fprintf( ofile, "#CSS precovery fields\n");
      fprintf( ofile, "%s", precovery_header_line);
      for( i = 0; i < 2; i++)
         {
         const int err_val = find_precovery_plates( obs, n_obs,
                           (i ? "css.idx" : "css_new.idx"),
                           ofile, orbit,
                           n_objects, epoch_jd, min_jd, max_jd,
                           ephemeris_mag_limit);

         if( !err_val)           /* at least one index worked */
            rval = 0;
         }
      fclose( ofile);
      return( rval);
      }
   if( planet_no < 0)      /* bad observatory code or satellite */
      return( -3);
   if( !abs_mag)
      abs_mag = atof( get_environment_ptr( "ABS_MAG"));
   if( ephem_type == OPTION_OBSERVABLES && !(options & OPTION_SHOW_SIGMAS))
      n_objects = 1;
   orbits_at_epoch = (double *)calloc( n_objects, 8 * sizeof( double));
   memcpy( orbits_at_epoch, orbit, n_objects * 6 * sizeof( double));
   stored_ra_decs = (DPT *)( orbits_at_epoch + 6 * n_objects);
   setvbuf( ofile, NULL, _IONBF, 0);
   switch( step_units)
      {
      case 'd':
         hh_mm = 0;
         date_format = FULL_CTIME_DATE_ONLY;
         break;
      case 'h':
         hh_mm = 1;
         date_format = FULL_CTIME_FORMAT_HH;
         break;
      case 'm':
         hh_mm = 2;
         date_format = FULL_CTIME_FORMAT_HH_MM;
         break;
      case 's':
      default:
         hh_mm = 3;
         date_format = FULL_CTIME_FORMAT_SECONDS;
         break;
      }
   date_format |= FULL_CTIME_YEAR_FIRST | FULL_CTIME_MONTH_DAY
            | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_LEADING_ZEROES
            | CALENDAR_JULIAN_GREGORIAN;
   date_format |= (n_step_digits << 4);
                     /* For ease of automated processing,  people may want */
   if( *override_date_format)     /* the time in some consistent format... */
      sscanf( override_date_format, "%x", &date_format);
   if( step && (options & OPTION_ROUND_TO_NEAREST_STEP))
       real_jd_start = round_to( jd_start - .5, step) + .5;
   if( ephem_type == OPTION_STATE_VECTOR_OUTPUT ||
       ephem_type == OPTION_POSITION_OUTPUT ||
       ephem_type == OPTION_MPCORB_OUTPUT ||
       ephem_type == OPTION_8_LINE_OUTPUT)
      {
      timescale = "y";        /* force TT output */
      fprintf( ofile, "%.5f %f %d %s %s\n", real_jd_start, step, n_steps,
                     get_environment_ptr( "VECTOR_OPTS"), note_text);
      }
   else if( ephem_type != OPTION_CLOSE_APPROACHES)
      {
      char hr_min_text[80], added_prec_text_ra[20], added_prec_text_dec[20];
      const char *pre_texts[4] = { "", "-HH", "-HH:MM", "-HH:MM:SS" };

      output_angle_to_buff( added_prec_text_ra, 0., ra_format);
      output_signed_angle_to_buff( added_prec_text_dec, 0., dec_format);
      remove_trailing_cr_lf( added_prec_text_ra);
      remove_trailing_cr_lf( added_prec_text_dec);
      memset( added_prec_text_ra,  '-', strlen( added_prec_text_ra));
      memset( added_prec_text_dec, '-', strlen( added_prec_text_dec));
      strcpy( hr_min_text, pre_texts[hh_mm]);
      if( n_step_digits)
         {
         strcat( hr_min_text, ".");
         for( i = n_step_digits; i; i--)
            {
            char tbuff[2];

            tbuff[0] = step_units;
            tbuff[1] = '\0';
            strcat( hr_min_text, tbuff);
            }
         }
      if( note_text)
         fprintf( ofile, "#%s\n", note_text);
      snprintf( buff, sizeof( buff), "Date %s%s  ",
                     (*timescale ? "(TT)"  : "(UTC)"), hr_min_text);
      if( !(options & OPTION_SUPPRESS_RA_DEC))
         snprintf_append( buff, sizeof( buff), "-RA%s   -Dec%s  ",
                                    added_prec_text_ra + 3, added_prec_text_dec + 4);
      if( !(options & OPTION_SUPPRESS_DELTA))
         snprintf_append( buff, sizeof( buff), "-delta ");
      if( !(options & OPTION_SUPPRESS_SOLAR_R))
         snprintf_append( buff, sizeof( buff), "-r---- ");
      if( !(options & OPTION_SUPPRESS_ELONG))
         snprintf_append( buff, sizeof( buff), "elong ");
      if( options & OPTION_VISIBILITY)
         snprintf_append( buff, sizeof( buff), "SM ");
      if( options & OPTION_SKY_BRIGHTNESS)
         snprintf_append( buff, sizeof( buff), "SkyBr ");
      if( options & OPTION_SNR)
         snprintf_append( buff, sizeof( buff), " -SNR ");
      if( options & OPTION_EXPOSURE_TIME)
         snprintf_append( buff, sizeof( buff), " ExpT ");
      if( options & OPTION_PHASE_ANGLE_OUTPUT)
         snprintf_append( buff, sizeof( buff), " ph_ang  ");
      if( options & OPTION_PHASE_ANGLE_BISECTOR)
         snprintf_append( buff, sizeof( buff), " ph_ang_bisector  ");
      if( options & OPTION_HELIO_ECLIPTIC)
         snprintf_append( buff, sizeof( buff), " helio_ecliptic   ");
      if( options & OPTION_TOPO_ECLIPTIC)
         snprintf_append( buff, sizeof( buff), " topo_ecliptic    ");
      if( options & OPTION_GALACTIC_COORDS)
         snprintf_append( buff, sizeof( buff), "Gal_Lon- Gal_Lat- ");
      if( options & OPTION_GALACTIC_CONFUSION)
         snprintf_append( buff, sizeof( buff), "GC ");
      if( options & OPTION_SUN_TARGET_PA)
         snprintf_append( buff, sizeof( buff), "-PsAng- ");
      if( options & OPTION_SUN_HELIO_VEL_PA)
         snprintf_append( buff, sizeof( buff), "-PsAMV- ");
      if( options & OPTION_ORBIT_PLANE_ANGLE)
         snprintf_append( buff, sizeof( buff), "-PlAng- ");
      if( abs_mag)
         snprintf_append( buff, sizeof( buff), " mag");

      if( options & OPTION_LUNAR_ELONGATION)
         snprintf_append( buff, sizeof( buff), "  LuElo");
      if( options & OPTION_MOTION_OUTPUT)
         snprintf_append( buff, sizeof( buff),
                (options & OPTION_SEPARATE_MOTIONS) ? " -RA--'/hr--dec-"
                                                    : " -'/hr-- --PA--");
      if( options & OPTION_ALT_AZ_OUTPUT)
         snprintf_append( buff, sizeof( buff), " alt -az");
      if( options & OPTION_SUN_ALT)
         snprintf_append( buff, sizeof( buff), " Sal");
      if( options & OPTION_SUN_AZ)
         snprintf_append( buff, sizeof( buff), " Saz");
      if( options & OPTION_MOON_ALT)
         snprintf_append( buff, sizeof( buff), " Mal");
      if( options & OPTION_MOON_AZ)
         snprintf_append( buff, sizeof( buff), " Maz");
      if( options & OPTION_RADIAL_VEL_OUTPUT)
         snprintf_append( buff, sizeof( buff), "  rvel-");
      if( show_radar_data)
         snprintf_append( buff, sizeof( buff), "  SNR");
      if( options & OPTION_GROUND_TRACK)
         snprintf_append( buff, sizeof( buff), " -lon---- -lat---- -alt-(km)-");
      if( options & OPTION_SPACE_VEL_OUTPUT)
         snprintf_append( buff, sizeof( buff), "  svel-");
      if( options & OPTION_SHOW_SIGMAS)
         snprintf_append( buff, sizeof( buff), " \"-sig-PA");
      if( ephem_type == OPTION_OBSERVABLES)
         {
         const char *cf_filename = "eph_json.txt";

         header = (char *)malloc( 1024);
         assert( header);
         strcpy( header, buff);
         computer_friendly_ofile = fopen_ext( cf_filename, is_default_ephem ? "tfcw+" : "fw+");
         assert( computer_friendly_ofile);
         }
      if( show_radar_data)
         exposure_config.min_alt = rdata.altitude_limit * 180. / PI;
      if( !computer_friendly)
         {
         if( show_radar_data)
            {
            fprintf( ofile,
                  "Assumes power=%.2f kW, Tsys=%.1f deg K, gain %.2f K/Jy\n",
                  rdata.power_in_watts / 1000., rdata.system_temp_deg_k,
                  rdata.gain);
            fprintf( ofile,
               "Assumed rotation period = %.2f hours, diameter %.1f meters\n",
                           guessed_rotation_period_in_hours( abs_mag),
                           diameter_from_abs_mag( abs_mag, optical_albedo));
            }
         for( i = 0; buff[i]; i++)
            if( buff[i] == '-')
               fprintf( ofile, " ");
            else
               {
               fprintf( ofile, "%c", buff[i]);
               if( buff[i] != ' ')
                  buff[i] = '-';
               }
         fprintf( ofile, "\n%s\n", buff);
         }
      }

   prev_r[0] = prev_r[1] = 0.;
   latlon.x = lon;
   for( i = 0; i < n_steps; i++)
      {
      unsigned obj_n;
      bool show_this_line = true;
      double ephemeris_t, utc;
      double obs_posn[3], obs_vel[3];
      double obs_posn_equatorial[3];
      double geo_posn[3], geo_vel[3];
      double delta_t;
      long rgb = 0;

      if( use_observation_times)
         {
         if( !i || strcmp( obs[i].mpc_code, obs[i - 1].mpc_code))
            {
            get_observer_data( obs[i].mpc_code, buff, &lon,
                                             &rho_cos_phi, &rho_sin_phi);
            latlon.x = lon;
            reset_lat_alt = true;
            }
         curr_jd = obs[i].jd;
         }
      else if( *stepsize == 't')
         curr_jd = list_of_ephem_times[i];
      else if( *stepsize == 'a')
         {
         if( i)
            curr_jd = find_next_auto_step( max_auto_step,
                                 stepsize[1] == '-', curr_jd);
         }
      else
         {
         curr_jd = jd_start + (double)i * step;
         if( options & OPTION_ROUND_TO_NEAREST_STEP)
            curr_jd = round_to( curr_jd - .5, step) + .5;
         }
      if( reset_lat_alt)
         parallax_to_lat_alt(
                  rho_cos_phi / planet_radius_in_au,
                  rho_sin_phi / planet_radius_in_au, &latlon.y,
                                      &ht_in_meters, planet_no);
      reset_lat_alt = false;
      delta_t = td_minus_utc( curr_jd) / seconds_per_day;
      if( use_observation_times)
         curr_jd -= delta_t;
      if( *timescale)                     /* we want a TT ephemeris */
         {
         ephemeris_t = curr_jd;
         utc = curr_jd - delta_t;
         }
      else                                /* "standard" UTC ephemeris */
         {
         ephemeris_t = curr_jd + delta_t;
         utc = curr_jd;
         }
      compute_observer_loc( ephemeris_t, planet_no, rho_cos_phi, rho_sin_phi,
                              lon, obs_posn);
      compute_observer_vel( ephemeris_t, planet_no, rho_cos_phi, rho_sin_phi,
                              lon, obs_vel);
      compute_observer_loc( ephemeris_t, planet_no, 0., 0., 0., geo_posn);
      compute_observer_vel( ephemeris_t, planet_no, 0., 0., 0., geo_vel);
                /* we need the observer position in equatorial coords too: */
      memcpy( obs_posn_equatorial, obs_posn, 3 * sizeof( double));
      ecliptic_to_equatorial( obs_posn_equatorial);
      strcpy( buff, "Nothing to see here... move along... uninteresting... who cares?...");
      for( obj_n = 0; obj_n < n_objects; obj_n++)
         {
         double *orbi = orbits_at_epoch + obj_n * 6;
         double radial_vel, v_dot_r;
         double topo[3], topo_vel[3], geo[3], r;
         double topo_ecliptic[3];
         double orbi_after_light_lag[6];
         OBSERVE temp_obs;
         int j;

         integrate_orbit( orbi, prev_ephem_t, ephemeris_t);
         for( j = 0; j < 3; j++)
            {
            topo[j] = orbi[j] - obs_posn[j];
            geo[j] = orbi[j] - geo_posn[j];
            topo_vel[j] = orbi[j + 3] - obs_vel[j];
            }
         r = vector3_length( topo);
                    /* for "ordinary ephemeris" (not state vectors or */
                    /* orbital elements),  include light-time lag:    */
         if( ephem_type == OPTION_OBSERVABLES)
            {
            double derivs[6], dt = -r / AU_PER_DAY;

            calc_derivatives( ephemeris_t, orbi, derivs, 0);
            for( j = 0; j < 6; j++)
               {
               const double diff = derivs[j] * dt;

               orbi_after_light_lag[j] = orbi[j] + diff;
               if( j < 3)
                  {
                  topo[j] += diff;
                  geo[j] += diff;
                  }
               else
                  topo_vel[j - 3] += diff;
               }
            }

         memset( &temp_obs, 0, sizeof( OBSERVE));
         temp_obs.r = vector3_length( topo);
         for( j = 0; j < 3; j++)
            {
            temp_obs.vect[j] = topo[j] / r;
            temp_obs.obs_vel[j] = -topo_vel[j];
            }
         memcpy( temp_obs.obs_posn, obs_posn, 3 * sizeof( double));
         for( j = 0; j < 3; j++)
            temp_obs.obj_posn[j] = temp_obs.obs_posn[j] + topo[j];
                     /* rotate topo from ecliptic to equatorial */
         memcpy( topo_ecliptic, topo, 3 * sizeof( double));
         ecliptic_to_equatorial( topo);                           /* mpc_obs.cpp */
         ecliptic_to_equatorial( geo);
         ecliptic_to_equatorial( topo_vel);
         v_dot_r = 0.;
         for( j = 0; j < 3; j++)
            v_dot_r += topo[j] * topo_vel[j];
         r = vector3_length( topo);
         radial_vel = v_dot_r / r;
         if( *stepsize == 'a')
            max_auto_step = fabs( atof( stepsize + 1)) * r / vector3_length( topo_vel);
         if( ephem_type == OPTION_STATE_VECTOR_OUTPUT ||
               ephem_type == OPTION_POSITION_OUTPUT)
            {
            int ecliptic_coords = 0;      /* default to equatorial J2000 */
            double posn_mult = 1., vel_mult = 1.;     /* default to AU & AU/day */
            double tval = 1.;
            char format_text[20];

            snprintf( buff, sizeof( buff), "%.5f", curr_jd);
            sscanf( get_environment_ptr( "VECTOR_OPTS"), "%d,%lf,%lf",
                        &ecliptic_coords, &posn_mult, &tval);
            assert( tval);
            assert( posn_mult);
            vel_mult = posn_mult / tval;
            if( ecliptic_coords)
               {
               equatorial_to_ecliptic( topo);
               equatorial_to_ecliptic( topo_vel);
               }
            for( j = 10, tval = posn_mult; tval > 1.2; j--)
               tval /= 10.;
            snprintf( format_text, sizeof( format_text), "%%18.%df", j + 3);
            for( j = 0; j < 3; j++)
               snprintf_append( buff, sizeof( buff), format_text,
                                 topo[j] * posn_mult);
            if( ephem_type == OPTION_STATE_VECTOR_OUTPUT)
               {
               strcat( buff, " ");
               for( j = 12, tval = vel_mult; tval > 1.2; j--)
                  tval /= 10.;
               snprintf( format_text, sizeof( format_text), "%%16.%df", j);
               for( j = 0; j < 3; j++)
                  snprintf_append( buff, sizeof( buff), format_text,
                                 topo_vel[j] * vel_mult);
               }
            }
         else if( ephem_type == OPTION_8_LINE_OUTPUT
                || ephem_type == OPTION_MPCORB_OUTPUT)
            {
            if( !obj_n)
               {
               FILE *ifile;
                              /* only show comment data on last line */
               const int output_options = // ELEM_OUT_HELIOCENTRIC_ONLY |
                            (i == n_steps - 1 ? 0 : ELEM_OUT_NO_COMMENT_DATA);

               write_out_elements_to_file( orbi, ephemeris_t, ephemeris_t,
                          obs, n_obs, "", 5, 0, output_options);
               ifile = fopen_ext( get_file_name( buff,
                        (ephem_type == OPTION_8_LINE_OUTPUT) ? elements_filename : "mpc_fmt.txt"),
                        "tfcrb");
               while( fgets( buff, sizeof( buff), ifile))
                  fputs( buff, ofile);
               fclose( ifile);
               }
            *buff = '\0';
            show_this_line = last_line_shown = false;
            }
         else if( ephem_type == OPTION_CLOSE_APPROACHES)
            {
            prev_r[2] = prev_r[1];
            prev_r[1] = prev_r[0];
            prev_r[0] = vector3_length( geo);
            if( i >= 2 && prev_r[1] < prev_r[2] && prev_r[1] < prev_r[0])
               {
               char date_buff[80];
               double dist;
               const double max_close_approach = 0.5;
               const double jd = find_closest_approach( orbi, ephemeris_t, planet_no,
                                    &dist, step, prev_r);

               if( dist < max_close_approach)
                  {
                  full_ctime( date_buff, jd,
                        FULL_CTIME_FORMAT_SECONDS
                      | FULL_CTIME_YEAR_FIRST | FULL_CTIME_MONTH_DAY
                      | FULL_CTIME_MONTHS_AS_DIGITS
                      | FULL_CTIME_LEADING_ZEROES);
                  snprintf( buff, sizeof( buff), "Close approach at %s: ",
                                 date_buff);
                  format_dist_in_buff( buff + strlen( buff), dist);
                  fprintf( ofile, "%s\n", buff);
                  }
               }
            *buff = '\0';
            last_line_shown = show_this_line = false;
            }
         else if( ephem_type == OPTION_OBSERVABLES)
            {
            DPT ra_dec;
            const bool fake_astrometry = ((options & 7) == OPTION_FAKE_ASTROMETRY);
            char tbuff[80], *alt_tptr;

            ra_dec.x = atan2( topo[1], topo[0]);
            ra_dec.y = asin( topo[2] / r);
            stored_ra_decs[obj_n] = ra_dec;
            if( n_objects > 1 && obj_n == n_objects - 1 && show_this_line)
               {
               double dist, posn_ang;
               int int_pa;

               if( n_objects == 2)
                  calc_dist_and_posn_ang( (const double *)&stored_ra_decs[0],
                                       (const double *)&ra_dec,
                                       &dist, &posn_ang);
               else
                  {
                  const char *offset_dir = get_environment_ptr( "OFFSET_FILES");
                  FILE *offset_ofile = NULL;

                  if( *offset_dir)
                     {
                     char date_buff[80];
                     static bool path_already_made = false;

                     strcpy( tbuff, offset_dir);
                     strcat( tbuff, "/");
                     if( !path_already_made)
                        make_path_available( tbuff);
                     path_already_made = true;
                     full_ctime( date_buff, curr_jd,
                           FULL_CTIME_FORMAT_HH_MM | FULL_CTIME_YMD
                         | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_NO_SPACES
                         | FULL_CTIME_NO_COLONS | FULL_CTIME_LEADING_ZEROES);
                     strcat( tbuff, date_buff);
                     strcat( tbuff, ".off");
                     offset_ofile = fopen( tbuff, "wb");
                     assert( offset_ofile);
                     full_ctime( date_buff, curr_jd,
                            FULL_CTIME_FORMAT_HH_MM | FULL_CTIME_YMD);
                     fprintf( offset_ofile, "# JD %f = %s\n", curr_jd, date_buff);
                     fprintf( offset_ofile, "# %s\n", obs->packed_id);
                     }
                  calc_sr_dist_and_posn_ang( stored_ra_decs, n_objects,
                                       &dist, &posn_ang, offset_ofile);
                  if( offset_ofile)
                     fclose( offset_ofile);
                  }
               tbuff[0] = ' ';
               int_pa = put_ephemeris_posn_angle_sigma( tbuff + 1, dist, posn_ang, false);
               strcat( buff, tbuff);
               snprintf_append( alt_buff, sizeof( alt_buff), " %8.3f %3d",
                                    dist * 3600. * 180. / PI, int_pa);
               }
            if( !obj_n)
               {
               const double dec = ra_dec.y * 180. / PI;
               double ra = ra_dec.x * 12. / PI;
               double lunar_elong = 99., dist_moon = 99.;  /* values are ignored */
               char fake_line[81];
               double cos_elong, solar_r, elong;
               double fraction_illum = 1.;    /* i.e.,  not in earth's shadow */
               double mags_per_arcsec2 = 99.99;    /* sky brightness */
               DPT alt_az[3], sun_ra_dec;
               double earth_r = 0., hour_angle[3];
               char ra_buff[80], dec_buff[80];
               double phase_ang, curr_mag, air_mass = 40.;
               bool output_ra_in_degrees = (ra_format >= 100 && ra_format < 110);
               char visibility_char = ' ';

               solar_r = vector3_length( orbi_after_light_lag);
               earth_r = vector3_length( obs_posn_equatorial);
               cos_elong = r * r + earth_r * earth_r - solar_r * solar_r;
               cos_elong /= 2. * earth_r * r;
               elong = acose( cos_elong);
               if( ra < 0.) ra += 24.;
               if( ra >= 24.) ra -= 24.;
               output_angle_to_buff( ra_buff, ra * (output_ra_in_degrees ? 15. : 1.),
                                               ra_format);
               remove_trailing_cr_lf( ra_buff);
               for( j = 0; j < 3; j++)    /* compute alt/azzes of object (j=0), */
                  {                       /* sun (j=1), and moon (j=2)          */
                  DPT obj_ra_dec = ra_dec;

                  if( j)
                     {
                     int k;
                     double vect[3], earth_loc[3];

                     if( j == 1)    /* solar posn */
                        for( k = 0; k < 3; k++)
                           vect[k] = -obs_posn_equatorial[k];
                     else
                        {           /* we want the lunar posn */
                        earth_lunar_posn( ephemeris_t, earth_loc, vect);
                        for( k = 0; k < 3; k++)
                           vect[k] -= earth_loc[k];
                        cos_elong = dot_product( earth_loc, vect)
                                 / (vector3_length( earth_loc) * vector3_length( vect));
                        lunar_elong = acose( -cos_elong);
                        ecliptic_to_equatorial( vect);   /* mpc_obs.cpp */
                        fraction_illum = shadow_check( earth_loc, orbi_after_light_lag,
                                    EARTH_MAJOR_AXIS_IN_AU);
                        cos_elong = dot_product( vect, geo)
                                 / (vector3_length( vect) * vector3_length( geo));
                        dist_moon = acose( cos_elong);
                        }
                     vector_to_polar( &obj_ra_dec.x, &obj_ra_dec.y, vect);
                     if( j == 1)
                        sun_ra_dec = obj_ra_dec;
                     }
                  obj_ra_dec.x = -obj_ra_dec.x;
                  full_ra_dec_to_alt_az( &obj_ra_dec, &alt_az[j], NULL, &latlon, utc,
                                              &hour_angle[j]);
                  alt_az[j].x = centralize_ang( alt_az[j].x + PI);
                  }
               if( is_under_horizon( alt_az[0].y * 180. / PI,
                                     alt_az[0].x * 180. / PI, &exposure_config))
                  {
                  visibility_char = 'a';
                  rgb = RGB_OUTSIDE_POINTING_LIMITS;
                  }
               else if( elong < exposure_config.min_elong * PI / 180. ||
                        elong > exposure_config.max_elong * PI / 180.)
                  {
                  visibility_char = 'e';
                  rgb = RGB_OUTSIDE_POINTING_LIMITS;
                  }
               else if( dec < exposure_config.min_dec ||
                        dec > exposure_config.max_dec)
                  {
                  visibility_char = 'd';
                  rgb = RGB_OUTSIDE_POINTING_LIMITS;
                  }
               else if( -hour_angle[0] < exposure_config.min_ha * PI / 180. ||
                        -hour_angle[0] > exposure_config.max_ha * PI / 180.)
                  {
                  visibility_char = 'h';
                  rgb = RGB_OUTSIDE_POINTING_LIMITS;
                  }
               if( alt_az[0].y < 0. && visibility_char != ' ')
                  {
                  visibility_char = 'B';
                  rgb = RGB_BELOW_HORIZON;
                  }
               if( alt_az[0].y > 0.)
                  {
                  BRIGHTNESS_DATA bdata;

                  bdata.ht_above_sea_in_meters = ht_in_meters;
                  bdata.latitude = latlon.y;
                  bdata.temperature_in_c = 20.;      /* centigrade */
                  bdata.moon_elongation = lunar_elong;
                  bdata.relative_humidity = 20.;  /* 20% */
                  bdata.year = 2000. + (ephemeris_t - J2000) / 365.25;
                  bdata.month = bdata.year * 12. - floor( bdata.year * 12.);
                  bdata.dist_moon = dist_moon;
                  bdata.dist_sun = elong;
                  bdata.zenith_angle    = PI / 2. - alt_az[0].y;
                  bdata.zenith_ang_sun  = PI / 2. - alt_az[1].y;
                  bdata.zenith_ang_moon = PI / 2. - alt_az[2].y;
                  bdata.mask = 31;
                  set_brightness_params( &bdata);
                  compute_sky_brightness( &bdata);
                  rgb = 0;
                  for( j = 0; j < 3; j++)
                     {
                     double component = log10( bdata.brightness[j + 1]);
                     static double dark[3] = { -13.5, -13, -12.4 };
                     static double light[3] = { -7, -7, -5 };

                     component = (component - dark[j]) * 255. / (light[j] - dark[j]);
                     if( component > 255.)
                        component = 255.;
                     if( component > 0.)
                        rgb |= (unsigned)component << (j * 8);
                     }
                  mags_per_arcsec2 = -2.5 * log10( bdata.brightness[3]) - 11.055;  /* R brightness */
                  air_mass = bdata.air_mass;
                  }
               output_signed_angle_to_buff( dec_buff, dec, dec_format);
               remove_trailing_cr_lf( dec_buff);
               if( fake_astrometry)
                  {
                  strcpy( fake_line, obs->packed_id);
                  strcpy( fake_line + 12, "  C");
                  full_ctime( fake_line + 15, curr_jd + 5e-7,
                           FULL_CTIME_MICRODAYS | FULL_CTIME_YMD
                         | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_LEADING_ZEROES);
                  strcat( fake_line, ra_buff);
                  strcat( fake_line, dec_buff);
                  strcat( fake_line, "        ");      /* columns 57 to 65 */
                  }
               snprintf( alt_buff, sizeof( alt_buff), "%15.7f ", curr_jd);
               if( computer_friendly)
                  strlcpy( buff, alt_buff, sizeof( buff));
               else
                  {
                  full_ctime( buff, curr_jd, date_format);
                  strcat( buff, " ");
                  }
               iso_time( alt_buff + strlen( alt_buff), curr_jd, 0);
               if( !(options & OPTION_SUPPRESS_RA_DEC))
                  {
                  if( computer_friendly)
                     snprintf_append( buff, sizeof( buff), " %12.8f %12.8f",
                                    ra * 15, dec);
                  else
                     {
                     strcat( buff, " ");
                     strcat( buff, ra_buff);
                     strcat( buff, "   ");
                     strcat( buff, dec_buff);
                     strcat( buff, " ");
                     }
                  text_search_and_replace( ra_buff, " ", "_");
                  text_search_and_replace( dec_buff, " ", "_");
                  snprintf_append( alt_buff, sizeof( alt_buff), " %12.8f %s %12.8f %s",
                                    ra * 15, ra_buff, dec, dec_buff);
                  }
               if( !(options & OPTION_SUPPRESS_DELTA))
                  {
                  alt_tptr = alt_buff + strlen( alt_buff);
                  snprintf_append( alt_buff, sizeof( alt_buff), " %14.9f", r);
                  if( computer_friendly)
                     strlcat( buff, alt_tptr, sizeof( buff));
                  else
                     {
                       /* the radar folks prefer the distance to be always in */
                       /* AU,  w/no switch to km for close approach objects: */
                     use_au_only = show_radar_data;
                     format_dist_in_buff( buff + strlen( buff), r);
                     use_au_only = false;
                     }
                  }
               if( !(options & OPTION_SUPPRESS_SOLAR_R))
                  {
                  alt_tptr = alt_buff + strlen( alt_buff);
                  snprintf_append( alt_buff, sizeof( alt_buff), " %14.9f", solar_r);
                  if( computer_friendly)
                     strlcat( buff, alt_tptr, sizeof( buff));
                  else
                     format_dist_in_buff( buff + strlen( buff), solar_r);
                  }
               if( !(options & OPTION_SUPPRESS_ELONG))
                  {
                  snprintf_append( buff, sizeof( buff), " %5.1f", elong * 180. / PI);
                  snprintf_append( alt_buff, sizeof( alt_buff), " %9.5f", elong * 180. / PI);
                  }

               if( options & OPTION_VISIBILITY)
                  {
                  tbuff[0] = tbuff[2] = ' ';
                  tbuff[1] = visibility_char;
                  if( alt_az[1].y > 0.)
                     tbuff[2] = '*';         /* daylight */
                  else if( alt_az[1].y > -6. * PI / 180.)
                     tbuff[2] = 'C';         /* civil twilight */
                  else if( alt_az[1].y > -12. * PI / 180.)
                     tbuff[2] = 'N';         /* civil twilight */
                  else if( alt_az[1].y > -18. * PI / 180.)
                     tbuff[2] = 'A';         /* civil twilight */
                  else if( alt_az[2].y > 0.)      /* moon's up */
                     tbuff[2] = (lunar_elong > PI / 2. ? 'M' : 'm');
                  tbuff[3] = '\0';
                  if( visibility_char == ' ')
                     {
                     int galact_conf = galactic_confusion( ra * 15, dec);

                     if( galact_conf > 15)
                        {
                        int band;

                        galact_conf = galact_conf * 100 / 256;       /* scale to 0-99 */
                        for( band = 0; band < 24; band += 8)
                           if( ((rgb >> band) & 0xff) < (long)galact_conf)
                              {
                              rgb &= ~(0xff << band);
                              rgb |= (unsigned long)galact_conf << band;
                              }
                        if( galact_conf > 25)
                           tbuff[1] = (galact_conf > 60 ? 'G' : 'g');
                        }
                     }
                  if( computer_friendly)
                     strcat( buff, tbuff);
                  else
                     snprintf_append( buff, sizeof( buff), "$%06lx%s", rgb, tbuff);
                  if( tbuff[1] == ' ' && tbuff[2] == ' ')
                     tbuff[1] = '-';
                  strcat( alt_buff, tbuff);
                  }
               if( options & OPTION_SKY_BRIGHTNESS)
                  {
                  if( mags_per_arcsec2 > 99.9)
                     mags_per_arcsec2 = 99.99;
                  snprintf( tbuff, sizeof( tbuff), " %5.2f", mags_per_arcsec2);
                  strcat( alt_buff, tbuff);
                  if( mags_per_arcsec2 > 99.9 && !computer_friendly)
                     strcpy( tbuff, " --.--");
                  strcat( buff, tbuff);
                  }

               curr_mag = abs_mag + calc_obs_magnitude(
                          solar_r, r, earth_r, &phase_ang);  /* elem_out.cpp */
               if( fraction_illum != 1.)
                  {
                  if( fraction_illum > 0.)
                     curr_mag -= 2.5 * log10( fraction_illum);
                  else
                     curr_mag = 999.;
                  }
               if( curr_mag > 999.)       /* avoid overflow for objects     */
                  curr_mag = 999.;        /* essentially at zero elongation */
               if( curr_mag > ephemeris_mag_limit)
                  show_this_line = false;

               if( fake_astrometry)
                  {
                  if( abs_mag && curr_mag < 99.)
                     snprintf_append( fake_line, sizeof( fake_line), "%4.1f V", curr_mag);
                  else
                     strcat( fake_line, "      ");
                  snprintf_append( fake_line, sizeof( fake_line), " Synth%.3s",
                                           note_text + 1);
                  }

               if( rho_cos_phi || rho_sin_phi)
                  {
                  exposure_config.sky_brightness = mags_per_arcsec2;
                  if( alt_az[0].y < 0.)
                     exposure_config.airmass = 1e+10;
                  else
                     exposure_config.airmass = air_mass;
                  }

               if( options & OPTION_SNR)
                  {
                  if( exposure_config.airmass > 1e+9)
                     {
                     strcat( buff, (computer_friendly ? " 99999" : " --.--"));
                     strcat( alt_buff, " 99999");
                     }
                  else
                     {
                     const double exp_time = atof( get_environment_ptr( "EXPTIME"));
                     const double snr = snr_from_mag_and_exposure( &exposure_config,
                                    curr_mag, (exp_time ? exp_time : 30.));
                     const char *fmt = (snr > 99. ? " %5.0f" : " %5.2f");

                     snprintf( tbuff, sizeof( tbuff), fmt, snr);
                     strcat( buff, tbuff);
                     strcat( alt_buff, tbuff);
                     }
                  }

               if( options & OPTION_EXPOSURE_TIME)
                  {
                  const double target_snr = atof( get_environment_ptr( "SNR"));
                  double exposure_time;

                  if( exposure_config.airmass > 1e+9)
                     exposure_time = 99999.;
                  else
                     exposure_time = exposure_from_snr_and_mag( &exposure_config,
                                  (target_snr ? target_snr : 4.), curr_mag);
                  if( exposure_time > 99999.)
                     exposure_time = 99999.;
                  snprintf_append( alt_buff, sizeof( alt_buff), " %.1f", exposure_time);
                  if( exposure_time > 99999. && !computer_friendly)
                     strcat( buff, " -----");
                  else
                     snprintf_append( buff, sizeof( buff),
                           (exposure_time < 999. ? " %5.1f" : " %5.0f"), exposure_time);
                  }

               *tbuff = '\0';
               if( options & OPTION_PHASE_ANGLE_OUTPUT)
                  snprintf( tbuff, sizeof( tbuff), " %8.4f", phase_ang * 180. / PI);


               if( options & OPTION_PHASE_ANGLE_BISECTOR)
                  {
                  double pab_vector[3], pab_lon, pab_lat;

                  for( j = 0; j < 3; j++)
                     pab_vector[j] = topo_ecliptic[j] / r
                                          + orbi_after_light_lag[j] / solar_r;
                  vector_to_polar( &pab_lon, &pab_lat, pab_vector);
                  add_lon_lat_to_ephem( tbuff, sizeof( tbuff),
                              pab_lon, pab_lat);
                  }

               if( options & OPTION_HELIO_ECLIPTIC)
                  {
                  double eclip_lon, eclip_lat;

                  vector_to_polar( &eclip_lon, &eclip_lat, orbi_after_light_lag);
                  add_lon_lat_to_ephem( tbuff, sizeof( tbuff),
                              eclip_lon, eclip_lat);
                  }

               if( options & OPTION_TOPO_ECLIPTIC)
                  {
                  double eclip_lon, eclip_lat;

                  vector_to_polar( &eclip_lon, &eclip_lat, topo_ecliptic);
                  add_lon_lat_to_ephem( tbuff, sizeof( tbuff),
                              eclip_lon, eclip_lat);
                  }
               if( options & OPTION_GALACTIC_COORDS)
                  {
                  double galactic_lat, galactic_lon;

                  ra_dec_to_galactic( ra_dec.x, ra_dec.y,
                                                  &galactic_lat, &galactic_lon);
                  if( galactic_lon < 0.)
                     galactic_lon += PI + PI;
                  add_lon_lat_to_ephem( tbuff, sizeof( tbuff),
                              galactic_lon, galactic_lat);
                  }
               if( options & OPTION_GALACTIC_CONFUSION)
                  {
                  const int galact_conf =
                               galactic_confusion( ra * 15, dec) * 100 / 256;

                  snprintf_append( tbuff, sizeof( buff), " %02d", galact_conf);
                  }
               if( options & OPTION_SUN_TARGET_PA)
                  {
                  double sun_target_pa, unused_elongation;

                  calc_dist_and_posn_ang( &ra_dec.x, &sun_ra_dec.x,
                           &unused_elongation, &sun_target_pa);
                  sun_target_pa = PI - sun_target_pa;
                  if( sun_target_pa < 0.)
                     sun_target_pa += PI + PI;
                  snprintf_append( tbuff, sizeof( buff), " %7.3f",
                                  sun_target_pa * 180. / PI);
                  }
               if( options & OPTION_SUN_HELIO_VEL_PA)
                  {
                  double unused_dist, pa, helio_vel[3];
                  DPT vel_pole;

                  memcpy( helio_vel, orbi_after_light_lag + 3, 3 * sizeof( double));
                  ecliptic_to_equatorial( helio_vel);
                  vector_to_polar( &vel_pole.x, &vel_pole.y, helio_vel);
                  calc_dist_and_posn_ang( &ra_dec.x, &vel_pole.x,
                                        &unused_dist, &pa);
                  pa = PI - pa;
                  if( pa < 0.)
                     pa += PI + PI;
                  snprintf_append( tbuff, sizeof( buff), " %7.3f",
                                 pa * 180. / PI);
                  }
               if( options & OPTION_ORBIT_PLANE_ANGLE)
                  {
                  double orbit_norm[3], unused_pa, dist;
                  DPT orbit_pole;

                  vector_cross_product( orbit_norm, orbi, orbi + 3);
                  ecliptic_to_equatorial( orbit_norm);
                  vector_to_polar( &orbit_pole.x, &orbit_pole.y, orbit_norm);
                  calc_dist_and_posn_ang( &ra_dec.x, &orbit_pole.x,
                                        &dist, &unused_pa);
                  snprintf_append( tbuff, sizeof( buff), " %7.3f",
                                    dist * 180. / PI - 90.);
                  }

               strcat( buff, tbuff);
               strcat( alt_buff, tbuff);
               if( abs_mag)           /* don't show a mag if you dunno how bright */
                  {                   /* the object really is! */
                  const bool two_place_mags =
                                   (*get_environment_ptr( "MAG_DIGITS") == '2');

                  if( two_place_mags)
                     {
                     const char *fmt = (curr_mag > 99.8 ? " %5.1f" : " %5.2f");

                     snprintf_append( buff, sizeof( buff), fmt, curr_mag + .005);
                     }
                  else if( fraction_illum == 0.)
                     strcat( buff, " Sha ");
                  else if( curr_mag < 99 && curr_mag > -9.9)
                     snprintf_append( buff, sizeof( buff), " %4.1f", curr_mag + .05);
                  else
                     snprintf_append( buff, sizeof( buff), " %3d ", (int)( curr_mag + .5));
                  if( phase_ang > PI * 2. / 3.)    /* over 120 degrees */
                     if( object_type != OBJECT_TYPE_COMET
                                 && !computer_friendly)
                        {
                        char *endptr = buff + strlen( buff);

                        endptr[-1] = '?';      /* signal doubtful magnitude */
                        if( endptr[-2] == '.')
                           endptr[-2] = '?';
                        }
                  snprintf_append( alt_buff, sizeof( alt_buff), " %6.3f", curr_mag + .0005);
                  }

               if( options & OPTION_LUNAR_ELONGATION)
                  {
                  snprintf_append( buff, sizeof( buff), " %6.1f", dist_moon * 180. / PI);
                  snprintf_append( alt_buff, sizeof( alt_buff), " %10.5f",
                                                             dist_moon * 180. / PI);
                  }

               if( options & OPTION_MOTION_OUTPUT)
                  {
                  MOTION_DETAILS m;
                  char *end_ptr = buff + strlen( buff);
                  char *alt_tptr = alt_buff + strlen( alt_buff);

                  compute_observation_motion_details( &temp_obs, &m);
                  strcat( buff, " ");
                  if( options & OPTION_SEPARATE_MOTIONS)
                     {
                     format_motion( end_ptr + 1, m.ra_motion);
                     format_motion( end_ptr + 9, m.dec_motion);
                     snprintf_append( alt_buff, sizeof( alt_buff),
                                 " %f %f", m.ra_motion, m.dec_motion);
                     }
                  else
                     {
                     format_motion( end_ptr + 1, m.total_motion);
                     snprintf( end_ptr + 9, 7, "%5.1f ",
                                     m.position_angle_of_motion);
                     snprintf_append( alt_buff, sizeof( alt_buff),
                                 " %f %.4f", m.total_motion, m.position_angle_of_motion);
                     }
                  end_ptr[8] = ' ';
                  if( computer_friendly)
                     strcpy( end_ptr, alt_tptr);
                  }
               for( j = 0; j < 3; j++)
                  {
                  bool show_alt, show_az;
                  const double alt = alt_az[j].y * 180. / PI;
                  const double az  = alt_az[j].x * 180. / PI;

                  if( j == 1)
                     {
                     show_alt = (options & OPTION_SUN_ALT);
                     show_az  = (options & OPTION_SUN_AZ);
                     }
                  else if( j == 2)
                     {
                     show_alt = (options & OPTION_MOON_ALT);
                     show_az  = (options & OPTION_MOON_AZ);
                     }
                  else
                     show_alt = show_az = (options & OPTION_ALT_AZ_OUTPUT);
                  alt_tptr = alt_buff + strlen( alt_buff);
                  *tbuff = '\0';
                  if( show_alt)
                     {
                     snprintf( tbuff, sizeof( tbuff), " %c%02d",
                                       (alt > 0. ? '+' : '-'),
                                       (int)( fabs( alt) + .5));
                     snprintf_append( alt_buff, sizeof( alt_buff), " %8.4f", alt);
                     }
                  if( show_az)
                     {
                     snprintf_append( tbuff, sizeof( tbuff), " %03d",
                                       (int)( az + .5));
                     snprintf_append( alt_buff, sizeof( alt_buff), " %8.4f", az);
                     }
                  if( !computer_friendly)
                     strcat( buff, tbuff);
                  else
                     strcat( buff, alt_tptr);
                  }
               if( options & OPTION_RADIAL_VEL_OUTPUT)
                  {
                  char *end_ptr = buff + strlen( buff);
                  const double rvel_in_km_per_sec =
                                           radial_vel * AU_IN_KM / seconds_per_day;

                  snprintf_append( alt_buff, sizeof( alt_buff),
                                      " %11.6f", rvel_in_km_per_sec);
                  format_velocity_in_buff( end_ptr, rvel_in_km_per_sec);
                  }
               if( show_radar_data)
                  {
                  if( alt_az[0].y < 0.)
                     strcpy( tbuff, "  n/a");
                  else
                     {
                     const double radar_albedo = 0.1;
                     double snr = radar_snr_per_day( &rdata, abs_mag,
                                    radar_albedo, r);

                     *tbuff = ' ';
                     show_packed_with_si_prefixes( tbuff + 1, snr);
                     }
                  strcat( alt_buff, tbuff);
                  strcat( buff, tbuff);
                  }
               if( options & OPTION_GROUND_TRACK)
                  {
                  double lat_lon[2], alt_in_meters;
                  const double meters_per_km = 1000.;

                  alt_in_meters = find_lat_lon_alt( utc, geo, planet_no, lat_lon,
                           *get_environment_ptr( "GEOMETRIC_GROUND_TRACK") == '1');
                  snprintf( tbuff, 30, "%9.4f %+08.4f %10.3f",
                        lat_lon[0] * 180. / PI,
                        lat_lon[1] * 180. / PI,
                        alt_in_meters / meters_per_km);
                  tbuff[30] = '\0';
                  strcat( alt_buff, tbuff);
                  strcat( buff, tbuff);
                  }

               if( options & OPTION_SPACE_VEL_OUTPUT)
                  {
                           /* get 'full' velocity; cvt AU/day to km/sec: */
                  const double total_vel =
                             vector3_length( topo_vel) * AU_IN_KM / seconds_per_day;

                  format_velocity_in_buff( tbuff, total_vel);
                  strcat( buff, tbuff);
                  snprintf_append( alt_buff, sizeof( alt_buff),  " %11.6f", total_vel);
                  }
               if( options & OPTION_SUPPRESS_UNOBSERVABLE)
                  if( show_this_line)
                     {
                     show_this_line = (visibility_char == ' ');
                     if( !show_radar_data && alt_az[1].y > 0.)
                        show_this_line = false;  /* optical data must be taken at night */
                     }

               if( fake_astrometry)
                  strcpy( buff, fake_line);
               if( !show_this_line)
                  {
                  if( last_line_shown)
                     strcpy( buff, "................");
                  else
                     *buff = '\0';
                  *alt_buff = '\0';
                  }
               last_line_shown = show_this_line;
               }
            }
         else        /* shouldn't happen */
            strcpy( buff, "DANGER!\n");
         if( !obj_n && (options & OPTION_MOIDS) && show_this_line)
            for( j = 1; j <= 8; j++)
               {
               double moid;
               ELEMENTS planet_elem, elem;
               const double GAUSS_K = .01720209895;  /* Gauss' grav const */
               const double SOLAR_GM = (GAUSS_K * GAUSS_K);

               elem.central_obj = 0;
               elem.gm = SOLAR_GM;
               elem.epoch = curr_jd;
               calc_classical_elements( &elem, orbi, curr_jd, 1);
               setup_planet_elem( &planet_elem, j, (curr_jd - J2000) / 36525.);
               moid = find_moid_full( &planet_elem, &elem, NULL);
               snprintf_append( buff, sizeof( buff), "%8.4f", moid);
               }
         if( !obj_n && show_this_line)
            n_lines_shown++;
         }
      if( *buff)
         fprintf( ofile, "%s\n", buff);
      if( computer_friendly_ofile && *alt_buff)
         fprintf( computer_friendly_ofile, "%s\n", alt_buff);
      prev_ephem_t = ephemeris_t;
      }
   free( orbits_at_epoch);
   if( ephem_type == OPTION_OBSERVABLES && !n_lines_shown)
      fprintf( ofile, "No ephemeris output.  Object was too faint,  or in daylight,\n"
                   "or below horizon for the specified times.  Check ephem options.\n");
   fclose( ofile);
   if( header && computer_friendly_ofile)
      {
      extern const char *combine_all_observations;

      ofile = open_json_file( buff, "JSON_EPHEM_NAME", "ephemeri.json",
                              obs->packed_id, "w+");
      fprintf( ofile, "{\n  \"ephemeris\":\n  {\n");
      fprintf( ofile, "    \"obscode\": \"%.3s\",\n", note_text + 1);
      if( combine_all_observations && *combine_all_observations)
         strcpy( buff, combine_all_observations);
      else
         strcpy( buff, obs->packed_id);
      text_search_and_replace( buff, " ", "");
      fprintf( ofile, "    \"packed\": \"%s\",\n", buff);
      fprintf( ofile, "    \"count\": %d,\n", n_lines_shown);
      fprintf( ofile, "    \"start\": %f,\n", real_jd_start);
      fprintf( ofile, "    \"step\": %.9f,\n", step);
      fprintf( ofile, "    \"n_steps\": %d,\n", n_steps);
      fprintf( ofile, "    \"start iso\": \"%s\",\n", iso_time( buff, real_jd_start, 3));
      fprintf( ofile, "    \"entries\":\n");
      fprintf( ofile, "    {\n");
      fseek( computer_friendly_ofile, 0L, SEEK_SET);
      create_json_ephemeris( ofile, computer_friendly_ofile, header, real_jd_start, step);
      fprintf( ofile, "    }\n");
      free( header);
      fprintf( ofile, "  }\n}\n");
      combine_json_elems_and_ephems( obs->packed_id, ofile);
      fclose( ofile);
      }
   if( computer_friendly_ofile)
      fclose( computer_friendly_ofile);
   if( list_of_ephem_times)
      {
      free( list_of_ephem_times);
      list_of_ephem_times = NULL;
      }
   return( 0);
}
/* The above ephemeris code may insert a six-character hexadecimal
color,  computed using sky background brightness data,  prefaced by
'$'.  This should be removed before the line is shown to a user.
The RGB value is returned,  and may (or may not) be made use of. */

int remove_rgb_code( char *buff)
{
   unsigned rval = (unsigned)-1;

   buff = strchr( buff, '$');
   if( buff && sscanf( buff + 1, "%06x", &rval) == 1)
      memmove( buff, buff + 7, strlen( buff + 6));
   return( (int)rval);
}

int save_ephemeris_file( const char *filename)
{
   int rval = -1;
   char buff[400];
   FILE *ifile = fopen_ext( get_file_name( buff, ephemeris_filename),
               is_default_ephem ? "tcr" : "r");

   if( ifile)
      {
      FILE *ofile = fopen_ext( filename, "fw");

      while( fgets( buff, sizeof( buff), ifile))
         {
         remove_rgb_code( buff);
         fputs( buff, ofile);
         }
      fclose( ofile);
      rval = 0;
      }
   return( rval);
}

/* "is_topocentric_mpc_code( )" is taken to mean "can you compute alt/az
ephems and/or visibility info from this station",  and is used in
figuring out which options are available for an ephemeris. */

bool is_topocentric_mpc_code( const char *mpc_code)
{
   char buff[100];
   double rho_cos_phi, rho_sin_phi, lon;
   const int planet_idx = get_observer_data( mpc_code,
                    buff, &lon, &rho_cos_phi, &rho_sin_phi);

   return( planet_idx >= 0 && (rho_cos_phi != 0. || rho_sin_phi != 0.));
}

static void get_scope_params( const char *mpc_code, expcalc_config_t *c)
{
   FILE *ifile = fopen_ext( "scope.json", "fclrb");

   find_expcalc_config_from_mpc_code( mpc_code, ifile, c);
   fclose( ifile);
}

int ephemeris_in_a_file_from_mpc_code( const char *filename,
         const double *orbit,
         OBSERVE *obs, const int n_obs,
         const double epoch_jd, const double jd_start, const char *stepsize,
         const int n_steps, const char *mpc_code,
         ephem_option_t options, const unsigned n_objects)
{
   double rho_cos_phi, rho_sin_phi, lon;
   char note_text[100], buff[100];
   int real_number_of_steps, rval;
   const int planet_no = get_observer_data( mpc_code, buff, &lon,
                                           &rho_cos_phi, &rho_sin_phi);

   assert( strlen( mpc_code) >= 3);
   strlcpy_err( ephem_mpc_code, mpc_code, sizeof( ephem_mpc_code));
   snprintf( note_text, sizeof( note_text),
                    "(%s) %s", mpc_code, mpc_station_name( buff));
   get_object_name( buff, obs->packed_id);
   snprintf_append( note_text, sizeof( note_text), ": %s", buff);
   get_scope_params( mpc_code, &exposure_config);
   if( !strcmp( stepsize, "Obs"))
      real_number_of_steps = n_obs;
   else if( *stepsize == 't')
      real_number_of_steps = get_ephem_times_from_file( stepsize + 1);
   else
      real_number_of_steps = n_steps;
   rval = ephemeris_in_a_file( filename, orbit, obs, n_obs, planet_no,
               epoch_jd, jd_start, stepsize, lon, rho_cos_phi, rho_sin_phi,
               real_number_of_steps,
               note_text, options, n_objects);
   free_expcalc_config_t( &exposure_config);
   return( rval);
}

static int64_t ten_to_the_nth( int n)
{
   int64_t rval = 1;

   while( n--)
      rval *= 10;
   return( rval);
}

/* See comments for get_ra_dec() in mpc_fmt.cpp for info on the meaning  */
/* of 'precision' in this function.                                      */

static void output_angle_to_buff( char *obuff, double angle, int precision)
{
   int n_digits_to_show = 0;
   int64_t power_mul, fraction;
   size_t i, full_len = 12;

   if( (precision >= 100 && precision <= 116) /* decimal quantity, dd.dd... */
        || ( precision >= 200 && precision <= 215)) /* decimal ddd.dd... */
      {
      const int two_digits = (precision <= 200);

      n_digits_to_show = precision % 100;
      power_mul = ten_to_the_nth( n_digits_to_show);
      fraction = (int64_t)( angle * (two_digits ? 1. : 15.) * (double)power_mul + .5);

      snprintf( obuff, 4, (two_digits ? "%02u" : "%03u"),
                  (int)( fraction / power_mul));
      fraction %= power_mul;
      }
   else
      switch( precision)
         {
         case -1:       /* hh mm,  integer minutes */
         case -2:       /* hh mm.m,  tenths of minutes */
         case -3:       /* hh mm.mm,  hundredths of minutes */
         case -4:       /* hh mm.mmm,  milliminutes */
         case -5:       /* hh mm.mmmm */
         case -6:       /* hh mm.mmmmm */
         case -7:       /* hh mm.mmmmmm */
            {
            n_digits_to_show = -1 - precision;
            power_mul = ten_to_the_nth( n_digits_to_show);
            fraction = (int64_t)( angle * 60. * (double)power_mul + .5);
            snprintf( obuff, 6, "%02u %02u", (unsigned)( fraction / ((int64_t)60 * power_mul)),
                                         (unsigned)( fraction / power_mul) % 60);
            fraction %= power_mul;
            }
            break;
         case 0:        /* hh mm ss,  integer seconds */
         case 1:        /* hh mm ss.s,  tenths of seconds */
         case 2:        /* hh mm ss.ss,  hundredths of seconds */
         case 3:        /* hh mm ss.sss,  thousands of seconds */
         case 4: case 5: case 6:    /* possible extra digits in ephems */
         case 7: case 8: case 9:
         case 307:      /* hhmmsss,  all packed together:  tenths */
         case 308:      /* hhmmssss,  all packed together: hundredths */
         case 309:      /* milliseconds (or milliarcseconds) */
         case 310:      /* formats 307-312 are the 'super-precise' formats */
         case 311:
         case 312:      /* microseconds (or microarcseconds) */
         case 400:      /* (RA) ddd mm ss,  integer arcseconds */
         case 401:      /* (RA) ddd mm ss.s,  tenths of arcseconds */
         case 402:      /* (RA) ddd mm ss.ss,  centiarcsec */
            {
            const char *format = "%02u %02u %02u";

            if( precision >= 400)
               {
               format = "%03u %02u %02u";
               precision -= 400;
               angle *= 15.;
               }
            n_digits_to_show = precision % 306;
            power_mul = ten_to_the_nth( n_digits_to_show);
            fraction = (int64_t)( angle * 3600. * (double)power_mul + .5);
            snprintf( obuff, 10, format,
                     (unsigned)( fraction / ((int64_t)3600 * power_mul)),
                     (unsigned)( fraction / ((int64_t)60 * power_mul)) % 60,
                     (unsigned)( fraction / power_mul) % 60);
            fraction %= power_mul;
            if( precision > 306)          /* remove spaces: */
               text_search_and_replace( obuff, " ", "");
            }
            break;
         default:                  /* try to show the angle,  but indicate */
            fraction = 0;   /* not really necessary;  evades nuisance GCC warning */
            if( angle > -1000. && angle < 1000.)   /* the format is weird  */
               snprintf( obuff, 10, "?%.5f", angle);
            else
               strcpy( obuff, "?");
            break;
         }
                  /* Formats not used in astrometry -- they don't fit the */
                  /* field size for punched-card data,  but are used in ephems */
   if( precision >= 4 && precision <= 11)    /* 'overlong' dd mm ss.ssss... */
      full_len += (size_t)( precision - 3);
   if( precision >= 209 && precision <= 215) /* 'overlong' ddd mm ss.sss... */
      full_len += (size_t)( precision - 208);
   if( precision >= 110 && precision <= 116) /* 'overlong' dd mm ss.sss... */
      full_len += (size_t)( precision - 109);
   if( n_digits_to_show)
      {
      char format[7];

      if( precision < 307 || precision > 312)   /* omit decimal point for */
         strcat( obuff, ".");                    /* super-precise formats */
      snprintf( format, sizeof( format), "%%0%dd", n_digits_to_show);
      snprintf_append( obuff, full_len + 1, format, fraction);
      }
   for( i = strlen( obuff); i < full_len; i++)
      obuff[i] = ' ';
   obuff[full_len] = '\0';
}

static void output_signed_angle_to_buff( char *obuff, const double angle,
                               const int precision)
{
   *obuff++ = (angle < 0. ? '-' : '+');
   output_angle_to_buff( obuff, fabs( angle), precision);
}

/* 'put_residual_into_text( )' expresses a residual,  from 0 to 180 degrees, */
/* such that the text starts with a space,  followed by four characters,   */
/* and a sign:  six bytes in all.  This can be in forms such as:           */
/*                                                                         */
/*  Err!+    (for values above 999 degrees)                                */
/*  179d-    (for values above 999 arcminutes = 16.6 degrees)              */
/*  314'+    (for values above 9999 arcsec but below 999 arcminutes)       */
/*  7821-    (for values below 9999 arcsec but above 99 arcsec)            */
/*  12.3+    (for values above one arcsec but below 99 arcsec)             */
/*  2.71-    (1" < value < 10",  if precise_residuals == 1)                */
/*   .87-    (for values under an arcsecond, usually)                      */
/*  .871-    (for values under an arcsecond, if precise_residuals == 1)    */
/*                                                                         */
/*     If precise_residuals == 2 ("super-precise" residuals) and the value */
/* is below a milliarcsecond,  some special formats such as '3.2u' or      */
/* '451u' ('u'=microarcseconds) or '4.5n' or '231n' (nanoarcsec),  or      */
/* even picoarcseconds,  are used.  (This because I was investigating      */
/* some roundoff/precision issues.  There would normally be no practical   */
/* reason to go to the picoarcsecond level of precision!)                  */

static void put_residual_into_text( char *text, const double resid,
                                 const int resid_format)
{
   double zval = fabs( resid);
   const int precise = resid_format
                & (RESIDUAL_FORMAT_OVERPRECISE | RESIDUAL_FORMAT_PRECISE);

   if( resid_format & RESIDUAL_FORMAT_COMPUTER_FRIENDLY)
      {                   /* resids in arcseconds at all times,  with */
      if( resid > -9.9999 && resid < 9.9999)
         snprintf( text, 11, " %+8.6f", resid);    /* some added precision */
      return;
      }
   if( zval > 999. * 3600.)      /* >999 degrees: error must have occurred */
      strcpy( text, " Err!");
   else if( zval > 59940.0)             /* >999': show integer degrees */
      snprintf( text, 6, "%4.0fd", zval / 3600.);
   else if( zval > 9999.9)              /* 999' > x > 9999": show ###' arcmin */
      snprintf( text, 6, "%4.0f'", zval / 60.);
   else if( zval > 99.9)
      snprintf( text, 6, "%5.0f", zval);
   else if( zval > .99 && zval < 9.99 && precise)
      snprintf( text, 6, "%5.2f", zval);
   else if( zval > .99)
      snprintf( text, 6, "%5.1f", zval);
   else if( (resid_format & RESIDUAL_FORMAT_OVERPRECISE) && zval < .00999)
      {          /* 'high-precision' residuals */
      unsigned i;
      const char *lower_si_prefixes = " munpfazy ";

      for( i = 0; zval < 0.99 && i < 9; i++)
         zval *= 1000.;
      snprintf( text, 6, (zval < 9.9 ? "%4.1f%c" : "%4.0f%c"),
                     zval, lower_si_prefixes[i]);
      }
   else
      {
      snprintf( text, 6, (precise ? "%5.3f" : "%5.2f"), zval);
      text[precise ? 0 : 1] = ' ';
      }
   if( !atof( text))
      text[5] = ' ';
   else
      text[5] = (resid > 0. ? '+' : '-');
   text[6] = '\0';
}

static void show_dd_hh_mm_ss_point_sss( char *text,
                          const double day, int precision)
{
   const int64_t milliseconds = (int64_t)( day * seconds_per_day * 1000. + .1);
   const int64_t ms_per_minute = (int64_t)( 60 * 1000);
   const int64_t ms_per_hour = 60 * ms_per_minute;
   const int64_t ms_per_day = 24 * ms_per_hour;

   snprintf( text, 15, "%02u %02u:%02u:%02u%03u",
            (unsigned)( milliseconds / ms_per_day),
            (unsigned)( milliseconds / ms_per_hour) % 24,
            (unsigned)( milliseconds / ms_per_minute) % 60,
            (unsigned)( milliseconds / 1000) % 60,
            (unsigned)( milliseconds % 1000));
   while( precision < 3)
      text[11 + precision++] = ' ';
}

static void put_mag_resid( char *output_text, const double obs_mag,
                           const double computed_mag, const char mag_band)
{

   INTENTIONALLY_UNUSED_PARAMETER( mag_band);

   if( obs_mag < BLANK_MAG && computed_mag)
      snprintf( output_text, 8, "%6.2f ",
               obs_mag - computed_mag);
//             obs_mag - computed_mag - mag_band_shift( mag_band);
   else
      strcpy( output_text, "------ ");
}

/* Output is a number in three bytes,  such as '3.1', ' 31', '314', ' 3k',
'31k', '.3M',  etc. */

static void show_number_in_three_bytes( char *buff, double ival)
{
   if( ival < 999)
      snprintf( buff, 4, (ival < 9.9 ? "%3.1f" : "%3.0f"), ival);
   else
      {
      int i, digit;

      ival /= 1000.;
      for( i = 0; ival > 99. && si_prefixes[i]; i++)
         ival /= 1000.;
      digit = (int)( ival * 10);
      if( !si_prefixes[i])
         strcpy( buff, "!!!");
      else
         {
         if( digit < 10)
            {
            *buff++ = '.';
            *buff++ = '0' + digit;
            }
         else
            {
            digit /= 10;
            *buff++ = '0' + digit / 10;
            *buff++ = '0' + digit % 10;
            }
         *buff++ = si_prefixes[i];
         *buff++ = '\0';
         }
      }
}

static void show_resid_in_sigmas( char *buff, const double sigmas)
{
   *buff++ = ' ';
   *buff++ = (sigmas > 0. ? '+' : '-');
   show_number_in_three_bytes( buff, fabs( sigmas));
   buff[3] = ' ';
   buff[4] = '\0';
}

/* format_observation( ) takes an observation and produces text for it,
   suitable for display on a console (findorb) or in a Windoze scroll
   box (FIND_ORB),  or for writing to a file.  */

void format_observation( const OBSERVE FAR *obs, char *text,
                                        const int resid_format)
{
   double angle;
   char xresid[30], yresid[30];
   int i;
   extern int apply_debiasing;
   const int base_format = (resid_format & 3);
   const int base_time_format = obs->time_precision / 10;
   const int n_time_digits = obs->time_precision % 10;
   const int four_digit_years =
                    (resid_format & RESIDUAL_FORMAT_FOUR_DIGIT_YEARS);
   int month;
   long year;
   double day, utc;
   char *original_text_ptr = text;
   MOTION_DETAILS m;

   utc = obs->jd - td_minus_utc( obs->jd) / seconds_per_day;
   day = decimal_day_to_dmy( utc, &year, &month, CALENDAR_JULIAN_GREGORIAN);

   if( base_format != RESIDUAL_FORMAT_SHORT)
      {
      switch( base_time_format)
         {
         case 2:        /* CYYMMDD:HHMMSSsss:  formats 20-23 */
         case 3:        /* CYYMMDD.ddddddddd:  formats 30-39 */
            snprintf( text, 6, "%c%02u%02u",     /* show century letter, 2digit yr, mo */
                     (char)( 'A' + year / 100 - 10), (unsigned)( year % 100),
                     (unsigned)month);
            if( base_time_format == 2)
               {
               show_dd_hh_mm_ss_point_sss( text + 5, day, n_time_digits);
               text[7] = ':';
               text[10] = text[11];     /* Turn dd hh:mm:ss into dd:hhmmss. This */
               text[11] = text[12];     /* corresponds to the somewhat weird way */
               text[12] = text[14];     /* in which you can specify precise times */
               text[13] = text[15];     /* within the limits of the MPC's 80-byte */
               text[14] = text[16];     /* format (a Find_Orb extension).         */
               text[15] = text[17];
               text[16] = text[18];
               text[17] = text[19];
               text[18] = '\0';
               }
            else
               snprintf( text + 5, 13, "%012.9f", day);
            break;
         case 4:        /* M056336.641592653: MJD formats 40-49 */
         case 1:        /* 2456336.641592653: JD formats 10-49 */
            {
            snprintf( text, 18, "%017.9f",
                           utc - (base_time_format == 4 ? 2400000.5 : 0.));
            if( base_time_format == 4)
               *text = 'M';
            }
            break;
         case 0:        /* "Standard" MPC YYYY MM DD.dddddd : formats 0-6 */
            {
            assert( n_time_digits <= 6);
            if( four_digit_years)
               snprintf( text, 9, "%04ld\t%02d\t", year, month);
            else
               snprintf( text, 7, "%02d\t%02d\t", abs( (int)year % 100), month);
            text += strlen( text);
            if( resid_format & RESIDUAL_FORMAT_HMS)
               show_dd_hh_mm_ss_point_sss( text, day, 0);
            else
               {
               const char *date_format_text[7] = { "%02.0f       ",
                                                   "%04.1f     ",
                                                   "%05.2f    ",
                                                   "%06.3f   ",
                                                   "%07.4f  ",
                                                   "%08.5f ",
                                                   "%09.6f" };

               snprintf( text, 10, date_format_text[n_time_digits], day);
               }
            }
            break;
         default:
            assert( true);
            break;
         }
      if( base_time_format == 3 || base_time_format == 1 || base_time_format == 4)
         for( i = n_time_digits + 8; i < 17; i++)     /* clear excess digits */
            text[i] = ' ';
      snprintf_append( text, 40, "\t%c\t%s\t",
                   (obs->is_included ? ' ' : 'X'), obs->mpc_code);
      angle = obs->ra * 12. / PI;
      if( apply_debiasing)
         angle += obs->ra_bias / (3600. * 15. * cos( obs->dec));
      angle = fmod( angle, 24.);
      if( angle < 0.) angle += 24.;
      output_angle_to_buff( text + strlen( text), angle, obs->ra_precision);
      strcat( text, (base_format == RESIDUAL_FORMAT_FULL_WITH_TABS) ?
                              "\t" : "\t ");
      }
   else        /* 'short' MPC format: */
      {
      if( four_digit_years)
         *text++ = int_to_mutant_hex_char( year / 100);
      snprintf( text, 11, "%02u%02u%02u %s", (unsigned)abs( (int)year % 100),
                  (unsigned)month, (unsigned)day, obs->mpc_code);
      }
   text += strlen( text);

   compute_observation_motion_details( obs, &m);        /* mpc_obs.cpp */
   if( obs->note2 == 'R')
      {
      RADAR_INFO rinfo;

      compute_radar_info( obs, &rinfo);
      if( !rinfo.rtt_obs)
         strcpy( xresid, " ---- ");
      else if( resid_format & RESIDUAL_FORMAT_TIME_RESIDS)
         show_resid_in_sigmas( xresid,
                   (rinfo.rtt_obs - rinfo.rtt_comp) / rinfo.rtt_sigma);
      else
         {
         const double time_resid_in_microseconds =
                   (rinfo.rtt_obs - rinfo.rtt_comp) * 1e+6;

         if( fabs( time_resid_in_microseconds) < 999.)
            {
            snprintf( xresid, sizeof( xresid), "%+05d ", (int)( time_resid_in_microseconds * 10.));
            xresid[5] = xresid[0];
            xresid[0] = ' ';
            }
         else
            strcpy( xresid, " HUGE ");
         }
      if( !rinfo.doppler_obs)
         strcpy( yresid, " ---- ");
      else
         {
         const double resid_in_hz = rinfo.doppler_obs - rinfo.doppler_comp;

         if( resid_format & RESIDUAL_FORMAT_TIME_RESIDS)
            show_resid_in_sigmas( yresid, resid_in_hz / rinfo.doppler_sigma);
         else if( fabs( resid_in_hz) < 999.)
            {
            snprintf( yresid, sizeof( yresid), "%+05d ", (int)( resid_in_hz * 10.));
            yresid[5] = yresid[0];
            yresid[0] = ' ';
            }
         else
            strcpy( yresid, " HUGE ");
         }
      }
   else if( resid_format & RESIDUAL_FORMAT_TIME_RESIDS)
      {
      const double abs_time_resid = fabs( m.time_residual);
      const char sign = (m.time_residual < 0. ? '-' : '+');

      if( abs_time_resid < .00094)               /* show as " -.4ms " */
         snprintf( xresid, sizeof( xresid), " %c.%01dms", sign,
                     (int)( abs_time_resid * 10000. + .5));
      else if( abs_time_resid < .099)            /* show as " -47ms " */
         snprintf( xresid, sizeof( xresid), " %c%02dms", sign,
                     (int)( abs_time_resid * 1000. + .5));
      else if( abs_time_resid < .994)            /* show as " +.31s " */
         snprintf( xresid, sizeof( xresid), " %c.%02ds", sign,
                     (int)( abs_time_resid * 100. + .5));
      else if( abs_time_resid < 9.9)             /* show as " -4.7s " */
         snprintf( xresid, sizeof( xresid), " %+4.1fs", m.time_residual);
      else if( abs_time_resid < 999.)            /* show as " -217s " */
         snprintf( xresid, sizeof( xresid), " %c%03ds", sign,
                     (int)( abs_time_resid + .5));
      else if( abs_time_resid / 60. < 999.)      /* show as " +133m " */
         snprintf( xresid, sizeof( xresid), " %c%03dm", sign,
                     (int)( abs_time_resid / 60. + .5));
      else if( abs_time_resid / 3600. < 9999.)   /* show as " +027h " */
         snprintf( xresid, sizeof( xresid), " %c%03dh", sign,
                     (int)( abs_time_resid / 3600. + .5));
      else                                   /* Give up after 1000 hours; */
         strcpy( xresid, " !!!! ");          /* show "it's a long time"   */
      put_residual_into_text( yresid, m.cross_residual, resid_format);
      }
   else
      {
      put_residual_into_text( xresid, m.xresid, resid_format);
      put_residual_into_text( yresid, m.yresid, resid_format);
      if( !strcmp( obs->mpc_code, "258"))
         {
         double resid1, resid2;

         get_residual_data( obs, &resid1, &resid2);
         show_resid_in_sigmas( xresid, resid1);
         show_resid_in_sigmas( yresid, resid2);
         }
      }
   if( base_format != RESIDUAL_FORMAT_SHORT)
      {
      const char *tab_separator =
             ((base_format == RESIDUAL_FORMAT_FULL_WITH_TABS) ? "\t" : "");

      angle = obs->dec * 180. / PI;
      if( apply_debiasing)
         angle += obs->dec_bias / 3600.;
      if( angle < 0.)
         {
         angle = -angle;
         *text++ = '-';
         if( angle < -99.)       /* error prevention */
            angle = -99.;
         }
      else
         {
         *text++ = '+';
         if( angle > 99.)       /* error prevention */
            angle = 99.;
         }
      output_angle_to_buff( text, angle, obs->dec_precision);

      snprintf_append( text, 50, "\t%s%s%s\t", xresid, tab_separator, yresid);
      format_dist_in_buff( xresid, obs->r);
      if( resid_format & RESIDUAL_FORMAT_MAG_RESIDS)
         put_mag_resid( yresid, obs->obs_mag, obs->computed_mag, obs->mag_band);
      else
         format_dist_in_buff( yresid, obs->solar_r);
      snprintf_append( text, 50, "%s%s%s", xresid, tab_separator, yresid);
      }
   else        /* 'short' MPC format */
      {
      if( resid_format & RESIDUAL_FORMAT_MAG_RESIDS)
         {
         put_mag_resid( yresid, obs->obs_mag, obs->computed_mag, obs->mag_band);
         put_residual_into_text( xresid, sqrt( m.xresid * m.xresid
                                             + m.yresid * m.yresid), resid_format);
         xresid[5] = ' ';        /* replace the '+' with a ' ' */
         }
      memcpy( text, xresid, 6);
      memcpy( text + 6, yresid, 6);
      text[0] = (obs->is_included ? ' ' : '(');
      text[12] = (obs->is_included ? ' ' : ')');
      text[13] = '\0';
      }
                       /* for all other formats, replace tabs w/spaces: */
   if( base_format != RESIDUAL_FORMAT_FULL_WITH_TABS)
      for( i = 0; original_text_ptr[i]; i++)
         if( original_text_ptr[i] == '\t')
            original_text_ptr[i] = ' ';

   if( (resid_format & RESIDUAL_FORMAT_EXTRA)
               && base_format != RESIDUAL_FORMAT_SHORT)
      {
      char tbuff[50];
      int tformat = (resid_format | RESIDUAL_FORMAT_SHORT)
            ^ RESIDUAL_FORMAT_TIME_RESIDS ^ RESIDUAL_FORMAT_EXTRA
            ^ RESIDUAL_FORMAT_FOUR_DIGIT_YEARS;

      format_observation( obs, tbuff, tformat);
      strcat( text, tbuff + 11);
      tformat &= ~RESIDUAL_FORMAT_TIME_RESIDS;
      tformat |= RESIDUAL_FORMAT_MAG_RESIDS;
      format_observation( obs, tbuff, tformat);
      strcat( text, tbuff + 11);
      }
}


/* The MPC report format goes to only six decimal places in time,
a "microday".  If the reported time is more precise than that -- as can
happen with video observations -- a workaround is to make use of the object
motion data to adjust the position by up to half a microday.  We only do
this if the time given is more than a nanoday away from an integer
microday.  That simply avoids processing in the (overwhelmingly likely)
case that the data falls exactly on a microday:  i.e.,  it was given
to six or fewer decimal places.

   NOTE:  this should now be obsolete,  at least for Find_Orb,  because
several workarounds for getting up to nine decimal places -- "nanoday"
precision,  or 86.4 microseconds -- are now available.  Hence the call
to this being commented out.  I've not removed the code,  because it's
possible we'll need to reformat astrometry in a manner that'll make MPC
reasonably happy.   */

#ifdef CURRENTLY_UNUSED_POSSIBLY_OBSOLETE
static inline void set_obs_to_microday( OBSERVE FAR *obs)
{
   const double utc = obs->jd - td_minus_utc( obs->jd) / seconds_per_day;
   double delta_jd = utc - floor( utc);

   delta_jd = 1e+6 * delta_jd + .5;
   delta_jd = (delta_jd - floor( delta_jd) - .5) * 1e-6;
   if( delta_jd > 1e-9 || delta_jd < -1e-9)
      {
      MOTION_DETAILS m;
      const double cvt_motions_to_radians_per_day =
                  (PI / 180.) * 24. / 60.;

      compute_observation_motion_details( obs, &m);
      obs->jd -= delta_jd;
                  /* motions are given in '/hr, a.k.a. "/min: */
      obs->ra -= m.ra_motion * delta_jd * cvt_motions_to_radians_per_day
                        / cos( obs->dec);
      obs->dec -= m.dec_motion * delta_jd * cvt_motions_to_radians_per_day;
      }
}
#endif

static void put_sigma( char *buff, const double val)
{
   char tbuff[15];
   const char *format = (val > 10. && val < 999. ? "%3.0f" : "%.2f");

   snprintf( tbuff, sizeof( tbuff), format, val);
   if( *tbuff == '0')    /* skip leading zero */
      memcpy( buff, tbuff + 1, 3);
   else
      memcpy( buff, tbuff, 3);
}


int sigmas_in_columns_57_to_65 = 0;
bool force_traditional_format = false;
      /* The above causes 'extended' time and position formats to
      be shown in punched-card format.  This can help if you want
      to create data for use in (most) software that doesn't support
      Find_Orb's extensions to the MPC format. */

void recreate_observation_line( char *obuff, const OBSERVE FAR *obs)
{
   char buff[100];
   int mag_digits_to_erase = 0;
   OBSERVE tobs = *obs;

   if( obs->note2 == 'R')     /* for radar obs,  we simply store the */
      {                       /* original observation line           */
      strcpy( obuff, obs->second_line + 81);
      return;
      }
// set_obs_to_microday( &tobs);
   if( force_traditional_format)
      {
      if( tobs.ra_precision > 3)
         tobs.ra_precision = 3;
      if( tobs.dec_precision > 2)
         tobs.dec_precision = 2;
      if( tobs.time_precision > 6)
         tobs.time_precision = 6;
      }
   format_observation( &tobs, buff, 4);
   memcpy( obuff, obs->packed_id, 12);
   obuff[12] = obs->discovery_asterisk;
   obuff[13] = obs->note1;
   obuff[14] = obs->note2;
   memcpy( obuff + 15, buff, 17);      /* date/time */
   memcpy( obuff + 32, buff + 24, 12);      /* RA */
   memcpy( obuff + 44, buff + 38, 13);      /* dec */
   snprintf( obuff + 57, 24, "%13.2f%c%c%s%s", obs->obs_mag,
              obs->mag_band, obs->astrometric_net_code, obs->reference, obs->mpc_code);
   if( obs->obs_mag == BLANK_MAG)        /* no mag given;  clean out that value */
      mag_digits_to_erase = 5;
   else
      mag_digits_to_erase = 2 - obs->mag_precision;
   memset( obuff + 70 - mag_digits_to_erase, ' ', mag_digits_to_erase);
   memcpy( obuff + 56, obs->columns_57_to_65, 9);
   if( sigmas_in_columns_57_to_65 &&
               !memcmp( obuff + 56, "         ", 9))
      {
      if( obs->posn_sigma_1 == obs->posn_sigma_2)
         put_sigma( obuff + 59, obs->posn_sigma_1);
      else
         {
         double multiplier = 1.;

         if( obs->posn_sigma_1 < .9 && obs->posn_sigma_2 < .9)
            {
            multiplier = 1000.;    /* show small sigmas as milliarcsec */
            obuff[64] = 'm';
            }
         put_sigma( obuff + 57, obs->posn_sigma_1 * multiplier);
         put_sigma( obuff + 61, obs->posn_sigma_2 * multiplier);
         }
      }
   if( !obs->is_included)
      obuff[64] = 'x';
   if( obs->flags & OBS_DONT_USE)
      obuff[64] = '!';
}

#ifdef NOT_QUITE_READY_YET
void recreate_second_observation_line( char *buff, const OBSERVE FAR *obs)
{
   int i;
   double vect[3];

   buff[32] = '0' + obs->satellite_obs;
   for( i = 0; i < 3; i++)
      vect[i] = obs->obs_posn[j] - ?; (gotta get earths loc somewhere...)
   ecliptic_to_equatorial( vect);
   for( i = 0; i < 3; i++)
      snprintf( buff + 33 + i * 12, 13, "%12.8f", vect[i]);
   buff[69] = ' ';
}
#endif

int process_count;

/* When running on multiple cores,  we need to keep the processes running
on each core from overwriting one another's files.  Thus,  a file such as
'elements.txt' can retain that name in the single-core case.  A second
core should use 'eleme1.txt',  a third 'eleme2.txt',  and so on.  */

char *get_file_name( char *filename, const char *template_file_name)
{
   if( !process_count)
      strcpy( filename, template_file_name);
   else
      {
      const char *tptr = strchr( template_file_name, '.');
      size_t count = tptr - template_file_name;

      assert( tptr);
      assert( process_count < 1000);
      if( count > 5)
         count = 5;
      memcpy( filename, template_file_name, count);
      snprintf( filename + count, 30, "%d%s", process_count, tptr);
      }
   return( filename);
}

void create_obs_file( const OBSERVE FAR *obs, int n_obs, const int append)
{
   char filename[81], curr_sigma_text[81];
   FILE *ofile = fopen_ext( get_file_name( filename, observe_filename),
                           append ? "tfcab" : "tfcwb");

   *curr_sigma_text = '\0';
   while( n_obs--)
      {
      char obuff[81];

      snprintf( obuff, sizeof( obuff), "COM Posn sigma %g", obs->posn_sigma_1);
      if( obs->posn_sigma_2 != obs->posn_sigma_1)  /* elliptical sigma */
         {
         snprintf_append( obuff, sizeof( obuff), " %g", obs->posn_sigma_2);
         if( obs->posn_sigma_theta)    /* uncertainty is a tilted ellipse */
            snprintf_append( obuff, sizeof( obuff), " tilt %.1f",
                        obs->posn_sigma_theta * 180. / PI);
         }
      if( obs->note2 != 'R' && strcmp( curr_sigma_text, obuff))
         {
         fprintf( ofile, "%s\n", obuff);
         strcpy( curr_sigma_text, obuff);
         }
      recreate_observation_line( obuff, obs);
      fprintf( ofile, "%s\n", obuff);
      if( obs->second_line)
         fprintf( ofile, "%s\n", obs->second_line);
      obs++;
      }
   fclose( ofile);
}

static void add_final_period( char *buff)
{
   if( *buff && buff[strlen( buff) - 1] != '.')
      strcat( buff, ".");
}

static void tack_on_names( char *list, const char *names)
{
   while( *names)
      {
      int i, len, already_in_list = 0;

      while( *names == ' ')
         names++;
      for( len = 0; names[len] && names[len] != ','; len++)
         ;
      for( i = 0; list[i]; i++)
         if( !i || (i > 1 && list[i - 2] == ','))
            if( !memcmp( list + i, names, len))
               if( list[i + len] == ',' || !list[i + len])
                  already_in_list = 1;
      if( !already_in_list)
         {
         char *lptr;

         if( *list)
            strcat( list, ", ");
         lptr = list + strlen( list);
         memcpy( lptr, names, len);
         lptr[len] = '\0';
         }
      names += len;
      if( *names == ',')
         names++;
      }
}

static bool got_obs_in_range( const OBSERVE *obs, int n_obs,
               const double jd_start, const double jd_end)
{
   bool rval = false;

   while( !rval && n_obs)
      {
      rval = (obs->jd > jd_start && obs->jd < jd_end);
      n_obs--;
      obs++;
      }
   return( rval);
}

/* For getting default observer/telescope details from 'details.txt' and/or
'scopes.txt',  we look through those files for the MPC observatory code
in question,  then look for OBS/MEA/TEL data.  In some cases,  observers or
measurers or telescopes may change.  A comment line such as

COM Valid 1993 Mar 4 - 2012 Dec 25

tells the program:  "Don't pay any attention to the following lines unless
observations were made at this code during this time span." */

static int get_observer_details( const char *observation_filename,
      const OBSERVE *obs, const int n_obs,
      const char *mpc_code, char *observers, char *measurers, char *scope)
{
   FILE *ifile = fopen_ext( observation_filename, "fclrb");
   int rval = 0, n_codes_found = 0;

   *observers = *measurers = *scope = '\0';
   if( ifile)
      {
      char buff[700];
      bool done = false;

      while( !done && fgets( buff, sizeof( buff), ifile))
         if( !memcmp( buff, "COD ", 4))
            {
            bool new_code_found = false, use_lines = true;

            *observers = *measurers = *scope = '\0';
            n_codes_found++;
            if( !memcmp( buff + 4, mpc_code, 3))
               while( !done && fgets_trimmed( buff, sizeof( buff), ifile) && !new_code_found)
                  {
                  if( !memcmp( buff, "COM Valid:", 10))
                     {
                     char *tptr = strchr( buff, '-');
                     double jd_start, jd_end;

                     assert( tptr);
                     *tptr = '\0';
                     jd_start = get_time_from_string( 0., buff + 10, 0, NULL);
                     jd_end = get_time_from_string( 0., tptr + 1, 0, NULL);
                     assert( jd_start > 2000000. && jd_start < 3000000.);
                     assert( jd_end > 2000000. && jd_end < 3000000.);
                     use_lines = got_obs_in_range( obs, n_obs, jd_start, jd_end);
                     }
                  if( use_lines && !memcmp( buff, "OBS ", 4))
                     tack_on_names( observers, buff + 4);
                  if( use_lines && !memcmp( buff, "MEA ", 4))
                     tack_on_names( measurers, buff + 4);
                  if( use_lines && !memcmp( buff, "TEL ", 4))
                     strcat( scope, buff + 4);
                  if( !memcmp( buff, "COD ", 4))
                     {
                     if( memcmp( buff + 4, mpc_code, 3))
                        {
                        new_code_found = true;
                        done = true;
                        }
                     else
                        *observers = *measurers = *scope = '\0';
                     }
                  }
            }
      fclose( ifile);
      add_final_period( observers);
      add_final_period( measurers);
      add_final_period( scope);
      if( !strcmp( observers, measurers))
         *measurers = '\0';
      }

   if( *observers)
      rval = 1;
   if( *measurers)
      rval |= 2;
   if( *scope)
      rval |= 4;
   if( !n_codes_found)        /* we can just ignore this file completely, */
      rval = -1;              /* even for other observatory codes */
   return( rval);
}

static int get_observer_details_from_obs( const OBSERVE *obs,
      size_t n_obs, const char *mpc_code, char *observers,
      char *measurers, char *scope)
{
   int rval = 0;
   size_t i;
   const char *tptr;

   *observers = *measurers = *scope = '\0';
   while( n_obs--)
      {
      if( obs->obs_details && !strcmp( mpc_code, obs->mpc_code))
         for( i = 0; (tptr = obs->obs_details[i]) != NULL; i++)
            {
            if( !memcmp( tptr, "OBS ", 4))
               tack_on_names( observers, tptr + 4);
            if( !memcmp( tptr, "MEA ", 4))
               tack_on_names( measurers, tptr + 4);
            if( !memcmp( tptr, "TEL ", 4))  /* allow for only one scope */
               tack_on_names( scope, tptr + 4);
            }
      obs++;
      }
   add_final_period( observers);
   add_final_period( measurers);
   add_final_period( scope);
   if( *observers)
      rval = 1;
   if( *measurers)
      rval |= 2;
   if( *scope)
      rval |= 4;
   if( !strcmp( observers, measurers))
      *measurers = '\0';
   return( rval);
}

#define REPLACEMENT_COLUMN 42

static void observer_link_substitutions( char *buff)
{
   FILE *ifile = fopen_ext( "observer.txt", "fcrb");

   if( ifile)
      {
      char line[200], *loc;

      while( fgets_trimmed( line, sizeof( line), ifile))
         if( *line != ';' && *line != '#')
            {
            line[REPLACEMENT_COLUMN - 1] = '\0';
            remove_trailing_cr_lf( line);
            loc = strstr( buff, line);
            if( loc)
               {
               const size_t len = strlen( line);
               size_t len2;

               if( loc[len] <= ' ' || loc[len] == '.' || loc[len] == ',')
                  {
                  len2 = strlen( line + REPLACEMENT_COLUMN);
                  memmove( loc + len2, loc + len, strlen( loc + len) + 2);
                  memcpy( loc, line + REPLACEMENT_COLUMN, len2);
                  }
               }
            }
      fclose( ifile);
      }
}

static unsigned get_list_of_stations( const unsigned n_obs,
               const OBSERVE FAR *obs_data, const unsigned max_n_stations,
               char stations[][5])
{
   unsigned n_stations = 0, i, j;

   for( i = 0; i < n_obs; i++)
      {
      int compare = 1;

      j = 0;
      while( j < n_stations &&
             (compare = strcmp( obs_data[i].mpc_code, stations[j])) > 0)
         j++;
      if( compare)         /* got a new one */
         {
         assert( strlen( obs_data[i].mpc_code) == 3);
         memmove( stations + j + 1, stations + j, (n_stations - j)  * sizeof( stations[0]));
         strcpy( stations[j], obs_data[i].mpc_code);
         n_stations++;
         assert( n_stations < max_n_stations);
         }
      }
   return( n_stations);
}

static int write_observer_data_to_file( FILE *ofile, const char *ast_filename,
                 const int n_obs, const OBSERVE FAR *obs_data)
{
   unsigned n_stations = 0, i, j;
   int try_ast_file = 1, try_details_file = 1, try_scope_file = 1;
   char stations[400][5];

   INTENTIONALLY_UNUSED_PARAMETER( ast_filename);
   n_stations = get_list_of_stations( n_obs, obs_data, 400, stations);
   for( i = 0; i < n_stations; i++)
      {
      char buff[200], tbuff[100];
      char details[4][310];
      int details_found = 0;
      size_t loc;

      FSTRCPY( tbuff, stations[i]);
      put_observer_data_in_text( tbuff, buff);
      snprintf( details[0], sizeof( details[0]), "(%s) %s.", tbuff, buff);

      if( try_ast_file)
         {
         details_found = get_observer_details_from_obs( obs_data, n_obs, tbuff,
                                 details[1], details[2], details[3]);
         if( details_found == -1)       /* file wasn't found,  or it has */
            {                           /* no observational details.  In */
            details_found = 0;          /* either case,  ignore it for   */
            try_ast_file = 0;           /* further MPC codes.            */
            }
         }

      if( !details_found && try_details_file)
         {
         details_found = get_observer_details( "details.txt", obs_data,
                   n_obs, tbuff, details[1], details[2], details[3]);
         if( details_found == -1)       /* file wasn't found,  or it has */
            {                           /* no observational details.  In */
            details_found = 0;          /* either case,  ignore it for   */
            try_details_file = 0;       /* further MPC codes.            */
            }
         }

      if( !details_found && try_scope_file)
         {
         details_found = get_observer_details( "scopes.txt", obs_data,
                   n_obs, tbuff, details[1], details[2], details[3]);
         if( details_found == -1)       /* file wasn't found,  or it has */
            {                           /* no observational details.  In */
      /*    details_found = 0;   */     /* either case,  ignore it for   */
            try_scope_file = 0;         /* further MPC codes.            */
            }
         }

      loc = 0;
      for( j = 0; j < 4; j++)
         if( *details[j])
            {
            char inserted_text[15], *outtext = details[j];

            if( j == 3)                      /* telescope */
               strcpy( inserted_text, " ");
            else if( j == 0)                 /* observatory name/location */
               *inserted_text = '\0';
            else              /* j=1: observer(s); j=2: measurer(s) */
               {
               strcpy( inserted_text, (j == 2) ? " Measurer" : "  Observer");
               if( strchr( outtext, ','))
                  strcat( inserted_text, "s");
               strcat( inserted_text, " ");
               }
            memmove( outtext + strlen( inserted_text), outtext,
                              strlen( outtext) + 1);
            memcpy( outtext, inserted_text, strlen( inserted_text));
            while( *outtext)
               {
               int len = (int)strlen( outtext);

               if( len > 78 - (int)loc)
                  {
                  int k = 78 - (int)loc;

                  while( k && outtext[k - 1] != ' ')
                     k--;
                  len = k;
                  while( k && outtext[k - 1] == ' ')
                     k--;
                  fprintf( ofile, "%.*s\n   ", k, outtext);
                  loc = 3;
                  }
              else
                  {
                  fprintf( ofile, "%s", outtext);
                  loc += len;
                  }
               outtext += len;
               }
            }
      fprintf( ofile, "\n");
      }
   return( 0);
}

bool residual_file_in_config_dir = true;

int write_residuals_to_file( const char *filename, const char *ast_filename,
       const int n_obs, const OBSERVE FAR *obs_data, const int resid_format)
{
   FILE *ofile = fopen_ext( filename,
               residual_file_in_config_dir ? "tfcw" : "fw");
   int rval = 0;

   if( ofile )
      {
      char buff[200];
      int number_lines = (n_obs + 2) / 3;
      int i;

      if( (resid_format & 3) == RESIDUAL_FORMAT_SHORT)
         for( i = 0; i < number_lines * 3; i++)
            {
            int num = (i % 3) * number_lines + i / 3;
            OBSERVE FAR *obs = ((OBSERVE FAR *)obs_data) + num;

            if( num < n_obs)
               format_observation( obs, buff, resid_format);
            else
               *buff = '\0';
            fprintf( ofile, "%s%s", buff, (i % 3 == 2) ? "\n" : "   ");
            }
      else
         for( i = 0; i < n_obs; i++)
            {
            format_observation( obs_data + i, buff, resid_format);
            fprintf( ofile, "%s\n", buff);
            }
      fprintf( ofile, "\nStation data:\n");
      write_observer_data_to_file( ofile, ast_filename, n_obs, obs_data);
      fclose( ofile);
      }
   else                    /* file not opened */
      rval = -1;
   return( rval);
}

#ifdef FUTURE_PROJECT_IN_WORKS

int create_residual_scattergram( const char *filename, const int n_obs,
                         const OBSERVE FAR *obs)
{
   const int tbl_height = 19, tbl_width = 71;
   const int xspacing = 12, yspacing = 5,
   short *remap_table = (short *)calloc( tbl_height * tbl_width, sizeof( short));
   int i, j;
   FILE *ofile = fopen_ext( filename, "tfcwb");

   for( i = 0; i < n_obs; i++, obs++)
      {
      const double yresid = 3600. * (180./pi) * (obs->dec - obs->computed_dec);
      const double xresid = 3600. * (180./pi) * (obs->ra - obs->computed_ra)
                                        * cos( obs->computed_dec);
      int xloc = (int)floor( xresid * (double)xspacing / .5)
      }

   fclose( ofile);
   free( remap_table);
}
#endif

void remove_trailing_cr_lf( char *buff)
{
   int i;

   for( i = 0; buff[i] && buff[i] != 13 && buff[i] != 10; i++)
      ;
   while( i && buff[i - 1] == ' ')
      i--;
   buff[i] = '\0';
}

/* MPC frowns upon redistribution of NEOCP astrometry.  So if an input
line is from NEOCP,  it's blacked out,  _unless_ it's from a station
that has given permission for republication.  Those stations are listed
in the GREENLIT and GREENLIT2 lines in 'environ.dat';  you can add your
own if desired.

In some cases,  only heliocentric observations have gotten the green light.
Those codes are followed by an asterisk in the GREENLIT/GREENLIT2 lines.

For private use,  you can turn NEOCP redaction off... just be sure that
if you do that,  you don't redistribute anything.        */

bool neocp_redaction_turned_on = true;

static bool is_neocp_line( const char *mpc_line)
{
   return( strlen( mpc_line) == 80 && !memcmp( mpc_line + 72, "NEOCP", 5));
}

static bool line_must_be_redacted( const char *mpc_line,
                                                const bool is_heliocentric)
{
   bool rval = false;

   if( is_neocp_line( mpc_line) && neocp_redaction_turned_on)
      {
      const char *greenlit = strstr( get_environment_ptr( "GREENLIT"),
                                       mpc_line + 77);

      rval = true;
      if( !greenlit)        /* try an alternative string w/additional codes */
         greenlit = strstr( get_environment_ptr( "GREENLIT2"),
                               mpc_line + 77);
      if( greenlit && (is_heliocentric || greenlit[3] != '*'))
         rval = false;
      }
   return( rval);
}

int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr)
{
   size_t ilen = FSTRLEN( str), rval = 0;
   const size_t oldlen = strlen( oldstr);
   const size_t newlen = strlen( newstr);

   while( ilen >= oldlen)
      if( !FMEMCMP( str, oldstr, oldlen))
         {
         FMEMMOVE( str + newlen, str + oldlen, ilen - oldlen + 1);
         FMEMCPY( str, newstr, newlen);
         str += newlen;
         ilen -= oldlen;
         rval++;
         }
      else
         {
         str++;
         ilen--;
         }
   return( (int)rval);
}

static long round_off( const double ival, const double prec)
{
   long rval = 0L, digit;
   int got_it = 0;
   double diff;

   for( digit = 10000000L; !got_it; digit /= 10)
      {
      rval = (((long)ival + digit / 2L) / digit) * digit;
      diff = fabs( (double)rval - ival);
      if( digit == 1 || diff < ival * prec)
         got_it = 1;
      }
   return( rval);
}

/* Astrometry from NEOCP is not to be redistributed,  and is shown in
pseudo-MPECs as being blacked out.  See

https://www.projectpluto.com/redacted.htm

   for some explanation of this,  and of the following code.  It figures
out where in the redacted text to put links to the above explanation so
people will understand why some data is removed in this seemingly
strange manner.  The redacted terms are placed randomly and sorted to be
in the correct order.  If they'd overlap,  or run off the left edge,  we
try again,  up to max_iteration times. */

static inline void redacted_locations( const char *terms[],
              const unsigned n_redacted_lines, unsigned *x, unsigned *y)
{
   unsigned i, j, n_terms = 0, iteration;
   const unsigned max_iteration = 100;
   bool success = false;

   if( !n_redacted_lines)
      return;
   while( terms[n_terms])
      n_terms++;
   srand( (unsigned)time( NULL));
   for( iteration = 0; iteration < max_iteration && !success; iteration++)
      {
      const unsigned max_column = 51;

      success = true;
      for( i = 0; i < n_terms; i++)
         {
         x[i] = rand( ) % (max_column - 4);
         y[i] = rand( ) % n_redacted_lines;
         for( j = 0; j < i; j++)
            if( y[j] > y[i] || (y[j] == y[i] && x[j] > x[i]))
               {
               unsigned temp = x[i];   x[i] = x[j];  x[j] = temp;
               temp = y[i];   y[i] = y[j];   y[j] = temp;
               }
         }
      for( i = 0; i < n_terms; i++)
         {
         const unsigned end_x = x[i] + (unsigned)strlen( terms[i]) + 1;

         if( end_x > max_column)
            success = false;
         if( i < n_terms - 1 && y[i] == y[i + 1] && end_x > x[i + 1])
            success = false;        /* runs into next word */
         }
      }
   if( !success)        /* can't place the text;  just omit it */
      for( i = 0; i < n_terms; i++)
         y[i] = -9;
}

char *mpec_error_message = NULL;

int make_pseudo_mpec( const char *mpec_filename, const char *obj_name)
{
   char buff[800], mpec_buff[7];
   const char *mpec_permits = (strchr( mpec_filename, '/') ? "fwb" : "tfcwb");
   FILE *ofile = fopen_ext( mpec_filename, mpec_permits);
   FILE *header_htm_ifile;
   FILE *residuals_ifile, *ephemeris_ifile, *observations_ifile;
   FILE *elements_file = fopen_ext( get_file_name( buff, elements_filename), "tfcrb");
   int line_no = 0, rval = 0, total_lines = 0;
   unsigned mpec_no = atoi( get_environment_ptr( "MPEC"));
   bool orbit_is_heliocentric = true, suppressed = false;
   extern char findorb_language;
   unsigned redacted_line_number = 0, n_redacted_lines = 0;
   unsigned n_neocp_lines = 0;
   static const char *explanations_url = "https://www.projectpluto.com/mpec_xpl.htm";

   assert( ofile);
   setvbuf( ofile, NULL, _IONBF, 0);
   if( elements_file)
      {
      unsigned i;

      for( i = 0; i < 10 && fgets( buff, sizeof( buff), elements_file); i++)
         if( i == 1 && !memcmp( buff, "   Perihel", 10))
            i = 10;     /* yes,  it's definitely heliocentric,  no need to go on */
         else if( *buff == 'P')
            orbit_is_heliocentric = false;
      }

   if( mpec_no)
      snprintf( mpec_buff, 4, "_%02x", mpec_no % 256);
   else
      *mpec_buff = '\0';
   snprintf( buff, 12, "%cheader.htm", findorb_language);
   header_htm_ifile = fopen_ext( buff, "crb");
   if( !header_htm_ifile)
      header_htm_ifile = fopen_ext( buff + 1, "fcrb");
   assert( header_htm_ifile);

   observations_ifile = fopen_ext( get_file_name( buff, observe_filename), "tfcrb");
   assert( observations_ifile);

                  /* Count number of redacted and (current) NEOCP lines : */
   while( fgets_trimmed( buff, sizeof( buff), observations_ifile))
      if( is_neocp_line( buff))
         {
         if( memcmp( buff + 56, "Removed", 7))
            n_neocp_lines++;
         if( line_must_be_redacted( buff, orbit_is_heliocentric))
            n_redacted_lines++;
         }

   if( header_htm_ifile)                 /* copy header data to pseudo-MPEC */
      {
      while( fgets( buff, sizeof( buff), header_htm_ifile))
         if( *buff != '#')
            {
            char *tptr = strstr( buff, "_xx");

            if( tptr)
               {
               memmove( tptr + strlen( mpec_buff), tptr + 3, strlen( tptr));
               memcpy( tptr, mpec_buff, strlen( mpec_buff));
               }
            if( !memcmp( buff, "$Error", 6))
               {
               if( !mpec_error_message)
                  *buff = '\0';
               else
                  {
                  strcpy( buff, "<p> <b>");
                  strcat( buff, mpec_error_message);
                  strcat( buff, "</b> </p>");
                  }
               }
            while( (tptr = strchr( buff, '$')) != NULL)
               {                       /* See comments in 'header.htm'.  */
               int got_it = 0, i;      /* code replaces text between $s  */
                                       /* in that file.                  */
               for( i = 1; tptr[i] && tptr[i] != '$'; i++)
                  ;        /* search for matching $ */
               if( i < 20 && tptr[i] == '$')
                  {
                  char search_str[80], replace_str[180];

                  memcpy( search_str, tptr, i);
                  search_str[i] = '\0';
                  if( elements_file)
                     {
                     char tbuff[300], *tptr2;

                     if( !strcmp( search_str, "$Tg"))
                        {
                        full_ctime( replace_str, current_jd( ),
                                       FULL_CTIME_YMD);
                        got_it = 1;
                        }
                     else if( !strcmp( search_str, "$Name"))
                        {
                        strcpy( replace_str, obj_name);
                        got_it = 1;
                        }
                     else if( !strcmp( search_str, "$SV"))   /* state vect */
                        {                           /* for Orbit Simulator */
                        extern double helio_ecliptic_j2000_vect[];
                        int year = (int)JD_TO_YEAR(
                                          helio_ecliptic_j2000_vect[6] + 182.6);

                        if( year < 1950)    /* Orbit Simul has precomputed */
                           year = 1950;     /* solar syst data from 1950   */
                        if( year > 2050)    /* to 2050 */
                           year = 2050;
                        snprintf( replace_str, sizeof( replace_str),
                                        "%d.html?sv,1,%s,%.2f", year, obj_name,
                                        helio_ecliptic_j2000_vect[6]);
                        text_search_and_replace( replace_str, " ", "%20");
                        for( i = 0; i < 6; i++)
                           {
                           double oval = helio_ecliptic_j2000_vect[i] * AU_IN_METERS;

                           if( i >= 3)    /* velocity,  in km/s,  is desired */
                              oval /= seconds_per_day;
                           snprintf_append( replace_str, sizeof( replace_str),
                                    ",%f", oval);
                           }
                        got_it = 1;
                        }

                     fseek( elements_file, 0L, SEEK_SET);
                     while( !got_it &&
                             fgets_trimmed( tbuff, sizeof( tbuff), elements_file))
                        if( (tptr2 = strstr( tbuff, search_str)) != NULL
                                    && tptr2[i] == '=')
                           {
                           tptr2 += i + 1;
                           for( i = 0; tptr2[i] > ' '; i++)
                              replace_str[i] = tptr2[i];
                           replace_str[i] = '\0';
                           got_it = 1;
                           }
                     strcat( search_str, "$");
                     if( got_it)
                        text_search_and_replace( buff, search_str, replace_str);
                     }
                  }
               else        /* no matching '$' found */
                  *tptr = '!';
               }
            if( !suppressed)
               fputs( buff, ofile);
            }
         else if( !memcmp( buff, "# helio_only", 12) && !orbit_is_heliocentric)
            suppressed = (buff[13] == '1');
         else if( !memcmp( buff, "# neocp_only", 12))
            if( !n_neocp_lines || strlen( obj_name) > 7)
               suppressed = (buff[13] == '1');

      fclose( header_htm_ifile);
      }

   if( mpec_no)
      {
      snprintf( buff, 4, "%u", mpec_no % 255 + 1);
      set_environment_ptr( "MPEC", buff);
      }

   if( observations_ifile)
      {
      unsigned max_term = 5, x[10], y[10];
      const char *terms[] = { "Astrometry", "redacted;",
                      "click", "here", "for", "explanation", NULL };

      redacted_locations( terms, n_redacted_lines, x, y);
      fseek( observations_ifile, 0L, SEEK_SET);
      while( fgets_trimmed( buff, sizeof( buff), observations_ifile))
         if( memcmp( buff, "COM ", 4))      /* skip comment/'sigma' lines */
            {
//          if( buff[14] == 's' || buff[14] == 'v' || buff[14] == 'r')
//             fprintf( ofile, "%s\n", buff);
//          else
               {
               char mpc_code[8];
               const bool redacted = line_must_be_redacted( buff,
                                             orbit_is_heliocentric);

               strcpy( mpc_code, buff + 77);
               buff[77] = '\0';
               if( buff[14] != 's' && buff[14] != 'v' && buff[14] != 'r')
                  total_lines++;
               fprintf( ofile, "<a name=\"o%s%03d\"></a><a href=\"#r%s%03d\">%.12s</a>",
                        mpec_buff, total_lines, mpec_buff, total_lines, buff);
               if( redacted)
                  {
                  int i;
                  const size_t start_of_redacted_text = 25;
                  const size_t length_of_redacted_text = 77 - start_of_redacted_text;
                  char *tptr = buff + start_of_redacted_text;

                  strcpy( tptr, "<code class=\"neocp\">");
                  tptr += strlen( tptr);
                  memset( tptr, '~', length_of_redacted_text);
                  strcpy( tptr + length_of_redacted_text, "</code>");
                  for( i = max_term; i >= 0; i--)
                     if( redacted_line_number == y[i])
                        {
                        char tbuff[180], *zptr;

                        strcpy( tbuff, "</code><a href='https://www.projectpluto.com/redacted.htm'>");
                        strcat( tbuff, terms[i]);
                        strcat( tbuff, "</a><code class=\"neocp\">");
                        zptr = tptr + 2 + x[i];
                        memcpy( zptr, terms[i], strlen( terms[i]));
                        text_search_and_replace( tptr, terms[i], tbuff);
                        }
                  for( i = 0; tptr[i]; i++)  /* replace all the tildes with */
                     if( tptr[i] == '~')     /* pseudorandom text, but skip */
                        {                    /* chars that can't be used in */
                        const char *forbidden = "~ <>\"&";          /* HTML */

                        while( strchr( forbidden, tptr[i]))
                           tptr[i] = (char)( ' ' + rand( ) % ('z' - ' '));
                        }
                  redacted_line_number++;
                  }
               else
                  {
                  text_search_and_replace( buff + 13, "&", "&amp;");
                  text_search_and_replace( buff + 13, "<", "&lt;");
                  text_search_and_replace( buff + 13, ">", "&gt;");
                  }
               text_search_and_replace( buff + 13, "JPLRS",
                     "<a href=\"http://ssd.jpl.nasa.gov/?radar\">JPLRS</a>");
               fprintf( ofile, "%s<a href=\"#stn_%s\">%s</a>\n",
                        buff + 12, mpc_code, mpc_code);
               }
            }
      fclose( observations_ifile);
      }

   residuals_ifile = fopen_ext( get_file_name( buff, residual_filename), "tfcrb");
   if( residuals_ifile)
      {
      FILE *obslinks_file = fopen_ext( "obslinks.htm", "fcrb");
      FILE *mpc_obslinks_file = fopen_ext( "ObsCodesF.html", "fcrb");
      long obslinks_header_len = 0, mpc_obslinks_header_len = 0;
      char url[200];
              /* In making a pseudo-MPEC,  we attempt to include links     */
              /* to the Web sites of the observatory codes mentioned.  The */
              /* files 'obslinks.htm' (my effort at a complete list of Web */
              /* sites of observatory codes) and 'ObsLinksF.html' (the     */
              /* MPC's list) are both searched,  in that order.            */

      while( fgets( url, sizeof( url), obslinks_file)
                     && memcmp( url, "<a name=\"0\">", 12))
         ;
      obslinks_header_len = ftell( obslinks_file);
      while( fgets( url, sizeof( url), mpc_obslinks_file)
                     && memcmp( url, "<pre>", 5))
         ;
      mpc_obslinks_header_len = ftell( mpc_obslinks_file);
      while( fgets( buff, sizeof( buff), residuals_ifile) && memcmp( buff, "Station", 7))
         ;
      fprintf( ofile, "<a name=\"stations\"></a>\n");
      fprintf( ofile, "<a href=\"%s#stns\"><b>%s</b></a>", explanations_url, buff);
      while( fgets_trimmed( buff, sizeof( buff), residuals_ifile))
         if( *buff == ' ')
            {
            observer_link_substitutions( buff);
            assert( strlen( buff) < sizeof( buff));
            fprintf( ofile, "%s\n", buff);
            }
         else
            {
            char tbuff[4];
            int compare = 1, url_index = 0;
            char *tptr, saved_char;
            bool got_lat_lon = false;

            memcpy( tbuff, buff + 1, 3);
            tbuff[3] = '\0';
            fprintf( ofile, "<a name=\"stn_%s\"></a>", tbuff);

            tptr = strstr( buff, "  (N");
            if( !tptr)
               tptr = strstr( buff, "  (S");
            if( tptr && strchr( tptr, ')'))
               got_lat_lon = true;
            if( !tptr)
               tptr = strstr( buff, ". ");
            if( !tptr)
               {
               tptr = buff + strlen( buff);
               if( tptr[-1] == '.')
                  tptr--;
               }
                        /* At this point,  tptr should point at the end */
                        /* of the observatory name. */

/*          if( compare)         */
               {
               char text_to_find[50], *tptr;

               snprintf( text_to_find, sizeof( text_to_find),
                                  "></a> %.3s  <", buff + 1);
               url[19] = '\0';
               fseek( obslinks_file, obslinks_header_len, SEEK_SET);
               while( (compare = memcmp( url + 13, text_to_find, 12)) != 0 &&
                                 fgets_trimmed( url, sizeof( url), obslinks_file))
                  ;
               tptr = strstr( url, "<br>");
               if( tptr)
                  *tptr = '\0';
               url_index = 23;   /* if there is a link,  it starts in byte 23 */
               }
            if( compare)   /* still don't have URL;  try ObsLinks.html */
               {
               *url = '\0';
               fseek( mpc_obslinks_file, mpc_obslinks_header_len, SEEK_SET);
               while( (compare = memcmp( url, buff + 1, 3)) != 0 &&
                                 fgets_trimmed( url, sizeof( url), mpc_obslinks_file))
                  ;
               url_index = 32;   /* if there is a link,  it starts in byte 32 */
               }

            saved_char = *tptr;
            *tptr = '\0';
            if( !compare)   /* we got a link to an observatory code */
               {
               buff[5] = '\0';
               fprintf( ofile, "%s%s", buff, url + url_index);
               }
            else
               fprintf( ofile, "%s", buff);
            *tptr = saved_char;

            if( got_lat_lon)
               {
               double lat, lon;
               char lat_sign, lon_sign;
               char *new_tptr = strchr( tptr, ')');
               int n_scanned;

               assert( new_tptr);
               n_scanned = sscanf( tptr + 3, "%c%lf %c%lf", &lat_sign, &lat, &lon_sign, &lon);
               if( n_scanned != 4)
                  printf( "%s\n", tptr);
               assert( n_scanned == 4);
               if( lat_sign == 'S')
                  lat = -lat;
               if( lon_sign == 'W')
                  lon = -lon;
               fprintf( ofile, "  (<a title=\"Click for map\"");
               fprintf( ofile, " href=\"http://maps.google.com/maps?q=%.5f,+%.5f\">",
                              lat, lon);
               *new_tptr = '\0';
               fprintf( ofile, "%s</a>)", tptr + 3);
               tptr = new_tptr + 1;
               }
            observer_link_substitutions( tptr);
            fprintf( ofile, "%s\n", tptr);
            }
      fclose( obslinks_file);
      fclose( mpc_obslinks_file);
      }
   else
      rval |= 2;

   if( elements_file)
      {
      bool in_comments = false;

      fseek( elements_file, 0L, SEEK_SET);
      fprintf( ofile, "<a name=\"elements%s\"></a>\n", mpec_buff);
      line_no = 0;
      while( fgets_trimmed( buff, sizeof( buff), elements_file))
         if( *buff != '#')
            {
            char *h_ptr = NULL;

            if( buff[19] == 'H')
               h_ptr = buff + 20;
            if( buff[27] == 'H')
               h_ptr = buff + 28;
            if( !line_no)
               fprintf( ofile, "<a href=\"%s#elems\"><b>%s</b></a>\n",
                                explanations_url, buff);
            else if( *buff == 'P' && h_ptr)
               {
               const double abs_mag = atof( h_ptr);
                        /* H=4 indicates 420 to 940 km,  so: */
               double upper_size = 940. * exp( (4. - abs_mag) * LOG_10 / 5.);
               const char *units = "km";
               const char *size_url =
                   "href=\"https://www.minorplanetcenter.net/iau/lists/Sizes.html\">";
               char title[50];

               h_ptr[-1] = '\0';
               if( upper_size < .004)   /* under four meters,  use cm as units: */
                  {
                  upper_size *= 1000. * 100.;
                  units = "cm";
                  }
               else if( upper_size < 4.)  /* under four km,  use meters: */
                  {
                  upper_size *= 1000.;
                  units = "meters";
                  }
               snprintf( title, sizeof( title),
                          "\"Size is probably %ld to %ld %s\"\n",
                          round_off( upper_size / sqrt( 5.), .1),
                          round_off( upper_size, .1), units);
               fprintf( ofile, "%s<a title=%s%sH</a>%s\n",
                           buff, title, size_url, h_ptr);
               }
            else
               {
               text_search_and_replace( buff, "<HUGE>", "&lt;HUGE&gt;");
               text_search_and_replace( buff, "m^2", "m<sup>2</sup>");
               text_search_and_replace( buff, "   Find_Orb",
                                "   <a href=\"https://www.projectpluto.com/find_orb.htm\">Find_Orb</a>");
               fprintf( ofile, "%s\n", buff);
               }
            line_no++;
            }
         else          /* put 'elements.txt' comments in the pseudo-mpec, */
            {                          /* still as comments */
            if( !in_comments)
               fprintf( ofile, "%s", "</pre> ");
            fprintf( ofile, "<!-- %s -->\n", buff + 2);
            in_comments = true;
            }
      fclose( elements_file);
      }
   else
      rval |= 4;

               /* _now_ write out residuals: */
   if( residuals_ifile)
      {
      fseek( residuals_ifile, 0L, SEEK_SET);
      fprintf( ofile, "<pre><b><a name=\"residuals%s\"></a>"
                      "<a href=\"%s#resids\">"
                      "Residuals in arcseconds:</a> </b>\n",
                           mpec_buff, explanations_url);
      line_no = 0;
      while( fgets( buff, sizeof( buff), residuals_ifile) && *buff > ' ')
         {
         int i, column_off = (total_lines + 2) / 3, line;

         line_no++;
         for( i = 0; i < 3; i++)
            if( (line = line_no + column_off * i) <= total_lines)
               {
               char *tptr = buff + i * 26, tbuff[20];

               tptr[6] = '\0';         /* put out the YYMMDD... */
               fprintf( ofile, "<a name=\"r%s%03d\"></a><a href=\"#o%s%03d\">%s</a>",
                        mpec_buff, line, mpec_buff, line, tptr);

               memcpy( tbuff, tptr + 7, 3);     /* ...then the obs code.. */
               tbuff[3] = '\0';
               fprintf( ofile, " <a href=\"#stn_%s\">%s</a>", tbuff, tbuff);

               tptr[23] = '\0';        /* ...and finally,  the residuals */
               strcpy( tbuff, tptr + 10);
               text_search_and_replace( tbuff, "u", "&#xb5;");
               fprintf( ofile, "%s   ", tbuff);
               }
         fprintf( ofile, "\n");
         }
      fclose( residuals_ifile);
      }

               /* ...and now,  the ephemeris: */
   ephemeris_ifile = fopen_ext( get_file_name( buff, ephemeris_filename),
               is_default_ephem ? "tcr" : "r");
   if( ephemeris_ifile && fgets_trimmed( buff, sizeof( buff), ephemeris_ifile))
      {
      fprintf( ofile, "\n<a name=\"eph%s\"></a>", mpec_buff);
      fprintf( ofile, "<a href=\"%s#ephems\">", explanations_url);
      if( *buff != '#')        /* non-observables ephemeris,  no MPC code */
         fprintf( ofile, "<b>Ephemerides:</b></a>\n");
      else if( !memcmp( buff + 2, "500", 3))
         fprintf( ofile, "<b>Ephemerides (geocentric):</b></a>\n");
      else
         fprintf( ofile, "<b>Ephemerides for %s:</b></a>\n", buff + 1);

      while( fgets( buff, sizeof( buff), ephemeris_ifile))
         {
         char *color = strchr( buff, '$');

         if( color)
            {
            unsigned rgb, sum_components;
            char replace[80];
            const char *format;

            sscanf( color + 1, "%6x", &rgb);
            sum_components = (rgb & 0xff) + ((rgb >> 8) & 0xff) + (rgb >> 16);
            if( sum_components < 0x180)      /* dark color = use white text */
               format = "<a class=\"whtext\" style=\"background-color:#%06x;\">";
            else           /* bright color:  normal (black) text  */
               format = "<a style=\"background-color:#%06x;\">";
            memmove( color + 1, color + 7, strlen( color + 6));  /* remove RGB */
            memmove( color + 8, color + 4, strlen( color + 9));
            memcpy( color + 4, "</a>", 4);                 /* insert end tag */
            snprintf( replace, sizeof( replace), format, rgb);
            text_search_and_replace( color, "$", replace);
            }
         fputs( buff, ofile);
         }
      fclose( ephemeris_ifile);
      }
   else
      rval |= 8;

   fprintf( ofile, "</pre></body></html>\n");
   fclose( ofile);
   return( rval);
}

