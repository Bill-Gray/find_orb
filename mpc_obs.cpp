/* mpc_obs.cpp: parsing/interpreting MPC and other observations

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

#ifdef __WATCOMC__
#include <io.h>          /* for unlink( ) prototype */
#endif
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#ifdef _WIN32
#include <windows.h>
#include <malloc.h>     /* for alloca() prototype */
#else
#ifdef __WATCOMC__
#include <sys/timeb.h>
#else
#include <sys/time.h>
#include <unistd.h>
#endif
#endif
#include <stdarg.h>
#include <assert.h>
#include "watdefs.h"
#include "comets.h"
#include "lunar.h"
#include "afuncs.h"
#include "mpc_obs.h"
#include "sigma.h"
#include "date.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_MINOR_AXIS 6356755.
#define EARTH_MAJOR_AXIS_IN_AU (EARTH_MAJOR_AXIS / AU_IN_METERS)
#define EARTH_MINOR_AXIS_IN_AU (EARTH_MINOR_AXIS / AU_IN_METERS)
#define J2000 2451545.

void set_up_observation( OBSERVE FAR *obs);                 /* mpc_obs.c */
static double observation_jd( const char *buff);
double centralize_ang( double ang);             /* elem_out.cpp */
int sort_obs_by_date_and_remove_duplicates( OBSERVE *obs, const int n_obs);
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
int lat_alt_to_parallax( const double lat, const double ht_in_meters,
             double *rho_cos_phi, double *rho_sin_phi, const int planet_idx);
int parallax_to_lat_alt( const double rho_cos_phi, const double rho_sin_phi,
       double *lat, double *ht_in_meters, const int planet_idx); /* ephem0.c */
void set_obs_vect( OBSERVE FAR *obs);        /* mpc_obs.h */
int planet_posn( const int planet_no, const double jd, double *vect_2000);
void remove_trailing_cr_lf( char *buff);            /* ephem0.cpp */
void format_dist_in_buff( char *buff, const double dist_in_au); /* ephem0.c */
char int_to_mutant_hex_char( const int ival);               /* mpc_obs.c */
double current_jd( void);                       /* elem_out.cpp */
void remove_insignificant_digits( char *tbuff);          /* monte0.c */
int compute_observer_loc( const double jde, const int planet_no,
             const double rho_cos_phi,           /* mpc_obs.cpp */
             const double rho_sin_phi, const double lon, double FAR *offset);
int compute_observer_vel( const double jde, const int planet_no,
             const double rho_cos_phi,           /* mpc_obs.cpp */
             const double rho_sin_phi, const double lon, double FAR *vel);
int get_residual_data( const OBSERVE *obs, double *xresid, double *yresid);
static int xref_designation( char *desig);
int64_t nanoseconds_since_1970( void);                      /* mpc_obs.c */
int debug_printf( const char *format, ...);                 /* mpc_obs.c */
char **load_file_into_memory( const char *filename, size_t *n_lines);
void shellsort_r( void *base, const size_t n_elements, const size_t esize,
         int (*compare)(const void *, const void *, void *), void *context);
int string_compare_for_sort( const void *a, const void *b, void *context);
int format_jpl_ephemeris_info( char *buff);                 /* pl_cache.c */
int set_tholen_style_sigmas( OBSERVE *obs, const char *buff);  /* mpc_obs.c */
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
#ifdef _MSC_VER
     /* Microsoft Visual C/C++ has no snprintf.  See 'ephem0.cpp'.  */
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif
int snprintf_append( char *string, const size_t max_len,      /* ephem0.cpp */
                                   const char *format, ...);
char *mpc_station_name( char *station_data);       /* mpc_obs.cpp */

int debug_printf( const char *format, ...)
{
   FILE *ofile = fopen_ext( "debug.txt", "ca");

   if( ofile)
      {
      va_list argptr;
      const time_t t0 = time( NULL);

      fprintf( ofile, "%02d:%02d:%02d ",
               ((int)t0 / 3600) % 24, ((int)t0 / 60) % 60, (int)t0 % 60);
      va_start( argptr, format);
      vfprintf( ofile, format, argptr);
      va_end( argptr);
      fclose( ofile);
      }
   return( 0);
}

/* 'nanoseconds_since_1970( )' returns something close to the result of   */
/* ctime( ),  except a billion times larger and with added precision.     */
/* Naturally,  DOS/Windows,  OS/X,  and POSIX systems each use different  */
/* methods.                                                               */
/*    As currently coded,  the actual precision provided is 10^-7 second  */
/* in Windows;  a millisecond with the WATCOM compiler;  and 10^-6 second */
/* in everything else.  But note that "true" nanosecond precision is      */
/* possible,  if actually desired (see the NOT_CURRENTLY_IN_USE code).    */
/*    The range of a 64-bit signed integer is large enough to enable      */
/* this function to work until Friday, 2262 Apr 11 23:47:16.  We can get  */
/* an addition 292 years by using unsigned integers,  but it may be wiser */
/* to switch to 128-bit integers.                                         */

#ifdef _WIN32
int64_t nanoseconds_since_1970( void)
{
   FILETIME ft;
   const uint64_t jd_1601 = 2305813;  /* actually 2305813.5 */
   const uint64_t jd_1970 = 2440587;  /* actually 2440587.5 */
   const uint64_t ten_million = 10000000;
   const uint64_t diff = (jd_1970 - jd_1601) * ten_million * seconds_per_day;
   uint64_t decimicroseconds_since_1970;   /* i.e.,  time in units of 1e-7 seconds */

   GetSystemTimeAsFileTime( &ft);
   decimicroseconds_since_1970 = ((uint64_t)ft.dwLowDateTime |
                                ((uint64_t)ft.dwHighDateTime << 32)) - diff;
   return( decimicroseconds_since_1970 * (int64_t)100);
}
#else
#ifdef __WATCOMC__
int64_t nanoseconds_since_1970( void)
{
   struct timeb t;
   const int64_t one_million = 1000000;
   int64_t millisec;

   ftime( &t);
   millisec = (int64_t)t.millitm + (int64_t)1000 * (int64_t)t.time;
   return( millisec * (int64_t)one_million);
}
#else      /* OS/X,  BSD,  and Linux */
int64_t nanoseconds_since_1970( void)
{
   struct timeval now;
   const int rv = gettimeofday( &now, NULL);
   int64_t rval;
   const int64_t one_billion = (int64_t)1000000000;

   if( !rv)
      rval = (int64_t)now.tv_sec * one_billion
           + (int64_t)now.tv_usec * (int64_t)1000;
   else
      rval = 0;
   return( rval);
}
#endif
#endif

/* At one time,  I was using the following in Linux.  It gives a
"real" precision of nanoseconds,  instead of getting microseconds
and multiplying by 1000 (or decimicroseconds and multiplying by 100).
However,  it does require the realtime library to be linked in...
I leave it here in case we someday need nanosecond precision.  */

#ifdef NOT_CURRENTLY_IN_USE
int64_t nanoseconds_since_1970( void)
{
   struct timespec t;

   clock_gettime( CLOCK_REALTIME, &t);
   return( t.tv_sec * (int64_t)1000000000 + t.tv_nsec);
}
#endif    /* NOT_CURRENTLY_IN_USE */

      /* If you save NEOCP astrometry as a 'Web page,  complete',  the */
      /* first line is prefaced with HTML tags.  The following removes */
      /* any HTML tags,  which should also allow you to load up        */
      /* astrometry from pseudo-MPECs and maybe other HTML-ized data.  */

static void remove_html_tags( char *buff)
{
   char *left_angle_bracket, *right_angle_bracket;

   while( (left_angle_bracket = strchr( buff, '<')) != NULL
        && (right_angle_bracket = strchr( left_angle_bracket, '>')) != NULL)
      {
      memmove( left_angle_bracket, right_angle_bracket + 1,
                     strlen( right_angle_bracket));
      }
}

      /* In some situations,  MPC lines end up with a bit of leading  */
      /* garbage.  My hope is that this function will fix most cases. */

static void fix_up_mpc_observation( char *buff)
{
   const size_t len = strlen( buff);

   if( len == 80 && observation_jd( buff))      /* doesn't need fixing */
      {                                      /* except maybe for desig */
      if( buff[7] == ' ' && buff[6] != ' ')
         {
         char tbuff[20];

         if( sscanf( buff, "%19s", tbuff) == 1 && strlen( tbuff) == 7)
            {
            memset( buff, ' ', 5);
            memcpy( buff + 5, tbuff, 7);
            }
         }
      return;
      }

   if( len < 90)     /* avoid buffer overruns */
      {
      char desig[100], year[10], month[10], day[10];
      int bytes_read;

      if( sscanf( buff, "%99s %9s %9s %9s%n", desig, year, month, day,
                    &bytes_read) == 4 && strlen( month) == 2)
         {
         const size_t desig_len = strlen( desig);
         const size_t year_len = strlen( year);
         const size_t day_len = strlen( day);

         if( desig_len < 10 && year_len >= 5 && year_len < 7)
            {
            char obuff[81];
            char tbuff[10], minutes[10], seconds[10];
            int tval;

            memset( obuff, ' ', 80);
            obuff[80] = '\0';

            if( desig_len == 7)     /* preliminary */
               memcpy( obuff + 5, desig, 7);
            else
               memcpy( obuff, desig, desig_len);
            memcpy( obuff + 19 - year_len, year, year_len);
            memcpy( obuff + 20, month, 2);
            memcpy( obuff + 23, day, day_len);
                     /* Normally,  there will be a space between the day */
                     /* and the RA hours fields.  But if the day is given */
                     /* to six places,  they'll merge. */
            if( day_len <= 8)
               {
               if( sscanf( buff + bytes_read, "%9s%n", tbuff, &tval) == 1
                        && strlen( tbuff) == 2)
                  {
                  memcpy( obuff + 32, tbuff, 2);
                  bytes_read += tval;
                  }
               else        /* formatting trouble */
                  return;
               }
            if( sscanf( buff + bytes_read, "%9s%9s%n", minutes, seconds,
                           &tval) != 2)
               return;
            bytes_read += tval;
            memcpy( obuff + 35, minutes, 2);  /* RA minutes */
            memcpy( obuff + 38, seconds, strlen( seconds));
                  /* Again,  it's possible for the RA seconds to run */
                  /* right into the declination degrees.  You can't count */
                  /* on a separator being there. */
            if( strlen( seconds) < 6)
               {
               if( sscanf( buff + bytes_read, "%9s%n", tbuff, &tval) == 1
                           && strlen( tbuff) == 3)
                  {
                  memcpy( obuff + 44, tbuff, 3);
                  bytes_read += tval;
                  }
               else        /* formatting trouble */
                  return;
               }
            if( sscanf( buff + bytes_read, "%9s%9s%n", minutes, seconds,
                           &tval) != 2)
               return;
            bytes_read += tval;
            memcpy( obuff + 48, minutes, 2);  /* dec minutes */
            memcpy( obuff + 51, seconds, strlen( seconds));
            if( sscanf( buff + bytes_read, "%9s%n", tbuff, &tval) == 1
                           && isdigit( tbuff[0]) && isdigit( tbuff[1])
                           && (!tbuff[2] || tbuff[2] == '.'))
               {                 /* got a magnitude value: */
               memcpy( obuff + 65, tbuff, strlen( tbuff));
               bytes_read += tval;
                        /* Might get a mag band: */
               obuff[70] = buff[bytes_read + 1];
               }
                     /* figure out mag bands later... */
            memcpy( obuff + 77, buff + strlen( buff) - 3, 3);
//          printf( "Line in:\n'%s'\n", buff);
            strcpy( buff, obuff);
//          printf( "Line out:\n'%s'\n", buff);
            }
         }
      }
}

/* get_ra_dec() looks at an RA or dec from an MPC report and returns     */
/* its precision.  It interprets the formats used by MPC,  plus a lot of */
/* "extended" formats that can be useful if your input data is in other  */
/* formats and/or has extra digits.  The return value has the following  */
/* meanings ('z' = 'hours or degrees',  i.e.,  this format can apply to  */
/* both RAs and decs.)                                                   */
/*    hh mm ss.sss       3    (MPC uses this for 'precise' RAs)          */
/*    zz mm ss.ss        2    (MPC uses for most RAs & 'precise' decs)   */
/*    zz mm ss.s         1    (MPC uses for most decs & low-precision RA)*/
/*    zz mm ss           0    (Used by MPC, only rarely)                 */
/*    zz mm             -1    (Used _very_ rarely by MPC)                */
/*    zz mm.m           -2    (Maybe used by MPC once or twice)          */
/*    zz mm.mm          -3       The following are for Find_Orb only     */
/*    zz mm.mmm         -4                                               */
/*    zz mm.mmmm        -5                                               */
/*    zz mm.mmmmm       -6                                               */
/*    zz mm.mmmmmm      -7                                               */
/*    zz.               100                                              */
/*    zz.z              101                                              */
/*    zz.zz             102                                              */
/*    zz.zzz            103                                              */
/*    zz.zzzz           104... can go up to nine places = 109 in RA,     */
/*                            or to 108 = eight places in dec            */
/*    ddd.              200  (used for RA only)                          */
/*    ddd.d             201                                              */
/*    ddd.dd            202... can go up to eight places = 208           */
/*    HHMMSSs           307 (RA to .1 second)                            */
/*    ddmmSSs           307 (dec to .1 arcsecond)                        */
/*    HHMMSSss          308 (RA to .01 second)                           */
/*    ddmmSSss          308 (dec to .01 arcsecond)                       */
/*            ... and so forth until...                                  */
/*    ddmmSSsssss       311 (dec to 10 microarcseconds)                  */
/*    HHMMSSssssss      312 (RA to one microsecond)                      */
/*    Undetermined      -99                                              */
/* Please note that (at least thus far) I've only seen the first six     */
/* cases used by MPC,  and they will probably balk at anything sent in   */
/* anything but the first four formats.                                  */
/*    The remaining formats have been quite useful in situations where I */
/* find data in a non-MPC format;  I can leave it in decimal degrees or  */
/* whatnot,  and can accommodate extra digits for super-high-accuracy    */
/* data.  That's why formats 307-312 were added;  they accommodate some  */
/* highly precise VLBA astrometry that really _is_ good to the tens of   */
/* microarcseconds level.  See                                           */
/*                                                                       */
/* http://iau-comm4.jpl.nasa.gov/plan-eph-data/vlbaobs.html              */
/*                                                                       */
/*    And Gaia,  at least in theory,  will be of a similar level of      */
/* accuracy;  we need to Be Prepared for that.                           */
/*    The precision is stored and used to "recreate" the RA/dec in the   */
/* original form.  (It could,  and probably should,  also be used to     */
/* weight the observations.)                                             */

#define BAD_RA_DEC_FMT           -99

static double get_ra_dec( const char *ibuff, int *format, double *precision)
{
   char buff[13];
   double rval;
   unsigned i = 0, n_digits = 0;
   const bool is_dec = (*ibuff == '-' || *ibuff == '+');
   const bool is_negative = (*ibuff == '-');

   *precision = 1.;   /* in arcseconds */
   if( is_dec)
      ibuff++;
   memcpy( buff, ibuff, 12);
   buff[12] = '\0';
   rval = atof( buff);
   while( isdigit( buff[i]))
      i++;
   if( i > 7)        /* "packed" highly precise RA/decs described above */
      {
      unsigned tval;
      double factor = 1. / 3600.;

      *format = 300 + i;
      n_digits = i - 6;
      buff[6] = '\0';
      tval = atoi( buff);
      rval = (double)( tval / 10000)
           + (double)( (tval / 100) % 100) / 60.
           + (double)( tval % 100) / 3600.;
      for( i = 6; i < 12 && isdigit( ibuff[i]); i++)
         {
         factor *= .1;
         rval += (double)( ibuff[i] - '0') * factor;
         }
      buff[6] = ibuff[6];
//    debug_printf( "Extended: '%s' = %.8f\n", ibuff, rval);
      }
   else if( buff[2] == '.')        /* decimal degrees or hours */
      {
      *precision = 3600.;
      for( i = 3; isdigit( buff[i]) && i < 12; i++)
         n_digits++;
      *format = 100 + n_digits;
      }
   else if( buff[3] == '.')        /* decimal degrees for RA,  ddd.ddd.... */
      {
      *precision = 3600. / 15.;
      rval /= 15.;
      for( i = 4; isdigit( buff[i]) && i < 12; i++)
         n_digits++;
      *format = 200 + n_digits;
      }
   else if( buff[2] == ' ')            /* zz mm(.mmm...) or zz mm ss(.sss...)  */
      {
      rval += atof( buff + 3) / 60.;
      if( buff[5] == ' ' && isdigit( buff[7]))  /* i.e., seconds are given */
         {
         rval += atof( buff + 6) / 3600.;
         for( i = 9; i < 12 && isdigit( buff[i]); i++)
            n_digits++;
         *format = n_digits;
         }
      else           /* minutes: */
         {
         for( i = 6; i < 12 && isdigit( buff[i]); i++)
            n_digits++;
         *format = -((int)n_digits + 1);
         *precision = 60.;
         }
      }
   else
      *format = BAD_RA_DEC_FMT;
   while( n_digits--)
      *precision *= .1;
   if( is_negative)
      rval = -rval;
   return( rval);
}

#ifndef _MSC_VER
   #define CONSOLE
      /* In the console version of Find_Orb,  the following two functions */
      /* get remapped to Curses functions.  In the non-interactive one,   */
      /* they're mapped to 'do-nothings'.  See fo.cpp & find_orb.cpp.     */
   void refresh_console( void);
   void move_add_nstr( const int col, const int row, const char *msg, const int n_bytes);
#endif

#ifdef CONSOLE
#define COLOR_DEFAULT_INQUIRY       9
int inquire( const char *prompt, char *buff, const int max_len,
                     const int color);
#endif

int set_tholen_style_sigmas( OBSERVE *obs, const char *buff)
{
   const int n_scanned = sscanf( buff, "%lf%lf%lf",
               &obs->posn_sigma_1, &obs->posn_sigma_2, &obs->posn_sigma_theta);

   if( n_scanned == 1)        /* just a circular error */
      obs->posn_sigma_2 = obs->posn_sigma_1;
   if( n_scanned != 3)        /* no position angle supplied (usual case) */
      obs->posn_sigma_theta = 0.;
   else
      obs->posn_sigma_theta *= PI / 180.;
   return( n_scanned);
}


static double grab_double( const char *iptr, const int field_size)
{
   char tbuff[20];

   memcpy( tbuff, iptr, field_size);
   tbuff[field_size] = '\0';
   return( atof( tbuff));
}

int generic_message_box( const char *message, const char *box_type)
{
   int rval = 0;
#ifdef CONSOLE
   inquire( message, NULL, 30, COLOR_DEFAULT_INQUIRY);
#else
   int box_flags = MB_YESNO;

   if( !strcmp( box_type, "o"))
      box_flags = MB_OK;
   rval = MessageBox( NULL, message, "Find_Orb", box_flags);
#endif
   debug_printf( "%s", message);
   return( rval);
}

static double center_around( double ival, const double center)
{
   while( ival > center + 180.)
      ival -= 360.;
   while( ival < center - 180.)
      ival += 360.;
   return( ival);
}

#define is_between( x, x1, x2) ((x < x1 && x > x2) || (x < x2 && x > x1))

static bool extract_region_data_for_mpc_station( char *buff,
            const double lat, const double lon)
{
   FILE *ifile = fopen_ext( "geo_rect.txt", "fcrb");
   const double lat_in_degrees = (180. / PI) * lat;
   const double lon_in_degrees = (180. / PI) * lon;

   *buff = '\0';
   if( ifile)
      {
      char tbuff[90];

      while( !*buff && fgets_trimmed( tbuff, sizeof( tbuff), ifile))
         if( *tbuff != '#')
            {
            double lon1 = center_around( atof( tbuff), lon);
            double lat1 = atof( tbuff + 10);
            double lon2 = center_around( atof( tbuff + 20), lon1);
            double lat2 = atof( tbuff + 30);

            if( is_between( lon_in_degrees, lon1, lon2) &&
                   is_between( lat_in_degrees, lat1, lat2))
               strcpy( buff, tbuff + 40);
            }
      fclose( ifile);
      }
   return( *buff ? true : false);
}

/* You can store locations in 'rovers.txt' in base-60 form,  with the lat */
/* and longitude smashed together;  e.g.,  19 13' 33.1" would be stored   */
/* as 191333.1.  The following code would take 191331.1 as input and return */
/* 19 + 13/60 + 33.1/3600 = 19.22586111.   */

static double convert_base_60_to_decimal( const double ival)
{
   const int deg = (int)( ival / 10000.);
   const int min = (int)( fmod( ival / 100., 100.));
   const double sec = fmod( ival, 100.);
   const double rval = (double)deg + (double)min / 60. + sec / 3600.;

   return( rval);
}

/* http://www.minorplanetcenter.net/iau/info/Astrometry.html#HowObsCode
suggests that you start out using "observatory code" (XXX),  whilst
including a comment such as

COM Long. 239 18 45 E, Lat. 33 54 11 N, Alt. 100m, Google Earth

   Find_Orb allows this,  but extends it to work with _any_ observatory
code,  existing or newly created,  or with code (XXX).  Thus,  a header
showing

COD Bow
COM Long. 69 54 00 W, Lat. 44 01 10 N, Alt. 50m, approximate

   would cause Find_Orb to add a new observatory code (Bow),  corresponding
to Bowdoinham,  Maine,  where the corporate headquarters of Project Pluto
are located.

   You _can_ use this to override the existing codes,  but the only case
I can think of where you'd do that would be if you thought there might
be an error in the MPC's ObsCodes.html listing.  (Or in the supplement
in 'rovers.txt'.)

   Up to MAX_XXX_CODES are kept track of,  using the following xxx_code_t
structure. */

typedef struct
{
   char mpc_code[4];
   double lon, lat, alt;
} xxx_code_t;

#define MAX_XXX_CODES 20

static xxx_code_t xxx_codes[MAX_XXX_CODES];
static char curr_xxx_code[8];

static inline void get_xxx_location( const char *buff)
{
   double deg, min, sec;
   char compass;
   const char *tptr;
   xxx_code_t * xxx;
   size_t i = 0;

   debug_printf( "Looking for code in '%s'\n", buff);
   while( i < MAX_XXX_CODES - 1 && xxx_codes[i].mpc_code[0]
                  && memcmp( curr_xxx_code, xxx_codes[i].mpc_code, 3))
      i++;
   xxx = &xxx_codes[i];
   memcpy( xxx->mpc_code, curr_xxx_code, 3);
   deg = min = sec = 0.;
   sscanf( buff + 9, "%lf %lf %lf %c", &deg, &min, &sec, &compass);
   xxx->lon = deg + min / 60. + sec / 3600.;
   if( xxx->lon > 180.)
      xxx->lon -= 360.;
   if( compass == 'W')
      xxx->lon *= -1.;
   tptr = strstr( buff, "Lat.");
   if( tptr)
      {
      deg = min = sec = 0.;
      sscanf( tptr + 4, "%lf %lf %lf %c", &deg, &min, &sec, &compass);
      xxx->lat = deg + min / 60. + sec / 3600.;
      if( compass == 'S' || compass == 's')
         xxx->lat = -xxx->lat;
      }
   tptr = strstr( buff, "Alt.");
   if( tptr)
      xxx->alt = atof( tptr + 4);
}

extern int debug_level;

/* A return value of -1 indicates malformed input.  Anything else indicates */
/* the number of the planet on which the station is located (earth=3,  the  */
/* usual value... but note that rovers.txt contains extraterrestrial "MPC   */
/* stations",  so values from 0 to 10 may occur.)                           */

static int extract_mpc_station_data( const char *buff, double *lon_in_radians,
                            double *rho_cos_phi, double *rho_sin_phi)
{
   int i, rval = -1;       /* assume no valid data */
   double tlon = 0.;

                    /* A valid station code line must have,  for starters, */
                    /* a three-character code followed by a space,  and at */
                    /* least 22 bytes:                                     */
   if( strlen( buff) < 22)
      return( -1);
   for( i = 0; i < 3; i++)
      if( buff[i] <= ' ' || buff[i] > 'z')
         return( -1);
   if( lon_in_radians && rho_cos_phi && rho_sin_phi)
      *lon_in_radians = *rho_cos_phi = *rho_sin_phi = 0.;
   if( buff[7] == '.' && buff[14] == '.' && buff[23] == '.')
      {              /* for-real MPC data,  with parallax constants: */
      tlon = grab_double( buff + 4, 9) * PI / 180.;
      if( rho_cos_phi && rho_sin_phi)
         {
         *rho_cos_phi = grab_double( buff + 13, 8) * EARTH_MAJOR_AXIS_IN_AU;
         *rho_sin_phi = grab_double( buff + 21, 9) * EARTH_MAJOR_AXIS_IN_AU;
         }
      rval = 3;
      }
   else if( !memcmp( buff + 3, "                           ", 27)
            && buff[30] > ' ')    /* satellite stations (SOHO, */
      rval = 3;                   /* WISE,  etc.:  all blanks) */
   else        /* latitude & altitude,  not parallaxes: */
      {
      double lat_in_radians, alt_in_meters;

      buff += 4;
      if( *buff == '!')       /* "long format" */
         buff++;
      if( sscanf( buff, "%lf %lf %lf", &tlon, &lat_in_radians,
                                     &alt_in_meters) == 3)
         {
         const char *tptr = strchr( buff, '@');

         rval = (tptr ? atoi( tptr + 1) : 3);
         if( fabs( tlon) > 361. || fabs( lat_in_radians) > 91.)
            {     /* lat/lon are actually in dddmmss.sss form: */
            tlon = convert_base_60_to_decimal( tlon);
            lat_in_radians = convert_base_60_to_decimal( lat_in_radians);
            }
         tlon *= PI / 180.;
         lat_in_radians *= PI / 180.;
         if( rho_cos_phi && rho_sin_phi)
            lat_alt_to_parallax( lat_in_radians, alt_in_meters,
                     rho_cos_phi, rho_sin_phi, rval);   /* ephem0.cpp */
         }
      }
   if( lon_in_radians)           /* keep longitude in -180 to +180 */
      {
      *lon_in_radians = centralize_ang( tlon);
      if( *lon_in_radians > PI)
         *lon_in_radians -= PI + PI;
      }
   return( rval);
}

/* On pass=0,  we just parse both files containing MPC stations:  one
provided by MPC -- usually 'ObsCodes.html' -- and a 'rovers.txt' file
containing some non-standard MPC codes,  to find out how many lines
we have and how much memory they'll consume.  At the end of pass 0,
we allocate that memory.
   On pass=1,  we actually read in and store those lines in the 'rval'
array of strings. */

static inline char **load_mpc_stations( int *n_stations)
{
   char **rval = NULL, *codes = NULL;
   int pass, loop;

   for( pass = 0; pass < 2; pass++)
      {
      size_t buff_loc = 0;

      *n_stations = 0;
      for( loop = 0; loop < 2; loop++)
         {
         FILE *ifile;

         if( loop)
            ifile = fopen_ext( "rovers.txt", "fcrb");
         else
            {
            ifile = fopen_ext( "ObsCodes.html", "crb");
            if( !ifile)
               ifile = fopen_ext( "ObsCodes.htm", "fcrb");
            }
         if( ifile)
            {
            char buff[200];

            while( fgets_trimmed( buff, sizeof( buff), ifile))
               {
               const int planet_idx = extract_mpc_station_data( buff, NULL,
                              NULL, NULL);

               if( planet_idx != -1)
                  {
                  if( rval)
                     {
                     rval[*n_stations] = codes + buff_loc;
                     strcpy( rval[*n_stations], buff);
                     }
                  buff_loc += strlen( buff) + 1;
                  (*n_stations)++;
                  }
               }
            fclose( ifile);
            }
         }
      if( !pass)     /* At end of first pass,  make the buffer: */
         {
         rval = (char **)calloc( (*n_stations + 1) * sizeof( char *)
                                       + buff_loc, 1);
         codes = (char *)( rval + *n_stations + 1);
         }
      }
   return( rval);
}

static int get_asteroid_observer_data( const char *mpc_code, char *buff)
{
   FILE *ifile = fopen_ext( "mu1.txt", "fcrb");
   char tbuff[100];
   int line_no = 0;

   assert( ifile);
   while( fgets_trimmed( tbuff, sizeof( tbuff), ifile))
      if( *tbuff != ';')
         {
         if( atoi( tbuff) == atoi( mpc_code + 3))
            {
            strcpy( buff + 30, tbuff + 19);
            memcpy( buff, mpc_code, 3);
            return( line_no + 100);
            }
         else
            line_no++;
         }
   return( 0);
}

/* For all MPC stations in ObsCodes.html,  the station name starts
in column 31.  If you look at 'rovers.txt',  you'll see that some
station data lines put an ! in column 5,  in which case the station
name starts in column 48,  allowing room for some extra digits of
precision.  (And also allowing,  eventually,  for four-character
MPC codes.)   */

char *mpc_station_name( char *station_data)
{
   return( station_data + (station_data[4] == '!' ? 47 : 30));
}

/* The following function paws through the STATIONS.TXT file (or the
   ObsCodes.html or .htm file),  looking for the observer code in
   question.  When it finds it,  it just copies the entire line into
   buff.  If lon_in_radians and the rho_xxx_phi values are non-NULL,
   they're extracted from the buffer.

   There are a few "supplemental" observers,  mostly satellite observers
   who don't have MPC codes.  These could be handled as roving observers
   (code 247),  but this can be a hassle;  it's better if they have their
   own codes.  These codes are put into 'rovers.txt',  and have designations
   that are the initials of the observer;  that way,  they don't look
   too much like "real,  official" MPC designations.  At present,  there
   are a few codes there for artificial satellite observers.

   Return value:

      -2:  stations.txt,  obscodes.htm,  obscodes.html not found (need
               any one of these)
       Other:  index of planet of MPC station (3=earth,  most common
            case;  0=sun, 1=mercury,  etc.)
*/

static double roving_lon, roving_lat, roving_ht_in_meters;
int n_obs_actually_loaded;

int get_observer_data( const char FAR *mpc_code, char *buff,
              double *lon_in_radians, double *rho_cos_phi, double *rho_sin_phi)
{
   static char *curr_station = NULL;
   static char **station_data = NULL;
   static int n_stations = 0;
   const char *blank_line = "!!!   0.0000 0.000000 0.000000Unknown Station Code";
   int rval = -1;
   size_t i;
   const char *format_string = NULL;
   double lat0 = 0., lon0 = 0., alt0 = 0.;

   if( !mpc_code)    /* freeing up resources */
      {
      free( station_data);
      station_data = NULL;
      curr_station = NULL;
      n_stations = 0;
      xref_designation( NULL);
      return( 0);
      }

   if( !n_stations)
      {
      int sort_column = 0;

      station_data = load_mpc_stations( &n_stations);
      shellsort_r( station_data, n_stations, sizeof( char *),
                     string_compare_for_sort, &sort_column);
      for( i = 1; i < (size_t)n_stations; i++)
         if( !memcmp( station_data[i], station_data[i - 1], 4))
            {              /* duplication found:  use the one from  */
            if( station_data[i][4] == '!')            /* rovers.txt */
               station_data[i - 1] = station_data[i];
            else
               station_data[i] = station_data[i - 1];
            }
      }
   if( lon_in_radians)
      *lon_in_radians = *rho_cos_phi = *rho_sin_phi = 0.;
   else           /* attempting to look up an MPC code from the station name */
      {
      int pass;

      for( pass = 0; pass < 2; pass++)
         for( i = 0; station_data[i]; i++)
            if( (!pass && !memcmp( buff, station_data[i] + 30, strlen( buff)))
                     || (pass && strstr( station_data[i] + 30, buff)))
               {
               strcpy( buff, station_data[i]);
               return( 0);
               }
      return( -1);
      }

   for( i = 0; i < MAX_XXX_CODES; i++)
      if( !strcmp( mpc_code, xxx_codes[i].mpc_code))
         {
         lat0 = xxx_codes[i].lat;
         lon0 = xxx_codes[i].lon;
         alt0 = xxx_codes[i].alt;
         format_string = "%11.5f%9.5f%7.1fTemporary MPC code";
         }

   if( !strcmp( mpc_code, "247"))
      {
      lat0 = roving_lat;
      lon0 = roving_lon;
      alt0 = roving_ht_in_meters;
      format_string = "%11.5f%9.5f%7.1fRoving observer";
      }

#ifdef TRY_THIS_SOME_OTHER_TIME
   if( strchr( "nsew", tolower( mpc_code[0])))
      {
      char sign2, sign1 = tolower( mpc_code[0]);

      debug_printf( "Looks like lat/lon: '%s'\n", mpc_code);
      if( sscanf( mpc_code + 1, "%lf%c%lf %lf", &lat0, &sign2, &lon0, &alt0) >= 3)
         {
         sign2 = tolower( sign2);

         if( sign1 == 'w' || sign1 == 's')
            lat0 = -lat0;
         if( sign2 == 'w' || sign2 == 's')
            lat0 = -lon0;
         if( sign1 == 'w' || sign1 == 'e')
            {                     /* actually gave longitude first; */
            lon0 += lat0;         /* swap lat & lon */
            lat0 = lon0 - lat0;
            lon0 -= lat0;
            }
         format_string = "%11.5%9.5%7.1User-supplied observer";
         }
      debug_printf( "lat %f sign2 '%c' lon %f\n", lon0, sign2, lat0);
      }
#endif

   if( format_string)
      {
      char tbuff[90];

      strcpy( tbuff, mpc_code);
      sprintf( tbuff + 3, format_string, lon0, lat0, alt0);
      if( buff)
         strcpy( buff, tbuff);
      extract_mpc_station_data( tbuff, lon_in_radians,
                                        rho_cos_phi, rho_sin_phi);
      return( 3);
      }

   if( !memcmp( mpc_code, "Ast", 3))
      {
      assert( buff);
      strcpy( buff, blank_line);
      return( get_asteroid_observer_data( mpc_code, buff));
      }

   if( !curr_station || memcmp( curr_station, mpc_code, 3))
      {
      int step, loc = -1, loc1;

      curr_station = NULL;
      for( step = 0x8000; step; step >>= 1)
         if( (loc1 = loc + step) < n_stations)
            {
            const int compare = memcmp( station_data[loc1], mpc_code, 3);

            if( compare <= 0)
               loc = loc1;
            if( !compare)
               curr_station = station_data[loc];
            }
      }
   if( !curr_station)
      {
      debug_printf( "Couldn't find MPC station '%s'\n", mpc_code);
      if( buff)
         {
         strcpy( buff, blank_line);
         memcpy( buff, mpc_code, 3);
         }
      }
   else
      {
      char *tptr = strchr( curr_station, '@');

      rval = (tptr ? atoi( tptr + 1) : 3);
      if( buff)
         strcpy( buff, curr_station);
      extract_mpc_station_data( curr_station, lon_in_radians,
                                        rho_cos_phi, rho_sin_phi);
      }
   return( rval);
}

/* As the function name suggests,  gets the lat/lon of an MPC station.
   Return value is the planet index (3=earth, 0=sun, 1=mercury,  etc.)
   or a negative value if the station doesn't exist,  or if there's no
   latitude/longitude (planet-centric case,  or spacecraft).  */

static int get_observer_data_latlon( const char FAR *mpc_code,
              char *buff, double *lon_in_radians, double *lat_in_radians,
              double *alt_in_meters)
{
   double rho_cos_phi, rho_sin_phi, alt = 0.;
   int rval;

   *lon_in_radians = *lat_in_radians = 0.;
   rval = get_observer_data( mpc_code, buff, lon_in_radians,
                                   &rho_cos_phi, &rho_sin_phi);
   if( rval >= 0 && (rho_cos_phi || rho_sin_phi))
      {       /* Cvt parallax data from AU back into earth-axis units: */
      rho_cos_phi /= EARTH_MAJOR_AXIS_IN_AU;
      rho_sin_phi /= EARTH_MAJOR_AXIS_IN_AU;
      parallax_to_lat_alt( rho_cos_phi, rho_sin_phi,
                                      lat_in_radians, &alt, rval);
      }
   else
      rval = -2;
   if( alt_in_meters)
      *alt_in_meters = alt;
   return( rval);
}

/*   The offset between a satellite observation and the earth or sun    */
/* is stored in a second line,  as described at                         */
/* http://www.minorplanetcenter.net/iau/info/SatelliteObs.html          */
/*    This format allows parallax type '1' in kilometers or type '2'    */
/* in AU.  If the input file contains the line '#relax_xyz',  Find_Orb  */
/* is less picky about where decimal points and +/- signs appear.       */
/* (It used to insist that the field be "filled out" with digits,  but  */
/* at least some records contain lower precision positions.  Thus far,  */
/* the ones I've seen still had enough digits to match the precision of */
/* the instrument,  so I'm letting "short" records slide if they're     */
/* only lacking three or fewer places.)                                 */

#define SATELL_COORD_BAD_SIGN       -1
#define SATELL_COORD_BAD_NUMBER     -2
#define SATELL_COORD_NO_DECIMAL     -3

static bool strict_sat_xyz_format = true;

inline double get_satellite_coordinate( const char *iptr, int *decimal_loc)
{
   char tbuff[12];
   const char sign_byte = *iptr;
   double rval = 0.;

   memcpy( tbuff, iptr, 11);
   tbuff[11] = '\0';
   if( !strict_sat_xyz_format)
      {
      rval = atof( tbuff + 1);
      if( sign_byte == '-')
         rval = -rval;
      *decimal_loc = 0;
      }
   else if( sign_byte != '+' && sign_byte != '-')
      *decimal_loc = -1;         /* signal bad sign */
   else
      {
      char *tptr;
      int n_bytes_read;

      if( sscanf( tbuff + 1, "%lf%n", &rval, &n_bytes_read) != 1
                     || n_bytes_read < 7)
         *decimal_loc = -2;
      else if( (tptr = strchr( tbuff, '.')) == NULL)
         *decimal_loc = -3;
      else
         *decimal_loc = (int)( tptr - tbuff);
      if( sign_byte == '-')
         rval = -rval;
      if( *decimal_loc <= 0)
         debug_printf( "decimal loc %d: '%s', n_bytes_read %d\n", *decimal_loc,
                           tbuff, n_bytes_read);
      }
   return( rval);
}

static int get_satellite_offset( const char *iline, double *xyz)
{
   unsigned i;
   int error_code = 0, decimal_loc;
   const int observation_units = (int)iline[32] - '0';

   for( i = 0; i < 3; i++)    /* in case of error,  use 0 offsets */
      xyz[i] = 0.;
   iline += 34;      /* this is where the offsets start */
   for( i = 0; !error_code && i < 3; i++, iline += 12)
      {
      xyz[i] = get_satellite_coordinate( iline, &decimal_loc);
      if( observation_units == 1)         /* offset given in km */
         {
         xyz[i] /= AU_IN_KM;
         if( strict_sat_xyz_format && decimal_loc != 6)
            error_code = -1;
         }
      else if( observation_units == 2)          /* offset in AU */
         {
         if( strict_sat_xyz_format && decimal_loc != 2)
            error_code = -2;
         }
      else      /* don't know about this sort of offset */
         error_code = -3;
      if( !error_code && xyz[i] == 0.)
         error_code = -4;
      }
   equatorial_to_ecliptic( xyz);
   return( error_code);
}

/* Used in part for sanity checks ("is the observed RA/dec above the
   horizon?  Is the sun _below_ the horizon at that time?")  Either
   alt/az can be NULL if you're only after one alt/az.

   Return value = 0 if successful,  nonzero otherwise.  (For the function
   to work,  the MPC station must be topocentric.  So you won't get alt/az
   values for a geocentric/planetocentric location,  nor from spacecraft.)

*/

static int get_obs_alt_azzes( const OBSERVE FAR *obs, DPT *sun_alt_az,
                                               DPT *object_alt_az)
{
   DPT latlon;
   int i;
   int rval = (get_observer_data_latlon( obs->mpc_code, NULL,
                                    &latlon.x, &latlon.y, NULL) != 3);

   if( !rval)
      {
      DPT ra_dec;
      const double utc = obs->jd - td_minus_utc( obs->jd) / seconds_per_day;

      for( i = 0; i < 2; i++)
         {
         DPT *alt_az = (i ? object_alt_az : sun_alt_az);

         if( alt_az)
            {
            if( !i)      /* compute solar alt/az */
               {
               double equat[3];

               memcpy( equat, obs->obs_posn, 3 * sizeof( double));
               ecliptic_to_equatorial( equat);
               ra_dec.x = atan2( equat[1], -equat[0]);
               ra_dec.y = -asin( equat[2] / vector3_length( equat));
               }
            else
               {
               if( obs->note2 == 'R')
                  {
                  ra_dec.x = -obs->computed_ra;
                  ra_dec.y = obs->computed_dec;
                  }
               else
                  {
                  ra_dec.x = -obs->ra;
                  ra_dec.y = obs->dec;
                  }
               }
            full_ra_dec_to_alt_az( &ra_dec, alt_az, NULL, &latlon, utc, NULL);
            }
         }
      }
   else if( obs->second_line && obs->second_line[14] == 's')
      {
      double xyz[3], len, cos_sun = 0., cos_obj = 0.;
      double observer_r = vector3_length( obs->obs_posn);

      get_satellite_offset( obs->second_line, xyz);
      len = vector3_length( xyz);
      for( i = 0; i < 3; i++)
         {
         cos_sun += xyz[i] * obs->obs_posn[i];
         cos_obj += xyz[i] * obs->vect[i];
         }
      object_alt_az->y = asin( cos_obj / len);
      sun_alt_az->y = asin( -cos_sun / (observer_r * len));
      object_alt_az->x = sun_alt_az->x = -99.;  /* flag azimuths as meaningless */
      rval = 0;
      }
   sun_alt_az->x *= 180. / PI;
   sun_alt_az->x += 180.;
   sun_alt_az->y *= 180. / PI;
   object_alt_az->x *= 180. / PI;
   object_alt_az->x += 180.;
   object_alt_az->y *= 180. / PI;

   return( rval);
}

/* "Mutant hex" is frequently used by MPC.  It uses the usual hex digits
   0123456789ABCDEF for numbers 0 to 15,  followed by G...Z for 16...35
   and a...z for 36...61.  */

static int mutant_hex_char_to_int( const char c)
{
   int rval;

   if( c >= '0' && c <= '9')
      rval = (int)c - '0';
   else if( c >= 'A' && c <= 'Z')
      rval = (int)c - 'A' + 10;
   else if( c >= 'a' && c <= 'z')
      rval = (int)c - 'a' + 36;
   else
      rval = -1;
   return( rval);
}

char int_to_mutant_hex_char( const int ival)
{
   int rval;

   if( ival < 0 || ival > 61)
      rval = '\0';
   else if( ival < 10)
      rval = '0';
   else if( ival < 36)
      rval = 'A' - 10;
   else
      rval = 'a' - 36;
   return( rval ? (char)( rval + ival) : '\0');
}

/* This will unpack a packed designation such as 'K04J42X' into '2004 JX42'. */

static int unpack_provisional_packed_desig( char *obuff, const char *ibuff)
{
   int rval = 0;

   if( *ibuff >= 'G' && *ibuff <= 'K' && isdigit( ibuff[1])
            && isdigit( ibuff[2]) && isupper( ibuff[3])
            && isdigit( ibuff[5]) && isalnum( ibuff[6]))
      {
      int output_no = mutant_hex_char_to_int( ibuff[4]);
//    char *orig_obuff = obuff;

      if( output_no == -1)
         rval = -2;

      *obuff++ = ((*ibuff >= 'K') ? '2' : '1');          /* millennium */
      *obuff++ = (char)( '0' + ((*ibuff - 'A') % 10));   /* century   */
      *obuff++ = ibuff[1];                               /* decade   */
      *obuff++ = ibuff[2];                               /* year    */
      *obuff++ = ' ';
      *obuff++ = ibuff[3];       /* half-month designator */
      if( isupper( ibuff[6]))    /* asteroid second letter */
         *obuff++ = ibuff[6];
      output_no = output_no * 10 + ibuff[5] - '0';
      if( !output_no)
         *obuff = '\0';
      else
         sprintf( obuff, "%d", output_no);
      if( islower( ibuff[6]))    /* comet fragment letter */
         sprintf( obuff + strlen( obuff), "%c", ibuff[6]);
//    debug_printf( "Resulting output: '%s'\n", orig_obuff);
      }
   else
      {
      *obuff = '\0';
      rval = -1;
      }
   return( rval);
}

char **load_file_into_memory( const char *filename, size_t *n_lines)
{
   FILE *ifile = fopen_ext( filename, "crb");
   char **rval = NULL;

   if( ifile)
      {
      size_t filesize = 0, lines_read = 0;
      char buff[400];

      while( fgets_trimmed( buff, sizeof( buff), ifile))
         {
         lines_read++;
         filesize += strlen( buff) + 1;
         }
      rval = (char **)malloc( filesize + (lines_read + 1) *  sizeof( char *));
      if( rval &&
             (rval[0] = (char *)( rval + lines_read + 1)) != NULL)
         {
         fseek( ifile, 0L, SEEK_SET);
         lines_read = 0;
         while( fgets_trimmed( buff, sizeof( buff), ifile))
            {
            strcpy( rval[lines_read], buff);
            rval[lines_read + 1] = rval[lines_read] + strlen( buff) + 1;
            lines_read++;
            }
         rval[lines_read] = NULL;
         }
      fclose( ifile);
      if( n_lines)
         *n_lines = lines_read;
      }
   return( rval);
}

static void try_adding_comet_name( const char *search, char *name)
{
   FILE *ifile = fopen_ext( "ELEMENTS.COMET", "crb");
   bool found_name = false;

   *name = '\0';
   if( ifile)
      {
      char buff[200];
      const size_t slen = strlen( search);

      while( !found_name && fgets_trimmed( buff, sizeof( buff), ifile))
         if( !memcmp( buff, search, slen))
            {
            found_name = true;
            if( buff[slen] != ' ')     /* got a real name */
               {
               buff[44] = '\0';
               remove_trailing_cr_lf( buff);
               *name++ = ' ';
               strcpy( name, buff + slen);
               }
            }
      fclose( ifile);
      }
}

/* An object with a name such as 1989-013A is probably an artsat,  and
probably has a NORAD designation within 'all_tle.txt' or similar file.
(Thus far,  only 'all_tle.txt' is searched.) The following will append
that designation if it's found. */

static bool try_artsat_xdesig( char *name, const char *filename)
{
   FILE *ifile = fopen_ext( "all_tle.txt", "crb");
   bool found_a_match = false;

   if( ifile)
      {
      char tbuff[100], xdesig[20];

      memset( xdesig, ' ', 10);
      xdesig[0] = name[2];      /* decade */
      xdesig[1] = name[3];      /* year */
      memcpy( xdesig + 2, name + 5, strlen( name + 5));
      while( !found_a_match && fgets( tbuff, sizeof( tbuff), ifile))
         if( *tbuff == '1' && !memcmp( tbuff + 9, xdesig, 9))
            {
            found_a_match = true;
            tbuff[7] = '\0';
            snprintf_append( name, 30, " = NORAD %s", tbuff + 2);
            }
      fclose( ifile);
      }
   return( found_a_match);
}

/* In an MPC astrometric report line,  the name can be stored in assorted
   highly scrambled and non-intuitive ways.  Those ways _are_ documented
   on the MPC Web site,  which is probably the best place to look to
   understand why the following bizarre code does what it does:

   http://www.minorplanetcenter.org/iau/info/PackedDes.html

   Return value is zero if it's an asteroid designation,  1 if it's a comet,
   -1 if it wasn't puzzled out.        */

int get_object_name( char *obuff, const char *packed_desig)
{
   int rval = -1;
   size_t i, gap;
   static size_t n_lines;
   static char **extra_names = NULL;
   char provisional_desig[40], xdesig[40];

   if( !obuff)       /* flag to free up internal memory */
      {
      if( extra_names)
         free( extra_names);
      extra_names = NULL;
      return( 0);
      }

   if( !extra_names)          /* see 'odd_name.txt' for comments on this */
      {
      int sort_column = 0;

      extra_names = load_file_into_memory( "odd_name.txt", &n_lines);
      shellsort_r( extra_names, n_lines, sizeof( char *),
                     string_compare_for_sort, &sort_column);
      }
   strcpy( xdesig, packed_desig);
   xref_designation( xdesig);

   for( i = 0, gap = 0x8000; gap; gap >>= 1)
      if( i + gap < n_lines)
         {
         const int compare = memcmp( xdesig, extra_names[i + gap], 12);

         if( compare >= 0)
            i += gap;
         if( !compare)
            {
            strcpy( obuff, extra_names[i] + 13);
            return( 0);
            }
         }

   if( xdesig[4] == 'S')   /* Possible natural satellite */
      {
      if( strchr( "MVEJSUNP", *xdesig) && isdigit( xdesig[1])
                  && isdigit( xdesig[2]) && isdigit( xdesig[3]) &&
                  !memcmp( xdesig + 5, "       ", 7))
         {
         const char *planet_names[8] = { "Venus", "Earth", "Mars", "Jupiter",
                  "Saturn", "Uranus", "Neptune", "Pluto" };
         const char *roman_digits[10] = { "", "I", "II", "III", "IV",
                     "V", "VI", "VII", "VIII", "IX" };
         const char *roman_tens[10] = { "", "X", "XX", "XXX", "XL",
                     "L", "LX", "LXX", "LXXX", "XL" };
         const char *roman_hundreds[10] = { "", "C", "CC", "CCC", "CD",
                     "D", "DC", "DCC", "DCCC", "CD" };
         const int obj_number = atoi( xdesig + 1);

         for( i = 0; i < 8; i++)
            if( planet_names[i][0] == *xdesig)
               {
               strcpy( obuff, planet_names[i]);
               strcat( obuff, " ");
               }
         if( obj_number / 100)
            strcat( obuff, roman_hundreds[obj_number / 100]);
         if( (obj_number / 10) % 10)
            strcat( obuff, roman_tens[(obj_number / 10) % 10]);
         if( obj_number % 10)
            strcat( obuff, roman_digits[obj_number % 10]);
         rval = 0;
         }
      else if( strchr( "MVEJSUNP", xdesig[8]) && isdigit( xdesig[6])
               && isdigit( xdesig[7]) && isdigit( xdesig[9])
               && isdigit( xdesig[10]) && xdesig[11] == '0'
               && !memcmp( xdesig + 1, "   ", 3) && xdesig[5] >= 'H'
               && xdesig[5] <= 'Z')
         {
         sprintf( obuff, "S/%d%c%c", 20 + xdesig[5] - 'K',
                                  xdesig[6], xdesig[7]);
         obuff[6] = ' ';
         obuff[7] = xdesig[8];       /* planet identifier */
         obuff[8] = ' ';
         if( xdesig[9] > '0')     /* double-digit ID (unlikely, */
            {                           /* but it _can_ happen)       */
            obuff[9] = xdesig[9];
            obuff[10] = xdesig[10];
            obuff[11] = '\0';
            }
         else                 /* more usual single-digit case: */
            {
            obuff[9] = xdesig[10];
            obuff[10] = '\0';
            }
         rval = 0;
         }
      }
   unpack_provisional_packed_desig( provisional_desig, xdesig + 5);
            /* For numbered asteroids or comets,  we require either that
            columns 6-12 be blank (usually the case) or have a valid
            provisional designation.  And,  of course,  the packed desig must
            start with an alphanumeric and be followed by three digits. */
   if( isalnum( xdesig[0]) && isdigit( xdesig[1]) &&  /* possible numbered */
       isdigit( xdesig[2]) && isdigit( xdesig[3]) &&  /* asteroid or comet */
               (*provisional_desig || !memcmp( xdesig + 5, "       ", 5)))
      {
      if( isdigit( xdesig[4]))
         {                                /* it's a numbered asteroid */
         int number = mutant_hex_char_to_int( *xdesig);

         if( number >= 0)
            {
            number = number * 10000L + atol( xdesig + 1);
            sprintf( obuff, "(%d)", number);
            rval = 0;
                /* Desig may be,  e.g., "U4330K06SL7X" : both the number */
                /* and the provisional desig.  We'd like to decipher this */
                /* as (for the example) "304330 = 2006 SX217".            */
            if( *provisional_desig)
               sprintf( obuff + strlen( obuff), " = %s", provisional_desig);
            }
         }
      else if( strchr( "PCDXA", xdesig[4]))   /* it's a numbered comet */
         {
         const char extra_suffix_char = xdesig[10];
         const char suffix_char = xdesig[11];
         char tbuff[20];

         *obuff++ = xdesig[4];
         *obuff++ = '/';
         i = 0;
         while( xdesig[i] == '0')         /* skip leading zeroes */
            i++;
         while( xdesig[i] >= '0' && xdesig[i] <= '9')   /* read digits... */
            *obuff++ = xdesig[i++];
                  /* possibly _two_ suffix letters... so far,  only the  */
                  /* components of 73P/Schwassmann-Wachmann 3 have this: */
         if( extra_suffix_char >= 'a' && extra_suffix_char <= 'z')
            *obuff++ = extra_suffix_char;
         if( suffix_char >= 'a' && suffix_char <= 'z')
            *obuff++ = suffix_char;
         sprintf( tbuff, "%3d%c/", atoi( xdesig), xdesig[4]);
         try_adding_comet_name( tbuff, obuff);
         rval = 1;
         }
      }

   if( rval == -1 && !memcmp( xdesig, "    ", 4)
            && strchr( " PCDXA", xdesig[4]))
      {
      if( xdesig[4] != ' ')
         {
         *obuff++ = xdesig[4];
         *obuff++ = '/';
         }
      if( !unpack_provisional_packed_desig( obuff, xdesig + 5))
         rval = (xdesig[4] != ' ');
      if( rval != -1 && xdesig[4] != ' ')
         {
         char tbuff[40];

         sprintf( tbuff, "   %s ", obuff - 2);
         try_adding_comet_name( tbuff, obuff + strlen( obuff));
         }
      }

#ifdef POSSIBLY_USEFUL_LATER
               /* Look for CSS/SSS/Mt. Lemmon/LAB desigs.  These should be  */
               /* a hex digit,  plus a half-month character,  plus five hex */
               /* more hex digits,  all uppercase.  The first of these is   */
               /* 0-3 for an object found at CSS,  4-7 for SSS, 8-B for Mt. */
               /* Lemmon,  C-F for LAB (which apparently is rare).          */
   if( rval == -1 && !memcmp( xdesig, "     ", 5))
      {
      int is_valid = 1;

      if( xdesig[6] < 'A' || xdesig[6] > 'Z')
         is_valid = 0;           /* not a valid half-month designation */
      for( i = 5; i < 12 && is_valid; i++)
         if( i != 6)       /* skip the half-month designation */
            if( !isdigit( xdesig[i]))
               if( xdesig[i] < 'A' || xdesig[i] > 'F')
                  is_valid = 0;
      if( is_valid)
         {
         const char *suffixes[4] = { " (CSS)", " (SSS)", " (MtL)", " (LAB)" };
         const int stn = mutant_hex_char_to_int( xdesig[7]) / 4;

         memcpy( obuff, xdesig + 5, 7);
         strcpy( obuff + 7, suffixes[stn]);
         rval = 0;
         }
      }
#endif

   if( rval == -1)    /* store the name "as is",  assuming no encoding */
      {               /* (except skip leading spaces) */
      for( i = 0; i < 12 && xdesig[i] == ' '; i++)
         ;
      memcpy( obuff, xdesig + i, 12 - i);
      obuff[12 - i] = '\0';
      remove_trailing_cr_lf( obuff);      /* ephem0.cpp */
      if( atoi( obuff) > 1956 && atoi( obuff) < 2100)    /* possible artsat */
         if( obuff[4] == '-' && atoi( obuff + 5))
            try_artsat_xdesig( obuff, "all_tle.txt");
      }
   return( rval);
}

void set_obs_vect( OBSERVE FAR *obs)
{
   obs->vect[0] = cos( obs->ra) * cos( obs->dec);
   obs->vect[1] = sin( obs->ra) * cos( obs->dec);
   obs->vect[2] = sin( obs->dec);
   equatorial_to_ecliptic( obs->vect);
}

/* Given a planet number (Sun=0, Mercury=1, ... Pluto=9) and a JD,
   the following code computes that planet's J2000 equatorial coordinates
   using the PS1996 theory (or,  for the moon,  ELP82). */

const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);

int earth_lunar_posn( const double jd, double FAR *earth_loc, double FAR *lunar_loc)
{
   double t_earth[3], t_lunar[3];
   int i;

   planet_posn( 3, jd, t_earth);
   planet_posn( 10, jd, t_lunar);
                     /* earth_loc is the E-M barycenter,  and lunar_loc */
                     /* is geocentric.  Modify earth_loc to be the earth */
                     /* geocenter loc,  and lunar_loc to be heliocentric: */
   for( i = 0; i < 3; i++)
      {
      const double earth_moon_barycenter_factor = 82.300679;

      t_earth[i] -= t_lunar[i] / earth_moon_barycenter_factor;
      t_lunar[i] += t_earth[i];
      if( earth_loc)
         earth_loc[i] = t_earth[i];
      if( lunar_loc)
         lunar_loc[i] = t_lunar[i];
      }
   return( 0);
}

/* Input time is a JD in UT.  Output offset is in equatorial    */
/* J2000, in AU.    Declared 'inline' because it's used exactly */
/* once,  in compute_observer_loc.                              */

static inline int compute_topocentric_offset( const double ut,
               const int planet_no,
               const double rho_cos_phi,
               const double rho_sin_phi, const double lon,
               double *offset, double *vel)
{
   double precess_matrix[9];
   int i;

   calc_planet_orientation( planet_no, 0, ut, precess_matrix);
   spin_matrix( precess_matrix, precess_matrix + 3, lon);
   for( i = 0; i < 3; i++)
      {
      const double omega = 2. * PI * 360.9856235 / 360.;
                   /* earth's rotation rate,  in radians/day */

      if( offset)
         offset[i] = (rho_cos_phi * precess_matrix[i]
                    + rho_sin_phi * precess_matrix[i + 6]);
      if( vel)
         vel[i] = -rho_cos_phi * precess_matrix[i + 3] * omega;
      }
   return( 0);
}


int compute_observer_loc( const double jde, const int planet_no,
               const double rho_cos_phi,
               const double rho_sin_phi, const double lon, double FAR *offset)
{
   if( planet_no != 3 && planet_no != 10)
      planet_posn( planet_no >= 0 ? planet_no : 12, jde, offset);
               /* planet_no == -1 means SS Barycenter */
   else
      earth_lunar_posn( jde, (planet_no == 3) ? offset : NULL,
                             (planet_no == 3) ? NULL : offset);

   if( rho_sin_phi || rho_cos_phi)
      {
      const double ut = jde - td_minus_ut( jde) / seconds_per_day;
      double geo_offset[3];
      int i;

      compute_topocentric_offset( ut, planet_no, rho_cos_phi, rho_sin_phi,
                                        lon, geo_offset, NULL);
      equatorial_to_ecliptic( geo_offset);
      for( i = 0; i < 3; i++)
         offset[i] += geo_offset[i];
      }
   return( 0);
}

int compute_observer_vel( const double jde, const int planet_no,
               const double rho_cos_phi,
               const double rho_sin_phi, const double lon, double FAR *vel)
{
   double loc1[3], loc2[3];
   const double delta_t = 10. / minutes_per_day;
   int i;

   compute_observer_loc( jde + delta_t, planet_no, 0., 0., 0., loc2);
   compute_observer_loc( jde - delta_t, planet_no, 0., 0., 0., loc1);
   for( i = 0; i < 3; i++)
      vel[i] = (loc2[i] - loc1[i]) / (2. * delta_t);
   if( rho_sin_phi || rho_cos_phi)
      {
      const double ut = jde - td_minus_ut( jde) / seconds_per_day;
      double geo_vel_offset[3];

      compute_topocentric_offset( ut, planet_no, rho_cos_phi, rho_sin_phi,
                                        lon, NULL, geo_vel_offset);
      equatorial_to_ecliptic( geo_vel_offset);
      for( i = 0; i < 3; i++)
         vel[i] += geo_vel_offset[i];
      }
   return( 0);
}

/* parse_observation( ) takes an MPC astrometric observation of the usual
   80-character variety,  and extracts all relevant data from it and puts
   it into the 'obs' structure.  It also computes the ecliptic J2000
   coordinates of the observer at the time of the observation,  using the
   PS1996 and ELP82 theories,  or a JPL ephemeris.  */

static double input_coordinate_epoch = 2000.;
static double override_time = 0.;
      /* 'override_time' allows one to set the observation time with the */
      /* #time keyword (see below).  If that's done,  the time from the  */
      /* next observation will be ignored.  */
static double observation_time_offset = 0.;
      /* 'time_offset' allows you to add,  or subtract,  a certain number */
      /* of days from the following observations.  This was added to support */
      /* meteor observations,  where one station might have its clock in */
      /* error by a second or more relative to another. */

static void comment_observation( OBSERVE FAR *obs, const char *comment)
{
   size_t i;

   for( i = 0; comment[i] && i < 9; i++)
      obs->columns_57_to_65[i] = comment[i];
}

void set_up_observation( OBSERVE FAR *obs)
{
   double rho_cos_phi = 0.;
   double rho_sin_phi = 0.;
   double lon = 0.;
   char tbuff[300];
   int observer_planet = get_observer_data( obs->mpc_code, tbuff,
                           &lon, &rho_cos_phi, &rho_sin_phi);

   if( observer_planet < 0)
      {
      static unsigned n_unfound = 0;
      static char unfound[10][4];
      unsigned i;

      obs->is_included = 0;
      obs->discovery_asterisk = '!';
      for( i = 0; FMEMCMP( unfound[i], obs->mpc_code, 3) && i < n_unfound; i++)
         ;
      if( i == n_unfound && n_unfound < 10)     /* got a new one! */
         {
         FMEMCPY( unfound[n_unfound++], obs->mpc_code, 3);
         sprintf( tbuff, "Didn't find observer %s\n", obs->mpc_code);
         strcat( tbuff, "Observation(s) will be excluded and treated as\n");
         strcat( tbuff, "geocentric.  You can fix this by downloading the\n");
         strcat( tbuff, "current list of MPC stations at\n\n");
         strcat( tbuff, "http://www.minorplanetcenter.net/iau/lists/ObsCodes.html\n\n");
         strcat( tbuff, "and saving it to the folder in which Find_Orb runs.\n");
         generic_message_box( tbuff, "o");
         comment_observation( obs, "? NoCode");
         }
      observer_planet = 3;    /* default to geocentric */
      }
   else if( obs->note2 != 'S' && strcmp( obs->mpc_code, "247")
               && !memcmp( tbuff + 3, "                          ", 26))
      {
      obs->is_included = 0;
      obs->flags |= OBS_NO_OFFSET;
      comment_observation( obs, "? offset");
      }
   compute_observer_loc( obs->jd, observer_planet,
               rho_cos_phi, rho_sin_phi, lon, obs->obs_posn);
   compute_observer_vel( obs->jd, observer_planet,
               rho_cos_phi, rho_sin_phi, lon, obs->obs_vel);
   set_obs_vect( obs);
}

static bool is_valid_mpc_code( const char *mpc_code)
{
   int i;

   for( i = 0; i < 3; i++)
      if( mpc_code[i] <= ' ' || mpc_code[i] > 'z')
         return( false);
   return( true);
}

static inline int get_two_digits( const char *iptr)
{
   return( (int)iptr[0] * 10 + (int)iptr[1] - (int)'0' * 11);
}

   /* By default,  Find_Orb will only handle arcs up to 200 years */
   /* long.  If the arc is longer than that,  observations will be */
   /* dropped to get an arc that fits.  The max arc length can be */
   /* adjusted in 'environ.dat'.                                  */
double maximum_observation_span = 200.;

double minimum_observation_year = -1e+9;   /* set in console Find_Orb's */
double maximum_observation_year =  1e+9;   /* command line      */

#define YEAR_TO_JD( year) (J2000 + (year - 2000.) * 365.25)
#define MINIMUM_OBSERVATION_JD YEAR_TO_JD( minimum_observation_year)
#define MAXIMUM_OBSERVATION_JD YEAR_TO_JD( maximum_observation_year)

/* The date/time in an 80-column MPC report line is stored in columns 16-32.
   MPC expects this to be in the form YYYY MM DD.dddddd,  with the date usually
   given to five digits = .864 second precision;  the sixth digit should only
   be given if you're really sure the time is good to that level.  In some
   cases (mostly older,  poorly timed observations),  fewer digits are given.
   A few ancient observations are just given to the nearest day.

   Find_Orb also supports some NON-STANDARD,  FIND_ORB ONLY (i.e.,  MPC will
   reject anything you do in these formats).  These formats support greater
   precision and spare you the need to convert dates from HH:MM:SS.sss format
   to decimal-day format,  and/or the need to convert from JD or MJD format
   to YYYY MM SS format.  To demonstrate,  the following lines give the same
   date/time in the various formats,  give or take a little rounding error :

               2013 02 13.141593     (MPC's expected format)
               2456336.641592653     (Julian Day format)
               M056336.141592653     (MJD = Modified Julian Day)
               K130213.141592653     (CYYMMDD.dddddd)
               K130213:032353605     (CYYMMDD HH:MM:SS.sss)

   The last two use the MPC "mutant hex" convention of K for 21st century
   years and J for twentieth-century years,  followed by a two-digit year.
   After the two-digit month and two-digit day of month,  one can have
   a decimal point and decimal day _or_ a colon and HHMMSS,  and optionally,
   up to millisecond precision.  (Though you can and should supply fewer
   digits if -- as will almost always be the case -- your timing is not
   really that exact.)  The last one is a little confusing,  given the
   placement of the ':':  it actually means "2013 02 13 03:23:53.605".

   Note that these non-standard formats allow precision up to 10^-9 day
   (86.4 microseconds) or,  for the last format,  one millisecond.  This
   precision should be good enough for anybody.  Indeed,  giving roundoff
   issues,  I doubt Find_Orb could really use extra digits anyway.

   For each of these formats,  a very little bit of format checking is
   done (make sure digits are in certain key places,  and that the full
   line is exactly 80 bytes).  Some malformed records _can_ slip past!

         49    M056336.641592653     (MJD, 10^-9 day)
         48    M056336.64159265      (MJD, 10^-8 day)
         47    M056336.6415926       (MJD, 10^-7 day)
         46    M056336.641592        (MJD, 10^-6 day)
         45    M056336.64159         (MJD, 10^-5 day)
         44    M056336.6415          (MJD, 10^-4 day)
         43    M056336.641           (MJD,  .001 day)
         42    M056336.64            (MJD,   .01 day)
         41    M056336.6             (MJD,    .1 day... not supported)
         40    M056336.              (MJD,     1 day... not supported)

         39    K130213.141592653     (High-prec, 10^-9 day)
         38    K130213.14159265      (High-prec, 10^-8 day)
         37    K130213.1415926       (High-prec, 10^-7 day)
         36    K130213.141592        (High-prec, 10^-6 day)
         35    K130213.14159         (High-prec, 10^-5 day)
         34    K130213.1415          (High-prec, 10^-4 day)
         33    K130213.141           (High-prec,  .001 day)
         32    K130213.14            (High-prec,   .01 day)
         31    K130213.1             (High-prec,    .1 day... not supported)
         30    K130213.              (High-prec,     1 day... not supported)

         23    K130213:032353605     (CYYMMDD HH:MM:SS.sss)
         22    K130213:03235360      (CYYMMDD HH:MM:SS.ss)
         21    K130213:0323536       (CYYMMDD HH:MM:SS.s)
         20    K130213:032353        (CYYMMDD HH:MM:SS)

         19    2456336.641592653     (Julian Day, 10^-9 day)
         18    2456336.64159265      (Julian Day, 10^-8 day)
         17    2456336.6415926       (Julian Day, 10^-7 day)
         16    2456336.641592        (Julian Day, 10^-6 day)
         15    2456336.64159         (Julian Day, 10^-5 day)
         14    2456336.6415          (Julian Day, 10^-4 day)
         13    2456336.641           (Julian Day, 10^-3 day)
         12    2456336.64            (Julian Day, 10^-2 day... not supported)
         11    2456336.6             (Julian Day, 10^-1 day... not supported)

          6    2013 02 13.141593     (MPC's expected format, 10^-6 day)
          5    2013 02 13.14159      (MPC's expected format, 10^-5 day)
          4    2013 02 13.1415       (MPC's expected format, 10^-4 day)
          3    2013 02 13.141        (MPC's expected format, 10^-3 day)
          2    2013 02 13.14         (MPC's expected format, 10^-2 day)
          1    2013 02 13.1          (MPC's expected format, 10^-1 day)
          0    2013 02 13.           (MPC's expected format, 10^-0 day) */

static double extract_date_from_mpc_report( const char *buff, unsigned *format)
{
   double rval = 0.;
   int year = 0, month = 0;
   size_t start_of_decimals = 0;
   unsigned format_found = 0;
   char tbuff[18];
   const size_t len = strlen( buff);
   unsigned i, bit, digits_mask = 0;

   if( len != 80)             /* check for correct length */
      return( 0.);
   if( buff[12] != ' ' && buff[12] != '*')
      return( 0.);
   if( !is_valid_mpc_code( buff + 77))
      return( 0.);
   memcpy( tbuff, buff + 15, 17);
   for( i = 0, bit = 1; i < 17; i++, bit <<= 1)
      if( isdigit( tbuff[i]))
         digits_mask |= bit;
   tbuff[17] = '\0';
   if( tbuff[4] == ' ')
      {                       /* standard YYYY MM DD.dddddd format */
      if( (digits_mask & 0x3ff) == 0x36f    /* i.e.,  'dddd dd dd' */
                            && tbuff[7] == ' ' && tbuff[10] == '.')
         {
         int divisor = 1;

         year = atoi( tbuff);
         month = atoi( tbuff + 5);
//       rval = atof( tbuff + 8);
                     /* atof( ) is a little slow,  so we use a little more */
         for( i = 11; i < 17 && tbuff[i] != ' '; i++)  /* code in exchange */
            divisor *= 10;                             /* for better speed */
         rval = (double)atoi( tbuff + 8) +
                           (double)atoi( tbuff + 11) / (double)divisor;
         format_found = 0;
         start_of_decimals = 11;
         }
      }
   else if( *tbuff >= 'H' && *tbuff <= 'K')  /* 18th through 21st century */
      {                                          /* CYYMMDD format */
      if( (tbuff[7] == '.' || tbuff[7] == ':')
               && (digits_mask & 0x3ff) == 0x37e)  /* i.e, 'Zdddddd.dd' */
         {
         year = (*tbuff - 'J') * 100 + 1900 +
                    get_two_digits( tbuff + 1);
         month = get_two_digits( tbuff + 3);
         rval = atof( tbuff + 5);
         if( tbuff[7] == ':')
            {
            rval += (double)get_two_digits( tbuff + 8) / hours_per_day
               + (double)get_two_digits( tbuff + 10) / minutes_per_day
               + (double)get_two_digits( tbuff + 12) / seconds_per_day;
            tbuff[13] = '.';
            rval += atof( tbuff + 13) / seconds_per_day;
            format_found = 20;      /* formats 20-23;  see above */
            start_of_decimals = 14;
            }
         else     /* decimal formats 32-40 */
            {
            format_found = 30;
            start_of_decimals = 8;
            }
         }
      }
   else if( tbuff[7] == '.')        /* MJD or JD format */
      {
      if( (digits_mask & 0x3fe) == 0x37e)   /* i.e., 'zdddddd.dd' */
         {
         if( *tbuff == 'M')    /* MJD */
            {
            format_found = 40;
            rval = 2400000.5 + atof( tbuff + 1);
            }
         else
            {
            format_found = 10;
            rval = atof( tbuff);      /* plain ol' JD */
            }
         start_of_decimals = 8;
         }
      }
   if( format)
      {
      if( start_of_decimals)
         while( isdigit( tbuff[start_of_decimals++]))
            format_found++;
      *format = format_found;
      }

   if( month >= 1 && month <= 12 && rval > 0. && rval < 99.)
      rval += (double)dmy_to_day( 0, month, year,
                                    CALENDAR_JULIAN_GREGORIAN) - .5;

   if( rval < MINIMUM_OBSERVATION_JD || rval > MAXIMUM_OBSERVATION_JD)
      rval = 0.;
             /* Radar obs are always given to the nearest UTC second. So  */
             /* some rounding is usually required with MPC microday data. */
   if( rval && (buff[14] == 'R' || buff[14] == 'r'))
      {
      const double time_of_day = rval - floor( rval);
      const double resolution = 1. / seconds_per_day;
      const double half = .5 / seconds_per_day;

      rval += half - fmod( time_of_day + half, resolution);
      }
   return( rval);
}

/* Some historical observations are provided in apparent coordinates of date.
When that happens,  they have to be adjusted for both aberration of light
and nutation and precession.  The following removes the aberration. */

static void adjust_for_aberration( OBSERVE FAR *obs)
{
   size_t i;

   set_up_observation( obs);
   for( i = 0; i < 3; i++)
      obs->vect[i] -= obs->obs_vel[i] / AU_PER_DAY;
   ecliptic_to_equatorial( obs->vect);
   obs->ra  = atan2( obs->vect[1], obs->vect[0]);
   obs->dec = asin( obs->vect[2] / vector3_length( obs->vect));
}


int apply_debiasing = 0;

int find_fcct_biases( const double ra, const double dec, const char catalog,
                 const double jd, double *bias_ra, double *bias_dec);

static int parse_observation( OBSERVE FAR *obs, const char *buff)
{
   unsigned time_format;
   double utc = extract_date_from_mpc_report( buff, &time_format);
   char tbuff[300];
   const bool is_radar_obs = (buff[14] == 'R' || buff[14] == 'r');
   static bool fcct_error_message_shown = false;
   const double saved_ra_bias = obs->ra_bias;
   const double saved_dec_bias = obs->dec_bias;
   double coord_epoch = input_coordinate_epoch;

   if( !utc)
      return( -1);
   assert( obs);
   assert( buff);
   memset( obs, 0, sizeof( OBSERVE));
   memcpy( obs->packed_id, buff, 12);
   obs->ra_bias = saved_ra_bias;
   obs->dec_bias = saved_dec_bias;
   obs->packed_id[12] = '\0';
   if( get_object_name( tbuff, obs->packed_id) == 1)
      obs->flags = OBS_IS_COMET;
   else
      obs->flags = 0;
   if( override_time)
      {
      utc = override_time;
      override_time = 0.;
      }
   utc += observation_time_offset;
   obs->jd = utc + td_minus_utc( utc) / seconds_per_day;

   obs->mag_band = buff[70];
   obs->mag_band2 = buff[71];
   obs->discovery_asterisk = buff[12];
   obs->note1 = buff[13];
   obs->note2 = buff[14];
   if( is_radar_obs || buff[14] == 's' || buff[14] == 'v')
      {     /* radar data and "second lines" have no RA/dec: */
      obs->ra = obs->dec = 0.;
      obs->flags |= OBS_DONT_USE;
      }
   else
      {
      obs->ra  = get_ra_dec( buff + 32, &obs->ra_precision, &obs->posn_sigma_1) * (PI / 12.);
      obs->posn_sigma_1 *= 15.;     /* cvt minutes/seconds to arcmin/arcsec */
      if( obs->ra_precision == BAD_RA_DEC_FMT)
         return( -2);
      obs->dec = get_ra_dec( buff + 44, &obs->dec_precision, &obs->posn_sigma_2) * (PI / 180.);
      if( obs->dec_precision == BAD_RA_DEC_FMT)
         return( -3);
      }

   obs->time_precision = time_format;
   if( is_radar_obs)                  /* Radar obs are treated as */
      obs->time_sigma = 0.;           /* being 'perfectly' timed  */
   else
      {
      obs->time_sigma = pow( .1, (double)( time_format % 10));
      if( time_format / 10 == 2)    /* CYYMMDD HH:MM:SS.ss.. formats */
         obs->time_sigma /= seconds_per_day;
      }

   obs->mag_precision = 2;         /* start out assuming mag to .01 mag */
   while( obs->mag_precision && buff[67 + obs->mag_precision] == ' ')
      obs->mag_precision--;
   obs->mag_sigma = pow( .1, (double)obs->mag_precision);
   if( buff[67] == ' ' && buff[66] >= '0')     /* mag given to integer value */
      obs->mag_precision = -1;

   obs->is_included = (buff[64] != 'x' && buff[12] != '-');
   FMEMCPY( obs->mpc_code, buff + 77, 3);
   obs->mpc_code[3] = '\0';
   FMEMCPY( obs->columns_57_to_65, buff + 56, 9);
   obs->columns_57_to_65[9] = '\0';
   if( !strcmp( obs->columns_57_to_65, "Apparent "))
      coord_epoch = -1.;
   else if( !strcmp( obs->columns_57_to_65, "Mean     "))
      coord_epoch = 0.;
   else if( !memcmp( obs->columns_57_to_65, "Epoch", 5))
      coord_epoch = atof( obs->columns_57_to_65 + 5);
   if( coord_epoch != 2000.)
      {
      double year = (obs->jd - J2000) / 365.25 + 2000.;

      if( coord_epoch == -1.)       /* true coords of date */
         {
         double matrix[9];

         setup_precession_with_nutation( matrix, year);
         precess_ra_dec( matrix, &obs->ra, &obs->ra, 1);
         adjust_for_aberration( obs);
         }
      else
         {
         if( coord_epoch)               /* specific epoch given, not just */
            year = coord_epoch;         /* mean coords of date */
         obs->ra *= -1.;
         precess_pt( (DPT DLLPTR *)&obs->ra,
                     (DPT DLLPTR *)&obs->ra, year, 2000.);
         obs->ra *= -1.;
         }
      obs->note2 = 'A';       /* mark as coordinates precessed */
      }
   if( obs->note2 == 'a')     /* input coords are alt/az,  not RA/dec */
      {
      DPT latlon, alt_az, ra_dec;

      if( get_observer_data_latlon( obs->mpc_code, NULL,
              &latlon.x, &latlon.y, NULL) == 3)
         {
         alt_az.x = obs->ra;     /* because input coords were really alt/az */
         alt_az.y = obs->dec;
         full_alt_az_to_ra_dec( &ra_dec, &alt_az, utc, &latlon);
         obs->ra = -ra_dec.x;
         obs->dec = ra_dec.y;
         }
      }
   if( isdigit( buff[66]) && obs->note2 != 'R' &&
               (buff[67] == '.' || buff[67] == ' '))
      obs->obs_mag = atof( buff + 65);
   else
      obs->obs_mag = BLANK_MAG;
   FMEMCPY( obs->reference, buff + 72, 5);
   obs->reference[5] = '\0';
   if( memcmp( buff + 72, ".rwo ", 5) &&
          find_fcct_biases( obs->ra, obs->dec, obs->mag_band2, obs->jd,
                                &obs->ra_bias, &obs->dec_bias) == -2)
      {        /* i.e.,  we tried to get FCCT14 debiasing and failed */
      if( !fcct_error_message_shown && apply_debiasing)
         {
         strcpy( tbuff, "Didn't find observational bias data file\n"
                        "To include the FCCT14 bias data,  Find_Orb needs\n"
                        "the file 'bias.dat'.  There is more information\n"
                        "about this at http://www.projectpluto.com/bias.htm\n"
                        "Biases will not be applied until the file is available.\n");
         generic_message_box( tbuff, "o");
         fcct_error_message_shown = true;
         }
      }
   set_up_observation( obs);
   return( 0);
}

int separate_periodic_comet_apparitions = 0;

           /* In theory,  an object should either have a permanent number, */
           /* or a temporary designation.  In practice,  some have both. */
           /* The 'separate_periodic_comet_apparitions' lets you tell    */
           /* Find_Orb which should be used when both are given;  if it's   */
           /* set to TRUE,  it'll zap the permanent designation data     */
           /* (i.e.,  set it to spaces.)  Otherwise,  it'll zap the      */
           /* temporary designation.                                     */
           /*    For natural satellites,  the situation is a little      */
           /* different.  There,  we zap the _permanent_ version.        */

static void reduce_designation( char *desig, const char *idesig)
{
   memcpy( desig, idesig, 12);
   desig[12] = '\0';
   if( desig[5] != ' ' && strlen( desig) > 9)
      if( isdigit( desig[1]) && isdigit( desig[2]) && isdigit( desig[3]))
         {
         if( strchr( "PCDXA", desig[4]))        /* it's a numbered comet... */
            if( isdigit( desig[0]) && isdigit( desig[1]))
               {
               if( separate_periodic_comet_apparitions)
                  memset( desig, ' ', 4);
               else
                  memset( desig + 5, ' ', 7);
               }

           /* Same problem exists for numbered asteroids,  sometimes: */
         if( isdigit( desig[0]) || isalpha( desig[0]))
            if( isdigit( desig[4]))
               {
               if( separate_periodic_comet_apparitions)
                  memset( desig, ' ', 5);  /* NB: _5_,  not _4_,  for these */
               else
                  memset( desig + 5, ' ', 7);
               }
         if( desig[4] == 'S')
            memset( desig, ' ', 4);
         }
}

/* If combine_all_observations is non-zero,  then all observations from */
/* the file are treated as if they're from a single object.  This can be */
/* useful if you're trying to link two arcs,  or if the designations from */
/* different sources aren't the same.  */

int combine_all_observations = 0;

static int compare_desigs( const char *desig1, const char *desig2)
{
   int rval;

   assert( *desig1);
   assert( *desig2);
   if( combine_all_observations)
      rval = 0;
   else
      {
      char reduced1[13], reduced2[13];

      reduce_designation( reduced1, desig1);
      reduce_designation( reduced2, desig2);
      rval = strcmp( reduced1, reduced2);
      }
   return( rval);
}

/* In load_observations(),  we need to know if a given "new" observation
is of an already-found object.  To do that,  we keep track of designations
in a hash table of OBJECT_INFO structures.  The following function either
finds the index of the desired designation within the hash table,  or it
finds a blank entry where that designation ought to go.   */

static unsigned find_in_hash_table( const OBJECT_INFO *objs, const char *desig,
                                    const unsigned table_size)
{
   if( combine_all_observations)
      return( 0);
   else
      {
      char reduced[13];
      unsigned loc = 42, i;
      uint32_t *tptr = (uint32_t *)reduced;  /* horrible type pun */

      reduce_designation( reduced, desig);
      for( i = 0; i < 3; i++)
         {
         loc ^= (unsigned)*tptr++;
         loc ^= loc << 13;
         loc ^= loc >> 8;
         }
      loc %= table_size;
      i = 1;
      assert( *desig);
      while( objs[loc].packed_desig[0]
                  && compare_desigs( objs[loc].packed_desig, desig))
         {
         loc = (loc + (i > table_size ? 1 : i)) % table_size;
         i++;
         }
      return( loc);
      }
}

int qsort_strcmp( const void *a, const void *b, void *ignored_context)
{
   return( strcmp( (const char *)a, (const char *)b));
}

static int xref_designation( char *desig)
{
   static char *xlate_table = NULL;
   static char prev_desig_in[12], prev_desig_out[12];
   static int n_lines = 0;
   char reduced_desig[13];
   int i, gap;

   if( !desig)                /* free up memory */
      {
      if( xlate_table)
         free( xlate_table);
      xlate_table = NULL;
      n_lines = 0;
      return( 0);
      }

   if( !xlate_table)
      {
      FILE *ifile = fopen_ext( "xdesig.txt", "fcrb");
      char buff[100];
      int j;

      assert( ifile);
      while( fgets( buff, sizeof( buff), ifile))
         if( *buff != ';' && *buff >= ' ')
            n_lines++;
      xlate_table = (char *)malloc( n_lines * 26);

      fseek( ifile, 0L, SEEK_SET);
      i = 0;
      while( fgets( buff, sizeof( buff), ifile))
         if( *buff != ';' && *buff >= ' ')
            {
            for( j = 0; buff[j] && buff[j] != 10; j++)
               ;
            while( j < 25)
               buff[j++] = ' ';
            buff[25] = '\0';
            strcpy( xlate_table + i * 26, buff);
            i++;
            }
      fclose( ifile);
      shellsort_r( xlate_table, n_lines, 26, qsort_strcmp, NULL);
      }

            /* Frequently,  this function is given the same designation */
            /* multiple times in a row.  You can get at least some speedup */
            /* by caching the previous input,  checking to see if it */
            /* matches the new input,  and if it does,  return the previous */
            /* output designation. */
   if( !memcmp( desig, prev_desig_in, 12))
      {
      memcpy( desig, prev_desig_out, 12);
      return( 0);
      }

   memcpy( prev_desig_in, desig, 12);
   reduce_designation( reduced_desig, desig);

   for( i = 0, gap = 0x8000; gap; gap >>= 1)
      if( i + gap < n_lines)
         {
         char *xlate_ptr = xlate_table + (i + gap) * 26;
         const int compare = memcmp( xlate_ptr, reduced_desig, 12);

         if( compare <= 0)
            i += gap;
         if( !compare)    /* it's a match */
            memcpy( desig, xlate_ptr + 13, 12);
         }
   memcpy( prev_desig_out, desig, 12);
   return( 0);
}

static int is_second_line( const char *buff)
{
   if( buff[14] == 's' || buff[14] == 'r' || buff[14] == 'v')
      return( buff[14]);
   else
      return( 0);
}

/* The following allows one to have Find_Orb accept only data from
given station(s),  or reject all data from given station(s).  For
example,  if you add

REJECTED_STATIONS=Z42,258,Q99

to 'environ.dat',  all observations from those three stations would
be ignored.  But make that

ACCEPTED_STATIONS=Z42,258,Q99

   and only data from those three would be used by Find_Orb.  I've
used,  for example,

ACCEPTED_STATIONS=253,251

   to make "radar-only" orbits.   */

static bool observatory_is_acceptable( const char *mpc_code)
{
   static const char *rejects = NULL;
   static const char *filter = NULL;

   if( !rejects)
      rejects = get_environment_ptr( "REJECTED_STATIONS");
   if( !filter)
      filter = get_environment_ptr( "ACCEPTED_STATIONS");
   if( *rejects)
      return( strstr( rejects, mpc_code) ? false : true);
   if( *filter)
      return( strstr( filter, mpc_code) ? true : false);
   return( true);
}

/*  The following does some simple error checking to 'buff' to see if it's
    really an observation,  or just a random line in the file.  If it's
    for-real data,  the observation JD is returned.  Otherwise,  you get 0.  */

static double observation_jd( const char *buff)
{
   const double utc = extract_date_from_mpc_report( buff, NULL);

   if( utc && is_valid_mpc_code( buff + 77))
      {
#ifdef STRICT_OBSERVATION_CHECKING
      static const char digits[] =
            { 15, 16, 17, 18, 21, 23, 24, 26, 0 };
      int i;

      for( i = 0; digits[i]; i++)
         if( buff[digits[i]] < '0' || buff[digits[i]] > '9')
            return( 0.);
#endif
      if( !observatory_is_acceptable( buff + 77))
         return( 0.);
      }
   return( utc);
}

/* Depending on how it was saved,  an NEOCP ephemeris can be plain text;
plain text modified in weird ways by the browser (at least Firefox) so
that bold items are "decorated" with asterisks;  or it might be real HTML
with real HTML <p><b> and </b> tags around the object designations. Example
of each of the three:

<p><b>VE82F84</b>
*VE82F84*
VE82F84
 */

#define NEOCP_FILE_TYPE_UNKNOWN           0
#define NEOCP_FILE_TYPE_HTML              1
#define NEOCP_FILE_TYPE_ASTERISKED        2

/* The code can handle ephemerides with either the default step size of
one hour,  or others with smaller step sizes.  The two look mostly
alike,  differing only in whether the minutes column contains spaces
or two digits:

Date       UT   *  R.A. (J2000) Decl.  Elong.  V        Motion     Object     Sun         Moon        Uncertainty
            h m                                      "/min   P.A.  Azi. Alt.  Alt.  Phase Dist. Alt.
2010 12 31 0540   06 39 55.0 +26 35 07 176.5  16.8   24.75  196.4  293  +63   -61    0.18  130  -60
2010 12 31 05     06 39 55.0 +26 35 07 176.5  16.8   24.75  196.4  293  +63   -61    0.18  130  -60
             ^^ Note difference here

   Hence the lines starting at the comment "ephem step size is in hours":
we're checking the first 26 bytes to see which have ASCII digits in them.  */

static int neocp_file_type;

static bool get_neocp_data( char *buff, char *desig, char *mpc_code)
{
   size_t len;
   bool rval = false;

   for( len = 0; buff[len] >= ' '; len++)
      ;
   if( debug_level > 9)
      debug_printf( "get_neocp_data len: %zd\n", len);
   while( len && buff[len - 1] == ' ')    /* eliminate trailing spaces */
      len--;
   if( len < 10 && len > 2 && *buff == '*' && buff[len - 1] == '*')
      {
      memcpy( desig, buff + 1, len - 2);
      desig[len - 2] = '\0';
      neocp_file_type = NEOCP_FILE_TYPE_ASTERISKED;
      }
   else if( !memcmp( buff, " <p><b>", 7) && !memcmp( buff + 14, "</b>", 4))
      {
      memcpy( desig, buff + 7, 7);
      desig[7] = '\0';
      neocp_file_type = NEOCP_FILE_TYPE_HTML;
      }
   else if( len && len < 10 && neocp_file_type == NEOCP_FILE_TYPE_UNKNOWN)
      {
      memcpy( desig, buff, len);
      desig[len] = '\0';
      }
   else if( !memcmp( buff, "Ephemerides are for ", 20))
      {
      static int already_warned = 0;

      if( !memcmp( buff + 20, "observatory code ", 17))
         memcpy( mpc_code, buff + 37, 3);
      else if( !already_warned)
         {
         char msg[500];

         already_warned = 1;
         strcpy( msg, "Find_Orb can use the NEOCP ephemeris you've given it.\n");
         strcat( msg, "However,  the ephemerides are geocentric.  You'll\n");
         strcat( msg, "get much better results if you select an observatory\n");
         strcat( msg, "code,  preferably one near the equator.  (781) Quito\n");
         strcat( msg, "is best for this purpose.\n");
         memcpy( mpc_code, "500", 3);
         generic_message_box( msg, "o");
         }
      }
   else if( !memcmp( buff, " observatory code ", 18))
      memcpy( mpc_code, buff + 18, 3);
   else if( len > 64 && buff[26] == '.' && *desig && *mpc_code)
      {
      int i;
      const long mask_if_hours   = 0x36c1b6f;
      const long mask_if_minutes = 0x36c7b6f;
      long mask = 0;

      for( i = 0; i < 26; i++)
         if( isdigit( buff[i]))
            mask |= (1L << i);
         else if( buff[i] != ' ')
            i = 99;                    /* mark as bad */
      if( i == 26 && (mask == mask_if_hours || mask == mask_if_minutes))
         {
         char obuff[82];
         int minutes;

         memset( obuff, ' ', 80);
         obuff[80] = '\0';
         strcpy( obuff + 5, desig);
         memcpy( obuff + 15, buff, 10);     /* year, month, day */
//       sprintf( obuff + 25, ".%05d", atoi( buff + 11) * 100000 / 24);
         if( mask == mask_if_hours)     /* ephem step size is in hours */
            minutes = atoi( buff + 11) * 60;
         else                       /* step size was in minutes */
            minutes = (atoi( buff + 11) / 100) * 60 +
                                          atoi( buff + 13);
         sprintf( obuff + 25, ".%06d", minutes * 1000000 / (int)minutes_per_day);
         memcpy( obuff + 32, buff + 18, 10);      /* RA */
         memcpy( obuff + 44, buff + 29, 9);       /* dec */
         memcpy( obuff + 65, buff + 46, 4);       /* mag */
         strcpy( obuff + 72, "neocp");
         memcpy( obuff + 77, mpc_code, 3);
         obuff[14] = 'C';                          /* assume a CCD obs */
         obuff[70] = 'V';                          /* mags are V */
         for( i = 0; i < 80; i++)
            if( !obuff[i])
               obuff[i] = ' ';
         strcpy( buff, obuff);
         rval = true;
         }
      }
   return( rval);
}

   /* In the RWO (NEODyS) format,  provisional designations are given */
   /* in the form ' YYYYLL(num)',  where YYYY is the year and LL are  */
   /* two letters A-Z,  and (num) is an optional number.              */
   /* Declared 'inline' mostly to emphasise that it's called in only  */
   /* one place,  from rwo_to_mpc.                                    */

static void inline reformat_rwo_designation_to_mpc( const char *buff, char *obuff)
{
   if( !isalpha( buff[6]))           /* numbered object */
      {
      const long ast_number = atoi( buff);
      const long leading_digit = ast_number / 10000L;

      *obuff = int_to_mutant_hex_char( leading_digit);
      if( !*obuff)      /* mutant hex fails past asteroid 619999 */
         *obuff = '?';  /* ...just put _something_ there         */
      sprintf( obuff + 1, "%04ld", ast_number % 10000L);
      obuff[5] = ' ';
      }
   else                       /* provisional designation */
      {
      obuff[5] = (char)( 'K' + (buff[1] - '2') * 10 + buff[2] - '0');
                                    /* obuff[5] is century marker */
      obuff[6] = buff[3];    /* decade */
      obuff[7] = buff[4];    /* year */
      obuff[8] = buff[5];    /* first letter */
      obuff[11] = buff[6];    /* 2nd   letter */
                  /* A provisional designation is a four-digit year, */
                  /* plus two letters,  plus either no number,  or a */
                  /* one, two,  or three-digit number. */
      if( buff[7] == ' ')              /* No number: say,  '2004RW' */
         obuff[9] = obuff[10] = '0';
      else if( buff[8] == ' ')         /* One-digit #: say,  '2004RW1' */
         {
         obuff[9] = '0';
         obuff[10] = buff[7];
         }
      else if( buff[9] == ' ')         /* 2-digit #: say,  '2004RW31' */
         {
         obuff[9] = buff[7];    /* 1st (tens) digit */
         obuff[10] = buff[8];    /* 2nd (units) digit */
         }
      else                             /* 3-digit #: say,  '2004RW314' */
         {
         const int tens = buff[7] * 10 + buff[8] - '0' * 11;

         obuff[9] = int_to_mutant_hex_char( tens);
         obuff[10] = buff[9];    /* 3rd (units) digit */
         }
      }
}

/* The AstDyS/NEODyS .rwo astrometry format provides almost all of
the data of the MPC's 80-column punched-card format.  However,  it
doesn't give the frequency for radar observations.  Fortunately,  since
the radar folks appear to have been fairly consistent in the frequencies
used,  we can (as of early 2014) determine the frequency based on year
and MPC code.  I figured this out based on the accumulated radar data at

http://ssd.jpl.nasa.gov/?grp=ast&fmt=html&radar=             */

static inline double get_radar_frequency( const int mpc_code, const int year)
{
   double freq_in_mhz;

   switch( mpc_code)
      {
      case 251:      /* Arecibo uses 2380 MHz almost all the time */
         if( year == 1975)    /* except its first two obs of Eros */
            freq_in_mhz = 430.;     /* on 1975 Jan 22 */
         else
            freq_in_mhz = 2380.;
         break;
      case 252:      /* Goldstone DSS 13, Fort Irwin */
      case 253:      /* Goldstone DSS 14, Fort Irwin */
      case 255:      /* Evpatoria */
      case 257:      /* Goldstone DSS 25, Fort Irwin */
         if( year > 1998)     /* obs from 1999 Mar 22 to present */
            freq_in_mhz = 8560;
         else if( year > 1990)    /* from 1991 jun 15 to 1998 aug 16 */
            freq_in_mhz = 8510;
         else if( year > 1974)    /* 1975 jan 23 - 1990 aug 3 */
            freq_in_mhz = 8495;
         else
            freq_in_mhz = 2388.;    /* Icarus observations in 1968 */
         break;
      case 254:          /* Haystack,  Westford */
         freq_in_mhz = 7840.;    /* two Icarus observations in 1968 */
         break;
      case 272:            /* so far,  they've made one observation, */
         freq_in_mhz = 929.;     /* of (387943) Duende = 2012 DA14  */
         break;
      default:                   /* unknown station;  frequency     */
         freq_in_mhz = 5000.;    /* chosen at random using fair die */
         break;
      }
   return( freq_in_mhz);
}

/* Code to convert an observation in the AstDyS and NeoDyS '.rwo' format
to the MPC's 80-byte format.  Return value indicates that the input
buffer _was_ in .rwo format (rval = 1) or was not (rval = 0),  in which
case the input buffer is unaltered.

12 Aug 2005:  revised to also accept the revised AstDyS format.

2010 Nov 2: revised to accept 'second line' data for satellite observations.
   In such observations,  the second line gives the xyz offset from the
   geocenter for the satellite.  Note that some pretty odd reformatting is
   necessary,  in particular to put the signs for the coordinates in the
   "proper" (according to MPC) columns.  Also fixed so that the observation
   note is copied over, and so that accuracies in RA and dec are respected
   (so that the trailing zero is dropped).

20 Jan 2011:  Revised to accept cases where the .rwo line includes the new
   catalog reference code.  When that happens,  it inserts four columns
   in the .rwo line.


   Example satellite data,  in RWO and MPC formats:

 53037     S s   2010 01 31.27280 1   5432.8041   3919.2531   1679.3707 C51
 2010SO16  S s   2010 09 17.94694 1    547.9162  -2632.9667  -6366.8349 C51
 2004GL32  S s   2010 06 16.32490 1   6826.5602   -950.7342    464.5735 C51
     K04G32L  s2010 06 16.32490 1 + 6826.5602 -  950.7342 +  464.5735   ~0J7nC51
     K10S16O  s2010 09 17.28554 1 +  520.1210 - 2628.1408 - 6370.9902   ~0NK7C51
 2004GL32  S s   2010 06 15.46487 1   6810.5652  -1047.8458    490.0909 C51
     K04G32L  s2010 06 15.46487 1 + 6810.5652 - 1047.8458 +  490.0909   ~0J7nC51

Further notes about .rwo and MPC radar data:  the former stores range and
range-rate data separately,  whereas MPC's format can stuff both into one
record.

NOTE that .rwo "distance" is half the round-trip distance,  and MPC
"distance" is round-trip travel time in microseconds;  hence some
factors of two to cause grief.  An example of a range and range-rate
observation in both formats,  and the conversion factors.  Note that
for the MPC's decimal day format,  we could instead use Find_Orb's
extension for HH:MM:SS format;  instead of storing "R1992 12 18.298611",
it would make sense to store "RJ921218:071000".

! Object   Obser ====== Date =======  ============ Radar range/range rate (km or km/d) ============= Station    Residual
! Design   K T N YYYY MM DD hh:mm:ss        Measure     Accuracy    rms    F      Bias       Resid   TRX RCX     Chi   S
 4179      R c   1992 12 18 07:10:00  10313006.57062   0.28180   0.28180 F     0.00000     0.06586 253 253      0.24 1
 4179      V c   1992 12 18 07:10:00    871504.67532   0.54787   0.54787 F     0.00000     0.00535 253 253      0.01 1
    RTDist 20625987 km = 68.80088959s; Dopp -20.17441496 km/s = -572677.09 Hz
 20625987 km = 68.80088959s * 299792.458 km/s;  halve to get measure in .rwo
 .2818 km = 0.00000188s * 299792.458 km/s / 2
 871504.67 km/d = 86400 s/d * 299792.458 km/s * 572657.34 Hz / 8510 MHz / 2
 .54787 km/d =    86400 s/d * 299792.458 km/s *       .36 Hz / 8510 MHz / 2

04179         R1992 12 18.298611   6880097411  -    57265734   8510 253 JPLRS253
04179         r1992 12 18.298611C         1880           036        253 JPLRS253

   One-way distance is 10313006.57062 km.  Double that and divide by
c=299792.458 km/s,  and the round-trip travel time (RTT) is
68.80097411 seconds.

   The .rwo "velocity" is 871504.67532 km/day.  Again double it,  divide
by 24*60*60 = 86400 to get 20.173719336 km/s,  then multiply by 8510 MHz
and divide by c=299792.458 km/s to get a shift of 572657.34 Hz.  The
sigmas are similarly converted.

*/

#define MINIMUM_RWO_LENGTH 117

static int rwo_to_mpc( char *buff, double *ra_bias, double *dec_bias,
                  double *posn_sigma_ra, double *posn_sigma_dec, double *mag_sigma)
{
   int rval = 0, i, j;
   const size_t line_len = strlen( buff);
   char obuff[82], second_radar_line[82];

   *second_radar_line = '\0';
   if( debug_level > 2)
      debug_printf( "rwo_to_mpc: Input: '%s'\n", buff);
   if( line_len == 75 && buff[27] == '.' && buff[54] == '.'
               && buff[13] == 's' && buff[12] == ' ')
      {                       /* satellite second line */
      if( debug_level > 2)
         debug_printf( "Input: '%s'\n", buff);
      memset( obuff, ' ', 80);
      obuff[14] = buff[13];   /* obs type */
      for( i = 15; i < 33; i++)     /* Year,  month,  day */
         obuff[i] = buff[i + 2];
      obuff[32] = buff[34];      /* units specifier */
      memcpy( obuff + 34, buff + 36, 36);  /* xyz */
      memcpy( obuff + 77, buff + 72, 3);   /* MPC code */
      for( i = 34; i < 59; i += 12)
         {
         for( j = i + 1; j < i + 12; j++)
            if( obuff[j] == '-' || obuff[j] == '+')
               {
               obuff[i] = obuff[j];
               obuff[j] = ' ';
               }
         if( obuff[i] == ' ')   /* if no sign given, make coordinate */
            obuff[i] = '+';     /* explicitly positive */
         }
      rval = 1;
      if( debug_level > 2)
         debug_printf( "Output: '%s'\n", obuff);
      }
   else if( line_len < MINIMUM_RWO_LENGTH)
      return( 0);

   if( line_len > 192 && *buff == ' ' && buff[16] == ' '
                && buff[21] == ' ' && buff[27] == '.' && buff[58] == '.'
                && (buff[131] == '.' || buff[122] == 'E'))
      {
      const int is_version_two = (buff[176] == ' ');
            /* "Version 2" .rwos have a code indicating the star catalog   */
            /* used,  with columns shifted by four bytes after that point. */

      memset( obuff, ' ', 80);
      obuff[14] = buff[13];   /* obs type */
      for( i = 15; i < 33; i++)     /* Year,  month,  day */
         obuff[i] = buff[i + 2];
      for( i = 32; i < 44; i++)     /* RA */
         obuff[i] = buff[i + 18];
      for( i = 44; i < 56; i++)     /* dec */
         obuff[i] = buff[i + 59];
      for( i = 65; i < 71; i++)     /* mag,  color */
         obuff[i] = buff[i + 91];
      for( i = 77; i < 80; i++)     /* MPC station code */
         obuff[i] = buff[i + (is_version_two ? 103 : 99)];
      if( buff[is_version_two ? 194 : 190] == '0')
         obuff[64] = 'x';           /* excluded observation */
      obuff[13] = buff[15];         /* note */
      if( is_version_two)
         obuff[71] = buff[178];     /* catalog code */
      if( !memcmp( buff + 64, "1.500E", 6)) /* RA is to full seconds,  or  */
         {                                  /* .1,  or .01 second          */
         int digits_to_drop = atoi( buff + 70) + 2;

         assert( digits_to_drop >= 0 && digits_to_drop <= 3);
         if( digits_to_drop == 3)      /* drop the decimal point too */
            digits_to_drop = 4;
         memset( obuff + 44 - digits_to_drop, ' ', digits_to_drop);
         }
      if( !memcmp( buff + 117, "1.000E", 6)) /* dec is to .1 arcsec;*/
         {
         const int digits_to_drop = atoi( buff + 123) + 2;

         assert( digits_to_drop >= 0 && digits_to_drop < 3);
         memset( obuff + 56 - digits_to_drop, ' ', digits_to_drop);
         }
      if( ra_bias)
         *ra_bias = atof( buff + 86);
      if( dec_bias)
         *dec_bias = atof( buff + 139);
      if( posn_sigma_ra)
         {
         if( buff[78] == '.')
            *posn_sigma_ra = atof( buff + 74);
         else
            *posn_sigma_ra = 1.2;
         }
      if( posn_sigma_dec)
         {
         if( buff[131] == '.')
            *posn_sigma_dec = atof( buff + 127);
         else
            *posn_sigma_dec = 1.2;
         }
      if( mag_sigma && buff[165] == '.')
         *mag_sigma = atof( buff + 164);
      if( line_len >= 203)                   /* reference supplied */
         {
         memcpy( obuff + 72, buff + 198, 5);
         memcpy( obuff + 58, " rwo", 4);
         }
      else              /* no reference... just mark it as a RWO */
         memcpy( obuff + 72, ".rwo ", 5);
      rval = 1;
      }
   else if( buff[30] == ':' && buff[33] == ':' &&
                                    strlen( buff) == 118)
      {                          /* Radar data: not really handled yet */
      double val1, val2;

      if( sscanf( buff + 36, "%lf %lf", &val1, &val2) == 2)
         {
         const int year = atoi( buff + 17);
         double freq_in_mhz;

         buff[36] = '\0';
         memset( obuff, ' ', 80);
         obuff[14] = 'R';
         obuff[80] = '\0';
         obuff[15] = (char)('A' + (year / 100) - 10);   /* century */

         obuff[16] = buff[19];   /* decade */
         obuff[17] = buff[20];   /* year */
         obuff[18] = buff[22];   /* tens of months */
         obuff[19] = buff[23];   /* months */
         obuff[20] = buff[25];   /* tens of days */
         obuff[21] = buff[26];   /* days */
         obuff[22] = ':';
         obuff[23] = buff[28];   /* tens of hours */
         obuff[24] = buff[29];   /* hours */
         obuff[25] = buff[31];   /* tens of minutes */
         obuff[26] = buff[32];   /* minutes */
         obuff[27] = buff[34];   /* tens of seconds */
         obuff[28] = buff[35];   /* seconds */
         memcpy( obuff + 68, buff + 99, 3);    /* transmitting MPC code */
         memcpy( obuff + 77, buff + 103, 3);    /* receiving MPC code */
         strcpy( second_radar_line, obuff);
         freq_in_mhz = get_radar_frequency( atoi( buff + 99), year);
         sprintf( obuff + 62, "%5.0f", freq_in_mhz);
         second_radar_line[14] = 'r';
         if( buff[11] == 'R')
            {
            val1 *= 2. / SPEED_OF_LIGHT;
            val2 *= 2. / SPEED_OF_LIGHT;
            if( debug_level > 2)
               debug_printf( "%s: round trip %f +/- %f microseconds\n",
                     buff, val1 * 1e+6, val2 * 1e+6);
            sprintf( obuff + 32, "%13.0f", val1 * 1e+8);
            sprintf( second_radar_line + 34, "%12.0f", val2 * 1e+9);
            rval = 1;
            }
         else if( buff[11] == 'V')
            {
            val1 *= 2. * freq_in_mhz / (SPEED_OF_LIGHT * seconds_per_day);
            val2 *= 2. * freq_in_mhz / (SPEED_OF_LIGHT * seconds_per_day);
            val1 = -val1;     /* .rwo & MPC have different sign convention */
            if( debug_level > 2)
               debug_printf( "%s: Doppler %f +/- %f Hz\n",
                     buff, val1, val2);
            sprintf( obuff + 48, "%13.0f", fabs( val1) * 1e+9);
            obuff[47] = (val1 > 0. ? '+' : '-');
            sprintf( second_radar_line + 48, "%13.0f", val2 * 1e+9);
            rval = 1;
            }
         if( !memcmp( obuff + 77, "251", 3) || !memcmp( obuff + 77, "253", 3))
            memcpy( obuff + 72, "JPLRS", 5);
         buff[36] = ' ';
         for( i = 0; i < 80; i++)
            if( !obuff[i])
               obuff[i] = ' ';
         for( i = 0; i < 80; i++)
            if( !second_radar_line[i])
               second_radar_line[i] = ' ';
         strcpy( second_radar_line + 68, obuff + 68);
         }
      }
   if( rval)
      {
      obuff[80] = '\0';
      reformat_rwo_designation_to_mpc( buff, obuff);
      strcpy( buff, obuff);
      if( *second_radar_line)
         {
         strcpy( buff + 81, second_radar_line);
         memcpy( buff + 81, buff, 13);
         }
      }
   return( rval);
}

#ifndef memicmp
#ifndef __WATCOMC__
static int memicmp( const char *s1, const char *s2, int n)
{
   int c1, c2;

   while( n--)
      {
      if( (c1 = tolower( *s1++)) != (c2 = tolower( *s2++)))
         return( c1 - c2);
      }
   return( 0);
}
#endif
#endif

   /* Satellite,  radar,  and roving observer observations are stored     */
   /* as two-line pairs in the MPC 80-column punched-card format.  Sadly, */
   /* you can't count on them being in the "proper" order.  (True of the  */
   /* satellite obs,  anyway;  I've not verified any errors with radar or */
   /* roving observations.  If there are any,  though,  this code should  */
   /* repair the errors.)                                                 */
   /*    To fix this,  when we encounter a satellite/radar/roving obs,    */
   /* we check to see if its matching line has already been read.  If we  */
   /* find it,  we return it in oline; else,  we store 'iline' in the     */
   /* stored_line list,  waiting for a mate.  The total number of stored  */
   /* lines is returned.                                                  */
   /*   At the end,  the total number of stored lines _should_ be zero;   */
   /* i.e.,  for each line,  a matching line should have been found.  If  */
   /* that's not the case,  Find_Orb shows a warning message to let you   */
   /* know those observations weren't used.                               */

typedef char mpc_line[81];

static int look_for_matching_line( char *iline, char *oline)
{
   static int n_stored = 0;
   static mpc_line *stored_lines = NULL;
   int i;

   if( !iline)                /* just checking for leftover lines */
      {
      const int rval = n_stored;

      if( n_stored)
         debug_printf( "%d unmatched satellite/roving observer lines:\n",
                        n_stored);
      for( i = 0; i < n_stored; i++)
         debug_printf( "%s\n", stored_lines[i]);
      n_stored = 0;
      if( stored_lines)
         {
         free( stored_lines);
         stored_lines = NULL;
         }
      return( rval);
      }
   *oline = '\0';    /* assume no match found */
   for( i = 0; i < n_stored; i++)
      {
              /* Lines should compare to byte 31,  except possibly for */
              /* a discovery asterisk in column 13 _or_ (added 2014 Oct */
              /* 17) a note in column 14:  */
      if( !memicmp( stored_lines[i], iline, 12) &&
                      !memicmp( stored_lines[i] + 14, iline + 14, 17))
         {        /* we have a match */
         memcpy( oline, stored_lines[i], sizeof( mpc_line));
         n_stored--;
               /* Move last line into place formerly used by the */
               /* now-matched line: */
         memcpy( stored_lines[i], stored_lines[n_stored], sizeof( mpc_line));
         }
      }
   if( !*oline)    /* we didn't find a match: */
      n_stored++;
            /* We've now either found a match,  and therefore removed a */
            /* line from 'stored_lines';  or didn't,  and are about to  */
            /* add a line.  Either way,  we need to realloc:            */
   stored_lines = (mpc_line *)realloc( stored_lines, n_stored * sizeof( mpc_line));
   if( !*oline)                /* store unmatched line: */
      memcpy( stored_lines[n_stored - 1], iline, sizeof( mpc_line));
   else if( iline[14] >= 'a')  /* lines are backward;       */
      {                        /* swap 'em to proper order  */
      char swap_buff[80];

      memcpy( swap_buff, iline, 80);
      memcpy( iline, oline, 80);
      memcpy( oline, swap_buff, 80);
      }
   return( *oline ? 1 : 0);
}

/* Occasionally,  observations are duplicates,  except that one has */
/* a blank reference or magnitude and the other doesn't.  In this   */
/* function,  we'll copy the reference that _is_ given over the blank, */
/* and the mag that's given over the zero mag.  Then we check to see */
/* if that's caused the observations to become the same.  If it has, */
/* we consider the difference to have been "corrected".              */

static void correct_differences( OBSERVE *obs1, const OBSERVE *obs2)
{
   if( obs1->reference[0] == ' ' && obs2->reference[0] != ' ')
      strcpy( obs1->reference, obs2->reference);
   if( obs1->obs_mag == BLANK_MAG && obs2->obs_mag != BLANK_MAG)
      {
      obs1->obs_mag = obs2->obs_mag;
      obs1->mag_precision = obs2->mag_precision;
      obs1->mag_band = obs2->mag_band;
      }
   if( obs1->mag_band2 == ' ' && obs2->mag_band2 != ' ')
      obs1->mag_band2 = obs2->mag_band2;
}

static int compare_observations( const void *a, const void *b, void *unused_context)
{
   const OBSERVE *obs1 = (const OBSERVE *)a;
   const OBSERVE *obs2 = (const OBSERVE *)b;
   int rval;

   if( obs1->jd < obs2->jd)
      rval = -1;
   else if( obs1->jd > obs2->jd)
      rval = 1;
   else
      rval = FMEMCMP( obs1->mpc_code, obs2->mpc_code, 3);
   if( !rval)
      rval = obs1->note1 - obs2->note1;
   if( !rval && obs1->ra != obs2->ra)
      rval = (obs1->ra > obs2->ra ? 1 : -1);
   if( !rval && obs1->dec != obs2->dec)
      rval = (obs1->dec > obs2->dec ? 1 : -1);
   return( rval);
}

/* The re-entrant version of qsort(),  qsort_r(),  is implemented in
subtly but infuriatingly different ways on different platforms (including
as 'qsort_s' on MS Windows).  You can "sort of" work around these,  as

https://github.com/noporpoise/sort_r

   shows.  But I found it easier to write the following.  At least right
now,  I'm sorting arrays small enough that Shellsort is more than adequate. */

void shellsort_r( void *base, const size_t n_elements, const size_t esize,
         int (*compare)(const void *, const void *, void *), void *context)
{
// printf( "Sort begins: %f, %u elements, size %u\n",
//                  (double)clock( ) / (double)CLOCKS_PER_SEC,
//                  (unsigned)n_elements, (unsigned)esize);
#if (defined _GNU_SOURCE || defined __GNU__ || defined __linux)
   qsort_r( base, n_elements, esize, compare, context);
#else
   size_t gap = 1;
   char *data = (char *)base;
   char *pivot = (char *)alloca( esize);

   while( gap < n_elements)
      gap = gap * 3 + 1;
   while( gap /= 3)
      {
      size_t j;
      const size_t spacing = esize * gap;

      for( j = gap; j < n_elements; j++)
         {
         char *tptr = data + j * esize;
         char *tptr2 = tptr - spacing;

         if( (compare)( tptr2, tptr, context) > 0)
            {
            memcpy( pivot, tptr, esize);
            memcpy( tptr, tptr2, esize);
            tptr = tptr2;
            tptr2 -= spacing;
            while( tptr2 >= base && (compare)( tptr2, pivot, context) > 0)
               {
               memcpy( tptr, tptr2, esize);
               tptr = tptr2;
               tptr2 -= spacing;
               }
            memcpy( tptr, pivot, esize);
            }
         }
      }
#endif
// printf( "Sort ends: %f\n", (double)clock( ) / (double)CLOCKS_PER_SEC);
}

/* Does what the function name suggests.  Return value is the number
of observations left after duplicates have been removed.   */

int sort_obs_by_date_and_remove_duplicates( OBSERVE *obs, const int n_obs)
{
   int i, j;

   if( !n_obs)
      return( 0);
   shellsort_r( obs, n_obs, sizeof( OBSERVE), compare_observations, NULL);
   if( debug_level)
      debug_printf( "%d obs sorted by date\n", n_obs);
   for( i = j = 1; i < n_obs; i++)
      if( memcmp( obs + i, obs + i - 1, sizeof( OBSERVE)))
         {
         OBSERVE temp1 = obs[i], temp2 = obs[i - 1];

         correct_differences( &temp1, &temp2);
         correct_differences( &temp2, &temp1);
         if( !memcmp( &temp1, &temp2, sizeof( OBSERVE)))
            {
            char buff[80];

            full_ctime( buff, temp1.jd, CALENDAR_JULIAN_GREGORIAN
                        | FULL_CTIME_FORMAT_DAY | FULL_CTIME_6_PLACES);
            debug_printf( "Observation %d from %s on %s is effectively a duplicate\n",
                                      j, temp1.mpc_code, buff);
            obs[j - 1] = temp1;
            }
         else
            {
            obs[j] = obs[i];
            j++;
            }
         }
      else
         comment_observation( obs + j - 1, "Duplicate");
   return( j);
}

static int fix_radar_obs( OBSERVE *obs, unsigned n_obs)
{
   unsigned i;

   for( i = 0; n_obs && i < n_obs - 1; i++, obs++)
      if( obs[0].note2 == 'R' && obs[1].note2 == 'R'
                 && obs[0].jd == obs[1].jd)
         {
         unsigned pass;

         assert( obs[0].second_line);
         assert( obs[1].second_line);
         for( pass = 0; pass < 4; pass++)
            {
            unsigned offset = ((pass & 1) ? 33 : 47);
            char *loc0 = obs[0].second_line + (pass >= 2 ? 81 : 0);
            char *loc1 = obs[1].second_line + (pass >= 2 ? 81 : 0);

            if( !memcmp( loc0 + offset, "              ", 14))
               memcpy( loc0 + offset, loc1 + offset, 14);
            }
         free( obs[1].second_line);
         n_obs--;
         memmove( obs + 1, obs + 2, (n_obs - i - 1) * sizeof( OBSERVE));
         }
   return( n_obs);
}

/* MPC gives radar data to six decimal places in days.  But this is  */
/* meant to correspond to the nearest UTC second.  Some legerdemain  */
/* is required to get around this.  The way I'm doing it,  at least  */
/* at present,  is to reformat the date/time to be in HHMMSS form.   */
/*    EDIT:  I decided this wasn't really necessary.  The time _is_  */
/* rounded to the nearest second (see 'extract_date_from_mpc_report' */
/* for details).  I suppose one could argue that,  when the data is  */
/* shown to the user,  it should be explicitly shown in HHMMSS form  */
/* so as to make it blatantly clear what's going on.  But people are */
/* used to the MPC's microday format;  it may be more confusing than */
/* helpful.  So at least for the nonce,  this is #ifdeffed out.      */

#ifdef POSSIBLY_NOT_A_REALLY_GOOD_IDEA
static void fix_radar_time( char *buff)
{
   char tbuff[70];
   const double jd = extract_date_from_mpc_report( buff, NULL);

   full_ctime( tbuff, jd, CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD
                 | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_LEADING_ZEROES);
   buff[15] = (char)( 'A' + (tbuff[0] - '1') * 10 + tbuff[1] - '0');  /* cen */
   buff[16] = tbuff[2];        /* decade */
   buff[17] = tbuff[3];        /* year */
   buff[18] = tbuff[5];        /* month, tens place */
   buff[19] = tbuff[6];        /* month, units */
   buff[20] = tbuff[8];        /* day, tens */
   buff[21] = tbuff[9];        /* day, units */
   buff[22] = ':';             /* signal HHMMSS format */
   buff[23] = tbuff[11];       /* hours   */
   buff[24] = tbuff[12];       /* hours   */
   buff[25] = tbuff[14];       /* minutes */
   buff[26] = tbuff[15];       /* minutes */
   buff[27] = tbuff[17];       /* seconds */
   buff[28] = tbuff[18];       /* seconds */
   buff[29] = buff[30] = buff[31] = ' ';
}
#endif

int unload_observations( OBSERVE FAR *obs, const int n_obs)
{
   int i;

   if( obs)
      {
      for( i = 0; i < n_obs; i++)
         if( obs[i].second_line)
            {
//          debug_printf( "Unloading %d: '%s'\n", i, obs[i].second_line);
            free( obs[i].second_line);
            }
      free( obs);
      }
   return( 0);
}

/* I had ideas of checking for discordant satellite coordinates
by fitting a quadratic to x, y, and z,  and looking for
outliers.  But at least thus far,  it looks like a lot of work
for catching not many errors. */

#ifdef FUTURE_SATELLITE_OBS_CHECKING

/* Evaluates the quadratic whose coefficients are determined by
the subsequent 'compute_quadratic_fits( )' function (q.v.). */

static double compute_quad( const double *coeffs, const double x)
{
   const double dx = x - coeffs[0];

   return( coeffs[3] + dx * (coeffs[2] + coeffs[1] * dx));
}

/* Computes the 'best-fit' quadratic

y = a(x-x0)^2 + b(x-x0) + c

   passing through the points

(x[0], y[0]), (x[1], y[1]), ... (x[n_pts-1], y[n_pts-1]).

x0 is the mean of the input x coordinates.  This (a) reduces roundoff error
in situations where the x values are of similar magnitude and (b) makes the
math slightly easier.  I essentially took this from the 'quadratic curve
fitting' section of Meeus' _Astronomical Algorithms_,  p.43,  except that
with this adjustment (and a temporary similar one in y),  P and T (mean
values of x and y) become zero.  */

static inline void compute_quadratic_fits( double *coeffs, const double *x,
        const double *y, const unsigned n_pts)
{
   double x0 = 0., y0 = 0., sum_x2 = 0., sum_x3 = 0., sum_x4 = 0;
   double quad, linear, constant;
   double sum_xy = 0., sum_x2y = 0.;
   double d;
   unsigned i;

   for( i = 0; i < n_pts; i++)
      {
      x0 += x[i];
      y0 += y[i];
      }
   x0 /= (double)n_pts;    /* compute mean values of x and y */
   y0 /= (double)n_pts;
   for( i = 0; i < n_pts; i++)
      {
      const double dx = x[i] - x0, dy = y[i] - y0, x2 = dx * dx, xy = dx * dy;

      sum_x2 += x2;
      sum_x3 += x2 * dx;
      sum_x4 += x2 * x2;
      sum_xy += xy;
      sum_x2y += xy * dx;
      }
   d = (double)n_pts * (sum_x2 * sum_x4 - sum_x3 * sum_x3)
                  - sum_x2 * sum_x2 * sum_x2;
   quad = (double)n_pts * (sum_x2 * sum_x2y - sum_x3 * sum_xy) / d;
   linear = (double)n_pts * (sum_x4 * sum_xy - sum_x3 * sum_x2y)
                  - sum_x2 * sum_x2 * sum_xy;
   linear /= d;
   constant = (sum_x2 * sum_x3 * sum_xy - sum_x2 * sum_x2 * sum_x2y) / d;
   coeffs[0] = x0;
   coeffs[1] = quad;
   coeffs[2] = linear;
   coeffs[3] = constant + y0;
}

inline bool is_satellite_obs( const OBSERVE *obs)
{
   return( obs->second_line && obs->second_line[14] == 's');
}

static void check_satellite_obs( const OBSERVE *obs, unsigned n_obs)
{
   while( n_obs && !is_satellite_obs( obs))
      {
      obs++;            /* skip ahead to a satellite observation */
      n_obs--;
      }
   if( n_obs)
      {
      unsigned i, n_sat;
      double *arrays[4];

      for( i = n_sat = 0; i < n_obs; i++)
         if( is_satellite_obs( obs + i))
            n_sat++;
      arrays[0] = (double *)calloc( n_sat * 4, sizeof( double));
      assert( arrays[0]);
      for( i = 1; i < 4; i++)
         arrays[i] = arrays[i - 1] + n_sat;
      for( i = n_sat = 0; i < n_obs; i++)
         if( is_satellite_obs( obs + i))
            {
            double xyz[3];

            get_satellite_offset( obs[i].second_line, xyz);
            arrays[0][n_sat] = obs[i].jd;
            arrays[1][n_sat] = xyz[0];
            arrays[2][n_sat] = xyz[1];
            arrays[3][n_sat] = xyz[2];
            n_sat++;
            }

      for( i = 1; i < 4; i++)
         {
         double coeffs[4];
         unsigned j;

         compute_quadratic_fits( coeffs, arrays[0], arrays[i], n_sat);
         debug_printf( "Axis %u\n", i);
         for( j = 0; j < n_sat; j++)
            {
            const double value = arrays[i][j];
            const double quad_value = compute_quad( coeffs, arrays[0][j]);
            const double multiplier = (obs->second_line[32] == '1' ?
                     AU_IN_KM : 1.);

            debug_printf( "   %3u: %15f %15f %15f\n", j,
                        value * multiplier, quad_value * multiplier,
                        (value - quad_value) * multiplier);
            }
         }
      free( arrays[0]);
      }
}
#endif    /* #ifdef FUTURE_SATELLITE_OBS_CHECKING */

/* (691) Spacewatch has a habit of providing near-duplicate observations,
one "hand-measured",  the other not,  like this :

     P10pKck  C2015 12 06.25527 04 35 11.84 +14 22 25.4          21.0 VoNEOCP691
     P10pKck HC2015 12 06.25527 04 35 11.85 +14 22 25.4          21.4 VoNEOCP691

   Normally,  two observations with identical times and observatories would
both be discarded.  In this case,  we want to toss out the first (non-hand
measured) observation only.         */

inline int spacewatch_duplication( OBSERVE FAR *obs)
{
   int rval = 0;

   if( !strcmp( obs->mpc_code, "691") && obs[0].note1 == ' '
                                      && obs[1].note1 == 'H')
      {
      rval = 1;
      obs->is_included = 0;
      obs->flags |= OBS_DONT_USE;
      comment_observation( obs, "Replaced");
      }
   return( rval);
}

int sanity_check_observations = 1;
bool use_sigmas = true;

OBSERVE FAR *load_observations( FILE *ifile, const char *packed_desig,
                           const int n_obs)
{
   const double days_per_year = 365.25;
   char buff[350], mpc_code_from_neocp[4], desig_from_neocp[15];
   char obj_name[80];
   OBSERVE FAR *rval;
   bool including_obs = true;
   int i = 0;
   unsigned line_no = 0;
   unsigned n_below_horizon = 0, n_in_sunlight = 0;
   unsigned lines_actually_read = 0;
   unsigned n_spurious_matches = 0;
   unsigned n_sat_obs_without_offsets = 0;
   bool obs_are_of_a_comet = false;
   double override_posn_sigma_1 = 0.;  /* in arcsec */
   double override_posn_sigma_2 = 0.;
   double override_posn_sigma_theta = 0.;
   double override_mag_sigma = 0.;   /* in mags */
   double override_time_sigma = 0.;  /* in seconds */
            /* We distinguish between observations that are complete clones */
            /* of each other,  and those with the same time, RA/dec, MPC    */
            /* code, and magnitude,  but which differ someplace else.       */
   unsigned n_duplicate_obs_found = 0;
   unsigned n_almost_duplicates_found = 0;
   unsigned n_parse_failures = 0;
   unsigned n_bad_satellite_offsets = 0;
   extern int monte_carlo_object_count;  /* we just want to zero this */
   extern int n_monte_carlo_impactors;   /* and this,  too */
   const bool fixing_trailing_and_leading_spaces =
               (*get_environment_ptr( "FIX_OBSERVATIONS") != '\0');

   *desig_from_neocp = '\0';
   strcpy( mpc_code_from_neocp, "500");   /* default is geocenter */
   neocp_file_type = NEOCP_FILE_TYPE_UNKNOWN;
   get_object_name( obj_name, packed_desig);
   rval = (OBSERVE FAR *)FCALLOC( n_obs + 1, sizeof( OBSERVE));
   if( !rval)
      return( NULL);
   input_coordinate_epoch = 2000.;
   while( fgets_trimmed( buff, sizeof( buff), ifile) && i != n_obs)
      {
      int is_rwo = 0;
      double rwo_posn_sigma_1 = 0., rwo_posn_sigma_2 = 0., rwo_mag_sigma = 0.;
      char original_packed_desig[13];
      size_t ilen = strlen( buff);

      line_no++;
      lines_actually_read++;

      if( *buff == '<')
         remove_html_tags( buff);
      if( !memcmp( buff, "COD ", 4) && ilen == 7)
         strcpy( curr_xxx_code, buff + 4);
      if( !memcmp( buff, "COM Long.", 9))
         get_xxx_location( buff);
      if( debug_level > 2)
         debug_printf( "Line %d: %s\n", line_no, buff);
      if( get_neocp_data( buff, desig_from_neocp, mpc_code_from_neocp))
         if( !i && debug_level)
            debug_printf( "Got NEOCP data\n");
      if( ilen == 75 || ilen >= MINIMUM_RWO_LENGTH)
         {
         is_rwo = rwo_to_mpc( buff, &rval[i].ra_bias, &rval[i].dec_bias,
                &rwo_posn_sigma_1, &rwo_posn_sigma_2, &rwo_mag_sigma);
         if( is_rwo && !i && debug_level)
            debug_printf( "Got .rwo data\n");
         }
      if( ilen >= 92 && buff[82] == '.' && buff[88] == '.')
         {                    /* Possible Tholen-style sigmas:  2F6.3 */
         double sig1, sig2;   /* right after the usual 80 columns */
         int bytes_read;

         if( sscanf( buff + 80, "%lf%lf%n", &sig1, &sig2, &bytes_read) >= 2
                     && bytes_read == 12)
            {
            ilen = 80;
            buff[80] = '\0';
            rwo_posn_sigma_1 = sig1;
            rwo_posn_sigma_2 = sig2;
            }
         }
      if( fixing_trailing_and_leading_spaces)
         fix_up_mpc_observation( buff);
      original_packed_desig[12] = '\0';
      memcpy( original_packed_desig, buff, 12);
      xref_designation( buff);
      if( observation_jd( buff) &&
                     !compare_desigs( packed_desig, buff))
         {
         const int error_code = parse_observation( rval + i, buff);

         strcpy( rval[i].packed_id, original_packed_desig);
         if( error_code)
            {
            n_parse_failures++;
            debug_printf( "Bad obs; error code %d:\n%s\n", error_code, buff);
            }
         else           /* Successfully-loaded observation: */
            {
            bool observation_is_good = true;
            char second_line[81];

            if( buff[14] == 's' || buff[14] == 'S'
                    || buff[14] == 'v' || buff[14] == 'V'
                    || buff[14] == 'r' || buff[14] == 'R')
               {
               if( buff[14] == 'R' && is_rwo)
                  strcpy( second_line, buff + 81);
               else if( !look_for_matching_line( buff, second_line))
                  observation_is_good = false;
               if( observation_is_good)
                  parse_observation( rval + i, buff);
               }

            if( buff[14] == 'S' && observation_is_good)
               {      /* we did find the "matching" line: */
               double vect[3];
               int j, error_code;

               lines_actually_read++;
               rval[i].satellite_obs = (char)(second_line[32] - '0');
               error_code = get_satellite_offset( second_line, vect);
               if( error_code)
                  {
                  rval[i].flags |= OBS_DONT_USE;
                  comment_observation( rval + i, "?off");
                  n_bad_satellite_offsets++;
                  debug_printf( "Error code %d; offending line was:\n%s\n",
                     error_code, second_line);
                  }
               for( j = 0; j < 3; j++)
                  rval[i].obs_posn[j] += vect[j];
               rval[i].second_line = (char *)malloc( 81);
               strcpy( rval[i].second_line, second_line);
               }
            else if( buff[14] == 'R' && observation_is_good)
               {      /* we did find the "matching" line: */
               lines_actually_read++;
//             fix_radar_time( buff);
               rval[i].second_line = (char *)malloc( 81 * 2);
               strcpy( rval[i].second_line, second_line);
               strcpy( rval[i].second_line + 81, buff);
               }
            else if( buff[14] == 'V' && observation_is_good)
               {
               double rho_sin_phi, rho_cos_phi;

               lines_actually_read++;
               xref_designation( second_line);
               roving_lon = atof( second_line + 34);
               roving_lat = atof( second_line + 45);
               roving_ht_in_meters = atof( second_line + 56);
               lat_alt_to_parallax( roving_lat * PI / 180., roving_ht_in_meters,
                                    &rho_cos_phi, &rho_sin_phi, 3);
               compute_observer_loc( rval[i].jd, 3, rho_cos_phi, rho_sin_phi,
                                      roving_lon * PI / 180., rval[i].obs_posn);
               compute_observer_vel( rval[i].jd, 3, rho_cos_phi, rho_sin_phi,
                                      roving_lon * PI / 180., rval[i].obs_vel);
               set_obs_vect( rval + i);
               rval[i].second_line = (char *)malloc( 81);
               strcpy( rval[i].second_line, second_line);
               }
            if( observation_is_good)
               {
               const double radians_per_arcsec = PI / (180. * 3600.);
               double mag_sigma = 0., time_sigma = 0.;
               double posn_sigma_1 = 0.0001, posn_sigma_2 = 0.0001;
               double posn_sigma_theta = 0.;
                    /* If we want all observations to have the same sigma, */
                    /* we use a nonexistent MPC code.  "Unknown" codes will */
                    /* just get the default sigma assigned at the end of    */
                    /* 'sigma.txt' (q.v.)                                   */
               const char *mpc_code_to_use = (use_sigmas ? rval[i].mpc_code : "***");

               posn_sigma_1 = get_observation_sigma( rval[i].jd,
                       (int)( rval[i].obs_mag * 10. + .001),
                       mpc_code_to_use, &mag_sigma, &time_sigma, rval[i].note1);
               posn_sigma_2 = posn_sigma_1;     /* ...for the nonce,  anyway */
                           /* The sigmas just determined may be overruled */
                           /* by keywords given in the observation file : */
               if( use_sigmas)
                  {
                  double ra_sigma, dec_sigma;
                  int bytes_read;

                  if( override_mag_sigma)
                     mag_sigma = override_mag_sigma;
                  if( override_time_sigma)
                     time_sigma = override_time_sigma;
                  if( override_posn_sigma_1)
                     posn_sigma_1 = override_posn_sigma_1;
                  if( override_posn_sigma_2)
                     posn_sigma_2 = override_posn_sigma_2;
                  if( override_posn_sigma_theta)
                     posn_sigma_theta = override_posn_sigma_theta;

                     /* AstDyS/NEODyS .rwo data has position and mag sigmas */
                     /* in it;  Dave Tholen puts position sigmas in columns */
                     /* 81-92.  If available,  we use these :               */
                  if( rwo_posn_sigma_1 && rwo_posn_sigma_2)
                     {
                     posn_sigma_1 = rwo_posn_sigma_1;
                     posn_sigma_2 = rwo_posn_sigma_2;
                     posn_sigma_theta = 0.;
                     }
                  if( rwo_mag_sigma)
                     mag_sigma = rwo_mag_sigma;
                  if( sscanf( rval[i].columns_57_to_65, "%lf %lf%n",
                              &ra_sigma, &dec_sigma, &bytes_read) >= 2
                              && bytes_read >= 8
                              && ra_sigma > 0 && dec_sigma > 0.)
                     {
                     const char end_char = rval[i].columns_57_to_65[8];
                     double units = 0.;

                     if( bytes_read == 9)
                        units = 1.;
                     if( bytes_read == 8)
                        {
                        if( end_char == 'u')
                           units = 1.e-6;
                        if( end_char == 'm')
                           units = 0.001;
                        if( end_char == '\'')
                           units = 60.;
                        if( end_char == 'd')
                           units = 3600.;
                        }
                     if( units)
                        {
                        posn_sigma_1 =  ra_sigma * units;
                        posn_sigma_2 = dec_sigma * units;
                        }
                     posn_sigma_theta = 0.;
                     }
                  }
                           /* The observation data's precision has already been */
                           /* used to figure out a minimum sigma;  e.g.,  a mag */
                           /* of '16.3' results in a sigma of 0.1.  If the above */
                           /* methods return smaller values,  we stick with the */
                           /* sigma determined from the # of places given. */
                           /* If we're enforcing uniform sigmas,  though,  we use */
                           /* those sigmas whether they "make sense" or not. */
               if( mag_sigma > rval[i].mag_sigma || !use_sigmas)
                  rval[i].mag_sigma = mag_sigma;
               if( time_sigma > rval[i].time_sigma || !use_sigmas)
                  rval[i].time_sigma = time_sigma;
               if( posn_sigma_1 > rval[i].posn_sigma_1 || !use_sigmas)
                  rval[i].posn_sigma_1 = posn_sigma_1;
               if( posn_sigma_2 > rval[i].posn_sigma_2 || !use_sigmas)
                  rval[i].posn_sigma_2 = posn_sigma_2;
               rval[i].posn_sigma_theta = posn_sigma_theta;
               if( !including_obs)
                  rval[i].is_included = 0;
               if( apply_debiasing)
                  {
                  rval[i].ra  -= rval[i].ra_bias  * radians_per_arcsec / cos( rval[i].dec);
                  rval[i].dec -= rval[i].dec_bias * radians_per_arcsec;
                  }
               if( rval[i].note2 == 'n')        /* video observations:  assume */
                  rval[i].time_sigma = 0.01 / seconds_per_day;  /* 10ms sigma */
               set_obs_vect( rval + i);
               if( is_rwo && !rval[i].is_included)
                  rval[i].flags |= OBS_DONT_USE;
               if( obs_are_of_a_comet)
                  rval[i].flags |= OBS_IS_COMET;
               i++;
               }
            }
         }
      override_time = 0.;
               /* For backwards compatibility,  we'll handle both 'weight' */
               /* and 'sigma' keywords,  with the former considered to     */
               /* indicate 1/sigma arcseconds.                             */
      if( *buff == '#')
         {
         if( !memcmp( buff, "#Weight ", 8))
            {
            override_posn_sigma_1 = override_posn_sigma_2 = 1. / atof( buff + 8);
            override_posn_sigma_theta = 0.;
            }
         else if( !memcmp( buff, "#Posn sigma ", 12))
            {
            OBSERVE tmp;

            if( set_tholen_style_sigmas( &tmp, buff + 12))
               {
               override_posn_sigma_1 = tmp.posn_sigma_1;
               override_posn_sigma_2 = tmp.posn_sigma_2;
               override_posn_sigma_theta = tmp.posn_sigma_theta;
               }
            }
         else if( !memcmp( buff, "#Mag sigma ", 11))
            override_mag_sigma = atof( buff + 11);
         else if( !memcmp( buff, "#Time sigma ", 12))
            override_time_sigma = atof( buff + 12) / seconds_per_day;
         else if( !memcmp( buff, "#override_weight", 15)
               || !memcmp( buff, "#override_sigma", 14))
            override_posn_sigma_1 = override_posn_sigma_2
            = override_posn_sigma_theta
            = override_time_sigma = override_mag_sigma = 0.;
         else if( !memcmp( buff, "#coord epoch ", 13))
            {
            if( buff[13] == 'a')       /* apparent coords */
               input_coordinate_epoch = -1.;
            else
               input_coordinate_epoch = atof( buff + 13);
            }
         else if( !memcmp( buff, "#suppress_obs", 13))
            including_obs = false;
         else if( !memcmp( buff, "#include_obs", 12))
            including_obs = true;
         else if( !memcmp( buff, "#toffset", 7))
            observation_time_offset = atof( buff + 8) / seconds_per_day;
         else if( !memcmp( buff, "#relax_xyz", 10))
            strict_sat_xyz_format = false;
         else if( !memcmp( buff, "#time ", 6))
            override_time = get_time_from_string( 0, buff + 6,
                              CALENDAR_JULIAN_GREGORIAN, NULL);
                  /* Above allows one to reset the time of the preceding obs */
         else if( !memcmp( buff, "#comet", 6))
            obs_are_of_a_comet = (atoi( buff + 6) != 0);
         }
      }

   n_obs_actually_loaded = i;
   if( debug_level)
      debug_printf( "%u obs found in file\n",  n_obs_actually_loaded);
   for( i = 0; i < n_obs_actually_loaded; i++)
      if( rval[i].note2 == 'X' || rval[i].note2 == 'x')     /* deleted observation */
         {
         comment_observation( rval + i, "Deleted");
         rval[i].flags |= OBS_DONT_USE;
         rval[i].is_included = 0;
         }
   for( i = n_obs_actually_loaded; i < n_obs; i++)
      rval[i].jd = 0.;

   i = sort_obs_by_date_and_remove_duplicates( rval, n_obs_actually_loaded);
   if( rval[i - 1].jd - rval[0].jd > maximum_observation_span * days_per_year)
      {
      sprintf( buff, "The observations span %.1f years,  greater than Find_Orb's\n"
                     "maximum span of %.1f years.  If you really want the whole\n"
                     "set of observations (which may not process correctly),  edit\n"
                     "'environ.dat' and look for the MAX_OBSERVATION_SPAN line.",
                     (rval[i - 1].jd - rval[0].jd) /  days_per_year,
                     maximum_observation_span);
      generic_message_box( buff, "o");
      while( rval[i - 1].jd - rval[0].jd > maximum_observation_span * days_per_year)
         i--;
      }

   n_duplicate_obs_found = n_obs_actually_loaded - i;
   n_obs_actually_loaded = fix_radar_obs( rval, i);

   monte_carlo_object_count = 0;
   n_monte_carlo_impactors = 0;
   if( look_for_matching_line( NULL, NULL))
      {
      strcpy( buff, "Not all satellite observations were read correctly.\n");
      strcat( buff, "This shouldn't happen.  Please send your observation\n");
      strcat( buff, "file to pluto@projectpluto.com so it can be fixed.\n");
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "o");
      }
   if( n_duplicate_obs_found)
      {
      sprintf( buff, "%u observations were duplicates.\n",
                     n_duplicate_obs_found);
      strcat( buff, "The duplicates have been removed.\n");
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "o");
      }
   if( n_bad_satellite_offsets)
      {
      debug_printf( "%u bad sat offsets\n", n_bad_satellite_offsets);
      snprintf( buff, sizeof( buff),
          "%u satellite observations had incorrect offset data.  'debug.txt'\n\
may have further information.  These observations will be excluded.\n",
                     n_bad_satellite_offsets);
      generic_message_box( buff, "o");
      }

   if( n_parse_failures)
      {
      sprintf( buff, "%u observations were incorrectly formatted.\n",
                     n_parse_failures);
      strcat( buff, "These observations will be ignored.\n");
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "o");
      }

   if( sanity_check_observations)
      for( i = 0; i < n_obs_actually_loaded; i++)
         if( rval[i].note2 != 'R')        /* no way to sanity-test */
            {                             /* radar obs yet         */
            DPT obj_alt_az, sun_alt_az;

            if( !get_obs_alt_azzes( rval + i, &sun_alt_az, &obj_alt_az)
                     && sun_alt_az.x > -90.)
               {                          /* I.e., not flagged as meaningless */
               if( sun_alt_az.y > 0.)
                  {
                  n_in_sunlight++;
                  rval[i].is_included = 0;
                  rval[i].flags |= OBS_DONT_USE;
                  comment_observation( rval + i, "Daylit");
                  }
               if( obj_alt_az.y < 0.)
                  {
                  n_below_horizon++;
                  rval[i].is_included = 0;
                  rval[i].flags |= OBS_DONT_USE;
                  comment_observation( rval + i, "Horizon");
                  }
               }
            }
   if( n_below_horizon || n_in_sunlight)
      {
      *buff = '\0';
      if( n_below_horizon)
         sprintf( buff, "%u observations are below the horizon.\n",
                        n_below_horizon);
      if( n_in_sunlight)
         sprintf( buff + strlen( buff),
                        "%u observations were taken in sunlight.\n",
                        n_in_sunlight);
      sprintf( buff + strlen( buff), "These observations will be ignored.\n");
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "o");
      }
   for( i = 1; i < n_obs_actually_loaded; i++)
      if( rval[i].jd == rval[i - 1].jd && !strcmp( rval[i].mpc_code, rval[i - 1].mpc_code))
         if( toupper( rval[i].note2) != 'X' && toupper( rval[i - 1].note2) != 'X')
            if( !spacewatch_duplication( rval + i - 1))
               {
               rval[i].is_included = 0;
               rval[i].flags |= OBS_DONT_USE;
               comment_observation( rval + i, "Duplicate");
               if( rval[i].ra == rval[i - 1].ra
                        && rval[i].dec == rval[i - 1].dec
                        && rval[i].obs_mag == rval[i - 1].obs_mag)
                  n_almost_duplicates_found++;
               else
                  {
                  comment_observation( rval + i - 1, "Duplicate");
                  n_spurious_matches++;
                  rval[i - 1].is_included = 0;
                  rval[i - 1].flags |= OBS_DONT_USE;
                  }
               }
   if( n_spurious_matches)
      {
      sprintf( buff, "%u observations match in date and observatory code,\n",
               n_spurious_matches);
      sprintf( buff + strlen( buff),
                "but not in other regards.  They will be ignored.\n");
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "o");
      }
   if( n_almost_duplicates_found)
      {
      sprintf( buff, "%u observations match in date, RA/dec, magnitude, and MPC\n",
               n_almost_duplicates_found);
      sprintf( buff + strlen( buff),
                "code but not in other regards.  They will be ignored.\n");
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "o");
      }
   for( i = 0; i < n_obs_actually_loaded; i++)
      if( rval[i].flags & OBS_NO_OFFSET)
         n_sat_obs_without_offsets++;
   if( n_sat_obs_without_offsets)
      {
      sprintf( buff, "%u observations are from spacecraft,  but aren't marked\n",
                        n_sat_obs_without_offsets);
      strcat( buff, "as such and lack the 'second line' offset data.  See\n\n");
      strcat( buff, "http://www.minorplanetcenter.net/iau/info/SatelliteObs.html\n");
      strcat( buff, "\nfor info on how to handle spacecraft-based observations.\n");
      generic_message_box( buff, "o");
      }
   if( rval[n_obs_actually_loaded - 1].jd > current_jd( ))
      generic_message_box(
            "The date/time for one or more observations is in the future.\n"
            "If that's not what you expected,  you should check the data.\n", "o");


#ifdef FUTURE_SATELLITE_OBS_CHECKING
   check_satellite_obs( rval, n_obs_actually_loaded);
#endif
   if( n_obs == 2)
      {
      if( rval[0].note2 != 'X' && rval[1].note2 == 'X')
         rval[1].flags |= OBS_DONT_USE;
      if( rval[0].note2 == 'X' && rval[1].note2 != 'X')
         rval[0].flags |= OBS_DONT_USE;
      }
   return( rval);
}

static int id_compare( const void *a, const void *b, void *context)
{
   const OBJECT_INFO *aptr = (const OBJECT_INFO *)a;
   const OBJECT_INFO *bptr = (const OBJECT_INFO *)b;
   int *comparison_method = (int *)context;
   int rval;

   switch( *comparison_method)
      {
      case OBJECT_INFO_COMPARE_PACKED:
         rval = strcmp( aptr->packed_desig, bptr->packed_desig);
         break;
      case OBJECT_INFO_COMPARE_LAST_OBSERVED:
         rval = (aptr->jd_end > bptr->jd_end ? 1 : -1);
         break;
      case OBJECT_INFO_COMPARE_NAME:
      default:
         rval = strcmp( aptr->obj_name, bptr->obj_name);
         break;
      }
   return( rval);
}

void sort_object_info( OBJECT_INFO *ids, const int n_ids,
                                          int object_info_compare_method)
{
   shellsort_r( ids, n_ids, sizeof( OBJECT_INFO), id_compare,
                                    &object_info_compare_method);
}

/* find_objects_in_file( ) reads through the file of MPC astrometric data
   specified by 'filename',  and figures out which objects appear in that
   file.  Those objects can then be listed on the console (findorb) or
   a scroll box (find_o32) or used to figure out what objects are going
   to be loaded (fo),  and the user can select one for which she wants
   to determine an orbit.

   We start out not knowing how many objects may be in the file,  or in what
   order they may appear.  So a hash table with room for 20 objects is
   allocated to start with.  Once that gets to be 75% full,  a new table,
   double in size,  is allocated;  everything is moved from the old table
   to the new one,  and the old one is freed.

   When we're done,  the new table is sorted by name (this puts any blank
   entries at the end of the table).   */

#define UNNEEDED_DEBUGGING_CODE

OBJECT_INFO *find_objects_in_file( const char *filename,
                                         int *n_found, const char *station)
{
   FILE *ifile = fopen( filename, "rb");
   OBJECT_INFO *rval;
   int i, n = 0, n_alloced = 20, prev_loc = -1;
   const int fixing_trailing_and_leading_spaces =
               *get_environment_ptr( "FIX_OBSERVATIONS");
   char buff[250], mpc_code_from_neocp[4], desig_from_neocp[15];
#ifdef CONSOLE
   const clock_t t0 = clock( );
   int next_output = 20000;
   long filesize;

   if( ifile)
      {
      fseek( ifile, 0L, SEEK_END);
      filesize = ftell( ifile);
      fseek( ifile, 0L, SEEK_SET);
      }
#endif

   if( !ifile)
      {
      *n_found = -1;
      return( NULL);
      }
   *mpc_code_from_neocp = '\0';
   *desig_from_neocp = '\0';
   strcpy( mpc_code_from_neocp, "500");   /* default is geocenter */
   neocp_file_type = NEOCP_FILE_TYPE_UNKNOWN;
   rval = (OBJECT_INFO *)calloc( n_alloced + 1, sizeof( OBJECT_INFO));
   if( debug_level > 8)
      debug_printf( "About to read input\n");
   while( fgets_trimmed( buff, sizeof( buff), ifile))
      {
      size_t iline_len = strlen( buff);
      bool is_neocp = false;
      double jd;

      if( debug_level > 8)
         debug_printf( "Input line len %zd\n", strlen( buff));
      if( *buff == '<')
         remove_html_tags( buff);
      if( !n || *mpc_code_from_neocp)
         is_neocp = get_neocp_data( buff, desig_from_neocp,
                                                 mpc_code_from_neocp);
      if( debug_level > 8)
         debug_printf( "After get_neocp_data\n");
      if( iline_len > MINIMUM_RWO_LENGTH)
         rwo_to_mpc( buff, NULL, NULL, NULL, NULL, NULL);
      if( iline_len >= 92 && buff[82] == '.' && buff[88] == '.')
         {
         buff[80] = '\0';     /* probably Dave Tholen-style sigmas */
         iline_len = 80;
         }
      if( fixing_trailing_and_leading_spaces)
         fix_up_mpc_observation( buff);
      if( debug_level > 8)
         debug_printf( "After fixup: %zd\n", strlen( buff));
      if( (jd = observation_jd( buff)) && !is_second_line( buff))
         if( !station || !memcmp( buff + 76, station, 3))
            {
            int loc;


            if( *buff == '#')
               *buff = ' ';           /* handle remarked-out lines,  too */
            xref_designation( buff);
            if( prev_loc >= 0 &&
                          !compare_desigs( rval[prev_loc].packed_desig, buff))
               loc = prev_loc;
            else
               loc = find_in_hash_table( rval, buff, n_alloced);
            prev_loc = loc;
            buff[46] = '\0';
            if( combine_all_observations && n)
               loc = 0;
            else if( !rval[loc].packed_desig[0])   /* it's a new one */
               {
               memcpy( rval[loc].packed_desig, buff, 12);
               rval[loc].packed_desig[12] = '\0';
               get_object_name( rval[loc].obj_name, rval[loc].packed_desig);
               rval[loc].n_obs = 0;
               rval[loc].jd_start = rval[loc].jd_end = jd;
               if( is_neocp)   /* for NEOCP obs,  we need to start at the */
                  rval[loc].file_offset = 0L;   /* beginning of the file  */
               else
                  {
                  rval[loc].file_offset = ftell( ifile) - (long)iline_len
                                       - 100;
                  if( (long)rval[loc].file_offset < 0)
                     rval[loc].file_offset = 0;
                  }
               n++;
#ifdef CONSOLE
               if( n == next_output)
                  {
                  char msg_buff[80];
                  const double t_elapsed =
                           (double)( clock( ) - t0) / (double)CLOCKS_PER_SEC;
                  const double fraction_file_read =
                           (double)ftell( ifile) / (double)filesize;

                  next_output += (int)((double)n / (t_elapsed + 1.));
                  sprintf( msg_buff, "%4.1f%% complete", fraction_file_read * 100.);
                  move_add_nstr( 3, 3, msg_buff, -1);
                  refresh_console( );
                  }
#endif
               }
            rval[loc].n_obs++;
            if( rval[loc].jd_start > jd)
               rval[loc].jd_start = jd;
            if( rval[loc].jd_end < jd)
               rval[loc].jd_end = jd;
            if( n == n_alloced - n_alloced / 4)  /* table is 75% full; */
               {                       /* reallocate & move everything */
               const unsigned new_size = n_alloced * 2 - 1;
               OBJECT_INFO *new_rval = (OBJECT_INFO *)calloc(
                                 new_size, sizeof( OBJECT_INFO));

               for( i = 0; i < n_alloced; i++)
                  if( rval[i].packed_desig[0])
                     {
                     const unsigned new_loc = find_in_hash_table( new_rval,
                                   rval[i].packed_desig, new_size);
                     new_rval[new_loc] = rval[i];
                     }
               free( rval);
               rval = new_rval;
               n_alloced = new_size;
               prev_loc = -1;
               }
            }
      }
   fclose( ifile);
   *n_found = n;
               /* The allocated hash table is,  at most,  80% full,  with */
               /* plenty of empty entries.  Sorting the entries will put  */
               /* the blank entries at the end of the table,  and the     */
               /* 'n' entries that were actually found at the start.      */
   for( i = n = 0; i < n_alloced; i++)
      if( rval[i].packed_desig[0])
         rval[n++] = rval[i];
   assert( n == *n_found);
   sort_object_info( rval, n, OBJECT_INFO_COMPARE_PACKED);
   return( rval);
}

/* put_observer_data_in_text( ) takes a 'station_no' and fills 'buff'
   with a little bit of text about that station,  as found from
   STATIONS.TXT:  bits such as the lat/lon and name of the station.
   In Find_Orb, these details are shown for the station that made
   the currently-selected observation. */

void put_observer_data_in_text( const char FAR *mpc_code, char *buff)
{
   double lon, rho_cos_phi, rho_sin_phi;
   const int planet_idx = get_observer_data( mpc_code, buff,
                                  &lon, &rho_cos_phi, &rho_sin_phi);

   if( planet_idx < 0)
      {
      char tbuff[4];

      FMEMCPY( tbuff, mpc_code, 4);
      sprintf( buff, "No information about station '%s'", tbuff);
      }
   else
      {
      char *name = mpc_station_name( buff);

      memmove( buff, name, strlen( name) + 1);
      if( rho_cos_phi)
         {
         double lat, unused_ht_in_meters;
         const char *output_format = "  (%c%.6f %c%.6f)";

                 /* Cvt parallax data from AU back into earth-axis units: */
         rho_cos_phi /= EARTH_MAJOR_AXIS_IN_AU;
         rho_sin_phi /= EARTH_MAJOR_AXIS_IN_AU;
         parallax_to_lat_alt( rho_cos_phi, rho_sin_phi, &lat,
                           &unused_ht_in_meters, planet_idx);
         snprintf_append( buff, 80, output_format,
                           (lat > 0. ? 'N' : 'S'), fabs( lat) * 180. / PI,
                           (lon > 0. ? 'E' : 'W'), fabs( lon) * 180. / PI);
         if( planet_idx == 3)
            extract_region_data_for_mpc_station( buff + strlen( buff),
                           lat, lon);
         }
      }
}

static const char *environ_dot_dat = "environ.dat";

const char *get_environment_ptr( const char *env_ptr)
{
   static char **edata;
   static size_t n_lines;
   size_t i, i1, step;
   const char *rval = "";
   size_t tbuff_len;
   char tbuff[80];
   int compare;

   if( !env_ptr)
      {
      if( edata)
         free( edata);
      edata = NULL;
      return( rval);
      }

   if( !edata)
      {
      int sort_column = 0;

      edata = load_file_into_memory( environ_dot_dat, NULL);
      if( !edata)       /* fall back on a default version: */
         edata = load_file_into_memory( "environ.def", NULL);
      assert( edata);
                  /* remove blank/commented-out lines : */
      for( i = n_lines = 0; edata[i]; i++)
         if( edata[i][0] > ' ' && edata[i][0] != '#')
            edata[n_lines++] = edata[i];
                  /* sort remaining lines : */
      shellsort_r( edata, n_lines, sizeof( char *),
                     string_compare_for_sort, &sort_column);
      }
   assert( edata);
   strcpy( tbuff, env_ptr);
   strcat( tbuff, "=");
   tbuff_len = strlen( tbuff);
   for( i = (size_t)-1, step = 0x8000; step; step >>= 1)
      if( (i1 = i + step) < n_lines &&
               (compare = memcmp( edata[i1], tbuff, tbuff_len)) <= 0)
         {
         i = i1;
         if( !compare)
            rval = edata[i] + tbuff_len;
         }
// if( !*rval)
//    debug_printf( "Didn't find '%s'\n", env_ptr);
   return( rval);
}

void set_environment_ptr( const char *env_ptr, const char *new_value)
{
   size_t i;
   char **text = load_file_into_memory( environ_dot_dat, NULL);
   FILE *ofile;
   bool found_it = false;
   const size_t env_ptr_len = strlen( env_ptr);

   assert( env_ptr);
   assert( new_value);
   if( !text)
      text = load_file_into_memory( "environ.def", NULL);
   assert( text);
   ofile = fopen_ext( environ_dot_dat, "fcwb");
   assert( ofile);
   for( i = 0; text[i]; i++)
      if( !memcmp( text[i], env_ptr, env_ptr_len) && text[i][env_ptr_len] == '=')
         {
         found_it = true;
         fprintf( ofile, "%s=%s\n", env_ptr, new_value);
         }
      else
         fprintf( ofile, "%s\n", text[i]);

   if( !found_it)       /* not a replacement,  but a new value: */
      fprintf( ofile, "%s=%s\n", env_ptr, new_value);
   fclose( ofile);
   free( text);
   get_environment_ptr( NULL);      /* force a reload of the env data */
}

static inline void compute_relative_velocity_vectors( const OBSERVE FAR *obs,
                                    double *vel)
{
   double j2000_vel[3], matrix[9], length;
   int i;

   for( i = 0; i < 3; i++)
      {
      j2000_vel[i] = obs->obj_vel[i] - obs->obs_vel[i];
//    matrix[i] = obs->vect[i];
      matrix[i] = (obs->obj_posn[i] - obs->obs_posn[i]) / obs->r;
      }
   ecliptic_to_equatorial( j2000_vel);
   ecliptic_to_equatorial( matrix);
   length = sqrt( matrix[0] * matrix[0] + matrix[1] * matrix[1]);
   matrix[3] =  matrix[1] / length;
   matrix[4] = -matrix[0] / length;
   matrix[5] = 0.;

   matrix[6] =  matrix[4] * matrix[2];
   matrix[7] = -matrix[3] * matrix[2];
   matrix[8] = length;

            /* Now we've got an orthonormal matrix,  matrix[012] pointing */
            /* in the direction of the observation,  matrix[345] at right */
            /* angles in the equatorial plane,  matrix[678] at right angles */
            /* to both.  So we can multiply: */
   for( i = 0; i < 9; i += 3)
      vel[i / 3] = matrix[  i  ] * j2000_vel[0]
                 + matrix[i + 1] * j2000_vel[1]
                 + matrix[i + 2] * j2000_vel[2];
}

static int get_year_and_mpc_half_month_letter( const double jd, char *letter)
{
   char buff[30];
   int half_month;

   full_ctime( buff, jd, FULL_CTIME_YMD | FULL_CTIME_LEADING_ZEROES
                  | FULL_CTIME_DATE_ONLY | FULL_CTIME_MONTHS_AS_DIGITS);
   half_month = atoi( buff + 5) * 2 - 2;
   if( atoi( buff + 8) >= 16)
      half_month++;
   *letter = (char)( 'A' + half_month);
   if( *letter >= 'I')
      *letter = (char)( *letter + 1);
   return( atoi( buff));
}

/* MPC references are stored in five characters in each 80-byte record,
   in columns 73-77,  in an highly packed manner described at

   http://www.minorplanetcenter.net/iau/info/References.html

   The following function converts those five bytes to human-readable form.

   MPC references are stored as five-digit numbers.  At the current (2012)
rate,  this should work until sometime around 2018,  when we'll have the
MPC 100K Bug to fix.

   MPS references are stored as 'a...z' plus four digits,  allowing 260K
references.  To get beyond that,  MPC used a tilde (~) followed by four
"mutant hex" digits (base 62),  to get us 260000+62^4 = 15 036 336
references.  That oughta last us at least a little while.  I assume
something similar will be done for MPC references when the time comes.

   MPEC references lack a year;  that's why they are often shown in the form
'MPEC ????-whatever'.  But if the observation was made in the last year, we
_can_ replace those four question marks with the correct year.  (In some
cases,  almost two years;  for example,  if today is 2012 Jan 19,  and the
MPEC is ????-C13,  and the observation was made on 2010 Feb 25,  then the
year must be 2011.  It couldn't be 2010,  because late February would
result in a D or later half-month specification.  And it couldn't be 2012,
because it would then have a B or A half-month specification.)   That's why
the observation JD is passed in. */

static void reference_to_text( char *obuff, const char *reference,
                                            const double jd)
{
   if( !strcmp( reference, "     "))       /* no reference given */
      *obuff = '\0';
   else if( !strcmp( reference, "neocp"))
      strcpy( obuff, "NEOCP");
   else if( *reference >= '0' && *reference <= '9')
      sprintf( obuff, "MPC %s", reference);
   else if( *reference >= 'a' && *reference <= 'z')
      sprintf( obuff, "MPS %d%s", *reference - 'a', reference + 1);
   else if( *reference == 'E')
      {
      int obs_year, curr_year;
      char obs_letter, curr_letter;

      sprintf( obuff, "MPEC ?  ?-%c%d", reference[1], atoi( reference + 2));
      obuff[6] = obuff[7] = '?';    /* attempt to evade trigraph oddities */
      obs_year = get_year_and_mpc_half_month_letter( jd, &obs_letter);
      curr_year = get_year_and_mpc_half_month_letter( current_jd( ), &curr_letter);
      if( curr_letter < reference[1])  /* observation is from last year or earlier */
         curr_year--;
      if( obs_letter > reference[1])
         obs_year++;
      if( curr_year == obs_year)    /* this reference can only mean one year: */
         {
         sprintf( obuff + 5, "%4d", curr_year);
         obuff[9] = '-';
         }
      }
   else if( *reference == 'D' && isdigit( reference[1]))
      sprintf( obuff, "DASO %d", atoi( reference + 1));
   else if( *reference == '~')      /* MPS number, packed as four digits */
      {                             /* in base 62 "mutant hex"           */
      unsigned i, mps_number = 0;

      for( i = 1; i < 5; i++)
         {
         mps_number *= 62;
         mps_number += mutant_hex_char_to_int( reference[i]);
         }
      sprintf( obuff, "MPS %u", 260000 + mps_number);
      }
   else           /* just copy it in,  but add a space */
      {
      while( *reference >= 'A')
         *obuff++ = *reference++;
      *obuff++ = ' ';
      while( *reference == '0')
         reference++;
      strcpy( obuff, reference);
      }
}

int compute_observation_motion_details( const OBSERVE FAR *obs,
               MOTION_DETAILS *m)
{
   double vel[3];

   compute_relative_velocity_vectors( obs, vel);
   if( !obs->r)
      m->ra_motion = m->dec_motion = m->position_angle_of_motion = 0.;
   else
      {
      m->ra_motion = vel[1] * (180. / PI) / obs->r;
      m->dec_motion = vel[2] * (180. / PI) / obs->r;
      m->ra_motion *= -60. / hours_per_day;       /* cvt to arcmin/hr, or */
      m->dec_motion *= 60. / hours_per_day;       /* arcsec/minute        */
      m->position_angle_of_motion =
                  180. + (180. / PI) * atan2( -m->ra_motion, -m->dec_motion);
      }
   m->total_motion = sqrt( m->ra_motion * m->ra_motion
                                        + m->dec_motion * m->dec_motion);
   m->xresid = (obs->ra - obs->computed_ra) * cos( obs->dec);
   m->yresid = obs->dec - obs->computed_dec;
               /* cvt xresid, yresid from radians to arcseconds: */
   m->xresid *= (180. / PI) * 3600.;
   m->yresid *= (180. / PI) * 3600.;
               /* time residual is in seconds */
   m->time_residual = m->xresid * m->ra_motion + m->yresid * m->dec_motion;
   m->time_residual *= 60. / (m->total_motion * m->total_motion);
   m->cross_residual = m->xresid * m->dec_motion - m->yresid * m->ra_motion;
   m->cross_residual /= m->total_motion;
   m->radial_vel = vel[0] * AU_IN_KM / seconds_per_day;
   return( 0);
}

         /* motion is in arcminutes/hour */
static void format_motion( char *obuff, const double motion)
{
#ifdef _WIN32
   static const char degree_symbol = (char)0xb0;
#else
   static const char degree_symbol = (char)0xf8;
#endif
   const double fabs_motion = fabs( motion);

   if( fabs_motion < 99.)
      sprintf( obuff, "%5.2f'/hr", motion);
   else if( fabs_motion < 999.)
      sprintf( obuff, "%5.1f'/hr", motion);
   else if( fabs_motion < 99999.)
      sprintf( obuff, "%5.0f'/hr", motion);
   else if( fabs_motion < 99999. * 60.)
      sprintf( obuff, "%5.0f%c/hr", motion / 60., degree_symbol);
   else if( fabs_motion < 99999. * 3600.)
      sprintf( obuff, "%5.0f%c/min", motion / 3600., degree_symbol);
   else if( fabs_motion < 99999. * 216000.)
      sprintf( obuff, "%5.0f%c/sec", motion / 216000., degree_symbol);
   else
      strcpy( obuff, "!!!!!");
}

static double relative_velocity( const double *loc1, const double *vel1,
                                 const double *loc2, const double *vel2)
{
   double dot_prod = 0, dist2 = 0., delta, rval;
   int i;

   for( i = 0; i < 3; i++)
      {
      delta = loc1[i] - loc2[i];
      dist2 += delta * delta;
      dot_prod += delta * (vel1[i] - vel2[i]);
      }
   rval = dot_prod / sqrt( dist2);
         /* Convert AU/day to km/s : */
   return( rval * AU_IN_KM / seconds_per_day);
}

/* MPC 80-column radar observations contain four floating-point values
   that are expressed in fifteen bytes.  Implicitly,  there's a decimal
   point between the eleventh and twelfth bytes.  The first byte may be
   a +/- sign (for the Doppler shift),  or an 'S' or 'C' (for the
   uncertainty in time delay).  There may also be trailing spaces,
   or the whole thing may be blank.  See
   http://www.minorplanetcenter.net/iau/info/RadarObs.html */

static double extract_radar_value( const char *ibuff)
{
   char buff[16];
   int i;
   double rval;

   memcpy( buff, ibuff, 15);
   buff[15] = '\0';
   for( i = 11; i < 15; i++)
      if( buff[i] == ' ')
         buff[i] = '0';
   if( buff[0] >= '0' && buff[0] <= '9')
      rval = atof( buff);
   else
      rval = atof( buff + 1);
   if( buff[0] == '-')
      rval *= -1;
   return( rval * 1.e-4);
}

int compute_radar_info( const OBSERVE *obs, RADAR_INFO *rinfo)
{
   double time_diff = obs->r * AU_IN_KM / SPEED_OF_LIGHT;
   const double jd = obs->jd - time_diff / seconds_per_day;
   double rho_cos_phi, rho_sin_phi, longitude;
   double xyz[3], vel[3];
   double v1, v2, doppler_factor;
   char tbuff[20];
   char *first_line = obs->second_line + 81;
   int iter;

   assert( obs);
   assert( obs->second_line);
   memcpy( tbuff, obs->second_line + 68, 3);
   tbuff[3] = '\0';
   get_observer_data( tbuff, NULL, &longitude, &rho_cos_phi, &rho_sin_phi);
   for( iter = 0; iter < 3; iter++)
      {
      int i;
      double delta[3];

      compute_observer_loc( jd - time_diff / seconds_per_day, 3,
               rho_cos_phi, rho_sin_phi, longitude, xyz);
      for( i = 0; i < 3; i++)
         delta[i] = xyz[i] - obs->obj_posn[i];
      time_diff = vector3_length( delta) * AU_IN_KM / SPEED_OF_LIGHT;
      }
   compute_observer_vel( jd - time_diff / seconds_per_day, 3,
               rho_cos_phi, rho_sin_phi, longitude, vel);
   time_diff += obs->r * AU_IN_KM / SPEED_OF_LIGHT;
   v1 = relative_velocity( obs->obj_posn, obs->obj_vel,
                           obs->obs_posn, obs->obs_vel);
   v2 = relative_velocity( obs->obj_posn, obs->obj_vel,
                                 xyz,           vel);
   v1 /= SPEED_OF_LIGHT;
   v2 /= SPEED_OF_LIGHT;
   doppler_factor = sqrt( (1. - v1) * (1. - v2) / ((1. + v1) * (1. + v2)));
            /* Transmitter frequency is split up with six digits in columns */
            /* 63-68 of the first record,  six more in same columns rec 2,  */
            /* with an implicit decimal after the first five digits.  And   */
            /* maybe nothing but spaces after the decimal point.            */
   memcpy( tbuff, first_line + 62, 5);
   tbuff[5] = '.';
   tbuff[6] = first_line[67];
   memcpy( tbuff + 7, obs->second_line + 62, 8);
   tbuff[15] = '\0';
   rinfo->freq_hz = atof( tbuff) * 1e+6;   /* cvt MHz to Hz */
   rinfo->rtt_obs = extract_radar_value( first_line + 32) * 1.e-6;
   rinfo->rtt_sigma = extract_radar_value( obs->second_line + 32) * 1e-6;
   rinfo->doppler_obs = extract_radar_value( first_line + 47);
   rinfo->doppler_sigma = extract_radar_value( obs->second_line + 47);
   rinfo->doppler_comp = (doppler_factor - 1.) * rinfo->freq_hz;
   rinfo->rtt_comp = time_diff;
   return( 0);
}

static void show_radar_info( char *buff, const OBSERVE *obs)
{
   RADAR_INFO rinfo;

   compute_radar_info( obs, &rinfo);
   sprintf( buff, "RTDist (C) %.8fs = %.3f km; Dopp %.8f km/s = %.2f Hz",
               rinfo.rtt_comp, rinfo.rtt_comp * SPEED_OF_LIGHT,
               rinfo.doppler_comp * SPEED_OF_LIGHT / rinfo.freq_hz,
               rinfo.doppler_comp);
}

static size_t strip_trailing_zeroes( char *buff)
{
   size_t i = strlen( buff);

   while( i && buff[i - 1] == '0')
      i--;
   if( i && buff[i - 1] == '.')
      i--;
   buff[i] = '\0';
   return( i);
}

/*  Generates text such as:

line 0: Elong 167.3    Phase   6.6    RA vel -0.88'/hr   decvel -0.37'/hr   dT=8.15 sec
   or.. Doppler -2.7182818 km/s = 4578.23 Hz
line 1: ang vel  0.95'/hr at PA 247.0   radial vel -7.716 km/s  cross -0.18
   or.. round trip 3141592.653 km = 10.4792251s
line 2: Delta= .92837  r= 1.9137  mag=17.40  mag (computed)=17.76   1997 Oct 15 11:48:11.52
line 3: Obj alt 4.9 az 271.2  Sun alt -17.4 az 89.1
line 4: (709) W & B Observatory, Cloudcroft  (N32.95580 E254.22882)
*/

static int generate_observation_text( const OBSERVE FAR *optr,
                       const int line_number, char *buff)
{
   const double earth_sun = vector3_length( optr->obs_posn);

   *buff = '\0';
   switch( line_number)
      {
      case 0:
         if( optr->note2 == 'R')
            {
            RADAR_INFO rinfo;

            compute_radar_info( optr, &rinfo);
            if( rinfo.rtt_obs)
               {
               sprintf( buff, "Time (obs) %.7f", rinfo.rtt_obs);
               buff += strip_trailing_zeroes( buff);
               sprintf( buff, " +/- %.1fus", rinfo.rtt_sigma * 1e+6);
               buff += strip_trailing_zeroes( buff);
               }
            if( rinfo.doppler_obs)
               {
               if( rinfo.rtt_obs)
                  {
                  strcpy( buff, "   ");
                  buff += 3;
                  }
               sprintf( buff, "Shift(obs) %f", rinfo.doppler_obs);
               buff += strip_trailing_zeroes( buff);
               sprintf( buff, " +/- %f", rinfo.doppler_sigma);
               buff += strip_trailing_zeroes( buff);
               strcpy( buff, " Hz");
//             buff += 3;
               }
            else
               {
               sprintf( buff, "  Dist (comp) %.9f = %.2f km", optr->r,
                           optr->r * AU_IN_KM);
//             buff += strlen( buff);
               }
//          strcpy( buff, optr->second_line);
            }
         else
            {
            const double cos_elong = (earth_sun * earth_sun + optr->r * optr->r
                                       - optr->solar_r * optr->solar_r)
                                       / (2. * earth_sun * optr->r);
            const double cos_phase = (optr->r * optr->r + optr->solar_r * optr->solar_r
                                       - earth_sun * earth_sun)
                                       / (2. * optr->solar_r * optr->r);
            MOTION_DETAILS m;
            char ra_motion_buff[15], dec_motion_buff[15];

            compute_observation_motion_details( optr, &m);
            sprintf( buff, "Elong %5.1f    Phase %5.1f    ",
                                        acose( cos_elong) * 180. / PI,
                                        acose( cos_phase) * 180. / PI);
            format_motion( ra_motion_buff, m.ra_motion);
            format_motion( dec_motion_buff, m.dec_motion);
            sprintf( buff + strlen( buff), "RA vel %s   decvel %s   dT=",
                                           ra_motion_buff, dec_motion_buff);
            buff += strlen( buff);

            if( fabs( m.time_residual) < .999)
               {
               sprintf( buff, "%.3f sec", fabs( m.time_residual));
               *buff = (m.time_residual > 0. ? '+' : '-');
               }
            else if( fabs( m.time_residual) < 99.9)
               sprintf( buff, "%.2f sec", m.time_residual);
            else if( fabs( m.time_residual / 60.) < 99.9)
               sprintf( buff, "%.2f min", m.time_residual / 60.);
            else if( fabs( m.time_residual / 60.) < 9999.)
               sprintf( buff, "%d min", (int)( m.time_residual / 60.));
            else if( fabs( m.time_residual / 3600.) < 9999.)
               sprintf( buff, "%d hr", (int)( m.time_residual / 3600.));
            else
               strcpy( buff, "!!!!");
            }
         break;
      case 1:
         if( optr->note2 == 'R')
            show_radar_info( buff, optr);
         else
            {
            MOTION_DETAILS m;
            char tbuff[15];
            double tdiff;

            compute_observation_motion_details( optr, &m);
            format_motion( tbuff, m.total_motion);
            sprintf( buff, "ang vel %s at PA %.1f", tbuff,
                      m.position_angle_of_motion);
            sprintf( buff + strlen( buff), "   radial vel %.3f km/s  cross ",
                                      m.radial_vel);
            if( fabs( m.cross_residual) < 9.9)
               sprintf( tbuff, "%.2f", m.cross_residual);
            else if( fabs( m.cross_residual) < 99.9)
               sprintf( tbuff, "%4.1f", m.cross_residual);
            else if( fabs( m.cross_residual) < 9999.)
               sprintf( tbuff, "%4d", (int)m.cross_residual);
            else
               strcpy( tbuff, "!!!!");
            strcat( buff, tbuff);
            tdiff = current_jd() - optr->jd;
            if( fabs( tdiff) < 1. / hours_per_day)     /* less than an hour ago */
               sprintf( tbuff, "%d min", (int)( tdiff * minutes_per_day));
            else if( fabs( tdiff) < 1.)
               sprintf( tbuff, "%.1f hr", tdiff * hours_per_day);
            else if( fabs( tdiff) < 100.)
               sprintf( tbuff, "%.1f days", tdiff);
            else
               *tbuff = '\0';
            if( *tbuff)
               sprintf( buff + strlen( buff), "  %s ago", tbuff);
            if( tdiff < 0.)
               strcat( buff, " <FUTURE!>");
            }
         break;
      case 2:
         {
         strcpy( buff, "Delta=");
         buff += strlen( buff);
         format_dist_in_buff( buff, optr->r);  /* ephem0.cpp */

         strcat( buff, "  r=");
         buff += strlen( buff);
         format_dist_in_buff( buff, optr->solar_r);  /* ephem0.cpp */
         strcat( buff, "  ");
         if( optr->obs_mag < BLANK_MAG)
            sprintf( buff + strlen( buff), "mag=%5.2f  ", optr->obs_mag);
         else
            strcat( buff, "           ");
         if( optr->computed_mag)
            sprintf( buff + strlen( buff), "mag (computed)=%5.2f   ",
                      optr->computed_mag);

         full_ctime( buff + strlen( buff),
                           optr->jd - td_minus_utc( optr->jd) / seconds_per_day,
                           FULL_CTIME_HUNDREDTH_SEC | FULL_CTIME_YMD
                            | CALENDAR_JULIAN_GREGORIAN);
         }
         break;
      case 3:
         {
         DPT sun_alt_az, object_alt_az;
         int i;
         static const char *net_codes[] = {
           /* http://www.minorplanetcenter.net/iau/info/CatalogueCodes.html
                G. V. Williams, 2012, ``Minor Planet Astrophotometry'', PhD
                thesis, Open University. [2012PhDT.........7W] */
                  "aUSNO-A1",
                  "bUSNO-SA1",
                  "cUSNO-A2",
                  "dUSNO-SA2",
                  "eUCAC-1",
                  "fTycho-1",
                  "gTycho-2",
                  "hGSC-1.0",
                  "iGSC-1.1",
                  "jGSC-1.2",
                  "kGSC-2.2",
                  "lACT",
                  "mGSC-ACT",
                  "nSDSS-DR8",       /* was TRC */
                  "oUSNO-B1",
                  "pPPM",
                  "qUCAC-4",     /* also UCAC2-beta for earlier obs */
                  "rUCAC-2",
                  "sUSNO-B2",
                  "tPPMXL",
                  "uUCAC-3",
                  "vNOMAD",
                  "wCMC-14",
                  "xHIP-2",
                  "yHIP",
                  "zGSC-1.x",
                  "AAC",
                  "BSAO 1984",
                  "CSAO",
                  "DAGK 3",
                  "EFK4",
                  "FACRS",
                  "GLick Gaspra Catalogue",
                  "HIda93 Catalogue",
                  "IPerth 70",
                  "JCOSMOS/UKST Southern Sky Catalogue",
                  "KYale",
                  "L2MASS", /* used for WISE & PanSTARRS astrometry */
                  "MGSC-2.3",
                  "NSDSS-DR7",
                  "OSST-RC1",
                  "PMPOSC3",
                  "QCMC-15",
                  "RSST-RC4",
                  "SURAT-1",
                  "TURAT-2",
                  "UGAIA-DR1",
                  "VGAIA-DR2",
                  NULL };

         if( optr->posn_sigma_1 > 1.01 || optr->posn_sigma_1 < .99
                  || optr->posn_sigma_1 != optr->posn_sigma_2)
            if( optr->note2 != 'R')
               {
               int tilt_angle = 0;

               strcpy( buff, "Sigma ");
               sprintf( buff + 6, "%.6f", optr->posn_sigma_1);
               remove_insignificant_digits( buff + 6);
               if( optr->posn_sigma_1 != optr->posn_sigma_2)
                  {
                  char *end_ptr = buff + strlen( buff);

                  sprintf( end_ptr, "x%.6f", optr->posn_sigma_2);
                  remove_insignificant_digits( end_ptr);
                  tilt_angle = (int)( optr->posn_sigma_theta * 180. / PI);
                  }
               strcat( buff, "\" ");
               if( tilt_angle)
                  sprintf( buff + strlen( buff), "%d ", tilt_angle);
               }
         buff += strlen( buff);
         reference_to_text( buff, optr->reference, optr->jd);
         if( *buff)
            strcat( buff, "  ");
         if( !get_obs_alt_azzes( optr, &sun_alt_az, &object_alt_az))
            {
            sprintf( buff + strlen( buff), "Obj alt %.1f",
                      object_alt_az.y);
            if( object_alt_az.x > -1.)
               sprintf( buff + strlen( buff), " az %.1f",
                      object_alt_az.x);
            if( optr->note2 == 'R')
               {
               double xresid, yresid;

               get_residual_data( optr, &xresid, &yresid);
               buff += strlen( buff);
               if( xresid && yresid)
                  sprintf( buff, "   %.2f,%.2f sigmas",
                           xresid, yresid);
               else
                  sprintf( buff, "   %.2f sigmas",
                           xresid + yresid);
               }
            else
               {
               sprintf( buff + strlen( buff), "  Sun alt %.1f",
                     sun_alt_az.y);
               if( sun_alt_az.x > -1.)
                  sprintf( buff + strlen( buff), " az %.1f",
                        sun_alt_az.x);
               }
#ifdef TEST_SPACECRAFT_LOC
            if( object_alt_az.x <= -1.)
               {        /* flagged as spacecraft obs */
               double xyz[3], dist, ra, dec;

               get_satellite_offset( optr->second_line, xyz);
               ecliptic_to_equatorial( xyz);
               ra = atan2( xyz[1], xyz[0]) + PI;
               dist = vector3_length( xyz);
               dec = asin( xyz[2] / dist);
               sprintf( buff + strlen( buff),
                          " PrimRA %.3f dec %.3f dist %.1f km",
                          ra * 180. / PI,
                          dec * 180. / PI, dist * AU_IN_KM);
               }
#endif
            }
         for( i = 0; net_codes[i]; i++)
            if( optr->mag_band2 == net_codes[i][0])
               {
               strcat( buff, "  ");
               strcat( buff, net_codes[i] + 1);
               }
         }
         break;
      case 4:
         if( optr->obs_mag < BLANK_MAG)
            sprintf( buff, "Mag sigma %g; ", optr->mag_sigma);
         else
            *buff = '\0';
         sprintf( buff + strlen( buff), "time sigma %g",
                          optr->time_sigma * seconds_per_day);
         if( optr->ra_bias || optr->dec_bias)
             sprintf( buff + strlen( buff), "   RA bias %.3f\" dec bias %.3f\"",
                           optr->ra_bias, optr->dec_bias);
         break;
      case 5:
         put_observer_data_in_text( optr->mpc_code, buff);
         break;
     }
   return( 0);
}

void add_version_and_de_text( char *buff)
{
   strcpy( buff, "Version ");
   strcat( buff, __DATE__);
   format_jpl_ephemeris_info( buff + strlen( buff));
}

int generate_obs_text( const OBSERVE FAR *obs, const int n_obs, char *buff)
{
   size_t i, n_selected = 0, first = -1, last = 0;
   int n_lines = 1;
   char *tptr = buff;

   for( i = 0; i < (size_t)n_obs; i++)
      if( obs[i].flags & OBS_IS_SELECTED)
         {
         if( !n_selected)
            first = i;
         last = i;
         n_selected++;
         }
   if( !n_selected)
      strcpy( buff, "(No observations selected)\n");
   else if( n_selected > 1)
      {
      double mean_xresid = 0., mean_yresid = 0.;
      double mean_xresid2 = 0., mean_yresid2 = 0.;

      sprintf( buff, "%d observations selected of %d\n", (int)n_selected, n_obs);
      for( i = 0; i < (size_t)n_obs; i++)
         if( obs[i].flags & OBS_IS_SELECTED)
            {
            MOTION_DETAILS m;

            compute_observation_motion_details( obs + i, &m);
            mean_xresid += m.xresid;
            mean_yresid += m.yresid;
            mean_xresid2 += m.xresid * m.xresid;
            mean_yresid2 += m.yresid * m.yresid;
            }
      mean_xresid /= (double)n_selected;
      mean_yresid /= (double)n_selected;
      mean_xresid2 /= (double)n_selected;
      mean_yresid2 /= (double)n_selected;
      sprintf( buff + strlen( buff), "Mean RA residual %.2f +/- %.2f; dec %.2f +/- %.2f\n",
             mean_xresid, sqrt( mean_xresid2 - mean_xresid * mean_xresid),
             mean_yresid, sqrt( mean_yresid2 - mean_yresid * mean_yresid));
      n_lines++;
      if( n_selected == 2)
         {
         double dist, posn_ang, delta_time;
         const OBSERVE FAR *optr1 = obs + first;
         const OBSERVE FAR *optr2 = obs + last;

         calc_dist_and_posn_ang( &optr1->ra, &optr2->ra, &dist, &posn_ang);
         dist *= 180. / PI;      /* cvt radians to degrees */
         delta_time = optr2->jd - optr1->jd;

         sprintf( buff + strlen( buff), "Observations are %.1f\" = %.2f' = %.3f degrees apart\n",
               dist * 3600., dist * 60., dist);
         if( fabs( delta_time) < 1.)
            sprintf( buff + strlen( buff), "Time diff: %.1f sec = %.2f min = %.3f hrs\n",
                     delta_time * seconds_per_day,
                     delta_time * minutes_per_day,
                     delta_time * hours_per_day);
         else
            sprintf( buff + strlen( buff), "Time diff: %.1f hrs = %.2f days\n",
                     delta_time * 24., delta_time);
         dist /= delta_time;     /* get motion in degrees/day */
         dist *= 60. / 24.;      /* then convert to '/hr */
                           /* Dunno how the PA got flipped,  but it did: */
         posn_ang = 2. * PI - posn_ang;
         sprintf( buff + strlen( buff), "Motion: %.1f'/hr in RA, %.1f'/hr in dec",
                  dist * sin( posn_ang), dist * cos( posn_ang));
         sprintf( buff + strlen( buff), " (total %.1f'/hr at PA %.1f)\n",
                  dist, posn_ang * 180. / PI);
         n_lines += 3;
         }
      make_date_range_text( buff + strlen( buff),
                                   obs[first].jd, obs[last].jd);
      strcat( buff, "\n");
      n_lines++;
      }
   else
      {
      n_lines = 0;
#ifdef _MSC_VER
      for( i = 0; i < 6; i++)
#else
      for( i = 0; i < 5; i++)
#endif
         {
         generate_observation_text( obs + first, (int)i, tptr);
         if( *tptr)
            {
            strcat( tptr, "\n");
            tptr += strlen( tptr);
            n_lines++;
            }
         }
#ifdef _MSC_VER
      recreate_observation_line( tptr, obs + first);
      strcat( tptr, "\n");
      n_lines++;
#endif
      }
   if( n_lines <= 4)    /* got room for version info */
      {
      add_version_and_de_text( buff + strlen( buff));
      n_lines++;
      }
   return( n_lines);
}

int sanity_test_observations( const char *filename)
{
   FILE *ifile = fopen( filename, "rb");
   FILE *ofile = NULL;
   long line_no = 0L;
   char buff[250];
   OBSERVE obs;
   int n_problems_found = 0;

   if( !ifile)
      return( -1);
   while( fgets_trimmed( buff, sizeof( buff), ifile))
      {
      line_no++;
      if( observation_jd( buff) && !is_second_line( buff))
         if( !parse_observation( &obs, buff))
            {
            DPT alt_az_sun, alt_az_obj;

            if( !get_obs_alt_azzes( &obs, &alt_az_sun, &alt_az_obj))
               if( alt_az_sun.y > 0. || alt_az_obj.y < 0.)
                  {
                  char tbuff[100];

                  if( !ofile)
                     ofile = fopen_ext( "sanity.txt", "fwb");
                  sprintf( tbuff, "Line %ld: Sun alt %.1f az %.1f; obj alt %.1f az %.1f\n",
                           line_no,
                           alt_az_sun.y, alt_az_sun.x,
                           alt_az_obj.y, alt_az_obj.x);
                  printf( "%s %s", tbuff, buff);
                  fprintf( ofile, "%s %s", tbuff, buff);
                  n_problems_found++;
                  }
            override_time = 0.;
            }
      }
   fclose( ifile);
   if( ofile)
      fclose( ofile);
   return( n_problems_found);
}

