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
#endif
#include <stdarg.h>
#include <assert.h>
#include <errno.h>
#include "watdefs.h"
#include "details.h"
#include "comets.h"
#include "lunar.h"
#include "afuncs.h"
#include "mpc_func.h"
#include "mpc_obs.h"
#include "stackall.h"
#include "stringex.h"
#include "sigma.h"
#include "date.h"
#include "pl_cache.h"
#include "constant.h"

int pattern_match(const char* pattern, const char* string);   /* miscell.c */
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
double utc_from_td( const double jdt, double *delta_t);     /* ephem0.cpp */
int apply_excluded_observations_file( OBSERVE *obs, const int n_obs);
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
void remove_trailing_cr_lf( char *buff);            /* ephem0.cpp */
void format_dist_in_buff( char *buff, const double dist_in_au); /* ephem0.c */
double current_jd( void);                       /* elem_out.cpp */
void remove_insignificant_digits( char *tbuff);          /* monte0.c */
int compute_observer_loc( const double jde, const int planet_no,
             const double rho_cos_phi,           /* mpc_obs.cpp */
             const double rho_sin_phi, const double lon, double FAR *offset);
int compute_observer_vel( const double jde, const int planet_no,
             const double rho_cos_phi,           /* mpc_obs.cpp */
             const double rho_sin_phi, const double lon, double FAR *vel);
int get_satellite_offset( const char *iline, double *xyz);  /* mpc_obs.cpp */
int get_residual_data( const OBSERVE *obs, double *xresid, double *yresid);
bool nighttime_only( const char *mpc_code);                 /* mpc_obs.cpp */
char *find_numbered_mp_info( const int number);             /* mpc_obs.cpp */
bool is_sungrazing_comet( const OBSERVE *obs, const int n_obs);  /* orb_func.c */
static int xref_designation( char *desig);
int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
char *make_config_dir_name( char *oname, const char *iname);  /* miscell.cpp */
int download_a_file( const char *ofilename, const char *url);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);
char **load_file_into_memory( const char *filename, size_t *n_lines,
                        const bool fail_if_not_found);      /* mpc_obs.cpp */
void shellsort_r( void *base, const size_t n_elements, const size_t esize,
         int (*compare)(const void *, const void *, void *), void *context);
void *bsearch_ext( const void *key, const void *base0,
      size_t nmemb, const size_t size,                   /* shellsor.cpp */
      int (*compar)(const void *, const void *), bool *found);
int string_compare_for_sort( const void *a, const void *b, void *context);
const char *get_find_orb_text( const int index);      /* elem_out.cpp */
static void reduce_designation( char *desig, const char *idesig);
int set_tholen_style_sigmas( OBSERVE *obs, const char *buff);  /* mpc_obs.c */
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
#ifdef _MSC_VER
     /* Microsoft Visual C/C++ has no strncasecmp.  strncmp will do.  */
#define strncasecmp strncmp
#endif
int compute_canned_object_state_vect( double *loc, const char *mpc_code,
                     const double jd);                 /* elem_out.cpp */
char *mpc_station_name( char *station_data);       /* mpc_obs.cpp */
int get_object_name( char *obuff, const char *packed_desig);   /* mpc_obs.c */
void compute_error_ellipse_adjusted_for_motion( double *sigma1, double *sigma2,
                  double *posn_angle, const OBSERVE *obs,
                  const MOTION_DETAILS *m);                  /* orb_func.cpp */
double n_nearby_obs( const OBSERVE FAR *obs, const unsigned n_obs,
          const unsigned idx, const double time_span);       /* orb_func.cpp */
void convert_ades_sigmas_to_error_ellipse( const double sig_ra,
         const double sig_dec, const double correl, double *major,
         double *minor, double *angle);                      /* errors.cpp */

#ifndef PATH_MAX
   #define PATH_MAX 256
#endif

static void *_ades_ids_stack = NULL;

int debug_printf( const char *format, ...)
{
   const char *debug_file_name = "debug.txt";
   FILE *ofile = fopen_ext( debug_file_name, "ca");

   if( ofile)
      {
      va_list argptr;
      const time_t t0 = time( NULL);
      const long max_debug_file_size = 10000000;  /* 10 MBytes should be enough */

      if( ftell( ofile) > max_debug_file_size)
         {
         fclose( ofile);
         ofile = fopen_ext( debug_file_name, "cw");
         }
      fprintf( ofile, "%02d:%02d:%02d ",
               ((int)t0 / 3600) % 24, ((int)t0 / 60) % 60, (int)t0 % 60);
      va_start( argptr, format);
      vfprintf( ofile, format, argptr);
      va_end( argptr);
      if( *format && format[strlen( format) - 1] != '\n')
         fprintf( ofile, "\n");     /* ensure a line break */
      fclose( ofile);
      }
   return( 0);
}

/* Quick and dirty check on MPC80 format compliance.  FAILS on
Find_Orb format extensions,  which is OK for its current
use in ensuring that lines would be accepted by MPC. */

static int quick_mpc80_check( const char *buff)
{
   const char *check = "nnnn nn nn.nnnnn?nn nn nn.nn??nn nn nn.n";
   size_t i;
   int rval = 0;

   buff += 15;
   for( i = 0; !rval && check[i]; i++)
      switch( check[i])
         {
         case 'n':
            if( buff[i] < '0' || buff[i] > '9')
               rval = -1;
            break;
         case ' ':
            if( buff[i] != ' ')
               rval = -1;
            break;
         default:
            break;
         }
   if( i >= 17 && strchr( "Rrvs", buff[-1]))
      rval = 0;   /* radar obs or 2nd line for rover or satellite obs;  */
                  /* the date will check out OK in such cases           */
   return( rval);
}

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

/* By default,  Find_Orb will ignore observations before 1100 AD or after
2300 AD.  This range can be reset at runtime (see the '-r' option in
findorb.cpp) or as a configuration parameter (see OBSERVATION_DATE_RANGE
in environ.dat).        */

double minimum_observation_jd = 0.;
double maximum_observation_jd = 0;

static bool is_in_range( const double jd)
{
   return( jd > minimum_observation_jd && jd < maximum_observation_jd);
}

/* Packed designations are routinely misaligned.  This fixes the most common
cases,  using the following rules :
   A seven- or eight-byte packed designation should always be right-aligned.
(Which will mean starting in column 6 or 5,  respectively.)
   If a five-byte desig starts in column 0,  it's probably a numbered object;
leave it alone.  Otherwise,  it's probably a temp desig and should start in
column 6.  (This gets a little tricky.  Five-byte desigs should always
start either in column 1,  if it's a numbered object,  or column 6... but
figuring out which way to go is difficult.)
   1-4 byte or six-byte desigs are definitely temporary and should always
start in column 6 and leave blank(s) at the end. */

static void check_packed_desig_alignment( char *buff)
{
   int i = 0, j = 11, new_i, len;

   while( buff[i] == ' ' && i < j)
      i++;
   while( buff[j] == ' ' && i < j)
      j--;
   len = j - i + 1;
   if( len > 8)            /* over-long */
      return;
   new_i = i;
   switch( len)
      {
      case 7:
      case 8:
         new_i = 12 - len;
         break;
      case 5:
         if( i)
            new_i = 5;
         break;
      default:
         new_i = 5;
      }
   if( i != new_i)
      {
//    debug_printf( "Was : '%s';  len %d, %d to %d\n", buff, len, i, new_i);
      memmove( buff + new_i, buff + i, len);
      memset( buff, ' ', new_i);
      memset( buff + new_i + len, ' ', 12 - new_i - len);
//    debug_printf( "Now : '%s'\n", buff);
      }
}

/* In some situations,  MPC observations end up with one or more
leading or trailing spaces.  Or people copy/paste observations and
leave off a space or three at the beginning,  or don't put the
designation exactly where MPC wants it to be.  This function attempts
to puzzle out 'what the observer really meant to say',  for the most
common sorts of errors I've seen.

   If the input line appears to be a malformed observation with fixable
errors, the return value will have bits set from the following values.
Otherwise,  zero will be returned. */

#define OBS_FORMAT_LEADING_SPACES               1
#define OBS_FORMAT_WRONG_DESIG_PLACEMENT        2
#define OBS_FORMAT_INCORRECT                    4

static int fix_up_mpc_observation( char *buff, double *jd)
{
   size_t len = strlen( buff);
   int rval = 0;
   char tchar, packed[13];

   while( len > 40 && buff[len - 1] <= ' ')
      len--;                  /* lop off trailing spaces */
   buff[len] = '\0';
   if( len <= 40 || !is_valid_mpc_code( buff + len - 3))
      return( 0);
   if( len != 80 && len > 70 && is_valid_mpc_code( buff + len - 3)
                  && !quick_mpc80_check( buff + len - 80))
      {
      if( len < 80)        /* insert missing leading spaces */
         {
         memmove( buff + 80 - len, buff, len + 1);
         memset( buff, ' ', 80 - len);
         }
      else                 /* remove spurious leading spaces */
         memmove( buff, buff + len - 80, 81);
      rval = OBS_FORMAT_LEADING_SPACES;
      len = 80;
      }
   tchar = buff[12];
   buff[12] = '\0';
   if( !create_mpc_packed_desig( packed, buff))
      memcpy( buff, packed, 12);
   else
      check_packed_desig_alignment( buff);
   buff[12] = tchar;
   if( len == 80)
      {
      double temp_jd = observation_jd( buff);

      if( jd)
         *jd = temp_jd;
      if( temp_jd)      /* doesn't need fixing */
         return( rval);
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
               else        /* formatting trouble;  giving up */
                  return( 0);
               }
            if( sscanf( buff + bytes_read, "%9s%9s%n", minutes, seconds,
                           &tval) != 2)
               return( 0);
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
                  return( 0);
               }
            if( sscanf( buff + bytes_read, "%9s%9s%n", minutes, seconds,
                           &tval) != 2)
               return( 0);
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
            while( len > 3 && buff[len - 1] <= ' ')
               len--;
            if( len == 81 && buff[80] == 'x')
               {               /* CSS folk sometimes add 'x' to a */
               len = 80;       /* line to mark it as deleted. */
               obuff[64] = 'x';
               }
            memcpy( obuff + 77, buff + len - 3, 3);
            obuff[80] = '\0';
            if( strcmp( buff, obuff))
               {
               rval |= OBS_FORMAT_INCORRECT;
               strlcpy_err( buff, obuff, 81);
               }
            }
         }
      }
   return( rval);
}

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

int generic_message_box( const char *message, const char *box_type)
{
   int rval, color = (*box_type == '!' ? COLOR_ATTENTION : COLOR_DEFAULT_INQUIRY);

   rval = inquire( message, NULL, 30, color);
   debug_printf( "%s", message);
   return( rval);
}

/* https://www.minorplanetcenter.net/iau/info/Astrometry.html#HowObsCode
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
are located.  Also (again Find_Orb only),  the degrees or minutes can be
expressed as decimals,  so the above example for Bowdoinham could be

COD Bow
COM Long. 69.9 0 0 W, Lat. 44.019444 0 0 N, Alt. 50m, approximate

   You _can_ use this to override the existing codes,  but the only case
I can think of where you'd do that would be if you thought there might
be an error in the MPC's ObsCodes.html listing.  (Or in the supplement
in 'rovers.txt'.)    */

static void *obs_details;

static inline int get_lat_lon_from_header( double *lat,
            double *lon, double *alt, const char *mpc_code,
            const char **name_from_header)
{
   const char **lines = get_code_details( obs_details, mpc_code);
   size_t i;
   int rval = -1;

   *name_from_header = NULL;
   for( i = 0; rval && lines && lines[i]; i++)
      {
      static bool warning_shown = false;
      mpc_code_t cinfo;

      rval = get_xxx_location_info( &cinfo, lines[i]);
      if( rval == -2 && !warning_shown)      /* malformed position line */
         {
         char tbuff[200];

         warning_shown = true;      /* just do this once */
         snprintf_err( tbuff, sizeof( tbuff),   /* see efindorb.txt */
                  get_find_orb_text( 2000), mpc_code);
         generic_message_box( tbuff, "!");
         }
      else if( !rval)
         {
         *lat = cinfo.lat * 180. / PI;
         *lon = cinfo.lon * 180. / PI;
         *alt = cinfo.alt;
         }
      }
   if( !rval && strlen( lines[0]) > 8)
      *name_from_header = lines[0] + 8;
               /* if observatory name is specified in header,  e.g., */
               /* COD Bow Generic Observatory,  Bowdoinham */
   return( rval);
}

extern int debug_level;

/* A return value of -1 indicates malformed input; -2 indicates a satellite */
/* (for which no position is available) .  Anything else indicates          */
/* the number of the planet on which the station is located (earth=3,  the  */
/* usual value... but note that rovers.txt contains extraterrestrial "MPC   */
/* stations",  so values from 0 to 10 may occur.)                           */

static int extract_mpc_station_data( const char *buff, mpc_code_t *cinfo)
{
   int rval;

   if( !cinfo)
      {
      mpc_code_t unused;

      rval = get_mpc_code_info( &unused, buff);
      }
   else                          /* keep longitude in -180 to +180 */
      {                          /* put parallaxes in AU */
      rval = get_mpc_code_info( cinfo, buff);
      if( rval >= 0)
         {
         const double scale = planet_radius_in_meters( rval) / AU_IN_METERS;

         cinfo->rho_cos_phi *= scale;
         cinfo->rho_sin_phi *= scale;
         if( cinfo->lon > PI)
             cinfo->lon -= PI + PI;
         }
      else
         {
         memset( cinfo, 0, sizeof( mpc_code_t));
         cinfo->planet = rval;
         }
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
            char buff[100];

            while( fgets_trimmed( buff, sizeof( buff), ifile))
               {
               const int planet_idx = extract_mpc_station_data( buff, NULL);

               if( planet_idx != -1)
                  {
                  if( rval)
                     {
                     rval[*n_stations] = codes + buff_loc;
                     strlcpy_err( rval[*n_stations], buff, sizeof( buff));
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
   int line_no = 0, rval = 0;

   assert( ifile);
   while( !rval && fgets_trimmed( tbuff, sizeof( tbuff), ifile))
      if( *tbuff != ';')
         {
         if( atoi( tbuff) == atoi( mpc_code + 3))
            {
            if( buff)
               {
               strlcpy_err( buff + 30, tbuff + 19, 15);
               memcpy( buff, mpc_code, 3);
               }
            rval = line_no + 100;
            }
         else
            line_no++;
         }
   fclose( ifile);
   return( rval);
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

static int mpc_code_cmp( const void *ptr1, const void *ptr2)
{
   const char **p1 = (const char **)ptr1;
   const char **p2 = (const char **)ptr2;

   return( memcmp( *p1, *p2, 4));
}

/* The first (247) roving observer retains that code.  If another
rover is found with a different lat/lon,  it is assigned (24a).
The 27th rover is assigned (24z).  Rovers 28 to 53 get codes
(24A) to (24Z).  After that,  the second character runs from
A to Z (for 247 codes) or a to z (270) codes,  and the third
is a mutant hex (0...9, A-Z, a-z),  for an additional
26*62 = 1612 codes.  1665 codes should be enough for anybody. */

static int get_rover_index( const char *obscode)
{
   int rval = -1;

   if( obscode[0] == '2' && (obscode[1] == '4' || obscode[1] == '7'))
      {
      if( obscode[1] == '4' && obscode[2] == '7')     /* 247 */
         rval = 0;
      else if( obscode[1] == '7' && obscode[2] == '0')     /* 270 */
         rval = 0;
      else if( obscode[2] >= 'a')
         rval = obscode[2] - 'a' + 1;
      else if( obscode[2] >= 'A')
         rval = obscode[2] - 'A' + 27;
      }
   if( obscode[0] == '2' && isalpha( obscode[1]))
      {
      if( obscode[1] > 'Z')
         rval = (obscode[1] - 'a') * 62;
      else
         rval = (obscode[1] - 'A') * 62;
      rval += 53 + mutant_hex_char_to_int( obscode[2]);
      }
   return( rval);
}

/* The following function paws through the ObsCodes.htm or ObsCodes.html
   file,  looking for the observer code in question.  If found,  the
   line is simply copied into 'buff'.  If lon_in_radians and the
   rho_xxx_phi values are non-NULL,  they're extracted from the buffer.

   There are a few "supplemental" observers,  mostly satellite observers
   who don't have MPC codes.  These could be handled as roving observers
   (code 247),  or as 'temporary' observers (code XXX).  However,  it can
   be better if they have their own codes.  These non-MPC-approved codes
   are put into 'rovers.txt'.  They are three characters,  but not of the
   uppercase letter and two digits sort;  that way,  they don't look
   too much like "real,  official" MPC designations.  Also,  some fixes
   for defective MPC code lat/lon/alt values are provided in the file.

   Return value:

      -2:  obscodes.htm,  obscodes.html not found (need one of these)
       Other:  index of planet of MPC station (3=earth,  most common
            case;  0=sun, 1=mercury,  etc.)
*/

typedef struct
{
   double lon, lat, alt;         /* alt is in meters */
} rover_t;

static rover_t *rovers = NULL;
int n_obs_actually_loaded, n_rovers = 0;

int get_observer_data( const char FAR *mpc_code, char *buff, mpc_code_t *cinfo)
{
   static char *curr_station = NULL;
   static char **station_data = NULL;
   static int n_stations = 0;
   const char *blank_line = "!!!   0.0000 0.000000 0.000000Unknown Station Code";
   int rval = -1, rover_idx;
   size_t i;
   const char *override_observatory_name = NULL;
   double lat0 = 0., lon0 = 0., alt0 = 0.;
   char temp_code[5];
   const size_t buffsize = 81;

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
   if( !memcmp( mpc_code, "@90", 3))      /* looking for dynamical point info */
      {
      for( i = 1; i < (size_t)n_stations; i++)
         if( !memcmp( station_data[i] + 41, mpc_code, 5))
            {
            strcpy( buff, station_data[i]);
            return( 0);
            }
      assert( 0);          /* should never get here */
      return( -1);
      }

   memcpy( temp_code, mpc_code, 4);
   if( 4 != strlen( mpc_code))
      temp_code[3] = ' ';
   temp_code[4] = '\0';

   if( !cinfo)   /* attempting to look up an MPC code from the station name */
      {
      int pass;

      for( pass = 0; pass < 2; pass++)
         for( i = 0; station_data[i]; i++)
            if( (!pass && !memcmp( buff, station_data[i] + 30, strlen( buff)))
                     || (pass && strstr( station_data[i] + 30, buff)))
               {
               strlcpy( buff, station_data[i], buffsize);
               return( 0);
               }
      return( -1);
      }
   memset( cinfo, 0, sizeof( mpc_code_t));
   cinfo->planet = 3;         /* default to geocentric */

   if( !get_lat_lon_from_header( &lat0, &lon0, &alt0, temp_code,
                                             &override_observatory_name))
      if( !override_observatory_name)
         override_observatory_name = "Temporary MPC code";

   rover_idx = get_rover_index( temp_code);
   if( rover_idx >= 0)
      {
      if( rover_idx < n_rovers)
         {
         lat0 = rovers[rover_idx].lat;
         lon0 = rovers[rover_idx].lon;
         alt0 = rovers[rover_idx].alt;
         }
      override_observatory_name = "Roving observer";
      }

   if( !get_lat_lon_info( cinfo, mpc_code))
      {
      mpc_code = "zzzz";
      lat0 = cinfo->lat * 180. / PI;
      lon0 = cinfo->lon * 180. / PI;
      alt0 = cinfo->alt;
      override_observatory_name = "User-supplied observer";
      }

   if( override_observatory_name)
      {
      char tbuff[90];

      strlcpy_error( tbuff, mpc_code);
      strlcat_error( tbuff, "   ");
      snprintf_err( tbuff + 4, sizeof( tbuff) - 4, "!%15.9f%13.9f%13.3f %s",
                    lon0, lat0, alt0, override_observatory_name);
      if( buff)
         strlcpy_err( buff, tbuff, sizeof( tbuff));
      rval = extract_mpc_station_data( tbuff, cinfo);
      return( rval);
      }

   if( !memcmp( mpc_code, "Ast", 3))
      {
      memset( cinfo, 0, sizeof( mpc_code_t));
      if( buff)
         strlcpy_err( buff, blank_line, buffsize);
      rval = get_asteroid_observer_data( mpc_code, buff);
      cinfo->planet = rval;
      return( rval);
      }

   mpc_code = temp_code;
   if( !curr_station || mpc_code_cmp( &curr_station, &mpc_code))
      {
      char **search = (char **)bsearch_ext( &mpc_code, station_data, n_stations,
                  sizeof( char *), mpc_code_cmp, NULL);

      curr_station = (search ? *search : NULL);
      }
   if( !curr_station)
      {
      const char *envar = "UPDATE_OBSCODES_HTML";
      const double curr_t = current_jd( );
      double last_t = 0., delay = 0.1;

      sscanf( get_environment_ptr( envar), "%lf %lf", &last_t, &delay);
      debug_printf( "Couldn't find MPC station '%s'\n", mpc_code);
      debug_printf( "Curr t %f; last t %f; delay %f\n",
                  curr_t, last_t, delay);
      if( curr_t > last_t + delay && isalnum( mpc_code[0])
                                  && isdigit( mpc_code[1]) && isdigit( mpc_code[2]))
         {
         char path[PATH_MAX];

         snprintf( path, sizeof( path), "%f %f", curr_t, delay);
         set_environment_ptr( envar, path);
         make_config_dir_name( path, "ObsCodes.htm");
         if( !download_a_file( path,
                  "http://www.minorplanetcenter.org/iau/lists/ObsCodes.html"))
            {
            get_observer_data( NULL, NULL, NULL);
            return( get_observer_data( mpc_code, buff, cinfo));
            }
         }
      if( buff)
         {
         strcpy( buff, blank_line);
         memcpy( buff, mpc_code, 3);
         }
      }
   else
      {
      if( buff)
         strcpy( buff, curr_station);
      rval = extract_mpc_station_data( curr_station, cinfo);
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
   mpc_code_t cinfo;
   int rval;

   rval = get_observer_data( mpc_code, buff, &cinfo);
   *lon_in_radians = cinfo.lon;
   *lat_in_radians = cinfo.lat;
   if( alt_in_meters)
      *alt_in_meters = cinfo.alt;
   return( rval);
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

   if( !memcmp( obs->mpc_code, "500", 3))
      rval = -1;
   if( !rval)
      {
      DPT ra_dec;
      const double ut1 = obs->jd - td_minus_ut( obs->jd) / seconds_per_day;

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
            full_ra_dec_to_alt_az( &ra_dec, alt_az, NULL, &latlon, ut1, NULL);
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
   if( !rval)
      {
      sun_alt_az->x *= 180. / PI;
      sun_alt_az->x += 180.;
      sun_alt_az->y *= 180. / PI;
      object_alt_az->x *= 180. / PI;
      object_alt_az->x += 180.;
      object_alt_az->y *= 180. / PI;
      }
   return( rval);
}

char **load_file_into_memory( const char *filename, size_t *n_lines,
                        const bool fail_if_not_found)
{
   FILE *ifile = fopen_ext( filename, (fail_if_not_found ? "fcrb" : "crb"));

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
               strcat( name, strchr( name, '/') ? " " : "/");
               strcat( name, buff + slen);
               }
            }
      fclose( ifile);
      }
}

/* To get data for a particular numbered object from the MPC's file at

https://www.minorplanetcenter.net/iau/lists/NumberedMPs.txt.gz

we do a search based on the fact that the average line is 88.6 bytes
long (including line feed at the end).  So we search to an estimated
location,  almost certainly in the middle of a line;  fgets() the
rest of the line;  then fgets() a full line.  The minor planet number
read from that line will either be slightly less than the target
number (in which case we just read a few lines),  or will let us
compute a more accurate place in the file to seek to.

At present,  this program uses the object name or provisional designation.
At some point,  we might show discovery date, location,  reference,
and/or discoverer name.       */

char *find_numbered_mp_info( const int number)
{
   long loc = (number - 1) * 886 / 10;
   int i, curr_no;
   bool got_it = false;
   const int tolerance = 10;        /* try to get within ten lines */
   static int cached_number = -1;
   static char *cached = NULL;
   char buff[200];
   FILE *ifile;

   if( !number)
      {
      if( cached)
         free( cached);
      cached = NULL;
      cached_number = -1;
      return( 0);
      }
   if( cached_number == -2)      /* haven't got the file */
      return( NULL);
   ifile = fopen_ext( "NumberedMPs.txt", "crb");
   if( !ifile)
      {
      cached_number = -2;
      return( NULL);
      }
   while( !got_it)
      {
      if( number < tolerance)
         curr_no = 0;
      else
         {
         fseek( ifile, loc, SEEK_SET);
         for( i = 0; i < 2; i++)
            if( !fgets( buff, sizeof( buff), ifile))
               {
               fclose( ifile);
               return( NULL);
               }
         loc += 190;
         for( i = 0; i < 7 && buff[i] != '('; i++)
            ;
         curr_no = atoi( buff + i + 1);
         }
      if( curr_no > number - tolerance && curr_no <= number)
         {
         for( i = number - curr_no; i; i--)
            if( !fgets( buff, 200, ifile))
               {
               fclose( ifile);
               return( NULL);
               }
         got_it = true;    /* success */
         }
      else
         loc += (number - curr_no - 4) * 89;
      }
   fclose( ifile);
   if( !cached)
      cached = (char *)malloc( 200);
   strlcpy_err( cached, buff, 200);
   cached_number = number;
   return( cached);
}

/* An object with a name such as 1989-013A is probably an artsat,  and
probably has a NORAD designation,  name,  and other info in 'satcat.html',
a master list of artsats available at

https://planet4589.org/space/gcat/data/cat/satcat.html

   The following code can turn,  for example,  '1966-092A' into
'1966-092A = NORAD 02501 = Molniya-1'.       */

static bool try_artsat_xdesig( char *name)
{
   FILE *ifile = fopen_ext( "satcat.html", "crb");
   bool found_a_match = false;

   if( ifile)
      {
      char tbuff[500];
      size_t slen;
      size_t max_out = 80;    /* max len of 'name' will be 80 bytes */

      while( *name == ' ')    /* skip leading spaces */
         {
         name++;
         max_out--;
         }
      remove_trailing_cr_lf( name);
      slen = strlen( name);
      if( max_out > slen)
         while( !found_a_match && fgets( tbuff, sizeof( tbuff), ifile))
            if( !memcmp( tbuff + 20, name, slen) && tbuff[slen + 20] == ' ')
               {
               found_a_match = true;
               snprintf_append( name, max_out - slen, " = NORAD %.5s = %.28s",
                        tbuff + 1, tbuff + 48);
               remove_trailing_cr_lf( name);
               }
      fclose( ifile);
      }
   return( found_a_match);
}

/* Looks for designations of the form YYYY-NNN(letter).  */

static int is_artsat_desig( const char *desig)
{
   size_t slen;

   for( slen = strlen( desig); slen > 8; desig++, slen--)
      if( desig[4] == '-' && atoi( desig) > 1956 && atoi( desig) < 2100
               && isdigit( desig[5]) && isdigit( desig[6])
               && isdigit( desig[7]) && isupper( desig[8]))
         return( 1);
   return( 0);
}

/* In an MPC astrometric report line,  the name can be stored in assorted
   highly scrambled and non-intuitive ways.  Those ways _are_ documented
   on the MPC Web site,  which is probably the best place to look to
   understand why the following bizarre code does what it does:

   https://www.minorplanetcenter.org/iau/info/PackedDes.html

   Return values are as follows.  Note that 'other' includes artsats and
   temporary designations.        */

typedef struct
{
   void *next;
   char *line;
} odd_name_t;

int get_object_name( char *obuff, const char *packed_desig)
{
   int rval;
   size_t i, gap;
   static size_t n_lines;
   static char **extra_names = NULL;
   static odd_name_t *added = NULL;
   char xdesig[40];

   if( !packed_desig && !obuff)   /* flag to free up internal memory */
      {
      if( extra_names)
         free( extra_names);
      extra_names = NULL;
      while( added)
         {
         odd_name_t *next = (odd_name_t *)added->next;
         free( added->line);
         free( added);
         added = next;
         }
      return( 0);
      }

   if( !packed_desig && obuff)      /* adding a new name using 'COM desig' */
      {
      odd_name_t *next = (odd_name_t *)calloc( 1, sizeof( odd_name_t));

      next->line = (char *)malloc( strlen( obuff) + 1);
      strcpy( next->line, obuff);
      xref_designation( next->line);
      next->next = (odd_name_t *)added;
      added = next;
      return( 0);
      }

   if( !extra_names)          /* see 'odd_name.txt' for comments on this */
      {
      int sort_column = 0;

      extra_names = load_file_into_memory( "odd_name.txt", &n_lines, true);
      shellsort_r( extra_names, n_lines, sizeof( char *),
                     string_compare_for_sort, &sort_column);
      }

   reduce_designation( xdesig, packed_desig);
   xref_designation( xdesig);

   if( obuff)
      for( i = 0, gap = 0x8000; gap; gap >>= 1)
         if( i + gap < n_lines)
            {
            const int compare = memcmp( xdesig, extra_names[i + gap], 12);

            if( compare >= 0)
               i += gap;
            if( !compare)
               {
               strcpy( obuff, extra_names[i] + 13);
               if( is_artsat_desig( obuff))
                  return( OBJ_DESIG_ARTSAT);
               return( get_object_name( NULL, packed_desig));
               }
            }

   if( obuff && added)
      {
      odd_name_t *tptr = added;

      while( tptr)
         {
         if( !memcmp( xdesig, tptr->line, 12))
            {
            strcpy( obuff, tptr->line + 13);
            if( is_artsat_desig( obuff))
               return( OBJ_DESIG_ARTSAT);
            return( get_object_name( NULL, packed_desig));
            }
         tptr = (odd_name_t *)tptr->next;
         }
      }

   rval = unpack_mpc_desig( obuff, xdesig);

   if( obuff)
      {
      if( rval == OBJ_DESIG_ARTSAT)
         try_artsat_xdesig( obuff);
      if( rval == OBJ_DESIG_COMET_NUMBERED)
         {
         snprintf_err( xdesig, sizeof( xdesig), "%3d%c/", atoi( obuff), xdesig[4]);
         try_adding_comet_name( xdesig, obuff);
         }
      if( rval == OBJ_DESIG_COMET_PROVISIONAL)
         {
         snprintf_err( xdesig, sizeof( xdesig), "   %s ", obuff);
         try_adding_comet_name( xdesig, obuff);
         }
      if( rval == OBJ_DESIG_ASTEROID_NUMBERED)
         {
         const int number = atoi( obuff + 1);
         char *info = find_numbered_mp_info( number);

         snprintf_err( obuff, 11, "(%d)", number);
         obuff += strlen( obuff);
         if( info)
            {
            if( info[29] != ' ')        /* prov ID */
               {
               strlcpy( obuff, " =", 3);
               strlcpy( obuff + 2, info + 28, 12);
               remove_trailing_cr_lf( obuff);
               }
            if( info[9] != ' ')              /* append name */
               strlcpy( obuff + strlen( obuff), info + 8, 19);
            remove_trailing_cr_lf( obuff);
            }
         }
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

int load_environment_file( const char *filename);          /* mpc_obs.cpp */
static int load_default_environment_file( void);           /* mpc_obs.cpp */
void update_environ_dot_dat( void);                        /* mpc_obs.cpp */

static int earth_lunar_posn_vel( const double jd, double FAR *earth_loc,
                             double FAR *lunar_loc, const int is_vel)
{
   const int vel_offset = (is_vel ? PLANET_POSN_VELOCITY_OFFSET : 0);

   if( earth_loc)
      planet_posn( PLANET_POSN_EARTH + vel_offset, jd, earth_loc);
   if( lunar_loc)
      planet_posn( PLANET_POSN_MOON  + vel_offset, jd, lunar_loc);
   return( 0);
}

int earth_lunar_posn( const double jd, double FAR *earth_loc, double FAR *lunar_loc)
{
   return( earth_lunar_posn_vel( jd, earth_loc, lunar_loc, 0));
}

int earth_lunar_vel( const double jd, double FAR *earth_loc, double FAR *lunar_loc)
{
   return( earth_lunar_posn_vel( jd, earth_loc, lunar_loc, 1));
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
      const double omega = planet_rotation_rate( planet_no, 0) * PI / 180.;
                   /* planet rotation rate,  in radians/day */

      if( !omega)
         generic_message_box( "'cospar.txt' not found!", "!");
      assert( omega);
      if( offset)
         offset[i] = (rho_cos_phi * precess_matrix[i]
                    + rho_sin_phi * precess_matrix[i + 6]);
      if( vel)
         vel[i] = -rho_cos_phi * precess_matrix[i + 3] * omega;
      }
   return( 0);
}

static int get_canned_object_posn_vel( double *posn, double *vel, const double jde)
{
   extern const char *mpc_code_for_ephems;
   double state_vect[6];
   const int rval = compute_canned_object_state_vect( state_vect,
                                    mpc_code_for_ephems, jde);

   if( posn)
      memcpy( posn, state_vect, 3 * sizeof( double));
   if( vel)
      memcpy( vel, state_vect + 3, 3 * sizeof( double));
   return( rval);
}

/* To compute the location for L4 or L5,  we take the state vector and
compute the cross-product of position and velocity to get a vector
normal to the plane.  We then rotate both position and velocity
around that by +/- 60 degrees,  depending on whether we're at L4
or L5. */

static void _rotate_state_vector( double *state_vect, const double angle)
{
   double normal[3], orthog[6];
   size_t i;
   const double cos_ang = cos( angle), sin_ang = sin( angle);

   vector_cross_product( normal, state_vect, state_vect + 3);
   normalize_vect3( normal);
   vector_cross_product( orthog, state_vect, normal);
   vector_cross_product( orthog + 3, state_vect + 3, normal);
   for( i = 0; i < 6; i++)
      state_vect[i] = cos_ang * state_vect[i] + sin_ang * orthog[i];
}

/* See remarks in 'rovers.txt' about how L1-5 are computed. */

static int _compute_lagrange_point( double *vect, const int point_number,
                      const double jde, const bool get_vel)
{
   static double state[6], cached_jde;
   static double multiplier;
   static int cached_point_number = -1, obj1 = 0, obj2 = 0;
   int i;

   if( cached_point_number != point_number)
      {
      char buff[100], mpc_code[7];

      snprintf_err( mpc_code, sizeof( mpc_code), "@%d", point_number);
      get_observer_data( mpc_code, buff, NULL);
      i = sscanf( buff + 5, "%d %d %lf", &obj1, &obj2, &multiplier);
      assert( i == 3);
      assert( obj1 >= 0 && obj1 <= 10);
      assert( obj2 >= 0 && obj2 <= 10);
      if( obj1 == 3 && obj2 == 10)     /* Earth-Moon case */
         {
         obj1 = PLANET_POSN_EARTH;     /* see 'pl_cache.h' */
         obj2 = PLANET_POSN_MOON;
         }
      cached_point_number = point_number;
      }
   if( cached_jde != jde)
      {
      const int lpoint = (point_number - 1) % 5 + 1;
      double state1[6], state2[6];

      planet_posn( obj1, jde, state1);
      planet_posn( obj1 + PLANET_POSN_VELOCITY_OFFSET, jde, state1 + 3);
      planet_posn( obj2, jde, state2);
      planet_posn( obj2 + PLANET_POSN_VELOCITY_OFFSET, jde, state2 + 3);
      for( i = 0; i < 6; i++)    /* make state2 = vector from primary to secondary */
         state2[i] -= state1[i];
      if( lpoint >= 1 && lpoint <= 3)     /* collinear cases */
         {
         for( i = 0; i < 6; i++)
            state[i] = state1[i] + multiplier * state2[i];
         }
      else                        /* L4, L5 : equilateral triangle cases */
         {
         _rotate_state_vector( state2, (lpoint == 5 ? PI / 3. : -PI / 3.));
         for( i = 0; i < 6; i++)
            state[i] = state1[i] + state2[i];
         }
      cached_jde = jde;
      }
   for( i = 0; i < 3; i++)
      vect[i] = state[i + (get_vel ? 3 : 0)];
   return( 0);
}

int compute_observer_loc( const double jde, const int planet_no,
               const double rho_cos_phi,
               const double rho_sin_phi, const double lon, double FAR *offset)
{

   if( planet_no == -2)
      return( get_canned_object_posn_vel( offset, NULL, jde));

   if( planet_no > 9000)
      _compute_lagrange_point( offset, planet_no, jde, false);
   else if( planet_no != 3 && planet_no != 10)
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

      assert( planet_no < 9000);
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
   int i;

   if( planet_no == -2)    /* spacecraft-based observation */
      return( get_canned_object_posn_vel( NULL, vel, jde));
   if( planet_no > 9000)
      _compute_lagrange_point( vel, planet_no, jde, true);
   else if( planet_no != 3 && planet_no != 10)
      planet_posn( (planet_no >= 0 ? planet_no : 12) + PLANET_POSN_VELOCITY_OFFSET,
                   jde, vel);        /* planet_no == -1 means SS Barycenter */
   else
      earth_lunar_vel( jde, (planet_no == 3) ? vel : NULL,
                            (planet_no == 3) ? NULL : vel);
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
static double override_time = 0., override_ra = -100., override_dec = -100.;
static double override_mag = -100.;
      /* 'override_time' allows one to set the observation time with the */
      /* #time keyword (see below) or with certain keywords.  Similarly  */
      /* for override_ra and override_dec,  primarily used for ADES  */
      /* observations where the 80-column field lacks precision.     */
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

/* Converts JPL's SPICE codes,  such as 399=geocenter,  599=Jupiter body
center,  etc.  to plain ol' Sun=0, Merc=1,  etc. codes used in JPL's
DE ephemerides.  This matters because in ADES,  satellite reference
locations are specified using SPICE index values,  and we need to know
which object to look up from the DE ephem.  (Almost all such offsets
are relative to the geocenter.  But Gaia,  for example,  is relative
to the solar system barycenter.)   */

static int jpl_code_to_planet_idx( const int spice_code)
{
   int rval = -1;

   if( !spice_code)        /* SSB */
      rval = 12;
   else if( spice_code == 10)    /* sun */
      rval = 0;
   else if( spice_code % 100 == 99)
      rval = spice_code / 100;
   assert( rval >= 0);
   return( rval);
}

void set_up_observation( OBSERVE FAR *obs)
{
   char tbuff[300];
   mpc_code_t cinfo;
   int observer_planet = get_observer_data( obs->mpc_code, tbuff, &cinfo);

   if( observer_planet == -1)
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
         int text_to_add;

         FMEMCPY( unfound[n_unfound++], obs->mpc_code, 3);
         snprintf_err( tbuff, sizeof( tbuff), get_find_orb_text( 2002),
                             obs->mpc_code);
         if( strcmp( obs->mpc_code, "XXX"))
            text_to_add = 2003;   /* See efindorb.txt.  These reference */
         else                     /* possible error messages.  */
            text_to_add = 2004;
         strlcat_err( tbuff, get_find_orb_text( text_to_add), sizeof( tbuff));
         generic_message_box( tbuff, "!");
         comment_observation( obs, "? NoCode");
         obs->is_included = 0;
         obs->flags |= OBS_DONT_USE;
         }
      observer_planet = 3;    /* default to geocentric */
      }
   else if( obs->note2 != 'S' && obs->note2 != 'V'
               && !memcmp( tbuff + 3, "                          ", 26))
      {
      obs->is_included = 0;
      obs->flags |= OBS_NO_OFFSET | OBS_DONT_USE;
      comment_observation( obs, "? offset");
      }
   if( observer_planet == -2)          /* satellite observation */
      observer_planet = jpl_code_to_planet_idx( obs->ref_center);
   compute_observer_loc( obs->jd, observer_planet,
               cinfo.rho_cos_phi, cinfo.rho_sin_phi, cinfo.lon, obs->obs_posn);
   compute_observer_vel( obs->jd, observer_planet,
               cinfo.rho_cos_phi, cinfo.rho_sin_phi, cinfo.lon, obs->obs_vel);
   set_obs_vect( obs);
}

static int set_data_from_obs_header( OBSERVE *obs);

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
int object_type;

int find_fcct_biases( const double ra, const double dec, const char catalog,
                 const double jd, double *bias_ra, double *bias_dec);

static int parse_observation( OBSERVE FAR *obs, const char *buff)
{
   unsigned time_format;
   double utc = extract_date_from_mpc_report( buff, &time_format);
   const bool is_radar_obs = (buff[14] == 'R' || buff[14] == 'r');
   static bool fcct_error_message_shown = false;
   const OBSERVE saved_obs = *obs;
   double coord_epoch = input_coordinate_epoch;
   int obj_desig_type;

   if( !utc)
      return( -1);
   assert( obs);
   assert( buff);
   memset( obs, 0, sizeof( OBSERVE));
   memcpy( obs->packed_id, buff, 12);
   obs->ra_bias = saved_obs.ra_bias;
   obs->dec_bias = saved_obs.dec_bias;
   obs->ref_center = saved_obs.ref_center;
   obs->packed_id[12] = '\0';
   obj_desig_type = get_object_name( NULL, obs->packed_id);
   if( obj_desig_type == OBJ_DESIG_OTHER)
      {
      size_t i = 7;

      while( i < 12 && isalnum( buff[i]))
         i++;
      while( i < 12 && buff[i] == ' ')
         i++;
      if( i != 12)
         {
         static bool warning_not_yet_shown = true;

         if( warning_not_yet_shown)
            {
            char tbuff[300];

            warning_not_yet_shown = false;
            snprintf_err( tbuff, sizeof( tbuff), get_find_orb_text( 2074),
                        obs->packed_id);
            generic_message_box( tbuff, "!");
            }
         }
      }
   if( obj_desig_type == OBJ_DESIG_COMET_PROVISIONAL
         || obj_desig_type == OBJ_DESIG_COMET_NUMBERED)
      object_type = OBJECT_TYPE_COMET;
   if( override_time)
      {
      utc = override_time;
      override_time = 0.;
      }
   utc += observation_time_offset;
   obs->jd = utc + td_minus_utc( utc) / seconds_per_day;

   if( !is_radar_obs)
      {
      obs->mag_band = buff[70];
      obs->astrometric_net_code = buff[71];
      }
   else
      obs->mag_band = obs->astrometric_net_code = ' ';
   obs->discovery_asterisk = buff[12];
   obs->note1 = buff[13];
   obs->note2 = buff[14];
   if( is_radar_obs || buff[14] == 's' || buff[14] == 'v')
      obs->ra = obs->dec = 0.;  /* radar data and "second lines" have no RA/dec */
   else
      {
      const int rval = get_ra_dec_from_mpc_report( buff,
                  &obs->ra_precision, &obs->ra, &obs->posn_sigma_1,
                  &obs->dec_precision, &obs->dec, &obs->posn_sigma_2);

      if( override_ra >= 0.)
         {
         obs->ra = override_ra * PI / 180.;
         override_ra = -100.;
         }
      if( override_dec >= -95.)
         {
         obs->dec = override_dec * PI / 180.;
         override_dec = -100.;
         }
      if( rval)
         return( rval);
      }

   obs->time_precision = time_format;
   if( is_radar_obs)                  /* Radar obs are treated as */
      obs->time_sigma = 0.;           /* being 'perfectly' timed  */
   else
      {
      obs->time_sigma = pow( .1, (double)( time_format % 10));
      if( time_format / 10 == 2)    /* CYYMMDD HH:MM:SS.ss.. formats */
         obs->time_sigma /= seconds_per_day;
      if( time_format / 10 == 5)    /* CYYMMDD HH:MM.mmm... formats */
         obs->time_sigma /= minutes_per_day;
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
      double year = JD_TO_YEAR( obs->jd);

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
      {
      obs->obs_mag = atof( buff + 65);
      if( override_mag >= -95.)
         {
         obs->obs_mag = override_mag;
         override_mag = -100.;
         }
      if( obs->mag_band == ' ' && strstr( "249 C49 C50", obs->mpc_code))
         obs->mag_band = 'V';   /* Sungrazer obs usually have a blank mag */
      }                         /* band, but are basically V */
   else
      obs->obs_mag = BLANK_MAG;
   FMEMCPY( obs->reference, buff + 72, 5);
   obs->reference[5] = '\0';
   set_data_from_obs_header( obs);
   if( memcmp( buff + 72, ".rwo ", 5) &&
          find_fcct_biases( obs->ra, obs->dec, obs->astrometric_net_code, obs->jd,
                                &obs->ra_bias, &obs->dec_bias) == -2)
      {        /* i.e.,  we tried to get FCCT14 debiasing and failed */
      if( !fcct_error_message_shown && apply_debiasing)
         {
         generic_message_box( get_find_orb_text( 2005), "!");
         fcct_error_message_shown = true;   /* see efindorb.txt */
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
   size_t i = 12;

   memcpy( desig, idesig, 12);
   desig[12] = '\0';
   while( i > 0 && desig[i - 1] == ' ')
      i--;
   if( desig[5] != ' ' && i > 9)
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

const char *combine_all_observations;

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
   INTENTIONALLY_UNUSED_PARAMETER( ignored_context);
   return( strcmp( (const char *)a, (const char *)b));
}

static const char *new_xdesig_indicator = "New xdesig";

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
      assert( n_lines);
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

   if( strlen( desig) > 34 && !strcmp( desig + 26, new_xdesig_indicator))
      {                /* Add new xdesignation.  Increase the table size by */
      char *tptr;      /* one;  binary-search to find where new entry goes; */
      int i1;          /* move everything following it;  insert new entry   */

      xlate_table = (char *)realloc( xlate_table, (n_lines + 1) * 26);
      assert( xlate_table);
      for( i = -1, gap = 0x8000; gap; gap >>= 1)
         if( (i1 = i + gap) < n_lines)
            if( memcmp( xlate_table + i1 * 26, desig, 12) <= 0)
               i = i1;
      i++;
      tptr = xlate_table + i * 26;
      memmove( tptr + 26, tptr, (n_lines - i) * 26);
      memcpy( tptr, desig, 25);
      tptr[25] = '\0';
      n_lines++;
      return( 0);
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

   for( i = -1, gap = 0x8000; gap; gap >>= 1)
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

<p><b>VE82F84</b>   (may have other bytes before it)
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
   char *tptr;

   for( len = 0; buff[len] >= ' '; len++)
      ;
   if( debug_level > 9)
      debug_printf( "get_neocp_data len: %d\n", (int)len);
   while( len && buff[len - 1] == ' ')    /* eliminate trailing spaces */
      len--;
   if( len < 10 && len > 2 && *buff == '*' && buff[len - 1] == '*')
      {
      memcpy( desig, buff + 1, len - 2);
      desig[len - 2] = '\0';
      neocp_file_type = NEOCP_FILE_TYPE_ASTERISKED;
      }
   else if( len > 15 && (tptr = strstr( buff, "<p><b>")) != NULL)
      {
      char *endptr = strstr( tptr, "</b>");

      tptr += 6;
      if( endptr > tptr && endptr <= tptr + 7)
         {
         memcpy( desig, tptr, endptr - tptr);
         desig[endptr - tptr] = '\0';
         neocp_file_type = NEOCP_FILE_TYPE_HTML;
         }
      }
   else if( len && len < 10 && neocp_file_type == NEOCP_FILE_TYPE_UNKNOWN)
      {
      memcpy( desig, buff, len);
      desig[len] = '\0';
      }
   else if( strstr( buff, "the geocenter"))
      {
      static int already_warned = 0;

      memcpy( mpc_code, "500", 3);
      if( !already_warned)
         {                             /* warn about geocentric ephems */
         already_warned = 1;           /* see 'efindorb.txt' for details */
         generic_message_box( get_find_orb_text( 2006), "!");
         }
      }
   else if( (tptr = strstr( buff, " observatory code ")) != NULL)
      memcpy( mpc_code, tptr + 18, 3);
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
         if( mask == mask_if_hours)     /* ephem step size is in hours */
            minutes = atoi( buff + 11) * 60;
         else                       /* step size was in minutes */
            minutes = (atoi( buff + 11) / 100) * 60 +
                                          atoi( buff + 13);
         snprintf_err( obuff + 25, sizeof( obuff) - 25,
                          ".%06d", minutes * 1000000 / (int)minutes_per_day);
         memcpy( obuff + 32, buff + 18, 10);      /* RA */
         memcpy( obuff + 44, buff + 29, 9);       /* dec */
         memcpy( obuff + 65, buff + 46, 4);       /* mag */
         strlcpy_err( obuff + 72, "neocp", sizeof( obuff) - 72);
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
   size_t i = 1;
   char obj_name[10], packed[13];

   assert( ' ' == buff[0]);
   while( isdigit( buff[i]))
      i++;
   assert( i < 10);
   if( buff[i] == ' ')           /* numbered object */
      {

      snprintf_err( obj_name, sizeof( obj_name), "(%ld)", atol( buff));
      create_mpc_packed_desig( packed, obj_name);
      memcpy( obuff, packed, 6);
      }
   else                       /* provisional designation */
      {
      assert( 5 == i);
      assert( isalpha( buff[5]));
      assert( isalpha( buff[6]));
      while( buff[i] != ' ')
         i++;
      memcpy( obj_name, buff + 1, i - 1);
      obj_name[i - 1] = '\0';
      create_mpc_packed_desig( packed, obj_name);
      memcpy( obuff, packed, 12);
      }
}

/* The AstDyS/NEODyS .rwo astrometry format provides almost all of
the data of the MPC's 80-column punched-card format.  However,  it
doesn't give the frequency for radar observations.  Fortunately,  since
the radar folks appear to have been fairly consistent in the frequencies
used,  we can (as of early 2014) determine the frequency based on year
and MPC code.  I figured this out based on the accumulated radar data at

http://ssd.jpl.nasa.gov/?grp=ast&fmt=html&radar=             */

static inline double get_radar_frequency( const char *mpc_code, const int year)
{
   double freq_in_mhz;
   int mpc_code_as_int = 100 * mutant_hex_char_to_int( *mpc_code)
                  + atoi( mpc_code + 1);

   switch( mpc_code_as_int)
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
      case 3947:     /* 'd47' = DSS 47 */
         freq_in_mhz = 7159.45;
         break;
      default:                   /* unknown station;  frequency     */
         freq_in_mhz = 5000.;    /* chosen at random using fair die */
         assert( 0);
         break;
      }
   return( freq_in_mhz);
}

/* In transferring/translating satellite offsets from the .rwo format
to the MPC's 80-column format,  there are nuances.  In what follows,
obuff[0] will be a plus or minus.  obuff[1...11] will be the offset,
in kilometers or AU.  ibuff must be set to point to the decimal point
in the .rwo's offset coordinate.  There's some extra weirdness because
if the output is in kilometers,  you have to pad with leading spaces
if there are four or fewer leading digits.   */

static inline void transfer_rwo_satellite_offset( char *obuff,
            const char units, const char *ibuff)
{
   size_t loc = 0, size_out = 11;

   assert( *ibuff == '.');
   ibuff--;
   while( isdigit( *ibuff))
      {
      ibuff--;
      loc++;
      }
   if( *ibuff == '-')
      *obuff++ = '-';
   else
      {
      assert( *ibuff == ' ');
      *obuff++ = '+';
      }
   ibuff++;
   if( units == '1' && loc < 5) /* for km,  must be at least five digits; */
      {                         /* leave some leading spaces */
      obuff += 5 - loc;
      size_out -= 5 - loc;
      }
   while( size_out-- && *ibuff != ' ')
      *obuff++ = *ibuff++;
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

static void xfer_rwo_time_to_mpc( char *obuff, const char *ibuff)
{
   if( ibuff[17] == ' ' || ibuff[16] == ' ') /* six or fewer decimals: */
      memcpy( obuff, ibuff, 17);    /* MPC 'standard' 80-column format */
   else                           /* >6 decimals,  won't fit 'standard'; */
      {                           /* use Find_Orb CYYMMDD.ddddddddd format */
      obuff[0] = 'A' + (ibuff[0] - '1') * 10 + ibuff[1] - '0';
      obuff[1] = ibuff[2];    /* decade */
      obuff[2] = ibuff[3];    /* year */
      obuff[3] = ibuff[5];    /* month, tens */
      obuff[4] = ibuff[6];    /* month, units */
      memcpy( obuff + 5, ibuff + 8, 12);     /* DD.ddddddddd */
      }
}

/* JPL gives DSS locations negative numbers.  These have to be converted
to three-character MPC codes for Find_Orb usage.  Note that four codes,
at present,  lack official MPC codes,  and have d## codes listed in
'rovers.txt'. */

static void _jpl_to_mpc_code( char *mpc_code, const int jpl_code)
{
   const int jpl_codes[] = { -1, -2, -9, -13, -14, -25, -35, -36, -38, -43, -47, -73 };
   const char *mpc_codes = "251 254 256 252 253 257 d35 d36 255 d43 d47 273";
   size_t i;

   for( i = 0; i < sizeof( jpl_codes) / sizeof( jpl_codes[0]); i++)
      if( jpl_code == jpl_codes[i])
         {
         memcpy( mpc_code, mpc_codes + i * 4, 3);
         return;
         }
   assert( 0);       /* shouldn't get here */
}

#define MINIMUM_RWO_LENGTH 117

/* Circa 2019,  AstDyS/NEODyS revised their spacecraft offset lines to
allow 24 bytes per coordinate instead of 12.  Both old and new formats can
be handled without too much trouble.   */

#define OLD_SECOND_LINE_FORMAT   1
#define NEW_SECOND_LINE_FORMAT   2

static int rwo_to_mpc( char *buff, double *ra_bias, double *dec_bias,
                  double *posn_sigma_ra, double *posn_sigma_dec, double *mag_sigma)
{
   int rval = 0, i, second_line = 0;
   const size_t line_len = strlen( buff);
   char obuff[82], second_radar_line[82];
   const size_t old_radar_length = 118, new_radar_length = 125;

   *second_radar_line = '\0';
   if( debug_level > 2)
      debug_printf( "rwo_to_mpc: Input: '%s'\n", buff);
   if( !memcmp( buff + 11, "S s", 3) && buff[27] == '.')
      {
      if( line_len == 75)
         if( buff[42] == '.' && buff[54] == '.' && buff[66] == '.')
            second_line = OLD_SECOND_LINE_FORMAT;
      if( line_len == 111)
         if( buff[46] == '.' && buff[70] == '.' && buff[94] == '.')
            second_line = NEW_SECOND_LINE_FORMAT;
      }
   if( second_line)
      {                       /* satellite second line */
      if( debug_level > 2)
         debug_printf( "Input: '%s'\n", buff);
      memset( obuff, ' ', 80);
      obuff[14] = buff[13];   /* obs type */
      xfer_rwo_time_to_mpc( obuff + 15, buff + 17);
      obuff[32] = buff[34];      /* units specifier */
      memcpy( obuff + 77, buff + line_len - 3, 3);   /* MPC code */
      for( i = 0; i < 3; i++)
         transfer_rwo_satellite_offset( obuff + 34 + i * 12, buff[34],
                     buff + 38 + (4 + i * 12) * second_line);
      rval = 1;
      obuff[80] = '\0';
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
      xfer_rwo_time_to_mpc( obuff + 15, buff + 17);
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
   else if( buff[30] == ':' && buff[33] == ':'
             && (buff[11] == 'R' || buff[11] == 'V')
             && (line_len == old_radar_length || line_len == new_radar_length))
      {                          /* Radar data:  */
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
         if( line_len == old_radar_length)    /* older,  shorter lines store the MPC codes */
            {
            memcpy( obuff + 68, buff + 99, 3);    /* transmitting MPC code */
            memcpy( obuff + 77, buff + 103, 3);    /* receiving MPC code */
            }
         else        /* newer,  ESA-style 125-byte lines store JPL codes */
            {
            _jpl_to_mpc_code( obuff + 68, atoi( buff + 102));
            _jpl_to_mpc_code( obuff + 77, atoi( buff + 107));
            }
         strlcpy_error( second_radar_line, obuff);
         freq_in_mhz = get_radar_frequency( buff + 99, year);
         snprintf_err( obuff + 62, sizeof( obuff) - 62,
                                          "%5.0f", freq_in_mhz);
         second_radar_line[14] = 'r';
         if( buff[11] == 'R')
            {
            val1 *= 2. / SPEED_OF_LIGHT;
            val2 *= 2. / SPEED_OF_LIGHT;
            if( debug_level > 2)
               debug_printf( "%s: round trip %f +/- %f microseconds\n",
                     buff, val1 * 1e+6, val2 * 1e+6);
            snprintf_err( obuff + 32, sizeof( obuff) - 32,
                                             "%13.0f", val1 * 1e+8);
            snprintf_err( second_radar_line + 34, sizeof( second_radar_line) - 34,
                                             "%12.0f", val2 * 1e+9);
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
            snprintf_err( obuff + 48, sizeof( obuff) - 48,
                                   "%13.0f", fabs( val1) * 1e+9);
            obuff[47] = (val1 > 0. ? '+' : '-');
            snprintf_err( second_radar_line + 48, sizeof( second_radar_line) - 48,
                                        "%13.0f", val2 * 1e+9);
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
         strlcpy_err( second_radar_line + 68, obuff + 68, sizeof( second_radar_line) - 68);
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

/* For the 'two-line' MPC formats (roving observers,  spacecraft-based
observations,  radar obs),  the two lines should match in many columns
and have a mismatch in case for column 14 ('V' vs. 'v', 'S' vs. 's',
or 'R' vs. 'r'.)        */

static inline bool matching_lines( const char *line1, const char *line2)
{
   bool rval = (!memcmp( line1, line2, 12)         /* desigs match */
             && !memcmp( line1 + 15, line2 + 15, 16)  /* times match */
             && !memcmp( line1 + 77, line2 + 77, 3)   /* MPC codes match */
             && line1[14] == (line2[14] ^ 32));

   return( rval);
}

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

static int look_for_matching_line( char *iline, char *oline, const size_t oline_size)
{
   static int n_stored = 0;
   static mpc_line *stored_lines = NULL;
   int i;

   *oline = '\0';    /* assume no match found */
   if( !iline)                /* just checking for leftover lines */
      {
      const int rval = n_stored;

      if( n_stored)          /* we _do_ have leftovers;  make an error msg */
         {
         const char *bytes = "VvRrSs";

         for( i = 0; bytes[i]; i++)
            {
            int j = 0, n_matches = 0;

            for( j = 0; j < n_stored; j++)
               if( stored_lines[j][14] == bytes[i])
                  n_matches++;
            if( n_matches)       /* we have this specific sort of error */
               snprintf_append( oline, oline_size, get_find_orb_text( 2040 + i),
                           n_matches);
            }
                    /* now add a general 'look here for info' message */
         strlcat( oline, get_find_orb_text( 2046), oline_size);
         debug_printf( "%d unmatched satellite/roving observer lines:\n",
                           n_stored);
         for( i = 0; i < n_stored; i++)
            debug_printf( "%s\n", stored_lines[i]);
         n_stored = 0;
         }
      if( stored_lines)
         {
         free( stored_lines);
         stored_lines = NULL;
         }
      return( rval);
      }

   for( i = 0; i < n_stored && !*oline; i++)
      if( matching_lines( stored_lines[i], iline))
         {        /* we have a match */
         memcpy( oline, stored_lines[i], sizeof( mpc_line));
         n_stored--;
               /* Move last line into place formerly used by the */
               /* now-matched line: */
         memcpy( stored_lines[i], stored_lines[n_stored], sizeof( mpc_line));
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

/* Occasionally,  observations are duplicates,  except that one has a
blank reference or magnitude or other datum and the other doesn't. In
this function,  we'll copy reference that _is_ given over the blank
one, and the mag that's given over the zero mag,  and so forth.  Then
we check to see if that's caused the observations to become the same.
If it has, we consider the difference to have been "corrected".    */

static void correct_differences( OBSERVE *obs1, const OBSERVE *obs2)
{
   if( obs1->reference[0] == ' ' && obs2->reference[0] != ' ')
      strcpy( obs1->reference, obs2->reference);
   if( obs1->obs_mag == BLANK_MAG && obs2->obs_mag != BLANK_MAG)
      {
      obs1->obs_mag = obs2->obs_mag;
      obs1->mag_sigma = obs2->mag_sigma;
      obs1->mag_precision = obs2->mag_precision;
      obs1->mag_band = obs2->mag_band;
      }
   if( obs1->astrometric_net_code == ' ' && obs2->astrometric_net_code != ' ')
      obs1->astrometric_net_code = obs2->astrometric_net_code;
   if( obs1->mag_band == ' ' && obs2->mag_band != ' ')
      obs1->mag_band = obs2->mag_band;
   if( obs1->note1 == ' ' && obs2->note1 != ' ')
      obs1->note1 = obs2->note1;
   if( obs1->discovery_asterisk == ' ' && obs2->discovery_asterisk == '*')
      obs1->discovery_asterisk = '*';
}

/* For certain forms of orbit determination,  having radar observations right
in the middle of the data can cause trouble.  Let's say you try to apply
the method of Gauss,  but one of your observations is a radar one.  To evade
that issue,  we sort observations to temporarily put all the radar data at
the end.  When we're done,  we re-sort "conventionally" (by date).  */

int compare_observations( const void *a, const void *b, void *context)
{
   const OBSERVE *obs1 = (const OBSERVE *)a;
   const OBSERVE *obs2 = (const OBSERVE *)b;
   int rval = FMEMCMP( obs1->mpc_code, obs2->mpc_code, 3);

   if( context)
      {
      int *sort_type = (int *)context;

      if( *sort_type == SORT_OBS_BY_CODE_THEN_DATE && rval)
         return( rval);
      if( *sort_type == SORT_OBS_RADAR_LAST)
         {
         if( obs1->note2 == 'R' && obs2->note2 != 'R')
            return( 1);
         if( obs1->note2 != 'R' && obs2->note2 == 'R')
            return( -1);
         if( (obs1->flags ^ obs2->flags) & OBS_DONT_USE)
            return( (obs1->flags & OBS_DONT_USE) ? 1 : -1);
         }
      }
   if( obs1->jd < obs2->jd)
      rval = -1;
   else if( obs1->jd > obs2->jd)
      rval = 1;
   if( !rval)
      rval = obs1->note1 - obs2->note1;
   if( !rval)
      rval = strcmp( obs1->reference, obs2->reference);
   if( !rval && obs1->ra != obs2->ra)
      rval = (obs1->ra > obs2->ra ? 1 : -1);
   if( !rval && obs1->dec != obs2->dec)
      rval = (obs1->dec > obs2->dec ? 1 : -1);
   return( rval);
}

/* If two observations are within a threshhold of each other in time,
they may be duplicates,  with one (say) recorded to five places and
the other to six.  Or one in decimal days (80-column format) and
the other in seconds (ADES).

The threshhold may depend on the type of observation.  Video,  for
example,  may provide a slew of observations within 1/30 second of
each other.  Such nuances aren't currently implemented,  but they
are why pointers to the full observation data are passed,  instead
of just the observed times. */

static bool times_very_close( OBSERVE *obs1, OBSERVE *obs2)
{
   const double threshhold = 3. / seconds_per_day;

   return( fabs( obs1->jd - obs2->jd) < threshhold);
}

/* If we have duplicate observations,  but their RA/decs are within
(say) two arcseconds,  they may plausibly work for initial orbit
determination.  Otherwise,  they'll probably just throw off IOD
and should be marked as OBS_DONT_USE and excluded. */

static void exclude_unusable_duplicate_obs( OBSERVE *obs, int n_obs)
{
   int i;

   while( n_obs)
      {
      double min_ra = 0., max_ra = 0., min_dec = 0., max_dec = 0.;
      const double tolerance = 2. * (PI / 180.) / 3600.;    /* two arcsec */
      int j;

      i = 1;
      while( i < n_obs && times_very_close( obs, obs + i)
                        && !strcmp( obs[i].mpc_code, obs->mpc_code))
         {
         const double d_ra = centralize_ang( obs[i].ra - obs->ra);
         const double d_dec = centralize_ang( obs[i].dec - obs->dec);

         if( min_ra > d_ra)
            min_ra = d_ra;
         if( max_ra < d_ra)
            max_ra = d_ra;
         if( min_dec > d_dec)
            min_dec = d_dec;
         if( max_dec < d_dec)
            max_dec = d_dec;
         i++;
         }
      if( max_ra - min_ra > tolerance || max_dec - min_dec > tolerance)
         for( j = 0; j < i; j++)
            {
            obs[j].flags |= OBS_DONT_USE;
            obs[j].is_included = 0;
            }
      obs += i;
      n_obs -= i;
      }
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
      if( !times_very_close( obs + i, obs + i - 1) || strcmp( obs[i].mpc_code,
                                                obs[i - 1].mpc_code))
         obs[j++] = obs[i];
      else if( memcmp( obs + i, obs + i - 1, sizeof( OBSERVE)))
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
            comment_observation( obs + j - 1, "NearDup  ");
            }
         else
            {                    /* Is still really a duplicate,  but the  */
            obs[j] = obs[i];     /* differences couldn't be fixed          */
            if( tolower( obs[i].note2) != 'x' && tolower( obs[j].note2) != 'x')
               {
               comment_observation( obs + j - 1, "FixDup   ");
               comment_observation( obs + j,     "FixDup   ");
               }
            j++;
            }
         }
      else
         comment_observation( obs + j - 1, "Duplicate");
   exclude_unusable_duplicate_obs( obs, j);
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
   extern int available_sigmas;

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
   if( obs_details)
      {
      free_observation_details( obs_details);
      obs_details = NULL;
      }
   if( _ades_ids_stack)
      {
      destroy_stack( _ades_ids_stack);
      _ades_ids_stack = NULL;
      }
   if( rovers)
      {
      free( rovers);
      rovers = NULL;
      n_rovers = 0;
      }
   available_sigmas = NO_SIGMAS_AVAILABLE;
   return( 0);
}

inline bool is_satellite_obs( const OBSERVE *obs)
{
   return( obs->second_line && obs->second_line[14] == 's');
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

/* Observation sigmas,  etc. can be set by inserting lines such as

#Posn sigma .3
#Mag sigma .23

   etc. (see below).  The following function allows one to use COM lines
such as

COM Posn sigma .3
COM Mag sigma .23

   instead.  MPC will complain about the former,  but will ignore the latter.
*/

static void convert_com_to_pound_sign( char *buff)
{
   if( !memcmp( buff, "COM ", 4))
      {
      *buff = '#';
      memmove( buff + 1, buff + 4, strlen( buff + 3));
      }
}

/* (691) Spacewatch has a habit of providing near-duplicate observations,
one "hand-measured",  the other not,  like this :

     P10pKck  C2015 12 06.25527 04 35 11.84 +14 22 25.4          21.0 VoNEOCP691
     P10pKck HC2015 12 06.25527 04 35 11.85 +14 22 25.4          21.4 VoNEOCP691

   Normally,  two observations with identical times and observatories would
both be discarded.  In this case,  we want to toss out the first (non-hand
measured) observation only.

   More recently,  there have been similar situations,  except one has no
note and the other has a 'K' ('stacked image').  In that case,  we're supposed
to ignore the stacked one and keep the blank-note one.  */

inline int spacewatch_duplication( OBSERVE FAR *obs)
{
   int rval = 0;

   if( !strcmp( obs->mpc_code, "691") && obs[0].note1 == ' ')
      {
      if( obs[1].note1 == 'H')
         {
         rval = 1;
         obs->is_included = 0;
         obs->flags |= OBS_DONT_USE;
         comment_observation( obs, "Replaced");
         }
      else if( obs[1].note1 == 'K')
         {
         rval = 1;
         obs[1].is_included = 0;
         obs[1].flags |= OBS_DONT_USE;
         comment_observation( obs + 1, "Replaced");
         }
      }
   return( rval);
}

static inline void check_for_star( const OBSERVE *obs, const int n_obs)
{
   double min_ra, max_ra, min_dec, max_dec;
   int i;
   const double tolerance = 2. * (PI / 180.) / 3600.;

   min_ra = max_ra = obs[0].ra;
   min_dec = max_dec = obs[0].dec;
   for( i = 1; i < n_obs; i++)
      {
      if( min_ra > obs[i].ra)
         min_ra = obs[i].ra;
      else if( max_ra < obs[i].ra)
         max_ra = obs[i].ra;
      if( min_dec > obs[i].dec)
         min_dec = obs[i].dec;
      else if( max_dec < obs[i].dec)
         max_dec = obs[i].dec;
      }
   if( max_ra - min_ra < tolerance && max_dec - min_dec < tolerance)
      generic_message_box( get_find_orb_text( 2001), "!");
}

/* ADES sigmas are converted into a punched-card compatible format such as

COM Sigmas 0.34x0.42,-0.15 m:0.22 t:1.3

   In the above case,  the RA uncertainty is 0.34 arcsec,  dec 0.42,  with a
correlation of -0.15.  The magnitude sigma is 0.22;  the time sigma, 1.3.
Not all of these must be present,  except that a single RA/dec sigma does
have to be given. */

static inline int extract_ades_sigmas( const char *buff,
         double *posn1, double *posn2, double *theta,
         double *mag_sigma, double *time_sigma, double *unc_time)
{
   int loc, rval = -1;

   if( sscanf( buff, "%lf%n", posn1, &loc) == 1)
      {
      const char *tptr;

      rval = 0;
      *theta = 0.;
      if( buff[loc] == 'x')
         {
         *posn2 = atof( buff + loc + 1);
         while( buff[loc] > ' ' && buff[loc] != ',')
            loc++;
         if( buff[loc] == ',')
            {
            const double correlation = atof( buff + loc + 1);

            convert_ades_sigmas_to_error_ellipse( *posn1, *posn2,
                              correlation, posn1, posn2, theta);
            *theta = PI / 2. + *theta;
            }
         }
      else        /* circular position error */
         *posn2 = *posn1;
      tptr = strstr( buff + loc, "m:");
      if( tptr)
         *mag_sigma = atof( tptr + 2);
      tptr = strstr( buff + loc, "t:");
      if( tptr)
         *time_sigma = atof( tptr + 2) / seconds_per_day;
      tptr = strstr( buff + loc, "u:");
      if( tptr)
         *unc_time   = atof( tptr + 2) / seconds_per_day;
      }
   if( rval)
      debug_printf( "Malformed ADES sigma: '%s'\n", buff);
   return( rval);
}

/* If photometry is provided,  we count both comet and asteroid obs
to decide if we ought to use a comet or asteroid magnitude formula.
(If it's not provided,  a default judgment is made based on the desig.) */

static void reset_object_type( const OBSERVE *obs, const int n_obs)
{
   int i, n_asteroid = 0, n_nuclear = 0, n_total = 0;

   for( i = 0; i < n_obs; i++)
      if( obs[i].obs_mag != BLANK_MAG && obs[i].is_included)
         {
         if( obs[i].mag_band == 'N')
            n_nuclear++;
         else if( obs[i].mag_band == 'T')
            n_total++;
         else
            n_asteroid++;
         }
   if( n_total + n_nuclear > n_asteroid)
      object_type = OBJECT_TYPE_COMET;
   else if( is_sungrazing_comet( obs, n_obs))
      object_type = OBJECT_TYPE_COMET;
   else if( n_asteroid)
      object_type = OBJECT_TYPE_ASTEROID;
   if( object_type == OBJECT_TYPE_COMET)
      {
      extern char default_comet_magnitude_type;

      if( default_comet_magnitude_type == 'N' && !n_nuclear && n_total)
         default_comet_magnitude_type = 'T';
      if( default_comet_magnitude_type == 'T' && n_nuclear && !n_total)
         default_comet_magnitude_type = 'N';
      }
}

bool nighttime_only( const char *mpc_code)
{
   return( strstr( get_environment_ptr( "DAYTIME_OBS_OK1"), mpc_code) == NULL);
}

#define INSUFFICIENT_PRECISION_MAG      1
#define INSUFFICIENT_PRECISION_TIME     2
#define INSUFFICIENT_PRECISION_POSN1    4
#define INSUFFICIENT_PRECISION_POSN2    8

static bool _is_synthetic_obs( const OBSERVE *obs)
{
   return( !strcmp( obs->reference, "Synth")
            || !strcmp( obs->reference, "neocp")
            || !strcmp( obs->reference, "Dummy"));
}

/* Uncertainties on time,  magnitude,  and the error ellipse all will have
default values.  If ADES or Tholen or .rwo uncertainties have been set,  or
uncertainties in columns 57-65,  we use them.  (Such uncertainties apply to
_only_ one observation.)  Global uncertainties can be set using COM Time
sigma,  COM Posn sigma or COM Mag sigma;  if observation-specific sigmas
haven't been set,  we use those global uncertainties.

   I'm calling the various flavors of single-observation uncertainties
"ADES sigmas",  because (a) the behavior for any of those methods is the
same -- the sigma(s) are applied to the next observation and then forgotten;
(b) ADES will probably be the most common case;  (c) variable names such as
single_observation_posn_sigma_1 looked too long to me.     */

#define SET_SIGMA( sigma, ades_sigma, override_sigma)   \
       { if( ades_sigma) sigma = ades_sigma;  else if( override_sigma) sigma = override_sigma; }


   /* By default,  Find_Orb will only handle arcs up to 200 years */
   /* long.  If the arc is longer than that,  observations will be */
   /* dropped to get an arc that fits.  The max arc length can be */
   /* adjusted in 'environ.dat'.                                  */
double maximum_observation_span = 200.;

int sanity_check_observations = 1;
bool use_sigmas = true;
extern int is_interstellar;
static double _overall_obj_alt_limit, _overall_sun_alt_limit;

OBSERVE FAR *load_observations( FILE *ifile, const char *packed_desig,
                           const int n_obs)
{
   const double days_per_year = 365.25;
   char buff[650], mpc_code_from_neocp[4], desig_from_neocp[15];
   char obj_name[80], curr_ades_ids[100];
   OBSERVE FAR *rval;
   bool including_obs = true;
   int i, n_fixes_made = 0, n_future_obs = 0;
   unsigned line_no = 0;
   unsigned n_below_horizon = 0, n_in_sunlight = 0;
   unsigned n_spurious_matches = 0;
   unsigned n_sat_obs_without_offsets = 0;
   double override_posn_sigma_1 = 0., ades_posn_sigma_1 = 0.;  /* in arcsec */
   double override_posn_sigma_2 = 0., ades_posn_sigma_2 = 0.;
   double override_posn_sigma_theta = 0., ades_posn_sigma_theta = 0.;
   double override_mag_sigma = 0., ades_mag_sigma = 0.;   /* in mags */
   double override_time_sigma = 0., ades_time_sigma = 0.;  /* in seconds */
   double unc_time = 0.;
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
   bool is_fcct14_or_vfcc17_data = false;
   void *ades_context;
   int spacecraft_offset_reference = 399;    /* default is geocenter */
   double spacecraft_vel[3];
   static int suppress_private_obs = -1;
   int count_satellite_coord_errors[N_SATELL_COORD_ERRORS];

   memset( count_satellite_coord_errors, 0, sizeof( count_satellite_coord_errors));
   move_add_nstr( 1, 2, "Loading observations", -1);
   refresh_console( );
   *desig_from_neocp = '\0';
   strcpy( mpc_code_from_neocp, "500");   /* default is geocenter */
   neocp_file_type = NEOCP_FILE_TYPE_UNKNOWN;
   get_object_name( obj_name, packed_desig);
   rval = (OBSERVE FAR *)FCALLOC( n_obs + 1, sizeof( OBSERVE));
   if( !rval)
      return( NULL);
   input_coordinate_epoch = 2000.;
            /* Start out assuming asteroid-type magnitudes.  The name */
            /* may tell you it's really a comet,  and the orbit may   */
            /* tell you it's an artsat.                               */
   object_type = OBJECT_TYPE_ASTEROID;
   is_interstellar = 0;
   n_rovers = 0;
   if( rovers)
      {
      free( rovers);
      rovers = NULL;
      }
   if( !obs_details)
      obs_details = init_observation_details( );
   ades_context = init_ades2mpc( );
   memset( buff, 0, 13);         /* suppress spurious Valgrind messages */
   *curr_ades_ids = '\0';
   for( i = 0; i < 3; i++)
      spacecraft_vel[i] = 0.;
   i = 0;
   while( fgets_with_ades_xlation( buff, sizeof( buff), ades_context, ifile)
                  && i != n_obs)
      {
      int is_rwo = 0, fixes_made = 0;
      char original_packed_desig[13];
      size_t ilen = strlen( buff);
      double jd;

      line_no++;
      if( *buff == '<')
         remove_html_tags( buff);
      if( !strncmp( buff, "errmod  = 'fcct14'", 18))
         is_fcct14_or_vfcc17_data = true;
      if( !strncmp( buff, "errmod  = 'vfcc17'", 18))
         is_fcct14_or_vfcc17_data = true;
      if( debug_level > 2)
         debug_printf( "Line %d: %s\n", line_no, buff);
      if( get_neocp_data( buff, desig_from_neocp, mpc_code_from_neocp))
         if( !i && debug_level)
            debug_printf( "Got NEOCP data\n");
      if( ilen == 75 || ilen == 111 || ilen >= MINIMUM_RWO_LENGTH)
         {
         is_rwo = rwo_to_mpc( buff, &rval[i].ra_bias, &rval[i].dec_bias,
                &ades_posn_sigma_1, &ades_posn_sigma_2, &ades_mag_sigma);
         if( is_rwo && !i && debug_level)
            debug_printf( "Got .rwo data\n");
         }
      if( fixing_trailing_and_leading_spaces)
         fixes_made = fix_up_mpc_observation( buff, NULL);
      original_packed_desig[12] = '\0';
      memcpy( original_packed_desig, buff, 12);
      xref_designation( buff);
      add_line_to_observation_details( obs_details, buff);
      jd = observation_jd( buff);
      if( is_in_range( jd) && !compare_desigs( packed_desig, buff))
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
               else if( !look_for_matching_line( buff, second_line, sizeof( second_line)))
                  observation_is_good = false;
               if( observation_is_good)
                  {
                  rval[i].ref_center = spacecraft_offset_reference;
                  parse_observation( rval + i, buff);
                  }
               }

            if( buff[14] == 'S' && observation_is_good)
               {      /* we did find the "matching" line: */
               double vect[3];
               int j, error_code;

               rval[i].satellite_obs = (char)(second_line[32] - '0');
               error_code = get_satellite_offset( second_line, vect);
               assert( error_code <= 0 && error_code > -N_SATELL_COORD_ERRORS);
               count_satellite_coord_errors[-error_code]++;
               if( error_code)
                  {
                  char tbuff[6];

                  rval[i].flags |= OBS_DONT_USE;
                  rval[i].is_included = 0;
                  strcpy( tbuff, "?off ");
                  tbuff[4] = (char)( '0' - error_code);
                  tbuff[5] = '\0';
                  comment_observation( rval + i, tbuff);
                  n_bad_satellite_offsets++;
                  debug_printf( "Error code %d; offending line was:\n%s\n",
                     error_code, second_line);
                  }
               if( !spacecraft_vel[0] && !spacecraft_vel[1]
                              && !spacecraft_vel[2])
                  rval[i].flags |= OBS_NO_VELOCITY;
               for( j = 0; j < 3; j++)
                  {
                  rval[i].obs_posn[j] += vect[j];
                  rval[i].obs_vel[j] += spacecraft_vel[j] * seconds_per_day / AU_IN_KM;
                  spacecraft_vel[j] = 0.;
                  }
               rval[i].second_line = (char *)malloc( 81);
               strcpy( rval[i].second_line, second_line);
               }
            else if( buff[14] == 'R' && observation_is_good)
               {      /* we did find the "matching" line: */
//             fix_radar_time( buff);
               rval[i].second_line = (char *)malloc( 81 * 2);
               strcpy( rval[i].second_line, second_line);
               strcpy( rval[i].second_line + 81, buff);
               }
            else if( buff[14] == 'V' && observation_is_good)
               {
               double rho_sin_phi, rho_cos_phi;
               const double rlon = atof( second_line + 34);
               const double rlat = atof( second_line + 45);
               const double ralt = atof( second_line + 56);
               const double dist_tol = 0.0003;  /* 3e-4 degrees = about 30 meters */
               int idx = 0;

               xref_designation( second_line);
               while( idx < n_rovers &&
                            (fabs( rlon - rovers[idx].lon) > dist_tol ||
                            fabs( rlat - rovers[idx].lat) > dist_tol ||
                            fabs( ralt - rovers[idx].alt) > 30.))
                  idx++;
               if( idx == n_rovers)    /* got a new rover */
                  {
                  n_rovers++;
                  assert( idx < 1665);    /* can't handle more at the mo */
                  rovers = (rover_t *)realloc( rovers, n_rovers * sizeof( rover_t));
                  rovers[idx].lat = rlat;
                  rovers[idx].lon = rlon;
                  rovers[idx].alt = ralt;
                  }
               if( idx)    /* not our first,  default (247) rover */
                  {
                  if( idx >= 53)
                     {
                     idx -= 53;
                     if( second_line[78] == '4')   /* it's a (247) code */
                        second_line[78] = 'A' + idx / 62;
                     else                          /* it's a (270) code */
                        second_line[78] = 'a' + idx / 62;
                     rval[i].mpc_code[1] = second_line[78];
                     second_line[79] = int_to_mutant_hex_char( idx % 62);
                     }
                  else if( idx < 27)
                     second_line[79] = 'a' + idx - 1;
                  else    /* if( idx < 53) */
                     second_line[79] = 'A' + idx - 27;
                  rval[i].mpc_code[2] = second_line[79];
                  }
               lat_alt_to_parallax( rlat * PI / 180., ralt,
                                    &rho_cos_phi, &rho_sin_phi, 3);
               compute_observer_loc( rval[i].jd, 3, rho_cos_phi, rho_sin_phi,
                                      rlon * PI / 180., rval[i].obs_posn);
               compute_observer_vel( rval[i].jd, 3, rho_cos_phi, rho_sin_phi,
                                      rlon * PI / 180., rval[i].obs_vel);
               set_obs_vect( rval + i);
               rval[i].second_line = (char *)malloc( 81);
               strcpy( rval[i].second_line, second_line);
               }
            if( rval[i].reference[0] == '!' && suppress_private_obs == -1)
               {
               const char *suppression = get_environment_ptr( "PRIVATE_OBS");

               if( *suppression)
                  suppress_private_obs = atoi( suppression);
               else
                  {
                  const int c = generic_message_box( get_find_orb_text( 2052), "o");

                  if( c == 'y' || c == 'Y')     /* first time we see 'private' obs, */
                     suppress_private_obs = 0;  /* ask if user wants them included */
                  else
                     suppress_private_obs = 1;
                  }
               }
            if( rval[i].reference[0] == '!' && suppress_private_obs)
               observation_is_good = false;
            if( observation_is_good)
               {
               const double radians_per_arcsec = PI / (180. * 3600.);
               double mag_sigma = 0., time_sigma = 0.;
               double posn_sigma_1, posn_sigma_2;
               double posn_sigma_theta = 0.;
                    /* If we want all observations to have the same sigma, */
                    /* we use a nonexistent MPC code.  "Unknown" codes will */
                    /* just get the default sigma assigned at the end of    */
                    /* 'sigma.txt' (q.v.)                                   */
               const char *mpc_code_to_use = (use_sigmas ? rval[i].mpc_code : "***");

               if( rval[i].jd < 2e+6 || rval[i].jd > 3e+6)
                  debug_printf( "Weird obs JD %f:\n%s\n", rval[i].jd, buff);
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

                  SET_SIGMA( mag_sigma, ades_mag_sigma, override_mag_sigma);
                  SET_SIGMA( time_sigma, ades_time_sigma, override_time_sigma);
                  SET_SIGMA( posn_sigma_1, ades_posn_sigma_1, override_posn_sigma_1);
                  SET_SIGMA( posn_sigma_2, ades_posn_sigma_2, override_posn_sigma_2);
                  SET_SIGMA( posn_sigma_theta, ades_posn_sigma_theta, override_posn_sigma_theta);

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
               if( *curr_ades_ids)
                  {
                  if( !_ades_ids_stack)
                     _ades_ids_stack = create_stack( 2000);
                  rval[i].ades_ids = (char *)stack_alloc( _ades_ids_stack,
                           strlen( curr_ades_ids) + 1);
                  strcpy( rval[i].ades_ids, curr_ades_ids);
                  *curr_ades_ids = '\0';
                  }

               if( !strcmp( rval[i].reference, "neocp"))
                  {
                  posn_sigma_1 = 1.5;     /* NEOCP ephems are given to */
                  posn_sigma_2 = 1.0;     /* 0.1s in RA, 1" in dec */
                  }
                           /* The observation data's precision has already been */
                           /* used to figure out a minimum sigma;  e.g.,  a mag */
                           /* of '16.3' results in a sigma of 0.1.  If the above */
                           /* methods return smaller values,  we stick with the */
                           /* sigma determined from the # of places given. */
                           /* If we're enforcing uniform sigmas,  though,  we use */
                           /* those sigmas whether they "make sense" or not. */
               if( mag_sigma * 1.1 > rval[i].mag_sigma || !use_sigmas)
                  rval[i].mag_sigma = mag_sigma;
//             else if( rval[i].obs_mag != BLANK_MAG)
//                insufficient_precision |= INSUFFICIENT_PRECISION_MAG;
               if( time_sigma * 1.1 > rval[i].time_sigma || !use_sigmas)
                  rval[i].time_sigma = time_sigma;
//             else if( rval[i].note2 != 'X')
//                insufficient_precision |= INSUFFICIENT_PRECISION_TIME;
               if( posn_sigma_1 * 1.1 > rval[i].posn_sigma_1 || !use_sigmas)
                  rval[i].posn_sigma_1 = posn_sigma_1;
//             else if( rval[i].note2 != 'X')
//                insufficient_precision |= INSUFFICIENT_PRECISION_POSN1;
               if( posn_sigma_2 * 1.1 > rval[i].posn_sigma_2 || !use_sigmas)
                  rval[i].posn_sigma_2 = posn_sigma_2;
//             else if( rval[i].note2 != 'X')
//                insufficient_precision |= INSUFFICIENT_PRECISION_POSN2;
               rval[i].posn_sigma_theta = posn_sigma_theta;
               rval[i].unc_time = unc_time;
               if( !including_obs)
                  rval[i].is_included = 0;
               if( apply_debiasing)
                  {
                  rval[i].ra  -= rval[i].ra_bias  * radians_per_arcsec / cos( rval[i].dec);
                  rval[i].dec -= rval[i].dec_bias * radians_per_arcsec;
                  }
               if( is_fcct14_or_vfcc17_data && is_rwo)
                  rval[i].flags |= OBS_ALREADY_CORRECTED_FOR_OVEROBSERVING;
               if( rval[i].note2 == 'n')        /* video observations:  assume */
                  rval[i].time_sigma = 0.01 / seconds_per_day;  /* 10ms sigma */
               set_obs_vect( rval + i);
               if( !rval[i].is_included)
                  rval[i].flags |= OBS_DONT_USE;
               rval[i].obs_details = get_code_details( obs_details, rval[i].mpc_code);
               if( fixes_made)
                  {
                  char comment[15];

                  snprintf_err( comment, sizeof( comment), "FixMe%d ", fixes_made);
                  comment_observation( rval + i, comment);
                  n_fixes_made++;
                  }
               if( rval[i].packed_id[4] == 'I' &&
                    rval[i].packed_id[0] == '0' && rval[i].packed_id[1] == '0')
                  is_interstellar = 1;
               i++;
               spacecraft_offset_reference = 399;  /* default to geocentric offsets */
               }
            }
         }
      if( is_in_range( jd) && buff[14] != 'S')     /* Sigmas from ADES or Dave Tholen apply to only */
         {                 /* one observation.  Zero 'em out after that use : */
         ades_posn_sigma_1 = ades_posn_sigma_2 = ades_mag_sigma = 0.;
         ades_posn_sigma_theta = ades_time_sigma = unc_time = 0.;
         }
      override_time = 0.;
      override_ra = override_dec = override_mag = -100.;
      convert_com_to_pound_sign( buff);
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
         else if( !memcmp( buff, "#RA/dec ", 8))
            {
            if( 3 == sscanf( buff + 8, "%lf %lf %lf", &override_ra,
                                    &override_dec, &override_time))
               override_time += J2000;
            }
         else if( !memcmp( buff, "#full mag ", 9))
            override_mag = atof( buff + 9);
         else if( !memcmp( buff, "#Sigmas ", 8))
            extract_ades_sigmas( buff + 8, &ades_posn_sigma_1,
                        &ades_posn_sigma_2,
                        &ades_posn_sigma_theta,
                        &ades_mag_sigma, &ades_time_sigma, &unc_time);
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
         else if( !memcmp( buff, "#sun_alt_limit ", 15))
            _overall_sun_alt_limit = atof( buff + 15);
         else if( !memcmp( buff, "#reset_debug ", 13))
             debug_level = atoi( buff + 13);
         else if( !memcmp( buff, "#obj_alt_limit ", 15))
            _overall_obj_alt_limit = atof( buff + 15);
         else if( !memcmp( buff, "#include_obs", 12))
            including_obs = true;
         else if( !memcmp( buff, "#toffset", 7))
            observation_time_offset = atof( buff + 8) / seconds_per_day;
         else if( !memcmp( buff, "#interstellar", 12))
            is_interstellar = 1;
         else if( !memcmp( buff, "#time ", 6))
            override_time = get_time_from_string( 0, buff + 6,
                              CALENDAR_JULIAN_GREGORIAN, NULL);
                  /* Above allows one to reset the time of the preceding obs */
         else if( !memcmp( buff, "#comet", 6))
            object_type = OBJECT_TYPE_COMET;
         else if( !strcmp( buff, "#ignore obs"))
            while( fgets_trimmed( buff, sizeof( buff), ifile)
                     && !strstr( buff, "end ignore obs"))
               ;     /* deliberately empty loop */
         else if( !memcmp( buff, "#Offset center ", 15))
            spacecraft_offset_reference = atoi( buff + 15);
         else if( !memcmp( buff, "#vel (km/s) ", 12))
            {
            int n_scanned = sscanf( buff + 30, "%lf %lf %lf", spacecraft_vel,
                        spacecraft_vel + 1, spacecraft_vel + 2);

            assert( 3 == n_scanned);
            equatorial_to_ecliptic( spacecraft_vel);
            }
         else if( !memcmp( buff, "#IDs ", 5))
            strlcpy_err( curr_ades_ids, buff + 5, sizeof( curr_ades_ids));
         else if( !memcmp( buff, "#'getradar' version", 19))
            {
            int j = 0;  /* we're replacing/updating the radar observations */
                        /* using JPL data.  So remove the old MPC radar. */
            while( j < i)
               if( rval[j].note2 == 'R')
                  rval[j] = rval[--i];    /* remove the radar observation */
               else
                  j++;
            }
         }
      }
   free_ades2mpc_context( ades_context);
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
   if( i && rval[i - 1].jd - rval[0].jd > maximum_observation_span * days_per_year)
      {
      int j = 0;

      snprintf_err( buff, sizeof( buff), get_find_orb_text( 2007),
                     (rval[i - 1].jd - rval[0].jd) /  days_per_year,
                     maximum_observation_span);
      generic_message_box( buff, "!");
      while( rval[i - 1].jd - rval[j].jd > maximum_observation_span * days_per_year)
         j++;
      i -= j;
      n_obs_actually_loaded -= j;
      memmove( rval, rval + j, i * sizeof( OBSERVE));
      }

   n_duplicate_obs_found = n_obs_actually_loaded - i;
   n_obs_actually_loaded = fix_radar_obs( rval, i);

   monte_carlo_object_count = 0;
   n_monte_carlo_impactors = 0;
   if( look_for_matching_line( NULL, buff, sizeof( buff)))
      generic_message_box( buff, "!");
   if( n_fixes_made)
      {
      snprintf_err( buff, sizeof( buff), get_find_orb_text( 2028),
                     n_fixes_made);
      generic_message_box( buff, "!");
      }
   if( n_duplicate_obs_found)
      {
      snprintf_err( buff, sizeof( buff), get_find_orb_text( 2009),
                     n_duplicate_obs_found);
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "!");
      }
   if( n_bad_satellite_offsets)
      {
      *buff = '\0';
      for( i = 0; i < N_SATELL_COORD_ERRORS; i++)
         if( count_satellite_coord_errors[i])
            snprintf_append( buff, sizeof( buff), get_find_orb_text( i + 2090),
                                   count_satellite_coord_errors[i]);
      strlcat_error( buff, get_find_orb_text( 2010));
      generic_message_box( buff, "!");
      }

   if( n_parse_failures)
      {
      snprintf_err( buff, sizeof( buff), get_find_orb_text( 2011),
                     n_parse_failures);
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "!");
      }

   if( sanity_check_observations)
      for( i = 0; i < n_obs_actually_loaded; i++)
         if( rval[i].note2 != 'R')        /* no way to sanity-test */
            {                             /* radar obs yet         */
            DPT obj_alt_az, sun_alt_az;
            const bool is_synthetic = _is_synthetic_obs( rval + i);

            if( !is_synthetic
                     && !get_obs_alt_azzes( rval + i, &sun_alt_az, &obj_alt_az)
                     && sun_alt_az.x > -90.)
               {                          /* I.e., not flagged as meaningless */

               if( sun_alt_az.y > _overall_sun_alt_limit
                                     && nighttime_only( rval[i].mpc_code))
                  {
                  rval[i].is_included = 0;
                  rval[i].flags |= OBS_DONT_USE;
                  n_in_sunlight++;
                  comment_observation( rval + i, "Daylit");
                  }
               if( obj_alt_az.y < _overall_obj_alt_limit)
                  {
                  rval[i].is_included = 0;
                  rval[i].flags |= OBS_DONT_USE;
                  n_below_horizon++;
                  comment_observation( rval + i, "Horizon");
                  }
               }
            }
   if( n_below_horizon || n_in_sunlight)
      {
      *buff = '\0';
      if( n_below_horizon)
         snprintf_err( buff, sizeof( buff), get_find_orb_text( 2012),
                        n_below_horizon);
      if( n_in_sunlight)
         snprintf_err( buff, sizeof( buff), get_find_orb_text( 2013),
                        n_in_sunlight);
      strlcat_error( buff, get_find_orb_text( 2014));
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "!");
      }
   for( i = 1; i < n_obs_actually_loaded; i++)
      if( times_very_close( rval + i, rval + i - 1) && !strcmp( rval[i].mpc_code, rval[i - 1].mpc_code))
         if( toupper( rval[i].note2) != 'X' && toupper( rval[i - 1].note2) != 'X')
            if( !spacewatch_duplication( rval + i - 1))
               {
               comment_observation( rval + i, "Duplicate");
               if( rval[i].ra == rval[i - 1].ra
                        && rval[i].dec == rval[i - 1].dec
                        && rval[i].obs_mag == rval[i - 1].obs_mag)
                  n_almost_duplicates_found++;
               else
                  {
                  comment_observation( rval + i - 1, "Duplicate");
                  n_spurious_matches++;
                  }
               }
   if( n_spurious_matches)
      {
      snprintf_err( buff, sizeof( buff), get_find_orb_text( 2015),
               n_spurious_matches);
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "!");
      }
   if( n_almost_duplicates_found)
      {
      snprintf_err( buff, sizeof( buff), get_find_orb_text( 2016),
               n_almost_duplicates_found);
      debug_printf( "%s:\n", rval->packed_id);
      generic_message_box( buff, "!");
      }
   for( i = 0; i < n_obs_actually_loaded; i++)
      if( rval[i].flags & OBS_NO_OFFSET)
         if( tolower( rval[i].note2) != 'x')
            n_sat_obs_without_offsets++;
   if( n_sat_obs_without_offsets)
      {
      snprintf_err( buff, sizeof( buff), get_find_orb_text( 2017),
                        n_sat_obs_without_offsets);
      generic_message_box( buff, "!");
      }
   i = n_obs_actually_loaded;
   while( i-- >= 0 && rval[i].jd > current_jd( ))
      if( !_is_synthetic_obs( rval + i))
         n_future_obs++;    /* warn of non-synthetic obs in the future */
   if( n_future_obs)
      {
      snprintf_err( buff, sizeof( buff), get_find_orb_text( 2018),
                     n_future_obs);
      generic_message_box( buff, "!");
      }

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
   if( n_obs > 1)
      {
      check_for_star( rval, n_obs);
      for( i = 0; i < n_obs_actually_loaded; i++)
         if( rval[i].note2 == 'x')     /* check for deleted satellite observation */
            {
            int j = -1;

            if( i > 0 && times_very_close( rval + i, rval + i - 1)
                       && !strcmp( rval[i].mpc_code, rval[i - 1].mpc_code))
               j = i-1;
            if( i < n_obs_actually_loaded - 1 && times_very_close( rval + i + 1, rval + i)
                       && !strcmp( rval[i].mpc_code, rval[i + 1].mpc_code))
               j = i+1;
            if( j >= 0)
               {
               memcpy( rval[i].obs_posn, rval[j].obs_posn, 3 * sizeof( double));
               memcpy( rval[i].obs_vel, rval[j].obs_vel, 3 * sizeof( double));
               }
            }
      reset_object_type( rval, n_obs_actually_loaded);
      apply_excluded_observations_file( rval, n_obs_actually_loaded);
      }
   snprintf_err( buff, sizeof( buff), "%d observations actually loaded", n_obs_actually_loaded);
   move_add_nstr( 1, 2, buff, -1);
   refresh_console( );
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
      case OBJECT_INFO_COMPARE_LAST_UPDATED:
         rval = (aptr->jd_updated > bptr->jd_updated ? 1 : -1);
         break;
      case OBJECT_INFO_COMPARE_NAME:
      default:
         rval = strcmp( aptr->obj_name, bptr->obj_name);
         break;
      }
   return( rval);
}

/* My code to download NEOCP observations (see my 'miscell' repository)
stores the time each observation was first found in the normally
unused bytes 60 to 64.  Byte 60 = '~',  61 = month,  62 = day,  63
= hour,  64 = minute,  in MPC-style 'mutant hex'.  If our current
observation is from NEOCP and byte 60 is a tilde,  we can get the
update time.         */

static double get_neocp_update_jd( const char *buff)
{
   double jd = 0.;

   if( !memcmp( buff + 72, "NEOCP", 5) && buff[59] == '~')
      {
      const int month = mutant_hex_char_to_int( buff[60]);
      const int day = mutant_hex_char_to_int( buff[61]);
      const int hour = mutant_hex_char_to_int( buff[62]);
      const int minute = mutant_hex_char_to_int( buff[63]);
      long year = atol( buff + 15);
      const int obs_month = atoi( buff + 18);

      if( obs_month > month)    /* observed last year, posted this one */
         year++;
      jd = (double)dmy_to_day( day, month, year, 0) - 0.5
                              + (double)hour / hours_per_day
                              + (double)minute / minutes_per_day;
      }
   return( jd);
}

void sort_object_info( OBJECT_INFO *ids, const int n_ids,
                                          int object_info_compare_method)
{
   shellsort_r( ids, n_ids, sizeof( OBJECT_INFO), id_compare,
                                    &object_info_compare_method);
}

const char *desig_pattern = NULL, *fullname_pattern = NULL;

/* One can specify,  e.g.,  '-N C314159,K22Ea4C' to process only
those two objects.  I may eventually implement something more
regular expression-ish. */

static bool desig_pattern_matched( OBJECT_INFO *id, const char *pattern)
{
   char packed_desig[13];
   size_t i = 0, j;

   strlcpy_error( packed_desig, id->packed_desig);
   text_search_and_replace( packed_desig, " ", "");
   while( pattern[i])
      {
      char tpattern[30];

      j = i;
      while( pattern[j] != ',' && pattern[j])
         j++;
      assert( j - i < sizeof( tpattern) - 1);
      strlcpy( tpattern, pattern + i, j - i + 1);
      if( pattern_match( tpattern, packed_desig))
         return( true);
      i = j;
      if( pattern[i] == ',')
         i++;
      }
   return( false);
}

static bool fullname_pattern_matched( const char *full_name, const char *pattern)
{
   size_t i = 0, j;

   while( pattern[i])
      {
      char tpattern[80];

      j = i;
      while( pattern[j] != ',' && pattern[j])
         j++;
      assert( j - i < sizeof( tpattern) - 1);
      strlcpy( tpattern, pattern + i, j - i + 1);
      if( pattern_match( tpattern, full_name))
         return( true);
      i = j;
      if( pattern[i] == ',')
         i++;
      }
   return( false);
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


OBJECT_INFO *find_objects_in_file( const char *filename,
                                         int *n_found, const char *station)
{
   static void *obj_name_stack;
   FILE *ifile = (filename ? fopen( filename, "rb") : NULL);
   char new_xdesig[80], new_name[90];
   OBJECT_INFO *rval;
   int i, n = 0, n_alloced = 20, prev_loc = -1;
   const int fixing_trailing_and_leading_spaces =
               *get_environment_ptr( "FIX_OBSERVATIONS");
   char buff[550], mpc_code_from_neocp[4], desig_from_neocp[15];
   void *ades_context;
   const clock_t t0 = clock( );
   int next_output = 2000, n_obs_read = 0;
   long filesize;

   if( obj_name_stack)
      {
      destroy_stack( obj_name_stack);
      obj_name_stack = NULL;
      }

   if( !ifile)
      {
      if( filename)
         debug_printf( "find_objects_in_file: error opening %s: %s\n",
                 filename, strerror( errno));
      if( n_found)
         *n_found = -1;
      return( NULL);
      }
   fseek( ifile, 0L, SEEK_END);
   filesize = ftell( ifile);
   fseek( ifile, 0L, SEEK_SET);
   *mpc_code_from_neocp = '\0';
   *desig_from_neocp = '\0';
   *new_xdesig = *new_name = '\0';
   strcpy( mpc_code_from_neocp, "500");   /* default is geocenter */
   neocp_file_type = NEOCP_FILE_TYPE_UNKNOWN;
   rval = (OBJECT_INFO *)calloc( n_alloced + 1, sizeof( OBJECT_INFO));
   obj_name_stack = create_stack( 2000);
   if( debug_level > 8)
      debug_printf( "About to read input\n");
   ades_context = init_ades2mpc( );
   while( fgets_with_ades_xlation( buff, sizeof( buff), ades_context, ifile))
      {
      size_t iline_len = strlen( buff);
      bool is_neocp = false;
      double jd = 0.;

      if( debug_level > 8)
         debug_printf( "Input line len %d\n", (int)strlen( buff));
      if( *buff == '<')
         remove_html_tags( buff);
      convert_com_to_pound_sign( buff);
      if( !strcmp( buff, "#Combine all"))
         combine_all_observations = "";
      if( !n || *mpc_code_from_neocp)
         is_neocp = get_neocp_data( buff, desig_from_neocp,
                                                 mpc_code_from_neocp);
      if( debug_level > 8)
         debug_printf( "After get_neocp_data\n");
      if( iline_len > MINIMUM_RWO_LENGTH)
         rwo_to_mpc( buff, NULL, NULL, NULL, NULL, NULL);
      if( fixing_trailing_and_leading_spaces)
         fix_up_mpc_observation( buff, &jd);
      if( debug_level > 8)
         debug_printf( "After fixup: %d\n", (int)strlen( buff));
      if( !jd)
         jd = observation_jd( buff);
      if( jd && *new_xdesig)   /* previous line was "COM = (xdesig)";   */
         {                    /* add a new cross-designation to the table */
         memcpy( new_xdesig, buff, 12);
         new_xdesig[12] = ' ';
         for( i = (int)strlen( new_xdesig); i < 26; i++)
            new_xdesig[i] = ' ';
         strcpy( new_xdesig + 26, new_xdesig_indicator);
         xref_designation( new_xdesig);
         *new_xdesig = '\0';
         }
      if( jd && *new_name)    /* 'odd name' sort of xdesig */
         {
         memcpy( new_name, buff, 12);
         new_name[12] = ' ';
         get_object_name( new_name, NULL);
         *new_name = '\0';
         }
      if( is_in_range( jd) && !is_second_line( buff))
         if( !station || !memcmp( buff + 76, station, 3))
            {
            int loc;
            char *tptr;

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
               memset( rval + loc, 0, sizeof( OBJECT_INFO));
               memcpy( rval[loc].packed_desig, buff, 12);
               rval[loc].jd_start = rval[loc].jd_end = jd;
               rval[loc].jd_updated = jd;    /* at minimum */
               if( !is_neocp)   /* for NEOCP obs,  we need to start at the */
                  {                             /* beginning of the file  */
                  rval[loc].file_offset = ftell( ifile) - (long)iline_len
                                       - 100;
                  if( (long)rval[loc].file_offset < 0)
                     rval[loc].file_offset = 0;
                  }
               n++;
               }
            rval[loc].n_obs++;
            if( buff[14] == 'x' || buff[14] == 'X')   /* deleted observation */
               rval[loc].solution_exists++;     /* ..flagged via soln exists */
            if( rval[loc].jd_start > jd)
               rval[loc].jd_start = jd;
            if( rval[loc].jd_end < jd)
               rval[loc].jd_end = jd;
            if( rval[loc].jd_updated < jd)
               rval[loc].jd_updated = jd;
            i = 0;
            tptr = rval[loc].mpc_codes;
            while( tptr[i] && memcmp( tptr + i, buff + 77, 3))
               i += 3;
            if( i < 15)     /* new obscode for this object */
               memcpy( tptr + i, buff + 77, 3);
            if( !memcmp( buff + 72, "NEOCP", 5))
               {
               const double jd_u = get_neocp_update_jd( buff);

               if( rval[loc].jd_updated < jd_u)
                  rval[loc].jd_updated = jd_u;
               }
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
            n_obs_read++;
            if( n_obs_read == next_output)
               {
               char msg_buff[80];
               const double t_elapsed =
                        (double)( clock( ) - t0) / (double)CLOCKS_PER_SEC;
               const double fraction_file_read =
                        (double)ftell( ifile) / (double)filesize;
               const double t_total = t_elapsed / (fraction_file_read + .01);

               next_output += (int)((double)n_obs_read / (t_elapsed + 1.)) / 3;
               snprintf_err( msg_buff, sizeof( msg_buff),
                       "%4.1f%% complete; %.0f seconds elapsed, %.0f remain",
                                    fraction_file_read * 100.,
                                    t_elapsed, t_total - t_elapsed);
               move_add_nstr( 3, 3, msg_buff, -1);
               snprintf_err( msg_buff, sizeof( msg_buff),
                        "%d observations of %d objects read thus far",
                        n_obs_read, n);
               move_add_nstr( 4, 3, msg_buff, -1);
               refresh_console( );
               }
            }
      if( *buff == '#')
         {
         i = 1;            /* check for CSS-style artsat cross-desig,  of */
         while( isdigit( buff[i]))        /* form COM NORAD = Int'l desig */
            i++;
         if( i > 1 && i < 7 && !memcmp( buff + i, "U = ", 4))
            {
            int j = i + 4;

            while( buff[j] > ' ' && buff[j] != ',' && buff[j] != ';')
               j++;
            buff[j] = '\0';
            memmove( buff + 1, buff + i + 2, strlen( buff + i + 1));
            }
         }
      if( !memcmp( buff, "#= ", 3))
         {
         *new_xdesig = '!';
         strlcpy( new_xdesig + 13, buff + 3, 20);
         strlcat( new_xdesig + 13, "            ", 14);
         check_packed_desig_alignment( new_xdesig + 13);
         }
      if( !memcmp( buff, "#fullname ", 10))
         {
         *new_name = '!';
         strlcpy_err( new_name + 13, buff + 10, sizeof( new_name) - 13);
         }
      if( !strcmp( buff, "#ignore obs"))
         while( fgets_trimmed( buff, sizeof( buff), ifile)
                    && !strstr( buff, "end ignore obs"))
            ;     /* deliberately empty loop */
      }
   fclose( ifile);
   free_ades2mpc_context( ades_context);
   *n_found = n;
               /* The allocated hash table is,  at most,  80% full,  with */
               /* plenty of empty entries.  Sorting the entries will put  */
               /* the blank entries at the end of the table,  and the     */
               /* 'n' entries that were actually found at the start.      */
   for( i = n = 0; i < n_alloced; i++)
      if( rval[i].packed_desig[0])
         rval[n++] = rval[i];
   assert( n == *n_found);
   if( desig_pattern)
      {
      i = 0;
      while( i < n)
         if( !desig_pattern_matched( rval + i, desig_pattern))
            {
            memmove( rval + i, rval + i + 1, (n - i - 1) * sizeof( OBJECT_INFO));
            n--;
            }
         else
            i++;
      }
   rval = (OBJECT_INFO *)realloc( rval, n * sizeof( OBJECT_INFO));
           /* Only now do we set the 'full' object names,  because  */
           /* some may have been added mid-file using COM fullname. */
   for( i = 0; i < n; i++)
      {
      get_object_name( buff, rval[i].packed_desig);
      if( fullname_pattern && !fullname_pattern_matched( buff, fullname_pattern))
         {
         memmove( rval + i, rval + i + 1, (n - i - 1) * sizeof( OBJECT_INFO));
         n--;
         i--;
         }
      else
         {
         rval[i].obj_name = (char *)stack_alloc( obj_name_stack, strlen( buff) + 1);
         strcpy( rval[i].obj_name, buff);
         }
      }
   *n_found = n;
   sort_object_info( rval, n, OBJECT_INFO_COMPARE_PACKED);
   return( rval);
}

/* put_observer_data_in_text( ) takes a 'station_no' and fills 'buff'
   with a little bit of text about that station,  as found from
   ObsCodes.htm:  bits such as the lat/lon and name of the station.
   In Find_Orb, these details are shown for the station that made
   the currently-selected observation. */

int put_observer_data_in_text( const char FAR *mpc_code, char *buff)
{
   double lon, lat, alt_in_meters;
   const int planet_idx = get_observer_data_latlon( mpc_code, buff,
                             &lon, &lat, &alt_in_meters);
   const size_t buffsize = 100;

   if( planet_idx == -1)
      {
      char tbuff[4];

      FMEMCPY( tbuff, mpc_code, 4);
      snprintf_err( buff, buffsize, "No information about station '%s'", tbuff);
      }
   else
      {
      const char *name = mpc_station_name( buff);

      memmove( buff, name, strlen( name) + 1);
      if( lon || lat)
         {
         lon *= 180. / PI;
         lat *= 180. / PI;
         if( memcmp( name, "User-supplied", 13))
            {
            const char *output_format = "  (%c%.6f %c%.6f)";

            snprintf_append( buff, buffsize, output_format,
                           (lat > 0. ? 'N' : 'S'), fabs( lat),
                           (lon > 0. ? 'E' : 'W'), fabs( lon));
            }
         if( planet_idx == 3)
            {
            FILE *ifile = fopen_ext( "geo_rect.txt", "fcrb");

            if( ifile)
               {
               extract_region_data_for_lat_lon( ifile,
                                       buff + strlen( buff), lat, lon);
               fclose( ifile);
               }
            }
         }
      }
   return( planet_idx);
}

static const char *environ_dot_dat = "environ.dat";
static char **edata = NULL;
static size_t n_lines = 0, n_lines_allocated = 0;
static bool is_default_environment = false;

int write_environment_pointers( void)
{
   FILE *ofile = fopen_ext( "env.txt", "fcw");
   size_t i;

   fprintf( ofile, "%d vars\n", (int)n_lines);
   for( i = 0; i < n_lines; i++)
      fprintf( ofile, "%s\n", edata[i]);
   fclose( ofile);
   return( (int)n_lines);
}

static size_t get_environment_ptr_index( const char *env_ptr, bool *got_it)
{
   size_t i = 0, n = n_lines;
   const size_t len = strlen( env_ptr);

   *got_it = false;
   while( !*got_it && n)
      {
      size_t j = 0, mid = i + n / 2;

      assert( edata[mid]);
      j = 0;
      while( j < len && env_ptr[j] == edata[mid][j])
         j++;
      if( j == len && edata[mid][j] == '=')
         {
         *got_it = true;
         i = mid;
         }
      else if( j < len && edata[mid][j] > env_ptr[j])
         n /= 2;
      else
         {
         n -= n / 2 + 1;
         i = mid + 1;
         }
      }
   return( i);
}

const char *get_environment_ptr( const char *env_ptr)
{
   size_t i;
   bool got_it;

   if( !env_ptr)
      {
      if( edata)
         {
         for( i = 0; i < n_lines; i++)
            free( edata[i]);
         free( edata);
         }
      edata = NULL;
      return( NULL);
      }
   if( !edata)
      load_default_environment_file( );
   i = get_environment_ptr_index( env_ptr, &got_it);
   if( !got_it)
      return( "");
   else
      return( edata[i] + strlen( env_ptr) + 1);
}

void set_environment_ptr( const char *env_ptr, const char *new_value)
{
   bool got_it;
   const size_t idx = get_environment_ptr_index( env_ptr, &got_it);

   if( !got_it && n_lines == n_lines_allocated)       /* need to expand array */
      {
      n_lines_allocated *= 2;
      if( !n_lines_allocated)
         n_lines_allocated = 8;
      edata = (char **)realloc( edata, n_lines_allocated * sizeof( char *));
      }
   if( !got_it)
      {
      memmove( edata + idx + 1, edata + idx, (n_lines - idx) * sizeof( edata[0]));
      n_lines++;
      edata[idx] = NULL;
      }
   edata[idx] = (char *)realloc( edata[idx],
                        strlen( env_ptr) + strlen( new_value) + 2);
   strcpy( edata[idx], env_ptr);
   strcat( edata[idx], "=");
   strcat( edata[idx], new_value);
}

static int load_json_environment_file( const char *buff)
{
   char key[300];
   const char *tptr;
   size_t depth = 0, starts[30];

   *key = '\0';
   starts[0] = 0;
   for( tptr = buff; *tptr; tptr++)
      if( *tptr == '{')
         {
         depth++;
         starts[depth] = strlen( key);
         }
      else if( *tptr == '}')
         {
         assert( depth);
         depth--;
         }
      else if( *tptr == '"')
         {
         const char *tptr2;
         char *kptr;

         tptr++;
         tptr2 = strchr( tptr, '"');
         assert( tptr2);
         kptr = key + starts[depth];
         if( kptr != key)
            *kptr++ = '_';
         memcpy( kptr, tptr, tptr2 - tptr);
         kptr[tptr2 - tptr] = '\0';
         tptr = tptr2 + 2;
         while( *tptr <= ' ' || *tptr == ':')
            tptr++;
         if( *tptr != '{')     /* we're setting a parameter */
            {
            const char *param = tptr;
            char value[100];

            if( *tptr == '"')
               {
               tptr++;
               param++;
               while( *tptr && *tptr != '"')
                  tptr++;
               }
            else
               while( *tptr > ' ' && *tptr != ',')
                  tptr++;
            memcpy( value, param, tptr - param);
            value[tptr - param] = '\0';
            set_environment_ptr( key, value);
            key[starts[depth]] = '\0';
            }
         else
            tptr--;
         }
   return( 0);
}

int load_environment_file( const char *filename)
{
   FILE *ifile = fopen_ext( filename, (is_default_environment ? "crb" : "rb"));
   char buff[300], *tptr;

   if( !ifile)
      return( -1);
   if( edata)
      is_default_environment = false;
   if( fgets_trimmed( buff, sizeof( buff), ifile) && *buff == '{')
      {
      size_t len, n_read;

      fseek( ifile, 0L, SEEK_END);
      len = (size_t)ftell( ifile);
      fseek( ifile, 0L, SEEK_SET);
      tptr = (char *)malloc( len + 1);
      n_read = fread( tptr, 1, (int)len, ifile);
      fclose( ifile);
      assert( n_read == len);
      if( n_read != len)
         exit( -2);
      tptr[len] = '\0';
      is_default_environment = false;
      load_json_environment_file( tptr);
      free( tptr);
      return( 0);
      }
   fseek( ifile, 0L, SEEK_SET);
   while( fgets_trimmed( buff, sizeof( buff), ifile))
      if( *buff != ' ' && (tptr = strchr( buff, '=')) != NULL)
         {
         *tptr = '\0';
         set_environment_ptr( buff, tptr + 1);
         }
   fclose( ifile);
   return( 0);
}

static int load_default_environment_file( void)
{
   int rval;

   is_default_environment = true;
   rval = load_environment_file( environ_dot_dat);
   if( rval)
      rval = load_environment_file( "environ.def");
   assert( !rval);
   return( rval);
}

void update_environ_dot_dat( void)
{
   if( is_default_environment)
      {
      size_t i, j;
      char **text = load_file_into_memory( environ_dot_dat, NULL, false);
      FILE *ofile;
      char *found = (char *)calloc( n_lines, sizeof( char));

      if( !text)
         text = load_file_into_memory( "environ.def", NULL, true);
      ofile = fopen_ext( environ_dot_dat, "fcwb");
      assert( ofile);
      for( i = 0; text[i]; i++)
         {
         bool updated = false;
         char *tptr = strchr( text[i], '=');

         if( text[i][0] != ' ' && tptr)
            for( j = 0; j < n_lines && !updated; j++)
               if( !strncmp( text[i], edata[j], tptr - text[i] + 1))
                  {
                  fprintf( ofile, "%s\n", edata[j]);
                  found[j] = 1;
                  updated = true;
                  }
         if( !updated)
            fprintf( ofile, "%s\n", text[i]);
         }
      for( i = 0; i < n_lines; i++)
         if( !found[i])
            fprintf( ofile, "%s\n", edata[i]);
      fclose( ofile);
      free( text);
      free( found);
      }
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
   length = hypot( matrix[0], matrix[1]);
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

   https://www.minorplanetcenter.net/iau/info/References.html

   The following function converts those five bytes to human-readable form.

   MPC references are stored as five-digit numbers.  The dreaded MPC 100K
bug was triggered on 2016 May 21,  when MPCs 99895-100318 were issued.  MPC
references past 99999 are now stored as '@' plus four digits.  Note that this
didn't buy us much;  we hit the MPC 110K bug sometime in mid-2018.

   MPC references 110000 and up are stored as '#' plus four "mutant hex"
(base 62) characters,  allowing us an additional 62^4 references.  At the
current rate of using about 5000 MPC references a year,  this should last
about 2900 years.

   MPS references are stored as 'a...z' plus four digits,  allowing 260K
references.  To get beyond that,  MPC used a tilde (~) followed by four
"mutant hex" digits (base 62),  to get us 260000+62^4 = 15 036 336
references.

   MPEC references lack a year;  that's why they are often shown in the form
'MPEC ????-whatever'.  But if the observation was made in the last year, we
_can_ replace those four question marks with the correct year.  (In some
cases,  almost two years;  for example,  if today is 2012 Jan 19,  and the
MPEC is ????-C13,  and the observation was made on 2010 Feb 25,  then the
year must be 2011.  It couldn't be 2010,  because late February would
result in a D or later half-month specification.  And it couldn't be 2012,
because it would then have a B or A half-month specification.)   That's why
the observation JD is passed in. */

static void reference_to_text( char *obuff, const size_t obuff_size,
                                 const char *reference, const double jd)
{
   assert( 5 == strlen( reference));
   if( !strcmp( reference, "     "))       /* no reference given */
      *obuff = '\0';
   else if( !strcmp( reference, "neocp"))
      strcpy( obuff, "NEOCP");
   else if( *reference >= '0' && *reference <= '9')
      snprintf_err( obuff, obuff_size, "MPC %s", reference);
   else if( *reference == '@')
      snprintf_err( obuff, obuff_size, "MPC 10%s", reference + 1);
   else if( *reference >= 'a' && *reference <= 'z')
      snprintf_err( obuff, obuff_size, "MPS %d%s", *reference - 'a', reference + 1);
   else if( *reference == 'E')
      {
      int obs_year, curr_year;
      char obs_letter, curr_letter;

      snprintf_err( obuff, obuff_size, "MPEC ?  ?-%c%d", reference[1], atoi( reference + 2));
      obuff[6] = obuff[7] = '?';    /* attempt to evade trigraph oddities */
      obs_year = get_year_and_mpc_half_month_letter( jd, &obs_letter);
      curr_year = get_year_and_mpc_half_month_letter( current_jd( ), &curr_letter);
      if( curr_letter < reference[1])  /* observation is from last year or earlier */
         curr_year--;
      if( obs_letter > reference[1])
         obs_year++;
      if( curr_year == obs_year)    /* this reference can only mean one year: */
         {
         snprintf_err( obuff + 5, obuff_size - 5, "%4d", curr_year);
         obuff[9] = '-';
         }
      }
   else if( *reference == 'D' && isdigit( reference[1]))
      snprintf_err( obuff, obuff_size, "DASO %d", atoi( reference + 1));
   else if( *reference == '~' || *reference == '#')   /* MPS or MPC number, */
      {                 /* packed as four "mutant hex" (base 62) digits */
      const unsigned number = get_mutant_hex_value( reference + 1, 4);

      snprintf_err( obuff, obuff_size, "MP%c %u", ((*reference == '~') ? 'S' : 'C'),
                        number + ((*reference == '~') ? 260000 : 110000));
      }
   else           /* just copy it in,  but add a space */
      {
      while( *reference >= 'A')
         *obuff++ = *reference++;
      *obuff++ = ' ';
      while( *reference == '0')
         reference++;
      strlcpy_err( obuff, reference, obuff_size);
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
      if( !m->ra_motion && !m->dec_motion)
         m->position_angle_of_motion = 0.;
      else
         m->position_angle_of_motion =
                  180. + (180. / PI) * atan2( -m->ra_motion, -m->dec_motion);
      }
   m->total_motion = hypot( m->ra_motion, m->dec_motion);
   m->xresid = (obs->ra - obs->computed_ra) * cos( obs->dec);
   m->yresid = obs->dec - obs->computed_dec;
               /* cvt xresid, yresid from radians to arcseconds: */
   m->xresid *= (180. / PI) * 3600.;
   m->yresid *= (180. / PI) * 3600.;
               /* time residual is in seconds */
   m->time_residual = m->xresid * m->ra_motion + m->yresid * m->dec_motion;
   m->cross_residual = m->xresid * m->dec_motion - m->yresid * m->ra_motion;
   if( m->total_motion)
      {
      m->cross_residual /= m->total_motion;
      m->time_residual *= 60. / (m->total_motion * m->total_motion);
      }
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
   const size_t obuff_size = 12;

   if( fabs_motion < 99.)
      snprintf_err( obuff, obuff_size, "%5.2f'/hr", motion);
   else if( fabs_motion < 999.)
      snprintf_err( obuff, obuff_size, "%5.1f'/hr", motion);
   else if( fabs_motion < 99999.)
      snprintf_err( obuff, obuff_size, "%5.0f'/hr", motion);
   else if( fabs_motion < 99999. * 60.)
      snprintf_err( obuff, obuff_size, "%5.0f%c/hr", motion / 60., degree_symbol);
   else if( fabs_motion < 99999. * 3600.)
      snprintf_err( obuff, obuff_size, "%5.0f%c/min", motion / 3600., degree_symbol);
   else if( fabs_motion < 99999. * 216000.)
      snprintf_err( obuff, obuff_size, "%5.0f%c/sec", motion / 216000., degree_symbol);
   else
      strlcpy_err( obuff, "!!!!!", obuff_size);
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
   https://www.minorplanetcenter.net/iau/info/RadarObs.html */

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
   mpc_code_t cinfo;
   double xyz[3], vel[3];
   double v1, v2, doppler_factor;
   char tbuff[20];
   char *first_line = obs->second_line + 81;
   int iter;

   assert( obs);
   assert( obs->second_line);
   memcpy( tbuff, obs->second_line + 68, 3);
   tbuff[3] = '\0';
   get_observer_data( tbuff, NULL, &cinfo);
   for( iter = 0; iter < 3; iter++)
      {
      int i;
      double delta[3];

      compute_observer_loc( jd - time_diff / seconds_per_day, 3,
               cinfo.rho_cos_phi, cinfo.rho_sin_phi, cinfo.lon, xyz);
      for( i = 0; i < 3; i++)
         delta[i] = xyz[i] - obs->obj_posn[i];
      time_diff = vector3_length( delta) * AU_IN_KM / SPEED_OF_LIGHT;
      }
   compute_observer_vel( jd - time_diff / seconds_per_day, 3,
               cinfo.rho_cos_phi, cinfo.rho_sin_phi, cinfo.lon, vel);
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
   snprintf_err( buff, 90, "RTDist (C) %.8fs = %.3f km; Dopp %.8f km/s = %.2f Hz",
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

/* get_net_used_from_obs_header( ) looks through the observation header
for the given MPC code for a line starting with 'NET '.  It then looks
through the above 'net_codes' array in hopes of finding a match to the
rest of the line. If it does, we've got our byte for column 72 of the
punched-card astrometry format. */

static int set_data_from_obs_header( OBSERVE *obs)
{
   const char **lines;
   int rval = 0;
   size_t i;

   if( obs->astrometric_net_code != ' ')
      rval = 1;
   if( obs->mag_band != ' ' || obs->obs_mag == BLANK_MAG)
      rval |= 2;
   if( obs_details && (lines = get_code_details( obs_details, obs->mpc_code)) != NULL)
      for( i = 0; lines[i] && rval != 3; i++)
         {
         if( !(rval & 1) && !memcmp( lines[i], "NET ", 4))
            {
            obs->astrometric_net_code = net_name_to_byte_code( lines[i] + 4);
            rval |= 1;
            }
         if( !(rval & 2) && !memcmp( lines[i], "BND ", 4))
            {
            obs->mag_band = lines[i][4];
            rval |= 2;
            }
         }
   return( rval);
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

int show_alt_info = 0;

static int generate_observation_text( const OBSERVE FAR *obs, const int idx,
                 const int n_obs, const int line_number, char *const buff,
                 size_t buffsize)
{
   const OBSERVE FAR *optr = obs + idx;
   const double earth_sun = vector3_length( optr->obs_posn);

   *buff = '\0';
   if( buffsize > 100)      /* no line should be larger than this, */
      buffsize = 100;       /* even if the buffer does have room for it */
   switch( line_number)
      {
      case 0:
         if( optr->note2 == 'R')
            {
            RADAR_INFO rinfo;

            compute_radar_info( optr, &rinfo);
            if( rinfo.rtt_obs)
               {
               snprintf_err( buff, buffsize, "Time (obs) %.7f", rinfo.rtt_obs);
               strip_trailing_zeroes( buff);
               snprintf_append( buff, buffsize, " +/- %.1fus", rinfo.rtt_sigma * 1e+6);
               strip_trailing_zeroes( buff);
               }
            if( rinfo.doppler_obs)
               {
               if( rinfo.rtt_obs)
                  strlcat_err( buff, "   ", buffsize);
               snprintf_append( buff, buffsize, "Shift(obs) %f", rinfo.doppler_obs);
               strip_trailing_zeroes( buff);
               snprintf_append( buff, buffsize, " +/- %f", rinfo.doppler_sigma);
               strip_trailing_zeroes( buff);
               strlcat_err( buff, " Hz", buffsize);
               }
            else
               {
               snprintf_append( buff, buffsize, "  Dist (comp) %.9f = %.2f km", optr->r,
                           optr->r * AU_IN_KM);
               }
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
            snprintf_err( buff, buffsize, "Elong %5.1f    Phase %5.1f    ",
                                        acose( cos_elong) * 180. / PI,
                                        acose( cos_phase) * 180. / PI);
            format_motion( ra_motion_buff, m.ra_motion);
            format_motion( dec_motion_buff, m.dec_motion);
            if( optr->flags & OBS_NO_VELOCITY)
               strlcat_err( buff, "No velocity information", buffsize);
            else
               snprintf_append( buff, buffsize, "RA vel %s   decvel %s   dT=",
                                           ra_motion_buff, dec_motion_buff);

            if( show_alt_info)
               {
               double sig1, sig2, tilt;

               compute_error_ellipse_adjusted_for_motion( &sig1, &sig2, &tilt,
                              optr, &m);
               snprintf_append( buff, buffsize, " %.3fx%.3f PA %.2f\n",
                         sig1, sig2, tilt * 180. / PI);
               }
            else if( !(optr->flags & OBS_NO_VELOCITY))
               {
               if( fabs( m.time_residual) < .999)
                  {
                  char *tptr = buff + strlen( buff);

                  snprintf_append( buff, buffsize, "%.3f sec", fabs( m.time_residual));
                  *tptr = (m.time_residual > 0. ? '+' : '-');
                  }
               else if( fabs( m.time_residual) < 99.9)
                  snprintf_append( buff, buffsize, "%.2f sec", m.time_residual);
               else if( fabs( m.time_residual / 60.) < 99.9)
                  snprintf_append( buff, buffsize, "%.2f min", m.time_residual / 60.);
               else if( fabs( m.time_residual / 60.) < 9999.)
                  snprintf_append( buff, buffsize, "%d min", (int)( m.time_residual / 60.));
               else if( fabs( m.time_residual / 3600.) < 9999.)
                  snprintf_append( buff, buffsize, "%d hr", (int)( m.time_residual / 3600.));
               else
                  strlcat_err( buff, "!!!!", buffsize);
               }
            }
         break;
      case 1:
         if( optr->note2 == 'R')
            show_radar_info( buff, optr);
         else if( optr->flags & OBS_NO_VELOCITY)
            strlcpy_err( buff,
                     "See https://www.projectpluto.com/no_vel.htm for info/fix",
                     buffsize);
         else
            {
            MOTION_DETAILS m;
            char tbuff[15];
            double tdiff;

            compute_observation_motion_details( optr, &m);
            format_motion( tbuff, m.total_motion);
            snprintf_err( buff, buffsize, "ang vel %s at PA %.1f", tbuff,
                      m.position_angle_of_motion);

            snprintf_append( buff, buffsize, "   radial vel %.3f km/s  cross ",
                                      m.radial_vel);
            if( fabs( m.cross_residual) < 9.9)
               snprintf_err( tbuff, sizeof( tbuff), "%.2f", m.cross_residual);
            else if( fabs( m.cross_residual) < 99.9)
               snprintf_err( tbuff, sizeof( tbuff), "%4.1f", m.cross_residual);
            else if( fabs( m.cross_residual) < 9999.)
               snprintf_err( tbuff, sizeof( tbuff), "%4d", (int)m.cross_residual);
            else
               strlcpy_error( tbuff, "!!!!");
            strlcat_err( buff, tbuff, buffsize);
            tdiff = current_jd() - optr->jd;
            if( fabs( tdiff) < 1. / hours_per_day)     /* less than an hour ago */
               snprintf_err( tbuff, sizeof( tbuff), "%d min", (int)( tdiff * minutes_per_day));
            else if( fabs( tdiff) < 1.)
               snprintf_err( tbuff, sizeof( tbuff), "%.1f hr", tdiff * hours_per_day);
            else if( fabs( tdiff) < 100.)
               snprintf_err( tbuff, sizeof( tbuff), "%.1f days", tdiff);
            else
               *tbuff = '\0';
            if( *tbuff)
               snprintf_append( buff, buffsize, "  %s ago", tbuff);
            if( tdiff < 0.)
               strlcat_err( buff, " <FUTURE!>", buffsize);
            }
         break;
      case 2:
         if( show_alt_info && optr->second_line)
            strlcpy_err( buff, optr->second_line, buffsize);
         else if( show_alt_info && optr->ades_ids)
            strlcpy_err( buff, optr->ades_ids, buffsize);
         else
            {
            strlcpy_err( buff, "Delta=", buffsize);
            format_dist_in_buff( buff + strlen( buff), optr->r);  /* ephem0.cpp */

            strlcat_err( buff, "  r=", buffsize);
            format_dist_in_buff( buff + strlen( buff), optr->solar_r);  /* ephem0.cpp */
            strlcat_err( buff, "  ", buffsize);
            if( optr->obs_mag < BLANK_MAG)
               snprintf_append( buff, buffsize, "mag=%5.2f  ", optr->obs_mag);
            else
               strlcat_err( buff, "           ", buffsize);
            if( optr->computed_mag)
               snprintf_append( buff, buffsize, "mag (computed)=%5.2f   ",
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
         const char *net_name = byte_code_to_net_name( optr->astrometric_net_code);
         char *end_ptr;

         if( optr->posn_sigma_1 > 1.01 || optr->posn_sigma_1 < .99
                  || optr->posn_sigma_1 != optr->posn_sigma_2)
            if( optr->note2 != 'R')
               {
               int tilt_angle = 0;
               char sig1_buff[20], sig2_buff[20];

               snprintf_err( sig1_buff, sizeof( sig1_buff), "%.6f", optr->posn_sigma_1);
               remove_insignificant_digits( sig1_buff);
               snprintf_err( sig2_buff, sizeof( sig2_buff), "%.6f", optr->posn_sigma_2);
               remove_insignificant_digits( sig2_buff);
               if( strcmp( sig1_buff, sig2_buff))
                  {
                  strlcat_err( sig1_buff, "x", sizeof( sig1_buff));
                  strlcat_err( sig1_buff, sig2_buff, sizeof( sig1_buff));
                  tilt_angle = (int)( optr->posn_sigma_theta * 180. / PI);
                  }
               snprintf_err( buff, buffsize, "Sigma %s\" ", sig1_buff);
               if( tilt_angle % 180)
                  snprintf_append( buff, buffsize, "%d ", tilt_angle);
               }
         end_ptr = buff + strlen( buff);
         reference_to_text( end_ptr, 15, optr->reference, optr->jd);
         if( *end_ptr)
            strlcat_err( buff, "  ", buffsize);
         if( !get_obs_alt_azzes( optr, &sun_alt_az, &object_alt_az))
            {
            snprintf_append( buff, buffsize, "Obj alt %.1f", object_alt_az.y);
            if( object_alt_az.x > -1.)
               snprintf_append( buff, buffsize, " az %.1f",  object_alt_az.x);
            if( optr->note2 == 'R')
               {
               double xresid, yresid;

               get_residual_data( optr, &xresid, &yresid);
               if( xresid && yresid)
                  snprintf_append( buff, buffsize, "   %.2f,%.2f sigmas",
                           xresid, yresid);
               else
                  snprintf_append( buff, buffsize, "   %.2f sigmas",
                           xresid + yresid);
               }
            else
               {
               snprintf_append( buff, buffsize,
                                 "  Sun alt %.1f", sun_alt_az.y);
               if( sun_alt_az.x > -1.)
                  snprintf_append( buff, buffsize,
                                 " az %.1f", sun_alt_az.x);
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
               snprintf_append( buff, buffsize,
                          " PrimRA %.3f dec %.3f dist %.1f km",
                          ra * 180. / PI,
                          dec * 180. / PI, dist * AU_IN_KM);
               }
#endif
            }
         if( net_name)
            {
            strlcat_err( buff, "  ", buffsize);
            strlcat_err( buff, net_name, buffsize);
            }
         }
         break;
      case 4:
         if( optr->obs_mag < BLANK_MAG)
            snprintf_err( buff, buffsize, "Mag sigma %g; ", optr->mag_sigma);
         else
            *buff = '\0';
         snprintf_append( buff, buffsize, "time sigma %g",
                          optr->time_sigma * seconds_per_day);
         if( optr->unc_time)
            snprintf_append( buff, buffsize, "  uncTime %g",
                          optr->unc_time * seconds_per_day);
         if( optr->ra_bias || optr->dec_bias)
             snprintf_append( buff, buffsize, "  %cRA bias %.3f\" dec bias %.3f\"%c",
                           (apply_debiasing ? ' ' : '['),
                           optr->ra_bias, optr->dec_bias,
                           (apply_debiasing ? ' ' : ']'));
         if( !strcmp( optr->reference, "NEOCP") && optr->columns_57_to_65[3] == '~'
                                       && !show_alt_info)
            {
            const int month = mutant_hex_char_to_int( optr->columns_57_to_65[4]);

            snprintf_append( buff, buffsize, "  TTag:%s %d %02d:%02d",
                     set_month_name( month, NULL),
                     mutant_hex_char_to_int( optr->columns_57_to_65[5]),
                     mutant_hex_char_to_int( optr->columns_57_to_65[6]),
                     mutant_hex_char_to_int( optr->columns_57_to_65[7]));
            }
         if( show_alt_info)
            {
            extern double overobserving_time_span;

            snprintf_append( buff, buffsize, " Nnear=%.3f",
                        n_nearby_obs( obs, n_obs, idx,
                                     overobserving_time_span));
            }
         else if( optr->ades_ids)
            {
            const char *tptr = strstr( optr->ades_ids, "obsID:");

            if( tptr)
               {
               const int ymd = get_mutant_hex_value( tptr + 6, 3);
               const int hms = get_mutant_hex_value( tptr + 9, 3);
               const int millisec = get_mutant_hex_value( tptr + 12, 2);
               const int year = ymd / (31 * 12) + 1800;
               const int month = (ymd / 31) % 12 + 1;
               const int day = ymd % 31 + 1;

               snprintf_append( buff, buffsize,
                           "  TTag %04d-%02d-%02d %02d:%02d:%02d.%03d",
                           year, month, day,
                           hms / 3600, (hms / 60) % 60, hms % 60,
                           millisec);
               }
            }
         break;
     }
   return( *buff ? 0 : -1);      /* indicate if anything was written */
}

static void _add_version_and_de_text( char *buff, const size_t buffsize)
{
   snprintf_err( buff, buffsize, "Version %s\n",
                        find_orb_version_jd( NULL));
   format_jpl_ephemeris_info( buff + strlen( buff));
}

int show_observational_details = 0;

int generate_obs_text( const OBSERVE FAR *obs, const int n_obs, char *buff,
                                               const size_t buffsize)
{
   size_t i, n_selected = 0, first = 0, last = 0;
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
   if( !n_selected)     /* "Click on an observation to..." */
      strlcpy_err( buff, get_find_orb_text( 2025), buffsize);
   else if( n_selected > 1)
      {
      double mean_xresid = 0., mean_yresid = 0., mean_mresid = 0.;
      double mean_xresid2 = 0., mean_yresid2 = 0., mean_mresid2 = 0.;
      double mean_tresid = 0., mean_cross_resid = 0.;
      double mean_tresid2 = 0., mean_cross_resid2 = 0.;
      int n_mags = 0;

      snprintf_err( buff, buffsize, "%d observations selected of %d\n",
                          (int)n_selected, n_obs);
      n_selected = 0;
      for( i = 0; i < (size_t)n_obs; i++)
         if( (obs[i].flags & OBS_IS_SELECTED) && obs[i].note2 != 'R')
            {
            MOTION_DETAILS m;

            compute_observation_motion_details( obs + i, &m);
            mean_xresid += m.xresid;
            mean_yresid += m.yresid;
            mean_tresid += m.time_residual;
            mean_cross_resid += m.cross_residual;
            mean_xresid2 += m.xresid * m.xresid;
            mean_yresid2 += m.yresid * m.yresid;
            mean_tresid2 += m.time_residual * m.time_residual;
            mean_cross_resid2 += m.cross_residual * m.cross_residual;
            n_selected++;
            if( obs[i].obs_mag && obs[i].obs_mag != BLANK_MAG)
               {
               const double mresid = obs[i].obs_mag - obs[i].computed_mag;

               mean_mresid +=   mresid;
               mean_mresid2 +=   mresid *   mresid;
               n_mags++;
               }
            }
      if( n_selected > 1)
         {
         mean_xresid /= (double)n_selected;
         mean_yresid /= (double)n_selected;
         mean_xresid2 /= (double)n_selected;
         mean_yresid2 /= (double)n_selected;
         mean_tresid /= (double)n_selected;
         mean_cross_resid /= (double)n_selected;
         mean_tresid2 /= (double)n_selected;
         mean_cross_resid2 /= (double)n_selected;
         n_lines++;
         snprintf_append( buff, buffsize,
                "Mean RA residual %.3f +/- %.3f; dec %.3f +/- %.3f\n",
                mean_xresid, sqrt( mean_xresid2 - mean_xresid * mean_xresid),
                mean_yresid, sqrt( mean_yresid2 - mean_yresid * mean_yresid));
         }
      if( n_mags > 1)
         {
         mean_mresid /= (double)n_mags;
         mean_mresid2 /= (double)n_mags;
         n_lines++;
         snprintf_append( buff, buffsize,
             "mean mag residual %.2f +/- %.2f\n",
             mean_mresid, sqrt( mean_mresid2 - mean_mresid * mean_mresid));
         }
      if( n_selected)
         snprintf_append( buff, buffsize,
             "Mean time residual %.3f +/- %.3fs; mean cross-track resid %.3f +/- %.3f\"\n",
             mean_tresid, sqrt( mean_tresid2 - mean_tresid * mean_tresid),
             mean_cross_resid, sqrt( mean_cross_resid2 - mean_cross_resid * mean_cross_resid));
      if( n_selected == 2)
         {
         double dist, posn_ang, delta_time;
         const OBSERVE FAR *optr1 = obs + first;
         const OBSERVE FAR *optr2 = obs + last;

         calc_dist_and_posn_ang( &optr1->ra, &optr2->ra, &dist, &posn_ang);
         dist *= 180. / PI;      /* cvt radians to degrees */
         delta_time = optr2->jd - optr1->jd;

         n_lines++;
         snprintf_append( buff, buffsize,
                   "Observations are %.2f\" = %.2f' = %.3f degrees apart\n",
                   dist * 3600., dist * 60., dist);
         n_lines++;
         if( fabs( delta_time) < 1.)
            snprintf_append( buff, buffsize,
                     "Time diff: %.2f sec = %.2f min = %.3f hrs\n",
                     delta_time * seconds_per_day,
                     delta_time * minutes_per_day,
                     delta_time * hours_per_day);
         else
            snprintf_append( buff, buffsize,
                     "Time diff: %.1f hrs = %.2f days\n",
                     delta_time * 24., delta_time);
         dist /= delta_time;     /* get motion in degrees/day */
         dist *= 60. / 24.;      /* then convert to '/hr */
                           /* Dunno how the PA got flipped,  but it did: */
         posn_ang = 2. * PI - posn_ang;
         n_lines++;
         snprintf_append( buff, buffsize,
                  "Motion: %.2f'/hr in RA, %.2f'/hr in dec",
                  dist * sin( posn_ang), dist * cos( posn_ang));
         n_lines++;
         snprintf_append( buff, buffsize,
                  " (total %.2f'/hr at PA %.1f)\n",
                  dist, posn_ang * 180. / PI);
         }
      make_date_range_text( buff + strlen( buff),
                                   obs[first].jd, obs[last].jd);
      strlcat_err( buff, "\n", buffsize);
      n_lines++;
      }
   else if( show_observational_details)
      {
      const char **lines = obs[first].obs_details;

      n_lines = 0;
      *buff = '\0';
      if( lines)
         while( lines[n_lines])
            {
            strlcat_err( buff, lines[n_lines], buffsize);
            strlcat_err( buff, "\n", buffsize);
            n_lines++;
            }
      else
         {
         FILE *ifile = fopen_ext( "details.txt", "fcrb");

         while( fgets( tptr, 100, ifile))
            if( !memcmp( tptr, "COD ", 4) && !memcmp( tptr + 4, obs[first].mpc_code, 3))
               {
               while( fgets( tptr, 100, ifile) && memcmp( tptr, "COD", 3) && n_lines < 8)
                  if( *tptr >= ' ')
                     {
                     tptr += strlen( tptr);
                     n_lines++;
                     }
               }
         *tptr = '\0';
         fclose( ifile);
         }
      if( !n_lines)
         {
         strlcpy_err( buff, get_find_orb_text( 2026), buffsize);
         n_lines = 1;            /* "No obs header" message */
         }
      }
   else        /* "standard",  computed details */
      {
      n_lines = 0;
      *buff = '\0';
      for( i = 0; i < 5; i++)
         {
         const size_t remaining_bytes = buffsize - strlen( buff);

         if( !generate_observation_text( obs, (int)first, n_obs, (int)i,
                           buff + strlen( buff), remaining_bytes))
            {
            strlcat_err( buff, "\n", buffsize);
            n_lines++;
            }
         }
      }
   if( n_lines <= 4)    /* got room for version info */
      {
      _add_version_and_de_text( buff + strlen( buff), buffsize - strlen( buff));
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
                  snprintf_err( tbuff, sizeof( tbuff),
                     "Line %ld: Sun alt %.1f az %.1f; obj alt %.1f az %.1f\n",
                           line_no,
                           alt_az_sun.y, alt_az_sun.x,
                           alt_az_obj.y, alt_az_obj.x);
                  printf( "%s %s", tbuff, buff);
                  fprintf( ofile, "%s %s", tbuff, buff);
                  n_problems_found++;
                  }
            override_time = 0.;
            override_ra = override_dec = override_mag = -100.;
            }
      }
   fclose( ifile);
   if( ofile)
      fclose( ofile);
   return( n_problems_found);
}

