/* elem_out.cpp: formatting elements into human-friendly form

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

#define __STDC_FORMAT_MACROS

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include "watdefs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "mpc_func.h"
#include "date.h"
#include "afuncs.h"
#include "lunar.h"
#include "monte0.h"     /* for put_double_in_buff() proto */
#include "showelem.h"

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>

bool findorb_already_running = false;
#endif

            /* Pretty much every platform I've run into supports */
            /* Unicode display,  except OpenWATCOM and early     */
            /* versions of MSVC.                                 */
#if !defined( __WATCOMC__)
   #if !defined( _MSC_VER) || (_MSC_VER > 1100)
      #define HAVE_UNICODE
   #endif
#endif

#define J2000 2451545.
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define JD_TO_YEAR(jd)  (2000. + ((jd)-J2000) / 365.25)
#define YEAR_TO_JD( year) (J2000 + (year - 2000.) * 365.25)


extern double optical_albedo;

#ifdef _MSC_VER   /* MSVC/C++ lacks snprintf.  See 'ephem0.cpp' for details. */
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif
// void elements_in_tle_format( char *buff, const ELEMENTS *elem);
int snprintf_append( char *string, const size_t max_len,      /* ephem0.cpp */
                                   const char *format, ...)
#ifdef __GNUC__
         __attribute__ (( format( printf, 3, 4)))
#endif
;
#ifdef NOT_CURRENTLY_IN_USE
#define ssnprintf_append( obuff, ...) snprintf_append( obuff, sizeof( obuff), __VA_ARGS__)
#define ssnprintf( obuff, ...) snprintf( obuff, sizeof( obuff), __VA_ARGS__)
#endif
int store_defaults( const ephem_option_t ephemeris_output_options,
         const int element_format, const int element_precision,
         const double max_residual_for_filtering,
         const double noise_in_arcseconds);           /* elem_out.cpp */
int get_defaults( ephem_option_t *ephemeris_output_options, int *element_format,
         int *element_precision, double *max_residual_for_filtering,
         double *noise_in_arcseconds);                /* elem_out.cpp */
static int elements_in_mpcorb_format( char *buff, const char *packed_desig,
                const char *full_desig, const ELEMENTS *elem,
                const OBSERVE FAR *obs, const int n_obs);   /* orb_func.c */
static int elements_in_guide_format( char *buff, const ELEMENTS *elem,
                     const char *obj_name, const OBSERVE *obs,
                     const unsigned n_obs);                /* orb_func.c */
FILE *open_json_file( char *filename, const char *env_ptr, const char *default_name,
                  const char *packed_desig, const char *permits); /* ephem0.cpp */
int find_worst_observation( const OBSERVE FAR *obs, const int n_obs);
double initial_orbit( OBSERVE FAR *obs, int n_obs, double *orbit);
int set_locs( const double *orbit, double t0, OBSERVE FAR *obs, int n_obs);
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
double calc_obs_magnitude( const double obj_sun,
          const double obj_earth, const double earth_sun, double *phase_ang);
int find_best_fit_planet( const double jd, const double *ivect,
                                 double *rel_vect);         /* runge.cpp */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
void remove_trailing_cr_lf( char *buff);      /* ephem0.cpp */
int write_tle_from_vector( char *buff, const double *state_vect,
        const double epoch, const char *norad_desig, const char *intl_desig);
int setup_planet_elem( ELEMENTS *elem, const int planet_idx,
                                          const double t_cen);   /* moid4.c */
void set_environment_ptr( const char *env_ptr, const char *new_value);
double find_collision_time( ELEMENTS *elem, double *latlon, const int is_impact);
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
int get_idx1_and_idx2( const int n_obs, const OBSERVE FAR *obs,
                                int *idx1, int *idx2);      /* elem_out.c */
double mag_band_shift( const char mag_band);                /* elem_out.c */
int get_jpl_ephemeris_info( int *de_version, double *jd_start, double *jd_end);
double *get_asteroid_mass( const int astnum);   /* bc405.cpp */
double current_jd( void);                       /* elem_out.cpp */
double centralize_ang( double ang);             /* elem_out.cpp */
char *get_file_name( char *filename, const char *template_file_name);
void get_relative_vector( const double jd, const double *ivect,
          double *relative_vect, const int planet_orbiting);  /* orb_func.c */
void push_orbit( const double epoch, const double *orbit);  /* orb_fun2.c */
int pop_orbit( double *epoch, double *orbit);               /* orb_fun2.c */
double get_planet_mass( const int planet_idx);                /* orb_func.c */
double observation_rms( const OBSERVE FAR *obs);            /* elem_out.cpp */
double find_epoch_shown( const OBSERVE *obs, const int n_obs); /* elem_out */
double evaluate_initial_orbit( const OBSERVE FAR *obs,      /* orb_func.c */
                              const int n_obs, const double *orbit);
double diameter_from_abs_mag( const double abs_mag,      /* ephem0.cpp */
                                     const double optical_albedo);
char **load_file_into_memory( const char *filename, size_t *n_lines,
                        const bool fail_if_not_found);      /* mpc_obs.cpp */
const char *get_find_orb_text( const int index);      /* elem_out.cpp */
int set_language( const int language);                      /* elem_out.cpp */
void get_find_orb_text_filename( char *filename);     /* elem_out.cpp */
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
static int names_compare( const char *name1, const char *name2);
static int get_uncertainty( const char *key, char *obuff, const bool in_km);
char *iso_time( char *buff, const double jd, const int precision); /* elem_out.cpp */
char *real_packed_desig( char *obuff, const char *packed_id);  /* ephem0.cpp */
extern int debug_level;
double asteroid_magnitude_slope_param = .15;
double comet_magnitude_slope_param = 10.;
char default_comet_magnitude_type = 'N';
const char *mpc_fmt_filename = "mpc_fmt.txt";
const char *sof_filename = "sof.txt";
const char *sofv_filename = "sofv.txt";
int force_model = 0;
extern int forced_central_body;
int get_planet_posn_vel( const double jd, const int planet_no,
                     double *posn, double *vel);         /* runge.cpp */
void compute_variant_orbit( double *variant, const double *ref_orbit,
                     const double n_sigmas);       /* orb_func.cpp */
void make_config_dir_name( char *oname, const char *iname);  /* miscell.cpp */
int put_elements_into_sof( char *obuff, const char *templat,
         const ELEMENTS *elem,
         const int n_obs, const OBSERVE *obs);                /* elem_ou2.cpp */
uint64_t parse_bit_string( const char *istr);                /* miscell.cpp */
const char *write_bit_string( char *ibuff, const uint64_t bits);

int debug_printf( const char *format, ...)                 /* runge.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;

bool is_inverse_square_force_model( void)
{
   return( (force_model & 0x10) != 0);
}

char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile)
{
   char *rval = fgets( buff, (int)max_bytes, ifile);

   if( rval)
      {
      int i;

      for( i = 0; buff[i] && buff[i] != 10 && buff[i] != 13; i++)
         ;
      buff[i] = '\0';
      }
   return( rval);
}

void get_first_and_last_included_obs( const OBSERVE *obs,
              const int n_obs, int *first, int *last)       /* elem_out.c */
{
   if( first)
      for( *first = 0; *first < n_obs - 1 && !obs[*first].is_included;
                                    (*first)++)
         ;
   if( last)
      for( *last = n_obs - 1; *last && !obs[*last].is_included; (*last)--)
         ;
}

void make_date_range_text( char *obuff, const double jd1, const double jd2)
{
   long year, year2;
   int month, month2;
   const int day1 = (int)decimal_day_to_dmy( jd1, &year,  &month,
                                    CALENDAR_JULIAN_GREGORIAN);
   const int day2 = (int)decimal_day_to_dmy( jd2, &year2, &month2,
                                    CALENDAR_JULIAN_GREGORIAN);
   static const char *month_names[] = { "Jan.", "Feb.", "Mar.", "Apr.", "May",
            "June", "July", "Aug.", "Sept.", "Oct.", "Nov.", "Dec." };

   if( year == year2)
      {
      sprintf( obuff, "%ld %s %d", year, month_names[month - 1], day1);
      obuff += strlen( obuff);
      if( month == month2 && day1 != day2)
         sprintf( obuff, "-%d", day2);
      else if( month != month2)
         sprintf( obuff, "-%s %d", month_names[month2 - 1], day2);
      }
   else              /* different years */
      sprintf( obuff, "%ld %s %d-%ld %s %d", year, month_names[month - 1],
                             day1, year2, month_names[month2 - 1], day2);

   obuff += strlen( obuff);
   if( jd2 - jd1 < 10. / seconds_per_day)  /* less than 10 seconds: show to .01 sec */
      sprintf( obuff, " (%.2f sec)", (jd2 - jd1) * seconds_per_day);
   else if( jd2 - jd1 < 100. / seconds_per_day) /* less than 100 seconds: show to .1 sec */
      sprintf( obuff, " (%.1f sec)", (jd2 - jd1) * seconds_per_day);
   else if( jd2 - jd1 < 100. / minutes_per_day)     /* less than 100 minutes: show in min */
      sprintf( obuff, " (%.1f min)", (jd2 - jd1) * minutes_per_day);
   else if( jd2 - jd1 < 2.)
      sprintf( obuff, " (%.1f hr)", (jd2 - jd1) * hours_per_day);
}

/* This is useful for abbreviating on-screen text;  say,  displaying */
/* all planet names chopped down to four characters.  With UTF-8 text, */
/* this may not mean four bytes. */

const char *find_nth_utf8_char( const char *itext, size_t n)
{
   while( *itext && n--)
      {
      switch( ((unsigned char)*itext) >> 4)
         {
         case 0xf:          /* four-byte token;  U+10000 to U+1FFFFF */
            itext += 4;
            break;
         case 0xe:          /* three-byte token; U+0800 to U+FFFF */
            itext += 3;
            break;
         case 0xc:          /* two-byte token: U+0080 to U+03FF */
         case 0xd:          /* two-byte token: U+0400 to U+07FF */
            itext += 2;
            break;
         default:          /* "ordinary" ASCII (U+0 to U+7F) */
            itext++;       /* single-byte token              */
            break;
         }
      }
   return( itext);
}

            /* String file name defaults to English,  but can be replaced */
            /* with ifindorb.txt (Italian), ffindorb.txt (French), etc.   */
char findorb_language = 'e';

void get_find_orb_text_filename( char *filename)
{
   strcpy( filename, "efindorb.txt");
   *filename = findorb_language;
}

         /* The following only works for Win1252,  and even there, */
         /* the part from 0x80 to 0x9f fails.  But we don't have   */
         /* Euro signs and such in Find_Orb at this point.         */
#ifndef HAVE_UNICODE
void utf8_to_win1252( char *text)
{
   char *optr = text;

   while( *text)
      if( (unsigned char)*text < 0x80)
         *optr++ = *text++;
      else
         {
         const unsigned char t0 = (unsigned char)text[0];
         const unsigned char t1 = (unsigned char)text[1];

         *optr++ = (char)( (t0 << 6) | (t1 & 0x3f));
         text += 2;
         }
   *optr = '\0';
}
#endif

/* Multi-line text in ?findorb.txt files is 'joined' by this function.
If two or more consecutive lines have the same index,  we concatenate
them, with a line feed in between them.  */

static void collapse_findorb_txt( char **text)
{
   size_t i, j;

   for( i = 0; text[i]; i++)
      {
      text_search_and_replace( text[i], "\\n", "\n");
      text_search_and_replace( text[i], "[email]", "pluto@p");
      }
   for( i = 0; text[i]; i++)
      {
      const int idx = atoi( text[i]);

      if( text[i + 1] && idx && idx == atoi( text[i + 1]))
         {
         strcat( text[i], "\n");
         memmove( text[i] + strlen( text[i]), text[i + 1] + 8,
                          strlen( text[i + 1] + 7));
         for( j = i + 1; text[j + 1]; j++)
            text[j] = text[j + 1];
         text[j] = NULL;
         i--;        /* may be more lines to concatenate */
         }
      }
}

const char *get_find_orb_text( const int index)
{
   static char **text = NULL, **default_text = NULL;
   size_t i;
   static char currently_loaded_language = '\0';
   char filename[20];

   if( !index)          /* clean up */
      {
      if( text)
         free( text);
      if( default_text)
         free( default_text);
      default_text = text = NULL;
      return( NULL);
      }
   get_find_orb_text_filename( filename);
   if( findorb_language != 'e' &&
               currently_loaded_language != findorb_language)
      {
      if( text)          /* 'text' = pointer to strings for current language */
         free( text);
      text = load_file_into_memory( filename, NULL, true);
      collapse_findorb_txt( text);
      currently_loaded_language = findorb_language;
      }
   if( !default_text)     /* 'default_text' = strings for English,  and as */
      {                   /* a backstop when something isn't translated    */
      *filename = 'e';
      default_text = load_file_into_memory( filename, NULL, true);
      collapse_findorb_txt( default_text);
      }
   if( text && findorb_language != 'e')  /* try non-default language... */
      {
      for( i = 0; text[i]; i++)
         if( atoi( text[i]) == index)
            return( text[i] + 8);
      }              /* ...if that fails,  search strings for English : */
   for( i = 0; default_text[i]; i++)
      if( atoi( default_text[i]) == index)
         return( default_text[i] + 8);
   assert( 1);             /* i.e.,  should never get here */
   return( NULL);
}

/* observation_summary_data( ) produces the final line in an MPC report,
   such as 'From 20 observations 1997 Oct. 20-22;  mean residual 0".257.   '
   Note that the arcsecond mark comes before the decimal point;  this
   oddity is handled using the text_search_and_replace() function.
*/

static void observation_summary_data( char *obuff, const OBSERVE FAR *obs,
                              const int n_obs, const int options)
{
   int i, n_included, first_idx, last_idx;

   get_first_and_last_included_obs( obs, n_obs, &first_idx, &last_idx);
   for( i = n_included = 0; i < n_obs; i++)
      n_included += obs[i].is_included;
   if( options == -1)      /* 'guide.txt' bare-bones format */
      sprintf( obuff, "%d of %d", n_included, n_obs);
   else if( (options & ELEM_OUT_ALTERNATIVE_FORMAT) && n_included != n_obs)
      sprintf( obuff, get_find_orb_text( 15), n_included, n_obs);
   else
      sprintf( obuff, get_find_orb_text( 16), n_included);
   if( options != -1 && n_included)
      {
      const double rms = compute_rms( obs, n_obs);
      char rms_buff[14];
      const char *rms_format = "%.2f";

      strcat( obuff, " ");
      obuff += strlen( obuff);
      make_date_range_text( obuff, obs[first_idx].jd, obs[last_idx].jd);
      obuff += strlen( obuff);
      if( options & ELEM_OUT_PRECISE_MEAN_RESIDS)
         rms_format = (rms > 0.003 ? "%.3f" : "%.1e");
      sprintf( rms_buff, rms_format, rms);
      text_search_and_replace( rms_buff, ".", "\".");
      sprintf( obuff, get_find_orb_text( 17), rms_buff);
      }                                 /* "; mean residual %s." */
}

double centralize_ang( double ang)
{
   ang = fmod( ang, PI + PI);
   if( ang < 0.)
      ang += PI + PI;
   return( ang);
}

void convert_elements( const double epoch_from, const double epoch_to,
      double *incl, double *asc_node, double *arg_per);     /* conv_ele.cpp */

   /* Packed MPC designations have leading and/or trailing spaces.  This */
   /* function lets you get the designation minus those spaces.          */

static void packed_desig_minus_spaces( char *obuff, const char *ibuff)
{
   while( *ibuff && *ibuff == ' ')
      ibuff++;
   while( *ibuff && *ibuff != ' ')
      *obuff++ = *ibuff++;
   *obuff = '\0';
}

double current_jd( void)
{
   static const double jan_1970 = 2440587.5;
   const double jd = jan_1970 + (double)time( NULL) / seconds_per_day;

   return( jd);
}

int n_clones_accepted = 0;

static int elements_in_mpcorb_format( char *buff, const char *packed_desig,
                const char *full_desig, const ELEMENTS *elem,
                const OBSERVE FAR *obs, const int n_obs)   /* orb_func.c */
{
   extern unsigned perturbers;
   int month, day, i, first_idx, last_idx, n_included_obs = 0;
   long year;
   const double rms_err = compute_rms( obs, n_obs);
   const unsigned hex_flags = 0;
            /* 'mpcorb' has four hexadecimal flags starting in column 162, */
            /* signifying if the object is in any of various classes such  */
            /* as Aten,  scattered-disk object,  PHA,  Jupiter Trojan,     */
            /*  etc.  None of those flags are set yet.                     */
   const int n_oppositions = 1;
            /* The above needs some work.  Problem is,  what constitutes   */
            /* an "opposition" for an NEO?  (It's more clearcut for MBOs.) */
            /* For the nonce,  we'll just say "one opposition".            */
   double arc_length;
   char packed_desig2[40];

   packed_desig_minus_spaces( packed_desig2, packed_desig);
   sprintf( buff, "%-8s%5.2f  %4.2f ", packed_desig2, elem->abs_mag,
                           asteroid_magnitude_slope_param);
   day = (int)( decimal_day_to_dmy( elem->epoch, &year,
                              &month, CALENDAR_JULIAN_GREGORIAN) + .0001);
   sprintf( buff + 20, "%c%02ld%X%c",
                  int_to_mutant_hex_char( year / 100),
                  year % 100L, month,
                  int_to_mutant_hex_char( day));
   sprintf( buff + 25, "%10.5f%11.5f%11.5f%11.5f%11.7f",
           centralize_ang( elem->mean_anomaly) * 180. / PI,
           centralize_ang( elem->arg_per) * 180. / PI,
           centralize_ang( elem->asc_node) * 180. / PI,
           centralize_ang( elem->incl) * 180. / PI,
           elem->ecc);
   sprintf( buff + 79, "%12.8f%12.7f",
            (180 / PI) / elem->t0,        /* n */
            elem->major_axis);
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         n_included_obs++;
   day = (int)( decimal_day_to_dmy( current_jd( ),
                         &year, &month, CALENDAR_JULIAN_GREGORIAN) + .0001);
   sprintf( buff + 103,
      "    FO %02d%02d%02d  %4d  %2d ****-**** ****         Find_Orb   %04x",
                  (int)( year % 100), month, (int)day,
                  n_included_obs, n_oppositions, hex_flags);
   get_first_and_last_included_obs( obs, n_obs, &first_idx, &last_idx);
   arc_length = obs[last_idx].jd - obs[first_idx].jd;
   if( arc_length < 99. / seconds_per_day)
      sprintf( buff + 127, "%4.1f sec ", arc_length * seconds_per_day);
   else if( arc_length < 99. / minutes_per_day)
      sprintf( buff + 127, "%4.1f min ", arc_length * minutes_per_day);
   else if( arc_length < 2.)
      sprintf( buff + 127, "%4.1f hrs ", arc_length * hours_per_day);
   else if( arc_length < 600.)
      sprintf( buff + 127, "%4d days", (int)arc_length + 1);
   else
      sprintf( buff + 127, "%4d-%4d",
                (int)JD_TO_YEAR( obs[first_idx].jd),
                (int)JD_TO_YEAR( obs[last_idx].jd));
   buff[136] = ' ';
   sprintf( buff + 165, " %-30s", full_desig);
   day = (int)( decimal_day_to_dmy( obs[last_idx].jd, &year,
                       &month, CALENDAR_JULIAN_GREGORIAN) + .0001);
   sprintf( buff + 194, "%04ld%02d%02d", year, month, day);
   if( rms_err < 9.9)
      sprintf( buff + 137, "%4.2f", rms_err);
   else if( rms_err < 99.9)
      sprintf( buff + 137, "%4.1f", rms_err);
   else if( rms_err < 9999.)
      sprintf( buff + 137, "%4.0f", rms_err);
   buff[141] = ' ';
   if( (perturbers & 0x1fe) == 0x1fe)
      {    /* we have Mercury through Neptune,  at least */
      const char *coarse_perturb, *precise_perturb;

      if( perturbers & 0x700000)    /* asteroids included */
         {
         precise_perturb = (perturbers & 0x400 ? "3E" : "38");
         coarse_perturb = "M-v";
         }
      else        /* non-asteroid case */
         {
         precise_perturb = (perturbers & 0x400 ? "06" : "00");
         coarse_perturb = (perturbers & 0x200 ? "M-P" : "M-N");
         }
      memcpy( buff + 142, coarse_perturb, 3);
      memcpy( buff + 146, precise_perturb, 2);
      }
   return( 0);
}

char *iso_time( char *buff, const double jd, const int precision)
{
   full_ctime( buff, jd, CALENDAR_JULIAN_GREGORIAN
                | FULL_CTIME_YMD | FULL_CTIME_MONTHS_AS_DIGITS
                | FULL_CTIME_LEADING_ZEROES | (precision << 4));
   buff[4] = buff[7] = '-';
   buff[10] = 'T';
   strcat( buff, "Z");
   return( buff);
}

/* This is ludicrously high,  but let's say that we'll never
report a linkage between 50 or more tracklets at once. */

#define MAX_LINKAGE_IDS 50

static int make_linkage_json( const int n_obs, const OBSERVE *obs, const ELEMENTS *elem)
{
   int i, j, n_ids = 0, idx[MAX_LINKAGE_IDS], n_designated = 0;
   FILE *ofile, *ifile;
   char buff[200];

   for( i = 0; i < n_obs && n_ids < MAX_LINKAGE_IDS; i++)
      {
      j = 0;
      while( j < n_ids && strcmp( obs[i].packed_id, obs[idx[j]].packed_id))
         j++;
      if( j == n_ids)
         {
         const int desig_type = unpack_mpc_desig( NULL, obs[i].packed_id);

         if( desig_type != OBJ_DESIG_OTHER && desig_type != OBJ_DESIG_ARTSAT)
            n_designated++;
         idx[n_ids++] = i;
         }
      }
   if( n_ids < 2 || n_ids >= MAX_LINKAGE_IDS)
      if( elem->central_obj != 3)
         return( n_ids);                  /* no ID to be made */
   ifile = fopen_ext( "link_hdr.json", "cr");
   if( !ifile)   /* no user-modified header; fall back to default hdr */
      ifile = fopen_ext( "link_def.json", "fcr");
   ofile = fopen_ext( (elem->central_obj == 3) ? "artsat.json" : "linkage.json", "fcw");
   while( fgets( buff, sizeof( buff), ifile))
      if( *buff != '#')
         {
         char *tptr = strstr( buff, "%t");
         char tbuff[200];

         if( tptr)
            {
            full_ctime( tbuff, current_jd( ), FULL_CTIME_YMD);
            text_search_and_replace( buff, "%t", tbuff);
            }
         if( elem->central_obj == 3)
            text_search_and_replace( buff, "\"comment\": \"",
                     "\"comment\": \"Identified as artsat. ");

         fputs( buff, ofile);
         }
   fclose( ifile);
   if( elem->central_obj == 3)
      fprintf( ofile, "      \"designations\": [\n"
                      "         \"ARTSAT\"\n"
                      "        ],\n");
   if( n_designated)
      {
      fprintf( ofile, "      \"designations\": [\n");
      for( i = j = 0; i < n_ids; i++)
         {
         const int desig_type = unpack_mpc_desig( NULL, obs[idx[i]].packed_id);

         if( desig_type != OBJ_DESIG_OTHER && desig_type != OBJ_DESIG_ARTSAT)
            {
            strcpy( buff, obs[idx[i]].packed_id);
            text_search_and_replace( buff, " ", "");
            j++;
            fprintf( ofile, "        \"%s\"%s\n", buff, (j == n_designated ? "" : ","));
            }
         }
      fprintf( ofile, "      ],\n");
      }
   if( n_ids != n_designated)
      {
      fprintf( ofile, "      \"trksubs\": [\n");
      for( i = j = 0; i < n_ids; i++)
         {
         const int desig_type = unpack_mpc_desig( NULL, obs[idx[i]].packed_id);

         if( desig_type == OBJ_DESIG_OTHER || desig_type == OBJ_DESIG_ARTSAT)
            {
            const double obs_utc_jd = obs[idx[i]].jd
                        - td_minus_utc( obs[idx[i]].jd) / seconds_per_day;

            strcpy( buff, obs[idx[i]].packed_id);
            text_search_and_replace( buff, " ", "");
            j++;
            fprintf( ofile, "        [\n");
            fprintf( ofile, "        \"%s\",\n", buff);
            full_ctime( buff, obs_utc_jd, FULL_CTIME_YMD | FULL_CTIME_LEADING_ZEROES
                           | FULL_CTIME_MONTHS_AS_DIGITS
                           | FULL_CTIME_FORMAT_DAY | FULL_CTIME_5_PLACES);
            fprintf( ofile, "        \"%s\",\n", buff);
            fprintf( ofile, "        \"%s\"\n", obs[idx[i]].mpc_code);
            fprintf( ofile, "        ]%s\n", (j == n_ids - n_designated ? "" : ","));
            }
         }
      fprintf( ofile, "      ],\n");
      }
   i = 0;
   while( i < n_obs && strcmp( obs[i].reference, "NEOCP"))
      i++;
   if( i < n_obs)         /* no NEOCP observations */
   fprintf( ofile, "      \"identification_type\": \"neocp\"%s\n",
                (elem->central_obj == 3) ? "" : ",");
   if( elem->central_obj != 3)
      {
      fprintf( ofile, "      \"orbit\": {\n");
      fprintf( ofile, "        \"arg_pericenter\": %f,\n",
                                    centralize_ang( elem->arg_per) * 180. / PI);
      fprintf( ofile, "        \"eccentricity\": %.8f,\n", elem->ecc);
      fprintf( ofile, "        \"epoch\": %f,\n", elem->epoch);
      fprintf( ofile, "        \"inclination\": %f,\n", elem->incl * 180. / PI);
      fprintf( ofile, "        \"lon_asc_node\": %f,\n",
                                   centralize_ang( elem->asc_node) * 180. / PI);
      fprintf( ofile, "        \"pericenter_distance\": %.15f,\n", elem->q);
      fprintf( ofile, "        \"pericenter_time\": %f\n", elem->perih_time);
      fprintf( ofile, "      }\n");
      }
   fprintf( ofile, "    }\n  }\n}\n");
   fclose( ofile);
   return( n_ids);
}

/* _Usually_,  you can tell what type of object you have from the
length of the packed designation,  and from that,  what the alignment
within the first twelve bytes of a punched-card record should be.
The only exception I know of is cases where somebody uses a five-byte
temporary identifier;  that can't be distinguished from a numbered
minor planet or comet so easily.  It may not really be a problem
for us,  though it's probably an issue for MPC.  */

static inline void align_packed_desig( char *obuff, const char *packed_desig)
{
   int loc;
   const size_t len = strlen( packed_desig);

   assert( len > 0 && len < 13);
   if( len == 8)     /* mostly comet desigs */
      loc = 4;
   else if( len == 5)         /* mostly permanent desigs */
      loc = 0;
   else if( len > 8)       /* artsat,  "overlong" non-standard desigs */
      loc = 0;
   else                    /* short,  usually temporary desigs */
      loc = 5;
   memset( obuff, ' ', 12);
   memcpy( obuff + loc, packed_desig, strlen( packed_desig));
   obuff[12] = '\0';
}

char *find_numbered_mp_info( const int number);             /* mpc_obs.cpp */

static char *object_name( char *buff, const int obj_index)
{
   if( !obj_index)
      strcpy( buff, "Sun");
   else if( obj_index > 0 && obj_index < 11) /* nine planets & moon */
      strcpy( buff, get_find_orb_text( 99107 + obj_index));
   else
      sprintf( buff, "Object_%d\n", obj_index);
   return( buff);
}

static int elements_in_json_format( FILE *ofile, const ELEMENTS *elem,
                     const char *obj_name, const OBSERVE *obs,
                     const unsigned n_obs, const double *moids,
                     const char *body_frame_note, const bool show_obs)
{
   extern int n_extra_params;
   double jd = current_jd( );
   double jd_first, jd_last;
   double q_sigma = 0., weighted_rms;
   char buff[180];
   int first, last, i, n_used;
   extern const char *combine_all_observations;
   const char *packed_id;
   extern double uncertainty_parameter;

   if( combine_all_observations && *combine_all_observations)
      {
      char aligned_desig[13];

      obj_name = buff;
      packed_id = combine_all_observations;
      align_packed_desig( aligned_desig, packed_id);
      get_object_name( buff, aligned_desig);
      }
   else
      packed_id = obs->packed_id;
   fprintf( ofile, "{\n  \"num\": 1,\n");
   fprintf( ofile, "  \"ids\":\n  [\n    \"%s\"\n", obj_name);
   fprintf( ofile, "  ],\n");
   fprintf( ofile, "  \"objects\":\n  {\n");
   fprintf( ofile, "    \"%s\":\n", obj_name);
   fprintf( ofile, "    {\n      \"object\": \"%s\",\n", obj_name);
   if( *obj_name == '(')      /* numbered obj */
      {
      const char *text = find_numbered_mp_info( atoi( obj_name + 1));

      if( text && *text)
         {
         const char *fields[] = { "name", "prov_desig", "disc_date", "disc_site",
                           "disc_ref", "discover" };
         const int offsets[] = { 9, 29, 41, 53, 72, 78, 170 };

         strcpy( buff, text);
         for( i = 0; i < 6; i++)
            if( buff[offsets[i]] != ' ')
               {
               buff[offsets[i + 1] - 1] = '\0';
               remove_trailing_cr_lf( buff + offsets[i]);
               fprintf( ofile, "      \"%s\": \"%s\",\n", fields[i],
                                       buff + offsets[i]);
               }
         }
      }
   fprintf( ofile, "      \"packed\": \"%s\",\n",
                  real_packed_desig( buff, packed_id));
   fprintf( ofile, "      \"created\": %.5f,\n", jd);
   fprintf( ofile, "      \"created iso\": \"%s\",\n", iso_time( buff, jd, 0));
   fprintf( ofile, "      \"elements\":\n      {\n");
   fprintf( ofile, "        \"central body\": \"%s\",\n", object_name( buff, elem->central_obj));
   strcpy( buff, body_frame_note + 1);
   i = (int)strlen( buff);
   buff[i - 1] = '\0';    /* strip trailing paren */
   fprintf( ofile, "        \"frame\": \"%s\",\n", buff);
   fprintf( ofile, "        \"epoch_iso\": \"%s\",\n", iso_time( buff, elem->epoch, 0));
   fprintf( ofile, "        \"epoch\": %17.8f,", elem->epoch);
   if( elem->ecc < 1.)
      {
      fprintf( ofile, "\n        \"P\": %17.13f,", 2. * PI * elem->t0);
      if( !get_uncertainty( "sigma_P:", buff, 0))
         fprintf( ofile, " \"P sigma\": %s,", buff);

      fprintf( ofile, "\n        \"M\": %17.13f,",
                        centralize_ang( elem->mean_anomaly) * 180. / PI);
      if( !get_uncertainty( "sigma_M", buff, 0))
         fprintf( ofile, " \"M sigma\": %s,", buff);
      }

   fprintf( ofile, "\n        \"n\": %17.13f,", (180 / PI) / elem->t0);
   if( !get_uncertainty( "sigma_n:", buff, 0))
      fprintf( ofile, " \"n sigma\": %s,", buff);

   fprintf( ofile, "\n        \"a\": %17.13f,", elem->major_axis);
   if( !get_uncertainty( "sigma_a:", buff, 0))
      fprintf( ofile, " \"a sigma\": %s,", buff);

   fprintf( ofile, "\n        \"e\": %17.13f,", elem->ecc);
   if( !get_uncertainty( "sigma_e", buff, 0))
      fprintf( ofile, " \"e sigma\": %s,", buff);

   fprintf( ofile, "\n        \"q\": %17.13f,", elem->q);
   if( !get_uncertainty( "sigma_q", buff, 0))
      {
      q_sigma = atof( buff);
      fprintf( ofile, " \"q sigma\": %s,", buff);
      }
   if( elem->ecc < 1.)
      {
      const double big_q = elem->q * (1. + elem->ecc) / (1. - elem->ecc);

      fprintf( ofile, "\n        \"Q\": %17.13f,", big_q);
      if( !get_uncertainty( "sigma_Q", buff, 0))
         fprintf( ofile, " \"Q sigma\": %s,", buff);
      }

   fprintf( ofile, "\n        \"i\": %17.13f,", elem->incl * 180. / PI);
   if( !get_uncertainty( "sigma_i", buff, 0))
      fprintf( ofile, " \"i sigma\": %s,", buff);

   fprintf( ofile, "\n        \"arg_per\":  %17.13f,",
                        centralize_ang( elem->arg_per) * 180. / PI);
   if( !get_uncertainty( "sigma_omega", buff, 0))
      fprintf( ofile, " \"arg_per sigma\":  %s,", buff);

   fprintf( ofile, "\n        \"asc_node\": %17.13f,",
                        centralize_ang( elem->asc_node) * 180. / PI);
   if( !get_uncertainty( "sigma_Omega", buff, 0))
      fprintf( ofile, " \"asc_node sigma\": %s,", buff);

   fprintf( ofile, "\n        \"Tp\": %16.8f,", elem->perih_time);
   if( !get_uncertainty( "sigma_Tp", buff, 0))
      fprintf( ofile, " \"Tp sigma\": %s,", buff);
   fprintf( ofile, "\n        \"Tp_iso\": \"%s\",", iso_time( buff, elem->perih_time, 3));
   for( i = 0; i < n_extra_params; i++)
      {
      char tbuff[20];
      extern double solar_pressure[];

      fprintf( ofile, "\n        \"A%d\": %.9g,", i + 1, solar_pressure[i]);
      snprintf( tbuff, sizeof( tbuff), "Sigma_A%d:", i + 1);
      if( !get_uncertainty( tbuff, buff, 0))
         fprintf( ofile, " \"sigma_A%d\": %s,", i + 1, buff);
      }
   if( elem->abs_mag)
      {
      fprintf( ofile, "\n        \"H\": %6.2f,", elem->abs_mag);
      if( !get_uncertainty( "sigma_H:", buff, 0))
         fprintf( ofile, " \"H sigma\": %s,", buff);
      fprintf( ofile, "\n        \"G\": %6.2f,", elem->slope_param);
      if( !get_uncertainty( "sigma_G:", buff, 0))
         fprintf( ofile, " \"G sigma\": %s,", buff);
      }
   fprintf( ofile, "\n        \"rms_residual\": %.5g,", compute_rms( obs, n_obs));
   weighted_rms = compute_weighted_rms( obs, n_obs, &n_used);
   fprintf( ofile, "\n        \"weighted_rms_residual\": %.4f,", weighted_rms);
   fprintf( ofile, "\n        \"n_resids\": %d,", n_used);
   if( uncertainty_parameter < 90.)
      fprintf( ofile, "\n        \"U\": %.4f,", uncertainty_parameter);
   if( q_sigma && !elem->central_obj)
      {
      double neo_score = erf( (1.3 - elem->q) / q_sigma) * .5 + .5;

      fprintf( ofile, "\n        \"p_NEO\": %.4f,", neo_score * 100.);
      }
/* if( moids[1])        */
      {
      char tbuff[80];

      fprintf( ofile, "\n        \"MOIDs\":");
      fprintf( ofile, "\n        {");
      for( i = 1; i <= 8; i++)
         fprintf( ofile, "\n          \"%s\" : %.6f%c",
                     object_name( tbuff, i), moids[i],
                     (i == 8) ? ' ' : ',');
      fprintf( ofile, "\n        }");
      }

   fprintf( ofile, "\n      },\n      \"observations\":\n      {");
   get_first_and_last_included_obs( obs, n_obs, &first, &last);
   for( i = n_used = 0; i < (int)n_obs; i++)
      n_used += (obs[i].is_included & 1);
   fprintf( ofile, "\n        \"count\": %u,", n_obs);
   fprintf( ofile, "\n        \"used\": %u,", n_used);

   jd_first = obs[0].jd - td_minus_utc( obs[0].jd) / seconds_per_day;
   jd_last  = obs[n_obs - 1].jd - td_minus_utc( obs[n_obs - 1].jd) / seconds_per_day;
   fprintf( ofile, "\n        \"earliest\": %16.8f,", jd_first);
   fprintf( ofile, "\n        \"latest\": %16.8f,", jd_last);
   fprintf( ofile, "\n        \"earliest iso\": \"%s\",", iso_time( buff, jd_first, 3));
   fprintf( ofile, "\n        \"latest iso\": \"%s\",", iso_time( buff, jd_last, 3));

   jd_first = obs[first].jd - td_minus_utc( obs[first].jd) / seconds_per_day;
   jd_last  = obs[last].jd - td_minus_utc( obs[last].jd) / seconds_per_day;
   fprintf( ofile, "\n        \"earliest_used\": %16.8f,", jd_first);
   fprintf( ofile, "\n        \"latest_used\": %16.8f,", jd_last);
   fprintf( ofile, "\n        \"earliest_used iso\": \"%s\",", iso_time( buff, jd_first, 3));
   fprintf( ofile, "\n        \"latest_used iso\": \"%s\",", iso_time( buff, jd_last, 3));

   fprintf( ofile, "\n        \"residuals\":\n        [");
   for( i = 0; show_obs && i < (int)n_obs; i++)
      {
      MOTION_DETAILS m;
      int n_digits = 3;
      double total_resid;

      jd = obs[i].jd - td_minus_utc( obs[i].jd) / seconds_per_day;
      compute_observation_motion_details( obs + i, &m);
      fprintf( ofile, "\n          {\"JD\": %.6f, \"iso date\": \"%s\", \"obscode\": \"%s\",",
                  jd, iso_time( buff, jd, 3), obs[i].mpc_code);
      fprintf( ofile, "\n                 \"RA\" : %.6f, \"Dec\": %.6f,",
                  obs[i].ra * 180. / PI, obs[i].dec * 180. / PI);
      total_resid = sqrt( m.xresid * m.xresid + m.yresid * m.yresid);
      while( n_digits < 12 && total_resid < .1)
         {
         total_resid *= 10.;
         n_digits++;
         }
      fprintf( ofile, "\n                 \"dRA\" : %.*f, \"dDec\": %.*f, \"dTime\": %.3f, \"cross\": %.*f,",
                    n_digits, m.xresid, n_digits, m.yresid, m.time_residual,
                    n_digits, m.cross_residual);
      fprintf( ofile, "\n                 \"reference\" : \"%s\",", obs[i].reference);
      fprintf( ofile, "  \"packed\" : \"%s\",", obs[i].packed_id);
      fprintf( ofile, "\n                 \"sigma_1\" : %f,", obs[i].posn_sigma_1);
      fprintf( ofile, "\n                 \"sigma_2\" : %f,", obs[i].posn_sigma_2);
      fprintf( ofile, "\n                 \"posn_sigma_theta\" : %f,",
                                                      obs[i].posn_sigma_theta * 180. / PI);
      fprintf( ofile, "\n                 \"RAMotion\" : %.4f,", m.ra_motion);
      fprintf( ofile, "\n                 \"decMotion\" : %.4f,", m.dec_motion);
      fprintf( ofile, "\n                 \"TotalMotion\" : %.4f,", m.total_motion);
      fprintf( ofile, "\n                 \"PAMotion\" : %.2f,", m.position_angle_of_motion);
      if( obs[i].obs_mag != BLANK_MAG)
         {
         fprintf( ofile, " \"dMag\" : %.2f,", obs[i].obs_mag - obs[i].computed_mag);
         fprintf( ofile, " \"MagObs\" : %.2f,", obs[i].obs_mag);
         fprintf( ofile, "\n                 \"MagBand\" : \"%c\",", obs[i].mag_band);
         fprintf( ofile, " \"MagSigma\" : %.4f,", obs[i].mag_sigma);
         }
      fprintf( ofile, "\n                 \"net\" : \"%c\",", obs[i].astrometric_net_code);
      if( obs[i].note1 == '\\' || obs[i].note1 == '\"')
         {
         buff[0] = '\\';
         buff[1] = obs[i].note1;
         buff[2] = '\0';
         }
      else
         {
         buff[0] = obs[i].note1;
         buff[1] = '\0';
         }
      fprintf( ofile, "\n                 \"note1\" : \"%s\",", buff);
      fprintf( ofile, "\n                 \"note2\" : \"%c\",", obs[i].note2);
      fprintf( ofile, "\n                 \"discovery_asterisk\" : \"%c\",",
                                       obs[i].discovery_asterisk);
      fprintf( ofile, "\n                 \"flags\" : %d,", obs[i].flags);
      fprintf( ofile, "\n                 \"incl\" : %d", obs[i].is_included);
      fprintf( ofile, " }%c", (i == (int)n_obs - 1 ? ' ' : ','));
      }
   fprintf( ofile, "\n        ]\n      }");
   fprintf( ofile, "\n    }\n  }\n}\n");
   return( 0);
}

static int elements_in_guide_format( char *buff, const ELEMENTS *elem,
                     const char *obj_name, const OBSERVE *obs,
                     const unsigned n_obs)
{
   int month;
   double day;
   long year;

   day = decimal_day_to_dmy( elem->perih_time, &year, &month,
                                              CALENDAR_JULIAN_GREGORIAN);
            /*      name day  mon yr MA      q      e */
   sprintf( buff, "%-43s%8.5f%3d%5ld Find_Orb %14.7f%12.7f%11.6f%12.6f%12.6f",
            obj_name, day, month, year,
            elem->q, elem->ecc,
            centralize_ang( elem->incl) * 180. / PI,
            centralize_ang( elem->arg_per) * 180. / PI,
            centralize_ang( elem->asc_node) * 180. / PI);
   if( elem->q < .01)
      {
      sprintf( buff + 71, "%12.10f", elem->q);
      buff[71] = buff[83] = ' ';
      }
   sprintf( buff + strlen( buff), " %9.1f%5.1f%5.1f %c",
            elem->epoch, elem->abs_mag,
            elem->slope_param * (elem->is_asteroid ? 1. : 0.4),
            (elem->is_asteroid ? 'A' : ' '));
   if( elem->central_obj)
      sprintf( buff + strlen( buff), "  Center: %d", elem->central_obj);
   strcat( buff, "  ");
   observation_summary_data( buff + strlen( buff), obs, n_obs, -1);
   return( 0);
}

static int is_cometary( const char *constraints)
{
   const char *ecc = strstr( constraints, "e=1");

   return( ecc && atof( ecc + 2) == 1.);
}

int monte_carlo_object_count = 0;
int n_monte_carlo_impactors = 0;
int append_elements_to_element_file = 0;
int using_sr = 0;
char orbit_summary_text[120];
double max_monte_rms;

void set_statistical_ranging( const int new_using_sr)
{
   using_sr = new_using_sr;
   n_monte_carlo_impactors = monte_carlo_object_count = 0;
}

static size_t space_pad_buffer( char *buff, const size_t length)
{
   const size_t rval = strlen( buff);

   if( rval < length)
      {
      memset( buff + rval, ' ', length - rval);
      buff[length] = '\0';
      }
   return( rval);
}

static int show_reference( char *buff)
{
   const size_t reference_loc = 62;
   const int rval = (space_pad_buffer( buff, reference_loc) <= reference_loc);

   if( rval)
      {
      const char *reference = get_environment_ptr( "REFERENCE");

      if( !*reference)
         reference = "Find_Orb";
      strcpy( buff + reference_loc, reference);
      }
   return( rval);
}

extern int available_sigmas;

int compute_available_sigmas_hash( const OBSERVE FAR *obs, const int n_obs,
         const double epoch, const unsigned perturbers, const int central_obj);

static int get_uncertainty( const char *key, char *obuff, const bool in_km)
{
   int rval = -1;
   FILE *ifile;
   const char *filenames[4] = { NULL, "covar.txt", "monte.txt", "monte.txt" };
   char buff[100];

   *obuff = '\0';
   if( available_sigmas && (ifile = fopen_ext(
                  get_file_name( buff, filenames[available_sigmas]), "tcrb")) != NULL)
      {
      const size_t keylen = strlen( key);

      while( rval && fgets( buff, sizeof( buff), ifile))
         if( !memcmp( buff, key, keylen) && buff[keylen] == ' ')
            {
            size_t loc = keylen;

            rval = 0;
            if( in_km)
               {
               char *tptr = strchr( buff, '(');

               if( tptr)
                  loc = tptr - buff + 1;
               }
            sscanf( buff + loc, "%20s", obuff);
            }
      fclose( ifile);
      }
   return( rval);
}


static void consider_replacing( char *buff, const char *search_text,
                                            const char *sigma_key)
{
   char *tptr = strstr( buff, search_text);

   if( tptr && available_sigmas)
      {
      int at_end_of_line;
      bool is_km = false;

      tptr += strlen( search_text);
      while( *tptr == ' ')    /* scan over spaces... */
         tptr++;
      while( *tptr && *tptr != ' ')   /* ...then over the actual value */
         tptr++;
      remove_trailing_cr_lf( tptr);
      at_end_of_line = (*tptr == '\0');
      if( tptr[-2] == 'k' && tptr[-1] == 'm')
         {
         is_km = true;
         tptr -= 2;
         }
      strcpy( tptr, " +/- ");
      tptr += 5;
      get_uncertainty( sigma_key, tptr, is_km);
      if( !at_end_of_line)
         tptr[strlen( tptr)] = ' ';
      }
}

/* For display purposes,  it can be useful to switch from the default
   J2000 ecliptic frame to a planet-centered frame.  In such a frame,
   the z-axis points in the direction of the planet's north pole;  the
   x-axis is the cross-product of z with the J2000 pole of the earth;
   and the y-axis is the cross-product of x and z. */

static void ecliptic_to_planetary_plane( const int planet_no,
               const double epoch_jd, double *state_vect)
{
   double planet_matrix[9], xform[3][3];
   int i;
   double tval;

   calc_planet_orientation( planet_no, 0, epoch_jd, planet_matrix);
         /* At this point,  planet_matrix[6, 7, 8] is a J2000 equatorial */
         /* vector pointing in the direction of the planet's north pole. */
         /* Copy that as our z-axis,  xform[2]: */
   memcpy( xform[2], planet_matrix + 6, 3 * sizeof( double));
         /* Our "X-axis" is perpendicular to "Z",  but in the plane */
         /* of the J2000 equator,  corresponding to the ascending node  */
         /* of the planet's equator relative to the J2000 equator:      */
   tval = sqrt( xform[2][0] * xform[2][0] + xform[2][1] * xform[2][1]);
   xform[0][0] = -xform[2][1] / tval;
   xform[0][1] =  xform[2][0] / tval;
   xform[0][2] = 0.;
         /* So we've got two of the three base vectors,  still in J2000 */
         /* equatorial.  Transform them to ecliptic... */
   equatorial_to_ecliptic( xform[0]);
   equatorial_to_ecliptic( xform[2]);
         /* ...and 'Y' is simply the cross-product of 'X' and 'Z': */
   vector_cross_product( xform[1], xform[2], xform[0]);
         /* OK.  Now we're ready to transform the state vector : */
   for( i = 0; i < 2; i++, state_vect += 3)
      {
      const double x = dot_product( state_vect, xform[0]);
      const double y = dot_product( state_vect, xform[1]);
      const double z = dot_product( state_vect, xform[2]);

      state_vect[0] = x;
      state_vect[1] = y;
      state_vect[2] = z;
      }
}

   /* The following will revise text such as
   "value=3.141e-009 +/- 2.34e-011" to read
   "value=3.141e-9 +/- 2.34e-11".         */

static void clobber_leading_zeroes_in_exponent( char *buff)
{
   while( *buff)
      {
      if( *buff == 'e' || *buff == 'E')
         if( buff[1] == '+' || buff[1] == '-')
            if( buff[2] == '0')
               {
               memmove( buff + 2, buff + 3, strlen( buff + 2));
               buff--;
               }
      buff++;
      }
}

#define MAX_SOF_LEN 400

char *get_file_name( char *filename, const char *template_file_name);

/* Write out the elements in SOF (Standard Orbit Format) at the end of an
existing file,  one in which the first line is the header.      */

static int add_sof_to_file( const char *filename,
             const ELEMENTS *elem,
             const int n_obs, const OBSERVE *obs)
{
   char templat[MAX_SOF_LEN], obuff[MAX_SOF_LEN];
   char output_filename[100];
   const char *obj_name = obs->packed_id;
   FILE *fp;
   int rval = -1, forking;

   while( *obj_name == ' ')
      obj_name++;
   fp = fopen_ext( filename, "cr+b");
   get_file_name( output_filename, filename);
   forking = strcmp( output_filename, filename);
   if( fp)
      {
      bool got_it = false;

      fseek( fp, 0L, SEEK_SET);
      if( !fgets( templat, sizeof( templat), fp))
         {
         fclose( fp);
         return( -1);
         }
      if( forking)
         {
         fclose( fp);
         fp = fopen_ext( output_filename, "cr+b");
         assert( fp);
         }
      fseek( fp, 0L, SEEK_SET);
      while( !got_it && fgets( obuff, sizeof( obuff), fp))
         if( !memcmp( obuff, obj_name, strlen( obj_name)) &&
                     obuff[strlen( obj_name)] == ' ')
            {
            got_it = true;
            fseek( fp, -(long)strlen( obuff), SEEK_CUR);
            }
      rval = put_elements_into_sof( obuff, templat, elem, n_obs, obs);
      fwrite( obuff, strlen( obuff), 1, fp);
      fclose( fp);
      }
   return( rval);
}

/* Code to figure out the delta-V required to rendezvous with an object
using Shoemaker and Helin's method,  1978,  Earth-approaching
asteroids as targets for exploration, NASA CP-2053, pp. 245-256.
This is basically a translation into C of some IDL code written by
Nick Moskovitz.          */
#ifdef NOT_READY_YET
static double shoemaker_helin_encounter_velocity( const ELEMENTS *elem)
{
   const double norm = 29.784;     /* Normalization factor = Earth's orbital speed [km/s] */
// const double U_0 = 7.94 / norm; /* Normalized orbital speed in LEO */
   const double U_0 = 7.727 / norm; /* Normalized orbital speed in LEO at 300km */
   const double sqrt_2 = 1.41421358;
   const double S = sqrt_2 * U_0   /* Normalized escape speed from Earth (OLD: 11.2 / norm) */
   const double Q = elem->major * (1. + elem->ecc);
   const double cos_half_incl = cos( elem->incl / 2.);
   double u_t2, u_c2, u_r2;

   if( elem->major > 1. && elem->q > .983)        /* Aten */
      {
      u_t2 = 2. - 2. * sqrt( 2. * Q - Q * Q) * cos_half_incl;
      u_c2 = 3. / Q - 1. - 2. * sqrt( 2. - Q) / Q;
      u_r2 = 3. / q - 1. / elem->major - 2. * sqrt( elem->major * (1.
      ))
      }
}

#endif

void get_periapsis_loc( double *ecliptic_lon, double *ecliptic_lat,
             const ELEMENTS *elem)
{
   *ecliptic_lon = elem->asc_node +
         atan2( cos( elem->incl) * sin( elem->arg_per), cos( elem->arg_per));
   *ecliptic_lon = centralize_ang( *ecliptic_lon);
   *ecliptic_lat = asin( sin( elem->incl) * sin( elem->arg_per));
}

/* From an e-mail from Alan Harris:

   "...the formula for encounter velocity (in FORTRAN), for a circular planet
orbit is:

DV=30.*SQRT(3.-A0/A-2.*SQRT(A*(1.-E**2)/A0)*COS(XI/57.2958))

Where A0 is the planet orbit SMA (1.0 for Earth), and A, E, XI are (a,e,i)
of the crossing body.  If DV is less than about 2.5 [km/s], the crossing body
can't make it to/from Mars or Venus no matter what direction it is going,
at greater than 2.5, it can.  That sort of provides the cutoff for natural
bodies crossing the Earth orbit.  There are a few exceptions, but they are
probably moon ejecta or old rocket cans." (Note: 57.2958 = 180/pi = degrees
per radian conversion,  not needed in C.)

   I've seen some inaccuracies in this formula,  though,  which I _think_
reflect the fact that it assumes the earth's orbit is circular.  Alan
points out that this really should only apply to crossing orbits (which,
with earth's orbit considered circular,  means q < 1 < Q.)  If the MOID
is non-zero,  it's hard to say exactly what the "encounter velocity" means.

   Note furthermore that the above can be reduced to

dv = v0 * sqrt(3. - tisserand)

   where v0 = orbital speed for the earth (roughly 30 km/s) and 'tisserand'
is the usual Tisserand criterion.           */

static double encounter_velocity( const ELEMENTS *elem, const double a0)
{
   const double a = elem->major_axis;
   const double tval = sqrt( a * (1. - elem->ecc * elem->ecc) / a0);
   const double tisserand = a0 / a + 2. * tval * cos( elem->incl);
   double rval;

   if( tisserand > 3.)    /* can happen if the orbits can't really */
      rval = 0.;          /* intersect (i.e.,  q > 1 or Q < 1) */
   else
      rval = 30. * sqrt( 3. - tisserand);
   return( rval);
}

/* Used to provide a link to Tony Dunn's Orbit Simulator when making
pseudo-MPECs.  Epoch is stored in the sixth element of the array. */

double helio_ecliptic_j2000_vect[7];

/* The results from write_out_elements_to_file() can be somewhat
varied.  The output for elliptical and parabolic/hyperbolic orbits are
very different.  Asteroids and comets differ in whether H and G are
given,  or M(N)/M(T) and K.  Or those fields can be blank if no
magnitude data is given.

   Constraints or Monte Carlo data can modify the perihelion line,  such as

   Perihelion 1998 Mar 16.601688 TT;  Constraint: a=2.5
   Perihelion 1998 Mar 16.601688 TT;  20.3% impact (203/1000)

   AMR data appears on the epoch line:

Epoch 1997 Oct 13.0 TT; AMR 0.034 m^2/kg

   MOIDs can carry on to next line,  wiping out P Q header and even
the (2000.0) text if there are enough extra MOIDs.

   Examples:

Orbital elements:
1996 XX1
   Perihelion 1998 Mar 8.969329 TT = 23:15:50 (JD 2450881.469329)
Epoch 1997 Oct 13.0 TT = JDT 2450734.5   Earth MOID: 0.0615   Ma: 0.0049
M 322.34260              (2000.0)            P               Q
n   0.25622619     Peri.   72.47318     -0.50277681     -0.86441429
a   2.45501241     Node    47.71084      0.79213697     -0.46159019
e   0.5741560      Incl.    0.14296      0.34602670     -0.19930484
P   3.85           H   15.8     U  8.4     q 1.04545221  Q 3.86457261
From 13 observations 1997 Oct. 12-22; mean residual 0".485.

Orbital elements:
2009 BD
   Perigee 2009 Jan 25.247396 TT; A1=5.76e-11, A2=-1.29e-12
Epoch 2009 Jan 22.0 TT = JDT 2454853.5                  Find_Orb
q685927.806km       (2000.0)            P               Q
H   28.6  G 0.15   Peri.   86.73833      0.01481057     -0.49410657
                   Node   119.50263     -0.32898297      0.81855881
e   2.9230182      Incl.   92.82529      0.94421970      0.29295079
From 173 observations 2009 Jan. 16-2011 June 20; mean residual 0".294.

Orbital elements:
1997 ZZ99
   Perilune 1997 Apr 22.543629 TT = 13:02:49 (JD 2450561.043629)
Epoch 1997 Apr 22.5 TT = JDT 2450561.0                  Find_Orb
q 24606.5028km           (2000.0)            P               Q
H   22.9  G 0.15   Peri.  111.16531      0.57204494     -0.81965341
                   Node   193.85513     -0.79189676     -0.54220560
e 659.1509995      Incl.    7.32769     -0.21369156     -0.18488202
From 35 observations 1997 Apr. 21-22; mean residual 1".128.

Possible new format,  allowing inclusion of uncertainties:

Orbital elements: 1996 XX1
   Perihelion 1998 Mar 8.977252 +/- 0.109 TT (JD 2450881.477252)
Epoch 1997 Oct 12.0 TT = JDT 2450733.5   Earth MOID: 0.0613   Ma: 0.0049
M     322.109164 +/-   0.0574     (More MOIDs or area/mass goes here)
n   0.256058518 +/- 0.000205      Peri.  72.530699 +/-   0.0861
a   2.456084084 +/- 0.00131       Node   47.657952 +/-   0.0576
e   0.57438774 +/- 0.000178       Incl.   0.143216 +/-   0.000294
H   16.6    G 0.15     U  7.1   P 3.85012 +/- .00021
q 1.045339477 +/- 0.000139      Q 3.866828691 +/- 0.0025
From 13 of 17 observations 1997 Oct. 12-22; mean residual 0".485.

   ...and the hyperbolic version:

Orbital elements: 1997 ZZ99
   Perilune 1997 Apr 22.543629 +/- 0.021 TT (JD 2450561.043629)
Epoch 1997 Apr 22.5 TT = JDT 2450561.0                  Find_Orb
q 24606.5028 +/- 3.345 km         (More MOIDs or area/mass here)
H   22.9  G 0.15                  Peri.  111.16531 +/- 0.0321
z 934.123 +/- 3.14159             Node   193.85513 +/- 0.145
e 659.1509995 +/- 3.456           Incl.    7.32769 +/- 0.0056
From 33 of 35 observations 1997 Apr. 21-22; mean residual 1".128.

   ....wherein Q, a, P, n, U, and M are omitted as no longer meaningful,
but we show z=1/a,  which is sometimes helpful with near-parabolic cases.
Actually,  we should always show z if the uncertainty in a is comparable
to a itself,  for all types of orbits.

   We basically retain the first three lines of the existing format,
except that the second line needs a little revision to contain the
periapsis uncertainty.  And we retain the "from x to y" line.

   We should also "prune" quantities so the number of digits in the sigma
matches the number of digits in the quantity itself.

   The program is theoretically capable of computing MOIDs relative to
N_MOIDS objects (currently 14:  eight planets,  six asteroids). However,
N_MOIDS_TO_SHOW = 8 at present (we only show planetary MOIDs.) */

#define N_MOIDS           14
#define N_MOIDS_TO_SHOW    8

double comet_total_magnitude = 0.;          /* a.k.a. "M1" */
double comet_nuclear_magnitude = 0.;        /* a.k.a. "M2" */
bool saving_elements_for_reuse = false;

#define ELEMENT_FRAME_DEFAULT                   0
#define ELEMENT_FRAME_J2000_ECLIPTIC            1
#define ELEMENT_FRAME_J2000_EQUATORIAL          2
#define ELEMENT_FRAME_BODY_FRAME                3

#define ORBIT_SUMMARY_Q                1
#define ORBIT_SUMMARY_H                2

int write_out_elements_to_file( const double *orbit,
            const double curr_epoch,
            const double epoch_shown,
            OBSERVE FAR *obs, const int n_obs, const char *constraints,
            const int precision, const int monte_carlo,
            const int options)
{
   char object_name[80], buff[260], more_moids[80];
   const char *file_permits = (append_elements_to_element_file ? "tfca" : "tfcw+");
   extern const char *elements_filename;
   FILE *ofile = fopen_ext( get_file_name( buff, elements_filename), file_permits);
   double rel_orbit[6], orbit2[6];
   int planet_orbiting, n_lines, i, bad_elements = 0;
   ELEMENTS elem, helio_elem;
   char *tptr, *tbuff;
   char impact_buff[80], elem_epoch_text[30];
   int n_more_moids = 0;
   int output_format = (precision | SHOWELEM_PERIH_TIME_MASK);
   extern int n_extra_params;
   extern unsigned perturbers;
   int reference_shown = 0;
   double moids[N_MOIDS + 1];
   double j2000_ecliptic_rel_orbit[6];
   double barbee_style_delta_v = 0.;   /* see 'moid4.cpp' */
   const char *monte_carlo_permits;
   const bool rms_ok = (compute_rms( obs, n_obs) < max_monte_rms);
   extern int available_sigmas;
   const char *body_frame_note = NULL;
   bool body_frame_note_shown = false;
   int showing_sigmas = available_sigmas;
   int elements_frame = atoi( get_environment_ptr( "ELEMENTS_FRAME"));
   const unsigned orbit_summary_options = atoi( get_environment_ptr( "ORBIT_SUMMARY_OPTIONS"));
   extern int is_interstellar;         /* orb_func.cpp */

   setvbuf( ofile, NULL, _IONBF, 0);
   if( default_comet_magnitude_type == 'N')
      output_format |= SHOWELEM_COMET_MAGS_NUCLEAR;
   if (options & ELEM_OUT_ALTERNATIVE_FORMAT)
      output_format |= SHOWELEM_OMIT_PQ_MASK;
   fprintf( ofile, "%s %s",
         get_find_orb_text( 99174),        /* "Orbital elements:" */
         options & ELEM_OUT_ALTERNATIVE_FORMAT ? " " : "\n");
   get_object_name( object_name, obs->packed_id);
   if( !curr_epoch || !epoch_shown)
      {
      fprintf( ofile, "%s\nNo elements available\n", object_name);
      observation_summary_data( buff, obs, n_obs, options);
      fprintf( ofile, "%s\n", buff);
      fclose( ofile);
      return( -1);
      }
   memcpy( orbit2, orbit, 6 * sizeof( double));
   integrate_orbit( orbit2, curr_epoch, epoch_shown);
   memcpy( helio_ecliptic_j2000_vect, orbit2, 6 * sizeof( double));
   helio_ecliptic_j2000_vect[6] = epoch_shown;
   if( options & ELEM_OUT_HELIOCENTRIC_ONLY)
      {
      planet_orbiting = forced_central_body;
      get_relative_vector( epoch_shown, orbit2, rel_orbit, planet_orbiting);
      }
   else
      planet_orbiting = find_best_fit_planet( epoch_shown, orbit2, rel_orbit);

   memcpy( j2000_ecliptic_rel_orbit, rel_orbit, 6 * sizeof( double));
            /* By default,  we use J2000 equatorial elements for geocentric
            elements,  J2000 ecliptic for everybody else. */
   if( elements_frame == ELEMENT_FRAME_DEFAULT)
      elements_frame = ((planet_orbiting == 3) ?
                  ELEMENT_FRAME_J2000_EQUATORIAL :
                  ELEMENT_FRAME_J2000_ECLIPTIC);

   if( elements_frame == ELEMENT_FRAME_J2000_ECLIPTIC)
      {
      const char *elem_epoch = get_environment_ptr( "ELEMENT_EPOCH");

      if( !*elem_epoch)
         body_frame_note = "(J2000 ecliptic)";
      else
         {
         double year = atof( elem_epoch + 1);
         double precession_matrix[9];

         body_frame_note = elem_epoch_text;
         if( !year)
             year = JD_TO_YEAR( epoch_shown);
         snprintf( elem_epoch_text, sizeof( elem_epoch_text), "(%s)", elem_epoch);
         setup_ecliptic_precession( precession_matrix, 2000., year);
         precess_vector( precession_matrix, rel_orbit, rel_orbit);
         precess_vector( precession_matrix, rel_orbit + 3, rel_orbit + 3);
         }
      }
   if( elements_frame == ELEMENT_FRAME_J2000_EQUATORIAL)
      {
      ecliptic_to_equatorial( rel_orbit);
      ecliptic_to_equatorial( rel_orbit + 3);
      body_frame_note = "(J2000 equator)";
      }
   if( elements_frame == ELEMENT_FRAME_BODY_FRAME)
      {
      ecliptic_to_planetary_plane(
                  (planet_orbiting == 100 ? -1 : planet_orbiting),
                  epoch_shown, rel_orbit);
      body_frame_note = "(body frame)";
      }

   if( !(options & ELEM_OUT_ALTERNATIVE_FORMAT))
      showing_sigmas = 0;
   if( showing_sigmas == COVARIANCE_AVAILABLE)
      {
      extern int available_sigmas_hash;

      if( available_sigmas_hash != compute_available_sigmas_hash( obs, n_obs,
                  epoch_shown, perturbers, planet_orbiting))
         showing_sigmas = 0;
      }

   elem.central_obj = planet_orbiting;
   elem.gm = get_planet_mass( planet_orbiting);
   elem.epoch = epoch_shown;
   calc_classical_elements( &elem, rel_orbit, epoch_shown, 1);
   if( elem.ecc < .9)
      sprintf( orbit_summary_text, "a=%.3f, ", elem.major_axis);
   else
      *orbit_summary_text = '\0';
   if( elem.ecc >= .9 || (orbit_summary_options & ORBIT_SUMMARY_Q))
      snprintf_append( orbit_summary_text, sizeof( orbit_summary_text),
                                   "q=%.3f, ", elem.q);
   snprintf_append( orbit_summary_text, sizeof( orbit_summary_text),
            "e=%.3f, i=%d", elem.ecc, (int)( elem.incl * 180. / PI + .5));
   elem.is_asteroid = (object_type == OBJECT_TYPE_ASTEROID);
   if( elem.is_asteroid)
      {
      elem.slope_param = asteroid_magnitude_slope_param;
      elem.abs_mag = calc_absolute_magnitude( obs, n_obs);
      if( orbit_summary_options & ORBIT_SUMMARY_H)
         snprintf_append( orbit_summary_text, sizeof( orbit_summary_text),
                  " H=%.1f", elem.abs_mag);
      }
   else              /* for comets, compute nuclear & total absolute mags */
      {
      elem.slope_param = comet_magnitude_slope_param;
      for( i = 0; i < 2; i++)
         {
         default_comet_magnitude_type = 'T' + 'N' - default_comet_magnitude_type;
         elem.abs_mag = calc_absolute_magnitude( obs, n_obs);
         if( default_comet_magnitude_type == 'N')
            comet_nuclear_magnitude = elem.abs_mag;
         else
            comet_total_magnitude = elem.abs_mag;
         }
      }

   add_sof_to_file( (n_extra_params >= 2 ? "cmt_sof.txt" : sof_filename),
                    &elem, n_obs, obs);            /* elem_ou2.cpp */
   if( saving_elements_for_reuse && available_sigmas == COVARIANCE_AVAILABLE)
      add_sof_to_file( "orbits.sof", &elem, n_obs, obs);
/* if( showing_sigmas == COVARIANCE_AVAILABLE)
*/    {
      ELEMENTS elem2 = elem;
      double rel_orbit2[6];

      if( showing_sigmas == COVARIANCE_AVAILABLE)
         {
         compute_variant_orbit( rel_orbit2, rel_orbit, 1.);    /* orb_func.cpp */
         calc_classical_elements( &elem2, rel_orbit2, epoch_shown, 1);
         }
      else        /* insert dummy elements */
         elem2.q = elem2.ecc = elem2.incl = elem2.arg_per = elem2.asc_node = 0.;
      add_sof_to_file( sofv_filename, &elem2, n_obs, obs);     /* elem_ou2.cpp */
      }
   helio_elem = elem;            /* Heliocentric J2000 ecliptic elems */
   helio_elem.central_obj = 0;
   helio_elem.gm = SOLAR_GM;
   calc_classical_elements( &helio_elem, orbit2, epoch_shown, 1);
   tbuff = (char *)malloc( 80 * 9);
   n_lines = elements_in_mpc_format( tbuff, &elem, object_name,
               is_cometary( constraints) && fabs( elem.ecc - 1.) < 1.e-6,
               output_format);
   fprintf( ofile, "%s\n", tbuff);
   tptr = tbuff + strlen( tbuff) + 1;
   *more_moids = '\0';
   for( i = 0; i < 9; i++)
      moids[i] = 0.;
   for( i = 1; i < n_lines && *tptr; i++)
      {
      char *tt_ptr;
      char sigma_buff[80];
      extern double solar_pressure[], uncertainty_parameter;
             /* "Solar radiation pressure at 1 AU",  in             */
             /* kg*AU^3 / (m^2*d^2),  from a private communication  */
             /* from Steve Chesley; see orb_func.cpp for details    */

      strcpy( sigma_buff, "+/- ");
      strcpy( buff, tptr);
      tt_ptr = strstr( buff, "TT") + 2;
      if( !memcmp( buff, "   Peri", 7))
         {
         assert( tt_ptr);
         if( *constraints)
            {
            if( constraints[0] == 'm' && constraints[1] == '(')
               {
               double *mass = get_asteroid_mass( atoi( constraints + 2));

               if( mass)
                  {
                  sprintf( tt_ptr, "; %s=%.4g", constraints, *mass);
                  consider_replacing( buff, constraints, "Sigma_mass:");
                  clobber_leading_zeroes_in_exponent( buff);
                  }
               else
                  sprintf( tt_ptr, "; bad '%s'", constraints);
               }
            else
               sprintf( tt_ptr, ";  Constraint: %s", constraints);
            }
         else if( n_monte_carlo_impactors && monte_carlo)
            sprintf( tt_ptr, ";  %.2f%% impact (%d/%d)",
                100. * (double)n_monte_carlo_impactors /
                       (double)monte_carlo_object_count,
                       n_monte_carlo_impactors, monte_carlo_object_count);
         if( showing_sigmas)
            if( !get_uncertainty( "sigma_Tp", sigma_buff + 4, false))
               {
               strcat( sigma_buff, " TT");
               text_search_and_replace( buff, "TT", sigma_buff);
               }

         if( n_extra_params >= 2 && n_extra_params <= 4)
            {
            char tbuff0[80], sig_name[20];
            int j = 0;

            if( !strcmp( constraints, "A=0") || !strcmp( constraints, "A1=0"))
               j = 1;      /* A1 constrained to be zero -> don't show it */
            strcat( tt_ptr, "\n");
            tt_ptr += strlen( tt_ptr);
            for( ; j < n_extra_params; j++)
               {
               put_double_in_buff( tbuff0, solar_pressure[j]);
               text_search_and_replace( tbuff0, " ", "");
               strcat( tbuff0, " ");
               sprintf( sig_name, "Sigma_A%d:", j + 1);
               if( showing_sigmas)
                  if( !get_uncertainty( sig_name, sigma_buff + 4, false))
                     strcat( tbuff0, sigma_buff);
/*             snprintf_append( tt_ptr, 180, "A%d: %s", j + 1, tbuff0); */
               if( j == 3)
                  snprintf_append( buff, sizeof( buff), "DT: %s", tbuff0);
               else
                  {
                  snprintf_append( buff, sizeof( buff), "A%d: %s", j + 1, tbuff0);
                  if( j == n_extra_params - 1 || j == 2)
                     {
                     strcat( tt_ptr, " AU/day^2");
                     if( is_inverse_square_force_model( ))
                        strcat( tt_ptr, " [1/r^2]");
                     }
                  }
               if( j < n_extra_params - 1)
                  strcat( tt_ptr, (strlen( tt_ptr) > 50) ? "\n" : "   ");
               }
            }
         assert( strlen( buff) < sizeof( buff) - 1);
         }

      else if( !memcmp( buff, "Epoch", 5))
         {
         int j;
         const double SRP1AU = 2.3e-7;

         if( n_extra_params == 1)
            {
            sprintf( tt_ptr, "; AMR %.5g",
                                 solar_pressure[0] * SOLAR_GM / SRP1AU);
            if( showing_sigmas)
               consider_replacing( tt_ptr, "AMR", "Sigma_AMR:");
            strcat( tt_ptr, " m^2/kg");
            }
         if( !planet_orbiting)
            for( j = 0; j < N_MOIDS_TO_SHOW; j++)
               {
               static const char moid_idx[N_MOIDS] = { 3, 5, 2, 1, 4, 6, 7, 8,
                                       10, 11, 12, 13, 14, 15 };
               moid_data_t mdata;
               double moid, moid_limit = .1;
               ELEMENTS planet_elem;
               const int forced_moid =
                   (atoi( get_environment_ptr( "MOIDS")) >> j) & 1;

               setup_planet_elem( &planet_elem, moid_idx[j],
                                (epoch_shown - J2000) / 36525.);
               moid = find_moid_full( &planet_elem, &helio_elem, &mdata);
               if( !j)     /* only get Barbee speed for the earth */
                  barbee_style_delta_v = mdata.barbee_speed * AU_IN_KM / seconds_per_day;
               if( j < 2)        /* Earth or Jupiter */
                  moid_limit = 1.;
               else if( j > 7)            /* asteroid */
                  moid_limit = .1;
               else if( j > 4)          /* Saturn,  Uranus,  Neptune */
                  moid_limit = 1.;
               moids[(int)moid_idx[j]] = moid;
               if( forced_moid || moid < moid_limit)
                  {
                  char addendum[30];
                  static const char *moid_text[N_MOIDS] = { "Earth MOID", "Ju",
                           "Ve", "Me", "Ma", "Sa", "Ur", "Ne",
                           "Ce", "Pa", "Vt", "(29)", "(16)", "(15)" };

                  sprintf( addendum, "   %s: %.4f", moid_text[j], moid);
                  if( strlen( addendum) + strlen( buff) < 79)
                     strcat( buff, addendum);
                  else
                     {
                     if( n_more_moids < 3)
                        strcat( more_moids, addendum);
                     n_more_moids++;
                     }
                  if( !j && moid < .5)
                     snprintf_append( orbit_summary_text,
                           sizeof( orbit_summary_text), " MOID %.3f", moid);
                  }
               }
         reference_shown = show_reference( buff);
         }
      else if( *more_moids)
         {
         space_pad_buffer( buff, 33);
         strcpy( buff + 33, more_moids);
         *more_moids = '\0';
         if( !reference_shown)
            reference_shown = show_reference( buff);
         }
      else if( !reference_shown && *buff == 'M' && buff[1] != '(' &&
                       !available_sigmas &&
                       (options & ELEM_OUT_ALTERNATIVE_FORMAT))
         {                        /* no ref shown yet: try fitting on next */
         buff[39] = '\0';         /* line, if using alternative format     */
         }
      if( !(options & ELEM_OUT_ALTERNATIVE_FORMAT))
         if( !memcmp( buff + 31, " G ", 3) && uncertainty_parameter < 90.)
            {                       /* replace slope w/uncertainty */
//          sprintf( buff + 32, "U%7.3f  ", uncertainty_parameter);
            sprintf( buff + 32, "U%5.1f  ", uncertainty_parameter);
            buff[40] = ' ';
            }
      if( showing_sigmas)
         {
         if( i >= 4 && i <= 6)        /* lines w/Peri., Node, & Incl. */
            {
            char *tptr = buff + (precision < 10 ? 19 : 9 + precision);

            memmove( tptr + 17, tptr, strlen( tptr) + 1);
            memset( tptr, ' ', 17);
            consider_replacing( buff, "Peri.", "sigma_omega");
            consider_replacing( buff, "Node", "sigma_Omega");
            consider_replacing( buff, "Incl.", "sigma_i");
            remove_trailing_cr_lf( buff); /* also removes trailing spaces */
            if( *buff == 'e')
               consider_replacing( buff, "e", "sigma_e");
            else if( *buff == 'a')
               consider_replacing( buff, "a", "sigma_a:");
            else if( *buff == 'n')
               consider_replacing( buff, "n", "sigma_n:");
            }
         switch( buff[0])
            {
            case 'M':
               if( buff[1] != '(')         /* don't do this if the 'M' */
                  {                        /* refers to a comet mag */
                  consider_replacing( buff, "M", "sigma_M");
                  if( !reference_shown)
                     reference_shown = show_reference( buff);
                  }
               break;
            case 'q':
               consider_replacing( buff, "q", "sigma_q");
               break;
            case 'P':
               {
               char *zptr = strchr( buff, 'Q');
               char tbuff[80];

               if( zptr)
                  {
                  strcpy( tbuff, zptr);
                  consider_replacing( tbuff, "Q", "sigma_Q");
                  zptr[-1] = '\0';
                  }
               else                 /* no aphelion given */
                  *tbuff = '\0';
               zptr = strchr( buff, 'q');
               if( zptr)
                  {
                  char phg_line[80];   /* contains P, H, G,  sometimes U text */

                  zptr[-1] = '\0';
                  if( buff[32] == 'G' && buff[19] == 'H')
                     {
                     assert( zptr - buff == 43);
                     memcpy( phg_line, buff, 19);
                     memcpy( phg_line + 19, "        H ", 10);
                     memcpy( phg_line + 29, buff + 23, 4);  /* move H */
                     memcpy( phg_line + 33, "   G ", 5);
                     memcpy( phg_line + 38, buff + 35, 5);  /* move G */
                     phg_line[43] = '\0';
                     }
                  else
                     strcpy( phg_line, buff);

                  fprintf( ofile, "%s", phg_line);
                  if( uncertainty_parameter < 90.)
                     fprintf( ofile, "   U%5.1f", uncertainty_parameter);
                  if( available_sigmas == 2)
                     fprintf( ofile, "  MC");
                  if( available_sigmas == 3)
                     fprintf( ofile, "  SR");
                  fprintf( ofile, "\n");
                  memmove( buff, zptr, strlen( zptr) + 1);
                  consider_replacing( buff, "q", "sigma_q");
                  strcat( buff, "    ");
                  strcat( buff, tbuff);
                  if( !reference_shown)
                     reference_shown = show_reference( buff);
                  }
               }
               break;
            default:
               break;
            }
         }
      assert( body_frame_note);
      if( !body_frame_note_shown)
         {
         size_t j = strlen( buff);

         if( j < 30)
            {
            while( j < 36)
               buff[j++] = ' ';
            strcpy( buff + j, body_frame_note);
            body_frame_note_shown = true;
            }
         else
            {
            j = 30;
            while( buff[j] == ' ')
               j++;
            if( j >= 60)    /* spaces to put the note in */
               {
               memcpy( buff + 36, body_frame_note, strlen( body_frame_note));
               body_frame_note_shown = true;
               }
            else if( i > 5 && (j = strlen( buff)) < 58)
               {
               memset( buff + j, ' ', 58 - j);
               strcpy( buff + 58, body_frame_note);
               body_frame_note_shown = true;
               }        /* above basically says,  "if we haven't gotten */
            }           /* the body frame note in the 'normal' places, try */
         }              /* wedging it in at the end" */
      fprintf( ofile, "%s\n", buff);
      tptr += strlen( tptr) + 1;
      }
   observation_summary_data( tbuff, obs, n_obs, options);
   fprintf( ofile, "%s\n", tbuff);
   if( elem.central_obj == 3 && elem.ecc < .99)
      {
      write_tle_from_vector( tbuff, rel_orbit, elem.epoch, NULL, NULL);
      tbuff[69] = tbuff[140] = '\0';
      fprintf( ofile, "# %s\n# %s\n", tbuff, tbuff + 71);
      }
   if( !(options & ELEM_OUT_NO_COMMENT_DATA))
      {
      double orb[6];
      const bool is_ecliptic =
                  (atoi( get_environment_ptr( "VECTOR_OPTS")) != 0);

      memcpy( orb, orbit2, 6 * sizeof( double));
      if( !is_ecliptic)
         {
         ecliptic_to_equatorial( orb);
         ecliptic_to_equatorial( orb + 3);
         ecliptic_to_equatorial( rel_orbit);
         ecliptic_to_equatorial( rel_orbit + 3);
         }
      fprintf( ofile, "# State vector (heliocentric %s J2000):\n",
               is_ecliptic ? "ecliptic" : "equatorial");
      fprintf( ofile, "# %+17.12f%+17.12f%+17.12f AU\n",
               orb[0], orb[1], orb[2]);
      fprintf( ofile, "# %+17.12f%+17.12f%+17.12f mAU/day\n",
               orb[3] * 1000., orb[4] * 1000., orb[5] * 1000.);
      if( planet_orbiting)
         {
         fprintf( ofile, "# State vector relative to central body:\n");
         fprintf( ofile, "# %17.12f%17.12f%17.12f AU\n",
               rel_orbit[0], rel_orbit[1], rel_orbit[2]);
         fprintf( ofile, "# %17.12f%17.12f%17.12f mAU/day\n",
               rel_orbit[3] * 1000., rel_orbit[4] * 1000., rel_orbit[5] * 1000.);
         }
      else   /* for heliocentric orbits,  show MOIDs: */
         {
         fprintf( ofile, "# MOIDs: Me%10.6f Ve%10.6f Ea%10.6f Ma%10.6f\n",
                  moids[1], moids[2], moids[3], moids[4]);
         fprintf( ofile, "# MOIDs: Ju%10.6f Sa%10.6f Ur%10.6f Ne%10.6f\n",
                  moids[5], moids[6], moids[7], moids[8]);
         }
      }
   if( monte_carlo)
      {
      if( !monte_carlo_object_count)
         n_clones_accepted = 0;
      if( rms_ok)
         n_clones_accepted++;
      }

   monte_carlo_permits = (n_clones_accepted == 1 ? "tfcwb" : "tfcab");
   if( monte_carlo && rms_ok)
      {
      FILE *ofile2 = fopen_ext( get_file_name( buff, "state.txt"), monte_carlo_permits);

      if( ofile2)
         {
         fseek( ofile, 0L, SEEK_SET);
         while( fgets( buff, sizeof( buff), ofile))
            fwrite( buff, strlen( buff), 1, ofile2);
         fclose( ofile2);
         fseek( ofile, 0L, SEEK_END);
         }
      }

   if( !(options & ELEM_OUT_NO_COMMENT_DATA))
      {
      const double jd = current_jd( );
      int jpl_de_version;

      full_ctime( buff, jd, CALENDAR_JULIAN_GREGORIAN);
      fprintf( ofile, "# Elements written: %s (JD %f)\n", buff, jd);
      make_date_range_text( buff, obs[0].jd, obs[n_obs - 1].jd);
      fprintf( ofile, "# Full range of obs: %s (%d observations)\n",
                              buff, n_obs);
      fprintf( ofile, "# Find_Orb ver: %s %s\n", __DATE__, __TIME__);
      fprintf( ofile, "# Perturbers: %08lx ", (unsigned long)perturbers);
      if( !perturbers)
         fprintf( ofile, "(unperturbed orbit)");
      else if( (perturbers & 0x3fe) == 0x3fe)
         fprintf( ofile, (perturbers & 0x400) ? "(Merc-Pluto plus Luna)" :
               "(Merc-Pluto, Earth & moon combined)");
      else if( perturbers == 0x408)
         fprintf( ofile, "(Sun/Earth/Moon)");
      get_jpl_ephemeris_info( &jpl_de_version, NULL, NULL);
      if( jpl_de_version)
         fprintf( ofile, ";  JPL DE-%d\n", jpl_de_version);
      else
         fprintf( ofile, ";  not using JPL DE\n");

      if( !elem.central_obj)
         {
         if( elem.ecc < 1.)
            {
            for( i = 0; i < 3; i++)  /* show Tisserand data for Ear, Jup & Nep */
               {                     /* if orbits come close to overlapping */
               const double semimajor_axes[3] = { 1., 5.2033, 30.069 };
               const char *names[3] = { "Earth", "Jupiter", "Neptune" };
               const double ratio =  semimajor_axes[i] / helio_elem.major_axis;
               const double tisserand = ratio
                  + 2. * sqrt( (1. - helio_elem.ecc * helio_elem.ecc) / ratio)
                  * cos( helio_elem.incl);
               const double aphelion = helio_elem.major_axis * 2. - helio_elem.q;

               if( helio_elem.q < semimajor_axes[i] / .7 && aphelion > semimajor_axes[i] * .7)
                  fprintf( ofile, "# Tisserand relative to %s: %.5f\n",
                        names[i], tisserand);
               }
            }
         if( helio_elem.q < 1.04)
            fprintf( ofile, "# Earth encounter velocity %.4f km/s\n",
                              encounter_velocity( &helio_elem, 1.));
         }
      if( !elem.central_obj || elem.central_obj == -1)
         if( elem.q < .15)  /* for both helio & barycentric orbits */
            {
            double ecliptic_lat, ecliptic_lon;

            get_periapsis_loc( &ecliptic_lon, &ecliptic_lat, &helio_elem);
            fprintf( ofile, "# Perihelion (%.3f, %.3f)\n",
                  ecliptic_lon * 180. / PI, ecliptic_lat * 180. / PI);
            }
      if( barbee_style_delta_v && helio_elem.q < 1.3)
         fprintf( ofile, "# Barbee-style encounter velocity: %.4f km/s\n",
                              barbee_style_delta_v);
      if( elem.abs_mag && elem.is_asteroid)
         {
         const double diam = diameter_from_abs_mag( elem.abs_mag, optical_albedo);

         fprintf( ofile, "# Diameter %.1f %s (assuming %.0f%% albedo)\n",
               (diam > 10000. ? diam / 1000. : diam),
               (diam > 10000. ? "km" : "meters"), optical_albedo * 100.);
         }

      fprintf( ofile, "# Score: %f\n", evaluate_initial_orbit( obs, n_obs, orbit));
      }

   *impact_buff = '\0';
   if( elem.central_obj < 15)
      {
      double latlon[2], t0;
      const int is_an_impact = (obs->jd < elem.perih_time);
                                         /* basically means,  "if we */
                                         /* observed the object after */
                                         /* periapsis, must be a launch; */
                                         /* otherwise,  must be impact." */
      ELEMENTS j2000_ecliptic_rel_elem = elem;

      calc_classical_elements( &j2000_ecliptic_rel_elem,
                                j2000_ecliptic_rel_orbit, epoch_shown, 1);

      t0 = find_collision_time( &j2000_ecliptic_rel_elem, latlon, is_an_impact);
      if( t0 < 1.)      /* t0 = 1 -> it was a miss after all */
         {
         char *end_ptr;
         const double lon = latlon[0] * 180. / PI;
         const double impact_time_td = elem.perih_time + t0;
         const double impact_time_utc = impact_time_td -
                        td_minus_utc( impact_time_td) / seconds_per_day;

         full_ctime( buff, impact_time_utc,
                       FULL_CTIME_HUNDREDTH_SEC | CALENDAR_JULIAN_GREGORIAN);
         snprintf( impact_buff, sizeof( impact_buff),
               " %.55s lat %+9.5f lon ", buff,
               latlon[1] * 180. / PI);
         end_ptr = impact_buff + strlen( impact_buff);
                     /* 0 < longitude < 360;  for Earth,  show this in */
                     /* "conventional" East/West 0-180 degree format:  */
         if( elem.central_obj == 3)
            {
            sprintf( end_ptr, "%c%.5f",
                  (lon < 180. ? 'E' : 'W'),
                  (lon < 180. ? lon : 360. - lon));
            fprintf( ofile, "%s at %s\n", (is_an_impact ? "IMPACT" : "LAUNCH"),
                                   impact_buff);
            }
                     /* Then show in 0-360 format,  for all other  */
                     /* planets, and for output to file:           */
         sprintf( end_ptr, "%9.5f", lon);
         if( elem.central_obj != 3)
            fprintf( ofile, "%s at %s\n", (is_an_impact ? "IMPACT" : "LAUNCH"),
                             impact_buff);
         }
      }
   if( *get_environment_ptr( "PLANET_STATES"))
      {
      fprintf( ofile, "# Planet states as of JD %.2f:\n", elem.epoch);
      for( i = 1; i < 11; i++)
         {
         double loc[3], vel[3];

         compute_observer_loc( elem.epoch, i, 0., 0., 0., loc);
         compute_observer_vel( elem.epoch, i, 0., 0., 0., vel);
         fprintf( ofile, "#%02d: %17.12f%17.12f%17.12f AU\n",
                   i, loc[0], loc[1], loc[2]);
         fprintf( ofile, "#    %17.12f%17.12f%17.12f mAU/day\n",
                   vel[0] * 1000., vel[1] * 1000., vel[2] * 1000.);
         }
      }
   if( !(options & ELEM_OUT_NO_COMMENT_DATA))
      {
      char time_buff[40];
      bool nongrav_sigmas_found = false;

      sprintf( tbuff, "#  $Name=%s", object_name);
      text_search_and_replace( tbuff + 4, " ", "%20");
               /* Epoch has to be given in YYYYMMDD.DDDDD format: */
      full_ctime( time_buff, helio_elem.perih_time,
               FULL_CTIME_YMD | FULL_CTIME_MONTHS_AS_DIGITS
               | FULL_CTIME_MICRODAYS | FULL_CTIME_LEADING_ZEROES);
      time_buff[4] = time_buff[7] = '\0';
      sprintf( tbuff + strlen( tbuff), "  $Ty=%s  $Tm=%s  $Td=%s",
               time_buff, time_buff + 5, time_buff + 8);
      sprintf( tbuff + strlen( tbuff), "  $MA=%.5f",
                  centralize_ang( helio_elem.mean_anomaly) * 180. / PI);
      fprintf( ofile, "%s\n", tbuff);

      sprintf( tbuff, "#  $ecc=%.7f  $Eqnx=2000.", helio_elem.ecc);
      fprintf( ofile, "%s\n", tbuff);

      sprintf( tbuff, "#  $a=%.7f  $Peri=%.5f  $Node=%.5f",
                  helio_elem.major_axis,
                  centralize_ang( helio_elem.arg_per) * 180. / PI,
                  centralize_ang( helio_elem.asc_node) * 180. / PI);
      sprintf( tbuff + strlen( tbuff), "  $Incl=%.5f",
                  helio_elem.incl * 180. / PI);
      fprintf( ofile, "%s\n", tbuff);

      sprintf( tbuff, "#  $EpJD=%.3f  $q=%.6f", helio_elem.epoch, helio_elem.q);
      sprintf( tbuff + strlen( tbuff), "  $T=%.6f  $H=%.1f",
               helio_elem.perih_time, helio_elem.abs_mag);
      fprintf( ofile, "%s\n", tbuff);
      fprintf( ofile, "# Sigmas avail: %d\n", available_sigmas);
      for( i = 1; i <= 3; i++)
         {
         char key[20], obuff[50];

         sprintf( key, "Sigma_A%d:", i);
         if( !get_uncertainty( key, obuff, false))
            {
            fprintf( ofile, (nongrav_sigmas_found ? "   %s %s" : "# %s %s"),
                                                   key, obuff);
            nongrav_sigmas_found = true;
            }
         }
      if( nongrav_sigmas_found)
         fprintf( ofile, "\n");
      }
   fclose( ofile);
         /* If the eccentricity is greater than about 1.2 for an heliocentric
            orbit,  and it's not marked as an interstellar object either by
            designation or a 'COM interstellar' line in the input,  the orbit
            is probably bogus and should be shown in,  say,  flashing red.  */

   if( !is_interstellar && !helio_elem.central_obj && helio_elem.ecc > 1.2)
      bad_elements = 1;
   if( helio_elem.q > 90.)
      bad_elements |= 2;

         /* Also,  write out elements in MPCORB-like format: */
   if( helio_elem.ecc < .999999)
      {
      const char *output_filename = (*elements_filename == 's' ?
                                             "mpc_sr.txt" : mpc_fmt_filename);

      ofile = fopen_ext( get_file_name( tbuff, output_filename), file_permits);
      elements_in_mpcorb_format( tbuff, obs->packed_id, object_name,
                              &helio_elem, obs, n_obs);
      fprintf( ofile, "%s\n", tbuff);
      fclose( ofile);
      }

   if( monte_carlo)
      monte_carlo_object_count++;
   if( monte_carlo && rms_ok)
      {
      const char *element_filename = get_environment_ptr( "MONTE_CARLO");
      char name_buff[48], virtual_full_desig[40];

      if( !*element_filename)
         element_filename = "mpcorb.dat";
      if( *impact_buff)
         n_monte_carlo_impactors++;
      sprintf( name_buff, "%05d", n_clones_accepted);
      packed_desig_minus_spaces( virtual_full_desig, obs->packed_id);
      sprintf( virtual_full_desig + strlen( virtual_full_desig), " [%d]",
                                  monte_carlo_object_count);
      if( elem.central_obj || elem.ecc > .999999)
         {
         ofile = fopen_ext( get_file_name( tbuff, "virtual.txt"), monte_carlo_permits);
         elements_in_guide_format( tbuff, &elem, virtual_full_desig, obs, n_obs);
         fprintf( ofile, "%s%s\n", tbuff, impact_buff);
         fclose( ofile);
         }

      if( helio_elem.ecc < .999999)
         {
         ofile = fopen_ext( get_file_name( tbuff, element_filename), monte_carlo_permits);
         if( !strcmp( monte_carlo_permits, "wb"))
            {        /* new file = write out a header for it */
            FILE *ifile = fopen_ext( "mpcorb.hdr", "fcrb");
            time_t t0 = time( NULL);

            fprintf( ofile, "Monte Carlo orbits from Find_Orb\nComputed %s", ctime( &t0));
            fprintf( ofile, "Find_Orb version %s %s\n", __DATE__, __TIME__);
            fprintf( ofile, (using_sr ? "Statistical Ranging\n" : "Full Monte Carlo\n"));
            if( ifile)
               {
               while( fgets( tbuff, 80 * 9, ifile))
                  fputs( tbuff, ofile);
               fclose( ifile);
               }
            }
         elements_in_mpcorb_format( tbuff, name_buff, virtual_full_desig,
                              &helio_elem, obs, n_obs);
         fprintf( ofile, "%s%s\n", tbuff, impact_buff);
         fclose( ofile);
         }
//    if( helio_elem.ecc > 2. || helio_elem.q > 90.)
//       set_statistical_ranging( 1);
      }

   if( (ofile = fopen_ext( get_file_name( tbuff, "guide.txt"), "tfcwb")) != NULL)
      {
      elements_in_guide_format( tbuff, &elem, object_name, obs, n_obs);
      fprintf( ofile, "%s%s\n", tbuff, impact_buff);
      fclose( ofile);
      }

   ofile = open_json_file( tbuff, "JSON_ELEMENTS_NAME", "elements.json", obs->packed_id,
                     "wb");
   assert( ofile);
   elements_in_json_format( ofile, &elem, object_name, obs, n_obs, moids,
                                             body_frame_note, true);
   fclose( ofile);
   ofile = open_json_file( tbuff, "JSON_SHORT_ELEMENTS", "elem_short.json", obs->packed_id,
                     "wb");
   assert( ofile);
   elements_in_json_format( ofile, &elem, object_name, obs, n_obs, moids,
                                             body_frame_note, false);
   fclose( ofile);
   make_linkage_json( n_obs, obs, &elem);
   free( tbuff);
   return( bad_elements);
}

void shellsort_r( void *base, const size_t n_elements, const size_t esize,
         int (*compare)(const void *, const void *, void *), void *context);

int string_compare_for_sort( const void *a, const void *b, void *context)
{
   const char **a1 = (const char **)a;
   const char **b1 = (const char **)b;
   const int *sort_column = (int *)context;

   if( *sort_column == -1)
      return( names_compare( a1[0], b1[0]));
   return( strcmp( a1[0] + *sort_column, b1[0] + *sort_column));
}

static const char *vector_data_file = "vectors.dat";

/* For each object,  we'd like to know if a solution has already been
computed for it and stored as a state vector in 'vectors.dat'.  We
do this by digging through the file and looking for lines with a
space,  epoch,  and object name,  such as

 2458127.5 2020 CD3

   We sort those names,  then go through each OBJECT_INFO structure
and binary-search through the sorted list to see if a match can be
found.  If one is,  then a 'stored solution' exists in 'vectors.dat'
for that object.   */

void set_solutions_found( OBJECT_INFO *ids, const int n_ids)
{
   size_t n_lines;
   char **ilines = load_file_into_memory( vector_data_file, &n_lines, false);
   int i, n_objs, sort_column = -1;

   for( i = 0; i < n_ids; i++)
      ids[i].solution_exists = 0;
   if( !ilines)
      return;
   for( i = n_objs = 0; (size_t)i < n_lines; i++)
      if( ilines[i][0] == ' ' && isdigit( ilines[i][1]))
         {
         const char *space = strchr( ilines[i] + 1, ' ');

         assert( space);     /* object name starts after a space */
         memmove( ilines[i], space + 1, strlen( space));
         ilines[n_objs++] = ilines[i];
         }
   shellsort_r( ilines, n_objs, sizeof( char *),
                           string_compare_for_sort, &sort_column);
   for( i = 0; i < n_ids; i++)
      {
      size_t loc = 0, loc1, step;

      for( step = 0x80000000; step; step >>= 1)
         if( (loc1 = loc + step) < (size_t)n_objs)
            {
            const int compare = names_compare( ilines[loc1], ids[i].obj_name);

            if( compare < 0)
               loc = loc1;
            if( !compare)
               ids[i].solution_exists = 1;
            }
      }
   free( ilines);
}

/* Can get some comet elements from

   http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET

   This is especially helpful for pre-2008 SOHO objects.  Without
ELEMENTS.COMET,  Find_Orb can sometimes flounder about a bit in its
efforts to determine an orbit. */

int extract_sof_data_ex( ELEMENTS *elem, const char *buff, const char *header,
                        double *extra_info);

static int get_orbit_from_mpcorb_sof( const char *filename,
                 const char *object_name, double *orbit, ELEMENTS *elems,
                 const double full_arc_len, double *max_resid)

{
   FILE *ifile = fopen_ext( filename, "crb");
   int got_vectors = 0;

   memset( elems, 0, sizeof( ELEMENTS));
   if( ifile)
      {
      char buff[300], header[300], tname[15];

      if( !fgets_trimmed( header, sizeof( header), ifile))
         {
         fprintf( stderr, "Error in %s header\n", filename);
         exit( -1);
         }
      while( *object_name == ' ')
         object_name++;
      if( *object_name == '(')         /* numbered asteroid */
         snprintf( tname, sizeof( tname), "%12d", atoi( object_name + 1));
      else if( *object_name == 'P' && object_name[1] == '/'
                  && atoi( object_name + 2) < 1000)      /* numbered comet */
         snprintf( tname, sizeof( tname), "%3uP        ",
                  atoi( object_name + 2));
      else
         snprintf( tname, sizeof( tname), "%-12.12s", object_name);
      while( !got_vectors && fgets_trimmed( buff, sizeof( buff), ifile))
         if( !memcmp( tname, buff, 12))
            {
            double extra_info[10];

            extract_sof_data_ex( elems, buff, header, extra_info);
            if( elems->epoch < 2400000.)
               printf( "JD %f\n", elems->epoch);
            assert( elems->epoch > 2400000.);
            if( extra_info[1] - extra_info[0] > full_arc_len / 2.)
               {
               derive_quantities( elems, get_planet_mass( elems->central_obj));
               comet_posn_and_vel( elems, elems->epoch, orbit, orbit + 3);
               if( elems->central_obj == 3)
                  {
                  equatorial_to_ecliptic( orbit);
                  equatorial_to_ecliptic( orbit + 3);
                  }
               if( elems->central_obj)
                  {
                  double planet_state[6];
                  int i;

                  get_planet_posn_vel( elems->epoch, elems->central_obj,
                                          planet_state, planet_state + 3);
                  for( i = 0; i < 6; i++)
                     orbit[i] += planet_state[i];
                  }
               got_vectors = 1;
               if( elems->ecc == 1.)     /* indicate parabolic-constraint orbit */
                  got_vectors = 2;
               *max_resid = extra_info[3];
               }
            }
      fclose( ifile);
      }
   return( got_vectors);
}

/* The following ensures that names starting with the same international
artsat designations compare as equal,  even if they diverge in irrelevant
ways after that.  For example,  2013-024B = NORAD 39169 = WGS 5 Rk would
compare as equal to 2013-024B or with 2013-024B = NORAD 39169.  */

static int names_compare( const char *name1, const char *name2)
{
   if( name1[4] == '-' && name1[9] == ' ')
      {
      unsigned mask = 0, i;

      for( i = 0; i < 9 && name1[i]; i++)
         if( name1[i] >= '0' && name1[i] <= '9')
            mask |= (1u << i);
      if( mask == 0xef)
         if( !memcmp( name1, name2, 9))
            return( 0);
      }
   return( strcmp( name1, name2));
}

#if !defined( _WIN32)
int _memicmp( const char *s1, const char *s2, int n)
{
   int c1, c2;

   while( n--)
      {
      if( (c1 = tolower( *s1++)) != (c2 = tolower( *s2++)))
         return( c1 - c2);
      }
   return( 0);
}

int _stricmp( const char *s1, const char *s2)
{
   int c1, c2;

   while( *s1 || *s2)
      {
      if( (c1 = tolower( *s1++)) != (c2 = tolower( *s2++)))
         return( c1 - c2);
      }
   return( 0);
}
#endif

/* You can specify a state vector on the command line (for fo and
Find_Orb) using the -v(state vect) switch.  The state vector must
have,  at minimum,  an epoch and six numbers (position and velocity).
Further options can be added after those numbers.  For example,

fo -oMade-up -v2020jan13,1.2,2.3,3.4,-.005,.002,.001,H=17,eq

would make a dummy object 'Made-up',  at epoch 2020 Jan 13,
position 1.2, 2.3, 3.4 (AU),  velocity -.005,-.002,.001 (AU/day),
with an absolute mag of 17,  in heliocentric J2000 equatorial
coordinates.  You can use km and seconds as units,  by adding
'km' and/or 's'; or meters with 'm';  or 'g' for geocentric.
'ep=2020.9' would set the coordinate plane epoch as that of
2020.9.  'UTC' would specify that the epoch is in UTC (some of
the artsat crowd give epochs on that time scale). 'TDB' would
specify a TDB epoch.

   To do : revise to allow barycentric and other-planet-centric
(and selenocentric) vectors,  and vectors in the body plane. */

static double extract_state_vect_from_text( const char *text,
            double *orbit, double *solar_pressure, double *abs_mag)
{
   char tbuff[81];
   size_t i = 0;
   double epoch = 0., coord_epoch = 2000.;
   int bytes_read, central_object = 0;

   while( text[i] && text[i] != ',' && i < 80)
      i++;
   if( text[i] == ',')
      {
      memcpy( tbuff, text, i);
      tbuff[i] = '\0';
      epoch = get_time_from_string( 0., tbuff, 0, NULL);
      }
   if( epoch)
      {
      text += i + 1;
      for( i = 0; i < 6 && sscanf( text, "%lf%n", orbit + i, &bytes_read) == 1; i++)
         text += bytes_read + 1;
      if( i != 6)       /* didn't read a whole state vect */
         epoch = 0.;
      }
   assert( epoch);
   *abs_mag = 18.;
   while( *text)
      {
      if( !memcmp( text, "eq", 2))
         {
         equatorial_to_ecliptic( orbit);
         equatorial_to_ecliptic( orbit + 3);
         text += 2;
         }
      else if( !memcmp( text, "ep=", 3))
         coord_epoch = atof( text + 3);
      else if( !memcmp( text, "H=", 2))
         *abs_mag = atof( text + 2);
      else if( !_memicmp( text, "UTC", 3))
         epoch += td_minus_utc( epoch) / seconds_per_day;
      else if( !_memicmp( text, "TDB", 3))
         {
         const double t_cen = (epoch - J2000) / 36525.;

         epoch -= tdb_minus_tdt( t_cen) / seconds_per_day;
         }
      else if( !memcmp( text, "km", 2))
         for( i = 0; i < 6; i++)
            orbit[i] /= AU_IN_KM;
      else if( *text == 'm')
         for( i = 0; i < 6; i++)
            orbit[i] /= AU_IN_METERS;
      else if( *text == 'g')
         central_object = 3;
      else if( *text == 's')
         for( i = 3; i < 6; i++)
            orbit[i] *= seconds_per_day;
      else if( *text == 'A' && text[1] >= '1' && text[1] <= '3'
                  && text[2] == '=')
         {
         extern int n_extra_params;

         n_extra_params = text[1] - '1';
         force_model = (n_extra_params | 0x10);   /* assume inverse square for the nonce */
         solar_pressure[n_extra_params - 1] = atof( text + 3);
         }
      while( *text != ',' && *text)
         text++;
      if( *text == ',')
         text++;
      }
   if( coord_epoch != 2000.)
      {
      double precess_matrix[9];

      setup_ecliptic_precession( precess_matrix, coord_epoch, 2000.);
      precess_vector( precess_matrix, orbit, orbit);
      precess_vector( precess_matrix, orbit + 3, orbit + 3);
      }
   if( central_object)
      {
      double vect[6];

      compute_observer_loc( epoch, central_object, 0., 0., 0., vect);
      compute_observer_vel( epoch, central_object, 0., 0., 0., vect + 3);
      for( i = 0; i < 6; i++)
         orbit[i] += vect[i];
      }
   return( epoch);
}

const char *state_vect_text = NULL;
int ignore_prev_solns;
bool take_first_soln = false, force_final_full_improvement = false;

static int fetch_previous_solution( OBSERVE *obs, const int n_obs, double *orbit,
               double *orbit_epoch, unsigned *perturbers)
{
   FILE *ifile = NULL;
   int got_vectors = 0, i;
   extern int n_extra_params;
   extern double solar_pressure[];
   char object_name[80];
   bool do_full_improvement = false;
   double abs_mag = 10.;         /* default value for 'dummy' use */

   get_object_name( object_name, obs->packed_id);
   for( i = 0; i < MAX_N_NONGRAV_PARAMS; i++)
      solar_pressure[i] = 0.;
   n_extra_params = 0;
   if( state_vect_text)    /* state vect supplied on cmd line w/ '-v' */
      {
      *orbit_epoch = extract_state_vect_from_text(
                  state_vect_text, orbit, solar_pressure, &abs_mag);
      got_vectors = (*orbit_epoch != 0.);
      }
   if( !got_vectors)
      ifile = (ignore_prev_solns ? NULL : fopen_ext( vector_data_file, "crb"));
   if( ifile)
      {
      char buff[120];
      double jd1 = 0., jd2 = 0.;
      double residual_filter_threshhold = 0.;

      while( fgets_trimmed( buff, sizeof( buff), ifile))
         {
         i = 1;
         while( buff[i] && buff[i] != ' ')
            i++;
         assert( buff[i] == ' ');
         if( !names_compare( object_name, buff + i + 1) || (take_first_soln && !got_vectors))
            {
            int n_read, center = 0, is_equatorial = 0;
            char *tptr;
            double dist_units = 1., time_units = .001;

            got_vectors = 1;
            *orbit_epoch = atof( buff);
            *perturbers = 0;
            fgets_trimmed( buff, sizeof( buff), ifile);
            tptr = strchr( buff, '#');
            if( tptr)
               {
               char *tptr2;

               *tptr++ = '\0';
               tptr2 = strstr( tptr, "Ctr");
               if( tptr2)
                  center = atoi( tptr2 + 3);
               if( strstr( tptr, "km"))
                  dist_units = AU_IN_KM;
               if( strstr( tptr, "eq"))
                  is_equatorial = 1;
               if( strstr( tptr, "sec"))
                  time_units = seconds_per_day;
               }
            for( i = 0; i < MAX_N_NONGRAV_PARAMS; i++)
                solar_pressure[i] = 0.;
            n_read = sscanf( buff, "%lf%lf%lf%x%lf %lf %lf",
                          &orbit[0], &orbit[1], &orbit[2], perturbers,
                          solar_pressure,
                          solar_pressure + 1,
                          solar_pressure + 2);
            assert( n_read >= 3 && n_read < 8);
            n_extra_params = n_read - 4;
            if( n_extra_params <= 0)
               n_extra_params = force_model = 0;
            else
               force_model = (n_extra_params | 0x10);  /* assume inv square for the nonce */
            fgets_trimmed( buff, sizeof( buff), ifile);
            sscanf( buff, "%lf%lf%lf%lf%lf%lf",
                          orbit + 3, orbit + 4, orbit + 5,
                          &residual_filter_threshhold, &jd1, &jd2);
            for( i = 3; i < 6; i++)
               orbit[i] /= dist_units / time_units;
            for( i = 0; i < 3; i++)
               orbit[i] /= dist_units;
            if( is_equatorial)
               {
               equatorial_to_ecliptic( orbit);
               equatorial_to_ecliptic( orbit + 3);
               }
            if( center)
               {
               double vect[3];

               compute_observer_loc( *orbit_epoch, center, 0., 0., 0., vect);
               for( i = 0; i < 3; i++)
                  orbit[i] += vect[i];
               compute_observer_vel( *orbit_epoch, center, 0., 0., 0., vect);
               for( i = 0; i < 3; i++)
                  orbit[i + 3] += vect[i];
               }
            for( i = 0; i < n_obs; i++)
               {
               obs[i].computed_ra  = obs[i].ra;
               obs[i].computed_dec = obs[i].dec;
               }
            }
         }
      if( got_vectors)
         {
         set_locs( orbit, *orbit_epoch, obs, n_obs);
         if( jd2)
            {
            for( i = 0; i < n_obs; i++)
               if( obs[i].jd < jd1 - .00001 || obs[i].jd > jd2 + .00001)
                  obs[i].is_included = 0;
               else
                  if( residual_filter_threshhold &&
                        observation_rms( obs + i) > residual_filter_threshhold)
                     obs[i].is_included = 0;
            }
         do_full_improvement = true;
         }
      fclose( ifile);
      }
   if( !got_vectors && !ignore_prev_solns)
      {
      ELEMENTS elems;
      const char *mpcorb_dot_sof_filename = get_environment_ptr( "MPCORB_SOF_FILENAME");
      const double full_arc_len = obs[n_obs - 1].jd - obs[0].jd;
      double rms_resid;

      got_vectors = get_orbit_from_mpcorb_sof( "orbits.sof", obs->packed_id,
                                          orbit, &elems, full_arc_len, &rms_resid);
      if( !got_vectors && *mpcorb_dot_sof_filename)
         got_vectors = get_orbit_from_mpcorb_sof( mpcorb_dot_sof_filename,
                             object_name, orbit, &elems, full_arc_len, &rms_resid);
      if( got_vectors)
         {
         *orbit_epoch = elems.epoch;
         if( got_vectors == 1)
            *perturbers = 0x7fe;    /* Merc-Pluto plus moon */
         set_locs( orbit, *orbit_epoch, obs, n_obs);
         filter_obs( obs, n_obs, rms_resid * 3., 1);
         abs_mag = elems.abs_mag;
         do_full_improvement = true;
         }
      }
   if( !got_vectors)
      {
      *perturbers = 0;
      *orbit_epoch = initial_orbit( obs, n_obs, orbit);
      if( force_final_full_improvement)
         do_full_improvement = true;
      }
   else if( n_obs == 1 && !strcmp( obs->reference, "Dummy"))
      {
      obs->jd = *orbit_epoch + td_minus_ut( *orbit_epoch) / seconds_per_day;
      *perturbers = 0x7fe;    /* Merc-Pluto plus moon */
      set_locs( orbit, *orbit_epoch, obs, n_obs);
      obs->ra = obs->computed_ra;
      obs->dec = obs->computed_dec;
      obs->obs_mag = abs_mag + calc_obs_magnitude( obs->solar_r,
                 obs->r, vector3_length( obs->obs_posn), NULL);
      obs->obs_mag = floor( obs->obs_mag * 10. + .5) * .1;
      obs->mag_precision = 1;         /* start out assuming mag to .1 mag */
      }
   if( do_full_improvement)
      {
      OBSERVE *saved_obs = (OBSERVE *)calloc( n_obs, sizeof( OBSERVE));
      double prev_score;

      if( !*get_environment_ptr( "KEEP_PREVIOUS_EPOCH"))
         {
         const double new_epoch = find_epoch_shown( obs, n_obs);

         integrate_orbit( orbit, *orbit_epoch, new_epoch);
         *orbit_epoch = new_epoch;
         }
      memcpy( saved_obs, obs, n_obs * sizeof( OBSERVE));
      push_orbit( *orbit_epoch, orbit);
      prev_score = evaluate_initial_orbit( obs, n_obs, orbit);
      full_improvement( obs, n_obs, orbit, *orbit_epoch, NULL,
                           ORBIT_SIGMAS_REQUESTED, *orbit_epoch);
      if( prev_score < evaluate_initial_orbit( obs, n_obs, orbit) - .001)
         {
         pop_orbit( orbit_epoch, orbit);    /* we were better off with the old orbit */
         memcpy( obs, saved_obs, n_obs * sizeof( OBSERVE));
         }
      else        /* throw out the saved orbit;  we've got something better */
         pop_orbit( NULL, NULL);
      free( saved_obs);
      }
   return( got_vectors);
}

double override_epoch_shown = 0.;

double find_epoch_shown( const OBSERVE *obs, const int n_obs)
{
   int first, last;
   double rval;

   if( override_epoch_shown)
      return( override_epoch_shown);
   get_first_and_last_included_obs( obs, n_obs, &first, &last);
   if( !last)      /* no observations included;  shouldn't actually happen */
      rval = floor( (obs[first].jd + obs[last].jd) / 2.) + .5;
   else
      {
      rval = floor( obs[last].jd) + .5;
      if( rval > obs[last].jd)
         if( rval - obs[last].jd > obs[first].jd - (rval - 1.))
            rval--;
      }
   return( rval);
}

   /* When solving orbits for objects orbiting very near to Jupiter and
      Saturn,  Find_Orb will default to including perturbations from the
      satellites of those planets.  That can be a problem if the object
      in question _is_ a satellite of that planet;  when computing an
      orbit for,  say,  Callisto,  you don't really want to include
      perturbations from Callisto.  Only the other three Galileans ought
      to be considered in that case.  (If you _do_ attempt to include
      "perturbations by Callisto on Callisto",  you ought to get a
      divide-by-zero error.  In practice,  you just get wacky results.)

      The following function looks at a packed designation and tells
      you which perturber it is.  It returns -1 in the (usual) case
      that you're not solving for the orbit of an inner satellite that
      also happens to be a perturbing body in Find_Orb.        */

static int obj_desig_to_perturber( const char *packed_desig)
{
   int rval = -1, i;
   extern unsigned excluded_perturbers;
   extern double object_mass;
   const char *outer_planets[11] = { "            ",
                                     "     Mercury",
                                     "     Venus  ",
                                     "     Earth  ",
                                     "     Mars   ",
                                     "     Jupiter",
                                     "     Saturn ",
                                     "     Uranus ",
                                     "     Neptune",
                                     "     Pluto  ", NULL };

                     /* The EXCLUDED environment variable provides a way */
                     /* to specifically exclude some perturbers,  either */
                     /* to speed things up or for debugging purposes.    */
   sscanf( get_environment_ptr( "EXCLUDED"), "%x", &excluded_perturbers);
   if( excluded_perturbers == (unsigned)-1)
      excluded_perturbers = 0;
// excluded_perturbers |= 512;      /* Pluto is _always_ excluded */

   if( !memcmp( packed_desig + 4, "S    ", 5) && packed_desig[1] == '0'
                                              && packed_desig[2] == '0')
      {
      if( *packed_desig == 'J')
         {
         if( packed_desig[3] >= '1' && packed_desig[3] <= '4')
            rval = 11 + packed_desig[3] - '1';  /* Io..Callisto = 11..14 */
         }
      else if( *packed_desig == 'S')
         {
         if( packed_desig[3] >= '3' && packed_desig[3] <= '6')
            rval = 15 + packed_desig[3] - '3';  /* Tethys...Titan = 15..18 */
         else if( packed_desig[3] == '8')
            rval = 19;                          /* Iapetus */
         }
      else if( *packed_desig == 'E' && packed_desig[3] == '1')
         rval = 10;           /* Earth's moon */
      }
   if( strstr( packed_desig, "D4340"))    /* Pluto as an asteroid */
      rval = 9;
   if( rval > 0)
      excluded_perturbers |= (1 << rval);
   object_mass = 0.;
   for( i = 1; outer_planets[i]; i++)
      if( !strcmp( packed_desig, outer_planets[i]))
         {
         rval = i;
         excluded_perturbers |= (1 << rval);
         if( rval == 6)   /* also exclude five Saturn sats */
            excluded_perturbers |= (0x1f << 15);
         if( rval == 5)   /* also exclude four (Galilean) Jupiter sats */
            excluded_perturbers |= (0xf << 11);
         object_mass = get_planet_mass( rval) / SOLAR_GM;
         debug_printf( "Excluded: %d (%x)\n", rval, excluded_perturbers);
         }
   return( rval);
}

OBSERVE FAR *load_object( FILE *ifile, OBJECT_INFO *id,
                       double *curr_epoch, double *epoch_shown, double *orbit)
{
   extern int n_obs_actually_loaded;
   extern int debug_level;
   extern unsigned perturbers;
   extern int excluded_asteroid_number;
   OBSERVE FAR *obs = load_observations( ifile, id->packed_desig,
                                                id->n_obs);

   if( debug_level || n_obs_actually_loaded != id->n_obs)
      {
      char buff[80];

      debug_printf( " %d observations loaded\n", n_obs_actually_loaded);
      make_date_range_text( buff, obs[0].jd,
                                    obs[n_obs_actually_loaded - 1].jd);
      id->n_obs = n_obs_actually_loaded;
      debug_printf( "%s\n", buff);
      }
   obj_desig_to_perturber( id->packed_desig);
   if( id->obj_name[0] == '(')    /* numbered asteroid:  shouldn't perturb */
      excluded_asteroid_number = atoi( id->obj_name + 1);      /* itself */
   else
      excluded_asteroid_number = 0;
   if( n_obs_actually_loaded > 0)
      {
      const int got_vector = fetch_previous_solution( obs, id->n_obs, orbit,
                            curr_epoch, &perturbers);

      *epoch_shown = *curr_epoch;
      if( got_vector <= 0)
         *epoch_shown = find_epoch_shown( obs, id->n_obs);
      }
   else                           /* indicate failure */
      *epoch_shown = *curr_epoch = 0.;
   return( obs);
}

double observation_rms( const OBSERVE FAR *obs);            /* elem_out.cpp */

int store_solution( const OBSERVE *obs, const int n_obs, const double *orbit,
       const double orbit_epoch, const int perturbers)
{
   FILE *ofile = fopen_ext( vector_data_file, "fcab");

   if( ofile)
      {
      char buff[80];
      int i, j, k;
      double max_resid_included_obs = 0.;

      get_object_name( buff, obs->packed_id);
      fprintf( ofile, "%10.1f %s\n", orbit_epoch, buff);
      for( i = 0; i < 3; i++)
         fprintf( ofile, "%21.16f", orbit[i]);
      if( perturbers)
         {
         extern int n_extra_params;

         fprintf( ofile, " %04x", perturbers);
         if( n_extra_params)
            {
            extern double solar_pressure[];

            for( i = 0; i < n_extra_params; i++)
               fprintf( ofile, " %.9g", solar_pressure[i]);
            }
         }
      get_first_and_last_included_obs( obs, n_obs, &i, &j);
      for( k = i; k <= j; k++)
         {
         const double resid = observation_rms( obs + k);

         if( obs[k].is_included && max_resid_included_obs < resid)
            max_resid_included_obs = resid;
         }
      fprintf( ofile, "\n%21.16f%21.16f%21.16f %.3f %.5f %.5f\n",
              orbit[3] * 1000., orbit[4] * 1000., orbit[5] * 1000.,
              max_resid_included_obs + .001, obs[i].jd, obs[j].jd);
      fclose( ofile);
      }
   return( ofile ? 0 : -1);
}


#define LOG_10 2.3025850929940456840179914546843642076011014886287729760333279009675726

double calc_obs_magnitude( const double obj_sun,
            const double obj_earth, const double earth_sun, double *phase_ang)
{
   double rval;
   double ph_ang = obj_sun * obj_sun +
                  obj_earth * obj_earth - earth_sun * earth_sun;

   ph_ang /= 2. * obj_earth * obj_sun;
   ph_ang = acose( ph_ang);
   if( phase_ang)
      *phase_ang = ph_ang;

   if( object_type == OBJECT_TYPE_COMET)
      rval = comet_magnitude_slope_param * log( obj_sun);
   else
      {
      double phi1, phi2, log_tan_half_phase;

      log_tan_half_phase = log( sin( ph_ang / 2.) / cos( ph_ang / 2.));
      phi1 = exp( -3.33 * exp( log_tan_half_phase * 0.63));
      phi2 = exp( -1.87 * exp( log_tan_half_phase * 1.22));
      rval = 5. * log( obj_sun) - 2.5 *
                  log( (1. - asteroid_magnitude_slope_param) * phi1
                + asteroid_magnitude_slope_param * phi2);
      }
   rval += 5. * log( obj_earth);
   rval /= LOG_10;         /* allow for common logs,  not naturals */
   return( rval);
}

/* The following function,  as the comment indicates,  assumes that
a "no band" case (obs->mag_band == ' ') must be an R mag.  That's
probably the best guess for most modern CCD observations.  MPC
assumes a B (photographic) magnitude,  which is probably the best
guess for older observations.  I suppose one would ideally look at
the second 'note'  which is C for CCD observations and P for
photographic observations.  The code could then assume a default of
R for CCD obs,  B for photographic,  and V for the admittedly rare
micrometer or encoder-based observations.

It may be more accurate to infer a mag_band code based on the catalog
being used:  'r' for any of the UCACs,  'G' for Gaia,  'V' for GSC-1.x
for declinations north of +3 and 'R' south of that,  etc.  The idea
is that if you used UCAC magnitudes for calibration,  you're reporting
an 'r' magnitude,  no matter what filter you had.  Problems are cases
where the photometry was done with a different catalog than the astrometry;
or the catalog had multiple magnitude bands (Ax.0,  B1.0);  or tricky
things were done to "adjust" the magnitudes to a different band,  using
either 2MASS colors to estimate a color correction or just ass uming
an average color correction.

'C' ('clear',  'unfiltered') magnitudes are something of a problem.
They never really should have existed.  As described above,  the key
thing isn't what filter you use (though that's a nice thing to
know);  it's what sort of reference photometry you used.  At
present,   I'm treating 'C' magnitudes as close enough to 'R' or
'blank'.  (I perhaps should treat them as 'probably wrong' and give
them minimal or no weight in the computation of H values.)

Estimates for the mag band shifts in R, B, I, and U are from a post by
Petr Pravec,  https://groups.io/g/mpml/message/24598
in turn derived from data in Shevchenko and Lupishko, Solar System
Research 32, 220-232, 1998 :

http://www.asu.cas.cz/~ppravec/ShevchenkoandLupishko1998.pdf

grizyw (Sloan magnitudes) come from table 1, p. 360, 2013 April,
Publications of the Astro Soc of Pacific, 2013 PASP, 125:357-395,
Denneau et. al., 'The Pan-STARRS Moving Object Process System'

2017 Jul 11:  David Tholen found V-G=0.28 (private communication at
present,  will give a "proper" reference once it's published...
also provides V-G as quadratic functions of V-I or V-R,  which
should be handy in other contexts) */

double mag_band_shift( const char mag_band)
{
   double rval = 0.;
   const char *bands = "R BIUCgrizywG";
   const double offsets[] = { .43, .43,
                     /* R and 'no band' are treated alike */
         -.77, .82, -1.16, .4,    /* B, I, U, C */
         -0.28, 0.23, 0.39, 0.37, 0.36, 0.16,     /* grizyw */
         0.28 };                                  /* G */
   const char *tptr = strchr( bands, mag_band);

   if( tptr)
      rval = offsets[tptr - bands];
   return( rval);
}

double calc_absolute_magnitude( OBSERVE FAR *obs, int n_obs)
{
   int obs_no;
   double n_mags = 0.;
   double rval = 0.;

   for( obs_no = 0; obs_no < n_obs; obs_no++)
      {
      obs->computed_mag = 0.;
      if( obs->r && obs->solar_r)
         {
         const double earth_sun = vector3_length( obs->obs_posn);

         if( earth_sun)
            {
            bool use_obs = true;

            if( object_type == OBJECT_TYPE_COMET
                            && obs->mag_band != default_comet_magnitude_type)
               use_obs = false;
            if( obs->obs_mag == BLANK_MAG || !obs->is_included)
               use_obs = false;
            obs->computed_mag = calc_obs_magnitude(
                  obs->solar_r, obs->r, earth_sun, NULL) - mag_band_shift( obs->mag_band);
            if( use_obs)
               {
               rval += (obs->obs_mag - obs->computed_mag) / obs->mag_sigma;
               n_mags += 1. / obs->mag_sigma;
               }
            }
         }
      obs++;
      }
   if( n_mags)
      rval /= n_mags;
   obs -= n_obs;
   for( obs_no = 0; obs_no < n_obs; obs_no++, obs++)
      if( n_mags)
         obs->computed_mag += rval;
      else
         obs->computed_mag = 0.;
   return( rval);
}

int find_worst_observation( const OBSERVE FAR *obs, const int n_obs)
{
   int i, rval = -1;
   double worst_rms = 0., rms;

   for( i = 0; i < n_obs; i++, obs++)
      if( obs->is_included)
         {
         rms = compute_rms( obs, 1);
         if( rms > worst_rms)
            {
            worst_rms = rms;
            rval = i;
            }
         }
   return( rval);
}

/* If you've got n_obs observations stored in the obs array,  the
   get_idx1_and_idx2( ) function will puzzle through them to find the first
   and last valid observation (those that haven't had their 'is_included'
   flags set to FALSE),  and will store the indices to them in *idx1 and
   *idx2.  These are shown near the top of the display,  and are used in
   the method of Herget.  Return value is the number of included obs.   */

int get_idx1_and_idx2( const int n_obs, const OBSERVE FAR *obs,
                                          int *idx1, int *idx2)
{
   int i, rval = 0;

   for( i = 0; i < n_obs && (!obs[i].is_included || !obs[i].r); i++)
      ;
   if( i == n_obs)
      *idx1 = *idx2 = 0;
   else
      {
      *idx1 = i;
      for( ; i < n_obs; i++)
         if( obs[i].is_included && obs[i].r)
            rval++;
      for( i = n_obs - 1; i && (!obs[i].is_included || !obs[i].r); i--)
         ;
      *idx2 = i;
      }
   return( rval);
}

int get_r1_and_r2( const int n_obs, const OBSERVE FAR *obs,
                           double *r1, double *r2)
{
   int idx1, idx2, rval = get_idx1_and_idx2( n_obs, obs, &idx1, &idx2);

   if( !rval)
      *r1 = *r2 = 0.;
   else
      {
      *r1 = obs[idx1].r;
      *r2 = obs[idx2].r;
      }
   return( rval);
}

extern char default_comet_magnitude_type;
extern double probability_of_blunder;
extern double overobserving_time_span;
extern unsigned overobserving_ceiling;
extern int use_blunder_method;
extern bool use_sigmas;
extern double ephemeris_mag_limit;
extern int apply_debiasing;
extern int sigmas_in_columns_57_to_65;

int store_defaults( const ephem_option_t ephemeris_output_options,
         const int element_format, const int element_precision,
         const double max_residual_for_filtering,
         const double noise_in_arcseconds)
{
   char buff[256];

   sprintf( buff, "%c,%d,%d,0,%f,%f",
               default_comet_magnitude_type,
               element_format, element_precision,
               max_residual_for_filtering, noise_in_arcseconds);
   set_environment_ptr( "SETTINGS", buff);
   set_environment_ptr( "EPHEM_OPTIONS",
                   write_bit_string( buff, ephemeris_output_options));
   sprintf( buff, "%.3f %f %u %d", probability_of_blunder * 100.,
              overobserving_time_span,
              overobserving_ceiling,
              use_blunder_method);
   set_environment_ptr( "FILTERING", buff);
   sprintf( buff, "%d %.2f %d %d %d", use_sigmas ? 1 : 0,
                                 ephemeris_mag_limit,
                                 sigmas_in_columns_57_to_65,
                                 forced_central_body,
                                 apply_debiasing);
   set_environment_ptr( "SETTINGS2", buff);

   buff[0] = findorb_language;
   buff[1] = '\0';
   set_environment_ptr( "LANGUAGE", buff);
   return( 0);
}

/* Modified from code at https://www.informit.com/articles/article.aspx?p=23618&seqNum=4 .
The first instance of Find_Orb,  fo,  or fo_serve will attempt to put a lock
on /tmp/fo_lock.  Subsequent instances will see that and fail to get a lock.
Basically,  an interprocess mutex.     */

#ifndef _WIN32
int check_for_other_processes( const int locking)
{
   const char *lock_filename = "/tmp/fo_lock";
   static int lock_fd;
   static struct flock lock;
   int rval;

   if( locking)
      {
      lock_fd = open( lock_filename, O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR);
      lock.l_type = F_WRLCK;
      rval = fcntl( lock_fd, F_SETLK, &lock);
      }
   else           /* unlocking */
      {
      lock.l_type = F_UNLCK;
      fcntl( lock_fd, F_SETLK, &lock);
      close( lock_fd);
      if( !findorb_already_running)
         unlink( lock_filename);
      rval = 0;
      }
   return( rval);
}
#endif

int set_language( const int language)
{
   int i;

   findorb_language = language;
   for( i = 0; i < 12; i++)
      {                          /* create abbreviated month names */
      static char month_names[12][17];
      char tbuff[100], *tptr;

      strcpy( tbuff, get_find_orb_text( i + 1000));
      tptr = (char *)find_nth_utf8_char( tbuff, 3);
      *tptr = '\0';        /* truncate month name at third char */
      assert( strlen( tbuff) < 16);
      strcpy( month_names[i], tbuff);
      set_month_name( i + 1, month_names[i]);
      }
   return( 0);
}

unsigned always_included_perturbers;

int get_defaults( ephem_option_t *ephemeris_output_options, int *element_format,
         int *element_precision, double *max_residual_for_filtering,
         double *noise_in_arcseconds)
{
   ephem_option_t unused_ephemeris_output_options;
   int unused_element_format;
   int unused_element_precision;
   double unused_max_residual_for_filtering;
   double unused_noise_in_arcseconds;
   const char *language = get_environment_ptr( "LANGUAGE");
   extern double minimum_jd, maximum_jd;
   extern double maximum_observation_span;
   extern int use_config_directory;
   extern double *sr_orbits;
   extern unsigned max_n_sr_orbits;
   int use_sigmas_int;
   unsigned long obsolete_ephem_output_options;
   const char *override_fcct14_filename = get_environment_ptr( "FCCT14_FILE");
   const char *ephem_bitstring = get_environment_ptr( "EPHEM_OPTIONS");
   const char *eop_filename = get_environment_ptr( "EOP_FILE");
   const char *albedo = get_environment_ptr( "OPTICAL_ALBEDO");

#ifndef _WIN32
   findorb_already_running = (check_for_other_processes( 1) != 0);
#endif
   if( *language)
      set_language( *language);
   if( *override_fcct14_filename)
      {
      extern const char *fcct14_bias_file_name;

      fcct14_bias_file_name = override_fcct14_filename;
      }
   if( eop_filename)
      load_earth_orientation_params( eop_filename, NULL);
   reset_td_minus_dt_string( get_environment_ptr( "DELTA_T"));
   sscanf( get_environment_ptr( "MAX_OBSERVATION_SPAN"), "%lf",
                                  &maximum_observation_span);
   if( !ephemeris_output_options)
      ephemeris_output_options = &unused_ephemeris_output_options;
   if( !element_format)
      element_format = &unused_element_format;
   if( !element_precision)
      element_precision = &unused_element_precision;
   if( !max_residual_for_filtering)
      max_residual_for_filtering = &unused_max_residual_for_filtering;
   if( !noise_in_arcseconds)
      noise_in_arcseconds = &unused_noise_in_arcseconds;
   if( sscanf( get_environment_ptr( "TIME_RANGE"), "%lf,%lf",
               &minimum_jd, &maximum_jd) == 2)
      {
      minimum_jd = YEAR_TO_JD( minimum_jd);
      maximum_jd = YEAR_TO_JD( maximum_jd);
      }
   *ephemeris_output_options = 0;
   sscanf( get_environment_ptr( "SETTINGS"), "%c,%d,%d,%lu,%lf,%lf",
               &default_comet_magnitude_type,
               element_format, element_precision,
               &obsolete_ephem_output_options,
               max_residual_for_filtering, noise_in_arcseconds);

   if( !*ephem_bitstring)     /* set defaults;  see 'environ.def' */
      ephem_bitstring = "10,16";
   *ephemeris_output_options = parse_bit_string( ephem_bitstring);

   sscanf( get_environment_ptr( "FILTERING"), "%lf %lf %u %d",
               &probability_of_blunder,
               &overobserving_time_span,
               &overobserving_ceiling,
               &use_blunder_method);
   probability_of_blunder /= 100.;        /* cvt from percent to raw prob */
   sscanf( get_environment_ptr( "SETTINGS2"), "%d %lf %d %d %d",
               &use_sigmas_int, &ephemeris_mag_limit,
               &sigmas_in_columns_57_to_65, &forced_central_body,
               &apply_debiasing);

   use_sigmas = (use_sigmas_int ? true : false);
   if( *get_environment_ptr( "COMBINE_ALL"))
      {
      extern const char *combine_all_observations;

      combine_all_observations = "";
      }
   if( use_config_directory)
      {
      char cospar_name[255];

      make_config_dir_name( cospar_name, "cospar.txt");
      load_cospar_file( cospar_name);
      }
   max_n_sr_orbits = atoi( get_environment_ptr( "MAX_SR_ORBITS"));
   if( !max_n_sr_orbits)
      max_n_sr_orbits = 500;
   sr_orbits = (double *)calloc( (size_t)max_n_sr_orbits, 7 * sizeof( double));
   assert( sr_orbits);
   sscanf( get_environment_ptr( "PERTURBERS"), "%x", &always_included_perturbers);
   if( *albedo)
      optical_albedo = atof( albedo);
   return( 0);
}

/* The following functions are used to "color" observations in the
console versions of Find_Orb.  The idea resembles that of the
four-color map problem,  except in this case, we'd like to show
observations in max_n_colors,  such that adjacent observations from
different MPC codes show up in different colors. You can't always do
this.  For example,  with observations from eight observatories,
mixed up so that each "pair" occurs,  you'd obviously need eight
different colors.  This code just tries a lot of possible colorings
and returns the one resulting in the fewest mismatches.

   To do this,  it uses an "annealing" sort of algorithm:  it first
sets the color for each MPC observatory at random (a value from zero
to max_n_colors - 1).  It then uses the improve_mpc_colors() to get
a better solution;  that function can make "obvious" improvements,
such as "if we change this MPC code to red,  there will be fewer
mismatches".  If the result has no mismatches,  we're home free and
stop looking for a "better" solution.  Otherwise,  we set a new set
of random colors and try to improve that... and repeat the procedure
for up to two seconds;  it's probably not worth spending much more time
on it than that.

   There are other ways to do this,  of course,  including a formal
pruned tree search among all possible color combinations.  But this
appears to work quite well,  and I thought it would result in simpler
code.  (I'm no longer so sure of that.  But I don't think I'll spend
the time to write a tree search version.)
*/

#define NO_MPC_COLOR_SET   -1

int find_mpc_color( const MPC_STATION *sdata, const char *mpc_code)
{
   int rval = NO_MPC_COLOR_SET;

   if( !mpc_code)       /* indicates 'just count colors */
      {
      rval = 0;
      while( sdata[rval].code[0])
         rval++;
      }
   else while( rval == NO_MPC_COLOR_SET && sdata->code[0])
      {
      if( mpc_code[0] == sdata->code[0] &&
             mpc_code[1] == sdata->code[1] &&
             mpc_code[2] == sdata->code[2])
         rval = sdata->color;
      sdata++;
      }
   return( rval);
}

static void set_mpc_colors_semirandomly( MPC_STATION *sdata,
               const int max_n_colors, unsigned long seed)
{
   int i;

   for( i = 0; sdata[i].code[0]; i++)
      sdata[i].color = (char)( i % max_n_colors);
   srand( seed);
   while( --i > 0)
      {
      const int j = (int)( rand( ) % i);
      const int color = sdata[j].color;

      sdata[j].color = sdata[i].color;
      sdata[i].color = color;
      }
}

/* After setting the colors at random,  we look for "simple" improvements:
for each MPC code,  we check the adjacent observations with different
MPC codes,  and see what colors they have.  We might see that (say)
there are four red neighbors,  three green,  and one blue.  In that case,
changing the color of the current MPC code to blue would result in only
one problem case,  instead of three or four.  (Ideally,  we'll find that
there are _no_ blue neighbors,  of course.)

   We keep trying this until no color changes are made.
*/

static void improve_mpc_colors( const int n_obs, const OBSERVE FAR *obs,
                   const int max_n_colors, MPC_STATION *sdata)
{
   int i, changes_made = 1, n_iterations = 0;

   while( changes_made)
      {
      changes_made = 0;
      for( i = 0; sdata[i].code[0]; i++)
         {
         int counts[20], j, color = sdata[i].color;

         assert( color >=0 && color < max_n_colors);
         for( j = 0; j < max_n_colors; j++)
            counts[j] = 0;
         for( j = 0; j < n_obs; j++)
            if( !strcmp( obs[j].mpc_code, sdata[i].code))
               {
               if( j > 0 && strcmp( obs[j - 1].mpc_code, sdata[i].code))
                  {
                  const int adjacent_color =
                          find_mpc_color( sdata, obs[j - 1].mpc_code);

                  assert( adjacent_color >=0 && adjacent_color < max_n_colors);
                  counts[adjacent_color]++;
                  }
               if( j < n_obs - 1 && strcmp( obs[j + 1].mpc_code, sdata[i].code))
                  {
                  const int adjacent_color =
                          find_mpc_color( sdata, obs[j + 1].mpc_code);

                  assert( adjacent_color >=0 && adjacent_color < max_n_colors);
                  counts[adjacent_color]++;
                  }
               }
         for( j = 0; j < max_n_colors; j++)
            if( counts[j] < counts[color])
               {
               color = j;
               sdata[i].color = (char)color;
               changes_made = 1;
               }
         assert( color >=0 && color < max_n_colors);
         }
      n_iterations++;
      assert( n_iterations < 10);
      }
}

extern int debug_level;

MPC_STATION *find_mpc_color_codes( const int n_obs, const OBSERVE FAR *obs,
                   const int max_n_colors)
{
   int n_codes = 0, i, j, n_alloced = 10;
   int best_score = 99999, best_seed = 0;
   clock_t t0;
   MPC_STATION *rval =
               (MPC_STATION *)calloc( n_alloced + 1, sizeof( MPC_STATION));

   for( i = 0; i < n_obs; i++)
      if( find_mpc_color( rval, obs[i].mpc_code) == NO_MPC_COLOR_SET)
         {
         int loc = 0;

         if( n_codes == n_alloced)
            {
            const int new_size = n_alloced * 2;
            MPC_STATION *new_array =
                   (MPC_STATION *)calloc( new_size + 1, sizeof( MPC_STATION));

            memcpy( new_array, rval, n_alloced * sizeof( MPC_STATION));
            n_alloced = new_size;
            free( rval);
            rval = new_array;
            }
         for( loc = 0; loc < n_codes
              && strcmp( rval[loc].code, obs[i].mpc_code) < 0; loc++)
            ;
                     /* move the rest of the array forward... */
         memmove( rval + loc + 1, rval + loc,
                          (n_codes - loc) * sizeof( MPC_STATION));
                     /* ...so we can copy in the new code: */
         strcpy( rval[loc].code, obs[i].mpc_code);
         n_codes++;
         }
   if( debug_level)
      {
      debug_printf( "%d obs, %d codes\n", n_obs, n_codes);
      for( i = 0; i < n_codes; i++)
         debug_printf( "%d: '%s'\n", i, rval[i].code);
      }
   t0 = clock( );
   for( i = 1; best_score && clock( ) < t0 + 2 * CLOCKS_PER_SEC; i++)
      {
      int score = 0;

      set_mpc_colors_semirandomly( rval, max_n_colors, (unsigned long)i);
      improve_mpc_colors( n_obs, obs, max_n_colors, rval);
      for( j = 0; j < n_obs - 1; j++)
         if( strcmp( obs[j].mpc_code, obs[j + 1].mpc_code))
            if( find_mpc_color( rval, obs[j].mpc_code) ==
               find_mpc_color( rval, obs[j + 1].mpc_code))
                  score++;
      if( score < best_score)   /* "lower" is "better" */
         {
         best_score = score;
         best_seed = i;
         }

      if( debug_level > 1 && (i < 10 || !(i % 100)))
          debug_printf( "Seed: %d, score %d, best %d\n", i, score, best_score);
      }
   if( debug_level)
      debug_printf( "Color setting: best score %d, i = %d\n", best_score, i);
   set_mpc_colors_semirandomly( rval, max_n_colors, (unsigned long)best_seed);
   improve_mpc_colors( n_obs, obs, max_n_colors, rval);
   return( rval);
}
