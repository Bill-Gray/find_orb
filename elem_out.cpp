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
#include <sys/stat.h>
#include "watdefs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "mpc_func.h"
#include "date.h"
#include "afuncs.h"
#include "lunar.h"
#include "monte0.h"     /* for put_double_in_buff() proto */
#include "showelem.h"
#include "stringex.h"
#include "constant.h"

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

static const char *_extras_filename = "hints.txt";
static const char *_default_extras_filename = "hints.def";
extern int available_sigmas;
extern double optical_albedo;
extern unsigned perturbers;

#ifdef NOT_CURRENTLY_IN_USE
#define ssnprintf_append( obuff, ...) snprintf_append( obuff, sizeof( obuff), __VA_ARGS__)
#define ssnprintf( obuff, ...) snprintf( obuff, sizeof( obuff), __VA_ARGS__)
#endif
int store_defaults( const ephem_option_t ephemeris_output_options,
         const int element_format, const int element_precision,
         const double max_residual_for_filtering,
         const double noise_in_sigmas);           /* elem_out.cpp */
int get_defaults( ephem_option_t *ephemeris_output_options, int *element_format,
         int *element_precision, double *max_residual_for_filtering,
         double *noise_in_sigmas);                /* elem_out.cpp */
int64_t nanoseconds_since_1970( void);                      /* mpc_obs.c */
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
void attempt_extensions( OBSERVE *obs, const int n_obs, double *orbit,
                  const double epoch);                  /* orb_func.cpp */
double calc_obs_magnitude( const double obj_sun,
          const double obj_earth, const double earth_sun, double *phase_ang);
int find_best_fit_planet( const double jd, const double *ivect,
                                 double *rel_vect);         /* runge.cpp */
void find_relative_state_vect( const double jd, const double *ivect,
               double *ovect, const int ref_planet);        /* runge.cpp */
void compute_effective_solar_multiplier( const char *constraints);   /* runge.c */
int get_orbit_from_mpcorb_sof( const char *object_name, double *orbit,
             ELEMENTS *elems, const double full_arc_len, double *max_resid);
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
double mag_band_shift( const char mag_band, int *err_code);   /* elem_out.c */
int get_jpl_ephemeris_info( int *de_version, double *jd_start, double *jd_end);
double *get_asteroid_mass( const int astnum);   /* bc405.cpp */
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
               const int n_obs, const double *orbit, const double epoch);
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
int compute_canned_object_state_vect( double *loc, const char *mpc_code,
                     const double jd);                /* elem_out.cpp */
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
char *make_config_dir_name( char *oname, const char *iname);  /* miscell.cpp */
int earth_lunar_posn( const double jd, double FAR *earth_loc, double FAR *lunar_loc);
double vect_diff2( const double *a, const double *b);
int get_residual_data( const OBSERVE *obs, double *xresid, double *yresid);
int put_elements_into_sof( char *obuff, const char *templat,
         const ELEMENTS *elem, const double *nongravs,
         const int n_obs, const OBSERVE *obs);                /* elem_ou2.cpp */
int qsort_strcmp( const void *a, const void *b, void *ignored_context);
uint64_t parse_bit_string( const char *istr);                /* miscell.cpp */
const char *write_bit_string( char *ibuff, const uint64_t bits,
                                          const size_t max_bitstring_len);
int save_ephemeris_settings( const ephem_option_t ephemeris_output_options,
      const int n_steps, const char *obscode, const char *step_size,
      const char *ephem_start, const char *config);      /* elem_out.cpp */
int load_ephemeris_settings( ephem_option_t *ephemeris_output_options,
      int *n_steps, char *obscode, char *step_size, char *ephem_start,
      const char *config);                               /* elem_out.cpp */
double generate_mc_variant_from_covariance( double *var_orbit,
                                    const double *orbit);    /* orb_func.c */
void rotate_state_vector_to_current_frame( double *state_vect,
                  const double epoch_shown, const int planet_orbiting,
                  char *body_frame_note);               /* elem_out.cpp */
void *bsearch_ext_r( const void *key, const void *base0, size_t nmemb,
      const size_t size, int (*compar)(const void *, const void *, void *),
      void *arg, bool *found);                           /* shellsor.cpp */

int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;

/* Old MSVCs and OpenWATCOM lack erf() and many other math functions: */

#if defined( _MSC_VER) && (_MSC_VER < 1800) || defined( __WATCOMC__)
double erf( double x);     /* orb_fun2.cpp */
#endif

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
   int i = 0;

   while( i < n_obs - 1 && !obs[i].is_included)
      i++;

   if( first)
      *first = i;
   if( last)
      for( *last = n_obs - 1; *last > i && !obs[*last].is_included; (*last)--)
         ;
}

void make_observatory_info_text( char *text, const size_t textlen,
             const OBSERVE *obs, int n_obs, const char *mpc_code)
{
   double jd_start = 0., jd_end = 0.;
   int n_found = 0, n_used = 0;
   char date_text[80], time_text[80];

   while( n_obs--)
      {
      if( !strcmp( obs->mpc_code, mpc_code))
         {
         if( !jd_start)
            jd_start = obs->jd;
         jd_end = obs->jd;
         n_found++;
         if( obs->is_included)
            n_used++;
         }
      obs++;
      }
   assert( n_found);
   snprintf_err( text, textlen, "%d observations from (%s)", n_found, mpc_code);
   if( n_used == n_found)
      strlcat_err( text, "\n", textlen);
   else
      snprintf_append( text, textlen, "; %d used\n", n_used);
   full_ctime( date_text, jd_start, CALENDAR_JULIAN_GREGORIAN
            | FULL_CTIME_YMD | FULL_CTIME_LEADING_ZEROES | FULL_CTIME_MICRODAYS);
   full_ctime( time_text, jd_start, FULL_CTIME_TIME_ONLY);
   snprintf_append( text, textlen, "First %s = %s\n", date_text, time_text);
   full_ctime( date_text, jd_end, CALENDAR_JULIAN_GREGORIAN
            | FULL_CTIME_YMD | FULL_CTIME_LEADING_ZEROES | FULL_CTIME_MICRODAYS);
   full_ctime( time_text, jd_end, FULL_CTIME_TIME_ONLY);
   snprintf_append( text, textlen, "Last  %s = %s", date_text, time_text);
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
   const size_t obuff_size = 37;
   const double dt = jd2 - jd1;

   if( year == year2)
      {
      snprintf_err( obuff, obuff_size, "%ld %s %d", year, month_names[month - 1], day1);
      obuff += strlen( obuff);
      if( month == month2 && day1 != day2)
         snprintf_err( obuff, obuff_size, "-%d", day2);
      else if( month != month2)
         snprintf_err( obuff, obuff_size, "-%s %d", month_names[month2 - 1], day2);
      }
   else              /* different years */
      snprintf_err( obuff, obuff_size, "%ld %s %d-%ld %s %d", year, month_names[month - 1],
                             day1, year2, month_names[month2 - 1], day2);

   if( dt < 10. / seconds_per_day)  /* less than 10 seconds: show to .01 sec */
      snprintf_append( obuff, obuff_size, " (%.2f sec)", dt * seconds_per_day);
   else if( dt < 100. / seconds_per_day) /* less than 100 seconds: show to .1 sec */
      snprintf_append( obuff, obuff_size, " (%.1f sec)", dt * seconds_per_day);
   else if( dt < 100. / minutes_per_day)     /* less than 100 minutes: show in min */
      snprintf_append( obuff, obuff_size, " (%.1f min)", dt * minutes_per_day);
   else if( dt < 2.)
      snprintf_append( obuff, obuff_size, " (%.1f hr)", dt * hours_per_day);
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
   debug_printf( "Requested index %d in language %c not found\n",
                     index, findorb_language);
   assert( 0);             /* i.e.,  should never get here */
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
   size_t obuff_size = 80;

   get_first_and_last_included_obs( obs, n_obs, &first_idx, &last_idx);
   for( i = n_included = 0; i < n_obs; i++)
      n_included += obs[i].is_included;
   if( options == -1)      /* 'guide.txt' bare-bones format */
      snprintf_err( obuff, obuff_size, "%d of %d", n_included, n_obs);
   else if( (options & ELEM_OUT_ALTERNATIVE_FORMAT) && n_included != n_obs)
      snprintf_err( obuff, obuff_size, get_find_orb_text( 15), n_included, n_obs);
   else
      snprintf_err( obuff, obuff_size, get_find_orb_text( 16), n_included);
   if( options != -1 && n_included)
      {
      const double rms = (options & ELEM_OUT_NORMALIZED_MEAN_RESID) ?
                  compute_weighted_rms( obs, n_obs, NULL) :
                  compute_rms( obs, n_obs);
      char rms_buff[24];
      const char *rms_format = "%.2f";

      strlcat_err( obuff, " ", obuff_size);
      obuff += strlen( obuff);
      make_date_range_text( obuff, obs[first_idx].jd, obs[last_idx].jd);
      obuff += strlen( obuff);
      if( options & ELEM_OUT_PRECISE_MEAN_RESIDS)
         rms_format = (rms > 0.003 ? "%.3f" : "%.1e");
      snprintf_err( rms_buff, sizeof( rms_buff), rms_format, rms);
      if( options & ELEM_OUT_NORMALIZED_MEAN_RESID)
         strlcat_err( rms_buff, " sigmas", sizeof( rms_buff));
      else
         text_search_and_replace( rms_buff, ".", "\".");
      snprintf_err( obuff, obuff_size, get_find_orb_text( 17), rms_buff);
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

/* For some time,  I only had second-hand info on MPC's definition of
an 'opposition'.  It sounded as if their definition was that 17 half-months
had elapsed between consecutive observations,  with the length varying
depending on the weirdnesses of the Gregorian calendar.

I ignored calendrical issues and just took it to be 237 days,  since
that's a prime number and is "close enough".  Margaret Pan at MPC tells me
the actual definition is :

-- If the arc is less than 238 days,  it's one opposition.
-- If it's greater than that,  the times of solar conjunction are computed.
If there are observations between consecutive conjunctions,  that's an
opposition at which the object was observed.  She notes that "for earth
co-orbital objects,  this can sometimes depend on the resolution of the
conjunction computation and the range of elongations allowed."

   I don't expect to implement this scheme soon.  Thus far,  all I've
done is to change my 237-day span to 238 days.  The result is usually
fairly close to the MPC value.         */

const double opposition_time = 238.;

bool opposition_break( const OBSERVE *obs)
{
   return( obs[1].jd - obs[0].jd > opposition_time);
}

static int _n_oppositions( const OBSERVE *obs, const int n)
{
   int i, j, rval = 1;

   for( i = j = 0; i < n; i++)
      if( obs[i].jd - obs[j].jd > opposition_time)
         {
         rval++;
         j = i;
         }
   return( rval);
}

int n_clones_accepted = 0;

static int elements_in_mpcorb_format( char *buff, const char *packed_desig,
                const char *full_desig, const ELEMENTS *elem,
                const OBSERVE FAR *obs, const int n_obs)   /* orb_func.c */
{
   int month, day, i, first_idx, last_idx, n_included_obs = 0;
   long year;
   const double rms_err = compute_rms( obs, n_obs);
   const unsigned hex_flags = 0;
            /* 'mpcorb' has four hexadecimal flags starting in column 162, */
            /* signifying if the object is in any of various classes such  */
            /* as Aten,  scattered-disk object,  PHA,  Jupiter Trojan,     */
            /*  etc.  None of those flags are set yet.                     */
   double arc_length;
   char packed_desig2[40];
   const size_t mpcorb_line_len = 203;
   int precision = 7;
   double tval;
   double abs_mag = elem->abs_mag;

   packed_desig_minus_spaces( packed_desig2, packed_desig);
   if( 12 == strlen( packed_desig2))  /* fix cases where number &  */
      packed_desig2[5] = '\0';        /* provisional ID are both set; */
   packed_desig2[8] = '\0';           /* prevent overrun otherwise */
   if( abs_mag > 99.9)
      abs_mag = 99.9;
   if( abs_mag < -9.9)
      abs_mag = -9.9;
   snprintf_err( buff, mpcorb_line_len, "%-8s%5.2f %5.2f ", packed_desig2, abs_mag,
                           asteroid_magnitude_slope_param);
   day = (int)( decimal_day_to_dmy( elem->epoch, &year,
                              &month, CALENDAR_JULIAN_GREGORIAN) + .0001);
   assert( 20 == strlen( buff));
   snprintf_append( buff, mpcorb_line_len, "%c%02ld%X%c",
                  (year >= 0 && year < 6200) ? int_to_mutant_hex_char( year / 100) : '~',
                  abs( year) % 100L, month,
                  int_to_mutant_hex_char( day));
   assert( 25 == strlen( buff));
   snprintf_append( buff, mpcorb_line_len, "%10.5f%11.5f%11.5f%11.5f%11.7f",
           centralize_ang( elem->mean_anomaly) * 180. / PI,
           centralize_ang( elem->arg_per) * 180. / PI,
           centralize_ang( elem->asc_node) * 180. / PI,
           centralize_ang( elem->incl) * 180. / PI,
           elem->ecc);
   assert( 79 == strlen( buff));
   tval = elem->major_axis;
   while( tval > 999.9999)
      {
      tval /= 10.;
      precision--;
      }
   snprintf_append( buff, mpcorb_line_len, "%12.8f%12.*f",
            (180 / PI) / elem->t0,        /* n */
            precision, elem->major_axis);
   if( 103 != strlen( buff))
      printf( "Weirdness '%s'\n", buff);
   assert( 103 == strlen( buff));
   for( i = 0; i < n_obs; i++)
      if( obs[i].is_included)
         n_included_obs++;
   day = (int)( decimal_day_to_dmy( current_jd( ),
                         &year, &month, CALENDAR_JULIAN_GREGORIAN) + .0001);
   snprintf_append( buff, mpcorb_line_len,
      "    FO %02d%02d%02d %5d %3d ****-**** ****         Find_Orb   %04x",
                  (int)( year % 100), month, (int)day,
                  n_included_obs, _n_oppositions( obs, n_obs), hex_flags);
   get_first_and_last_included_obs( obs, n_obs, &first_idx, &last_idx);
   arc_length = obs[last_idx].jd - obs[first_idx].jd;
   assert( strlen( buff) == 165);
   if( arc_length < 99. / seconds_per_day)
      snprintf_err( buff + 127, 10, "%4.1f sec ", arc_length * seconds_per_day);
   else if( arc_length < 99. / minutes_per_day)
      snprintf_err( buff + 127, 10, "%4.1f min ", arc_length * minutes_per_day);
   else if( arc_length < 2.)
      snprintf_err( buff + 127, 10, "%4.1f hrs ", arc_length * hours_per_day);
   else if( arc_length < 600.)
      snprintf_err( buff + 127, 10, "%4d days", (int)arc_length + 1);
   else
      snprintf_err( buff + 127, mpcorb_line_len, "%4d-%4d",
                (int)JD_TO_YEAR( obs[first_idx].jd),
                (int)JD_TO_YEAR( obs[last_idx].jd));
   buff[136] = ' ';
   assert( 165 == strlen( buff));
   if( strlen( full_desig) > 28)
      strlcat( buff, full_desig, 195);
   else
      snprintf_append( buff, mpcorb_line_len, " %-28s", full_desig);
   day = (int)( decimal_day_to_dmy( obs[last_idx].jd, &year,
                       &month, CALENDAR_JULIAN_GREGORIAN) + .0001);
   assert( 194 == strlen( buff));
   snprintf_append( buff, mpcorb_line_len, "%04ld%02d%02d", year, month, day);
   if( rms_err < 9.9)
      snprintf_err( buff + 137, 5, "%4.2f", rms_err);
   else if( rms_err < 99.9)
      snprintf_err( buff + 137, 5, "%4.1f", rms_err);
   else if( rms_err < 9999.)
      snprintf_err( buff + 137, 5, "%4.0f", rms_err);
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

static int _unpack_desig_for_linkage( const char *packed_id, char *reduced)
{
   size_t i = 0;

   strlcpy_err( reduced, packed_id, 13);
   while( i < 12 && reduced[i] != ' ')
      i++;
   if( i == 12)                  /* both permanent and provisional desigs; */
      memset( reduced + 5, ' ', 7);          /* just use the permanent one */
   return( unpack_mpc_desig( NULL, reduced));
}

double utc_from_td( const double jdt, double *delta_t);     /* ephem0.cpp */

/* This is ludicrously high,  but let's say that we'll never
report a linkage between 50 or more tracklets at once. */

#define MAX_LINKAGE_IDS 50

static int make_linkage_json( const int n_obs, const OBSERVE *obs, const ELEMENTS *elem)
{
   int i, j, n_ids = 0, idx[MAX_LINKAGE_IDS], n_designated = 0;
   FILE *ofile, *ifile;
   char buff[200], packed_id2[13];

   for( i = 0; i < n_obs && n_ids < MAX_LINKAGE_IDS; i++)
      if( !i || strcmp( obs[i].packed_id, obs[i - 1].packed_id))
         {
         const int desig_type = _unpack_desig_for_linkage( obs[i].packed_id, packed_id2);

         for( j = 0; j < n_ids; j++)
            {
            char packed_id[13];

            _unpack_desig_for_linkage( obs[idx[j]].packed_id, packed_id);
            if( !strcmp( packed_id, packed_id2))
               break;
            }
         if( j == n_ids)
            {
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
   ofile = fopen_ext( (elem->central_obj == 3) ? "artsat.json" : "linkage.json", "tfcw");
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
         tptr = strstr( buff, "%v");
         if( tptr)
            text_search_and_replace( buff, "%v", find_orb_version_jd( NULL));
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
         const int desig_type = _unpack_desig_for_linkage( obs[idx[i]].packed_id, packed_id2);

         if( desig_type != OBJ_DESIG_OTHER && desig_type != OBJ_DESIG_ARTSAT)
            {
            strlcpy_error( buff, packed_id2);
            text_search_and_replace( buff, " ", "");
            j++;
            fprintf( ofile, "        \"%s\"%s\n", buff, (j == n_designated ? "" : ","));
            }
         }
      if( n_ids == n_designated)
         fprintf( ofile, "      ]\n");
      else
         fprintf( ofile, "      ],\n");
      }
   if( n_ids != n_designated)
      {
      fprintf( ofile, "      \"trksubs\": [\n");
      for( i = j = 0; i < n_ids; i++)
         {
         const int desig_type = _unpack_desig_for_linkage( obs[idx[i]].packed_id, packed_id2);

         if( desig_type == OBJ_DESIG_OTHER || desig_type == OBJ_DESIG_ARTSAT)
            {
            const double obs_utc_jd = utc_from_td( obs[idx[i]].jd, NULL);

            strlcpy_error( buff, packed_id2);
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
      fprintf( ofile, "      ]");
      }
   i = 0;
   while( i < n_obs && strcmp( obs[i].reference, "NEOCP"))
      i++;
   if( i < n_obs)         /* no NEOCP observations */
      {
      fprintf( ofile, ",\n");
      fprintf( ofile, "      \"identification_type\": \"neocp\"");
      }
   if( elem->central_obj != 3)
      {
      fprintf( ofile, ",\n");
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
      fprintf( ofile, "      }");
      }
   fprintf( ofile, "\n    }\n  }\n}\n");
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
      snprintf_err( buff, 12, "Object_%d", obj_index);
   return( buff);
}

double vect3_squared( const double *v)
{
   return( v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

int geo_score( const double *orbit, const double epoch)
{
   extern unsigned n_sr_orbits;
   unsigned i, n_geo, n_orbits;
   const double gm = get_planet_mass( 3);
   extern int n_orbit_params;

   if( n_orbit_params > 6) /* don't bother with a geo score for objects */
      return( -1);         /* so well-determined that we have non-gravs */
   if( available_sigmas == COVARIANCE_AVAILABLE)
      n_orbits = 100;
   else if( available_sigmas == SR_SIGMAS_AVAILABLE)
      n_orbits = n_sr_orbits;
   else
      return( -1);
   for( i = n_geo = 0; i < n_orbits; i++)
      {
      double variant_orbit[MAX_N_PARAMS], earth_centric_orbit[MAX_N_PARAMS];
      double r, vel_squared;
      extern double *sr_orbits;

      if( available_sigmas == COVARIANCE_AVAILABLE)
         generate_mc_variant_from_covariance( variant_orbit, orbit);

      else if( available_sigmas == SR_SIGMAS_AVAILABLE)
         memcpy( variant_orbit, sr_orbits + 6 * i, 6 * sizeof( double));
      find_relative_state_vect( epoch, variant_orbit, earth_centric_orbit, 3);
      r = vector3_length( earth_centric_orbit);
      vel_squared = vect3_squared( earth_centric_orbit + 3);
      if( vel_squared < gm / r)
         n_geo++;
      }
   return( n_geo * 100 / n_orbits);
}

static int elements_in_json_format( FILE *ofile, const ELEMENTS *elem,
                     const double *orbit,
                     const char *obj_name, const OBSERVE *obs,
                     const unsigned n_obs, const double *moids,
                     const char *body_frame_note, const bool show_obs,
                     const int geocentric_score)
{
   extern int n_orbit_params;
   double jd = current_jd( );
   double jd_first, jd_last;
   double q_sigma = 0., weighted_rms;
   char buff[180];
   int first, last, i, n_used;
   extern const char *combine_all_observations;
   const char *packed_id;
   extern double uncertainty_parameter;
   const char *reference = get_environment_ptr( "REFERENCE");

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
         for( i = 5; i >= 0; i--)
            {
            const char *fields[] = { "name", "prov_desig", "disc_date", "disc_site",
                           "disc_ref", "discover" };
            const int offsets[] = { 9, 29, 41, 53, 72, 78 };
            char *tptr = (char *)find_nth_utf8_char( text, offsets[i]);

            remove_trailing_cr_lf( tptr);
            if( *tptr && *tptr != ' ')
               fprintf( ofile, "      \"%s\": \"%s\",\n", fields[i], tptr);
            *tptr = '\0';
            }
         }
      }
   fprintf( ofile, "      \"packed\": \"%s\",\n",
                  real_packed_desig( buff, packed_id));
   fprintf( ofile, "      \"created\": %.5f,\n", jd);
   fprintf( ofile, "      \"created iso\": \"%s\",\n", iso_time( buff, jd, 0));
   find_orb_version_jd( &jd);
   fprintf( ofile, "      \"Find_Orb_version\": %.5f,\n", jd);
   fprintf( ofile, "      \"Find_Orb_version_iso\": \"%s\",\n", iso_time( buff, jd, 0));
   fprintf( ofile, "      \"elements\":\n      {\n");
   fprintf( ofile, "        \"central body\": \"%s\",\n", object_name( buff, elem->central_obj));
   strlcpy_error( buff, body_frame_note + 1);
   i = (int)strlen( buff);
   buff[i - 1] = '\0';    /* strip trailing paren */
   fprintf( ofile, "        \"frame\": \"%s\",\n", buff);
   if( *reference)
      fprintf( ofile, "        \"reference\": \"%s\",\n", reference);
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
   for( i = 6; i < n_orbit_params; i++)
      {
      char tbuff[20];

      fprintf( ofile, "\n        \"A%d\": %.9g,", i - 5, orbit[i]);
      snprintf( tbuff, sizeof( tbuff), "Sigma_A%d:", i - 5);
      if( !get_uncertainty( tbuff, buff, 0))
         fprintf( ofile, " \"sigma_A%d\": %s,", i - 5, buff);
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
   if( geocentric_score > 0)
      fprintf( ofile, "\n        \"geo_score\": %d,", geocentric_score);

   fprintf( ofile, "\n        \"MOIDs\":");
   fprintf( ofile, "\n        {");
   for( i = 1; i <= 8; i++)
      fprintf( ofile, "\n          \"%s\" : %.6f%c",
                  object_name( buff, i), moids[i],
                  (i == 8) ? ' ' : ',');
   fprintf( ofile, "\n        }");

   fprintf( ofile, "\n      },\n      \"observations\":\n      {");
   get_first_and_last_included_obs( obs, n_obs, &first, &last);
   for( i = n_used = 0; i < (int)n_obs; i++)
      n_used += (obs[i].is_included & 1);
   fprintf( ofile, "\n        \"count\": %u,", n_obs);
   fprintf( ofile, "\n        \"used\": %u,", n_used);

   jd_first = utc_from_td( obs[0].jd, NULL);
   jd_last  = utc_from_td( obs[n_obs - 1].jd, NULL);
   fprintf( ofile, "\n        \"earliest\": %16.8f,", jd_first);
   fprintf( ofile, "\n        \"latest\": %16.8f,", jd_last);
   fprintf( ofile, "\n        \"earliest iso\": \"%s\",", iso_time( buff, jd_first, 3));
   fprintf( ofile, "\n        \"latest iso\": \"%s\",", iso_time( buff, jd_last, 3));

   jd_first = utc_from_td( obs[first].jd, NULL);
   jd_last  = utc_from_td( obs[last].jd, NULL);
   fprintf( ofile, "\n        \"earliest_used\": %16.8f,", jd_first);
   fprintf( ofile, "\n        \"latest_used\": %16.8f,", jd_last);
   fprintf( ofile, "\n        \"earliest_used iso\": \"%s\",", iso_time( buff, jd_first, 3));
   fprintf( ofile, "\n        \"latest_used iso\": \"%s\",", iso_time( buff, jd_last, 3));

   first = 0;
   while( first < (int)n_obs - 1 && (obs[first].flags & OBS_DONT_USE))
      first++;
   last = n_obs - 1;
   while( last > first && (obs[last].flags & OBS_DONT_USE))
      last--;
   jd_first = utc_from_td( obs[first].jd, NULL);
   jd_last  = utc_from_td( obs[last].jd, NULL);
   fprintf( ofile, "\n        \"earliest_unbanned\": %16.8f,", jd_first);
   fprintf( ofile, "\n        \"latest_unbanned\": %16.8f,", jd_last);
   fprintf( ofile, "\n        \"earliest_unbanned iso\": \"%s\",", iso_time( buff, jd_first, 3));
   fprintf( ofile, "\n        \"latest_unbanned iso\": \"%s\",", iso_time( buff, jd_last, 3));

   fprintf( ofile, "\n        \"residuals\":\n        [");
   for( i = 0; show_obs && i < (int)n_obs; i++)
      {
      MOTION_DETAILS m;
      int n_digits = 3;
      double total_resid, normalized_xresid, normalized_yresid;
      double ecliptic_lon, ecliptic_lat;

      jd = utc_from_td( obs[i].jd, NULL);
      compute_observation_motion_details( obs + i, &m);
      fprintf( ofile, "\n          {\"JD\": %.6f, \"iso date\": \"%s\", \"obscode\": \"%s\",",
                  jd, iso_time( buff, jd, 3), obs[i].mpc_code);
      fprintf( ofile, "\n                 \"RA\" : %.6f, \"Dec\": %.6f,",
                  obs[i].ra * 180. / PI, obs[i].dec * 180. / PI);
      total_resid = hypot( m.xresid, m.yresid);
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
      get_residual_data( obs + i, &normalized_xresid, &normalized_yresid);
      fprintf( ofile, "\n                 \"normalized_resid_1\" : %f,",
                                                normalized_xresid);
      fprintf( ofile, "\n                 \"normalized_resid_2\" : %f,",
                                                normalized_yresid);
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
      fprintf( ofile, "\n                 \"incl\" : %d,", obs[i].is_included);
      ecliptic_lon = atan2( obs[i].vect[1], obs[i].vect[0]);
      ecliptic_lat = asine( obs[i].vect[2]);
      fprintf( ofile, "\n                 \"ecliptic_lon\" : %f,",
                                    centralize_ang( ecliptic_lon) * 180. / PI);
      fprintf( ofile, "\n                 \"ecliptic_lat\" : %f,",
                                                   ecliptic_lat * 180. / PI);
      ecliptic_lon = 100.46435 + (obs[i].jd - J2000) * 0.9856091005;
      ecliptic_lon *= PI / 180.;
      fprintf( ofile, "\n                 \"earth_ecliptic_lon\" : %f",
                                    centralize_ang( ecliptic_lon) * 180. / PI);
      fprintf( ofile, " }%c", (i == (int)n_obs - 1 ? ' ' : ','));
      }
   fprintf( ofile, "\n        ]\n      }");
   fprintf( ofile, "\n    }\n  }\n}\n");
   return( 0);
}

/* First stab at writing out orbital elements in a form that can be
readily imported to JPL Horizon's "plain text" format.  That format
is noted as preliminary/in progress/unreliable,  so this may change. */

static int write_horizons_elems( const char *filename, const ELEMENTS *elem, const double *orbit)
{
   FILE *ofile = fopen( filename, "wb");
   extern int force_model;
   int i;

   assert( ofile);
   fprintf( ofile, "Object: From_Find_Orb\n");
   fprintf( ofile, "H= %.3f    G= %.2f\n", elem->abs_mag, elem->slope_param);
   fprintf( ofile, "EPOCH= %.3f\n", elem->epoch);
   fprintf( ofile, "EC= %.11f     QR= %.11f   TP= %.8f\n",
                        elem->ecc, elem->q, elem->perih_time);
   fprintf( ofile, "OM= %.7f     W= %.7f    IN= %.7f\n",
                     elem->asc_node * 180. / PI,
                     elem->arg_per * 180. / PI,
                     elem->incl * 180. / PI);
   for( i = 0; i < (force_model & 0xf); i++)
      fprintf( ofile, "A%d= %.6e    ", i + 1, orbit[i + 6]);
   fclose( ofile);
   return 0;
}

static int elements_in_guide_format( char *buff, const ELEMENTS *elem,
                     const char *obj_name, const OBSERVE *obs,
                     const unsigned n_obs)
{
   int month, q_prec = 10, e_prec = 8;
   double day, tval;
   long year;
   const size_t guide_line_len = 170;

   if( elem->q > 1e+10 || elem->ecc > 1e+8)
      return( -1);            /* can't format it;  bogus elems */
   day = decimal_day_to_dmy( elem->perih_time, &year, &month,
                                              CALENDAR_JULIAN_GREGORIAN);
            /*      name day  mon yr MA      q      e */
   tval = elem->q;
   while( tval > 9.999)
      {
      tval /= 10;
      q_prec--;
      }
   tval = elem->ecc;
   while( tval > 9.999)
      {
      tval /= 10;
      e_prec--;
      assert( e_prec);
      }
   strlcpy( buff, obj_name, 43);
   memset( buff + strlen( buff), ' ', 43);
   snprintf_err( buff + 43, guide_line_len,
            "%8.5f%3d%5ld Find_Orb    %12.*f %10.*f%11.6f %11.6f %11.6f",
            day, month, year,
            q_prec, elem->q, e_prec, elem->ecc,
            centralize_ang( elem->incl) * 180. / PI,
            centralize_ang( elem->arg_per) * 180. / PI,
            centralize_ang( elem->asc_node) * 180. / PI);
   assert( strlen( buff) == 130);
   snprintf_append( buff, guide_line_len,  " %9.1f%5.1f%5.1f %c",
            elem->epoch, elem->abs_mag,
            elem->slope_param * (elem->is_asteroid ? 1. : 0.4),
            (elem->is_asteroid ? 'A' : ' '));
   if( elem->central_obj)
      snprintf_append( buff, guide_line_len, "  Center: %d", elem->central_obj);
   strlcat_err( buff, "  ", guide_line_len);
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
      strcat( buff, "   ");
      strcat( buff, reference);
      }
   return( rval);
}

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
   tval = hypot( xform[2][0], xform[2][1]);
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
existing file,  one in which the first line is the header.  If it turns
out the file doesn't actually exist yet,  we use the 'fallback_filename'
to get the template line. */

static int add_sof_to_file( const char *filename,
             const ELEMENTS *elem, const double *nongravs,
             const int n_obs, const OBSERVE *obs, const char *fallback_filename)
{
   char templat[MAX_SOF_LEN], obuff[MAX_SOF_LEN];
   char output_filename[100];
   const char *obj_name = obs->packed_id;
   FILE *fp;
   int rval = -1, forking;

   while( *obj_name == ' ')
      obj_name++;
   get_file_name( output_filename, filename);
   forking = strcmp( output_filename, filename);
   fp = fopen_ext( output_filename, (forking ? "tr+b" : "cr+b"));
   if( !fp && fallback_filename)
      {
      fp = fopen_ext( fallback_filename, "fcr+b");
      if( !fgets( templat, sizeof( templat), fp))
         assert( 1);          /* should never happen */
      fclose( fp);
      fp = fopen_ext( output_filename, (forking ? "tfw+b" : "fcw+b"));
      fputs( templat, fp);
      }
   assert( fp);
   if( fp)
      {
      fseek( fp, 0L, SEEK_SET);
      if( !fgets( templat, sizeof( templat), fp))
         {
         fclose( fp);
         return( -1);
         }
      if( forking)                   /* simply append */
         fseek( fp, 0L, SEEK_END);
      else           /* we may be replacing an existing entry */
         {
         bool got_it = false;

         fseek( fp, 0L, SEEK_SET);
         while( !got_it && fgets( obuff, sizeof( obuff), fp))
            if( !memcmp( obuff, obj_name, strlen( obj_name)) &&
                        obuff[strlen( obj_name)] == ' ')
               {
               got_it = true;
               fseek( fp, -(long)strlen( obuff), SEEK_CUR);
               }
         }
      rval = put_elements_into_sof( obuff, templat, elem, nongravs, n_obs, obs);
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

static void _store_extra_orbit_info( const char *packed_id,
               const unsigned perturbers, const int n_extra_params,
               const double *solar_pressure, const char *constraints)
{
   extern int force_model;
   size_t n_lines;
   char **lines = load_file_into_memory( _extras_filename, &n_lines, false);
   FILE *ofile;
   int i, line_no = -1;
   char buff[200];

   if( !lines)
      lines = load_file_into_memory( _default_extras_filename, &n_lines, true);
   assert( strlen( packed_id) == 12);
   assert( (force_model && n_extra_params) || (!force_model && !n_extra_params));
   for( i = 0; i < (int)n_lines; i++)
      if( !memcmp( lines[i], packed_id, 12))
         line_no = i;
   if( (perturbers & ~0x7ff) || n_extra_params || *constraints)
      {
      strlcpy_error( buff, packed_id);
      if( perturbers & ~0x7ff)
         snprintf_append( buff, sizeof( buff), " p=%x", perturbers);
      if( force_model)
         snprintf_append( buff, sizeof( buff), " model=%x", force_model);
      for( i = 0; i < n_extra_params; i++)
         snprintf_append( buff, sizeof( buff), " A%d=%e", i + 1, solar_pressure[i]);
      if( *constraints)
         snprintf_append( buff, sizeof( buff), " Constraint=%s", constraints);
      if( line_no >= 0 && !strcmp( buff, lines[line_no]))
         {              /* already had this hint;  it didn't change; */
         free( lines);           /* we can bug out early */
         return;
         }
      }
   else
      {
      if( line_no == -1)      /* didn't have a hint before,  doesn't have */
         {                    /* have one now;  we can bug out early */
         free( lines);
         return;
         }
      *buff = '\0';
      }
   ofile = fopen_ext( _extras_filename, "fcw");
   for( i = 0; i < (int)n_lines; i++)
      if( i != line_no)
         fprintf( ofile, "%s\n", lines[i]);
   if( *buff)
      fprintf( ofile, "%s\n", buff);
   free( lines);
   fclose( ofile);
}

/* Packed designations may be (frequently are) misaligned.  This will
compare them 'accurately' even if the alignments are different,  ignoring
the different numbers of leading spaces. */

static int _compare_misaligned_packed_desigs( const char *desig1,
                  const char *desig2)
{
   while( *desig1 == ' ')
      desig1++;
   while( *desig2 == ' ')
      desig2++;
   while( *desig1 && *desig1 != ' ' && *desig1 == *desig2)
      {
      desig1++;
      desig2++;
      }
   if( !*desig1 || *desig1 == ' ')
      return( 0);
   return( (int)*desig1 - (int)*desig2);
}

static int _get_extra_orbit_info( const char *packed_id,
                unsigned *perturbers, int *n_extra_params,
                double *solar_pressure, char *constraints)
{
   extern int force_model;
   size_t n_lines;
   char **lines = load_file_into_memory( _extras_filename, &n_lines, false);
   int i, j, rval = 0;

   if( !lines)
      lines = load_file_into_memory( _default_extras_filename, &n_lines, true);
   assert( strlen( packed_id) == 12);
   for( i = 0; i < (int)n_lines; i++)
      if( !_compare_misaligned_packed_desigs( lines[i], packed_id))
         {
         const char *tptr;

         if( (tptr = strstr( lines[i], " p=")) != NULL)
            sscanf( tptr + 3, "%x", perturbers);
         if( (tptr = strstr( lines[i], " model=")) != NULL)
            sscanf( tptr + 7, "%x", (unsigned *)&force_model);
         for( j = 0; j < 5; j++)
            {
            char tbuff[6];

            snprintf( tbuff, sizeof( tbuff), " A%d=", j + 1);
            if( (tptr = strstr( lines[i], tbuff)) != NULL)
               {
               sscanf( tptr + 4, "%lf", solar_pressure + j);
               *n_extra_params = j + 1;
               }
            }
         if( (tptr = strstr( lines[i], " Constraint=")) != NULL)
            if( constraints)
               sscanf( tptr + 12, "%s", constraints);
         if( strstr( lines[i], "ignore"))
            rval = -1;
         break;
         }
   free( lines);
   return( rval);
}

/* see 'environ.def' for an explanation of how/why this works. */

static bool _include_comment( const char *keyval)
{
   return( !strstr( get_environment_ptr( "COMMENT_SHUTOFF"), keyval));
}


static void format_value_with_sigma( char *obuff, const double value, const char *sigma_text)
{
   const char *eptr = strchr( sigma_text, 'e');

   if( !eptr)
      {
      put_double_in_buff( obuff, value);
      strlcat_err( obuff, " +/- ", 80);
      strlcat_err( obuff, sigma_text, 80);
      }
   else
      {
      const char *decimal_ptr = strchr( sigma_text, '.');
      const int exponent = atoi( eptr + 1);
      int places = 0;
      char format[20], output_value[30];

      if( decimal_ptr)
         places = (int)( eptr - decimal_ptr) - 1;
      assert( places >= 0 && places < 10);
      strlcpy_error( format, "(%.?f +/- %s)e%s");
      format[3] = '0' + places;
      memcpy( output_value, sigma_text, eptr - sigma_text);
      output_value[eptr - sigma_text] = '\0';
      snprintf_err( obuff, 80, format, value * pow( 0.1, (double)exponent),
                                 output_value, eptr + 1);
      }
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

#define ELEMENT_FRAME_DEFAULT                   0
#define ELEMENT_FRAME_J2000_ECLIPTIC            1
#define ELEMENT_FRAME_J2000_EQUATORIAL          2
#define ELEMENT_FRAME_BODY_FRAME                3

void rotate_state_vector_to_current_frame( double *state_vect,
                  const double epoch_shown, const int planet_orbiting,
                  char *body_frame_note)
{
   int elements_frame = atoi( get_environment_ptr( "ELEMENTS_FRAME"));
   const char *frame_str = NULL;
   const size_t body_frame_note_len = 30;

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
         frame_str = "(J2000 ecliptic)";
      else
         {
         double year = atof( elem_epoch + 1);
         double precession_matrix[9];

         if( !year)
             year = JD_TO_YEAR( epoch_shown);
         if( body_frame_note)
            snprintf_err( body_frame_note, body_frame_note_len, "(%s)", elem_epoch);
         setup_ecliptic_precession( precession_matrix, 2000., year);
         precess_vector( precession_matrix, state_vect, state_vect);
         precess_vector( precession_matrix, state_vect + 3, state_vect + 3);
         }
      }
   if( elements_frame == ELEMENT_FRAME_J2000_EQUATORIAL)
      {
      ecliptic_to_equatorial( state_vect);
      ecliptic_to_equatorial( state_vect + 3);
      frame_str = "(J2000 equator)";
      }
   if( elements_frame == ELEMENT_FRAME_BODY_FRAME)
      {
      ecliptic_to_planetary_plane(
                  (planet_orbiting == 100 ? -1 : planet_orbiting),
                  epoch_shown, state_vect);
      frame_str = "(body frame)";
      }
   if( body_frame_note && frame_str)
      strlcpy_err( body_frame_note, frame_str, body_frame_note_len);
}

/* Used to provide a link to Tony Dunn's Orbit Simulator when making
pseudo-MPECs.  First six array elements are the actual state vector,
in the heliocentric J2000 ecliptic frame.  H and G are the sixth
and seventh array elements.  Epoch is stored in the eighth.   */

double helio_ecliptic_j2000_vect[9];

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
N_MOIDS objects (currently 14:  eight planets,  six asteroids).  We
default to 8 (planetary MOIDs only),  but can add six asteroid MOIDS
by setting N_MOIDS=14 in 'environ.dat'.  */

#define N_MOIDS           14

double comet_total_magnitude = 0.;          /* a.k.a. "M1" */
double comet_nuclear_magnitude = 0.;        /* a.k.a. "M2" */
bool saving_elements_for_reuse = false;

#define ORBIT_SUMMARY_Q                1
#define ORBIT_SUMMARY_H                2

int write_out_elements_to_file( const double *orbit,
            const double curr_epoch,
            const double epoch_shown,
            OBSERVE FAR *obs, const int n_obs, const char *constraints,
            const int precision, const int monte_carlo,
            const int options)
{
   char object_name[80], buff[320], more_moids[80];
   const char *file_permits = (append_elements_to_element_file ? "tfca" : "tfcw+");
   extern const char *elements_filename;
   FILE *ofile = fopen_ext( get_file_name( buff, elements_filename), file_permits);
   double rel_orbit[MAX_N_PARAMS], orbit2[MAX_N_PARAMS];
   int planet_orbiting, n_lines, i, bad_elements = 0;
   ELEMENTS elem, helio_elem;
   char *tptr, *tbuff;
   char impact_buff[80];
   int n_more_moids = 0;
   int output_format = (precision | SHOWELEM_PERIH_TIME_MASK);
   extern int n_orbit_params;
   int reference_shown = 0;
   double moids[N_MOIDS + 3];
   double j2000_ecliptic_rel_orbit[MAX_N_PARAMS];
   double barbee_style_delta_v = 0.;   /* see 'moid4.cpp' */
   const char *monte_carlo_permits;
   const bool rms_ok = (compute_rms( obs, n_obs) < max_monte_rms);
   extern int available_sigmas;
   int geocentric_score = -1;
   char body_frame_note[30];
   bool body_frame_note_shown = false;
   int showing_sigmas, saved_available_sigmas = available_sigmas;
   const unsigned orbit_summary_options = atoi( get_environment_ptr( "ORBIT_SUMMARY_OPTIONS"));
   extern int is_interstellar;         /* orb_func.cpp */
   const size_t tbuff_size = 80 * 9;
   const char *horizons_elems = get_environment_ptr( "HORIZONS_ELEMS");
   static const char *moid_text[N_MOIDS] = { "Earth MOID", "Ju",
                           "Ve", "Me", "Ma", "Sa", "Ur", "Ne",
                           "Ce", "Pa", "Vt", "(29)", "(16)", "(15)" };


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
   memcpy( orbit2, orbit, n_orbit_params * sizeof( double));
   integrate_orbit( orbit2, curr_epoch, epoch_shown);
   if( forced_central_body != ORBIT_CENTER_AUTO)
      {
      planet_orbiting = forced_central_body;
      get_relative_vector( epoch_shown, orbit2, rel_orbit, planet_orbiting);
      }
   else
      planet_orbiting = find_best_fit_planet( epoch_shown, orbit2, rel_orbit);

   memcpy( j2000_ecliptic_rel_orbit, rel_orbit, n_orbit_params * sizeof( double));
   rotate_state_vector_to_current_frame( rel_orbit, epoch_shown, planet_orbiting,
                       body_frame_note);

   if( available_sigmas == COVARIANCE_AVAILABLE)
      {
      extern int available_sigmas_hash;

      if( available_sigmas_hash != compute_available_sigmas_hash( obs, n_obs,
                  epoch_shown, perturbers, planet_orbiting))
         available_sigmas = 0;
      }
   if( !(options & ELEM_OUT_ALTERNATIVE_FORMAT))
      showing_sigmas = 0;
   else
      showing_sigmas = available_sigmas;

   elem.central_obj = planet_orbiting;
   elem.gm = get_planet_mass( planet_orbiting);
   elem.epoch = epoch_shown;
   calc_classical_elements( &elem, rel_orbit, epoch_shown, 1);
   if( elem.ecc < .9)
      snprintf_err( orbit_summary_text, sizeof( orbit_summary_text),
                                   "a=%.3f, ", elem.major_axis);
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
   memcpy( helio_ecliptic_j2000_vect, orbit2, 6 * sizeof( double));
   helio_ecliptic_j2000_vect[6] = elem.abs_mag;
   helio_ecliptic_j2000_vect[7] = elem.slope_param;
   helio_ecliptic_j2000_vect[8] = epoch_shown;

   add_sof_to_file( (n_orbit_params >= 8 ? "cmt_sof.txt" : sof_filename),
                    &elem, orbit + 6, n_obs, obs,
                    (n_orbit_params >= 8 ? "cometdef.sof" : "orbitdef.sof"));
   if( saving_elements_for_reuse && available_sigmas == COVARIANCE_AVAILABLE)
      {
      add_sof_to_file( "orbits.sof", &elem, orbit + 6, n_obs, obs, "orbitdef.sof");
      _store_extra_orbit_info( obs->packed_id, perturbers, n_orbit_params - 6,
                     orbit + 6, constraints);
      }
/* if( showing_sigmas == COVARIANCE_AVAILABLE)
*/    {
      ELEMENTS elem2 = elem;
      double rel_orbit2[MAX_N_PARAMS];

      if( showing_sigmas == COVARIANCE_AVAILABLE)
         {
         compute_variant_orbit( rel_orbit2, rel_orbit, 1.);    /* orb_func.cpp */
         calc_classical_elements( &elem2, rel_orbit2, epoch_shown, 1);
         }
      else        /* insert dummy elements */
         elem2.q = elem2.ecc = elem2.incl = elem2.arg_per = elem2.asc_node = 0.;
      add_sof_to_file( sofv_filename, &elem2, orbit + 6, n_obs, obs, "orbitdef.sof"); /* elem_ou2.cpp */
      }
   helio_elem = elem;            /* Heliocentric J2000 ecliptic elems */
   helio_elem.central_obj = 0;
   helio_elem.gm = SOLAR_GM;
   calc_classical_elements( &helio_elem, orbit2, epoch_shown, 1);
   tbuff = (char *)malloc( tbuff_size);
   n_lines = elements_in_mpc_format( tbuff, tbuff_size, &elem, object_name,
               is_cometary( constraints) && fabs( elem.ecc - 1.) < 1.e-6,
               output_format);
   fprintf( ofile, "%s\n", tbuff);
   tptr = tbuff + strlen( tbuff) + 1;
   *more_moids = '\0';
   for( i = 0; (size_t)i < sizeof( moids) / sizeof( moids[0]); i++)
      moids[i] = 0.;
   for( i = 1; i < n_lines && *tptr; i++)
      {
      char *tt_ptr;
      char sigma_buff[80];
      extern double uncertainty_parameter;
             /* "Solar radiation pressure at 1 AU",  in             */
             /* kg*AU^3 / (m^2*d^2),  from a private communication  */
             /* from Steve Chesley; see orb_func.cpp for details    */

      strlcpy_error( sigma_buff, "+/- ");
      strlcpy_error( buff, tptr);
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
                  snprintf_err( tt_ptr, 80, "; %s=%.4g", constraints, *mass);
                  consider_replacing( buff, constraints, "Sigma_mass:");
                  clobber_leading_zeroes_in_exponent( buff);
                  }
               else
                  snprintf_err( tt_ptr, 80, "; bad '%s'", constraints);
               }
            else
               snprintf_err( tt_ptr, 80, ";  Constraint: %s", constraints);
            }
         else if( n_monte_carlo_impactors && monte_carlo)
            snprintf_err( tt_ptr, 80, ";  %.2f%% impact (%d/%d)",
                100. * (double)n_monte_carlo_impactors /
                       (double)monte_carlo_object_count,
                       n_monte_carlo_impactors, monte_carlo_object_count);
         if( showing_sigmas)
            if( !get_uncertainty( "sigma_Tp", sigma_buff + 4, false))
               {
               strlcat_error( sigma_buff, " TT");
               text_search_and_replace( buff, "TT", sigma_buff);
               }

         if( (n_orbit_params >= 8 && n_orbit_params <= 10) || force_model == FORCE_MODEL_YARKO_A2)
            {
            int j = 0;

            if( !strcmp( constraints, "A=0") || !strcmp( constraints, "A1=0"))
               j = 1;      /* A1 constrained to be zero -> don't show it */
            strcat( tt_ptr, "\n");
            tt_ptr += strlen( tt_ptr);
            for( ; j < n_orbit_params - 6; j++)
               {
               char addenda[80];
               char tbuff0[80], sig_name[20];

               put_double_in_buff( tbuff0, orbit[j + 6]);
               text_search_and_replace( tbuff0, " ", "");
               strlcat_error( tbuff0, " ");
               snprintf_err( sig_name, sizeof( sig_name), "Sigma_A%d:", j + 1);
               if( showing_sigmas)
                  if( !get_uncertainty( sig_name, sigma_buff + 4, false))
                     format_value_with_sigma( tbuff0, orbit[j + 6], sigma_buff + 4);
               if( force_model == FORCE_MODEL_DELTA_V)
                  {
                  strlcpy_error( addenda, "T0: ");
                  if( j < 3)              /* dx, dy, or dz */
                     {
                     addenda[0] = 'd';
                     addenda[1] = 'x' + j;
                     }
                  else
                     {
                     full_ctime( tbuff0, orbit[9], CALENDAR_JULIAN_GREGORIAN
                               | FULL_CTIME_YMD | FULL_CTIME_LEADING_ZEROES | FULL_CTIME_MILLIDAYS);
                     snprintf_append( tbuff0, sizeof( tbuff0), " +/- %s", sigma_buff + 4);
                     }
                  strlcat_error( addenda, tbuff0);
                  if( j == 2)
                     strlcat_error( addenda, " m/s");
                  }
               else if( j == 3)
                  snprintf_err( addenda, sizeof( addenda), "DT: %s", tbuff0);
               else
                  {
                  const int n_to_use = (force_model == FORCE_MODEL_YARKO_A2 ? 2 : j + 1);

                  snprintf_err( addenda, sizeof( addenda), "A%d: %s", n_to_use, tbuff0);
                  if( j == n_orbit_params - 7 || j == 2)
                     {
                     strlcat_error( addenda, " AU/day^2");
                     if( is_inverse_square_force_model( ))
                        strlcat_error( addenda, " [1/r^2]");
                     }
                  }
               if( strlen( addenda) + strlen( tt_ptr) > 74)
                  {
                  strlcat_err( tt_ptr, "\n", 79);
                  tt_ptr += strlen( tt_ptr);
                  }
               else if( *tt_ptr)
                  strlcat_err( tt_ptr, "   ", 79);
               strlcat_err( tt_ptr, addenda, 79);
               }
            }
         assert( strlen( buff) < sizeof( buff) - 1);
         }

      else if( !memcmp( buff, "Epoch", 5))
         {
         int j;
         int n_moids_to_show = atoi( get_environment_ptr( "N_MOIDS"));

         if( !n_moids_to_show)      /* by default,  just show */
            n_moids_to_show = 8;    /* planetary MOIDs */
         if( force_model == FORCE_MODEL_SRP)
            {
            snprintf_err( tt_ptr, 80, "; AMR %.5g",
                                 orbit[6] * SOLAR_GM / SRP1AU);
            if( showing_sigmas)
               consider_replacing( tt_ptr, "AMR", "Sigma_AMR:");
            strlcat_err( tt_ptr, " m^2/kg", 80);
            }
         for( j = 0; j < n_moids_to_show; j++)
            {
            static const char moid_idx[N_MOIDS] = { 3, 5, 2, 1, 4, 6, 7, 8,
                                       10, 11, 12, 13, 14, 15 };
            const int planet = moid_idx[j];

            if( !planet_orbiting || planet_orbiting == planet)
               {
               moid_data_t mdata;
               double moid, moid_limit = .1;
               ELEMENTS planet_elem;
               const int forced_moid =
                   (atoi( get_environment_ptr( "MOIDS")) >> j) & 1;

               setup_planet_elem( &planet_elem, planet,
                                (epoch_shown - J2000) / 36525.);
               moid = find_moid_full( &planet_elem, &helio_elem, &mdata);
               if( !j)     /* only get Barbee speed for the earth */
                  {
                  double lunar_loc[3], earth_loc[3], lunar_r, earth_r;
                  double v2;

                  earth_lunar_posn( epoch_shown, earth_loc, lunar_loc);
                  lunar_r = sqrt( vect_diff2( orbit2, lunar_loc));
                  earth_r = sqrt( vect_diff2( orbit2, earth_loc));
                  v2 = mdata.barbee_speed * mdata.barbee_speed
                                    - 2. * get_planet_mass( 3) / earth_r
                                    - 2. * get_planet_mass( 10) / lunar_r;
                  if( v2 > 0.)
                     barbee_style_delta_v = sqrt( v2);
                  else
                     barbee_style_delta_v = -sqrt( -v2);
                  barbee_style_delta_v *= AU_IN_KM / seconds_per_day;
                  }
               if( j < 2)        /* Earth or Jupiter */
                  moid_limit = 1.;
               else if( j > 7)            /* asteroid */
                  moid_limit = .1;
               else if( j > 4)          /* Saturn,  Uranus,  Neptune */
                  moid_limit = 1.;
               moids[planet] = moid;
               if( !planet_orbiting && (forced_moid || moid < moid_limit))
                  {
                  char addendum[30];
                  snprintf_err( addendum, sizeof( addendum),
                                 "   %s: %.4f", moid_text[j], moid);
                  if( strlen( addendum) + strlen( buff) < 79)
                     strlcat_error( buff, addendum);
                  else
                     {
                     if( n_more_moids < 3)
                        strlcat_error( more_moids, addendum);
                     n_more_moids++;
                     }
                  if( !j && moid < .5)
                     snprintf_append( orbit_summary_text,
                           sizeof( orbit_summary_text), " MOID %.3f", moid);
                  }
               }
            }
         }
      else if( *more_moids)
         {
         space_pad_buffer( buff, 33);
         strlcpy_err( buff + 33, more_moids, sizeof( buff) - 33);
         *more_moids = '\0';
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
            snprintf_err( buff + 32, 9, "U%5.1f  ", uncertainty_parameter);
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
            else if( *buff == 'z')
               consider_replacing( buff, "z", "sigma_1/a");
            }
         switch( buff[0])
            {
            case 'M':
               if( buff[1] != '(')  /* if the 'M' doesn't refer to a comet mag */
                  consider_replacing( buff, "M", "sigma_M");
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
                  strlcpy_error( tbuff, zptr);
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
                     memcpy( phg_line + 29, buff + 23, 5);  /* move H */
                     memcpy( phg_line + 34, "  G ", 5);
                     memcpy( phg_line + 38, buff + 35, 5);  /* move G */
                     phg_line[43] = '\0';
                     }
                  else
                     strlcpy_error( phg_line, buff);

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
                  strlcat_error( buff, "    ");
                  strlcat_error( buff, tbuff);
                  }
               }
               break;
            default:
               break;
            }
         }
      if( !reference_shown)
         reference_shown = show_reference( buff);
      assert( *body_frame_note);
      if( !body_frame_note_shown)
         {
         size_t j = strlen( buff);

         if( j < 30)
            {
            while( j < 36)
               buff[j++] = ' ';
            strlcpy_err( buff + j, body_frame_note, sizeof( buff) - j);
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
               strlcpy_err( buff + 58, body_frame_note, sizeof( buff) - 58);
               body_frame_note_shown = true;
               }        /* above basically says,  "if we haven't gotten */
            }           /* the body frame note in the 'normal' places, try */
         }              /* wedging it in at the end" */
      fprintf( ofile, "%s\n", buff);
      tptr += strlen( tptr) + 1;
      }
   observation_summary_data( tbuff, obs, n_obs, options);
   fprintf( ofile, "%s\n", tbuff);
   if( elem.central_obj == 3 && elem.ecc < .99 && _include_comment( "TLE"))
      if( !write_tle_from_vector( tbuff, rel_orbit, elem.epoch, NULL, NULL))
         {
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
      if( _include_comment( "State"))
         {
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
         }
      if( !planet_orbiting && _include_comment( "MOID"))
         {
         fprintf( ofile, "# MOIDs: Me%10.6f Ve%10.6f Ea%10.6f Ma%10.6f\n",
                  moids[1], moids[2], moids[3], moids[4]);
         fprintf( ofile, "# MOIDs: Ju%10.6f Sa%10.6f Ur%10.6f Ne%10.6f\n",
                  moids[5], moids[6], moids[7], moids[8]);
         if( moids[10])
            {
            fprintf( ofile, "# MOIDs:");
            for( i = 8; i < N_MOIDS && moids[i + 2]; i++)
               fprintf( ofile, " %s%10.6f", moid_text[i], moids[i + 2]);
            fprintf( ofile, "\n");
            }
         }

      if( !planet_orbiting && _include_comment( "OPP"))
         fprintf( ofile, "# %d oppositions\n", _n_oppositions( obs, n_obs));

      geocentric_score = geo_score( orbit, curr_epoch);
      if( geocentric_score > 0)
         fprintf( ofile, "# GEO score %d\n", geocentric_score);
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
      if( _include_comment( "Twrit"))
         fprintf( ofile, "# Elements written: %s (JD %f)\n", buff, jd);
      make_date_range_text( buff, obs[0].jd, obs[n_obs - 1].jd);
      if( _include_comment( "ObsRange"))
         fprintf( ofile, "# Full range of obs: %s (%d observations)\n",
                              buff, n_obs);
      if( _include_comment( "Ver"))
         fprintf( ofile, "# Find_Orb ver: %s\n", find_orb_version_jd( NULL));
      if( _include_comment( "Per"))
         {
         fprintf( ofile, "# Perturbers: %08lx ", (unsigned long)perturbers);
         if( !perturbers)
            fprintf( ofile, "(unperturbed orbit)");
         else if( (perturbers & 0x3fe) == 0x3fe)
            fprintf( ofile, (perturbers & 0x400) ? "(Merc-Pluto plus Luna)" :
                  "(Merc-Pluto, Earth & moon combined)");
         else if( perturbers == 0x408)
            fprintf( ofile, "(Sun/Earth/Moon)");
         }
      get_jpl_ephemeris_info( &jpl_de_version, NULL, NULL);
      if( _include_comment( "DE"))
         {
         if( jpl_de_version)
            fprintf( ofile, ";  JPL DE-%d\n", jpl_de_version);
         else
            fprintf( ofile, ";  not using JPL DE\n");
         }

      if( !elem.central_obj)
         {
         if( elem.ecc < 1. && _include_comment( "Tis"))
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
         if( helio_elem.q < 1.04 && _include_comment( "Vin"))
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
      if( barbee_style_delta_v && helio_elem.q < 1.3 && _include_comment( "Vin"))
         fprintf( ofile, "# Barbee-style encounter velocity: %.4f km/s\n",
                              barbee_style_delta_v);
      if( elem.abs_mag && elem.is_asteroid)
         {
         const double diam = diameter_from_abs_mag( elem.abs_mag, optical_albedo);

         fprintf( ofile, "# Diameter %.1f %s (assuming %.0f%% albedo)\n",
               (diam > 10000. ? diam / 1000. : diam),
               (diam > 10000. ? "km" : "meters"), optical_albedo * 100.);
         }

      if( _include_comment( "Sco"))
         fprintf( ofile, "# Score: %f\n", evaluate_initial_orbit( obs, n_obs, orbit, curr_epoch));
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
         const double impact_time_utc = utc_from_td( impact_time_td, NULL);

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
            snprintf_err( end_ptr, 11, "%c%.5f",
                  (lon < 180. ? 'E' : 'W'),
                  (lon < 180. ? lon : 360. - lon));
            fprintf( ofile, "%s at %s\n", (is_an_impact ? "IMPACT" : "LAUNCH"),
                                   impact_buff);
            }
                     /* Then show in 0-360 format,  for all other  */
                     /* planets, and for output to file:           */
         snprintf_err( end_ptr, 10, "%9.5f", lon);
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
   if( !(options & ELEM_OUT_NO_COMMENT_DATA) && _include_comment( "Pse"))
      {
      char time_buff[40];
      bool nongrav_sigmas_found = false;
      double output_a;

      snprintf_err( tbuff, 80, "#  $Name=%s", object_name);
      text_search_and_replace( tbuff + 4, " ", "%20");
               /* Epoch has to be given in YYYYMMDD.DDDDD format: */
      full_ctime( time_buff, helio_elem.perih_time,
               FULL_CTIME_YMD | FULL_CTIME_MONTHS_AS_DIGITS
               | FULL_CTIME_MICRODAYS | FULL_CTIME_LEADING_ZEROES);
      time_buff[4] = time_buff[7] = '\0';
      snprintf_append( tbuff, tbuff_size, "  $Ty=%s  $Tm=%s  $Td=%s",
               time_buff, time_buff + 5, time_buff + 8);
      snprintf_append( tbuff, tbuff_size, "  $MA=%.5f",
                  centralize_ang( helio_elem.mean_anomaly) * 180. / PI);
      fprintf( ofile, "%s\n", tbuff);

      snprintf_err( tbuff, 80, "#  $ecc=%.7f  $Eqnx=2000.", helio_elem.ecc);
      fprintf( ofile, "%s\n", tbuff);

      output_a = helio_elem.major_axis;
      if( fabs( output_a) > 1000000.)
         output_a = 1000000.;       /* near-parabolic case */
      snprintf_err( tbuff, 80, "#  $a=%.7f  $Peri=%.5f  $Node=%.5f",
                  output_a,
                  centralize_ang( helio_elem.arg_per) * 180. / PI,
                  centralize_ang( helio_elem.asc_node) * 180. / PI);
      snprintf_append( tbuff, tbuff_size, "  $Incl=%.5f",
                  helio_elem.incl * 180. / PI);
      fprintf( ofile, "%s\n", tbuff);

      snprintf_err( tbuff, 80, "#  $EpJD=%.3f  $q=%.6f", helio_elem.epoch, helio_elem.q);
      snprintf_append( tbuff, tbuff_size, "  $T=%.6f  $H=%.1f",
               helio_elem.perih_time, helio_elem.abs_mag);
      fprintf( ofile, "%s\n", tbuff);
      fprintf( ofile, "# Sigmas avail: %d\n", available_sigmas);
      for( i = 1; i <= 3; i++)
         {
         char key[20], obuff[50];

         snprintf_err( key, sizeof( key), "Sigma_A%d:", i);
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

   if( *horizons_elems)
      write_horizons_elems( horizons_elems, &helio_elem, orbit);
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
      snprintf_err( name_buff, sizeof( name_buff), "%05d", n_clones_accepted);
      packed_desig_minus_spaces( virtual_full_desig, obs->packed_id);
      snprintf_append( virtual_full_desig, sizeof( virtual_full_desig), " [%d]",
                                  monte_carlo_object_count);
      if( elem.central_obj || elem.ecc > .999999)
         if( !elements_in_guide_format( tbuff, &elem, virtual_full_desig, obs, n_obs))
            {
            ofile = fopen_ext( get_file_name( buff, "virtual.txt"), monte_carlo_permits);
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

            fprintf( ofile, "Monte Carlo orbits from Find_Orb\nComputed %.24s\n",
                              asctime( gmtime( &t0)));
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

   if( !elements_in_guide_format( tbuff, &elem, object_name, obs, n_obs)
         && (ofile = fopen_ext( get_file_name( buff, "guide.txt"), "tfcwb")) != NULL)
      {
      fprintf( ofile, "%s%s\n", tbuff, impact_buff);
      fclose( ofile);
      }

   ofile = open_json_file( tbuff, "JSON_ELEMENTS_NAME", "elements.json", obs->packed_id,
                     "wb");
   assert( ofile);
   elements_in_json_format( ofile, &elem, orbit, object_name, obs, n_obs, moids,
                                             body_frame_note, true, geocentric_score);
   fclose( ofile);
   ofile = open_json_file( tbuff, "JSON_SHORT_ELEMENTS", "elem_short.json", obs->packed_id,
                     "wb");
   assert( ofile);
   elements_in_json_format( ofile, &elem, orbit, object_name, obs, n_obs, moids,
                                             body_frame_note, false, geocentric_score);
   fclose( ofile);
   make_linkage_json( n_obs, obs, &elem);
   free( tbuff);
   available_sigmas = saved_available_sigmas;
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
   size_t n_lines, name_len = 0, i;
   char **ilines = load_file_into_memory( "orbits.sof", &n_lines, false);
   int sort_column = 0;

   for( i = 0; i < (size_t)n_ids; i++)
      ids[i].solution_exists = 0;
   if( !ilines)
      return;
   while( ilines[0][name_len] && ilines[0][name_len] != '|')
      name_len++;
   assert( ilines[0][name_len] == '|');
   for( i = 0; i < (size_t)n_lines; i++)
      {
      ilines[i][name_len] = '\0';
      remove_trailing_cr_lf( ilines[i]);
      if( strlen( ilines[i]) == 12)    /* combined perm & provisional ID */
         ilines[i][5] = '\0';
      }
   shellsort_r( ilines, n_lines, sizeof( char *),
                           string_compare_for_sort, &sort_column);
   for( i = 0; i < (size_t)n_ids; i++)
      {
      char tname[15], *tptr = tname;

      assert( strlen( ids[i].packed_desig) == 12);
      strcpy( tname, ids[i].packed_desig);
      text_search_and_replace( tname, " ", "");
      if( strlen( tname) == 12)    /* combined perm & provisional ID */
         tname[5] = '\0';
      if( bsearch_ext_r( &tptr, ilines, n_lines, sizeof( char *),
                        string_compare_for_sort, &sort_column, NULL))
         ids[i].solution_exists = 1;
      }
   free( ilines);
}

static void compute_two_body_state_vect( ELEMENTS *elems, double *orbit, const double jd)
{
   derive_quantities( elems, get_planet_mass( elems->central_obj));
   comet_posn_and_vel( elems, jd, orbit, orbit + 3);
   if( elems->central_obj == 3)
      {
      equatorial_to_ecliptic( orbit);
      equatorial_to_ecliptic( orbit + 3);
      }
   if( elems->central_obj)
      {
      double planet_state[6];
      int i;

      get_planet_posn_vel( jd, elems->central_obj,
                              planet_state, planet_state + 3);
      for( i = 0; i < 6; i++)
         orbit[i] += planet_state[i];
      }
}

/* Can get some comet elements from

   http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET

   This is especially helpful for pre-2008 SOHO objects.  Without
ELEMENTS.COMET,  Find_Orb can sometimes flounder about a bit in its
efforts to determine an orbit. */

static int get_orbit_from_sof( const char *filename,
                 const char *object_name, double *orbit, ELEMENTS *elems,
                 const double full_arc_len, double *max_resid)
{
   FILE *ifile = fopen_ext( filename, "crb");
   int got_vectors = 0;

   memset( elems, 0, sizeof( ELEMENTS));
   if( ifile)
      {
      char buff[300], header[300], tname[15];
      int n_extra_params = 0;

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
         if( !memcmp( tname, buff, 12) && !_get_extra_orbit_info(
                       tname, &perturbers, &n_extra_params,
                       orbit + 6, NULL))
            {
            double extra_info[10];
            extern int n_orbit_params, force_model;

            extract_sof_data_ex( elems, buff, header, extra_info);
            if( elems->epoch < 2400000.)
               printf( "JD %f\n", elems->epoch);
            assert( elems->epoch > 2400000.);
            perturbers = 0x7fe;    /* Merc-Pluto plus moon */
            n_orbit_params = n_extra_params + 6;
            compute_two_body_state_vect( elems, orbit, elems->epoch);
            push_orbit( elems->epoch, orbit);
            if( extra_info[1] - extra_info[0] > full_arc_len / 2.)
               {
               got_vectors = 1;
               if( elems->ecc == 1.)     /* indicate parabolic-constraint orbit */
                  got_vectors = 2;
               *max_resid = extra_info[3];
               }
            else
               {
               force_model = 0;
               n_orbit_params = 6;
               }
            }
      fclose( ifile);
      }
   return( got_vectors);
}

int compute_canned_object_state_vect( double *loc, const char *mpc_code,
                     const double jd)
{
   char buff[300], header[300], *tptr;
   FILE *ifile;
   size_t tp_offset;
   double prev_jd_read = -1e+99;   /* 'infinitely long ago' */
   int rval = -2;

   snprintf( buff, sizeof( buff), "eph_%s.txt", mpc_code);
   ifile = fopen_ext( buff, "fcrb");
   if( !ifile)
      debug_printf( "File '%s' not found\n", buff);
   assert( ifile);
   if( !fgets_trimmed( header, sizeof( header), ifile))
      {
      fprintf( stderr, "Error in %s header\n", buff);
      exit( -1);
      }
   tptr = strstr( header, "|Te");
   assert( tptr);
   tp_offset = tptr - header + 1;
   while( rval && fgets( buff, sizeof( buff), ifile))
      {
      const double curr_jd_read = extract_yyyymmdd_to_jd( buff + tp_offset);

      if( curr_jd_read > jd)
         {
         rval = 0;
         if( jd < (curr_jd_read + prev_jd_read) / 2.)
            {            /* should have used previous line */
            fseek( ifile, -2 * (long)strlen( buff), SEEK_CUR);
            if( !fgets( buff, sizeof( buff), ifile))
               exit( -1);
            }
         }
      else if( curr_jd_read > 0.)
         rval = -1;
      prev_jd_read = curr_jd_read;
      }
   if( rval == -1)      /* read to end of file,  extrapolating */
      rval = 0;
   if( !rval)
      {
      ELEMENTS elems;

      memset( &elems, 0, sizeof( ELEMENTS));
      extract_sof_data_ex( &elems, buff, header, NULL);
      if( elems.epoch < 2400000.)
         debug_printf( "JD %f\n", elems.epoch);
      assert( elems.epoch > 2400000.);
      compute_two_body_state_vect( &elems, loc, jd);
      }
   fclose( ifile);
   return( rval);
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

#define FOUND_ECCENTRICITY     1
#define FOUND_TPERIH           2
#define FOUND_Q                4
#define FOUND_A                8
#define FOUND_MEAN_ANOMALY    16

static double extract_state_vect_from_text( const char *text,
            double *orbit, double *abs_mag)
{
   char tbuff[81];
   size_t i = 0;
   double epoch = 0., coord_epoch = 2000.;
   double dist_units = 1., time_units = 1.;
   int bytes_read, central_object = 0, quantities_found = 0;
   bool is_state_vector, is_ecliptic = true;
   ELEMENTS elem;

   while( text[i] && text[i] != ',' && i < 80)
      i++;
   if( text[i] == ',')
      {
      memcpy( tbuff, text, i);
      tbuff[i] = '\0';
      epoch = get_time_from_string( 0., tbuff, CALENDAR_JULIAN_GREGORIAN, NULL);
      }
   assert( epoch);
   text += i + 1;
   for( i = 0; i < 6 && sscanf( text, "%lf%n", orbit + i, &bytes_read) == 1; i++)
      text += bytes_read + 1;
   assert( i == 6 || i == 0);
   is_state_vector = (i == 6);
   memset( &elem, 0, sizeof( ELEMENTS));
   *abs_mag = 18.;
   while( *text)
      {
      if( !memcmp( text, "eq", 2))
         is_ecliptic = false;
      else if( !memcmp( text, "ep=", 3))
         coord_epoch = atof( text + 3);
      else if( !memcmp( text, "H=", 2))
         *abs_mag = atof( text + 2);
      else if( !memcmp( text, "Om=", 3))
         elem.asc_node = atof( text + 3) * PI / 180.;
      else if( !memcmp( text, "om=", 3))
         elem.arg_per = atof( text + 3) * PI / 180.;
      else if( !memcmp( text, "i=", 2))
         elem.incl = atof( text + 2) * PI / 180.;
      else if( !memcmp( text, "Tp=", 3))
         {
         i = 3;
         while( i < 70 && text[i] != ',' && text[i])
            i++;
         assert( i < 70);
         memcpy( tbuff, text + 3, i - 3);
         tbuff[i - 3] = '\0';
         elem.perih_time = get_time_from_string( 0., tbuff, CALENDAR_JULIAN_GREGORIAN, NULL);
         quantities_found |= FOUND_TPERIH;
         }
      else if( !memcmp( text, "M=", 2))
         {
         elem.mean_anomaly = atof( text + 2) * PI / 180.;
         quantities_found |= FOUND_MEAN_ANOMALY;
         }
      else if( !memcmp( text, "a=", 2))
         {
         elem.major_axis = atof( text + 2);
         quantities_found |= FOUND_A;
         }
      else if( !memcmp( text, "q=", 2))
         {
         elem.q = atof( text + 2);
         quantities_found |= FOUND_Q;
         }
      else if( !memcmp( text, "e=", 2))
         {
         elem.ecc = atof( text + 2);
         quantities_found |= FOUND_ECCENTRICITY;
         }
      else if( !_memicmp( text, "UTC", 3))
         epoch += td_minus_utc( epoch) / seconds_per_day;
      else if( !_memicmp( text, "TDB", 3))
         {
         const double t_cen = (epoch - J2000) / 36525.;

         epoch -= tdb_minus_tdt( t_cen) / seconds_per_day;
         }
      else if( !memcmp( text, "km", 2))
         dist_units = AU_IN_KM;
      else if( *text == 'm')
         dist_units = AU_IN_METERS;
      else if( *text == 'g')
         central_object = 3;
      else if( *text == 's')
         time_units = seconds_per_day;
      else if( *text == 'A' && text[1] >= '1' && text[1] <= '3'
                  && text[2] == '=')
         {
         extern int n_orbit_params;

         n_orbit_params = 6 + text[1] - '1';
         force_model = ((n_orbit_params - 6) | 0x10);   /* assume inverse square for the nonce */
         orbit[n_orbit_params - 7] = atof( text + 3);
         }
      while( *text != ',' && *text)
         text++;
      if( *text == ',')
         text++;
      }
   if( !is_state_vector)
      {
      elem.epoch = epoch;
      if( !(quantities_found & FOUND_TPERIH))
         elem.perih_time = elem.epoch;
      if( quantities_found & FOUND_A)
         elem.q = elem.major_axis * (1. - elem.ecc);
      elem.q /= dist_units;
      elem.gm = get_planet_mass( central_object);
      derive_quantities( &elem, elem.gm);
      if( quantities_found & FOUND_MEAN_ANOMALY)
         elem.perih_time = elem.epoch - elem.mean_anomaly * elem.t0;
      comet_posn_and_vel( &elem, elem.epoch, orbit, orbit + 3);
      }
   else
      {
      for( i = 0; i < 6; i++)
         orbit[i] /= dist_units;
      for( i = 3; i < 6; i++)
         orbit[i] *= time_units;
      }
   if( !is_ecliptic)
      {
      equatorial_to_ecliptic( orbit);
      equatorial_to_ecliptic( orbit + 3);
      text += 2;
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

int find_first_and_last_obs_idx( const OBSERVE *obs, const int n_obs,
         int *last)
{
   int i;

   if( last)
      {
      for( i = n_obs - 1; i > 0 && !obs[i].is_included; i--)
         ;
      *last = i;
      }
   for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
      ;
   return( i);
}

int get_orbit_from_mpcorb_sof( const char *object_name, double *orbit,
             ELEMENTS *elems, const double full_arc_len, double *max_resid)
{
   const char *mpcorb_dot_sof_filename = get_environment_ptr( "MPCORB_SOF_FILENAME");
   int rval = 0;

   if( *mpcorb_dot_sof_filename)
      rval = get_orbit_from_sof( mpcorb_dot_sof_filename,
                             object_name, orbit, elems, full_arc_len, max_resid);
   return( rval);
}

/* When doing (for example) a full six-parameter fit to an orbit,  it can
be helpful to use an epoch that is at the mid-point of the arc of
observations that is being fitted.  This improves stability.  */

double mid_epoch_of_arc( const OBSERVE *obs, const int n_obs)
{
   int first, last;

   first = find_first_and_last_obs_idx( obs, n_obs, &last);
   return( (obs[first].jd + obs[last].jd) / 2.);
}

static void _reset_sr_orbits( void)
{
   extern double *sr_orbits;
   extern unsigned max_n_sr_orbits;

   if( sr_orbits)
      free( sr_orbits);
   sr_orbits = (double *)calloc( (size_t)max_n_sr_orbits, 7 * sizeof( double));
   assert( sr_orbits);
}

const char *state_vect_text = NULL;
int ignore_prev_solns;
bool take_first_soln = false, force_final_full_improvement = false;
int n_extra_full_steps = 0;

static int fetch_previous_solution( OBSERVE *obs, const int n_obs, double *orbit,
               double *orbit_epoch, double *epoch_shown)
{
   int got_vectors = 0, i;
   extern int n_orbit_params;
   char object_name[80];
   bool do_full_improvement = false;
   double abs_mag = 10.;         /* default value for 'dummy' use */

   assert( n_obs > 0);
   get_object_name( object_name, obs->packed_id);
   n_orbit_params = 6;
   if( state_vect_text)    /* state vect supplied on cmd line w/ '-v' */
      {
      *orbit_epoch = extract_state_vect_from_text(
                  state_vect_text, orbit, &abs_mag);
      got_vectors = (*orbit_epoch != 0.);
      if( got_vectors)
         push_orbit( *orbit_epoch, orbit);
      }
   if( !got_vectors && !ignore_prev_solns)
      {
      ELEMENTS elems;
      const double full_arc_len = obs[n_obs - 1].jd - obs[0].jd;
      double rms_resid;

      got_vectors = get_orbit_from_sof( "orbits.sof", obs->packed_id,
                                          orbit, &elems, full_arc_len, &rms_resid);
#ifdef FALLING_BACK_ON_MPCORB          /* MPCORB isn't entirely reliable */
      if( !got_vectors)
         got_vectors = get_orbit_from_mpcorb_sof( object_name, orbit, &elems,
                                  full_arc_len, &rms_resid);
#endif
      if( got_vectors)
         {
         move_add_nstr( 2, 2, "Stored soln loaded", -1);
         refresh_console( );
         *orbit_epoch = elems.epoch;
         if( got_vectors == 1)
            perturbers = 0x7fe;    /* Merc-Pluto plus moon */
         if( n_obs > 1)
            {
            set_locs( orbit, *orbit_epoch, obs, n_obs);
            do_full_improvement = true;
            }
         abs_mag = elems.abs_mag;
         }
      }
   if( !got_vectors)
      {
      perturbers = 0;
      *orbit_epoch = initial_orbit( obs, n_obs, orbit);
      if( force_final_full_improvement)
         do_full_improvement = true;
      }
   else if( n_obs == 1 && !strcmp( obs->reference, "Dummy"))
      {
      obs->jd = *orbit_epoch + td_minus_ut( *orbit_epoch) / seconds_per_day;
      perturbers = 0x7fe;    /* Merc-Pluto plus moon */
      set_locs( orbit, *orbit_epoch, obs, n_obs);
      obs->ra = obs->computed_ra;
      obs->dec = obs->computed_dec;
      obs->obs_mag = abs_mag + calc_obs_magnitude( obs->solar_r,
                 obs->r, vector3_length( obs->obs_posn), NULL);
      obs->obs_mag = floor( obs->obs_mag * 10. + .5) * .1;
      obs->mag_precision = 1;         /* start out assuming mag to .1 mag */
      }
   for( i = 0; i < n_obs; i++)
      if( obs[i].note2 == 'R')
         do_full_improvement = true;
   *epoch_shown = *orbit_epoch;
   assert( n_obs > 0 && n_obs < 100000);
   if( do_full_improvement || available_sigmas == NO_SIGMAS_AVAILABLE)
      {
      extern double automatic_outlier_rejection_limit;
      OBSERVE *saved_obs = (OBSERVE *)calloc( n_obs, sizeof( OBSERVE));
      double prev_score = 0.;
      int pass;

      if( !*get_environment_ptr( "KEEP_PREVIOUS_EPOCH") && !got_vectors)
         *epoch_shown = find_epoch_shown( obs, n_obs);
      memcpy( saved_obs, obs, n_obs * sizeof( OBSERVE));
      filter_obs( obs, n_obs, automatic_outlier_rejection_limit, 0);
      for( i = 0; i < n_obs; i++)
         if( obs[i].note2 == 'R')
            obs[i].is_included = 1;
      if( got_vectors)
         attempt_extensions( obs, n_obs, orbit, *orbit_epoch);
      for( pass = 0; pass < 2; pass++)
         {
         const int64_t t0 = nanoseconds_since_1970( );
         const int64_t QUARTER_SECOND = 250000000;
         const double mid_epoch = mid_epoch_of_arc( obs, n_obs);

         push_orbit( *orbit_epoch, orbit);
         integrate_orbit( orbit, *orbit_epoch, mid_epoch);
         *orbit_epoch = mid_epoch;
         if( !pass)
            prev_score = evaluate_initial_orbit( obs, n_obs, orbit, *orbit_epoch);
         for( i = 0; i < 4 && (nanoseconds_since_1970( ) - t0) < QUARTER_SECOND; i++)
            {
            if( i)
               filter_obs( obs, n_obs, automatic_outlier_rejection_limit, 0);
            full_improvement( obs, n_obs, orbit, *orbit_epoch,
                           (pass ? "e=1" : NULL),
                           ORBIT_SIGMAS_REQUESTED, *epoch_shown);
            }
         if( !pass)
            for( i = n_extra_full_steps; i; i--)
               {
               if( i > 1)
                  filter_obs( obs, n_obs, automatic_outlier_rejection_limit, 0);
               full_improvement( obs, n_obs, orbit, *orbit_epoch,
                              (pass ? "e=1" : NULL),
                              ORBIT_SIGMAS_REQUESTED, *epoch_shown);

               }
         if( prev_score < evaluate_initial_orbit( obs, n_obs, orbit, *orbit_epoch) - .001)
            {
            pop_orbit( orbit_epoch, orbit);    /* we were better off with the old orbit */
            memcpy( obs, saved_obs, n_obs * sizeof( OBSERVE));
            }
         else        /* throw out the saved orbit;  we've got something better */
            {
            pass = 100;
            pop_orbit( NULL, NULL);
            }
         }
      free( saved_obs);
      }
               /* if a stored solution failed (i.e.,  didn't get sigmas), */
               /* we try again,  ignoring the stored solution.            */
   if( got_vectors && available_sigmas == NO_SIGMAS_AVAILABLE && !ignore_prev_solns)
      {
      move_add_nstr( 3, 2, "Stored soln FAILED", -1);
      ignore_prev_solns = 1;
      n_orbit_params = 6;
      force_model = 0;
      for( i = 0; i < n_obs; i++)
         obs[i].is_included = !(obs[i].flags & OBS_DONT_USE);
      got_vectors = fetch_previous_solution( obs, n_obs, orbit, orbit_epoch,
                        epoch_shown);
      ignore_prev_solns = 0;
      }
   else if( !got_vectors && available_sigmas == NO_SIGMAS_AVAILABLE && !ignore_prev_solns)
      {                          /* failed to get sigmas;  try more SR orbits */
      extern unsigned max_n_sr_orbits;

      ignore_prev_solns = 1;
      debug_printf( "%s failed; trying more SR orbits\n", obs->packed_id);
      max_n_sr_orbits *= 5;
      _reset_sr_orbits( );
      got_vectors = fetch_previous_solution( obs, n_obs, orbit, orbit_epoch,
                        epoch_shown);
      max_n_sr_orbits /= 5;
      _reset_sr_orbits( );
      ignore_prev_solns = 0;
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

double get_max_included_resid( const OBSERVE *obs, int n_obs)
{
   double rval = 0.;

   while( n_obs--)
      {
      if( obs->is_included)
         {
         const double resid = observation_rms( obs);

         if( rval < resid)
            rval = resid;
         }
      obs++;
      }
   return( rval);
}

static void _log_problems( const OBJECT_INFO *id, const OBSERVE FAR *obs)
{
   int i, n_obs_used = 0, n_real_obs = 0;
   double first_jd = -1., last_jd = obs->jd, arc_used, full_arc;
   const double rms_err = compute_weighted_rms( obs, id->n_obs, NULL);
   char error_message[200];

   error_message[0] = '\0';
   for( i = 0; i < id->n_obs; i++, obs++)
      if( obs->is_included)
         {
         if( first_jd < 0.)
            first_jd = obs->jd;
         last_jd = obs->jd;
         n_obs_used++;
         if( !(obs->flags & OBS_DONT_USE))
            n_real_obs++;
         }
   if( n_real_obs < 2)        /* skip singleton cases */
      return;
   obs -= id->n_obs;
   if( n_obs_used <= n_real_obs * 2 / 3)
      snprintf_append( error_message, sizeof( error_message),
               " %d of %d observations used", n_obs_used, n_real_obs);
   arc_used = last_jd - first_jd;
   full_arc = obs[id->n_obs - 1].jd - obs[0].jd;
   if(  arc_used < full_arc / 2.)
      snprintf_append( error_message, sizeof( error_message),
               " %f/%f of arclength used", arc_used, full_arc);
   if( available_sigmas == NO_SIGMAS_AVAILABLE)
      snprintf_append( error_message, sizeof( error_message),
               " no sigmas");
   if( rms_err > 3.)
      snprintf_append( error_message, sizeof( error_message),
               " rms=%f", rms_err);
   if( *error_message)
      {
      FILE *ofile = fopen_ext( "errors.txt", "ca");

      if( ofile)
         {
         fprintf( ofile, "%s : %s\n", id->packed_desig, error_message);
         fclose( ofile);
         }
      }
}

OBSERVE FAR *load_object( FILE *ifile, OBJECT_INFO *id,
                       double *curr_epoch, double *epoch_shown, double *orbit)
{
   extern int n_obs_actually_loaded;
   extern int debug_level;
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
   force_model = 0;
   obj_desig_to_perturber( id->packed_desig);
   if( id->obj_name[0] == '(')    /* numbered asteroid:  shouldn't perturb */
      excluded_asteroid_number = atoi( id->obj_name + 1);      /* itself */
   else
      excluded_asteroid_number = 0;
   compute_effective_solar_multiplier( NULL);
   if( n_obs_actually_loaded > 0)
      {
      const int got_vector = fetch_previous_solution( obs, id->n_obs, orbit,
                            curr_epoch, epoch_shown);

      if( got_vector <= 0)
         *epoch_shown = find_epoch_shown( obs, id->n_obs);
      _log_problems( id, obs);
      }
   else                           /* indicate failure */
      *epoch_shown = *curr_epoch = 0.;
   return( obs);
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
      double phi1, phi2, log_tan_half_phase, tval;

      if( ph_ang < 1e-10)   /* "close enough" to zero */
         phi1 = phi2 = 1.;
      else
         {
         log_tan_half_phase = log( sin( ph_ang / 2.) / cos( ph_ang / 2.));
         phi1 = exp( -3.33 * exp( log_tan_half_phase * 0.63));
         phi2 = exp( -1.87 * exp( log_tan_half_phase * 1.22));
         }
      tval = (1. - asteroid_magnitude_slope_param) * phi1
                + asteroid_magnitude_slope_param * phi2;
      if( tval < 1e-50)    /* avoid domain errors;  results in _very_ */
         rval = 9999.;     /* faint magnitudes */
      else
         rval = 5. * log( obj_sun) - 2.5 * log( tval);
      }
   rval += 5. * log( obj_earth);
   rval /= LOG_10;         /* allow for common logs,  not naturals */
   return( rval);
}

int generic_message_box( const char *message, const char *box_type);

/* Almost all mag band shifts,  as of 2021 Sep 19,  come from

https://minorplanetcenter.net/iau/info/BandConversion.txt

Previous versions drew from a variety of sources;  see old versions
for comments (many of which would probably apply here as well,  but
the MPC table lacks references.)

The following function,  as the comment indicates,  assumes that
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
either 2MASS colors to estimate a color correction or just assuming
an average color correction.

'C' ('clear',  'unfiltered') magnitudes are something of a problem.
They never really should have existed.  As described above,  the key
thing isn't what filter you use (though that's a nice thing to
know);  it's what sort of reference photometry you used.  At
present,  'C' magnitudes are treated as close enough to 'R',  per
MPC practice.  (I perhaps should treat them as 'probably wrong' and
give them minimal or no weight in the computation of H values.)

Note that along with the magnitude band shift,  the mag sigma
ought to be adjusted,  with some quantity added in quadrature.
I _think_ that quantity would be small for transforming,  say,
R to V,  but would be large for transforming a "less V-like"
band (such as the IR bands JHK) to V.     */

typedef struct
{
   char band;
   double v_mag_correction;
} band_shift_t;

double mag_band_shift( const char mag_band, int *err_code)
{
   size_t i;
   static bool unknown_band_warning_shown = false;
   static const band_shift_t bands[] = {
             {' ', 0.43},   {'B', -0.8 },   {'c', -0.05 },
             {'C', 0.4 },   {'g', -0.35 },  {'G', 0.28 },
             {'H', 1.4 },   {'i', 0.32 },   {'I', 0.8 },
             {'J', 1.2 },   {'K', 1.7 },    {'L', 0.2 },
             {'o', 0.33 },  {'r', 0.14 },   {'R', 0.4 }, {'u', 2.5 },
             {'U', -1.3 },  {'v', 0},       {'W', 0.4 },
             {'w', -0.13 }, {'y', 0.32 },   {'Y', 0.7 },
             {'z', 0.26, },    {'\0', 0. } };

   if( strchr( "VNT", mag_band))
      return( 0.);
   for( i = 0; bands[i].band; i++)
      if( mag_band == bands[i].band)
         return( bands[i].v_mag_correction);
   if( !unknown_band_warning_shown)    /* show warning only once */
      {
      char buff[200];

      unknown_band_warning_shown = true;
      snprintf( buff, sizeof( buff),
                "Band '%c' is unknown to Find_Orb.  Photometry in this\n"
                "band will be unadjusted.  If this really is a legitimate\n"
                "photometric band,  please contact Project Pluto.\n",  mag_band);
      if( *get_environment_ptr( "BAND_WARNING"))
         generic_message_box( buff, "o");
      else
         debug_printf( "%s", buff);
      if( err_code)
         *err_code = -1;
      }
   return( 0.);
}

double override_abs_mag = 0.;

static double _calc_absolute_magnitude_internal( OBSERVE FAR *obs, int n_obs)
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
            int err = 0;

            if( object_type == OBJECT_TYPE_COMET
                            && obs->mag_band != default_comet_magnitude_type)
               use_obs = false;
            if( obs->obs_mag == BLANK_MAG || !obs->is_included)
               use_obs = false;
            obs->computed_mag = calc_obs_magnitude(
                  obs->solar_r, obs->r, earth_sun, NULL) - mag_band_shift( obs->mag_band, &err);
            if( use_obs)
               {
               rval += (obs->obs_mag - obs->computed_mag) / obs->mag_sigma;
               n_mags += 1. / obs->mag_sigma;
               }
            if( err)
               {
               char buff[90];

               recreate_observation_line( buff, obs, 0);
               debug_printf( "Unknown band '%s'\n", buff);
               }
            }
         }
      obs++;
      }
   if( n_mags)
      rval /= n_mags;
   if( override_abs_mag)
      rval = override_abs_mag;
   obs -= n_obs;
   for( obs_no = 0; obs_no < n_obs; obs_no++, obs++)
      if( rval)
         obs->computed_mag += rval;
      else
         obs->computed_mag = 0.;
   return( rval);
}

/* Computing an absolute magnitude using the above function may fail,
simply because none of the included observations has a magnitude given.
If so,  we try again with all observations temporarily turned on.

_That_ may fail because _none_ of the observations has a magnitude.
If so,  and if the DEFAULT_V_MAG is set,  we try again with all
observations temporarily set to that magnitude.  */

double calc_absolute_magnitude( OBSERVE FAR *obs, const int n_obs)
{
   double rval = _calc_absolute_magnitude_internal( obs, n_obs);

   if( !rval)        /* no mag computed;  try first w/all obs on */
      {
      const double default_v_mag = atof( get_environment_ptr( "DEFAULT_V_MAG"));
      OBSERVE *temp_obs = (OBSERVE *)calloc( n_obs, sizeof( OBSERVE));
      int obs_no;

      memcpy( temp_obs, obs, n_obs * sizeof( OBSERVE));
      for( obs_no = 0; obs_no < n_obs; obs_no++)
         temp_obs[obs_no].is_included = 1;
      rval = _calc_absolute_magnitude_internal( temp_obs, n_obs);
      if( !rval && object_type == OBJECT_TYPE_COMET)
         {
         obs_no = 0;
         while( obs_no < n_obs && obs[obs_no].mag_band != 'T'
                           && obs[obs_no].mag_band != 'N')
            obs_no++;
         if( obs_no == n_obs)    /* ah,  no obs have a 'comet' mag */
            {
            for( obs_no = 0; obs_no < n_obs; obs_no++)
               temp_obs[obs_no].mag_band = default_comet_magnitude_type;
            rval = calc_absolute_magnitude( temp_obs, n_obs);
            }
         }
      if( !rval && default_v_mag)   /* see 'environ.def' for an explanation of this */
         {
         for( obs_no = 0; obs_no < n_obs; obs_no++)
            {
            temp_obs[obs_no].obs_mag = default_v_mag;
            temp_obs[obs_no].mag_band = 'V';
            temp_obs[obs_no].mag_sigma = 5.;
            if( object_type == OBJECT_TYPE_COMET)
               temp_obs[obs_no].mag_band = default_comet_magnitude_type;
            }
         rval = calc_absolute_magnitude( temp_obs, n_obs);
         }
      if( rval)
         for( obs_no = 0; obs_no < n_obs; obs_no++)
            obs[obs_no].computed_mag = temp_obs[obs_no].computed_mag;
      free( temp_obs);
      }
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
         const double noise_in_sigmas)
{
   char buff[256];

   snprintf_err( buff, sizeof( buff), "%c,%d,%d,0,%f,%f",
               default_comet_magnitude_type,
               element_format, element_precision,
               max_residual_for_filtering, noise_in_sigmas);
   set_environment_ptr( "SETTINGS", buff);
   set_environment_ptr( "EPHEM_OPTIONS",
                   write_bit_string( buff, ephemeris_output_options, sizeof( buff)));
   snprintf_err( buff, sizeof( buff), "%.3f %f %u %d", probability_of_blunder * 100.,
              overobserving_time_span,
              overobserving_ceiling,
              use_blunder_method);
   set_environment_ptr( "FILTERING", buff);
   snprintf_err( buff, sizeof( buff), "%d %.2f %d %d %d", use_sigmas ? 1 : 0,
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

#if !defined( _WIN32) && !defined( __WATCOMC__)
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
      rval = findorb_already_running;
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

      strlcpy_error( tbuff, get_find_orb_text( i + 1000));
      tptr = (char *)find_nth_utf8_char( tbuff, 3);
      *tptr = '\0';        /* truncate month name at third char */
      assert( strlen( tbuff) < 16);
      strlcpy_error( month_names[i], tbuff);
      set_month_name( i + 1, month_names[i]);
      }
   return( 0);
}

int load_ephemeris_settings( ephem_option_t *ephemeris_output_options,
      int *n_steps, char *obscode, char *step_size, char *ephem_start,
      const char *config)
{
   char buff[100];
   const char *settings;
   int n_found, bytes_read;

   snprintf_err( buff, sizeof( buff), "STORED_EPHEM_%s", config);
   settings = get_environment_ptr( buff);
   if( !*settings && !strcmp( config, "F"))  /* factory default: geocentric, */
      settings = "500 20 10,16 1h +0";    /* 20 one-hour steps starting now  */
   n_found = sscanf( settings, "%s %d %s %s %n", obscode, n_steps, buff,
                                                step_size, &bytes_read);
   assert( n_found == 4 || !n_found);
   if( n_found < 4)
      return( -1);
   *ephemeris_output_options = parse_bit_string( buff);
   strlcpy_err( ephem_start, settings + bytes_read, 50);
   return( 0);
}

int save_ephemeris_settings( const ephem_option_t ephemeris_output_options,
      const int n_steps, const char *obscode, const char *step_size,
      const char *ephem_start, const char *config)
{
   char buff[150], eph_option_string[256], env_var[25];

   write_bit_string( eph_option_string, ephemeris_output_options,
                        sizeof( eph_option_string));
   snprintf_err( env_var, sizeof( env_var), "STORED_EPHEM_%s", config);
   snprintf_err( buff, sizeof( buff), "%s %d %s %s %s", obscode, n_steps,
                     eph_option_string, step_size, ephem_start);
   set_environment_ptr( env_var, buff);
   return( 0);
}

unsigned always_included_perturbers;

int get_defaults( ephem_option_t *ephemeris_output_options, int *element_format,
         int *element_precision, double *max_residual_for_filtering,
         double *noise_in_sigmas)
{
   ephem_option_t unused_ephemeris_output_options;
   int unused_element_format;
   int unused_element_precision;
   double unused_max_residual_for_filtering;
   double unused_noise_in_sigmas;
   const char *language = get_environment_ptr( "LANGUAGE");
   extern double minimum_jd, maximum_jd;
   extern double maximum_observation_span;
   extern int use_config_directory;
   extern double automatic_outlier_rejection_limit;
   extern double default_automatic_outlier_rejection_limit;
   extern unsigned max_n_sr_orbits;
   int use_sigmas_int;
   unsigned long obsolete_ephem_output_options;
   const char *override_fcct14_filename = get_environment_ptr( "FCCT14_FILE");
   const char *ephem_bitstring = get_environment_ptr( "EPHEM_OPTIONS");
   const char *eop_filename = get_environment_ptr( "EOP_FILE");
   const char *albedo = get_environment_ptr( "OPTICAL_ALBEDO");
   extern double minimum_observation_jd;
   extern double maximum_observation_jd;
   const char *outlier_limit = get_environment_ptr( "OUTLIER_REJECTION_LIMIT");
   const char *output_dir = get_environment_ptr( "OUTPUT_DIR");
   extern const char *output_directory;

   if( !output_directory && *output_dir)
      output_directory = output_dir;
#if !defined( _WIN32) && !defined( __WATCOMC__)
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
   if( !noise_in_sigmas)
      noise_in_sigmas = &unused_noise_in_sigmas;
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
               max_residual_for_filtering, noise_in_sigmas);

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
   _reset_sr_orbits( );
   sscanf( get_environment_ptr( "PERTURBERS"), "%x", &always_included_perturbers);
   if( *albedo)
      optical_albedo = atof( albedo);
   if( !maximum_observation_jd)     /* hasn't already been set elsewhere */
      {
      const char *obs_range = get_environment_ptr( "OBSERVATION_DATE_RANGE");
      int n_found;

      if( !*obs_range)
         obs_range = "1100,2300";
      n_found = sscanf( obs_range, "%lf,%lf",
                    &minimum_observation_jd,
                    &maximum_observation_jd);         /* see environ.def */
      assert( n_found == 2);
      minimum_observation_jd = YEAR_TO_JD( minimum_observation_jd);
      maximum_observation_jd = YEAR_TO_JD( maximum_observation_jd);
      }
   if( *outlier_limit)
      default_automatic_outlier_rejection_limit = atof( outlier_limit);
   automatic_outlier_rejection_limit = default_automatic_outlier_rejection_limit;
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
            if( counts[j] < counts[color] ||
                     (counts[j] == counts[color] && !(rand( ) % 3)))
               {
               if( counts[j] < counts[color])
                  changes_made = 1;
               color = j;
               sdata[i].color = (char)color;
               }
         assert( color >=0 && color < max_n_colors);
         }
      n_iterations++;
      assert( n_iterations < 50);
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
         strlcpy_error( rval[loc].code, obs[i].mpc_code);
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
