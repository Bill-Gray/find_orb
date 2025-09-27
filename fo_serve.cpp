/* fo_serve.cpp: Find_Orb for a server

Copyright (C) 2012, Project Pluto

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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>
#include "watdefs.h"
#include "stringex.h"
#include "sigma.h"
#include "afuncs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "date.h"
#include "monte0.h"
#include "cgi_func.h"

extern int debug_level;

int debug_level = 0;

uint64_t parse_bit_string( const char *istr);                /* miscell.cpp */
int fetch_astrometry_from_mpc( FILE *ofile, const char *desig);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);
char *get_file_name( char *filename, const char *template_file_name);
int sanity_test_observations( const char *filename);
int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
int get_defaults( ephem_option_t *ephemeris_output_options, int *element_format,
         int *element_precision, double *max_residual_for_filtering,
         double *noise_in_sigmas);                /* elem_out.cpp */
int make_pseudo_mpec( const char *mpec_filename, const char *obj_name);
                                               /* ephem0.cpp */
int inquire( const char *prompt, char *buff, const int max_len,
                     const int color);          /* fo_serve.cpp */
void refresh_console( void);                    /* fo_serve.cpp */
void move_add_nstr( const int col, const int row, const char *msg,
                     const int n_bytes);        /* fo_serve.cpp */

/* In the server flavor of Find_Orb,  warning messages such as "3
observations were made in daylight" or "couldn't find thus-and-such file"
are stored/added to the 'mpec_error_message' text.  This is then shown
near the top of the pseudo-MPEC (see 'ephem0.cpp'). */

int inquire( const char *prompt, char *buff, const int max_len,
                     const int color)
{
   extern char *mpec_error_message;

   INTENTIONALLY_UNUSED_PARAMETER( buff);
   INTENTIONALLY_UNUSED_PARAMETER( max_len);
   INTENTIONALLY_UNUSED_PARAMETER( color);
   if( !mpec_error_message)
      {
      const size_t len = strlen( prompt) + 1;

      mpec_error_message = (char *)malloc( len);
      strlcpy_err( mpec_error_message, prompt, len);
      }
   else
      {
      const size_t new_msg_len = strlen( mpec_error_message) + strlen( prompt) + 1;
      char *new_message = (char *)malloc( new_msg_len);

      strlcpy_err( new_message, mpec_error_message, new_msg_len);
      strlcat_err( new_message, prompt, new_msg_len);
      free( mpec_error_message);
      mpec_error_message = new_message;
      }
   return( 0);
}

/* In the (interactive) console Find_Orb,  these allow some functions in
orb_func.cpp to show info as orbits are being computed,  or to let you
abort processing by hitting a key.  In this non-interactive code,
they're mapped to do nothing. */

void refresh_console( void)
{
}

void move_add_nstr( const int col, const int row, const char *msg, const int n_bytes)
{
   INTENTIONALLY_UNUSED_PARAMETER( col);
   INTENTIONALLY_UNUSED_PARAMETER( row);
   INTENTIONALLY_UNUSED_PARAMETER( msg);
   INTENTIONALLY_UNUSED_PARAMETER( n_bytes);
}

int curses_kbhit_without_mouse( )
{
   return( 0);
}

static void show_problem_message( void)
{
   FILE *ifile = fopen( "problem.htm", "rb");

   if( !ifile)
      {
      printf( "<h2> Internal error:  no 'problem.htm' </h2>\n");
      printf( "<p> Please report this to Project Pluto;  if this happens,\n"
              "a fix <i>really</i> needs to be made!</p>\n");
      }
   else
      {
      char buff[200];

      while( fgets( buff, sizeof( buff), ifile))
         printf( "%s", buff);
      fclose( ifile);
      }
}

double current_jd( void);                       /* elem_out.cpp */
void compute_variant_orbit( double *variant, const double *ref_orbit,
                     const double n_sigmas);       /* orb_func.cpp */

int main( const int argc, const char **argv)
{
   const size_t max_buff_size = 400000;       /* room for 5000 obs */
   char *buff = (char *)malloc( max_buff_size);
   char mpec_name[100];
   char ephemeris_step_size[80], mpc_code[20];
   int n_ids, n_ephem_steps = 0, cgi_status;
   char field[30];
   OBJECT_INFO *ids;
   FILE *ifile;
   FILE *lock_file = fopen( "lock.txt", "w");
   extern const char *combine_all_observations;
   extern char *temp_obs_filename;     /* miscell.cpp */
   extern int forced_central_body;
   double jd_start = 0., jd_end = 0., user_selected_epoch = 0.;
   const double min_jd = 2378527.5;    /* 1800 feb 1 */
   const double max_jd = 2488083.5;    /* 2099 dec 1 */
   ephem_option_t ephemeris_output_options = OPTION_ROUND_TO_NEAREST_STEP;
   int residual_format = RESIDUAL_FORMAT_SHORT;
   double mag_limit = 99.;
   extern double ephemeris_mag_limit;
   size_t i, bytes_written = 0;
   int element_format = ELEM_OUT_ALTERNATIVE_FORMAT;
   int file_no = 0;
   const char *file_names[4] = { "mpec.htm", "combined.json", "elements.json", "ephemeri.json" };
   extern bool neocp_redaction_turned_on;
   int center_object = -2;
#ifndef _WIN32
   extern char **environ;
   extern bool findorb_already_running;

   avoid_runaway_process( 90);
#endif         /* _WIN32 */
   INTENTIONALLY_UNUSED_PARAMETER( argv);
   INTENTIONALLY_UNUSED_PARAMETER( argc);
   setvbuf( lock_file, NULL, _IONBF, 0);
   neocp_redaction_turned_on = false;
   fprintf( lock_file, "We're in\n");
   combine_all_observations = "";
   *ephemeris_step_size = '\0';
#ifndef _WIN32
   for( i = 0; environ[i]; i++)
      fprintf( lock_file, "%s\n", environ[i]);
#endif

   fprintf( lock_file, "setrlimit called\n");
   strlcpy_error( mpc_code, "500");

               /* get_defaults( ) collects a lot of data that's for the  */
               /* interactive find_orb program.  But it also sets some   */
               /* important internal values for blunder detection,  etc. */
               /* So we still call it:                                   */
   get_defaults( NULL, NULL, NULL, NULL, NULL);
#ifndef _WIN32
   if( findorb_already_running)
      {
      printf( "Content-type: text/html\n\n");
      printf( "<h1> Server is busy.  Try again in a minute or two. </h1>");
      printf( "<h3> Your orbit is very important to us! </h3>");
      printf( "<p> I've had some problems with this server suddenly getting a few hundred\n");
      printf( "sets of astrometry thrown at it at once.  It doesn't have the computational\n");
      printf( "power to keep up.  It's now throttled to handling one orbit at a time. </p>\n");
      printf( "<p> The above 'try again in a minute' is about right:  try again,  and\n");
      printf( "it'll probably work. </p>\n");
      return( 0);
      }
#endif
   fprintf( lock_file, "defaults loaded\n");
   cgi_status = initialize_cgi_reading( );
   if( cgi_status <= 0)
      {
      printf( "Content-type: text/html\n\n");
      printf( "<p> <b> CGI data reading failed : error %d </b>", cgi_status);
      printf( "This isn't supposed to happen.</p>\n");
      return( 0);
      }
   while( !get_cgi_data( field, buff, NULL, max_buff_size))
      {
      if( !strcmp( field, "TextArea") || !strcmp( field, "upfile"))
         {
         if( strlen( buff) > 70)
            {
            char *tptr = strstr( buff, "Debug level:");
            FILE *ofile = fopen( temp_obs_filename,
                               (bytes_written ? "ab" : "wb"));

            assert( ofile);
            for( i = 1; buff[i + 1]; i++)     /* cvt Mac-style CR endings */
               if( buff[i] == 13 && buff[i + 1] != 10)   /* to LF endings */
                  buff[i] = 10;
            bytes_written += fwrite( buff, 1, strlen( buff), ofile);
            fprintf( ofile, "\n");
            fclose( ofile);
            if( tptr)
               debug_level = atoi( tptr + 12);
            }
         }
      else if( !strcmp( field, "obj_name") && *buff)
         {
         FILE *ofile = fopen( temp_obs_filename,
                               (bytes_written ? "ab" : "wb"));

         assert( ofile);
         bytes_written += fetch_astrometry_from_mpc( ofile, buff);
         fclose( ofile);
         }
      else if( !strcmp( field, "min_obs_t") && *buff)
         {
         extern double minimum_observation_jd;  /* default is 1     */

         minimum_observation_jd = get_time_from_string( 0., buff,
                          FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, NULL);
         }
      else if( !strcmp( field, "max_obs_t") && *buff)
         {
         extern double maximum_observation_jd;  /* default is +1e+9. */

         maximum_observation_jd = get_time_from_string( 0., buff,
                          FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, NULL);
         }
      else if( !strcmp( field, "year") || !strcmp( field, "t_end"))
         {
         const double jd = get_time_from_string( current_jd( ), buff,
             CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD | FULL_CTIME_TWO_DIGIT_YEAR, NULL);

         if( jd < min_jd || jd > max_jd)
            {
            printf( "Content-type: text/html\n\n");
            printf( "<b>Ephemeris date out of range</b>\n");
            printf( "<p>'%s' parsed as JD %f\n", buff, jd);
            printf( "The ephemeris starting date must be between JD %.1f and %.1f.</p>\n",
                              min_jd, max_jd);
            return( 0);
            }
         if( *field == 'y')
            jd_start = jd;
         else
            jd_end = jd;
         }
      else if( !strcmp( field, "epoch"))
         {
         if( strcmp( buff, "default"))
            user_selected_epoch = get_time_from_string( 0., buff,
                CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD | FULL_CTIME_TWO_DIGIT_YEAR, NULL);
         if( user_selected_epoch
                  && (user_selected_epoch < min_jd || user_selected_epoch > max_jd))
            {
            printf( "Content-type: text/html\n\n");
            printf( "Epoch of elements out of range\n");
            printf( "'%s' parsed as JD %f\n", buff, user_selected_epoch);
            printf( "The epoch must be between JD %.1f and %.1f.\n",
                              min_jd, max_jd);
            return( 0);
            }
         }
      else if( !strcmp( field, "n_steps"))
         n_ephem_steps = atoi( buff);
      else if( !strcmp( field, "stepsize"))
         strlcpy_err( ephemeris_step_size, buff, sizeof( ephemeris_step_size));
      else if( !strcmp( field, "mpc_code"))
         strlcpy_err( mpc_code, buff, sizeof( mpc_code));
      else if( !strcmp( field, "faint_limit"))
         mag_limit = atof( buff);
      else if( !strcmp( field, "element_center"))
         center_object = atoi( buff);
      else if( !strcmp( field, "resids"))
         {
         if( buff[0] == '1')
            residual_format |= RESIDUAL_FORMAT_PRECISE;
         if( buff[0] == '2')
            residual_format |= RESIDUAL_FORMAT_TIME_RESIDS;
         }
      else if( !strcmp( field, "constell"))
         ephemeris_output_options |= OPTION_CONSTELLATION;
      else if( !strcmp( field, "total_motion"))
         ephemeris_output_options |= OPTION_MOTION_OUTPUT;
      else if( !strcmp( field, "separate_motions"))
         ephemeris_output_options |= OPTION_SEPARATE_MOTIONS;
      else if( !strcmp( field, "alt_az"))
         ephemeris_output_options |= OPTION_ALT_AZ_OUTPUT;
      else if( !strcmp( field, "radial"))
         ephemeris_output_options |= OPTION_RADIAL_VEL_OUTPUT;
      else if( !strcmp( field, "phase"))
         ephemeris_output_options |= OPTION_PHASE_ANGLE_OUTPUT;
      else if( !strcmp( field, "pab"))
         ephemeris_output_options |= OPTION_PHASE_ANGLE_BISECTOR;
      else if( !strcmp( field, "hel_ec"))
         ephemeris_output_options |= OPTION_HELIO_ECLIPTIC;
      else if( !strcmp( field, "ground"))
         ephemeris_output_options |= OPTION_GROUND_TRACK;
      else if( !strcmp( field, "visib"))
         ephemeris_output_options |= OPTION_VISIBILITY;
      else if( !strcmp( field, "top_ec"))
         ephemeris_output_options |= OPTION_TOPO_ECLIPTIC;
      else if( !strcmp( field, "unobs"))
         ephemeris_output_options |= OPTION_SUPPRESS_UNOBSERVABLE;
      else if( !strcmp( field, "sigmas"))
         ephemeris_output_options |= OPTION_SHOW_SIGMAS;
      else if( !strcmp( field, "comp_fr"))
         ephemeris_output_options |= OPTION_COMPUTER_FRIENDLY;
      else if( !strcmp( field, "lun_elong"))
         ephemeris_output_options |= OPTION_LUNAR_ELONGATION;
      else if( !strcmp( field, "lun_alt"))
         ephemeris_output_options |= OPTION_MOON_ALT;
      else if( !strcmp( field, "lun_az"))
         ephemeris_output_options |= OPTION_MOON_AZ;
      else if( !strcmp( field, "sun_alt"))
         ephemeris_output_options |= OPTION_SUN_ALT;
      else if( !strcmp( field, "sun_az"))
         ephemeris_output_options |= OPTION_SUN_AZ;
      else if( !strcmp( field, "sky_br"))
         ephemeris_output_options |= OPTION_SKY_BRIGHTNESS;
      else if( !strcmp( field, "redact_neocp"))
         neocp_redaction_turned_on = true;
      else if( !strcmp( field, "ephem_type"))
         ephemeris_output_options += atoi( buff);
      else if( !strcmp( field, "file_no"))
         {
         file_no = atoi( buff);
         if( file_no < 0 || file_no > 3)
            file_no = 0;
         }
      else if( !strcmp( field, "ephem_opts"))
         ephemeris_output_options = parse_bit_string( buff);
      else if( *field == 'e' && isdigit( field[1]))
         ephemeris_output_options |= ((ephem_option_t)1 << atoi( field + 1));
      else if( !strcmp( field, "language"))
         {
         extern char findorb_language;

         findorb_language = buff[0];
         }
      }
   fprintf( lock_file, "Options read and parsed\n");
   if( !bytes_written)
      {
      printf( "Content-type: text/html\n\n");
      show_problem_message( );
      return( 0);
      }
   ephemeris_mag_limit = mag_limit;
   forced_central_body = center_object;

   load_up_sigma_records( "sigma.txt");
   fprintf( lock_file, "Default uncertainties table read\n");

   ids = find_objects_in_file( temp_obs_filename, &n_ids, NULL);
   fprintf( lock_file, "%d objects found in file\n", n_ids);
   if( n_ids <= 0)
      {
      printf( "Content-type: text/html\n\n");
      printf( n_ids == -1 ? "Couldn't open observation file\n" :
                            "No valid observations found in that data\n");
      return( 0);
      }

   const char *orbit_constraints = "";
   OBSERVE FAR *obs;
   const int n_obs = ids[0].n_obs;
   extern int n_obs_actually_loaded;
   const int element_precision = 5;
   double epoch_shown, curr_epoch, orbit[2 * MAX_N_PARAMS];
   double *orbits_to_use = orbit;
   extern const char *ephemeris_filename;
   extern const char *residual_filename;
   extern int available_sigmas;
   unsigned n_orbits_in_ephem = 1;

   ifile = fopen( temp_obs_filename, "rb");
   fprintf( lock_file, "'%s' %sopened\n", temp_obs_filename, ifile ? "" : "not ");
   obs = load_object( ifile, ids, &curr_epoch, &epoch_shown, orbit);
   fclose( ifile);
   fprintf( lock_file, "%d obs %sloaded\n", n_obs, (obs ? "" : "not "));
   if( available_sigmas == COVARIANCE_AVAILABLE)
      {
      extern unsigned perturbers;

      perturbers = 0x7fe;     /* all planets plus the moon */
      if( user_selected_epoch)
         epoch_shown = user_selected_epoch;
      full_improvement( obs, n_obs, orbit, curr_epoch, "",
                    ORBIT_SIGMAS_REQUESTED, epoch_shown);
      }
   fprintf( lock_file, "avail sigmas %d\n", available_sigmas);

   if( n_obs_actually_loaded > 1 && curr_epoch > 0.)
      {
      write_out_elements_to_file( orbit, curr_epoch, epoch_shown,
            obs, n_obs, orbit_constraints, element_precision,
            0, element_format);
      }
   else
      {
      show_problem_message( );
      return( 0);
      }

   create_obs_file( obs, n_obs_actually_loaded, 0, 0);
   if( available_sigmas == COVARIANCE_AVAILABLE)
      {
      extern int n_orbit_params;

      n_orbits_in_ephem = 2;
      compute_variant_orbit( orbit + n_orbit_params, orbit, 1.);
      }
   if( available_sigmas == SR_SIGMAS_AVAILABLE)
      {
      extern double *sr_orbits;
      extern unsigned n_sr_orbits;

      orbits_to_use = sr_orbits;
      n_orbits_in_ephem = n_sr_orbits;
      }
   if( jd_end > min_jd)
      {
      if( !n_ephem_steps)
         n_ephem_steps = 1;
      if( !*ephemeris_step_size)
         snprintf_err( ephemeris_step_size, sizeof( ephemeris_step_size),
                  "%.6fd", (jd_end - jd_start) / (double)n_ephem_steps);
      }
   fprintf( lock_file, "ready for ephem\n");
   if( ephemeris_in_a_file_from_mpc_code( ephemeris_filename,
               orbits_to_use, obs, n_obs_actually_loaded,
               curr_epoch, jd_start, ephemeris_step_size,
               n_ephem_steps, mpc_code,
               ephemeris_output_options,
               n_orbits_in_ephem))
      {
      printf( "Ephem generation failed\n");
      return( 0);
      }
   fprintf( lock_file, "Ephem generated\n");
   write_residuals_to_file( residual_filename, temp_obs_filename,
                  n_obs_actually_loaded, obs, residual_format);
   fprintf( lock_file, "Resids written\n");

   strlcpy_error( mpec_name, get_environment_ptr( "MPEC_NAME"));
   if( !*mpec_name)
      strlcpy_error( mpec_name, "mpec.htm");
   make_pseudo_mpec( mpec_name, ids[0].obj_name);
   fprintf( lock_file, "pseudo-MPEC made\n");
   unload_observations( obs, n_obs);
   fprintf( lock_file, "Obs unloaded\n");
   if( !file_no)
      {
      printf( "Content-type: text/html\n\n");
      ifile = fopen( mpec_name, "rb");
      }
   else
      {
      printf( "Content-type: application/json\n\n");
      ifile = fopen( file_names[file_no], "rb");
      }
   while( fgets( buff, max_buff_size, ifile))
      printf( "%s", buff);
   fflush( stdout);
   fclose( ifile);
   if( (i = strlen( mpec_name)) > 8)
      {
      const int counter = atoi( mpec_name + i - 7);

      if( counter)
         {
         snprintf_err( mpec_name + i - 7, 8, "%03d.htm", counter % 999 + 1);
         set_environment_ptr( "MPEC_NAME", mpec_name);
         }
      }
   fprintf( lock_file, "Cleaning memory\n");
   clean_up_find_orb_memory( );
   fprintf( lock_file, "Done!\n");
   return( 0);
}
