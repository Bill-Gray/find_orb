/* fo.cpp: main driver for non-interactive Find_Orb

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

/* Under *nix,  this program can distribute the work of determining
orbits among multiple processes.  If there are N objects for which
orbits need to be determined and Np processes (Np can be specified
on the command line),  each process gets N/Np objects to determine.
At the end,  the "original" process merges the results.  But at least
at present,  this only works on *nix systems... FORKING is undefined
for Windows and other non-*nix systems. */

#if defined( __linux) || defined( __unix__) || defined( __APPLE__)
#define FORKING
#endif

#ifdef FORKING
   #include <unistd.h>     /* Symbolic Constants */
   #include <sys/types.h>  /* Primitive System Data Types */
   #include <errno.h>      /* Errors */
   #include <stdio.h>      /* Input/Output */
   #include <sys/wait.h>   /* Wait for Process Termination  */
            /* above basically allows for forking so we can */
            /* run different objects on different cores     */
   #include <sys/time.h>         /* these allow resource limiting */
   #include <sys/resource.h>     /* see '-r' command switch below */
#endif
#ifdef __WATCOMC__
   #include <io.h>
#endif
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <time.h>
#include "watdefs.h"
#include "sigma.h"
#include "afuncs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "date.h"
#include "monte0.h"

extern int debug_level;

#ifndef strlcpy
size_t strlcpy(char *dest, const char *src, size_t size);   /* miscell.cpp */
#endif

   /* MSVC/C++ lacks snprintf.  See 'ephem0.cpp' for details. */
#if defined(_MSC_VER) && _MSC_VER < 1900
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

int debug_level = 0;
extern const char *sof_filename, *sofv_filename;

char *get_file_name( char *filename, const char *template_file_name);
int sanity_test_observations( const char *filename);
int debug_printf( const char *format, ...)                 /* runge.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
int get_defaults( ephem_option_t *ephemeris_output_options, int *element_format,
         int *element_precision, double *max_residual_for_filtering,
         double *noise_in_arcseconds);                /* elem_out.cpp */
int inquire( const char *prompt, char *buff, const int max_len,
                     const int color);                /* fo.cpp */
void refresh_console( void);                          /* fo.cpp */
void move_add_nstr( const int col, const int row, const char *msg,
                     const int n_bytes);              /* fo.cpp */
double current_jd( void);                       /* elem_out.cpp */
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
int make_pseudo_mpec( const char *mpec_filename, const char *obj_name);
                                               /* ephem0.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
int reset_astrometry_filename( int *argc, const char **argv);
uint64_t parse_bit_string( const char *istr);                /* miscell.cpp */
FILE *open_json_file( char *filename, const char *env_ptr, const char *default_name,
                  const char *packed_desig, const char *permits); /* ephem0.cpp */

/* In this non-interactive version of Find_Orb,  we just print out warning
messages such as "3 observations were made in daylight" or "couldn't find
thus-and-such file".  These will also be logged in 'debug.txt'.  We then
proceed as if nothing had happened: */

int inquire( const char *prompt, char *buff, const int max_len,
                     const int color)
{
   INTENTIONALLY_UNUSED_PARAMETER( buff);
   INTENTIONALLY_UNUSED_PARAMETER( max_len);
   INTENTIONALLY_UNUSED_PARAMETER( color);
   printf( "\n%s\n", prompt);
   return( 0);
}

static void object_comment_text( char *buff, const OBJECT_INFO *id)
{
   sprintf( buff, "%d observations; ", id->n_obs);
   make_date_range_text( buff + strlen( buff), id->jd_start, id->jd_end);
}

/* In the (interactive) console Find_Orb,  these allow some functions in
orb_func.cpp to show info as orbits are being computed.  In this
non-interactive code,  they're mapped to do nothing. */

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

static double curr_jd( void)
{
   const double jd_1970 = 2440587.5;

   return( jd_1970 + (double)time( NULL) / seconds_per_day);
}

static int remove_single_observation_objects( OBJECT_INFO *ids, const int n_ids)
{
   int i, j;

   for( i = j = 0; i < n_ids; i++)
      if( ids[i].n_obs > 1)
         ids[j++] = ids[i];
   return( j);
}

void make_config_dir_name( char *oname, const char *iname);

#ifdef FORKING
static int unlink_config_file( const char *filename)
{
   char buff[255];
   int err_code;
   extern int use_config_directory;          /* miscell.c */

   get_file_name( buff, filename);
   if( use_config_directory)
      {
      char cpath[255];

      make_config_dir_name( cpath, buff);
#ifdef _WIN32                /* MS is different. */
      err_code = _unlink( cpath);
#else
      err_code = unlink( cpath);
#endif
      }
   else
#ifdef _WIN32
      err_code = _unlink( buff);
#else
      err_code = unlink( buff);
#endif
   return( err_code);
}

static bool unlink_partial_files = true;

static void combine_element_files( const char *filename, const int n_processes)
{
   FILE **input_files = (FILE **)calloc( (size_t)n_processes, sizeof( FILE *));
   char buff[400];
   int i, quit;
   FILE *ofile;
   extern int process_count;

   if( !strcmp( filename, sof_filename) || !strcmp( filename, sofv_filename))
      {
      FILE *ifile = fopen_ext( filename, "fclr");

      assert( ifile);
      if( !fgets( buff, sizeof( buff), ifile))
         {        /* get header line */
         printf( "Failed to read sof.txt\n");
         exit( -1);
         }
      fclose( ifile);
      ofile = fopen( filename, "w");
      fwrite( buff, strlen( buff), 1, ofile);
      }           /* above writes header line back out */
   else
      ofile = fopen( filename, "w");
   assert( ofile);
   for( i = 0; i < n_processes; i++)
      {
      process_count = i + 1;
      input_files[i] = fopen_ext( get_file_name( buff, filename), "tfclr");
      }
   for( i = quit = 0; !quit; i = (i + 1) % n_processes)
      if( !fgets( buff, sizeof( buff), input_files[i]))
         quit = 1;
      else
         {
         extern const char *elements_filename;

         fputs( buff, ofile);
         if( !strcmp( filename, elements_filename))
            {
            while( fgets( buff, sizeof( buff), input_files[i]) &&
                     memcmp( buff, "# Sigmas avail", 14))
               fputs( buff, ofile);
            fputs( buff, ofile);
            }
         }
   for( i = 0; i < n_processes; i++)
      {
      int err_code;

      fclose( input_files[i]);
      process_count = i + 1;
      if( unlink_partial_files)
         {
         err_code = unlink_config_file( filename);
         if( err_code)
            perror( buff);
         assert( !err_code);
         }
      }
   free( input_files);
   fclose( ofile);
}
#endif

void compute_variant_orbit( double *variant, const double *ref_orbit,
                     const double n_sigmas);       /* orb_func.cpp */

static void extract_value( char *obuff, const char *ibuff,
                        const unsigned n_digits)
{
   unsigned i;

   while( *ibuff && *ibuff != '.')
      *obuff++ = *ibuff++;
   i = 0;
   while( i++ < n_digits + 1 && *ibuff > ' ')
      *obuff++ = *ibuff++;
   while( isdigit( *ibuff))
      ibuff++;
   if( !memcmp( ibuff, " +/- ", 5))
      {
      memcpy( obuff, ibuff, 5);
      ibuff += 5;
      obuff += 5;
      while( *ibuff > ' ')
         *obuff++ = *ibuff++;
      }
}

static size_t summ_sort_column = 0;

int summ_compare( const void *a, const void *b)
{
   const char *a1 = *(const char **)a;
   const char *b1 = *(const char **)b;

   a1 = strstr( a1, ".htm");
   b1 = strstr( b1, ".htm");
   assert( a1);
   assert( b1);
   if( strlen( a1) < summ_sort_column)
      return( -1);
   if( strlen( b1) < summ_sort_column)
      return( 1);
   a1 += summ_sort_column;
   b1 += summ_sort_column;
   if( summ_sort_column == 123)     /* sorting ephemeris unc: */
      {                             /* comparison is more complex */
      if( a1[4] == 'd' && b1[4] != 'd')
         return(  1);
      if( a1[4] != 'd' && b1[4] == 'd')
         return( -1);
      if( a1[4] == '\'' && b1[4] != '\'')
         return(  1);
      if( a1[4] != '\'' && b1[4] == '\'')
         return( -1);
      return( atof( a1) < atof( b1) ? -1 : 1);
      }
   return( strcmp( a1, b1));
}

static void get_summary_info( char *buff, const char *mpec_filename)
{
   FILE *ifile = fopen_ext( mpec_filename, "frb");
   char ibuff[400], *tptr;
   unsigned i;

   memset( buff, ' ', 80);
   buff[80] = '\0';
   while( fgets( ibuff, sizeof( ibuff), ifile))
      {
      if( (tptr = strstr( ibuff, "Pseudo-MPEC for")) != NULL)
         {
         tptr += 16;
         for( i = 0; *tptr && *tptr != '<' && i < 30; i++)
            buff[i] = *tptr++;
         }
      else if( ibuff[0] == 'a' && ibuff[1] == ' ')
         extract_value( buff + 17, ibuff + 2, 3);
      else if( ibuff[0] == 'e' && ibuff[1] == ' ')
         extract_value( buff + 34, ibuff + 2, 3);
      if( (tptr = strstr( ibuff, "Incl.")) != NULL)
         extract_value( buff + 51, tptr + 5, 1);
      else if( (tptr = strstr( ibuff, " Ea ")) != NULL)
         memcpy( buff + 68, tptr + 4, 7);
      }
   fclose( ifile);
}

#define VT100_RED       '1'
#define VT100_GREEN     '2'
#define VT100_YELLOW    '3'
#define VT100_BLUE      '4'
#define VT100_PURPLE    '5'
#define VT100_CYAN      '6'
#define VT100_GRAY      '7'

/* This inserts VT100 color codes so that text such as,  say,

...e=1.005...

can be transformed into the following,  which will cause the eccentricity
to be highlighted in red :

...\033[33me=1.005\033[0m...           */

static void add_vt100_colors( char *text, size_t nbytes, const char color)
{
   size_t len;

   memmove( text + 5, text, strlen( text) + 1);
   memcpy( text, "\033[3xm", 5);
   text[3] = color;
   text += 5;
   len = strlen( text);
   text += (len < nbytes ? len : nbytes);
   memmove( text + 4, text, strlen( text) + 1);
   memcpy( text, "\033[0m", 4);
}

static void colorize_text( char *text)
{
   char *tptr = strstr( text, "a=");

   if( tptr && atof( tptr + 2) > 7.)
      add_vt100_colors( tptr, 9, VT100_GREEN);
   tptr = strstr( text, "e=");
   if( tptr && atof( tptr + 2) > .9)
      add_vt100_colors( tptr, 8, VT100_YELLOW);
   tptr = strstr( text, "i=");
   if( tptr && atof( tptr + 2) > 60.)
      add_vt100_colors( tptr, 5, VT100_CYAN);
   tptr = strstr( text, "MOID ");
   if( tptr && atof( tptr + 5) < 0.0105)
      add_vt100_colors( tptr, 10, VT100_RED);
}

static int create_combined_json_header( const OBJECT_INFO *ids,
         const unsigned n_objs, const char *ofilename)
{
   char buff[200];
   FILE *ofile = fopen_ext( get_file_name( buff, ofilename), "tfcwb");
   unsigned i;

   fprintf( ofile, "{\n  \"num\": %u,\n", n_objs);
   fprintf( ofile, "  \"ids\":\n  [\n");
   for( i = 0; i < n_objs; i++)
      fprintf( ofile, "    \"%s\"%c\n", ids[i].obj_name,
               (i == n_objs - 1 ? ' ' : ','));
   fprintf( ofile, "  ],\n");
   fprintf( ofile, "  \"objects\":\n  {\n");
   fclose( ofile);
   return( 0);
}

static int add_json_data( const char *ofilename, const bool have_json_ephem,
            const char *packed_desig, const bool is_last_call)
{
   char buff[200];
   FILE *ofile = fopen_ext( get_file_name( buff, ofilename), "tfcab");
   FILE *ifile;
   bool found_start = false, found_end = false;

   if( have_json_ephem)
      ifile = open_json_file( buff, "JSON_COMBINED_NAME", "combined.json", packed_desig, "rb");
   else
      ifile = open_json_file( buff, "JSON_ELEMENTS_NAME", "elements.json", packed_desig, "rb");
   while( !found_start && fgets_trimmed( buff, sizeof( buff), ifile))
      if( !strcmp( buff, "  {"))
         found_start = true;
   assert( found_start);

   while( !found_end && fgets_trimmed( buff, sizeof( buff), ifile))
      {
      if( !strcmp( buff, "    }") && !is_last_call)
         {
         strcat( buff, ",");
         found_end = true;
         }
      fprintf( ofile, "%s\n", buff);
      }
   fclose( ifile);
   fclose( ofile);
   return( 0);
}


/* I really should use getopt() or a portable variant.  However,  this has
been sufficiently effective thus far... */

static const char *get_arg( const int argc, const char **argv, const int idx)
{
   if( argv[idx][2] || idx == argc - 1)
      return( argv[idx] + 2);
   else
      return( argv[idx + 1]);
}

int main( int argc, const char **argv)
{
   char tbuff[300], mpc_code[20];
   char **summary_lines = NULL;
   const char *separate_residual_file_name = NULL;
   const char *mpec_path = NULL;
   int n_ids, i, starting_object = 0;
   int n_processes = 1;
   OBJECT_INFO *ids;
   int total_objects = 0;
   FILE *ifile;
   extern int process_count;
   int n_lines_written = 0;
   FILE *summary_ofile = NULL;
   extern int forced_central_body;
   extern int use_config_directory;          /* miscell.c */
   int element_precision = 5;
   bool all_heliocentric = true;
   bool use_colors = true;
   bool show_processing_steps = true;
   ephem_option_t ephemeris_output_options
               = OPTION_SHOW_SIGMAS | OPTION_ROUND_TO_NEAREST_STEP;
   time_t update_time, t0;
   double ephem_end_jd = 0.;
   extern bool is_default_ephem;
   bool drop_single_obs = true;
   const char *ephem_option_string = NULL;
#ifdef FORKING
   int child_status;
#endif

   if( !strcmp( argv[0], "fo"))
      use_config_directory = true;
   else
      use_config_directory = false;
   *mpc_code = '\0';
   if( reset_astrometry_filename( &argc, argv))
      drop_single_obs = false;

   for( i = 1; i < argc; i++)       /* check to see if we're debugging: */
      if( argv[i][0] == '-')
         {
         const char *arg = get_arg( argc, argv, i);

         switch( argv[i][1])
            {
            case 'a':
               {
               extern int separate_periodic_comet_apparitions;

               separate_periodic_comet_apparitions ^= 1;
               }
               break;
            case 'b':
               separate_residual_file_name = arg;
               break;
            case 'c':
               {
               extern int combine_all_observations;

               combine_all_observations = 1;
               }
               break;
            case 'C':
               strlcpy( mpc_code, arg, sizeof( mpc_code));
               break;
            case 'd':
               debug_level = atoi( arg);
               if( !debug_level)
                  debug_level = 1;
               debug_printf( "fo: debug_level = %d; %s %s\n",
                           debug_level, __DATE__, __TIME__);
               break;
            case 'e':
               {
               extern const char *ephemeris_filename;

               ephemeris_filename = arg;
               is_default_ephem = false;
               }
               break;
            case 'E':
               ephem_option_string = arg;
               break;
            case 'f':                     /* obj desig specified;  fall through */
               break;
            case 'h':                     /* show planet-centric orbits */
               all_heliocentric = false;
               break;
#ifdef FORKING
            case 'k':
               unlink_partial_files = false;
               break;
#endif
            case 'i':
               {
               extern int ignore_prev_solns;

               ignore_prev_solns = 1;
               }
               break;
            case 'm':
               mpec_path = arg;
               break;
            case 'M':
               {
               extern const char *mpcorb_dot_sof_filename;

               mpcorb_dot_sof_filename = arg;
               }
               break;
            case 'n':
               starting_object = atoi( arg);
               break;
            case 'O':          /* write output files to specified dir */
               {
               extern const char *output_directory;

               output_directory = arg;
               }
               break;
            case 'o':            /* obj designation / ephemeris from orbital */
               break;            /* elems:  fall through, handle below */
            case 'p':
               {
               FILE *ifile = fopen_ext( "dummy.txt", "tfcw");

               fclose( ifile);
               n_processes = atoi( arg);
               }
               break;
            case 'q':            /* "quiet" */
               show_processing_steps = false;
               break;
#ifdef FORKING
            case 'r':
               {                    /* set 'soft' & 'hard' limits for CPU */
               struct rlimit r;     /* run time,  in seconds,  to avoid   */
               int soft_limit, hard_limit;      /* runaway processes */

               if( sscanf( arg, "%d,%d", &soft_limit, &hard_limit) == 2)
                  {
                  r.rlim_cur = (rlim_t)soft_limit;
                  r.rlim_max = (rlim_t)hard_limit;
                  setrlimit( RLIMIT_CPU, &r);
                  }
               }
               break;
#endif
            case 's':
               sanity_test_observations( argv[1]);
               printf( "Sanity check complete\n");
               return( 0);
            case 'S':
               {
               char curr_time[50];

               full_ctime( curr_time, current_jd( ), FULL_CTIME_YMD);
               summary_ofile = fopen( arg, "wb");
               assert( summary_ofile);
               ifile = fopen_ext( "summ.htm", "fcrb");
               assert( ifile);
               while( fgets( tbuff, sizeof( tbuff), ifile))
                  {
                  text_search_and_replace( tbuff, "%TIME%", curr_time);
                  fputs( tbuff, summary_ofile);
                  }
               fclose( ifile);
               }
               break;
            case 't':
               if( argv[i][2] == 'e' || argv[i][2] == 'E')
                  {
                  double jd = get_time_from_string( curr_jd( ), argv[i] + 3,
                           CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD
                           | FULL_CTIME_TWO_DIGIT_YEAR, NULL);

                  if( argv[i][2] == 'e')
                     ephem_end_jd = jd;
                  else
                     {
                     extern double override_epoch_shown;

                     override_epoch_shown = jd;
                     }
                  }
               else
                  total_objects = atoi( arg);
               break;
            case 'v':
               use_colors = false;
               break;
            case 'x':
               {
               extern const char *alt_config_directory;

               use_config_directory = true;
               alt_config_directory = arg;
               }
               break;
            case 'z':
               use_config_directory = true;
               break;
            default:
               printf( "Unknown command-line option '%s'\n", argv[i]);
               return( -1);
            }
         }
   for( i = 1; i < argc; i++)
      {
      const char *tptr = strchr( argv[i], '=');

      if( tptr && argv[i][0] != '-')
         {
         const size_t len = tptr - argv[i];

         memcpy( tbuff, argv[i], len);
         tbuff[len] = '\0';
         set_environment_ptr( tbuff, argv[i] + len + 1);
         }
      }

               /* get_defaults( ) collects a lot of data that's for the  */
               /* interactive find_orb program.  But it also sets some   */
               /* important internal values for blunder detection,  etc. */
               /* So we still call it:                                   */
   get_defaults( &ephemeris_output_options,
                         NULL, &element_precision, NULL, NULL);
   if( all_heliocentric)
      forced_central_body = 0;
   if( ephem_option_string)
      ephemeris_output_options = parse_bit_string( ephem_option_string);

   load_up_sigma_records( "sigma.txt");
   if( debug_level)
      debug_printf( "%d sigma recs read\n", i);

   if( argc < 2)
      {
      printf( "'fo' needs the name of an input file of MPC-formatted\n");
      printf( "astrometry as a command-line argument.\n");
      return( -2);
      }

   ids = find_objects_in_file( argv[1], &n_ids, NULL);
   if( n_ids <= 0)
      {        /* no objects found,  or file not found */
      const char *err_msg;

      if( n_ids == -1)
         err_msg = "Couldn't locate the file";
      else
         err_msg = "No objects found in file";
      printf( "%s '%s'\n", err_msg, argv[1]);
      return( -1);
      }

   if( drop_single_obs)
      n_ids = remove_single_observation_objects( ids, n_ids);
   if( show_processing_steps)
      printf( "Processing %d objects\n", n_ids);
   if( !total_objects)
      total_objects = n_ids;

   create_combined_json_header( ids, total_objects, "total.json");
   t0 = update_time = time( NULL);
#ifdef FORKING
   while( process_count < n_processes - 1)
      {
      const pid_t childpid = fork( );

      if( childpid == -1)      /* fork( ) returns -1 on failure */
         {
         perror("fork"); /* display error message */
         exit(0);
         }
      else if( childpid == 0)     /* we're a child process */
         {
//       printf( "Hi!  I'm child %d.  My PID is %d; parent's is %d\n",
//                process_count, getpid( ), getppid( ));
         }
      else
         break;       /* break out of loop,  signalling we're a parent */
      process_count++;
      }
   if( n_processes > 1)
      process_count++;
   if( show_processing_steps)
      printf( "Process count %d\n", process_count);
#endif

   if( summary_ofile)
      summary_lines = (char **)calloc( n_ids - starting_object + 1,
                                                    sizeof( char *));
   ifile = fopen( argv[1], "rb");
   if( total_objects > n_ids - starting_object)
      total_objects = n_ids - starting_object;
   for( i = starting_object; i < starting_object + total_objects; i++)
      if( n_processes == 1 || i % n_processes == process_count - 1)
         {
         const char *orbit_constraints = "";
         OBSERVE FAR *obs;
         const int n_obs = ids[i].n_obs;

         if( n_processes == 1 && show_processing_steps)
            printf( "%d: %s", i + 1, ids[i].obj_name);
         if( n_obs < 2 && drop_single_obs)
            printf( "; skipping\n");
         else
            {
            extern int append_elements_to_element_file;
            extern int n_obs_actually_loaded;
            extern char orbit_summary_text[];
            long file_offset = ids[i].file_offset - 40L;
            int element_options = ELEM_OUT_ALTERNATIVE_FORMAT;
            double epoch_shown, curr_epoch, orbit[12];
            bool have_json_ephem = false;

            if( all_heliocentric)
               element_options |= ELEM_OUT_HELIOCENTRIC_ONLY;
                /* Start a bit ahead of the actual data,  just in case */
                /* there's a #Sigma: or similar command in there: */
            if( file_offset < 0L)
               file_offset = 0L;
            fseek( ifile, file_offset, SEEK_SET);
            obs = load_object( ifile, ids + i, &curr_epoch, &epoch_shown, orbit);

            if( (n_obs_actually_loaded > 1 || !drop_single_obs) && curr_epoch > 0.)
               {
               int n_obs_included = 0;
               unsigned j = 0;

               write_out_elements_to_file( orbit, curr_epoch, epoch_shown,
                     obs, n_obs_actually_loaded, orbit_constraints, element_precision,
                     0, element_options);
               strcpy( tbuff, orbit_summary_text);
               if( use_colors)
                  colorize_text( tbuff);
               if( show_processing_steps)
                  {
                  if( n_processes > 1)
                     {
                     if( process_count == 1 && time( NULL) != update_time)
                        {
                        int elapsed, n_done = i - starting_object + 1;

                        update_time = time( NULL);
                        elapsed = (int)update_time - (int)t0;
                        printf( "%d seconds elapsed, %d remain\n", elapsed,
                               elapsed * (total_objects - n_done) / n_done);
                        }
                     printf( "(%d) %d: %s", process_count, i + 1, ids[i].obj_name);
                     }
                  printf( "; %s ", tbuff);
                  }
               if( separate_residual_file_name)
                  {
                  extern bool residual_file_in_config_dir;

                  residual_file_in_config_dir = false;
                  write_residuals_to_file( separate_residual_file_name, argv[1],
                               n_obs_actually_loaded, obs, RESIDUAL_FORMAT_PRECISE
                               | RESIDUAL_FORMAT_COMPUTER_FRIENDLY
                               | RESIDUAL_FORMAT_FOUR_DIGIT_YEARS
                               | RESIDUAL_FORMAT_EXTRA);
                  residual_file_in_config_dir = true;
                  }
               if( !mpec_path)
                  append_elements_to_element_file = 1;
               if( mpec_path || !is_default_ephem)
                  {
                  int n_orbits_in_ephem = 1;
                  int n_ephemeris_steps = 50;
                  char ephemeris_step_size[20];
                  extern const char *ephemeris_filename;
                  extern const char *residual_filename;
                  extern int available_sigmas;
                  extern double ephemeris_mag_limit;
                  double *orbits_to_use = orbit;
                  const double jd_start = get_time_from_string( curr_jd( ),
                           get_environment_ptr( "EPHEM_START"),
                           CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD
                           | FULL_CTIME_TWO_DIGIT_YEAR, NULL);

                  sscanf( get_environment_ptr( "EPHEM_STEPS"), "%d %9s",
                         &n_ephemeris_steps, ephemeris_step_size);
                  if( ephem_end_jd)
                     {
                     n_ephemeris_steps = 1;
                     snprintf( ephemeris_step_size,
                           sizeof( ephemeris_step_size),
                           "%f", ephem_end_jd - jd_start);
                     }
                  if( !*mpc_code)
                     sscanf( get_environment_ptr( "CONSOLE_OPTS"), "%9s",
                                 mpc_code);
                  create_obs_file( obs, n_obs_actually_loaded, 0);
                  ephemeris_mag_limit = 999.;
                  if( available_sigmas == COVARIANCE_AVAILABLE)
                     {
                     n_orbits_in_ephem = 2;
                     compute_variant_orbit( orbit + 6, orbit, 1.);
                     }
                  if( available_sigmas == SR_SIGMAS_AVAILABLE)
                     {
                     extern double *sr_orbits;
                     extern unsigned n_sr_orbits;

                     orbits_to_use = sr_orbits;
                     n_orbits_in_ephem = n_sr_orbits;
                     }
                  if( !ephemeris_in_a_file_from_mpc_code( ephemeris_filename,
                              orbits_to_use, obs, n_obs_actually_loaded,
                              curr_epoch, jd_start, ephemeris_step_size,
                              n_ephemeris_steps, mpc_code,
                              ephemeris_output_options,
                              n_orbits_in_ephem))
                     {
                     write_residuals_to_file( residual_filename, argv[1],
                                    n_obs_actually_loaded, obs, RESIDUAL_FORMAT_SHORT);
                     if( mpec_path)
                        {
                        char fullpath[100];

                        sprintf( fullpath, "%s/%s.htm", mpec_path, ids[i].packed_desig);
                        text_search_and_replace( fullpath, " ", "");

                        make_pseudo_mpec( fullpath, ids[i].obj_name);
                        get_summary_info( tbuff, fullpath);
                        if( summary_ofile)
                           {
                           FILE *ephemeris_ifile = fopen_ext( ephemeris_filename, "tfcrb");
                           char new_line[420];

                           tbuff[14] = '\0';
                           snprintf( new_line, sizeof( new_line), "<a href=\"%s\">%s</a>%s",
                                    fullpath, tbuff, tbuff + 15);
                           memset( tbuff, 0, sizeof( tbuff));
                           while( j < 4 && fgets_trimmed( tbuff, sizeof( tbuff),
                                                         ephemeris_ifile))
                              j++;
                           if( j == 4)
                              {
                              tbuff[23] = tbuff[39] = tbuff[73] = '\0';
                              sprintf( new_line + strlen( new_line), "  %s  %s  %s",
                                       tbuff + 15, tbuff + 30, tbuff + 57);
                                             /* now add sigma from end of ephem: */
                              while( fgets_trimmed( tbuff, sizeof( tbuff), ephemeris_ifile))
                                 ;
                              tbuff[73] = '\0';
                              sprintf( new_line + strlen( new_line), "%s", tbuff + 68);
                              }
                           fclose( ephemeris_ifile);
                           summary_lines[n_lines_written]
                                      = (char *)malloc( strlen( new_line) + 1);
                           strcpy( summary_lines[n_lines_written], new_line);
                           n_lines_written++;
                           }
                        }
                     }
                  }

               for( j = 0; j < (unsigned)n_obs_actually_loaded; j++)
                  if( obs[j].is_included)
                     n_obs_included++;
               if( n_obs_included != n_obs_actually_loaded
                              && show_processing_steps)
                  {
                  if( use_colors)        /* reverse colors to draw attn */
                     printf( "\033[30;47m");
                  printf( " %d /", n_obs_included);
                  }
               }
            else
               printf( "; not enough observations\n");
            unload_observations( obs, n_obs_actually_loaded);
            if( ephemeris_output_options & OPTION_COMPUTER_FRIENDLY)
               if( mpec_path || !is_default_ephem)
                  have_json_ephem = true;
            add_json_data( "total.json", have_json_ephem, obs->packed_id,
                     i == starting_object + total_objects - 1);
            }
         object_comment_text( tbuff, ids + i);
                  /* Abbreviate 'observations:' to 'obs:' */
         text_search_and_replace( tbuff, "ervations", "");
         if( show_processing_steps)
            printf( "  %s\n", tbuff);
         if( use_colors)
            printf( "\033[0m");
         }
   free( ids);
   if( summary_ofile)
      {
      int pass;

      for( pass = 0; pass < 4; pass++)
         {
         const char *headers[4] = {
               "<h3><a name=\"desig\">Sorted by designation</a></h3>",
               "<h3><a name=\"unc\">Sorted by current ephemeris uncertainty</a></h3>",
               "<h3><a name=\"moid\">Sorted by Earth MOID</a></h3>",
               "<h3><a name=\"mag\">Sorted by magnitude</a></h3>" };
         const char *ephem_header =
             "Object            Semimajor axis    Eccentricity        Incl         "
             "MOID           RA        dec      Elong mag  Sig1  Sig2";

         fprintf( summary_ofile, "\n%s\n", headers[pass]);
         fprintf( summary_ofile, "%s\n\n", ephem_header);
         if( pass)
            {
            static int columns[] = { 0, 123, 77, 118 };

            summ_sort_column = columns[pass];
            qsort( summary_lines, n_lines_written, sizeof( char *), summ_compare);
            }
         for( i = 0; summary_lines[i]; i++)
            {
            fprintf( summary_ofile, "%s%s", summary_lines[i],
                ( i % 5 == 4) ? "\n\n" : "\n");
            }
         }
      for( i = 0; summary_lines[i]; i++)
         free( summary_lines[i]);
      free( summary_lines);
      fprintf( summary_ofile, "</body> </html>\n");
      fclose( summary_ofile);
      }
   fclose( ifile);
#ifdef FORKING
   if( show_processing_steps)
      printf( "Process %d is done\n", process_count);
   wait( &child_status); /* wait for child to exit, and store its status */
   if( process_count == 1)
      {
      extern const char *elements_filename;
      extern const char *mpc_fmt_filename;

      combine_element_files( elements_filename, n_processes);
      combine_element_files( mpc_fmt_filename, n_processes);
      combine_element_files( sof_filename, n_processes);
      combine_element_files( sofv_filename, n_processes);
      for( i = 0; i < n_processes; i++)
         {                             /* clean up temp files: */
         process_count = i + 1;
         unlink_config_file( "monte.txt");
         unlink_config_file( "guide.txt");
         unlink_config_file( "covar.txt");
         unlink_config_file( "sr_elems.txt");
         unlink_config_file( "mpc_sr.txt");
//       unlink_config_file( sof_filename);
         }
      }
#ifdef TEST_PLANET_CACHING_HASH_FUNCTION
   if( process_count == 0)
      {
      extern int64_t planet_ns;
      extern long total_n_searches, total_n_probes, max_probes_required;

      printf( "  tp:%ld.%09ld\n",
                  (long)( planet_ns / (int64_t)1000000000),
                  (long)( planet_ns % (int64_t)1000000000));
      printf( "%ld searches; %ld probes; %ld max probes\n",
                  total_n_searches, total_n_probes, max_probes_required);
      }
#endif
#endif
   clean_up_find_orb_memory( );
   return( 0);
}
