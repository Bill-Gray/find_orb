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

int debug_level = 0;

char *get_file_name( char *filename, const char *template_file_name);
int sanity_test_observations( const char *filename);
int debug_printf( const char *format, ...);                /* runge.cpp */
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
int get_defaults( int *ephemeris_output_options, int *element_format,
         int *element_precision, double *max_residual_for_filtering,
         double *noise_in_arcseconds);                /* elem_out.cpp */
int inquire( const char *prompt, char *buff, const int max_len,
                     const int color);                /* fo.cpp */
void refresh_console( void);                          /* fo.cpp */
void move_add_nstr( const int col, const int row, const char *msg,
                     const int n_bytes);              /* fo.cpp */
double current_jd( void);                       /* elem_out.cpp */
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
int make_pseudo_mpec( const char *mpec_filename, const char *obj_name);
                                               /* ephem0.cpp */

/* In this non-interactive version of Find_Orb,  we just print out warning
messages such as "3 observations were made in daylight" or "couldn't find
thus-and-such file".  These will also be logged in 'debug.txt'.  We then
proceed as if nothing had happened: */

int inquire( const char *prompt, char *buff, const int max_len,
                     const int color)
{
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
}

static int remove_single_observation_objects( OBJECT_INFO *ids, const int n_ids)
{
   int i, j;

   for( i = j = 0; i < n_ids; i++)
      if( ids[i].n_obs > 1)
         ids[j++] = ids[i];
   return( j);
}

#ifdef FORKING
static void combine_element_files( const char *filename, const int n_processes)
{
   FILE **input_files = (FILE **)calloc( (size_t)n_processes, sizeof( FILE *));
   FILE *ofile = fopen( filename, "w");
   char buff[400];
   int i, quit;
   extern int process_count;

   assert( ofile);
   for( i = 0; i < n_processes; i++)
      {
      process_count = i + 1;
      input_files[i] = fopen( get_file_name( buff, filename), "r");
      assert( input_files[i]);
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
      err_code = unlink( get_file_name( buff, filename));
      assert( !err_code);
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
   FILE *ifile = fopen( mpec_filename, "rb");
   char ibuff[400], *tptr;
   unsigned i;

   assert( ifile);
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

int main( const int argc, const char **argv)
{
   char tbuff[300];
   char **summary_lines = NULL;
   const char *mpec_path = NULL;
   int n_ids, i, starting_object = 0;
   int n_processes = 1;
   OBJECT_INFO *ids;
   int heliocentric_only = 1, total_objects = 0;
   FILE *ifile;
   extern int process_count;
   int n_lines_written = 0;
   FILE *summary_ofile = NULL;
   extern int forced_central_body;
#ifdef FORKING
   int child_status;
#endif

               /* get_defaults( ) collects a lot of data that's for the  */
               /* interactive find_orb program.  But it also sets some   */
               /* important internal values for blunder detection,  etc. */
               /* So we still call it:                                   */
   get_defaults( NULL, NULL, NULL, NULL, NULL);
   forced_central_body = 0;
   for( i = 1; i < argc; i++)       /* check to see if we're debugging: */
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'a':
               {
               extern int separate_periodic_comet_apparitions;

               separate_periodic_comet_apparitions ^= 1;
               }
            case 'c':
               {
               extern int combine_all_observations;

               combine_all_observations = 1;
               }
               break;
            case 'd':
               debug_level = atoi( argv[i] + 2);
               if( !debug_level)
                  debug_level = 1;
               debug_printf( "fo: debug_level = %d; %s %s\n",
                           debug_level, __DATE__, __TIME__);
               break;
            case 'h':                     /* show planet-centric orbits */
               if( argv[i][2])
                  forced_central_body = atoi( argv[i] + 2);
               else
                  heliocentric_only = 0;
               break;
            case 'i':
               {
               extern int ignore_prev_solns;

               ignore_prev_solns = 1;
               }
               break;
            case 'm':
               mpec_path = argv[i] + 2;
               break;
            case 'n':
               starting_object = atoi( argv[i] + 2);
               break;
            case 'p':
               n_processes = atoi( argv[i] + 2);
               break;
#ifdef FORKING
            case 'r':
               {                    /* set 'soft' & 'hard' limits for CPU */
               struct rlimit r;     /* run time,  in seconds,  to avoid   */
               int soft_limit, hard_limit;      /* runaway processes */

               if( sscanf( argv[i] + 2, "%d,%d", &soft_limit, &hard_limit) == 2)
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
               summary_ofile = fopen( argv[i] + 2, "wb");
               assert( summary_ofile);
               ifile = fopen( "summ.htm", "rb");
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
               total_objects = atoi( argv[i] + 2);
               break;
            case 'z':
               {
               extern int alt_mpcorb;

               alt_mpcorb = 1;
               }
               break;
            default:
               printf( "Unknown command-line option '%s'\n", argv[i]);
               return( -1);
            }

   load_up_sigma_records( "sigma.txt");
   if( debug_level)
      debug_printf( "Default uncertainties table read\n");

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

   n_ids = remove_single_observation_objects( ids, n_ids);
   printf( "Processing %d objects\n", n_ids);
   if( !total_objects)
      total_objects = n_ids;

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
   printf( "Process count %d\n", process_count);
#endif

   if( summary_ofile)
      summary_lines = (char **)calloc( n_ids - starting_object + 1,
                                                    sizeof( char *));
   ifile = fopen( argv[1], "rb");
   for( i = starting_object; i < n_ids && i < starting_object + total_objects; i++)
      if( n_processes == 1 || i % n_processes == process_count - 1)
         {
         const char *orbit_constraints = "";
         OBSERVE FAR *obs;
         const int n_obs = ids[i].n_obs;

         if( n_processes > 1)
            printf( "(%d) ", process_count);
         printf( "%d: %s", i + 1, ids[i].obj_name);
         if( n_obs < 2)
            printf( "; skipping\n");
         else
            {
            extern int append_elements_to_element_file;
            extern int n_obs_actually_loaded;
            extern char orbit_summary_text[];
            long file_offset = ids[i].file_offset - 40L;
            const int element_precision = 5;
            const int element_options = ELEM_OUT_ALTERNATIVE_FORMAT
                       | (mpec_path ? 0 : heliocentric_only);
            double epoch_shown, curr_epoch, orbit[12];

                /* Start a bit ahead of the actual data,  just in case */
                /* there's a #Sigma: or similar command in there: */
            if( file_offset < 0L)
               file_offset = 0L;
            fseek( ifile, file_offset, SEEK_SET);
            obs = load_object( ifile, ids + i, &curr_epoch, &epoch_shown, orbit);

            if( n_obs_actually_loaded > 1 && curr_epoch > 0.)
               {
               write_out_elements_to_file( orbit, curr_epoch, epoch_shown,
                     obs, n_obs, orbit_constraints, element_precision,
                     0, element_options);
               printf( "; %s ", orbit_summary_text);
               if( !mpec_path)
                  append_elements_to_element_file = 1;
               else
                  {
                  char fullpath[100];
                  int n_orbits_in_ephem = 1;
                  extern const char *ephemeris_filename;
                  extern const char *residual_filename;
                  extern int available_sigmas;
                  extern double ephemeris_mag_limit;
                  double *orbits_to_use = orbit;
                  const double jd_start = get_time_from_string( 0., "now",
                           CALENDAR_JULIAN_GREGORIAN | FULL_CTIME_YMD
                           | FULL_CTIME_TWO_DIGIT_YEAR, NULL);

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
                  sprintf( fullpath, "%s/%s.htm", mpec_path, ids[i].packed_desig);
                  text_search_and_replace( fullpath, " ", "");
                  if( !ephemeris_in_a_file_from_mpc_code( ephemeris_filename,
                              orbits_to_use, obs, n_obs_actually_loaded,
                              curr_epoch, jd_start, "1h", 50, "500",
                              OPTION_SHOW_SIGMAS | OPTION_ROUND_TO_NEAREST_STEP,
                              n_orbits_in_ephem))
                     {
                     write_residuals_to_file( residual_filename, argv[1],
                                    n_obs_actually_loaded, obs, RESIDUAL_FORMAT_SHORT);
                     make_pseudo_mpec( fullpath, ids[i].obj_name);
                     get_summary_info( tbuff, fullpath);
                     if( summary_ofile)
                        {
                        FILE *ephemeris_ifile = fopen( ephemeris_filename, "rb");
                        char new_line[300];
                        unsigned j = 0;

                        tbuff[14] = '\0';
                        sprintf( new_line, "<a href=\"%s\">%s</a>%s",
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
            else
               printf( "; not enough observations\n");
            unload_observations( obs, n_obs);
            }
         object_comment_text( tbuff, ids + i);
                  /* Abbreviate 'observations:' to 'obs:' */
         text_search_and_replace( tbuff, "ervations", "");
         printf( "  %s\n", tbuff);
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
   printf( "Process %d is done\n", process_count);
   wait( &child_status); /* wait for child to exit, and store its status */
   if( process_count == 1)
      {
      extern const char *elements_filename;
      extern const char *mpc_fmt_filename;

      combine_element_files( elements_filename, n_processes);
      combine_element_files( mpc_fmt_filename, n_processes);
      for( i = 0; i < n_processes; i++)
         {                             /* clean up temp files: */
         process_count = i + 1;
         unlink( get_file_name( tbuff, "monte.txt"));
         unlink( get_file_name( tbuff, "guide.txt"));
         unlink( get_file_name( tbuff, "covar.txt"));
         unlink( get_file_name( tbuff, "sr_elems.txt"));
         unlink( get_file_name( tbuff, "mpc_sr.txt"));
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
