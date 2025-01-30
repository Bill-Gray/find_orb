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

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <assert.h>
#include <errno.h>
#include "stringex.h"
#include "mpc_func.h"

int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;

#ifdef CONFIG_DIR_AUTOCOPY
#include <string>

/* Support older (pre-10.15 Catalina) versions of macOS via cpp-filesystem
   shim. From https://github.com/gulrak/filesystem#using-it-as-single-file-header
*/
#ifdef __APPLE__
#include <Availability.h> // for deployment target to support pre-catalina targets without std::fs
#endif
#if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || (defined(__cplusplus) && __cplusplus >= 201703L)) && defined(__has_include)
#if __has_include(<filesystem>) && (!defined(__MAC_OS_X_VERSION_MIN_REQUIRED) || __MAC_OS_X_VERSION_MIN_REQUIRED >= 101500)
#define GHC_USE_STD_FS
#include <filesystem>
namespace fs = std::filesystem;
#endif
#endif
#ifndef GHC_USE_STD_FS
#include <ghc/filesystem.hpp>
namespace fs = ghc::filesystem;
#endif

#endif

#if defined( _WIN32) || defined( __WATCOMC__)
   #include <direct.h>        /* for _mkdir() definition */
#else
   #include <sys/stat.h>
   #include <sys/types.h>
   #include <unistd.h>
#endif

const char *get_find_orb_text( const int index);

/* This function allows one to put the following options in front of
the 'permits' string :

t  File is temporary (doesn't actually have an effect yet).
f  If the file isn't opened,  a fatal error message is shown.
c  Configuration file:  look for it in ~/.find_orb (for Linux) or in
   the directory in which Find_Orb is running (Windows).
l  Try opening locally;  if that fails,  try the config directory.

   One can combine these.  For example, 'permits' of tfcw would
tell the function that the file is a temporary one,  should be
opened within the configuration directory for writing,  and that
if it can't be opened, it's a fatal error.  'cl' = try config,
then local; 'lc' = try local,  then config;  'c' = try config
only.

   I'm still working out the degree to which a separate directory
for these files will be needed.  In the past,  find_orb,  fo,
and fo_serve.cgi simply read what they needed from the current
folder.  At present,  the only case in which use_config_directory
is 'true' is for console Find_Orb in Linux,  and even there,  you
can turn it back to 'false'.
*/

#ifndef _WIN32
int get_temp_dir( char *name, const size_t max_len);      /* miscell.cpp */
#endif
int fetch_astrometry_from_mpc( FILE *ofile, const char *desig);
int download_a_file( const char *ofilename, const char *url);
int generic_message_box( const char *message, const char *box_type);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
char *make_config_dir_name( char *oname, const char *iname);  /* miscell.cpp */
int reset_astrometry_filename( int *argc, const char **argv);
uint64_t parse_bit_string( const char *istr);                /* miscell.cpp */
const char *write_bit_string( char *ibuff, const uint64_t bits,
                                          const size_t max_bitstring_len);
int use_config_directory = false;
const char *alt_config_directory;

#ifndef CONFIG_DIR_AUTOCOPY

// make this an empty function if not building for conda
void ensure_config_directory_exists()
{
}

#else

// HACK: create the default config directory from a central copy
// Going forward it'd be good to reqork Find_Orb to search for data in
// default directories if local copies don't exist
#include "prefix.h"
static char PREFIX_STATIC[500] = PREFIX;
void ensure_config_directory_exists()
{
   if (!use_config_directory)
      return;

   // The c_str() magic in the next line allows conda-build's prefix
   // replacer to work as expected.
   // See https://github.com/conda/conda-build/issues/1674 for details.
   // Modified on 2023-03-17 to add PREFIX_STATIC as compilers have
   // gotten too clever & the existing workaround stopped working.
   std::string prefix = std::string(PREFIX_STATIC).c_str();

   if (prefix == "~") {
      // backwards compatibility; do nothing.
      return;
   }

   const char *home = getenv("HOME");
   if (home == NULL) {
      // home unknown; give up
      return;
   }
   auto dest = fs::path(home) / ".find_orb";
   if (fs::exists(dest)) {
      // ~/.find_orb already exists; nothing to do
      return;
   }

   // copy from $PREFIX/share/findorb/data to ~/.find_orb
   auto src = fs::path(prefix) / "share" / "findorb" / "data";
   fs::copy(src, dest, fs::copy_options::recursive);

   // symlink the ephemerides files from system-wide
   // location, if there are any
   auto path = fs::path(prefix) / "share" / "findorb" / "jpl_eph";
   if (fs::exists(path)) {
      for (const auto & entry : fs::directory_iterator(path))
      {
         auto lnk = dest / entry.path().filename();
         fs::create_symlink(entry, lnk);
      }
   }
}
#endif

/* Users may specify files such as ~/this/that.txt on non-Windows boxes.
The following function replaces ~ with the home directory. */

static FILE *fopen_tilde( const char *filename, const char *permits)
{
#ifdef _WIN32
   return( fopen( filename, permits));
#else
   if( *filename != '~' || filename[1] != '/')
      return( fopen( filename, permits));
   else
      {
      char fullname[255];

      strlcpy_err( fullname, getenv( "HOME"), sizeof( fullname));
      strlcat_err( fullname, filename + 1, sizeof( fullname));
      return( fopen( fullname, permits));
      }
#endif
}

char *make_config_dir_name( char *oname, const char *iname)
{
#ifndef _WIN32
   char *home_ptr = getenv( "HOME");
#endif

   if( alt_config_directory && *alt_config_directory)
      {
      strcpy( oname, alt_config_directory);
      strcat( oname, iname);
      return( oname);
      }
#ifdef _WIN32
   strcpy( oname, iname);
#else
   if( home_ptr)
      {
      strcpy( oname, home_ptr);
      strcat( oname, "/.find_orb/");
      }
   else
      *oname = '\0';
   strcat( oname, iname);
#endif
   return( oname);
}

const char *output_directory = NULL;

#ifndef _WIN32
int get_temp_dir( char *name, const size_t max_len)
{
   static int process_id = 0;

   if( output_directory)
      strlcpy_err( name, output_directory, max_len);
   else
      {
      const bool first_time = (process_id == 0);

      if( first_time)
#if defined( _WIN32) || defined( __WATCOMC__)
         process_id = 1;
#else
         process_id = getpid( );
#endif
      snprintf_err( name, max_len, "/tmp/find_orb%d", process_id);
      if( first_time)
#if defined( _WIN32) || defined( __WATCOMC__)
         _mkdir( name);
#else
         mkdir( name, 0777);
#endif
      }
   return( process_id);
}
#endif

/* We use a lock file to determine if Find_Orb is already running,  and
therefore putting some temporary files (ephemerides,  elements,  etc.)
into the config directory ~/.find_orb.   If that's happening,  we ought
to put our temporary files elsewhere,  in a directories of the form
/tmp/find_orb(process ID).

   We can also deliberately put output files to a desired 'output_directory',
via command-line options.    */

FILE *fopen_ext( const char *filename, const char *permits)
{
   FILE *rval = NULL;
   bool is_fatal = false;
   bool try_local = true;
   bool is_temporary = false;

   if( *permits == 't')
      {
#if !defined( _WIN32) && !defined( __WATCOMC__)
      extern bool findorb_already_running;

      is_temporary = findorb_already_running || (output_directory != NULL);
#endif
      permits++;
      }
   if( *permits == 'f')
      {
      is_fatal = true;
      permits++;
      }
   if( !use_config_directory || is_temporary)
      while( *permits == 'l' || *permits == 'c')
         permits++;
   if( *permits == 'c' && strchr( filename, '/'))
      permits++;              /* not really in the config directory */
   if( permits[0] == 'l' && permits[1] == 'c')
      {     /* try local,  then config version */
      try_local = false;
      permits++;
      rval = fopen_tilde( filename, permits + 1);
      }
#ifndef _WIN32
   if( is_temporary)
      {
      char tname[255];

      get_temp_dir( tname, sizeof( tname));
      snprintf_append( tname, sizeof( tname), "/%s",  filename);
      rval = fopen( tname, permits);
      }
#endif
   if( !rval && *permits == 'c' && !is_temporary)
      {
      char tname[255];

      make_config_dir_name( tname, filename);
      permits++;
      if( *permits == 'l')       /* permits are 'cl' = check both */
         permits++;
      else                       /* check config version only */
         try_local = false;
      rval = fopen_tilde( tname, permits);
      }
   if( try_local && !rval && !is_temporary)
      rval = fopen_tilde( filename, permits);
   if( !rval && is_fatal)
      {
      char buff[300];

      snprintf( buff, sizeof( buff), "Error opening %s: %s",
                 filename, strerror( errno));
      generic_message_box( buff, "o");
      exit( -1);
      }
   return( rval);
}

static int desig_matches( const char *iline, const char *desig)
{
   const size_t len = strlen( desig);
   int rval = 0;

   while( *iline == ' ')
      iline++;
   if( !memcmp( iline, desig, len))
      if( iline[len] == ' ' || iline[len] == '*')
         rval = 1;
   return( rval);
}

static int is_neocp_desig( const char *buff)
{
   const size_t len = strlen( buff);
   int rval = 0;

   if( len > 2 && len < 8 && !strchr( buff, ' '))
      {
      int unused, n_bytes;

      rval = 1;
      if( sscanf( buff, "%d%n", &unused, &n_bytes) == 1
                  && n_bytes == (int)len)    /* oops!  numbered object */
         rval = 0;
      }
   return( rval);
}

   /* Code to look for astrometry for a particular object,  usually in
      'neocp.txt' or 'neocp.old' (files of current or previous NEOCP
      astrometry,  respectively.)         */

static int look_for_astrometry_in_file( FILE *ofile, const char *ifilename, const char *desig)
{
   int bytes_written = 0;
   FILE *ifile;

   if( ifilename && *ifilename && (ifile = fopen( ifilename, "rb")) != NULL)
      {
      char tbuff[100];

      while( fgets( tbuff, sizeof( tbuff), ifile))
         if( desig_matches( tbuff, desig))
            bytes_written += (int)fwrite( tbuff, 1, strlen( tbuff), ofile);
      fclose( ifile);
      }
   return( bytes_written);
}

/* If astrometry is desired for a particular designation,  we hunt for it
in several possible locations.  If it looks like a temporary designation,
we look through up to five files specified by NEOCP_FILE_NAME0 through
NEOCP_FILE_NAME4.  For both the on-line Find_Orb and my own machine,
the first two correspond to a file containing astrometry for objects
currently on NEOCP and an 'neocp.old'-file for objects that have been
removed from NEOCP.  The others allow me to put oddball objects into
files for the on-line Find_Orb.

   If none of those files supplies astrometry for the object,  we use
the 'grab_mpc' program :

https://github.com/Bill-Gray/miscell/blob/master/grab_mpc.c

to download astrometry from MPC.  As described in 'grab_mpc.c',  some
steps are taken to cache downloads to avoid hammering MPC servers. The
temp*.ast files are time-stamped,  and we only get new data if three
hours have elapsed;  see 'grab_mpc.c' for details.

   (On Windows(R),  we may be able to use the URLDownloadToFile()
function (q.v.) to handle such downloads.)

   Note that your average user probably hasn't set up these various
files to provide NEOCP astrometry nor 'grab_mpc',  which will cause
this code to always return 0 (i.e.,  no astrometry fetched.)
*/

#define ERR_CODE_NO_DATA_AVAILABLE  (int)0xd100

int fetch_astrometry_from_mpc( FILE *ofile, const char *desig)
{
   char tbuff[300];
   int bytes_written = 0, pass;
   const char *grab_program = get_environment_ptr( "MPC_GRAB_PROGRAM");

   if( !*grab_program)
      grab_program = "grab_mpc";
   assert( ofile);
   if( is_neocp_desig( desig))
      for( pass = 0; pass < 5 && !bytes_written; pass++)
         {
         snprintf( tbuff, sizeof( tbuff), "NEOCP_FILE_NAME%d", pass);
         bytes_written = look_for_astrometry_in_file( ofile,
                    get_environment_ptr( tbuff), desig);
         }

   if( !bytes_written && *grab_program)
      {                 /* no NEOCP data;  maybe MPC has astrometry */
      unsigned j = 0;
      int err_code;
      size_t i;
      char filename[40], temp_desig[10];

      if( isalpha( desig[0]) && isdigit( desig[1]) && isdigit( desig[2])
                             && isdigit( desig[3]) && isdigit( desig[4]) && !desig[5])
         {
         snprintf_err( temp_desig, sizeof( temp_desig), "%d%s",
                  mutant_hex_char_to_int( *desig), desig + 1);
         desig = temp_desig;     /* unpacks G3141 -> 163141,  etc. */
         }
      for( i = 0; desig[i]; i++)
         j = j * 314159u + (unsigned)desig[i];
#if defined( _WIN32) || defined( __WATCOMC__)
      snprintf_err( filename, sizeof( filename),      "temp%02u.ast", j % 100);
#else
      snprintf_err( filename, sizeof( filename), "/tmp/temp%02u.ast", j % 100);
#endif
      snprintf_err( tbuff, sizeof( tbuff), "%s %s %s", grab_program,
                                                filename, desig);
      err_code = system( tbuff);
      if( !err_code)
         {
         FILE *ifile = fopen( filename, "rb");

         assert( ifile);
         while( fgets( tbuff, sizeof( tbuff), ifile))
            bytes_written += (int)fwrite( tbuff, 1, strlen( tbuff), ofile);
         fclose( ifile);
         }
      else if( err_code == ERR_CODE_NO_DATA_AVAILABLE)
         {
         snprintf_err( tbuff, sizeof( tbuff), get_find_orb_text( 2080), desig);
         generic_message_box( tbuff, "o");
         }
      else
         {
         debug_printf( "grab_mpc error code %d (%x) %s\n", err_code, err_code,
                     strerror( err_code));
         generic_message_box( get_find_orb_text( 2058), "o");
         }
      }
   return( bytes_written);
}

int download_a_file( const char *ofilename, const char *url)
{
   const char *grab_program = get_environment_ptr( "MPC_GRAB_PROGRAM");
   char *tbuff;
   size_t buffsize;
   int err_code;

   if( !*grab_program)
      grab_program = "grab_mpc";
   buffsize = strlen( grab_program) + strlen( ofilename) + strlen( url) + 4;
   tbuff = (char *)malloc( buffsize);
   snprintf_err( tbuff, (int)buffsize, "%s %s %s", grab_program, ofilename, url);
   err_code = system( tbuff);
   free( tbuff);
   return( err_code);
}

/* Code to write a single valid,  but completely meaningless observation
to a temporary file for a specified object.  This can then be read in
and used to get around the fact that Find_Orb essentially assumes that
there's astrometry available for the object(s) under consideration.
(If you search for 'Dummy' in elem_ou2.cpp,  you'll see that the date/time
is replaced with the epoch of the elements,  and the RA/dec is replaced
with that for the object at the epoch.  The values given below are just
placeholders to ensure that one valid observation gets read in.) */

static int make_fake_astrometry( const char *obj_name, const char *filename)
{
   FILE *ofile = fopen( filename, "wb");
   int rval = -1;

   assert( ofile);
   if( ofile)
      {
      const char *dummy_line =
        "  C2000 01 01.00000 00 00 00.00 +00 00 00.0               V Dummy500";
      char packed_desig[20];

      create_mpc_packed_desig( packed_desig, obj_name);
      fprintf( ofile, "%s%s\n", packed_desig, dummy_line);
      fclose( ofile);
      rval = 0;
      }
   return( rval);
}

/* The various flavors of Find_Orb have some similar logic in place
to allow the programs to either download astrometry,  or to create
a dummy file to generate ephemerides from stored orbital elements.

A command-line argument '-f' followed by an MPC designation causes the
desired astrometry to be fetched from MPC,  and an orbit is computed.
For a command-line argument '-o' followed by an MPC designation,  the
orbital elements are _read_ (not computed) from 'mpcorb.sof',  and a
single dummy observation is made that reflects the object's position
at the epoch (see above 'make_fake_astrometry()' function).  */

char *temp_obs_filename;

int reset_astrometry_filename( int *argc, const char **argv)
{
   int rval = 0;

   temp_obs_filename = (char *)malloc( 260);
   make_config_dir_name( temp_obs_filename, "temp_obs.txt");
   if( *argc > 1 && argv[1][0] == '-'
                    && (argv[1][1] == 'o' || argv[1][1] == 'f'))
      {
      char obj_name[50];
      int i, j;

      if( argv[1][2])
         strcpy( obj_name, argv[1] + 2);
      else
         *obj_name = '\0';
      for( i = 2; i < *argc && argv[i][0] != '-' && !strchr( argv[i], '='); i++)
         {
         if( *obj_name)
            strcat( obj_name, " ");
         strcat( obj_name, argv[i]);
         argv[i] = "";
         }
      if( argv[1][1] == 'f')
         {
         FILE *ofile = fopen_ext( temp_obs_filename, "fwb");

         assert( ofile);
         fetch_astrometry_from_mpc( ofile, obj_name);
         fclose( ofile);
         }
      else
         {
         const size_t len = strlen( obj_name);
         char unpacked[50], aligned[13];

         if( len > 4 && len < 9 && !strchr( obj_name, ' '))
            {
            memset( aligned, ' ', 12);
            aligned[12] = '\0';
            if( len == 5)
               memcpy( aligned, obj_name, 5);
            else
               memcpy( aligned + 12 - len, obj_name, len);
            if( unpack_mpc_desig( unpacked, aligned) != OBJ_DESIG_OTHER)
               strlcpy_error( obj_name, unpacked);
            }
         make_fake_astrometry( obj_name, temp_obs_filename);
         }
      argv[1] = temp_obs_filename;
      for( j = 2; i < *argc; j++, i++)
         argv[j] = argv[i];
      argv[j] = NULL;
      *argc = j;
      rval = 1;
      }
   return( rval);
}

/* Reads a string such as,  say,  "1,12,3-6,9" and returns a 64-bit
integer with (in this case) bits 1, 12,  3-6,  and 9 turned on:
hex 127A = binary 0001 0010 0111 1010.  Note that bits are toggled,
so "4-9,6" would turn bits 4-9 on,  then turn off bit 6.  This also
means that the order is unimportant;  "6,4-9" would give the same result. */

uint64_t parse_bit_string( const char *istr)
{
   uint64_t rval = 0;
   int bytes_scanned, bit1;
   int prev_bit1 = -1;

   while( sscanf( istr, "%d%n", &bit1, &bytes_scanned) == 1)
      {
      assert( bit1 >= 0 && bit1 < 64);
      if( prev_bit1 >= 0)     /* we're doing a range here */
         rval ^= (((uint64_t)2 << bit1) - ((uint64_t)2 << prev_bit1));
      else
         rval ^= (uint64_t)1 << bit1;
      istr += bytes_scanned;
      prev_bit1 = -1;
      if( *istr == ',')
         istr++;
      else if( *istr == '-')
         {
         prev_bit1 = bit1;
         istr++;
         }
      }
   return( rval);
}

/* Inverse of the above function.  I _think_ the longest possible bitstring
length is 122 bytes :

0-1,3-4,6-7,9,11-12,14-15,17-18,20-21,23-24,26-27,29-30,32-33,35-36,38-40,41-42,44-45,47-48,50-51,53-54,56-57,59-60,62-63

   but I lack a formal proof of this.  I'm sure 256 bytes will be more
than enough.   */

const char *write_bit_string( char *ibuff, const uint64_t bits, const size_t max_bitstring_len)
{
   int i, j;

   *ibuff = '\0';
   for( i = 0; i < 64; i++)
      if( (bits >> i) & 1)
         {
         j = i + 1;
         while( (bits >> j) & 1)
            j++;
         if( j > i + 1)
            snprintf_append( ibuff, max_bitstring_len, "%d-%d,", i, j - 1);
         else
            snprintf_append( ibuff, max_bitstring_len, "%d,", i);
         i = j;
         }
   if( *ibuff)             /* remove trailing comma */
      ibuff[strlen( ibuff) - 1] = '\0';
   return( ibuff);
}

/* Very simple pattern matcher,  using ? and * only */

int pattern_match(const char* pattern, const char* string)
{
   while( pattern[0])
      if( pattern[0] == '*')
            return pattern_match(pattern+1, string)
                     || (string[0] && pattern_match(pattern, string+1));
      else
         {
         if( pattern[0] == '?' && !string[0])
            return 0;
         if( pattern[0] != '?' && pattern[0] != string[0])
            return 0;
         pattern++;
         string++;
         }
   return !string[0];
}

const char *find_orb_version_jd( double *jd)
{
    if( jd)
      *jd = 2460705.5;
    return( "2025 Jan 30");
}
