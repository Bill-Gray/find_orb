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
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <errno.h>
#include "mpc_func.h"

#ifndef _WIN32
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

   /* MSVC/C++ lacks snprintf.  See 'ephem0.cpp' for details. */
#if defined(_MSC_VER) && _MSC_VER < 1900
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

int snprintf_append( char *string, const size_t max_len,      /* ephem0.cpp */
                                   const char *format, ...)
#ifdef __GNUC__
         __attribute__ (( format( printf, 3, 4)))
#endif
;

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

int generic_message_box( const char *message, const char *box_type);
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
void make_config_dir_name( char *oname, const char *iname);  /* miscell.cpp */
int reset_astrometry_filename( const int argc, const char **argv);
uint64_t parse_bit_string( const char *istr);                /* miscell.cpp */

int use_config_directory = false;
const char *alt_config_directory;

void make_config_dir_name( char *oname, const char *iname)
{
#ifndef _WIN32
   char *home_ptr = getenv( "HOME");
#endif

   if( alt_config_directory && *alt_config_directory)
      {
      strcpy( oname, alt_config_directory);
      strcat( oname, iname);
      return;
      }
#ifdef _WIN32
   strcpy( oname, iname);
#else
   if( home_ptr)
      {
      strcpy( oname, getenv( "HOME"));
      strcat( oname, "/.find_orb/");
      }
   else
      *oname = '\0';
   strcat( oname, iname);
#endif
}

/* We use a lock file to determine if Find_Orb is already running,  and
therefore putting some temporary files (ephemerides,  elements,  etc.)
into the config directory ~/.find_orb.   If that's happening,  we ought
to put our temporary files elsewhere,  in a directories of the form
/tmp/find_orb(process ID).

   We can also deliberately put output files to a desired 'output_directory',
via command-line options.    */

const char *output_directory = NULL;

FILE *fopen_ext( const char *filename, const char *permits)
{
   FILE *rval = NULL;
   bool is_fatal = false;
   bool try_local = true;
   bool is_temporary = false;

   if( *permits == 't')
      {
#ifndef _WIN32
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
   if( permits[0] == 'l' && permits[1] == 'c')
      {     /* try local,  then config version */
      try_local = false;
      permits++;
      rval = fopen( filename, permits + 1);
      }
#ifndef _WIN32
   if( is_temporary)
      {
      char tname[255];
      static int process_id = 0;
      const bool first_time = (process_id == 0);

      if( first_time)
         process_id = getpid( );
      if( output_directory)
         strcpy( tname, output_directory);
      else
         snprintf( tname, sizeof( tname), "/tmp/find_orb%d", process_id);
      if( first_time)
         mkdir( tname, 0777);
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
      rval = fopen( tname, permits);
      }
   if( try_local && !rval && !is_temporary)
      rval = fopen( filename, permits);
   if( !rval && is_fatal)
      {
      char buff[300];

      sprintf( buff, "Error opening %s: %s",
                 filename, strerror( errno));
      generic_message_box( buff, "o");
      exit( -1);
      }
   return( rval);
}

#ifndef strlcpy

      /* Some systems (BSD, Solaris,  others) offer this handy function. */
      /* Linux and Windows don't.                                        */

size_t strlcpy(char *dest, const char *src, size_t size)
{
   size_t i;

   for( i = 0; i < size && *src; i++)
      *dest++ = *src++;
   if( i == size)
      dest--;
   *dest = '\0';
   return( i);
}
#endif

static int desig_matches( const char *iline, const char *desig)
{
   const size_t len = strlen( desig);
   int rval = 0;

   while( *iline == ' ')
      iline++;
   if( !memcmp( iline, desig, len) && iline[len] == ' ')
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


/* If astrometry is desired for a particular designation,  we hunt for it
in several possible locations.  If it looks like a temporary designation,
we look through 'neocp.txt'.  If that doesn't turn up anything,  we use
the 'grab_mpc' program :

https://github.com/Bill-Gray/miscell/blob/master/grab_mpc.c

to download astrometry from MPC.  As described in 'grab_mpc.c',  some
steps are taken to cache downloads to avoid hammering MPC servers. The
temp*.ast files are time-stamped,  and we only get new data if three
hours have elapsed;  see 'grab_mpc.c' for details.       */

static int fetch_astrometry_from_mpc( FILE *ofile, const char *desig)
{
   char tbuff[100];
   int bytes_written = 0;

   assert( ofile);
   if( is_neocp_desig( desig))
      {
      FILE *ifile = fopen( "../../neocp2/neocp.txt", "rb");

      if( ifile)
         {
         while( fgets( tbuff, sizeof( tbuff), ifile))
            if( desig_matches( tbuff, desig))
               bytes_written += (int)fwrite( tbuff, 1, strlen( tbuff), ofile);
         fclose( ifile);
         }
      }
   if( !bytes_written)    /* no NEOCP data;  maybe MPC has astrometry */
      {
      unsigned j = 0;
      size_t i;
      char filename[40];

      for( i = 0; desig[i]; i++)
         j = j * 314159u + (unsigned)desig[i];
      snprintf( filename, sizeof( filename), "/tmp/temp%02u.ast", j % 100);
      snprintf( tbuff, sizeof( tbuff), "grab_mpc %s %s", filename, desig);
      if( !system( tbuff))
         {
         FILE *ifile = fopen( filename, "rb");

         assert( ifile);
         while( fgets( tbuff, sizeof( tbuff), ifile))
            bytes_written += (int)fwrite( tbuff, 1, strlen( tbuff), ofile);
         fclose( ifile);
         }
      }
   return( bytes_written);
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
        "  C2000 01 01.00000 00 00 00.00 +00 00 00.0                 Dummy500";
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
a dummy file to generate ephemerides from stored orbital elements. */

int reset_astrometry_filename( const int argc, const char **argv)
{
   int rval = 0;

   if( argv[1][0] == '-' && (argv[1][1] == 'o' || argv[1][1] == 'f'))
      {
      const char *filename = "/tmp/temp_obs.txt";
      char obj_name[50];
      int i;

      if( argv[1][2])
         strcpy( obj_name, argv[1] + 2);
      else
         *obj_name = '\0';
      for( i = 2; i < argc && argv[i][0] != '-' && !strchr( argv[i], '='); i++)
         {
         if( *obj_name)
            strcat( obj_name, " ");
         strcat( obj_name, argv[i]);
         argv[i] = "";
         }
      if( argv[1][1] == 'f')
         {
         FILE *ofile = fopen( filename, "wb");

         assert( ofile);
         fetch_astrometry_from_mpc( ofile, obj_name);
         fclose( ofile);
         }
      else
         make_fake_astrometry( obj_name, filename);
      argv[1] = filename;
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
