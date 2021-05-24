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

const char *get_find_orb_text( const int index);
#endif

size_t strlcpy( char *dst, const char *src, size_t dsize);   /* miscell.c */
size_t strlcat( char *dst, const char *src, size_t dsize);   /* miscell.c */
size_t strlcpy_err( char *dst, const char *src, size_t dsize); /* miscell.c */
size_t strlcat_err( char *dst, const char *src, size_t dsize); /* miscell.c */

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

int fetch_astrometry_from_mpc( FILE *ofile, const char *desig);
int generic_message_box( const char *message, const char *box_type);
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
void make_config_dir_name( char *oname, const char *iname);  /* miscell.cpp */
int reset_astrometry_filename( int *argc, const char **argv);
uint64_t parse_bit_string( const char *istr);                /* miscell.cpp */
const char *write_bit_string( char *ibuff, const uint64_t bits);

int use_config_directory = false;
const char *alt_config_directory;

/* Users may specify files such as ~/this/that.txt on non-Windows boxes.
The following function replaces ~ with the home directory. */

static FILE *fopen_tilde( const char *filename, const char *permits)
{
#ifdef _WIN32
   return( fopen( filename, permits));
#else
   if( *filename != '~' && filename[1] != '/')
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
      strcpy( oname, home_ptr);
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
      rval = fopen_tilde( filename, permits + 1);
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
      rval = fopen_tilde( tname, permits);
      }
   if( try_local && !rval && !is_temporary)
      rval = fopen_tilde( filename, permits);
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

int fetch_astrometry_from_mpc( FILE *ofile, const char *desig)
{
   char tbuff[100];
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
      size_t i;
      char filename[40];

      for( i = 0; desig[i]; i++)
         j = j * 314159u + (unsigned)desig[i];
#ifdef _WIN32
      snprintf( filename, sizeof( filename),      "temp%02u.ast", j % 100);
#else
      snprintf( filename, sizeof( filename), "/tmp/temp%02u.ast", j % 100);
#endif
      snprintf( tbuff, sizeof( tbuff), "%s %s %s", grab_program,
                                                filename, desig);
      if( !system( tbuff))
         {
         FILE *ifile = fopen( filename, "rb");

         assert( ifile);
         while( fgets( tbuff, sizeof( tbuff), ifile))
            bytes_written += (int)fwrite( tbuff, 1, strlen( tbuff), ofile);
         fclose( ifile);
         }
#ifndef _WIN32
      else
         generic_message_box( get_find_orb_text( 2058), "o");
#endif
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

const char *temp_obs_filename = "temp_obs.txt";

int reset_astrometry_filename( int *argc, const char **argv)
{
   int rval = 0;

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
         FILE *ofile = fopen_ext( temp_obs_filename, "twb");

         assert( ofile);
         fetch_astrometry_from_mpc( ofile, obj_name);
         fclose( ofile);
         }
      else
         make_fake_astrometry( obj_name, temp_obs_filename);
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

const char *write_bit_string( char *ibuff, const uint64_t bits)
{
   int i, j;
   const size_t max_bitstring_len = 256;

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

/* strlcat() and strlcpy() appear in some BSDs,  but I don't think they
appear anyplace else (though in my opinion,  they should).   The following
is from http://www.openbsd.org/cgi-bin/cvsweb/src/lib/libc/string/ . */

/* $OpenBSD: strlcpy.c,v 1.16 2019/01/25 00:19:25 millert Exp $   */

/*
 * Copyright (c) 1998, 2015 Todd C. Miller <millert@openbsd.org>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * Copy string src to buffer dst of size dsize.  At most dsize-1
 * chars will be copied.  Always NUL terminates (unless dsize == 0).
 * Returns strlen(src); if retval >= dsize, truncation occurred.
 */
size_t strlcpy( char *dst, const char *src, size_t dsize)
{
   const char *osrc = src;
   size_t nleft = dsize;

   /* Copy as many bytes as will fit. */
   if (nleft != 0) {
      while (--nleft != 0) {
         if ((*dst++ = *src++) == '\0')
            break;
      }
   }

   /* Not enough room in dst, add NUL and traverse rest of src. */
   if (nleft == 0) {
      if (dsize != 0)
         *dst = '\0';      /* NUL-terminate dst */
      while (*src++)
         ;
   }

   return(src - osrc - 1); /* count does not include NUL */
}


/*
 * Appends src to string dst of size dsize (unlike strncat, dsize is the
 * full size of dst, not space left).  At most dsize-1 characters
 * will be copied.  Always NUL terminates (unless dsize <= strlen(dst)).
 * Returns strlen(src) + MIN(dsize, strlen(initial dst)).
 * If retval >= dsize, truncation occurred.
 */

size_t strlcat( char *dst, const char *src, size_t dsize)
{
   const char *odst = dst;
   const char *osrc = src;
   size_t n = dsize;
   size_t dlen;

   /* Find the end of dst and adjust bytes left but don't go past end. */
   while (n-- != 0 && *dst != '\0')
      dst++;
   dlen = dst - odst;
   n = dsize - dlen;

   if (n-- == 0)
      return(dlen + strlen(src));
   while (*src != '\0') {
      if (n != 0) {
         *dst++ = *src;
         n--;
      }
      src++;
   }
   *dst = '\0';

   return(dlen + (src - osrc));  /* count does not include NUL */
}

/* Same as strlcpy() and strlcat(),  except that if truncation occurs,
we abort.  strlcpy()/strlcat() should only be used in situations where
truncation is entirely to be expected;  if truncation indicates a bug,
use the _err() versions instead.  */

size_t strlcpy_err( char *dst, const char *src, size_t dsize)
{
   const size_t rval = strlcpy( dst, src, dsize);

   if( rval >= dsize)
      {
      fprintf( stderr, "strlcpy overflow: dsize = %ld, rval %ld, '%s'\n",
                     (long)dsize, (long)rval, src);
      exit( -1);
      }
   return( rval);
}

size_t strlcat_err( char *dst, const char *src, size_t dsize)
{
   const size_t rval = strlcat( dst, src, dsize);

   if( rval >= dsize)
      {
      fprintf( stderr, "strlcat overflow: dsize = %ld, rval %ld, '%s'\n",
                     (long)dsize, (long)rval, src);
      exit( -1);
      }
   return( rval);
}
