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
#include <assert.h>
#include <errno.h>

   /* MSVC/C++ lacks snprintf.  See 'ephem0.cpp' for details. */
#if defined(_MSC_VER) && _MSC_VER < 1900
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

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
int fetch_astrometry_from_mpc( FILE *ofile, const char *desig); /* miscell.c */

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

FILE *fopen_ext( const char *filename, const char *permits)
{
   FILE *rval = NULL;
   bool is_fatal = false;
   bool try_local = true;

   if( *permits == 't')
      permits++;
   if( *permits == 'f')
      {
      is_fatal = true;
      permits++;
      }
   if( !use_config_directory)
      while( *permits == 'l' || *permits == 'c')
         permits++;
   if( permits[0] == 'l' && permits[1] == 'c')
      {     /* try local,  then config version */
      try_local = false;
      permits++;
      rval = fopen( filename, permits + 1);
      }
   if( !rval && *permits == 'c')
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
   if( try_local && !rval)
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

int fetch_astrometry_from_mpc( FILE *ofile, const char *desig)
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
