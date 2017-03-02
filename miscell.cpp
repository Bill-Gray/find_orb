#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <errno.h>

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

int use_config_directory = false;

void make_config_dir_name( char *oname, const char *iname)
{
   strcpy( oname, getenv( "HOME"));
   strcat( oname, "/.find_orb/");
   strcat( oname, iname);
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
