#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/* create_mpc_packed_desig( ) takes a "normal" name for a comet/asteroid,
such as P/1999 Q1a or 2005 FF351,  and turns it into the 12-byte packed
format used in MPC reports and element files.  Documentation of this format
is given on the MPC Web site.  A 'test main' at the end of this file
shows the usage of this function.
   This should handle all "normal" asteroid and comet designations.
It doesn't handle natural satellites.  */

/* "Mutant hex" uses the usual hex digits 0123456789ABCDEF for numbers
0 to 15,  followed by G...Z for 16...35 and a...z for 36...61.  MPC stores
epochs and certain other numbers using this scheme to save space.  */

int create_mpc_packed_desig( char *packed_desig, const char *obj_name);

static inline char mutant_hex( const int number)
{
   int rval;

   if( number < 10)
      rval = '0';
   else if( number < 36)
      rval = 'A' - 10;
   else
      rval = 'a' - 36;
   return( (char)( rval + number));
}

int create_mpc_packed_desig( char *packed_desig, const char *obj_name)
{
   int i, j, rval = 0;
   char buff[20], comet_desig = 0;

               /* Check for comet-style desigs such as 'P/1995 O1' */
               /* and such.  Leading character can be P, C, X, D, or A. */
   if( strchr( "PCXDA", *obj_name) && obj_name[1] == '/')
      {
      comet_desig = *obj_name;
      obj_name += 2;
      }

               /* Create a version of the name with all spaces removed: */
   for( i = j = 0; obj_name[i] && j < 19; i++)
      if( obj_name[i] != ' ')
         buff[j++] = obj_name[i];
   buff[j] = '\0';

   for( i = 0; isdigit( buff[i]); i++)
      ;
               /* If the name starts with four digits followed by an */
               /* uppercase letter,  it's a provisional designation: */
   if( buff[i] && i == 4 && isupper( buff[4]))
      {
      const int year = atoi( buff);
      int sub_designator;

      memset( packed_desig, ' ', 12);
      for( i = 0; i < 4; i++)
         {
         const char *surveys[4] = { "P-L", "T-1", "T-2", "T-3" };

         if( !strcmp( buff + 4, surveys[i]))
            {
            const char *surveys_packed[4] = {
                     "PLS", "T1S", "T2S", "T3S" };

            memcpy( packed_desig + 8, buff, 4);
            memcpy( packed_desig + 5, surveys_packed[i], 3);
            return( rval);
            }
         }

      sprintf( packed_desig + 5, "%c%02d",
                  'A' - 10 + year / 100, year % 100);
      packed_desig[6] = buff[2];    /* decade */
      packed_desig[7] = buff[3];    /* year */

      packed_desig[8] = (char)toupper( buff[4]);    /* prelim desigs are */
      i = 5;                                        /* _very_ scrambled  */
      if( isupper( buff[i]))                        /* when packed:      */
         {
         packed_desig[11] = buff[i];
         i++;
         }
      else
         packed_desig[11] = '0';

      sub_designator = atoi( buff + i);
      packed_desig[10] = (char)( '0' + sub_designator % 10);
      if( sub_designator < 100)
         packed_desig[9] = (char)( '0' + sub_designator / 10);
      else if( sub_designator < 360)
         packed_desig[9] = (char)( 'A' + sub_designator / 10 - 10);
      else
         packed_desig[9] = (char)( 'a' + sub_designator / 10 - 36);
      if( comet_desig)
         {
         packed_desig[4] = comet_desig;
         while( isdigit( buff[i]))
            i++;
         if( buff[i] >= 'a' && buff[i] <= 'z')
            packed_desig[11] = buff[i];
         }
      }
   else if( !buff[i])       /* simple numbered asteroid or comet */
      {
      const int number = atoi( buff);

      if( comet_desig)
         sprintf( packed_desig, "%04d%c       ", number, comet_desig);
      else
         sprintf( packed_desig, "%c%04d       ", mutant_hex( number / 10000),
               number % 10000);
      }
   else                 /* strange ID that isn't decipherable.  For this, */
      {                 /* we just copy the first twelve non-space bytes, */
      while( j < 12)                              /* padding with spaces. */
         buff[j++] = ' ';
      buff[12] = '\0';
      strcpy( packed_desig, buff);
      rval = -1;
      }
   return( rval);
}

int main( const int argc, const char **argv)
{
   if( argc >= 2)
      {
      char obuff[80];
      const int rval = create_mpc_packed_desig( obuff, argv[1]);

      obuff[12] = '\0';
      printf( "'%s'\n", obuff);
      if( rval)
         printf( "Error code %d\n", rval);
      }
   return( 0);
}
