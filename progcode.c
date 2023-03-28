/* Utility to read the 'ProgramCodes.txt' file from

https://minorplanetcenter.net/iau/lists/ProgramCodes.txt

and produce a file 'progcode.txt' suitable for Find_Orb to extract
observer info from.  Note that some editing of the input (or output)
is required; some observers appear twice,  with modifications due to
diacritical marks.  When doing updates,  I compare to the previous
version to make sure such problems don't creep back in.

   Compile with :

gcc -Wall -pedantic -Wextra -O3 progcode.c -o progcode

   Extension of the table provided at the above URL.  Some inference
is involved here,  and I'm not totally certain of the extended part.

10 !       20 +      30 `     40 ?      50 I      60 S
11 "       21 ,      31 {     41 @      51 J      61 T
12 #       22 -      32 |     42 A      52 K      62 U
13 $       23 .      33 }     43 B      53 L      63 V
14 %       24 /      34 ~     44 C      54 M
15 &       25 [      35 :     45 D      55 N
16 '       26 \      36 ;     46 E      56 O
17 (       27 ]      37 <     47 F      57 P
18 )       28 ^      38 =     48 G      58 Q
19 *       29 _      39 >     49 H      59 R
*/

#include <stdio.h>
#include <string.h>

/* Returns 0 if the line was good,  or -1 if it should be expunged.
The original 'ProgramCodes.txt' has several blunders in it. */

static int correct_line( char *buff)
{
   size_t i = strlen( buff);
   char *tptr;
   const char *removers[6] = {
               "568 1 D, C. Jewitt",
               "658 4 JJ Kavelaars",
               "G40 9 B. Luetkenhoener",
               "Q62 } S. K\\\"urti",
               "U69 # B. Luetkenhoener",
               NULL };

   while( i && buff[i - 1] <= ' ')
      i--;
   buff[i] = '\n';      /* remove trailing spaces */
   buff[i + 1] = '\0';
                  /* ad hoc fixes to blunders in original text : */
   for( i = 0; removers[i]; i++)
      if( strstr( buff, removers[i]))
         return( -1);
   tptr = strstr( buff, "JJ Kavelaars");
   if( tptr)
      strcpy( tptr, "J. J. Kavelaars\n");
   if( strstr( buff, "D, C. Jewitt"))
      return( -1);
   return( 0);
}

int main( const int argc, const char **argv)
{
   const char *filename = (argc == 1 ? "ProgramCodes.txt" : argv[1]);
   FILE *ifile = fopen( filename, "rb"), *ofile;
   char buff[200], prev_code[5];

   if( !ifile)
      {
      fprintf( stderr, "Couldn't open %s\n", filename);
      return( -1);
      }
   memset( prev_code, 0, sizeof( prev_code));
   ofile = fopen( "progcode.txt", "wb");
   while( fgets( buff, sizeof( buff), ifile))
      if( buff[0] > ' ' && buff[1] > ' ' && buff[2] > ' '
                        && buff[3] == ' ' && buff[4] > ' '
                        && buff[5] == ' ' && buff[6] > ' '
                        && !correct_line( buff))
         {
         if( memcmp( prev_code, buff, 5))
            {
            fprintf( ofile, "\nCOD %.5s\n", buff);
            memcpy( prev_code, buff, 5);
            }
         fprintf( ofile, "OBS %s", buff + 6);
         }
   fclose( ifile);
   fclose( ofile);
   return( 0);
}
