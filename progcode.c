/* Utility to read the 'ProgramCodes.txt' file from

https://minorplanetcenter.net/iau/lists/ProgramCodes.txt

   (should switch over to the JSON form at

https://minorplanetcenter.net/static/downloadable-files/program_codes.json

and produce a file 'progcode.txt' suitable for Find_Orb to extract
observer info from.  Note that some editing of the input (or output)
is required; some observers appear twice,  with modifications due to
diacritical marks.  That's because MPC uses the file at the above URL
to determine the program code associated with a given name,  and
people may submit as (say) "G. Borisov",  then as "G. V. Borisov".
I'm using the file to find the (canonical) name for a given program
code,  and therefore drop the duplicates.

When doing updates,  I compare to the previous version to make sure
such problems don't creep back in.

   Compile with :

gcc -Wall -pedantic -Wextra -O3 progcode.c -o progcode

   Extension of the table provided at the above URL.  The top line
goes one character further to include the pound (currency) symbol,
but I'd expect a non-7-bit-ASCII character to raise hell.

10 !    20 +    30 `    40 ?    50 I    60 S    70 c    80 m    90 w
11 "    21 ,    31 {    41 @    51 J    61 T    71 d    81 n    91 x
12 #    22 -    32 |    42 A    52 K    62 U    72 e    82 o    92 y
13 $    23 .    33 }    43 B    53 L    63 V    73 f    83 p    93 z
14 %    24 /    34 ~    44 C    54 M    64 W    74 g    84 q
15 &    25 [    35 :    45 D    55 N    65 X    75 h    85 r
16 '    26 \    36 ;    46 E    56 O    66 Y    76 i    86 s
17 (    27 ]    37 <    47 F    57 P    67 Z    77 j    87 t
18 )    28 ^    38 =    48 G    58 Q    68 a    78 k    88 u
19 *    29 _    39 >    49 H    59 R    69 b    79 l    89 v
*/

#include <stdio.h>
#include <string.h>

/* Returns 0 if the line was good,  or -1 if it should be expunged.
The original 'ProgramCodes.txt' has several blunders in it. */

static int correct_line( char *buff)
{
   size_t i = strlen( buff);
   const char *removers[] = {
               "095 2 V. Rumyantsev",
               "095 4 G. Borisov",
               "266 3 J. Kavelaars",
               "267 3 J. Kavelaars",
               "568 1 D. Jewitt",
               "658 4 JJ Kavelaars",
               "703 1 E. J. Christensen",
               "G40 9 B. Luetkenhoener",
               "G96 1 E. J. Christensen",
               "Q62 } S. K\\\"urti",
               "U69 # B. Luetkenhoener",
               "W84 \" J. Pena Z.",
               NULL };
   const char *replacers[] = {
            "JJ Kavelaars\n",     "J. J. Kavelaars\n",
#ifdef LATEXIZE_NAMES
            "F. Ocana\n",         "F. Oca\\~na\n",
            "B. Lutkenhoner\n",   "B. L\\\"utkenh\\\"oner\n",
            "S. Kurti\n",         "S. K\\\"urti",
            "M. Juric\n",         "M. Juri\\'c",
#endif
            "D, C. Jewitt\n",     "D. C. Jewitt\n",
            };

   while( i && buff[i - 1] <= ' ')
      i--;
   buff[i] = '\n';      /* remove trailing spaces */
   buff[i + 1] = '\0';
                  /* ad hoc fixes to blunders in original text : */
   for( i = 0; removers[i]; i++)
      if( strstr( buff, removers[i]))
         return( -1);
   for( i = 0; replacers[i]; i += 2)
      if( !strcmp( buff + 6, replacers[i]))
         strcpy( buff + 6, replacers[i + 1]);
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
