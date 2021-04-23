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
                        && buff[5] == ' ' && buff[6] > ' ')
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
