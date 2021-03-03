/* This reads the 'program.txt' file from

https://minorplanetcenter.net/iau/lists/ProgramCodes.txt

and produces a file 'progcode.txt' suitable for Find_Orb to extract
observer info from.  Note that some editing of the input is required;
some observers appear twice. */

#include <stdio.h>
#include <string.h>

int main( void)
{
   FILE *ifile = fopen( "program.txt", "rb"), *ofile;
   char buff[200], prev_code[5];

   if( !ifile)
      {
      printf( "Couldn't open program.txt\n");
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
