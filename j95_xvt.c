#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "watdefs.h"
#include "date.h"

/*    First few lines from input pointing log :
"File","ExposureStart","RAHMS","DecDMS","Exposure","RollAngle"
"M51_20020526_234229_60.fit",26/5/2002 23:42:29,"13 28 20.1 ","+47 12 50",60.00,66.95
"C2000WM1_20020530_233548_60.fit",30/5/2002 23:35:48,"17 18 35.4 ","+35 09 40",60.00,-4.85
"C2000WM1_20020530_233838_10.fit",30/5/2002 23:38:38,"17 18 35.7 ","+35 09 38",10.00,-4.87
"C2001N2_20020530_234556_30.fit",30/5/2002 23:45:56,"19 18 42.9 ","+22 00 07",30.00,-4.89
"C2001N2_20020530_234845_30.fit",30/5/2002 23:48:45,"19 18 42.5 ","+22 00 04",30.00,-4.88
"C2001N2_20020530_235102_20.fit",30/5/2002 23:51:02,"19 18 42.7 ","+22 00 02",20.00,-4.88
"C2001K5_20020530_235546_30.fit",30/5/2002 23:55:46,"16 44 49.0 ","+14 32 11",30.00,-4.83
*/

const double seconds_per_day = 24. * 60. * 60.;

static double get_angle( const char *istr)
{
   double deg, min, sec;
   int is_neg = 0, n_scanned;

   if( *istr == '-')
      {
      is_neg = 1;
      istr++;
      }
   n_scanned = sscanf( istr, "%lf %lf %lf", &deg, &min, &sec);
   assert( n_scanned == 3);
   deg += min / 60. + sec / 3600.;
   if( is_neg)
      deg = -deg;
   return( deg);
}

int main( const int argc, const char **argv)
{
   FILE *ifile = fopen( argv[1], "rb");
   char buff[200], results[6][100];
   char exp_text[100];
   double roll_ang = 0.;

   assert( ifile);
   *exp_text = '\0';
   printf( "# (J95) Great Shefford logs,  processed with j95_xvt.c (q.v.)\n");
#ifdef __TIMESTAMP__
   printf( "# Source file date %s\n", __TIMESTAMP__);
#else
   printf( "# Version %s %s\n", __DATE__, __TIME__);
#endif
   printf( "# Input file '%s'\n", argv[1]);
   while( fgets( buff, sizeof( buff), ifile))
      {
      size_t i, j, k;
      double jd;
      char *time_str = results[1];

      for( i = j = 0; buff[i]; i++)    /* remove quote marks */
         if( buff[i] != '"')
            buff[j++] = buff[i];
      buff[j] = '\0';
      for( i = j = 0; i < 6; i++)
         {
         k = j;
         while( buff[k] != ',' && buff[k] >= ' ')
            k++;
         assert( k > j);
         assert( buff[k] == ',' || i == 5);
         memcpy( results[i], buff + j, k - j);
         results[i][k - j] = '\0';
         j = k + 1;        /* skip the comma */
         }
      jd = get_time_from_string( 0., time_str, FULL_CTIME_DMY, NULL);
      if( jd > 2e+6)
         {
         const double exposure_time = atof( results[4]) / seconds_per_day;
         const double ra = get_angle( results[2]) * 15.;
         const double dec = get_angle( results[3]);
         const double new_roll_ang = atof( results[5]);

         if( strcmp( exp_text, results[4]))
            {
            strcpy( exp_text, results[4]);
            printf( "# Exposure: %s s\n", exp_text);
            }
         if( fabs( roll_ang - new_roll_ang) > .5)
            {
            roll_ang = new_roll_ang;
            printf( "# Tilt: %.1f\n", roll_ang);
            }
         jd += exposure_time / 2.;
         full_ctime( time_str, jd, FULL_CTIME_YMD | FULL_CTIME_MILLISECS
                        | FULL_CTIME_MONTHS_AS_DIGITS | FULL_CTIME_LEADING_ZEROES);
         time_str[4] = time_str[7] = '-';
         time_str[10] = 'T';
         printf( "%.3f,%.3f,%s,J95,%s\n",
                              ra, dec, time_str, results[0]);
         }
      }
   fclose( ifile);
   return( 0);
}
