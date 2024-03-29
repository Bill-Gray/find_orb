#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"

double geo_potential_in_au( const double x, const double y, const double z,
                 double *derivs, const int n_terms);    /* geo_pot.c */
extern int _starting_term;

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
/* Compile with
gcc -Wall -Wextra -DTEST_MAIN -o geo_test geo_test.cpp geo_pot.o -lm
for test code to show potential and its derivatives (i.e.,  acceleration)
for a given lat/lon/radius from the geocenter,  and the geoid height. */

static double jn_potential( const double x, const double y, const double z,
                         const double j2,  const double j3, const double j4)
{
   const double r2 = x * x + y * y + z * z;
   const double r = sqrt( r2);
   const double mu = z / r;
   const double mu2 = mu * mu;
   const double p3 = mu * (2.5 * mu2 - 1.5);
   const double p4 = (35. * mu2 * mu2 - 30. * mu2 + 3.) / 8.;
   const double p2 = 1.5 * mu2 - .5;     /* Danby, p. 115 */

   return( (j2 * p2 + j3 * p3 / r + j4 * p4 / r2) / (r * r2));
// return( (          j3 * p3 / r + j4 * p4 / r2) / (r * r2));
}

#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_R (EARTH_MAJOR_AXIS / AU_IN_METERS)
#define EARTH_R2 (EARTH_R * EARTH_R)
#define EARTH_R3 (EARTH_R * EARTH_R2)
#define J2_IN_EARTH_UNITS (1.0826355254e-3)
#define J3_IN_EARTH_UNITS (-2.5325505977e-6)
#define J4_IN_EARTH_UNITS (-1.62870177e-6)
#define EARTH_J2 (J2_IN_EARTH_UNITS * EARTH_R2)
#define EARTH_J3 (J3_IN_EARTH_UNITS * EARTH_R3)
#define EARTH_J4 (J4_IN_EARTH_UNITS * EARTH_R3 * EARTH_R)

int main( const int argc, const char **argv)
{
   const double lat = atof( argv[1]) * PI / 180.;
   const double lon = atof( argv[2]) * PI / 180.;
   int n_terms = 99999;
   size_t i;
   double oval;
   long double loc[3];
   double derivs[3], dpot;
   const double earth_gm_mks = 0.3986004415E+15;   /* GGM03 value, in m^3/s^2 */
   const double earth_gm_aud = earth_gm_mks * seconds_per_day * seconds_per_day
               / (AU_IN_METERS * AU_IN_METERS * AU_IN_METERS);
   const double omega = 7292115e-11;           /* in radians/second */
   const double earth_minor_axis = 6356755.;    /* polar radius */
   const double centrifugal_pot_at_equator =
         EARTH_MAJOR_AXIS * EARTH_MAJOR_AXIS * EARTH_MAJOR_AXIS * omega * omega / earth_gm_mks;
   const double r = atof( argv[3]);
   const double delta = r * .0001 * EARTH_MAJOR_AXIS / AU_IN_METERS;

   assert( argc > 3);
   assert( r > 0.);
   if( argc > 4)
      n_terms = atoi( argv[4]);
   loc[0] = (long double)( cos( lon) * cos( lat));
   loc[1] = (long double)( sin( lon) * cos( lat));
   loc[2] = (long double)  sin( lat);
   for( i = 0; i < 3; i++)       /* cvt distances to AU */
      loc[i] *= EARTH_R;
   _starting_term = 0;
   oval = (double)geo_potential_in_au( loc[0], loc[1],
                  loc[2] * earth_minor_axis / EARTH_MAJOR_AXIS, derivs, n_terms);
   oval = 1 + oval * EARTH_R / earth_gm_aud;
   oval -= centrifugal_pot_at_equator * cos( lat) * cos( lat) / 2.;
   printf( "%e  (geoid height %f meters)\n", oval,
                 -14497.68 - oval * EARTH_MAJOR_AXIS);
                     /* cvt radius from units of earth radii to AU: */
   printf( "  %e\n", jn_potential( loc[0] * AU_IN_METERS / EARTH_MAJOR_AXIS,
                                   loc[1] * AU_IN_METERS / EARTH_MAJOR_AXIS,
                                   loc[2] * AU_IN_METERS / EARTH_MAJOR_AXIS,
                                   EARTH_J2, EARTH_J3, EARTH_J4));
   for( i = 0; i < 3; i++)
      loc[i] *= r;
   printf( "  %e\n",
              geo_potential_in_au( loc[0], loc[1], loc[2], derivs, n_terms));
   printf( "  Derivs: %.10e %.10e %.10e (analytical)\n", derivs[0], derivs[1], derivs[2]);
   for( size_t i = 0; i < 3; i++)
      derivs[i] *= AU_IN_METERS / ((double)seconds_per_day * (double)seconds_per_day);
   printf( "  Derivs: %.10e %.10e %.10e (analytical, m/s^2)\n", derivs[0], derivs[1], derivs[2]);

   dpot = geo_potential_in_au( loc[0] + delta, loc[1], loc[2], NULL, n_terms)
        - geo_potential_in_au( loc[0] - delta, loc[1], loc[2], NULL, n_terms);
   derivs[0] = .5 * dpot / delta;
   dpot = geo_potential_in_au( loc[0], loc[1] + delta, loc[2], NULL, n_terms)
        - geo_potential_in_au( loc[0], loc[1] - delta, loc[2], NULL, n_terms);
   derivs[1] = .5 * dpot / delta;
   dpot = geo_potential_in_au( loc[0], loc[1], loc[2] + delta, NULL, n_terms)
        - geo_potential_in_au( loc[0], loc[1], loc[2] - delta, NULL, n_terms);
   derivs[2] = .5 * dpot / delta;
   printf( "  Derivs: %.10e %.10e %.10e (numerical, AU/day^2)\n", derivs[0], derivs[1], derivs[2]);
   printf( "Lat %c%f Long %c%f\n",
         (lat > 0. ? 'N' : 'S'), fabs( lat) * 180. / PI,
         (lon > 0. ? 'E' : 'W'), fabs( lon) * 180. / PI);
   return( 0);
}

/*
   Comparisons made to

http://geographiclib.sourceforge.net/cgi-bin/GeoidEval?input=44+-69.9&option=Submit

0 0:       17.226  16.927145
45 0:      46.767  47.289086
90 0:      14.898  14.898461
0 120:     59.82   58.211502
0 100:     -6.92   -6.904056
0 80:    -102.59 -102.701770
0 60:     -62.59  -62.740079
10 125:    63.91   64.616471
44 -69.9: -26.70  -26.649304
-10 -125: -10.25  -10.678747
lat lon
         online   geo_pot.c

   Note that rearrangement of the order of computation has changed the results
by a millimeter or so...

phred@phred:~/find_orb$ ./geo_pot 0 0 1 19
  -2.085725e-05
  Derivs: 4.8973461530e-01 2.3583555353e-06 -5.7623158275e-06 (analytical)
  Derivs: 9.8142864339e+00 4.7261467768e-05 -1.1547686499e-04 (analytical, m/s^2)
  Derivs: 4.8973462022e-01 2.3583549381e-06 -5.7623145439e-06 (numerical, AU/day^2)
phred@phred:~/find_orb$ ./geo_pot 31.4 159 1 19
  -2.084797e-05
  Derivs: -3.8937807981e-01 1.4946019510e-01 2.5541973048e-01 (analytical)
  Derivs: -7.8031404906e+00 2.9951837573e+00 5.1186138726e+00 (analytical, m/s^2)
  Derivs: -3.8937808013e-01 1.4946019321e-01 2.5541972839e-01 (numerical, AU/day^2)
phred@phred:~/find_orb$ ./geo_pot 31.4 159 .8 19
  -2.606101e-05
  Derivs: -6.0807926146e-01 2.3343379352e-01 3.9961972747e-01 (analytical)
  Derivs: -1.2185914289e+01 4.6780154831e+00 8.0083832089e+00 (analytical, m/s^2)
  Derivs: -6.0807926193e-01 2.3343379057e-01 3.9961972424e-01 (numerical, AU/day^2) */
