/* Compiles with

g++ -Wall -Wextra -pedantic -DFIND_MAX -o geo_max geo_max.cpp geo_pot.cpp

We pick points randomly over the unit sphere (or scaled sphere,  if
its radius is set on the command line) and compute the magnitude
of the acceleration provided by the terms in each degree l=1...N_TERMS-1.
This allows us to ignore small terms;  if we know,  for example,
that the l=5 term contributes a maximum acceleration of 1.793e-4 m/s^2
on the earth's surface (see table below),  and we know that that acceleration
will decrease by r^(l+1) = r^6 with distance,  we can very quickly
determine the maximum acceleration due to the l=5 terms and (if that
value is less than a tolerance value) not bother computing those terms.

At least in theory,  we really only need run this once and save the
results... except we may want to do the same thing for planets/satellites
other than the earth.  (Changes in Earth gravity models will happen,  but
will probably not affect the following very significantly.)

The bit about resetting the radius value is basically a sanity check to
make sure that if,  for example,  we go to two Earth radii,  then the l=2
term drops eightfold,  the l=3 terms drops 16-fold,  and so on.  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define N_TERMS 50
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#define is_power_of_two( X)   (!((X) & ((X) - 1)))
#define EARTH_MAJOR_AXIS 6378140.
#define AU_IN_KM 1.495978707e+8
#define AU_IN_METERS (AU_IN_KM * 1000.)
#define EARTH_R (EARTH_MAJOR_AXIS / AU_IN_METERS)

const double seconds_per_day = 24. * 60. * 60.;

double geo_potential_in_au( const double x, const double y, const double z,
                 double *derivs, const int n_terms);    /* geo_pot.c */

int main( const int argc, const char **argv)
{
   int l, iter, max_iter = (argc > 2 ? atoi( argv[2]) : INT_MAX);
   double maxima[N_TERMS + 1], r = EARTH_R;

   for( l = 0; l <= N_TERMS; l++)
      maxima[l] = 0.;
   if( argc > 1)
      r *= atof( argv[1]);
   for( iter = 0; iter < max_iter; iter++)
      {
      const double golden_ratio = 0.618033988749894848204586834365638117720309179805762862135448;
      const double theta = (double)iter * golden_ratio * 2. * PI;
      const double phi = (double)rand( ) * PI / (double)RAND_MAX - PI / 2.;
      const double x = cos( theta) * cos( phi) * r;
      const double y = sin( theta) * cos( phi) * r;
      const double z =               sin( phi) * r;
      double derivs[3];

      for( l = 0; l < N_TERMS; l++)
         {
         double magn;
         extern int _starting_term;    /* see geo_pot.cpp */

         _starting_term = l;
         geo_potential_in_au( x, y, z, derivs, l + 1);
         magn = sqrt( derivs[0] * derivs[0] + derivs[1] * derivs[1]
                                            + derivs[2] * derivs[2]);
         magn *= AU_IN_METERS / (seconds_per_day * seconds_per_day);
         if( maxima[l] < magn)
            maxima[l] = magn;
         if( iter > 8 && is_power_of_two( iter))
            {
            if( !l)
               printf( "%d iterations\n", iter);
            printf( "%d %e\n", l, maxima[l]);
            }
         }
      }
   return( 0);
}

/*    Example output after letting the above run a few hours,  slightly annotated :
67108864 iterations
0 0.000000e+00
1 0.000000e+00
2 3.182389e-02   max accel in m/s^2 at one Earth radius from this degree
3 2.565965e-04
4 1.793323e-04
5 1.735515e-04
6 1.728922e-04
7 1.640899e-04
8 1.204194e-04
9 1.270214e-04
10 1.180861e-04
11 1.066728e-04
12 6.365494e-05
13 9.537023e-05
14 6.969121e-05
15 6.828519e-05
16 7.836637e-05
17 7.223730e-05
18 8.435522e-05
19 7.640814e-05
20 7.170008e-05
21 7.496947e-05
22 7.833297e-05
23 6.913263e-05
24 6.512307e-05
25 7.156668e-05
26 6.328121e-05
27 4.923338e-05
28 6.474922e-05
29 7.520465e-05
30 6.661961e-05
31 7.355495e-05
32 7.175292e-05
33 7.857323e-05
34 8.297337e-05
35 8.565782e-05
36 7.856592e-05
37 7.883188e-05
38 6.901209e-05
39 9.012344e-05
40 7.331633e-05
41 6.597861e-05
42 8.001026e-05
43 8.139275e-05
44 7.368858e-05
45 7.691823e-05
46 8.372209e-05
47 8.640496e-05
48 8.220453e-05
49 7.372291e-05
50 7.849325e-05 */

