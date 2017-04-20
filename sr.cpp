/* sr.cpp: computes range limits for statistical ranging

Copyright (C) 2010, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */

/*  First off,  be advised that after deriving all of this,  I found
that much of it is explained a bit in Andrea Milani and Giovanni
Gronchi's book _The Theory of Orbit Determination_:

https://www.projectpluto.com/books.htm#milani

   ...the relevant part appeared originally in:

A. Milani, G.-F. Gronchi, M. de' Michieli Vitturi \& Z. Kne\v
zevi\'c : {\it Orbit Determination with Very Short Arcs. I Admissible
Regions}, Celestial Mechanics, 90, 57--85, 2004.

   You may like their explanation more than mine,  and there are other
good reasons to get the _Theory_ book.

   One problem we'll run into with SR is figuring out what the range of
possible distances to the first observation (r1) might be.  For example,
you'll usually find that if you assume r1=100 AU,  any orbit that gets
you to the second observation has to be absurdly hyperbolic.  And a
nearby orbit is (in theory) always possible for a very short arc.
It would be nice to know what's possible in between.  In other words:
For what range(s) of r1 are elliptical orbits possible?

   Rephrased again:  for what range(s) of r1 is it possible to get to
the second observation without exceeding the solar escape speed?

   With the latter phrasing,  it becomes possible to solve the problem.
(It seems to be common to _not_ solve it,  and just ask the user to
enter the range of r1 to be searched.  But I'd like to avoid that.
The user will either select too great a range,  in which case lots of
orbits get computed that can't possibly be elliptical;  or selects too
small a range,  in which case perfectly valid orbits will be overlooked.)

   At a distance r from the sun,  the square of the escape speed is
2GM/r.  With our distance from the first observation r1,  we can find
the distance to the ray defined by the second observation.  If that
distance is greater than the escape speed multiplied by (t2-t1),  then
there's no way to get to that second observation in an elliptical orbit.

   So.  Assume observations made with the observer at vectors q1 and q2,
with the object observed at unit vectors p1, p2 where

x-component of p1 = cos( ra1) * cos( dec1)
y-component of p1 = sin( ra1) * cos( dec1)
z-component of p1 = sin( dec1)

   and similarly for p2.  Then the actual location of the object at
time t1 would be

loc1 = q1 + r1 * p1

   and the escape velocity v is given by

v^2 = 2GM / |loc1|

   The shortest vector connecting loc1 to the ray defined by the second
observation is

dist = loc1 - q2 - p2( p2.(loc1 - q2))
     = q1 + r1 * p1 - q2 - p2( p2.(q1 + r1 * p1 - q2))

   Define dq = q1 - q2,  then

dist = dq + r1 * p1 - p2( p2.(dq + r1 * p1))
     = dq + r1 * p1 - p2( p2.dq + r1(p2.p1))
     = dq - p2(p2.dq) + r1(p1 - p2(p2.p1))

   If we set the vectors a = dq - p2(p2.dq) and b = p1 - p2(p2.p1),  then

dist = a + r1 * b

   and the square of that distance is

dist^2 = a.a + 2(a.b)r1 + (b.b)r1^2

   So.  |dist| is the _minimum_ distance it would take to get to the second
observation,  and v * (t2-t1) is the _maximum_ distance we can go while
still staying below the heliocentric escape speed.  So acceptable values
of r1 have to satisfy

v * (t2 - t1) = v * dt < |dist|

   Squaring both sides,  we get

2GM * dt^2 / |loc1| < a.a + 2(a.b)r1 + (b.b)r1^2

2GM * dt^2 < |loc1| * (a.a + 2(a.b)r1 + (b.b)r1^2)

   Now we gotta square again,  and replace |loc1|^2 with (q1 + r1 * p2) dotted
by itself:

4 * GM^2 * dt^4 < (q1.q1 + 2(p1.q1)r1 + r1^2)(a.a + 2(a.b)r1 + (b.b)r1^2) ^ 2

   ...which is a sixth-degree polynomial in r1.  This can have zero,  two,
four,  or six real roots,  corresponding to no,  one,  two,  or three
possible ranges for r1.  "Real" data ought to have at least one range
that encloses r1 = 0.  Fortunately,  I already have a good function to find
real roots of a real polynomial,  something written for Gauss' method.

   In theory,  at least,  one or more ranges could be entirely negative.
(I've never seen this happen with real-world data,  but it's trivial to
construct such a case:  if you reverse p1 and p2 -- as if the observations
were made 180 degrees from where they really were -- then the sign of
all the roots will be reversed.)  So at the end of find_sr_ranges( ),
such ranges are dropped.  Also,  if the lower limit of the range is
less than zero (always true for at least one range),  then that lower
limit is bumped up to zero.   */

int find_real_polynomial_roots( const double *poly, int poly_degree,
                                double *real_roots);        /* roots.cpp */

double dot_product( const double *v1, const double *v2)
{
   return( v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
}

/* multiply_poly() takes two polynomials of degree degree_a and degree_b,
with coefficients in a[] and b[],  and computes the product polynomial
of degree degree_a + degree_b,  putting the coefficients in poly[].  This
is convenient because the sixth-degree polynomial described above is the
product of three quadratics.  Two of them are multiplied to make a quartic
equation,  then that's multiplied by the remaining quadratic to make the
coefficients of the sixth-degree poly.   */

static void multiply_poly( double *poly, const double *a, const int degree_a,
                                         const double *b, const int degree_b)
{
   int i, j;

   for( i = 0; i <= degree_a + degree_b; i++)
      poly[i] = 0;
   for( i = 0; i <= degree_a; i++)
      for( j = 0; j <= degree_b; j++)
         poly[i + j] += a[i] * b[j];
}

int find_sr_ranges( double *ranges, const double *q1, const double *p1,
                                    const double *q2, const double *p2,
                                    const double gm, const double dt)
{
   int i, n_roots;
   const double q1_dot_q1 = dot_product( q1, q1);
   const double p1_dot_q1 = dot_product( p1, q1);
   const double p1_dot_p2 = dot_product( p1, p2);
   double a_dot_a = 0., b_dot_b = 0., a_dot_b = 0., dq[3];
   const double gm_times_dt_squared = gm * dt * dt;
   double p2_dot_dq;
   double quad[3], quart[5], poly[7];

   for( i = 0; i < 3; i++)
      dq[i] = q1[i] - q2[i];
   p2_dot_dq = dot_product( p2, dq);
   for( i = 0; i < 3; i++)
      {
      double a, b;

      a = dq[i] - p2[i] * p2_dot_dq;
      b = p1[i] - p2[i] * p1_dot_p2;
      a_dot_a += a * a;
      a_dot_b += a * b;
      b_dot_b += b * b;
      }
   quad[0] = a_dot_a;
   quad[1] = 2 * a_dot_b;
   quad[2] = b_dot_b;
   multiply_poly( quart, quad, 2, quad, 2);
   quad[0] = q1_dot_q1;
   quad[1] = 2. * p1_dot_q1;
   quad[2] = 1.;
   multiply_poly( poly, quad, 2, quart, 4);
   poly[0] -= 4. * gm_times_dt_squared * gm_times_dt_squared;
   n_roots = find_real_polynomial_roots( poly, 6, ranges);
   for( i = 0; i < n_roots; i += 2)
      if( ranges[i + 1] < 0.)     /* entire range is less than zero; */
         {                       /* delete this range */
         int j;

         for( j = i; j + 2 < n_roots; j++)
            ranges[j] = ranges[j + 2];
         i -= 2;
         n_roots -= 2;
         }
      else if( ranges[i] < 0.)    /* part of the range is less than zero; */
         ranges[i] = 0.;          /* correct this */
   return( n_roots / 2);
}

#ifdef TEST_MAIN

/* Use

gcc -DTEST_MAIN -o sr sr.cpp roots.cpp       */

#include <stdio.h>

#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)

int main( const int argc, const char **argv)
{
   const double q1[3] = { -.7256942252134, -0.7001830768909, .0000542781309 };
   const double q2[3] = { -.1516775381631, -1.0041053311688, .0000563616383 };
   const double p1[3] = { -.1609547406418, -0.8467946127346, .5069836834038 };
   const double p2[3] = { -.3530129352031, -0.8511896328040, .3884045269079 };
   const double dt = 38.9386800001375;
   double ranges[101];
   const double poly[7] = { -0.6269244422348,    0.06095061759976,
                             0.005371809737067, -0.07844593351442,
             0.061886725027, -0.01861368496397, 0.002531828864788 };
   int i, n_roots;

   setvbuf( stdout, NULL, _IONBF, 0);
   n_roots = find_real_polynomial_roots( poly, 6, ranges);
   printf( "%d roots found\n", n_roots);
   for( i = 0; i < n_roots; i++)
      printf( "Root %d: %.13lg\n", i, ranges[i]);
   printf( "Before:\n");
   printf( "q1: %.13lf %.13lf %.13lf\n", q1[0], q1[1], q1[2]);
   printf( "q2: %.13lf %.13lf %.13lf\n", q2[0], q2[1], q2[2]);
   printf( "p1: %.13lf %.13lf %.13lf\n", p1[0], p1[1], p1[2]);
   printf( "p2: %.13lf %.13lf %.13lf\n", p2[0], p2[1], p2[2]);
   printf( "Calling sr routine\n");
   n_roots = find_sr_ranges( ranges, q1, p1, q2, p2, SOLAR_GM, dt);
   for( i = 0; i < n_roots; i++)
      printf( "Root %d: %.13lg\n", i, ranges[i]);
   return( 0);
}
#endif
