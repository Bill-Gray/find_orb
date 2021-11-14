/* roots.cpp: finds all real roots of a real-coeff polynomial

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

/* Code to find all real roots of a polynomial with real coefficients.
Used in 'gauss.cpp' to get solutions to the equation of Lagrange:

   x^8 + a * x^6 + b * x ^ 3 + c = 0

   and in 'sr.cpp' to solve a particular sixth-degree polynomial used
in figuring out how far away an object might be.  But it could be used
for any polynomial with real coefficients.  It's awkward code,  but
has the virtue of finding _all_ real roots, including double roots and
closely-paired roots that might be neglected by other methods.  (Such
close pairs of solutions are common in the equation of Lagrange.)

   Its other virtue is that finding all real roots is guaranteed,  at
least in mathematical theory and in my experience (I can't rule out
the possibility that some sort of roundoff issue,  for example,  might
come out and savage you.)  For my purposes,  this was very important;
I needed "it always works,  no matter what" more than "it's very fast,
but sometimes fails."

   It's currently limited to polynomials of degree less than ten.
Redefine MAX_POLY_ORDER if this is a problem.

   It works via a recursive process:  it finds the real roots of the
_derivative_ of the polynomial.  (Therefore,  if you call this function
with an eighth-degree polynomial,  it will internally generate the
seventh-degree derivative polynomial,  find the real roots of that,
which will mean finding the roots of a sixth-degree polynomial, and so
on. As you'll see,  this eventually works down to a first-degree
polynomial, a.k.a. a linear equation,  which is trivially solved without
recursion.)

   Anyway,  the "real roots of the derivative of the polynomial" are also
known as the minima/maxima of the polynomial.  We can know as a fact that
all real roots of the polynomial lie between these maxima/minima.  So the
code evaluates the polynomial at each maximum/minimum;  if its value goes
through a sign change between consecutive maxima/minima, we've got a root
bracketed and can use find_poly_root_between() to nail it down precisely
with the method of Laguerre.

   One odd thing has to be done here,  though.  There are maxima/minima at
+/- infinity,  too.  So if we find,  say,  three minima at x1 <= x2 <= x3,
then there may be roots between (x1, x2) and (x2, x3), plus roots between
-inf and x1 and between x3 and +inf.  Since an infinite limit is awkward,
we find a limit on the possible magnitude of all roots, due to Cauchy
(see below).  (Plus another limit I came up with,  which turned out to be
a "rediscovery" of a limit from Donald Knuth's _The Art of Programming_.
Knuth's limit is frequently smaller than Cauchy's... but not always,  so
it pays to check them both.)  In general, if we find N minima/maxima,
then add in the two determined via Cauchy or Knuth,  there will be N+1
ranges between them to search.  If the polynomial has all real roots,
then every one of those ranges will bracket a root.

   One other small improvement has been made:  zero roots are factored out
up front.  At least for use here,  such roots are extremely common, since
the equation of Lagrange has lots of zero coefficients (five out of nine).

   Another small improvement would be to note that,  if the polynomial is
of odd degree, we can bracket and remove a root quite readily,  then work
on a polynomial of lesser degree.  If it's of even degree,  and the
constant and leading terms are of opposite sign,  we have _two_ real roots
that can be easily found and removed,  one positive and one negative.
Since we know that f(0) is the constant term and f(1) and f(-1) can be
easily determined via addition,  we may be able to bracket and remove
still more roots with little trouble.  I've started on this.  I may not
bother finishing it, because the code is plenty fast enough for my needs
at present;  the slight loss of precision inherent in removing roots this
way is more worrisome than the slight gain in speed. */

#include <math.h>       /* used for fabs( ) prototype */
#include <assert.h>

int find_real_polynomial_roots( const double *poly, int poly_degree,
                                double *real_roots);        /* roots.cpp */

static double evaluate_poly( const double *poly, int poly_degree, const double x)
{
   double rval = 0., power = 1.;

   poly_degree++;
   while( poly_degree--)
      {
      rval += (*poly++) * power;
      power *= x;
      }
   return( rval);
}

/* For the method of Laguerre,  we need the value of the polynomial for a
given x,  _and_ the values of the first and second derivatives of said
poly.  Evaluating all three at once speeds matters up slightly.

   As explained below,  calls with x=0 are reasonably frequent.  In such
cases,  we can skip the loop.
*/

static inline double evaluate_poly_and_derivs( const double *poly,
              int poly_degree, const double x, double *deriv, double *deriv2)
{
   double rval = poly[2], power = x;
   int i;

   *deriv = *deriv2 = 2 * poly[2];
   if( x)
      for( i = 3; i <= poly_degree; i++)
         {
         const double term = poly[i] * power;
         const double term2 = (double)i * term;

         *deriv2 += term2 * (double)(i - 1);
         *deriv += term2;
         rval += term;
         power *= x;
         }
   *deriv = poly[1] + (*deriv * x);
   return( poly[0] + x * (poly[1] + rval * x));
}

#define SEARCH_TOLERANCE 1.e-10

/* find_poly_root_between() looks for a zero of a polynomial between
two bracketing points.  It starts with a linear interpolation between
the initial two points,  then uses the method of Laguerre to find the
actual root.  If the Laguerre iteration fails due to a negative square
root,  we do a Newton-Raphson step (we've already computed the first
derivative anyway).  If the Laguerre or N-R step would put us outside
the bracketing points,  we reject it and do bisection instead.

   Laguerre is noted for almost always finding a root briskly.  If it
goes past 'max_iterations',  it's probably failing to converge,  and
we switch to bisecting on every other step.  A few bisections will
probably nudge us off the problem region,  resulting in speedy
convergence.  If it doesn't reinstate convergence,  it'll still
converge with worst-case linear performance.  I've never actually seen
Laguerre fail,  except when testing with a low 'max_iterations' value,
but it's theoretically possible.  (In particular,  a double root will
converge slowly,  and a triple or higher root will converge _really_
slowly.  The use of this bisection code ensures decent worst case
behavior,  even in such cases.)

   If x1 < 0 and x2 > 0,  then the bracket contains x=0,  and we can
"compute" the polynomial at that point y(0) = poly[0] trivially, and
then adjust the brackets accordingly.  The following code uses this
little trick to get a slight boost in performance.

   It's almost as simple to compute the value of the polynomial at x=1
(just add all the coefficients) and at x=-1 (add all even coeffs,  subtract
all the odd ones),  and I gave this a try.  It does help a little,  and
might be worthwhile in some situations,  but I ended up commenting it out
with #ifdef PROBABLY_NOT_WORTHWHILE.

   There are three reasons for bisecting,  enumerated below.  */

#define BISECTING_BECAUSE_OTHER_METHODS_FAILED       1
#define BISECTING_BECAUSE_STEPPED_OUT_OF_BRACKETS    2
#define BISECTING_BECAUSE_TOO_MANY_ITERATIONS        3

static double find_poly_root_between( const double *poly, const int poly_degree,
                               double x1, double y1, double x2, double y2)
{
   double delta, x;
   int iteration = 0;

   if( x1 < 0. && x2 > 0.)            /* brackets span zero;  move   */
      {
      if( y1 * poly[0] > 0.)          /* one end to x=0, y = poly[0] */
         {
         x1 = 0.;
         y1 = poly[0];
         }
      else
         {
         x2 = 0.;
         y2 = poly[0];
         }
      }
#ifdef PROBABLY_NOT_WORTHWHILE            /* see comments above */
   if( (x1 < 1. && x2 > 1.) || (x1 < -1. && x2 > -1.))
      {
      double sum_odd_coeffs = 0., sum_even_coeffs = 0.;
      double new_y, new_x;

      for( unsigned i = 0; i <= poly_degree; i += 2)
         sum_even_coeffs += poly[i];
      for( unsigned i = 1; i <= poly_degree; i += 2)
         sum_odd_coeffs += poly[i];
      if( x1 < 1. && x2 > 1.)          /* brackets span x=1:  move */
         {                             /* one end to x=1, y=poly(1) */
         new_y = sum_even_coeffs + sum_odd_coeffs;
         new_x = 1.;
         }
      else                           /* brackets span x=-1:  move */
         {                           /* one end to x=-1, y=poly(-1) */
         new_y = sum_even_coeffs - sum_odd_coeffs;
         new_x = -1.;
         }

      if( y1 * new_y > 0.)       /* move lower bound up */
         {
         x1 = new_x;
         y1 = new_y;
         }
      else                       /* move upper bound down */
         {
         x2 = new_x;
         y2 = new_y;
         }
      }
#endif            /* #ifdef PROBABLY_NOT_WORTHWHILE */
   if( y1 == y2)
      x = (x1 + x2) * .5;
   else
      x = x1 + (x2 - x1) * y1 / (y1 - y2);    /* linear interpolation */
   do
      {
      double y, deriv, deriv2;
      int use_bisection = 0;
      const int max_iterations = 9;

      iteration++;
      y = evaluate_poly_and_derivs( poly, poly_degree, x, &deriv, &deriv2);
      if( (y1 < 0. && y < 0.) || (y1 > 0. && y > 0.))
         {
         x1 = x;
         y1 = y;
         }
      else
         {
         x2 = x;
      /* y2 = y; */     /* We never actually access y2 from this point */
         }              /* forward,  so there's no need to update it   */
      if( !y)
         delta = 0.;
      else if( iteration > max_iterations && (iteration % 2))
         use_bisection = BISECTING_BECAUSE_TOO_MANY_ITERATIONS;
                   /* Laguerre has probably failed; bisect instead */
      else
         {
         const double n = poly_degree;
         const double big_g = deriv / y;
         const double big_h = big_g * big_g - deriv2 / y;
         const double discr = (n - 1) * (n * big_h - big_g * big_g);

         if( discr >= 0.)     /* take Laguerre step */
            {
            if( big_g > 0.)
               delta = -n / (big_g + sqrt( discr));
            else
               delta = -n / (big_g - sqrt( discr));
            }
         else if( deriv)     /* if Laguerre fails,  try N-R: */
            delta = -y / deriv;
         else                 /* Laguerre _and_ Newton-Raphson fail;  */
            {
            delta = 0.;
            use_bisection = BISECTING_BECAUSE_OTHER_METHODS_FAILED;
            }
         x += delta;
         if( x < x1 || x > x2)
            use_bisection = BISECTING_BECAUSE_STEPPED_OUT_OF_BRACKETS;
         }
      if( use_bisection)
         {
         x =     (x1 + x2) * .5;
         delta = (x1 - x2) * .5;
         }
      }
      while( delta > SEARCH_TOLERANCE || delta < -SEARCH_TOLERANCE);
   return( x);
}

/* Cauchy came up with a useful bound to the absolute value of the roots of
a polynomial:  for the polynomial y = a0 + a1 * x + a2 * x^2 + ... + an * x^n,
the roots will all have absolute values less than

c = 1. + max | ai / an |,  with i=0...n-1

   For comments/proof,  see http://en.wikipedia.org/wiki/Sturm_chain (search
for 'Cauchy' within the text) or:
http://fermatslasttheorem.blogspot.com/2009/02/cauchys-bound-for-real-roots.html
*/

static inline double cauchy_upper_root_bound( const double *poly,
                               int poly_degree)
{
   double max = 0., new_max;

   while( poly_degree--)
      if( (new_max = fabs( *poly++)) > max)
         max = new_max;
   return( 1. + max / fabs( *poly));
}

static inline double knuth_upper_root_bound( const double *poly,
                               int poly_degree)
{
   double rval = 0.;
   int i;

   for( i = 0; i < poly_degree; i++)
      {
      const double new_ratio =
                   log( fabs( poly[i] / poly[poly_degree])) / (double)(poly_degree - i);

      if( !i || new_ratio > rval)
         rval = new_ratio;
      }
   return( 2. * exp( rval));
}

static double best_upper_root_bound( const double *poly, int poly_degree)
{
   const double cauchy_bound = cauchy_upper_root_bound( poly, poly_degree);
   const double knuth_bound  =  knuth_upper_root_bound( poly, poly_degree);

   return( cauchy_bound < knuth_bound ? cauchy_bound : knuth_bound);
}

#define MAX_POLY_ORDER 10

int find_real_polynomial_roots( const double *poly, int poly_degree,
                                double *real_roots)
{
   int i, n_roots_found = 0;

   while( poly_degree > 1 && poly[0] == 0.)   /* get rid of zero roots quickly. */
      {                                  /* Not strictly necessary, but it */
      real_roots[n_roots_found++] = 0.;  /* speeds things up if there are  */
      poly_degree--;                    /* lots of zero roots... which there */
      poly++;                          /* are when solving Legendre eqn */
      }

#ifdef NOT_COMPLETED_YET
            /* If this is an odd-degree polynomial,  or an even-degree one */
            /* where the leading and constant coefficients have opposite   */
            /* sign,  we can bracket a root trivially and factor it out.   */
            /* There might be a slight loss of precision when doing this,  */
            /* and 'poly' would have to become non-const.                  */
            /*   At present,  it's really important that this routine work */
            /* correctly for every case and be simple to debug.  If speed  */
            /* became an issue,  this sort of bracketing might be desirable. */
            /* But for the nonce... it's best not to bother.                 */
   while( poly_degree > 1 && ((poly_root & 1) || poly[0] * poly[degree] < 0.))
      {
      double bound = best_upper_root_bound( poly, poly_degree);
      double root, ybound;

      if( poly[0] * poly[degree] > 0.)
         bound = -bound;
      ybound = evaluate_poly( poly, poly_degree, bound);
      if( bound < 0.)
         root = find_poly_root_between( poly, poly_degree,
                           bound, ybound, 0, poly[0]);
      else
         root = find_poly_root_between( poly, poly_degree,
                           0, poly[0], bound, ybound);
            /* ...then synthetic-divide by root,  then... */
      real_roots[n_roots_found++] = root;
      poly_degree--;
      poly++;
      }
#endif

   assert( poly_degree >= 1);
   if( poly_degree == 1)      /* simple linear case */
      real_roots[n_roots_found++] = -poly[0] / poly[1];
   else
      {
      double slope_poly[MAX_POLY_ORDER], minmax[MAX_POLY_ORDER];
      double x1, y1, x2, y2;
      int n_minmax, i;
      const double bound = best_upper_root_bound( poly, poly_degree);

      for( i = 0; i < poly_degree; i++)    /* find derivative of poly: */
         slope_poly[i] = (double)( i + 1) * poly[i + 1];
      n_minmax = find_real_polynomial_roots( slope_poly, poly_degree - 1,
                                                minmax);

      x1 = -bound;
      y1 = evaluate_poly( poly, poly_degree, x1);
      for( i = -1; i < n_minmax; i++)
         {
         if( i != n_minmax - 1)
            x2 = minmax[i + 1];
         else
            x2 = bound;
         y2 = evaluate_poly( poly, poly_degree, x2);
/*       printf( "Range %lf %lf (%lf %lf)\n", x1, x2, y1, y2);    */
               /* Make sure there is root searching to do (i.e.,  there is */
               /* a range to search and a sign change within that range): */
         if( y1 * y2 <= 0.)
            {
            real_roots[n_roots_found++] = find_poly_root_between( poly,
                                  poly_degree, x1, y1, x2, y2);
/*          printf( "Root = %lf\n", real_roots[n_roots_found - 1]);  */
            }
         if( i != n_minmax - 1)
            {
            x1 = minmax[i + 1];
            y1 = y2;
            }
         }
      }
                  /* If zero roots were found,  then our list of roots may */
                  /* not be in numerical order.  The following ugly sort   */
                  /* (well-suited to almost-sorted lists) will fix that.   */
   for( i = 0; i < n_roots_found - 1; i++)
      if( real_roots[i + 1] < real_roots[i])
         {
         const double tval = real_roots[i + 1];

         real_roots[i + 1] = real_roots[i];
         real_roots[i] = tval;
         if( i)
            i -= 2;
         }
   return( n_roots_found);
}

#ifdef TEST_CODE

/* If compiled with TEST_CODE defined,  you get the following little
test routine which can be run with coefficients on the command line.  The
roots are then printed out on the console. */

#include <stdio.h>
#include <stdlib.h>

int main( const int argc, const char **argv)
{
   int degree = argc - 2, n_roots;
   double poly[10], roots[10];

   for( int i = 0; i <= degree; i++)
      poly[i] = atof( argv[i + 1]);
   n_roots = find_real_polynomial_roots( poly, degree, roots);
   for( int i = 0; i < n_roots; i++)
      printf( "%lf   ", roots[i]);
   return( 0);
}
#endif
