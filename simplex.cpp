/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

#include <string.h>

/* Slightly modified Nelder-Mead simplex method,  using parameters from

http://www.webpages.uidaho.edu/~fuchang/res/ANMS.pdf

   Note that this boils down to plain old Nelder-Mead for two-dimensional
(n=2) problems.  The paper argues that their modifications have some
theoretical and empirical evidence for improving the speed of convergence
(and,  in some cases,  make convergence happen at all.)  Set alpha = 1.,
beta = 2., gamma = delta = 0.5 if you want the original N-M parameters.
*/

void init_simplex( double **vects, double *fvals,
         double (*f)( void *context, const double *vect),
               void *context, const int n);        /* simplex.c */
int simplex_step( double **vects, double *fvals,
         double (*f)( void *context, const double *vect),
               void *context, const int n);        /* simplex.c */

static void sort_simplex( double *fvals, double **vects, const int n)
{
   int i, j;

   for( i = 1; i < n; i++)
      {
      const double tval = fvals[i];
      double *tptr = vects[i];

      for( j = i; j && fvals[j - 1] > tval; j--)
         {
         fvals[j] = fvals[j - 1];
         vects[j] = vects[j - 1];
         }
      fvals[j] = tval;
      vects[j] = tptr;
      }
}

void init_simplex( double **vects, double *fvals,
         double (*f)( void *context, const double *vect),
               void *context, const int n)
{
   int i;

   for( i = 0; i <= n; i++)
      fvals[i] = f( context, vects[i]);
}

#define MAX_DIM 20

static inline void vector_adjust( double *ovect, const double *a,
            double *b, const double amount, const int n)
{
   int i;

   for( i = 0; i < n; i++)
      ovect[i] = a[i] + amount * (b[i] - a[i]);
}

/* This creates a reflected,  expanded,  or 1-D contracted point.  If
the objective function is lower at that point,  we accept that step. */

static inline void try_improvement( double *fval, double *curr,
            const double *centroid, const double factor,
            double (*f)( void *context, const double *vect),
            void *context, const int n)
{
   double new_point[MAX_DIM], new_fval;

   vector_adjust( new_point, centroid, curr, factor, n);
   new_fval = f( context, new_point);
   if( *fval > new_fval)    /* it's an improvement;  take the step */
      {
      *fval = new_fval;
      memcpy( curr, new_point, n * sizeof( double));
      }
}

/* Return value is the number of function evaluations made.  We try various
steps to get the "highest",  worst point to have a lower objective function
value than the next-to-highest step.  (If we can accomplish that,  we can
repeat simplex_step() and have it work with that point instead.  Otherwise,
we could get stuck in a loop where we keep playing with the same high point.)
However,  after trying two possible improvements,  we give up and do a "full
contraction" around the lowest point.

Thus,  the returned number of evaluations can be :

   rval = 1 if we reflected the high point through the centroid and
      got an improvement,  and did nothing more.
   rval = 2 if the reflection worked out so well that it was better
      than our lowest point,  so we did an "expansion" step,  resulting
      in a second evaluation of the objective function.
   rval = 2 if it _didn't_ work out so well,  and we tried a one-D
      contraction step,  and then _did_ get an improvement.
   rval = n + 2 if the contraction step didn't get fvals[n] < fvals[n-1],
      so we contracted around the lowest point.
*/

int simplex_step( double **vects, double *fvals,
               double (*f)( void *context, const double *vect),
               void *context, const int n)
{
   double cent[MAX_DIM];
   int i, j, n_evals = 1;
#ifndef ANMS    /* "accelerated" Nelder-Mead params */
   const double alpha = 1.;
   const double beta = 1. + 2. / (double)n;
   double gamma = 0.75 - 0.5 / (double)n;
   const double delta = 1. - 1. / (double)n;
#else       /* original params from Nelder-Mead */
   const double alpha = 1.;
   const double beta = 2.;
   double gamma = 0.5;
   const double delta = 0.5;
#endif

   sort_simplex( fvals, vects, n + 1);
   for( i = 0; i < n; i++)   /* find centroid of first n points */
      {
      cent[i] = 0.;
      for( j = 0; j < n; j++)
         cent[i] += vects[j][i];
      cent[i] /= (double)n;
      }
   try_improvement( fvals + n, vects[n], cent, -alpha, f, context, n);

   if( fvals[n] < fvals[0])   /* best value yet;  try extrapolate factor of 2 */
      {
//    printf( "Expanding\n");
      try_improvement( fvals + n, vects[n], cent, beta, f, context, n);
      n_evals = 2;
      }
   else if( fvals[n] >= fvals[n - 1])
      {                       /* try a one-dimensional contraction */
//    printf( "1-d contract\n");
      try_improvement( fvals + n, vects[n], cent, gamma, f, context, n);
      n_evals = 2;
      }
   if( fvals[n - 1] < fvals[n])  /* Haven't been able to eliminate the high */
      {                         /* point,  so contract around lowest point. */
//    printf( "Full contraction\n");
      for( i = 1; i <= n; i++)
         {
         vector_adjust( vects[i], vects[0], vects[i], delta, n);
         fvals[i] = f( context, vects[i]);
         }
      n_evals = n + 2;
      }
   sort_simplex( fvals, vects, n + 1);
   return( n_evals);
}
