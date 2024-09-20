/* lsquare.cpp: least-squares computations

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

/* 2013 October:  Revised to use long doubles for matrix inversion,
and using a better choice for partial pivoting (see comments
below above the pivot_value( ) function).  I did this thinking
it would reduce roundoff errors in matrix inversion,  resulting
in more stable results for short-arc solutions.  In practice,  it
doesn't really seem to help or hinder,  but it's still a good
idea for nearly-singular matrices.  */

// #define LSQUARE_ERROR_DEBUGGING

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "lsquare.h"

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

/*    Matrix inversion offers plenty of opportunity to lose precision.  So
      I decided to use GNU C's __float128 type where possible (seems to be
      versions 4.6 and up) and long doubles otherwise.  Plain old long
      doubles usually offer 64 bits of precision;  __float128 offers
      112 bits of precision.  Either is better than "ordinary" doubles
      with a mere 52 bits of precision.   */

#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406 && !defined( __arm__) && !defined( __aarch64_)
   #define ldouble    __float128
#else
   #define ldouble    long double
#endif

#define LSQUARE struct lsquare

LSQUARE
   {
   int n_params, n_obs;
   ldouble *wtw, *uw;
   };

void *lsquare_init( const int n_params)
{
   LSQUARE *rval = (LSQUARE *)calloc( 1, sizeof( LSQUARE));

   if( !rval)
      return( (void *)rval);
   rval->n_params = n_params;
   rval->n_obs = 0;
   rval->uw = (ldouble *)calloc( (n_params + 1) * n_params, sizeof( ldouble));
   rval->wtw = rval->uw + n_params;
   return(( void *)rval);
}

#ifdef LSQUARE_ERROR_DEBUGGING
static inline ldouble lfabs( ldouble ival)
{
   return( ival > 0. ? ival : -ival);
}

#include <stdio.h>

static FILE *debug_file;

static void dump_matrix( const ldouble *matrix, const int xsize, const int ysize)
{
   int i;
   ldouble largest_element = 0., sum = 0.;

   for( i = 0; i < xsize * ysize; i++)
      {
      ldouble curr = lfabs( *matrix);

      if( largest_element < curr)
         largest_element = curr;
      sum += *matrix * *matrix;
      fprintf( debug_file, "%11.2e%s", (double)*matrix++, !((i + 1) % xsize) ? "\n" : "");
      }
   sum /= (ldouble)( ysize * xsize);
   fprintf( debug_file, "Largest element: %11.2e; RMS: %11.2e\n",
               (double)largest_element, sqrt( (double)sum));
}
#endif

double levenberg_marquardt_lambda = 0.;      /* damping factor */

int lsquare_add_observation( void *lsquare, const double residual,
                                  const double weight, const double *obs)
{
   LSQUARE *lsq = (LSQUARE *)lsquare;
   int i, j;
   const int n_params = lsq->n_params;

   for( i = 0; i < n_params; i++)
      {
      const ldouble w2_obs_i = (ldouble)( weight * weight * obs[i]);

      lsq->uw[i] += (ldouble)residual * w2_obs_i;
      for( j = 0; j < n_params; j++)
         lsq->wtw[i + j * n_params] += w2_obs_i * (ldouble)obs[j];
      lsq->wtw[i + i * n_params] +=
                  w2_obs_i * (ldouble)( obs[i] * levenberg_marquardt_lambda);
      }
   lsq->n_obs++;
   return( lsq->n_obs);
}

ldouble lsquare_determinant;

   /* A simple Gauss-Jordan matrix inverter,  with partial pivoting.  It
      first extends the size x size square matrix into a size-high by
      (2*size) wide one,  with the expanded space on the right side filled
      with an identity matrix.  It then zeroes out the lower left triangle
      of the original size x size matrix.  Some row-swapping is done in
      this procedure (that's the "partial pivoting" aspect;  see _Numerical
      Recipes_,  chap. 2.1 for details.)   */


/* Commonly,  in partial pivoting schemes,  you just use whichever
line has the greatest absolute value.  So,  when starting with a matrix
such as

87   32  -12
99 1e+30 -4e+29
.7   .3   -.2

   you'd pivot on the second line.  The problem is that,  while
99 is indeed the largest of the three values in the first column,
you'll get some horrendous loss of precision because of the remaining
values in the second line.

   This "standard" scheme would correspond to the pivot_value function
below returning fabs( *line).  Instead,  it returns *line^2 divided
by the sum of the squares of all values in the row (i.e.,  the square
of the vector norm for that row).  This is somewhat similar to the
'scaled pivoting' scheme described at

https://en.wikipedia.org/wiki/Pivot_element
*/

static ldouble pivot_value( const ldouble *line, const unsigned line_size)
{
   ldouble sum_squares = line[0] * line[0], square0 = sum_squares;
   unsigned i;

   for( i = 1; i < line_size; i++)
      sum_squares += line[i] * line[i];
   return( sum_squares ? square0 / sum_squares : 0.);
}

#define swap_ldoubles( a, b)   { const ldouble __temp = a;  a = b;  b = __temp; }

static ldouble *calc_inverse( const ldouble *src, const int size)
{
   ldouble *rval;
   ldouble *temp = (ldouble *)calloc( 2 * size * size, sizeof( ldouble));
   ldouble *tptr = temp;
   ldouble *tptr1, *tptr2;
   int i, j, k;
   const int dsize = 2 * size;

   lsquare_determinant = 1.;
   if( !temp)
      return( NULL);
               /* copy 'src' to 'temp',  doubling its width and putting */
               /* an identity matrix in the right half: */
   for( i = 0; i < size; i++)
      {
      for( j = 0; j < size; j++)
         *tptr++ = (ldouble)*src++;
      for( j = 0; j < size; j++)
         *tptr++ = ((i == j) ? 1. : 0.);
      }
/* src -= size * size;    Restores src to its original value.  But we */
/*         aren't actually using it anyway,  so we'll comment it out. */

   tptr1 = temp;
   for( i = 0; i < size; i++, tptr1 += dsize)
      {
      int pivot = -1;
      ldouble best_val = 0.;

      tptr = tptr1;
      for( j = i; j < size; j++, tptr += dsize)
         {
         const ldouble zval = pivot_value( tptr + i, size - i);

         if( zval > best_val)
            {
            best_val = zval;
            pivot = j;
            }
         }

      if( pivot == -1)     /* un-invertable matrix:  return NULL */
         {
         free( temp);
         return( NULL);
         }

      if( pivot != i)                  /* swap rows */
         {
         tptr2 = temp + dsize * pivot;
         for( j = i; j < dsize; j++)
            swap_ldoubles( tptr1[j], tptr2[j]);
         }

      for( j = i + 1; j < size; j++)
         {
         ldouble tval;

         tptr2 = temp + dsize * j;
         tval = tptr2[i] / tptr1[i];
         for( k = i; k < dsize; k++)
            tptr2[k] -= tptr1[k] * tval;
         }
      }
                  /* the lower left triangle is now cleared;  time to */
                  /* zero out the upper right triangle: */

   for( i = size - 1; i >= 0; i--)
      {
      tptr1 = temp + i * dsize;
      for( j = size; j < dsize; j++)
         {
         lsquare_determinant /= tptr1[i];
         tptr1[j] /= tptr1[i];
         }
      tptr2 = temp;
      for( k = 0; k < i; k++, tptr2 += dsize)
         for( j = size; j < dsize; j++)
            tptr2[j] -= tptr2[i] * tptr1[j];
      }

   rval = (ldouble *)calloc( size * size, sizeof( ldouble));
   if( rval)            /* copy the right-hand half of 'temp',  which */
      {                 /* now has the inverse we wanted              */
      tptr1 = rval;
      for( i = 0; i < size; i++)
         for( j = 0; j < size; j++)
            *tptr1++ = temp[(i * 2 + 1) * size + j];
      }
   free( temp);
   return( rval);
}

#ifdef CHOLESKY_INVERSION

   /* Cholesky decomposition,  as described in _Numerical Recipes_ 2.9.
      Cholesky decomposition appears to be widely recommended for inverting
      covariance matrices,  and I may end up going that route.  I started
      out with Gauss-Jackson,  though,  and it appears to work Just Fine.
      (Admittedly,  with tweaks to improve pivot selection and the
      calc_improved_matrix() trick given below to "polish" an initial
      inversion.)  I may take the time,  at some point,  to try to get
      Cholesky inversion to work... though I don't think it'll actually
      get me anything at this point.  At least for the nonce,  I'm leaving
      it ifdeffed out. */

static int cholesky_decomposition( ldouble *a, const int size, ldouble *diag)
{
   int i, j;

   for( i = 0; i < size; i++)
      for( j = i; j < size; j++)
         {
         ldouble sum = a[i * size + j];
         int k;

         for( k = i - 1; k >= 0; k--)
            sum -= a[i * size + k] * a[j * size + k];
         if( i == j)
            {
            if( sum < 0.)       /* not actually positive definite after all */
               return( -1);
            else
               diag[i] = sqrt( sum);
            }
         else
            a[j * size + i] = sum / diag[i];
         }
   return( 0);
}

   /* Inverting a symmetric positive-definite matrix,  again as described */
   /* in _Numerical Recipes_,  with slight modifications.  In their       */
   /* version,  the lower triangle is computed;  I added a few lines to   */
   /* get the upper.  NOTE: does not work yet.                            */

static int cholesky_inversion( ldouble *a, const int size)
{
   ldouble diag[30];
   int i, j, k, rval = cholesky_decomposition( a, size, diag);

   assert( !rval);
   if( !rval)
      {
      for( i = 0; i < size; i++)
         {
         a[i * (size + 1)] = 1. / diag[i];
         for( j = i + 1; j < size; j++)
            {
            ldouble sum = 0;

            for( k = i; k < j; k++)
               sum -= a[j * size + k] * a[k * size + i];
            a[j * size + i] = sum / diag[j];
            }
         }
      for( i = 1; i < size; i++)          /* copy lower diag to upper */
         for( j = 0; j < i; j++)
            a[j * size + i] = a[i * size + j];
      }
   return( rval);
}

static ldouble *invert_symmetric_positive_definite_matrix( const ldouble *a,
                           const int size)
{
   ldouble *rval = (ldouble *)calloc( size * size, sizeof( ldouble));
   int cholesky_rval;

   memcpy( rval, a, size * size * sizeof( ldouble));
   cholesky_rval = cholesky_inversion( rval, size);
   if( cholesky_rval)
      {
      free( rval);
      rval = NULL;
      }
   return( rval);
}
#endif

static void mult_matrices( ldouble *prod, const ldouble *a, const int awidth,
                  const int aheight, const ldouble *b, const int bwidth)
{
   int i, j;

   for( j = 0; j < aheight; j++)
      for( i = 0; i < bwidth; i++, prod++)
         {
         int k;
         const ldouble *aptr = a + j * awidth, *bptr = b + i;

         *prod = 0.;
         for( k = awidth; k; k--, bptr += bwidth)
            *prod += *aptr++ * (*bptr);
         }
}

/* calc_inverse_improved() computes a matrix inverse using the simpler
   'calc_inverse()' function,  a plain ol' Gauss-Jordan inverter (see above).
   It then uses a rather simple trick from _Numerical Recipes_,  chap. 2.5,
   "Iterative Improvement of a Solution to Linear Equations",  to "polish"
   the result.  The idea is this.  Suppose you've inverted A to create a matrix
   B,  with (inevitably) some error in it,  so that instead of AB=I,  you have
   AB=I+E,  where E is an "error matrix" with lots of small values.  In that
   case,  A(B - BE) = AB - (AB)E = I+E - (I+E)E = I+E - E - E^2 = I - E^2.
   In other words,  B-BE is a better approximation to the inverse of A.

      I added this step when I had some concerns that my least squares
   solutions weren't what they ought to be.  The problem lay elsewhere.
   But this _should_ ensure that matrix inversion is more accurate than it
   otherwise would be,  at almost no computational cost.  */

static ldouble *calc_inverse_improved( const ldouble *src, const int size)
{
#ifdef LSQUARE_ERROR_DEBUGGING
   debug_file = fopen( "lsquare.dat", "ab");
#endif

   ldouble *inverse = calc_inverse( src, size);

#ifdef LSQUARE_ERROR_DEBUGGING
   fprintf( debug_file, "Inverse:\n");
   dump_matrix( inverse, size, size);
#endif
   if( inverse)
      {
      ldouble *err_mat = (ldouble *)calloc( 2 * size * size, sizeof( ldouble));
      ldouble *b_times_delta = err_mat + size * size;
      int i;

      mult_matrices( err_mat, src, size, size, inverse, size);
      for( i = 0; i < size; i++)
         err_mat[i * (size + 1)] -= 1.;
      mult_matrices( b_times_delta, inverse, size, size, err_mat, size);
      for( i = 0; i < size * size; i++)
         inverse[i] -= b_times_delta[i];
#ifdef LSQUARE_ERROR_DEBUGGING
      fprintf( debug_file, "%d-square matrix delta (= AB - I):\n", size);
      dump_matrix( err_mat, size, size);
      fprintf( debug_file, "Error matrix after adjustment:\n");
      mult_matrices( err_mat, src, size, size, inverse, size);
      for( i = 0; i < size; i++)
         err_mat[i * (size + 1)] -= 1.;
      dump_matrix( err_mat, size, size);
      fclose( debug_file);
#endif
      free( err_mat);
      }
   return( inverse);
}

int lsquare_solve( const void *lsquare, double *result)
{
   const LSQUARE *lsq = (const LSQUARE *)lsquare;
   int i, j, n_params = lsq->n_params;
   ldouble *inverse;

   if( n_params > lsq->n_obs)       /* not enough observations yet */
      return( -1);

// inverse = invert_symmetric_positive_definite_matrix( lsq->wtw, n_params);
   inverse = calc_inverse_improved( lsq->wtw, n_params);
   if( !inverse)
      return( -2);            /* couldn't invert matrix */

   for( i = 0; i < n_params; i++)
      {
      ldouble temp_result = 0;

      for( j = 0; j < n_params; j++)
         temp_result += inverse[i + j * n_params] * lsq->uw[j];
      result[i] = (double)temp_result;
      }

   free( inverse);
   return( 0);
}

static double *convert_ldouble_matrix_to_double( const ldouble *matrix,
                              const int size)
{
   double *rval = (double *)calloc( size * size, sizeof( double));
   int i;

   for( i = 0; i < size * size; i++)
      rval[i] = (double)matrix[i];
   return( rval);
}

double *lsquare_covariance_matrix( const void *lsquare)
{
   const LSQUARE *lsq = (const LSQUARE *)lsquare;
   ldouble *lrval = NULL;

   if( lsq->n_params <= lsq->n_obs)       /* got enough observations */
      lrval = calc_inverse_improved( lsq->wtw, lsq->n_params);
//    lrval = invert_symmetric_positive_definite_matrix( lsq->wtw, lsq->n_params);
   if( lrval)
      {
      double *rval = convert_ldouble_matrix_to_double( lrval, lsq->n_params);

      free( lrval);
      return( rval);
      }
   else
      return( NULL);
}

double *lsquare_wtw_matrix( const void *lsquare)
{
   const LSQUARE *lsq = (const LSQUARE *)lsquare;

   return( convert_ldouble_matrix_to_double( lsq->wtw, lsq->n_params));
}

void lsquare_free( void *lsquare)
{
   const LSQUARE *lsq = (const LSQUARE *)lsquare;

   free( lsq->uw);
   free( lsquare);
}

#ifdef TEST_CODE
#include <stdio.h>

void main( int argc, char **argv)
{
   FILE *ifile = fopen( "imatrix", "rb");
   int size, i, j;
   double *matrix, *inv;

   fscanf( ifile, "%d", &size);
   matrix = (double *)calloc( size * size, sizeof( double));
   for( i = 0; i < size * size; i++)
      fscanf( ifile, "%lf", matrix + i);
   inv = calc_inverse_improved( matrix, size);
   for( i = 0; i < size; i++)
      {
      printf( "\n");
      for( j = 0; j < size; j++)
         printf( "%10.5lf", *inv++);
      }
   free( inv);
}
#endif
