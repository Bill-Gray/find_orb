/* eigen.cpp: compute eigenvectors/values for REAL, SYMMETRIC matrices

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

#include <math.h>
#include <assert.h>

#define MAX_MATRIX_SIZE 40

void jacobi_eigenvalues( double *a, const int size, double *eigenvals,
                        double *eigenvects);       /* eigen.cpp */

#ifdef TEST_PROGRAM
#include <stdio.h>

static void show_matrix( const double *a, const int size)
{
   int i;

   for( i = 0; i < size * size; i++)
      printf( "%20.13lg%s", a[i], (i % size == size - 1 ? "\n" : " "));
}
#endif

/* Closely based on the eigenvector/value finding routine in
_Numerical Recipes in C_ (2nd edition),  p. 467-468,  with a few ideas
borrowed from

      http://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm

   NOTE THAT THIS ALGORITHM IS FOR SYMMETRIC MATRICES ONLY!

   I made several modifications,  all of which I'd consider to be
improvements:

   - Instead of a true two-dimensional array,  this assumes an array
a[0...size * size - 1].

   - The NR routine is based on Fortran-like 1 to N indexing,  instead
of C-like 0 to N-1.  The following uses the C convention.

   - In this function,  you can call with eigenvects == NULL if you have
no need for the eigenvectors.

   - Instead of doing various sweeps that look for "pretty big" values
and zeroing them,  this code just looks for the largest off-diagonal
value and zeroes it and repeats.  (Except if it's very small compared
to the corresponding eigenvalues;  in that case,  we skip it.)
Eventually,  all off-diagonal values have been zeroed, or are close
enough to zero,  and we're done.

   - NR points out that doing this sort of search "makes each Jacobi
rotation a process of order N^2 instead of N",  and they are quite
right about that.  If you used this for large N,  this could absorb
a significant amount of CPU effort.  But it would have to be pretty
large N,  and NR itself mentions that you probably shouldn't be using
the Jacobi algorithm for really large matrices anyway.  (There is
actually a rather slick way around this problem, described at
http://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm .  It finds
the largest off-diagonal value,  but doesn't go to N^2 time to do it.
But unless one had a really large matrix to solve,  and insisted on
using the Jacobi algorithm despite its not being well suited to large
matrices,  I wouldn't bother.)

   - In addition to making the code simpler,  this should mean that
fewer actual math operations,  with associated roundoff error,  will
be needed. (This may,  in fact,  make the code run a little
_faster_... possibly enough so to compensate for the extra time
spent searching for the largest off-diagonal value.  I've not
checked the speed,  though.  For my purposes,  the code is more than
fast enough; I'm _much_ more concerned about accuracy than speed!)

   - NR supplied a routine to sort the eigenvalues/vectors at the end,
using insertion sort.  I used selection sort,  which should provide
an admittedly _very_ minor speed improvement.  This is mostly an
aesthetic point:  if you're sorting large objects and therefore
attempting to minimize the number of swaps you do,  selection sort
is a better choice.

   - NR points out that half of the input matrix above the diagonal is
destroyed,  but the diagonal itself and the other half of the matrix
is left intact.  At the end of this routine,  use is made of that fact
to "repair" the input matrix.  (Turns out the aforementioned Wikipedia
page on the Jacobi eigenvalue algorithm describes the same idea.)  */

static inline void rotate( double *a1, double *a2, const double tau, const double sine)
{
   const double g = *a1;

   *a1 -= sine * (*a2 + tau * g);
   *a2 += sine * (g - tau * *a2);
}

void jacobi_eigenvalues( double *a, const int size, double *eigenvals,
                        double *eigenvects)
{
   double b[MAX_MATRIX_SIZE], z[MAX_MATRIX_SIZE];
   double largest_value = 1.;
   const double thresh = 1e-19;
   int i, j;

   assert( size < MAX_MATRIX_SIZE);
   if( eigenvects)
      for( i = 0; i < size * size; i++)
         eigenvects[i] = 0.;
   for( i = 0; i < size; i++)
      {
      j = i * (size + 1);
      if( eigenvects)
         eigenvects[j] = 1.;
      eigenvals[i] = b[i] = a[j];
      z[i] = 0.;
      }
#ifdef TEST_PROGRAM
   printf( "Initial eigenvects:\n");
   show_matrix( eigenvects, size);
#endif
   while( largest_value)
      {
      double tval;
      int iq = 0, ip = 0;

      largest_value = 0.;
      for( i = 0; i < size; i++)
         for( j = i + 1; j < size; j++)
            if( (tval = fabs( a[i + j * size])) > largest_value)
               {
               const double ej = fabs( eigenvals[j]);
               const double ei = fabs( eigenvals[i]);

               if( tval > (ei > ej ? ej : ei) * thresh)
                  {
                  largest_value = tval;
                  ip = i;
                  iq = j;
                  }
               }
      if( largest_value)
         {
         const int idx = ip + size * iq;
         const double theta = .5 * (eigenvals[iq] - eigenvals[ip]) / a[idx];
         const double tangent = (theta > 0. ? 1. : -1.)
                        / (fabs( theta) + sqrt( 1. + theta * theta));
         const double cosine = 1. / sqrt( 1. + tangent * tangent);
         const double sine = tangent * cosine;
         const double h = tangent * a[idx];
         const double tau = sine / (1. + cosine);

         z[ip] -= h;
         z[iq] += h;
         eigenvals[ip] -= h;
         eigenvals[iq] += h;
         a[idx] = 0.;
         for( j = 0; j < size; j++)
            if( j != ip && j != iq)
               {
               const int idx1 = (j < ip ? j + ip * size : ip + j * size);
               const int idx2 = (j < iq ? j + iq * size : iq + j * size);

               rotate( a + idx1, a + idx2, tau, sine);
               }
         if( eigenvects)
            for( j = 0; j < size; j++)
               rotate( eigenvects + j + ip * size,
                       eigenvects + j + iq * size, tau, sine);
         for( i = 0; i < size; i++)
            {
            b[i] += z[i];
            eigenvals[i] = b[i];
            z[i] = 0.;
            }
#ifdef TEST_PROGRAM
         printf( "ip = %d; iq = %d\n", ip, iq);
         printf( "New eigenvects:\n");
         show_matrix( eigenvects, size);
         printf( "New a:\n");
         show_matrix( a, size);
         for( i = 0; i < size; i++)
            printf( "   %f\n", eigenvals[i]);
#endif
         }
      }
         /* Selection-sort the eigenvalues & vectors.  In theory, this  */
         /* is a slow and stupid way of doing things.  But the matrices */
         /* will usually not be big;  if they are,  the behavior will   */
         /* be outweighed by the actual eigenvect/value finding;  and   */
         /* for this sort,  comparisons are cheap and actual moving may */
         /* not be (even more so for larger matrices).  So this is      */
         /* actually a case where selection sort makes sense.           */
   for( i = 0; i < size; i++)
      {
      int best_idx = i;

      for( j = i + 1; j < size; j++)
         if( eigenvals[j] < eigenvals[best_idx])
            best_idx = j;
      if( best_idx != i)
         {
         double tval = eigenvals[best_idx];

         eigenvals[best_idx] = eigenvals[i];
         eigenvals[i] = tval;
         if( eigenvects)
            for( j = 0; j < size; j++)
               {
               tval = eigenvects[j + i * size];
               eigenvects[j + i * size] = eigenvects[j + best_idx * size];
               eigenvects[j + best_idx * size] = tval;
               }
         }
      }
            /* After all this,  half the input matrix has been scrambled.   */
            /* But it can be restored by copying from the unscrambled half. */
            /* Thus,  we leave with the input matrix a unaltered.           */
   for( i = 0; i < size; i++)
      for( j = i + 1; j < size; j++)
         a[i + j * size] = a[j + i * size];
#ifdef TEST_PROGRAM
   printf( "Final eigenvects,  after sorting:\n");
   show_matrix( eigenvects, size);
   for( i = 0; i < size; i++)
      printf( "   %f\n", eigenvals[i]);
   printf( "Input matrix restored:\n");
   show_matrix( a, size);
#endif
}

#ifdef TEST_PROGRAM

      /* Example case from

      http://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm

   Output ought to be:

Final eigenvects,  after sorting:
     0.79261      0.45192      0.32242      0.25216
    -0.58208       0.3705      0.50958      0.51405
    -0.17919      0.74192     -0.10023     -0.63828
    0.029193     -0.32871      0.79141     -0.51455
Final eigenvals,  after sorting:
   0.166643
   1.478055
   37.101491
   2585.253811                  */

int main( const int argc, const char **argv)
{
// double test[9] = { 1., .9794, .6967,
//                   .9794, 1., .7742,
//                   .6967, .7742, 1. };
   double test[16] = { 4, -30, 60, -35,
                      -30, 300, -675, 420,
                      60, -675, 1620, -1050,
                     -35, 420, -1050, 700 };
   double eigenvects[16], eigenvals[4];

   jacobi_eigenvalues( test, 4, eigenvals, eigenvects);
   return( 0);
}
#endif
