/* The following function performs Lagrangian interpolation through n_pts
points with x=0, 1, 2, ... n_pts-1 and y=y[0], y[1],... y[n_pts - 1],
determining the y corresponding to a given x.  It does this with the
usual Lagrangian polynomial,  of order n_pts-1,

      n_pts-1
      ____
      \
y(x) = >    y[i] poly( i, x)
      /___
       i=0


              __
poly( i, x) = || (x - j) / (i - j)
            j=0..n_pts-1,  j != i

   To see how this works,  consider what y(x) will be for x=0...n_pts-1.
poly(i, x) will be equal to zero for all i not equal to x,  and will
be equal to 1 for i=x.  So we've found an n_pts-1 order polynomial that
fits though our n_pts points.

   At first glance,  it may look as if this is an order n^2 process,
but there's a bit of trickery we can do to make it order n.  First,
if we compute

     __
C =  || (x - j)
   j=0..n_pts-1

   (note the removal of the j != i condition),  we can pull C in front
of the sigma and get

        n_pts-1
         ____
         \
y(x) = C  > y[i] / (t(i) * (x - i))
         /___
          i=0

   where

       __
t(i) = || (i - j)
      j=0..n_pts-1,  j != i

   It still looks as if we're at order n^2.  We can do either of two
things.  We can pre-compute t(i) and use it over and over,  recomputing
only when n_pts changes.  Or we can compute t(0) and then make use of
recurrence:

t(i + 1) = t(i) * i / (i - n_pts)

   ...which is the approach used below.

   Also,  note that when we computed C,  we left ourselves open to a
divide-by-zero error if x was on an abscissa.  That happens iff C=0.
So if C=0,  we just return with the abscissa value y[(int)x]. */

static double interpolate( const double *y, const double x, const int n_pts)
{
   double t = 1., c = 1., rval;
   int i;

   for( i = 0; i < n_pts; i++)
      {
      c *= x - (double)i;
      if( i)
         t *= -(double)i;
      }
   if( !c)        /* we're on an abscissa */
      rval = y[(int)( x + .5)];
   else
      {
      rval = y[0] / (t * x);
      for( i = 1; i < n_pts; i++)
         {
         t *= (double)i / (double)( i - n_pts);
         rval += y[i] / (t * (x - (double)i));
         }
      rval *= c;
      }
   return( rval);
}
 /* A little test program.  Run it with no command line arguments,
and it just produces a table showing how the interpolation works with
five simulated data points. */

#include <stdlib.h>
#include <graph.h>
#include <conio.h>
#include <stdio.h>

void main( int argc, char **argv)
{
   double t;
   double grid_points[6] = { -2., -1.2, 0.2, 1.4, 2.8, 4.0 };
   int i, j, n_pts;

   if( argc == 1)
      n_pts = 5;
   else
      n_pts = atoi( argv[1]);
   _setvideomode( _VRES16COLOR);
   for( t = -1.; t < 5.1; t += .02)
      {
      double y = interpolate( grid_points, t, n_pts);

      _setpixel( (short)( (t + 1.) * 100.),
                 (short)( 240. -  y * 30.));
      }
   _setcolor( 3);
   for( i = 0; i < 7; i++)
      for( j = -9; j < 10; j++)
         _setpixel( (short)( i * 100), (short)( 240 - j * 30));
   for( i = 0; i < n_pts; i++)
      {
      int x = (int)( (i + 1) * 100);
      int y = (int)( 240. - grid_points[i] * 30.);

      _moveto( (short)( x - 3), (short)y);
      _lineto( (short)( x + 3), (short)y);
      _moveto( (short)x, (short)( y - 3));
      _lineto( (short)x, (short)( y + 3));
      }
   getch( );
   _setvideomode( _DEFAULTMODE);
   for( t = -1.; t < 5.1; t += .2)
      printf( "%5.1lf %7.3lf\n", t, interpolate( grid_points, t, n_pts));
}
