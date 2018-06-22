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

#include <math.h>
#include <assert.h>
#ifdef TEST_CODE
#include <stdio.h>
#include <stdlib.h>
#endif          /* #ifdef TEST_CODE */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

/* Code to convert a "standard" astrometric error ellipse into one that has
been adjusted for timing error.  As part of this,  code is provided to convert
an error ellipse (major axis,  minor axis,  position angle) into a covariance
matrix/quadratic form,  and vice versa.

   The simplest case of dealing with timing errors would be one in which you
have astrometric measurements with a simple Gaussian circular error,  but also
think you have a Gaussian-distributed error in timing.  In that situation,
you would expect the timing error to "smear" your measured position along
the line of motion.  If,  say,  your astrometry was good to 0.7 arcseconds,
your timing was good to about six seconds,  and the object was moving about
0.4" per second,  you might say:  'That timing error alone translates into
a positional uncertainty of 6 * 0.4 = 2.4".  The object's uncertainty ellipse
is therefore probably still only about 0.7" perpendicular to the direction
of motion,  but must be at least 2.4" along the direction of motion.'

   This turns out to be correct,  except that along the direction of motion,
you have to combine the 0.7" error due to astrometric measurement with the
2.4" error due to timing,  and do so in quadrature:  the actual error is
sqrt( 2.4^2 + 0.7^2) = 2.5" (values carefully chosen to make exact math).

   So correcting observations with circular uncertainties for timing error
is really quite straightforward,  and for some time,  it was the only way
in which I handled errors.  However,  some sources provide astrometric data
with elliptical uncertainties (AstDyS and NEODyS,  for example),  and
everybody probably should do this.  (Though one could argue that in many
cases,  the errors are close enough to circular that the small deviation
between major and minor axes could be ignored.)  For some time,  I dealt
with this by taking the geometric mean:  if AstDyS or NEODyS told me an
observation had a sigma of 0.4" in RA and of 0.9" in declination,  I
assigned that observation a uniform sigma of 0.6" on both axes.

   But I still felt that this was stupid,  or at least not a very good idea.
Occasionally,  I'd find datasets where the uncertainty was highly elliptical,
and I knew I wasn't treating the data correctly.  Fortunately,  it turns
out that dealing with elliptical uncertainties,  and correcting them for
timing error,  is not as difficult as I'd envisioned.  Here goes:

   An error ellipse with major axis a,  minor axis b,  at angle theta,
results in a probability density function

(1)  phi0(x, y) = k1 * exp( -(Ax^2 + 2Bxy + Cy^2) / 2)
                = k1 * exp( -((u/a)^2 + (v/b)^2) / 2)

   where u = x * cos( theta) + y * sin( theta) = (x, y) dotted with a unit
   vector pointing at angle theta,  and v = -x * sin( theta) + y * cos( theta)
   = (x, y) dotted with a unit vector pointing at right angles to that.  Run
   the math through,  and you can express A,  B,  and C in terms of a, b,
   theta,  as is done in convert_error_ellipse_to_quadratic_form( ) below.
   (The normalization constant k1,  along with similar normalization constants
   that crop up after adjusting for timing errors, turns out to be irrelevant
   to our current humble purposes,  which is why I'm not bothering to evaluate
   them.)  In covariance matrix form,

                              /x\  / A  B \
(2)    phi0(x, y) = k1 * exp(-| |  |      | (x y) / 2)
                              \y/  \ B  C /

   ...which,  if you plug through the matrix multiplies,  gets you the same
   form as shown above in (1).

   Finding the eigenvalues and eigenvectors of this matrix formulation allows
   us to reverse the process,  determining the error ellipse axes and angles
   from the quadratic form variables A, B, C,  so we can have functions to
   convert easily between error ellipse form and quadratic form.

   After adjusting for timing errors,  one has converted the original
   error ellipse/covariance errors into somewhat different,  "adjusted"
   errors.  To see this,  consider the effect if the object is imaged at
   nominal time t=0,  with timing error sigma_t,  while moving at velocity
   (v_x, v_y).  The adjusted probability density function phi is given as

        inf
        (
phi = k2 \  phi0( x + v_x * t, y + v_y * t) * exp( -(t / sigma_t)^2 / 2) dt
          )
       -inf

   where again,  k2 is a normalization constant that doesn't matter here.
Substitute u = t / sigma_t, vx = v_x * sigma_t, vy = v_y * sigma_t,  and we
get a somewhat simpler

            inf
            (
(3) phi = k3 \  phi0( x + vx * u, y + vy * u) * exp( -u^2 / 2) du
              )
           -inf

   Using (1) above for the definition of phi0,  we can get...

(4) phi0( x + vx * u, y + vy * u) =
           k1 * exp( -(Ax^2 + 2Axvx * u + Avx^2 * u^2
                      +2Bxy  + 2B(yvx + xvy) * u + 2Bvxvy * u^2
                      +Cy^2 + 2Cyvy * u + Cvy^2 * u^2) / 2)

   Combining terms leads us to a simpler form for the integrand in (3):

            inf
            (
(5) phi = k3 \  exp( -(au^2 + 2bu + c) / 2) du
              )
           -inf

   where    a = Avx^2 + 2Bvxvy + Cvy^2 - 1
            b = Axvx + B(xvy + yvx) + Cyvy
              = (Avx + Bvy)x + (Bvx + Cvy)y
            c = Ax^2 + 2Bxy + Cy^2

   The above integral can be found in standard texts.   (Note that c is the
same quadratic used in phi0 in (1),  which is how we get to the last part
of the following.)

(6) phi = k4 * exp( -(c - b^2/a) / 2) = k5 * exp( -(b^2/a) / 2) * phi0(x, y)


   If we set Fx = (Avx + Bvy) and Fy = (Bvx + Cvy),  then b = xFx + yFy,
   b^2 = Fx^2 * x^2 + 2(FxFy) * xy + Fy^2 * y^2,  and

(7) phi = k5 * exp( -((A - Fx * Fx / a)x^2
                   + 2(B - Fx * Fy / a)xy
                   +  (C - Fy * Fy / a)y^2) / 2)

   Set A1 = A - Fx * Fx / a, B1 = B - Fx * Fy / a, C1 = C - Fy * Fy / a,  and
   this all becomes...

(8) phi = k5 * exp( -(A1 x^2 + 2B1xy + C1y^2) / 2)

   ...i.e.,  we have reduced it to a new quadratic form/covariance matrix.
*/

void adjust_error_ellipse_for_timing_error( double *sigma_a, double *sigma_b,
         double *angle, const double vx, const double vy);   /* errors.cpp */
void convert_ades_sigmas_to_error_ellipse( const double sig_ra,
         const double sig_dec, const double correl, double *major,
         double *minor, double *angle);                      /* errors.cpp */

static void adjust_quadratic_form_for_timing_error( const double A,
         const double B, const double C, const double vx, const double vy,
         double *A1, double *B1, double *C1)
{
   const double E = A * vx * vx + 2. * B * vx * vy + C * vy * vy - 1.;
   const double Fx = A * vx + B * vy;
   const double Fy = C * vy + B * vx;

   *A1 = A - Fx * Fx / E;
   *B1 = B - Fx * Fy / E;
   *C1 = C - Fy * Fy / E;
#ifdef TEST_CODE
   printf( "E = %f; Fx = %f; Fy = %f\n", E, Fx, Fy);
   printf( "Bits: %f %f %f\n", Fx * Fx / E, Fx * Fy / E, Fy * Fy / E);
   printf( "Results: %f %f %f\n", *A1, *B1, *C1);
#endif          /* #ifdef TEST_CODE */
}

static void convert_error_ellipse_to_quadratic_form( const double a,
            const double b, const double angle,
            double *A, double *B, double *C)
{
   const double cos_ang = cos( angle), sin_ang = sin( angle);
   const double a2 = a * a, b2 = b * b;

   *A = -cos_ang * cos_ang / a2 - sin_ang * sin_ang / b2;
   *B = -sin_ang * cos_ang * (1. / a2 - 1. / b2);
   *C = -sin_ang * sin_ang / a2 - cos_ang * cos_ang / b2;
}

/* The eigenvalues of / A  B \
                      \ B  C /  are the roots of the quadratic
(lambda-A)(lambda-C) - B^2.  'eigenval2' is the larger eigenvalue,
corresponding to the minor axis.  We get 'eigenval1' by dividing the
constant term of the quadratic by 'eigenval2';  this avoids possible
precision problems that can crop up when you're taking the difference of
two similar quantities.  (Though it may not avoid such problems,  if
AC is close to B^2.... not much we can do about that,  though.)  */

static void convert_quadratic_form_to_error_ellipse( const double A,
         const double B, const double C, double *a, double *b,
         double *angle)
{
   const double tval = sqrt( (A - C) * (A - C) + 4. * B * B);
   const double eigenval2 = (A + C - tval) * .5;
   const double eigenval1 = (A * C - B * B) / eigenval2;

#ifdef TEST_CODE
   printf( "Eigenvals %f %f\n", eigenval1, eigenval2);
#endif          /* #ifdef TEST_CODE */
   assert( eigenval1 < 0.);
   assert( eigenval2 < 0.);
   *a = 1. / sqrt( -eigenval1);
   *b = 1. / sqrt( -eigenval2);
   *angle = atan2( eigenval1 - A, B);
}

/* ADES gives uncertainties in RA and dec,  plus their correlation
(between -1 and +1).  That becomes the correlation matrix

/ a  b \    a = -sig_ra^2     b = -correl * sig_ra * sig_dec
|      |
\ b  c /    c = -sig_dec^2

   Then the quadratic form we want is the inverse of the above
matrix,

/ A B \       / a b \      1    / c b \
|     | = inv |     | =  ------ |     |
\ B C /       \ b c /    ac-b^2 \ b a /

   and we can feed said quadratic form through the above function
to get the error ellipse,  which is what Find_Orb actually wants. */

void convert_ades_sigmas_to_error_ellipse( const double sig_ra,
         const double sig_dec, const double correl, double *major,
         double *minor, double *angle)
{
   const double a = -sig_ra * sig_ra;
   const double b = -sig_ra * sig_dec * correl;
   const double c = -sig_dec * sig_dec;
   const double det = a * c - b * b;
   const double A = c / det;
   const double B = b / det;
   const double C = a / det;

   convert_quadratic_form_to_error_ellipse( A, B, C, major, minor, angle);
}

/* adjust_error_ellipse_for_timing_error( ) puts the above pieces
together : given the estimated error ellipse and the uncertainty
vector from timing,  it computes an adjusted error ellipse "stretched
out" in the direction of motion.   */

void adjust_error_ellipse_for_timing_error( double *sigma_a, double *sigma_b,
         double *angle, const double vx, const double vy)
{
   double A, B, C;
   double A1, B1, C1;

   convert_error_ellipse_to_quadratic_form( *sigma_a,
               *sigma_b, *angle, &A, &B, &C);

   adjust_quadratic_form_for_timing_error( A, B, C, vx, vy,
                  &A1, &B1, &C1);
   convert_quadratic_form_to_error_ellipse( A1, B1, C1,
                 sigma_a, sigma_b, angle);
}


#ifdef TEST_CODE
int main( const int argc, const char **argv)
{
   const double sigma_a = atof( argv[1]);
   const double sigma_b = atof( argv[2]);
   const double theta = atof( argv[3]) * PI / 180.;
   double A, B, C, a, b, angle;

   convert_error_ellipse_to_quadratic_form( sigma_a,
               sigma_b, theta, &A, &B, &C);
   printf( "Quad form: %f %f %f\n", A, B, C);
   convert_quadratic_form_to_error_ellipse( A, B, C,
               &a, &b, &angle);
   printf( "Converted back: %f %f at angle %f\n",
               a, b, angle * 180. / PI);
   if( argc >= 6)
      {
      const double vx = atof( argv[4]);
      const double vy = atof( argv[5]);
      double A1, B1, C1;

      adjust_quadratic_form_for_timing_error( A, B, C, vx, vy,
                  &A1, &B1, &C1);
      printf( "Quad form after adjustment: %f %f %f\n", A1, B1, C1);
      convert_quadratic_form_to_error_ellipse( A1, B1, C1,
                 &a, &b, &angle);
      printf( "After timing error: %f %f at angle %f\n",
               a, b, angle * 180. / PI);
      }
   return( 0);
}
#endif          /* #ifdef TEST_CODE */
