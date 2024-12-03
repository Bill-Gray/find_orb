/* Code to compute the intensity of sunlight as a function of radius
on the solar disk (limb darkening),  and to compute the integrated
intensity of the sunlight within a given radius.  This has two uses
within Find_Orb :

   -- Find_Orb adjusts the brightness of an asteroid or artsat downward
as it passes through the earth's shadow.  Originally,  this was done
solely according to the portion of the sun's area blocked by the earth,
as if the sun's disk was uniformly bright.  It actually darkens as one
approaches the limb :

https://en.wikipedia.org/wiki/Limb_darkening

   -- The same fraction of illumination during partial eclipses is used
for computing non-gravitational effects on artsats.  For some lightweight
objects,  this is a surprisingly strong effect;  eclipses matter.  It
_may_ be the case that the non-uniform solar disk also matters,  though
I've not seen such a case yet.

   Find_Orb now splits the sun's disk into a series of rings and
computes,  for each ring,  the fraction of the sun's intensity within
that ring and the fraction of that ring which is uneclipsed.  Summing
up those ring contributions gets you a value between 0% (total eclipse)
and 100% (no eclipse).

   A formula on page 216 of
https://ui.adsabs.harvard.edu/abs/2018MPBu...45..215B/abstract gives
a decent approximation of disk brightness as a function of radius,
as found by S. Youles :

I = 0.436 + 0.72 * mu - 0.16 mu^2

   where I = intensity and mu^2 = 1 - (r/R)^2,  r = distance from
the center of the solar disk and R = radius of the solar disk.  It
refers to a no-longer-available source,

Youles, S. (2017). "PH3010: Solar Limb Darkening."
https://twiki.ph.rhul.ac.uk/twiki/pub/Public/SamYoules/PH3010_Solar_limb_darkening.pdf

   Two similar polynomials that include terms to mu^5 are given in

https://www.physics.hmc.edu/faculty/esin/a101/limbdarkening.pdf

   The polynomials are from Pierce & Slaughter (1977) and Neckel & Labs
(1994).   Reference for the latter:

http://articles.adsabs.harvard.edu/pdf/1994SoPh..153...91N

    While the coefficients of the two quintics are very different,  they
give similar results until you get close to the edge,  and even there,
the maximum difference is about 0.02.  Comparing intensities as a function
of radius gets you this.  Note that they're normalized to have intensity=1
at the center (with some roundoff occurring).

 r/R  Youles  PS1977  NL1994
0.00: 0.99600 1.00001 1.00000
0.05: 0.99550 0.99942 0.99937
0.10: 0.99399 0.99763 0.99748
0.15: 0.99145 0.99462 0.99432
0.20: 0.98785 0.99035 0.98985
0.25: 0.98314 0.98475 0.98404
0.30: 0.97724 0.97775 0.97684
0.35: 0.97006 0.96923 0.96815
0.40: 0.96149 0.95907 0.95787
0.45: 0.95138 0.94712 0.94586
0.50: 0.93954 0.93317 0.93192
0.55: 0.92572 0.91698 0.91578
0.60: 0.90960 0.89822 0.89709
0.65: 0.89075 0.87646 0.87540
0.70: 0.86858 0.85109 0.85003
0.75: 0.84224 0.82119 0.82005
0.80: 0.81040 0.78527 0.78399
0.85: 0.77088 0.74070 0.73934
0.90: 0.71944 0.68197 0.68091
0.95: 0.64522 0.59397 0.59390
0.96: 0.62506 0.56908 0.56919
0.97: 0.60158 0.53944 0.53955
0.98: 0.57294 0.50224 0.50184
0.99: 0.53438 0.45013 0.44766
1.00: 0.43600 0.30505 0.28392

   The PS and NL quintics agree fairly well.  The quadratic from Youles
is a decent approximation,  but not quite as good. Somewhat at random,
I've gone with NL (with PS as an #ifdeffed-out option.)

   Normalizing R=1,  the total intensity within a radius r0 is (you may
need fixed-width fonts to read the following)

              r0                 r0
             (                  (
TI(r0) = 2pi  \   rI(r) dr = 2pi \  r (a0 + a1 * mu + a2 * mu^2 + ...) dr
               )                  )
             r=0                r=0

   where,  as above,  mu = sqrt(1-r^2).  Conveniently,  though,
dmu/dr = -r/mu,  dr = -(mu/r) dmu,  and we can change variables :

             1
            (
TI(r0) = 2pi \  (a0 * mu + a1 * mu2 + a2 * mu^3 + ...) dmu
              )
            mu=mu0

   which is trivially integrable analytically.  In the following function,
the 2pi factor is ignored;  all we really care about is the fractional
intensity (0 if no disk visible,  1 if r=R : the 'normalized' intensity
I(r) = TI(r) / TI(r0).)

   Compiles with

g++ -Wall -Wextra -pedantic -oshadow shadow.cpp       */

#include <math.h>

#define YOULES                      0
#define PIERCE_SLAUGHTER            1
#define NECKEL_LABS                 2
#define UNIFORM_DISK                3

static const double pierce_slaughter_coeffs[6] = {
            0.30505, 1.13123, -0.78604, 0.40560,  0.02297, -0.07880 };
static const double neckel_labs_coeffs[6] = {
            0.28392, 1.36896, -1.75998, 2.22154, -1.56074, 0.44630 };
static const double youles_coeffs[6] = {
            0.436, 0.72, -0.16, 0., 0., 0.};
static const double *coeffs[4] = { youles_coeffs,
          pierce_slaughter_coeffs, neckel_labs_coeffs };

static double intensity( const double r, const int method)
{
   const double mu2 = 1. - r * r;
   double mu, rval = 0, power = 1.;
   size_t i;

   if( mu2 > 0.)        /* with roundoff,  we can be slightly negative */
      mu = sqrt( mu2);
   else
      mu = 0.;
   for( i = 0; i < 6; i++, power *= mu)
      rval += coeffs[method][i] * power;
   return( rval);
}

static double total_intensity( const double r, const int method)
{
   const double mu2 = 1. - r * r;
   double mu, rval = 0, power = mu2;
   size_t i;
   static double i0;
   static int curr_method = -1;

   if( method == UNIFORM_DISK)
      return( r * r);
   if( method != curr_method)    /* recalculate i0 */
      {
      i0 = 0.;
      curr_method = method;
      for( i = 0; i < 6; i++)
         i0 += coeffs[method][i] / (double)( i + 2);
      }
   if( mu2 > 0.)        /* with roundoff,  we can be slightly negative */
      mu = sqrt( mu2);
   else
      mu = 0.;
   for( i = 0; i < 6; i++, power *= mu)
      rval += coeffs[method][i] * power / (double)( i + 2);
   return( (i0 - rval) / i0);
}

#include <stdio.h>

#ifdef OLD_TEST_CODE
      /* I used this to verify that the numerical derivative of the
      'total_intensity()' function matches,  after normalization,
      the results from the 'intensity()' function.  Shouldn't have
      to do that,  but I didn't entirely trust myself on the above
      derivation... */

int main( const int argc, const char **argv)
{
   double r;
   const int method = (argc < 2 ? PIERCE_SLAUGHTER : atoi( argv[1]));

   printf( "I0 = %.16f\n", total_intensity( 0, method));
   for( r = 0.; r < 1.005; r += (r < 0.94 ? 0.05 : 0.01))
      {
      const double h = 0.001;
      double numeric_deriv =
               (total_intensity( r + h, method) - total_intensity( r - h, method)) / (2. * h);

      if( r)
         numeric_deriv /= r;
      printf( "%.2f: %.9f %.5f\n", r, total_intensity( r, method), numeric_deriv);
      }
   return( 0);
}

#endif         /* #ifdef OLD_TEST_CODE */

/* Default : generate table in the comments at top. */

int main( const int argc, const char **argv)
{
   double r;
   int method;

   printf( " r/R  Youles  PS1977  NL1994\n");
   for( r = 0.; r < 1.005; r += (r < 0.94 ? 0.05 : 0.01))
      {
      printf( "%.2f:", r);
      for( method = 0; method < 3; method++)
         printf( " %.5f", intensity( r, method));
      printf( "\n");
      }
   return( 0);
}
