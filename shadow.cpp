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
of radius,  with the intensity at the center being 1,  gets you this :

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
intensity (0 if no disk visible,  1 if r=R.  I0 = TI(1) / 2pi.)     */

#include <math.h>

#define NECKEL_LABS_COEFFS
/* #define YOULES_COEFFS   */

static double total_intensity( const double r)
{
#ifdef PIERCE_SLAUGHTER_COEFFS
   const double I0 = 0.4067828571428571;   /* above integral for r=0, mu=1 */
#endif
#ifdef NECKEL_LABS_COEFFS
   const double I0 = 0.4062268095238096;   /* above integral for r=0, mu=1 */
#endif
#ifdef YOULES_COEFFS
   const double I0 = 0.418;                /* above integral for r=0, mu=1 */
#endif
   const double mu2 = 1. - r * r;
   double mu, rval = 0, power = mu2;
   size_t i;
#ifdef PIERCE_SLAUGHTER_COEFFS
   static const double coeffs[6] = { 0.30505, 1.13123, -0.78604,
                                    0.40560,  0.02297, -0.07880 };
#endif
#ifdef NECKEL_LABS_COEFFS
   static const double coeffs[6] = { 0.28392, 1.36896, -1.75998,
                                    2.22154, -1.56074, 0.44630 };
#endif
#ifdef YOULES_COEFFS
   static const double coeffs[3] = { 0.436,   0.72,    -0.16 };
#endif

   if( mu2 > 0.)        /* with roundoff,  we can be slightly negative */
      mu = sqrt( mu2);
   else
      mu = 0.;
   for( i = 0; i < sizeof( coeffs) / sizeof( coeffs[0]); i++, power *= mu)
      rval += coeffs[i] * power / (double)( i + 2);
/* return( r * r);                  for a uniformly intense solar disc */
   return( (I0 - rval) / I0);
}

#include <stdio.h>

int main( const int argc, const char **argv)
{
   double r;

   printf( "I0 = %.16f\n", total_intensity( 0));
   for( r = 0.; r < 1.005; r += (r < 0.94 ? 0.05 : 0.01))
      {
      const double h = 0.001;
      double numeric_deriv =
               (total_intensity( r + h) - total_intensity( r - h)) / (2. * h);

      if( r)
         numeric_deriv /= r;
      printf( "%.2f: %.9f %.5f\n", r, total_intensity( r), numeric_deriv);
      }
   return( 0);
}

#ifdef OLD_CODE_USED_FOR_TABLE

/* I used this to make the table in the comments at top. */

int main( const int argc, const char **argv)
{
   double r;

   printf( " r/R  Youles  PS1977  NL1994\n");
   for( r = 0.; r < 1.005; r += (r < 0.94 ? 0.05 : 0.01))
      {
      double mu2 = 1. - r * r;
      double mu, youles, ps, nl;

      if( mu2 > 0.)
         mu = sqrt( 1. - r * r);
      else
         mu = 0.;
      youles = 0.436 + mu * (0.72 - 0.16 * mu);
      ps = 0.30505 + mu * (1.13123 + mu * (-0.78604 + mu * (0.40560 + mu * ( 0.02297 + mu * -0.07880))));
      nl = 0.28392 + mu * (1.36896 + mu * (-1.75998 + mu * (2.22154 + mu * (-1.56074 + mu *  0.44630))));
      printf( "%.2f: %.5f %.5f %.5f\n", r, youles, ps, nl);
      }
   return( 0);
}
#endif         /* #ifdef OLD_CODE_NOT_IN_USE */
