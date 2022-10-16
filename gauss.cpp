/* gauss.cpp: computes orbits using method of Gauss

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

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "watdefs.h"
#include "mpc_obs.h"
#include "afuncs.h"

FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
double gauss_method( const OBSERVE FAR *obs1, const OBSERVE FAR *obs2,
                     const OBSERVE FAR *obs3, double *orbit, const double mu,
                     const int desired_soln);               /* gauss.cpp */
double convenient_gauss( const OBSERVE FAR *obs, int n_obs, double *orbit,
                  const double mu, const int desired_soln); /* gauss.cpp */
int find_real_polynomial_roots( const double *poly, int poly_order,
                                double *real_roots);        /* roots.cpp */

/* References are to Boulet,  _Methods of Orbit Determination_. */

/* Higher-order terms in the f and g series are given in Boulet,
section 3.6.1 (program FANDG),  starting at line 16010. */

#define INCLUDE_FIFTH
#define INCLUDE_SIXTH

static double dot_prod3( const double FAR *a, const double FAR *b)
{
   return( a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

/* See p. 415, (10.13) and (10.14): we evaluate ten "(A x B) dot C"
   expressions,  so the following function is convenient: */

static double cross_then_dot( const double FAR *a, const double FAR *b,
                              const double FAR *c)
{
   double xprod[3];

   vector_cross_product( xprod, a, b);
   return( dot_prod3( xprod, c));
}

#define GAUSS_K .01720209895
#define P_MAX 10000.

double gauss_method( const OBSERVE FAR *obs1, const OBSERVE FAR *obs2,
                     const OBSERVE FAR *obs3, double *orbit, const double mu,
                     const int desired_soln)
{
   double tau1 = GAUSS_K * (obs1->jd - obs2->jd);
   double tau3 = GAUSS_K * (obs3->jd - obs2->jd);
   double tau  = tau3 - tau1, rval = 0.;
   const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923;

               /* See p. 415, (10.13): */
   const double d11 = cross_then_dot( obs2->vect, obs1->obs_posn, obs3->vect);
   const double d12 = cross_then_dot( obs2->vect, obs2->obs_posn, obs3->vect);
   const double d13 = cross_then_dot( obs2->vect, obs3->obs_posn, obs3->vect);

   const double d21 = cross_then_dot( obs1->obs_posn, obs1->vect, obs3->vect);
   const double d22 = cross_then_dot( obs2->obs_posn, obs1->vect, obs3->vect);
   const double d23 = cross_then_dot( obs3->obs_posn, obs1->vect, obs3->vect);

   const double d31 = cross_then_dot( obs1->obs_posn, obs2->vect, obs1->vect);
   const double d32 = cross_then_dot( obs2->obs_posn, obs2->vect, obs1->vect);
   const double d33 = cross_then_dot( obs3->obs_posn, obs2->vect, obs1->vect);
               /* See p. 415, (10.14): */

   const double d0 = cross_then_dot( obs1->vect, obs2->vect, obs3->vect);

               /* See p. 417, (10.20): */
   const double a1 = tau3 / tau, a3 = -tau1 / tau;
   const double b1 = a1 * (tau * tau - tau3 * tau3) / 6.;
   const double b3 = a3 * (tau * tau - tau1 * tau1) / 6.;


               /* See p. 417, (10.22): */
   const double big_a = -(a1 * d21 - d22 + a3 * d23) / d0;
   const double big_b = -(b1 * d21 + b3 * d23) / d0;

               /* See p. 417, (10.24): */
   const double big_e = 2. * dot_prod3( obs2->obs_posn, obs2->vect);
   const double big_f = dot_prod3( obs2->obs_posn, obs2->obs_posn);
   double f1 = 0., g1 = 0., f3 = 0., g3 = 0.;

               /* Set up the equation of Lagrange,  a degree-8 polynomial.
                  See p. 417, (10.26). */

   double poly[9], roots[8], r2;
   int i, n_roots, iteration, keep_iterating = 1;
   extern int process_count;
   FILE *ofile = (process_count? NULL : fopen_ext( "gauss.out", "tfcwb"));

   if( ofile)
      {
      setvbuf( ofile, NULL, _IONBF, 0);
      fprintf( ofile, "jds: %f %f %f; d0 = %g\n", obs1->jd, obs2->jd, obs3->jd, d0);
      fprintf( ofile, "Loc1: %f %f\n", obs1->ra * 180. / pi, obs1->dec * 180. / pi);
      fprintf( ofile, "Loc2: %f %f\n", obs2->ra * 180. / pi, obs2->dec * 180. / pi);
      fprintf( ofile, "Loc3: %f %f\n", obs3->ra * 180. / pi, obs3->dec * 180. / pi);
      }
   if( fabs( d0) < 1e-10 || !tau1 || !tau3)
      {
      if( ofile)
         fclose( ofile);
      return( 0.);
      }
   for( i = 0; i < 8; i++)
      poly[i] = 0.;
   poly[8] = 1.;
   poly[6] = -(big_a * (big_a + big_e) + big_f);
   poly[3] = -mu * big_b * (2. * big_a + big_e);
   poly[0] = -mu * mu * big_b * big_b;
   if( ofile)
      fprintf( ofile, "poly: x^8 + %f x^6 + %f x^3 + %f\n",
            poly[6], poly[3], poly[0]);

   n_roots = find_real_polynomial_roots( poly, 8, roots);    /* roots.cpp */
   if( ofile)
      {
      fprintf( ofile, "\n%d roots found:\n", n_roots);
      for( i = 0; i < n_roots; i++)
         fprintf( ofile, "Root %d: %f\n", i, roots[i]);
      }
           /* Toss negative roots and those leading to negative distances: */
   for( i = 0; i < n_roots; i++)
      if( roots[i] < 0.
               || big_a + mu * big_b / (roots[i] * roots[i] * roots[i]) < 0.)
         {
         roots[i] = roots[n_roots - 1];
         i--;
         n_roots--;
         }
   if( ofile)
      {
      fprintf( ofile, "%d roots survive\n", n_roots);
      for( i = 0; i < n_roots; i++)
         fprintf( ofile, "'surviving' root %d: %f\n", i, roots[i]);
      }

   if( desired_soln >= n_roots)
      {
      if( ofile)
         fclose( ofile);
      return( 0.);
      }
   if( desired_soln == -1)    /* we're supplying r2 */
      r2 = orbit[0];
   else
      r2 = roots[desired_soln];
   if( ofile)
      fprintf( ofile, "Initial rho %f\n", r2);

   for( iteration = 0; keep_iterating && iteration < 96; iteration++)
      {
            /* See p. 419, (10.28), (10.27): */
      const double u2 = mu / (r2 * r2 * r2);
      const double tau1_squared = tau1 * tau1, tau3_squared = tau3 * tau3;
      double f1_new = 1. - u2 * tau1_squared / 2.;
      double f3_new = 1. - u2 * tau3_squared / 2.;
      double g1_new = tau1 * (1. - u2 * tau1_squared / 6.);
      double g3_new = tau3 * (1. - u2 * tau3_squared / 6.);
      double f1_g3_minus_f3_g1, c1, c2, c3, p1, p2, p3, d1, d3;

      if( iteration)       /* we can refine the above; see (10.15), (3.80) */
         {
         const double z2 = dot_prod3( orbit, orbit + 3) / (r2 * r2);
         const double zs = z2 * z2;
         const double u2_z2_tau1_cubed = u2 * z2 * tau1_squared * tau1;
         const double u2_z2_tau3_cubed = u2 * z2 * tau3_squared * tau3;
         const double q = dot_prod3( orbit + 3, orbit + 3) / (r2 * r2) - u2;
         const double quartic = u2 * (3 * q - 15. * zs + u2) / 24.;
#ifdef INCLUDE_FIFTH
         const double quint_f = u2 * (3 * q +  7. * zs - u2) * z2 / 8.;
         const double quint_g = u2 * (9 * q - 45. * zs + u2) / 120.;
#ifdef INCLUDE_SIXTH
         const double sixth_f = u2 * (zs * (630 * q + 210 * u2 - 945 * zs)
               - u2 * (24 * q + u2) - 45 * q * q) / 720.;
         const double sixth_g = u2 * z2 * (210 * zs - 90 * q - 15 * u2) / 360.;
#endif
#endif
         const double tau1_fourth = tau1_squared * tau1_squared;
         const double tau3_fourth = tau3_squared * tau3_squared;

         f1_new += u2_z2_tau1_cubed / 2. + tau1_fourth * quartic;
         f3_new += u2_z2_tau3_cubed / 2. + tau3_fourth * quartic;
         g1_new += u2_z2_tau1_cubed * tau1 / 4.;
         g3_new += u2_z2_tau3_cubed * tau3 / 4.;
#ifdef INCLUDE_FIFTH
         f1_new += quint_f * tau1_fourth * tau1;
         f3_new += quint_f * tau3_fourth * tau3;
         g1_new += quint_g * tau1_fourth * tau1;
         g3_new += quint_g * tau3_fourth * tau3;
#ifdef INCLUDE_SIXTH
         f1_new += sixth_f * tau1_fourth * tau1_squared;
         f3_new += sixth_f * tau3_fourth * tau3_squared;
         g1_new += sixth_g * tau1_fourth * tau1_squared;
         g3_new += sixth_g * tau3_fourth * tau3_squared;
#endif
#endif
         }
      if( iteration && iteration < 4)
         {
         f1 = (f1 + f1_new) / 2.;      /* for stability reasons,  average */
         f3 = (f3 + f3_new) / 2.;      /* the new f1, f3, g1, and g3 with */
         g1 = (g1 + g1_new) / 2.;      /* values from previous iteration  */
         g3 = (g3 + g3_new) / 2.;
         }
      else        /* on first iteration,  we don't average f1, f3, g1, g3 */
         {        /* with the values from previous iterations... because  */
         f1 = f1_new;            /* there _were_ no previous iterations   */
         f3 = f3_new;
         g1 = g1_new;
         g3 = g3_new;
         }

               /* See p. 419, (10.29): */
      f1_g3_minus_f3_g1 = f1 * g3 - f3 * g1;
      c1 =  g3 / f1_g3_minus_f3_g1;
      c2 = -1.;
      c3 = -g1 / f1_g3_minus_f3_g1;

               /* See p. 419, (10.30): */
      p1 = (c1 * d11 + c2 * d12 + c3 * d13) / (c1 * d0);
      p2 = (c1 * d21 + c2 * d22 + c3 * d23) / (c2 * d0);
      p3 = (c1 * d31 + c2 * d32 + c3 * d33) / (c3 * d0);

               /* See p. 420, (10.33): */
      d1 = -f3 / f1_g3_minus_f3_g1;
      d3 =  f1 / f1_g3_minus_f3_g1;


      if( ofile)
         fprintf( ofile, "%d: r2 = %f\np1 = %f   p2 = %f    p3 = %f\n",
                  iteration, r2, p1, p2, p3);
      p1 = fabs( p1);
      p2 = fabs( p2);
      p3 = fabs( p3);
      if( p1 < 0. || p2 < 0. || p3 < 0.
               || p1 > P_MAX || p2 > P_MAX || p3 > P_MAX)
         {
         rval = 0.;        /* iterations diverging */
         keep_iterating = 0;
         }
      else
         {
         double r2_new_squared = 0., new_r2;
         const double tolerance = 1e-5;

         for( i = 0; i < 3; i++)
            {
            const double r1_component = obs1->obs_posn[i] + p1 * obs1->vect[i];
            const double r2_component = obs2->obs_posn[i] + p2 * obs2->vect[i];
            const double r3_component = obs3->obs_posn[i] + p3 * obs3->vect[i];

            orbit[i] = r2_component;
            orbit[i + 3] = d1 * r1_component + d3 * r3_component;
            rval = p2;
            r2_new_squared += r2_component * r2_component;
            }
         new_r2 = sqrt( r2_new_squared);
         if( fabs( new_r2 - r2) < tolerance)
            keep_iterating = 0;
         r2 = new_r2;
         }
               /* compute better tau values,  including light time lag: */

      tau1 = GAUSS_K * (obs1->jd - obs2->jd - (p1 - p2) / AU_PER_DAY);
      tau3 = GAUSS_K * (obs3->jd - obs2->jd - (p3 - p2) / AU_PER_DAY);
/*    tau  = tau3 - tau1;   We don't actually use the recomputed value */
      if( ofile)
         fprintf( ofile, "state vect: %f %f %f %f %f %f\n",
                           orbit[0], orbit[1], orbit[2],
                           orbit[3], orbit[4], orbit[5]);
      }

   for( i = 3; i < 6; i++)
      orbit[i] *= GAUSS_K;
   if( ofile)
      {
      fprintf( ofile, "rval = %f\n", rval);
      fclose( ofile);
      }
   return( rval);
}

/* If you have more than three observations,  you have the puzzle of
figuring out which three should be used in the method of Gauss.  The
usual reasoning is that they should have the greatest deviation from
a great circle,  i.e.,  the value of d0 in the above function should
be maximized.

   You _could_ do this by looking at all possible subsets of three
observations.  If you have (say) a hundred observations,  you'll be
looking at 100 * 99 * 98 / 6 = 161700 possible triplets,  which seems
like a bit much.  At present,  I just assume the first and last should
be used,  and look for the intervening observation that maximizes
d0 among them.  It might actually make slightly more sense to pick
a few pairs at pseudo-random,  then look for a third observation
that maximizes d0... you'd have a bit more processing time,  but also
some chance of improving d0.  */

static double find_best_obs_for_gauss( const OBSERVE FAR *obs,
               const unsigned n_obs, unsigned *idx)
{
   double max_d0 = 0.;
   unsigned i;
   double xprod[3];

   idx[0] = 0;
   idx[1] = 1;    /* mostly just to avoid an uninitialized var warning */
   idx[2] = n_obs - 1;
   vector_cross_product( xprod, obs[idx[0]].vect, obs[idx[2]].vect);
   for( i = 1; i < n_obs - 1; i++)
      if( obs[i].is_included)
         {
         const double d0 = fabs( dot_prod3( xprod, obs[i].vect));

         if( max_d0 < d0)
            {
            max_d0 = d0;
            idx[1] = i;
            }
         }
   return( max_d0);
}

/* The following function is,  as the name suggests,  "just for the sake
of convenience."  Given an array of observations,  it finds the first
and last valid observations.  Then it looks for the valid observation
that comes closest to the mid-time of the end points, and calls the
above function with the start,  end,  and middle observations. */

double convenient_gauss( const OBSERVE FAR *obs, int n_obs, double *orbit,
                  const double mu, const int desired_soln)
{
   double p2, old_orbit[6], d0;
   unsigned obs_idx[3], i;
   extern double uncertainty_parameter;

   uncertainty_parameter = 99.;
   while( n_obs && !obs->is_included)
      {
      obs++;
      n_obs--;
      }
   while( n_obs && !obs[n_obs - 1].is_included)
      n_obs--;
   if( n_obs < 3)
      return( 0.);
   d0 = find_best_obs_for_gauss( obs, n_obs, obs_idx);
   if( !d0)             /* either there weren't three valid obs,  or   */
      return( 0.);      /* there weren't three that were non-coplanar  */
   for( i = 0; i < 6; i++)
      old_orbit[i] = orbit[i];
   p2 = gauss_method( obs + obs_idx[0], obs + obs_idx[1],
                      obs + obs_idx[2], orbit, mu, desired_soln);
   if( p2 < .001)             /* no valid soln found */
      {
      for( i = 0; i < 6; i++)
         orbit[i] = old_orbit[i];
      return( 0.);
      }
   return( obs[obs_idx[1]].jd - p2 / AU_PER_DAY);
}

