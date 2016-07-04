/* conv_ele.cpp: convert orbital elements from epoch to epoch

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
#include "watdefs.h"
#include "afuncs.h"

const double PI = 3.141592653589793238462643383279502884197169399375;

static void compute_ecliptic_precession_angles( const double epoch_from,
            const double epoch_to, double *eta, double *pi, double *p)
{
   const double big_t = (epoch_from - 2000.) / 100.;
   const double big_t2 = big_t * big_t;
   const double t = (epoch_to - epoch_from) / 100.;
   const double t2 = t * t, t3 = t2 * t;

               /* See Meeus,  _Astro Algorithms_, p 128: */
   *eta = (47.0029 - .06603 * big_t + .000598 * big_t2) * t
                            + (-.03302 + .000598 * big_t) * t2
                            + .000060 * t3;
   *pi =  (174.876384 * 3600.) + 3289.4789 * big_t + .60622 * big_t2
                  - (869.8089 + .50491 * big_t) * t + .03536 * t2;
   *p = (5029.0966 + 2.22226 * big_t - .000042 * big_t2) * t
                            + (1.11113 - .000042 * big_t) * t2
                            - .000006 * t3;
            /* cvt arcseconds to radians: */
   *eta *= PI / (180. * 3600.);
   *pi *= PI / (180. * 3600.);
   *p *= PI / (180. * 3600.);
}

void convert_elements( const double epoch_from, const double epoch_to,
               double *incl, double *asc_node, double *arg_per)
{
   double pi, eta, p, phi, x, y, z;
   const double sin_incl = sin( *incl), cos_incl = cos( *incl);
   double sin_asc_node_minus_pi, cos_asc_node_minus_pi;

   compute_ecliptic_precession_angles( epoch_from, epoch_to, &eta, &pi, &p);
               /* See Meeus,  _Astro Algorithms_, p 147: */
   phi = pi + p;
   cos_asc_node_minus_pi = cos( *asc_node - pi);
   sin_asc_node_minus_pi = sin( *asc_node - pi);
   z = cos_incl * cos( eta) + sin_incl * sin( eta) *  cos_asc_node_minus_pi;
   y = sin_incl * sin_asc_node_minus_pi;
   x = -sin( eta) * cos_incl + cos( eta) * sin_incl * cos_asc_node_minus_pi;
   *asc_node = phi + atan2( y, x);
         /* avoid ill-conditioned cases: */
   if( z > .5 || z < -.5)
      *incl = acos( z);
   else
      {
      *incl = asine( sqrt( x * x + y * y));
      if( z < 0.)                   /* added 11 Sep 2006 to fix some */
         *incl = PI - *incl;        /* retrograde orbit errors       */
      }

   y = -sin( eta) * sin_asc_node_minus_pi;
   x = sin_incl * cos( eta) - cos_incl * sin( eta) * cos_asc_node_minus_pi;
   *arg_per += atan2( y, x);
}

#ifdef TEST_CODE

#include <stdio.h>

void main( int argc, char **argv)
{
   const double PI = 3.141592653589793238462643383279502884197169399375;
            /* test case from Meeus,  _Astro Algorithms_, p 148: */
   double incl = 47.1220 * PI / 180.;
   double omega = 151.4486 * PI / 180.;
   double Omega = 45.7481 * PI / 180.;

   convert_elements( 1744., 1950., &incl, &Omega, &omega);
   printf( "Results:   %f %f %f\n", incl * 180. / PI,
                omega * 180. / PI, Omega * 180. / PI);
   printf( "Should be: 47.1380   151.4782   48.6037\n");
   convert_elements( 1950., 1744., &incl, &Omega, &omega);
   printf( "Converted back: %f %f %f\n", incl * 180. / PI,
                omega * 180. / PI, Omega * 180. / PI);
   printf( "Should be:      47.1220   151.4486   45.7481\n");
}
#endif
