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

/*  This can be compiled with -DTEST_MAIN to get a simple command-line
test version,  or with -DCGI_VERSION to get code used for an on-line
converter.  See https://www.projectpluto.com/cvt_elem.htm for an
example of the latter usage.

gcc -Wall -Wextra -pedantic -I../include -DCGI_VERSION -o cvt_elem conv_ele.cpp -L../lib -l lunar -lm
*/


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

static void show_angles( const double incl, const double Omega, const double omega)
{
   printf( "incl %11.7f   asc %11.7f   argper %11.7f\n",
               incl * 180. / PI, Omega * 180. / PI, omega * 180. / PI);
}

int main( const int argc, const char **argv)
{
            /* test case from Meeus,  _Astro Algorithms_, p 148: */

   if( argc == 6)
      {
      double incl = atof( argv[1]) * PI / 180.;
      double omega = atof( argv[2]) * PI / 180.;
      double Omega = atof( argv[3]) * PI / 180.;
      double epoch_from = atof( argv[4]);
      double epoch_to   = atof( argv[5]);

      show_angles( incl, omega, Omega);
      convert_elements( epoch_from, epoch_to, &incl, &Omega, &omega);
      printf( "After conversion from %f to %f\n", epoch_from, epoch_to);
      show_angles( incl, omega, Omega);
      }
   else           /* test case from Meeus' _Astronomical Algorithms_ */
      {
      double incl = 47.1220 * PI / 180.;
      double omega = 151.4486 * PI / 180.;
      double Omega = 45.7481 * PI / 180.;

      show_angles( incl, omega, Omega);
      convert_elements( 1744., 1950., &incl, &Omega, &omega);
      printf( "After conversion from 1744 to 1950 :\n");
      show_angles( incl, omega, Omega);
      printf( "Should be: 47.1380   151.4782   48.6037\n");
      convert_elements( 1950., 1744., &incl, &Omega, &omega);
      printf( "Converted back to 1744 :\n");
      show_angles( incl, omega, Omega);
      printf( "Should be:      47.1220   151.4486   45.7481\n");
      }
   return( 0);
}
#endif

#ifdef CGI_VERSION

#include <string.h>
#include <stdio.h>
#ifdef __has_include
   #if __has_include("cgi_func.h")
       #include "cgi_func.h"
   #else
       #error   \
         'cgi_func.h' not found.  This project depends on the 'lunar'\
         library.  See www.github.com/Bill-Gray/lunar .\
         Clone that repository,  'make'  and 'make install' it.
       #ifdef __GNUC__
         #include <stop_compiling_here>
            /* Above line suppresses cascading errors. */
       #endif
   #endif
#else
   #include "cgi_func.h"
#endif

int main( void)
{
   char field[30], buff[100];
   int rval;
   double incl = 0., asc_node = 0., arg_per = 0.;
   double start_epoch = 1950., end_epoch = 2000.;

   printf( "Content-type: text/html\n\n");
   printf( "<pre>");
   avoid_runaway_process( 300);
   rval = initialize_cgi_reading( );
   if( rval <= 0)
      {
      printf( "<p> <b> CGI data reading failed : error %d </b>", rval);
      printf( "This isn't supposed to happen.</p>\n");
      return( 0);
      }
   while( !get_cgi_data( field, buff, NULL, sizeof( buff)))
      {
      if( !strcmp( field, "iota"))
         incl = atof( buff) * PI / 180.;
      else if( !strcmp( field, "arg_per"))
         arg_per = atof( buff) * PI / 180.;
      else if( !strcmp( field, "asc"))
         asc_node = atof( buff) * PI / 180.;
      else if( !strcmp( field, "epoch1"))
         start_epoch = atof( buff);
      else if( !strcmp( field, "epoch2"))
         end_epoch = atof( buff);
      else
         printf( "Unidentified field '%s'\n", field);
      }
   convert_elements( start_epoch, end_epoch, &incl, &asc_node, &arg_per);
   printf( "After conversion to epoch %.3f:\n", end_epoch);
   printf( "Ascending node (&Omega;) = %.6f\n", asc_node * 180. / PI);
   printf( "Argument of periapsis (&omega;) = %.6f\n", arg_per * 180. / PI);
   printf( "Inclination (&iota;) = %.6f\n", incl * 180. / PI);
   return( 0);
}
#endif
