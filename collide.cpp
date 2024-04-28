/* collide.cpp: computes when/where an object hits a planet or moon

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
#include <stdlib.h>
#include "watdefs.h"
#include "comets.h"
#include "afuncs.h"
#include "lunar.h"

#define PI 3.141592653589793238462643383279502884197169399375

int parallax_to_lat_alt( const double rho_cos_phi, const double rho_sin_phi,
       double *lat, double *ht_in_meters, const int planet_idx); /* ephem0.c */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
double find_collision_time( ELEMENTS *elem, double *latlon, const int is_impact);
double find_lat_lon_alt( const double ut, const double *ivect,
                  const int planet_no, double *lat_lon, const bool geometric);

int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;

double find_lat_lon_alt( const double ut, const double *ivect,
                  const int planet_no, double *lat_lon, const bool geometric)
{
   double planet_matrix[9];
   double loc[3], loc2k[3], alt_in_meters;
   const double planet_radius_in_au =
          planet_radius_in_meters( planet_no) / AU_IN_METERS;
   int i;
   double rho_sin_phi, rho_cos_phi;

   loc2k[0] = ivect[0];
   loc2k[1] = ivect[1];
   loc2k[2] = ivect[2];
   calc_planet_orientation( planet_no, 0, ut, planet_matrix);
               /* cvt J2000 to planet-centric coords: */
   precess_vector( planet_matrix, loc2k, loc);
   lat_lon[0] = -atan2( loc[1], loc[0]);
   if( lat_lon[0] < 0.)           /* keep in 0-360 range */
      lat_lon[0] += PI + PI;
   for( i = 0; i < 3; i++)          /* then cvt from AU to planet radii: */
      loc[i] /= planet_radius_in_au;
   rho_cos_phi = sqrt( loc[0] * loc[0] + loc[1] * loc[1]);
   rho_sin_phi = loc[2];
   if( geometric)
      {
      const double planet_radius_in_meters =
               planet_radius_in_au * AU_IN_METERS;

      lat_lon[1] = atan( rho_sin_phi / rho_cos_phi);
      alt_in_meters = vector3_length( loc) * planet_radius_in_meters;
      }
   else
      parallax_to_lat_alt( rho_cos_phi, rho_sin_phi,
                       lat_lon + 1, &alt_in_meters, planet_no);
   return( alt_in_meters);
}

double find_collision_time( ELEMENTS *elem, double *latlon, const int is_impact)
{
   double t_low = (is_impact ? -2. / 24. : 2. / 24.);
   double t_high = 0.;    /* assume impact within 2 hrs */
   double t0 = t_low / 2.;
   const double alt_0 = (elem->central_obj == 3 ?
                   atof( get_environment_ptr( "COLLISION_ALTITUDE")) : 0.);
   const double planet_radius_in_au =
          planet_radius_in_meters( elem->central_obj) / AU_IN_METERS;
   int iter = 25;

// debug_printf( "elem->q = %f, center %d\n", elem->q, elem->central_obj);
   if( elem->q > planet_radius_in_au + alt_0 / AU_IN_METERS)
      t0 = 1.;        /* no collision possible */
   else while( iter--)
      {
      double loc2k[4], vel2k[3];
      double alt_in_meters, ut, jd;

      comet_posn_and_vel( elem, elem->perih_time + t0, loc2k, vel2k);
      ecliptic_to_equatorial( loc2k);
      jd = elem->perih_time + t0;
      ut = jd - td_minus_ut( jd) / seconds_per_day;
      alt_in_meters = find_lat_lon_alt( ut, loc2k, elem->central_obj,
                                    latlon, false);
      if( alt_in_meters < alt_0)
         t_high = t0;
      else
         t_low = t0;
      t0 = (t_low + t_high) / 2.;
      }
   return( t0);
}
