/* orb_func.cpp: computes when/where an object hits a planet or moon

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
#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_MINOR_AXIS 6356755.
#define EARTH_AXIS_RATIO (EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS)

int parallax_to_lat_alt( const double rho_cos_phi, const double rho_sin_phi,
       double *lat, double *ht_in_meters, const int planet_idx); /* ephem0.c */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
double find_collision_time( ELEMENTS *elem, double *latlon, const int is_impact);
double planet_axis_ratio( const int planet_idx);            /* collide.cpp */
double planet_radius_in_meters( const int planet_idx);      /* collide.cpp */
int find_lat_lon_alt( const double jd, const double *ivect,  /* collide.cpp */
                                   double *lat_lon_alt, const bool geometric);
void calc_approx_planet_orientation( const int planet,        /* runge.cpp */
         const int system_number, const double jde, double *matrix);

int debug_printf( const char *format, ...);                /* runge.cpp */

double planet_radius_in_meters( const int planet_idx)
{
   static const double planet_radius[15] = { 695992000., 2439000., 6051000.,
/* earth-saturn */         6378140, 3393000., 71492000., 60268000.,
/* uranus-pluto, luna */   25559000., 24764000., 1195000., 1737400,
/* jupiter moons */        1821300, 1565000., 2634000., 2403000. };

   return( planet_radius[planet_idx]);
}

#define N_FLATTENINGS 9

double planet_axis_ratio( const int planet_idx)
{
   static const double flattenings[N_FLATTENINGS] = {
            0., 0., 0., .00335364,              /* Sun Mer Ven Earth */
            .00647630, .0647630, .0979624,      /* Mars Jup Satu */
            .0229273, .0171 };                  /* Uran Nep */

   return( planet_idx >= N_FLATTENINGS ? 1. : 1. - flattenings[planet_idx]);
}

static double find_geodetic_lat_lon_alt( const double ut, const double *ivect,
                     const int planet_no, double *lat_lon)
{
   double planet_matrix[9];
   double loc[3], loc2k[3], alt_in_meters;
   const double planet_radius_in_au =
          planet_radius_in_meters( planet_no) / AU_IN_METERS;
   int i;

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
   parallax_to_lat_alt( sqrt( loc[0] * loc[0] + loc[1] * loc[1]), loc[2],
                       lat_lon + 1, &alt_in_meters, planet_no);
   return( alt_in_meters);
}

double find_collision_time( ELEMENTS *elem, double *latlon, const int is_impact)
{
   double t_low = (is_impact ? -2. / 24. : 2. / 24.);
   double t_high = 0.;    /* assume impact within 2 hrs */
   double t0 = t_low / 2.;
   const double alt_0 = atof( get_environment_ptr( "COLLISION_ALTITUDE"));
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
                /* Geocentric elems are already in J2000 equatorial coords. */
                /* Others are in ecliptic,  so cvt them to equatorial:      */
      if( elem->central_obj != 3)
         ecliptic_to_equatorial( loc2k);
      jd = elem->perih_time + t0;
      ut = jd - td_minus_ut( jd) / seconds_per_day;
      alt_in_meters = find_geodetic_lat_lon_alt( ut, loc2k, elem->central_obj,
                                    latlon);
      if( alt_in_meters < alt_0)
         t_high = t0;
      else
         t_low = t0;
      t0 = (t_low + t_high) / 2.;
      }
   return( t0);
}

/* 30 Jan 2009:  Rob Matson asked if I could provide ephemerides in
   geodetic lat/lon/alt for 2008 TC3,  in aid of running the resulting
   vector through the atmosphere. Input vector is assumed to be in
   geocentric equatorial J2000 and in AU;  time is assumed to be TT.
   Note shameless cannibalisation of code from above. */

int find_lat_lon_alt( const double jd, const double *ivect,
                                  double *lat_lon_alt, const bool geometric)
{
   double planet_matrix[9];
   double loc[3], loc2k[3];
   double rho_cos_phi, rho_sin_phi;
   const double ut = jd - td_minus_ut( jd) / seconds_per_day;

   loc2k[0] = ivect[0];
   loc2k[1] = ivect[1];
   loc2k[2] = ivect[2];
   calc_planet_orientation( 3, 0, ut, planet_matrix);
               /* cvt J2000 to planet-centric coords: */
   precess_vector( planet_matrix, loc2k, loc);
   lat_lon_alt[0] = -atan2( loc[1], loc[0]);
   if( lat_lon_alt[0] < 0.)           /* keep in 0-360 range */
      lat_lon_alt[0] += PI + PI;

   rho_cos_phi = sqrt( loc[0] * loc[0] + loc[1] * loc[1]);
   rho_sin_phi = loc[2];
   if( geometric)
      {
      lat_lon_alt[1] = atan( rho_sin_phi / rho_cos_phi);
      lat_lon_alt[2] = vector3_length( loc);
      }
   else
      {
      const double planet_radius_in_au =
               planet_radius_in_meters( 3) / AU_IN_METERS;
      double ht_in_meters;

      rho_cos_phi /= planet_radius_in_au;
      rho_sin_phi /= planet_radius_in_au;
      parallax_to_lat_alt( rho_cos_phi, rho_sin_phi,
               lat_lon_alt + 1,  &ht_in_meters, 3);
      lat_lon_alt[2] = ht_in_meters / AU_IN_METERS;
      }
   return( 0);
}
