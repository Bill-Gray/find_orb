/* pl_cache.h: computes and caches planetary positions

SEE PL_CACHE.TXT (and 'pl_cache.cpp') FOR A DISCUSSION OF WHAT
THIS DOES.  It probably won't make much sense to you if you don't.

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

int planet_posn( const int planet_no, const double jd, double *vect_2000);
int format_jpl_ephemeris_info( char *buff);           /* pl_cache.cpp */
int get_jpl_ephemeris_info( int *de_version, double *jd_start, double *jd_end);

#define PLANET_POSN_VELOCITY_OFFSET 1000

/* Requesting planet 3 gets you the Earth-Moon barycenter.  Requesting
planet 10 gets the vector between Earth and Moon.  Those are the values
stored in the JPL ephemerides,  but it can be useful to say,  "just gimme
the vector from the sun to the earth" or "from the sun to the moon".  So
these values can be fed to planet_posn() in such cases.   */

#define PLANET_POSN_EARTH        20
#define PLANET_POSN_MOON         21
