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

/* details.h : code to store/access 'observational details' (header)
data for 80-column MPC-formatted astrometry,  as documented at

https://www.minorplanetcenter.net/iau/info/ObsDetails.html */

void *init_observation_details( void);
int add_line_to_observation_details( void *obs_details, const char *iline);
const char **get_code_details( const void *obs_details, const char *mpc_code);
void free_observation_details( void *obs_details);

/* Return codes for add_line_to_observation_details : */

#define OBS_DETAILS_IRRELEVANT_LINE             0
#define OBS_DETAILS_MPC_80_COLUMN_LINE          1
#define OBS_DETAILS_HEADER_LINE                 2
