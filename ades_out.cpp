/* findorb.cpp: main driver for console Find_Orb

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
02110-1301, USA.

   This is a start/rough hack at writing out astrometry in the XML ADES
format.  In particular,  the 'obsContext' data is essentially bogus at
present.  The output is that for a submission,  but no effort yet is
made to break up the input astrometry into a series of submissions,
each from a separate telescope/MPC observatory/program code.  The
current code is suitable for one very specific instance where I've been
tasked to provide ADES output;  a solution to the more general case will
have to wait for another day.       */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "mpc_obs.h"
#include "mpc_func.h"
#include "afuncs.h"

char *get_file_name( char *filename, const char *template_file_name);
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
char *iso_time( char *buff, const double jd);         /* elem_out.cpp */
int get_satellite_offset( const char *iline, double *xyz);  /* mpc_obs.cpp */
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

void create_ades_file( const char *filename, const OBSERVE FAR *obs, int n_obs)
{
   char buff[200];
   FILE *ofile = fopen_ext( get_file_name( buff, filename), "tfcwb");

   fprintf( ofile, "<?xml version=\"1.0\" ?>\n");
   fprintf( ofile, "<ades version=\"2017\">\n");
   fprintf( ofile, "  <obsBlock>\n");
   fprintf( ofile, "    <obsContext>\n");
   fprintf( ofile, "      <observatory>\n");
   fprintf( ofile, "        <mpcCode>%s</mpcCode>\n", obs->mpc_code);
   fprintf( ofile, "      </observatory>\n");
   fprintf( ofile, "      <submitter>\n");
   fprintf( ofile, "        <name>A. U. Submitter</name>\n");
   fprintf( ofile, "      </submitter>\n");
   fprintf( ofile, "      <observers>\n");
   fprintf( ofile, "        <name>Z. Z. Observer</name>\n");
   fprintf( ofile, "      </observers>\n");
   fprintf( ofile, "      <measurers>\n");
   fprintf( ofile, "        <name>Q. Q. Measurer</name>\n");
   fprintf( ofile, "      </measurers>\n");
   fprintf( ofile, "      <telescope>\n");
   fprintf( ofile, "        <design>Coronagraph</design>\n");
   fprintf( ofile, "        <aperture>0.11</aperture>\n");
   fprintf( ofile, "        <detector>CCD</detector>\n");
   fprintf( ofile, "        <name>C2</name>\n");
   fprintf( ofile, "      </telescope>\n");
   fprintf( ofile, "    </obsContext>\n");
   fprintf( ofile, "    <obsData>\n");
   while( n_obs--)
      {
      const double correlation = 0.;
      const char *catalogue = byte_code_to_net_name( obs->astrometric_net_code);

      fprintf( ofile, "      <optical>\n");
      strcpy( buff, obs->packed_id);
      text_search_and_replace( buff, " ", "");
      fprintf( ofile, "        <trkSub>%s</trkSub>\n", buff);
      fprintf( ofile, "        <mode>%s</mode>\n", "CCD");
      fprintf( ofile, "        <stn>%s</stn>\n", obs->mpc_code);
      if( obs->note2 == 'S')
         {
         int i, is_au = 0;
         double posn[3];

         get_satellite_offset( obs->second_line, posn);
         ecliptic_to_equatorial( posn);
         for( i = 0; i < 3; i++)
            if( fabs( posn[i]) > 999999. / AU_IN_KM)
               is_au = 1;
         fprintf( ofile, "        <sys>ICRF_%s</sys>\n",
                                            is_au ? "AU" : "KM");
         fprintf( ofile, "        <ctr>399</ctr>\n");
         for( i = 0; i < 3; i++)
            fprintf( ofile, "        <pos%d>%.*f</pos%d>\n",
                     i + 1, (is_au ? 13 : 4), posn[i] * (is_au ? 1. : AU_IN_KM), i + 1);
         }
      fprintf( ofile, "        <obsTime>%s</obsTime>\n",
               iso_time( buff, obs->jd));         /* elem_out.cpp */
      fprintf( ofile, "        <ra>%.6f</ra>\n", obs->ra * 180. / PI);
      fprintf( ofile, "        <dec>%.6f</dec>\n", obs->dec * 180. / PI);
      fprintf( ofile, "        <rmsRA>%.4f</rmsRA>\n", obs->posn_sigma_1);
      fprintf( ofile, "        <rmsDec>%.4f</rmsDec>\n", obs->posn_sigma_2);
      fprintf( ofile, "        <rmsCorr>%.4f</rmsCorr>\n", correlation);
      if( !catalogue)
         catalogue = "?";
      strcpy( buff, catalogue);
      text_search_and_replace( buff, "-", "");
      text_search_and_replace( buff, " ", "");
      fprintf( ofile, "        <astCat>%s</astCat>\n", buff);
      if( obs->obs_mag != BLANK_MAG)
         fprintf( ofile, "        <mag>%.*f</mag>\n",
                                 obs->mag_precision, obs->obs_mag);
      fprintf( ofile, "      </optical>\n");
      obs++;
      }
   fprintf( ofile, "    </obsData>\n");
   fprintf( ofile, "  </obsBlock>\n");
   fprintf( ofile, "</ades>\n");
   fclose( ofile);
}
