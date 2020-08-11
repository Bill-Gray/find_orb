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
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include "watdefs.h"
#include "mpc_obs.h"
#include "mpc_func.h"
#include "details.h"
#include "afuncs.h"

char *get_file_name( char *filename, const char *template_file_name);
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
char *iso_time( char *buff, const double jd, const int precision);   /* elem_out.c */
int get_satellite_offset( const char *iline, double *xyz);  /* mpc_obs.cpp */
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

static void output_names( FILE *ofile, const OBSERVE FAR *obs, const char *target)
{
   int i;
   size_t tlen = strlen( target);
   const char *tptr, *tptr2;

   if( obs->obs_details)
      for( i = 0; (tptr = obs->obs_details[i]) != NULL; i++)
         if( !strncmp( tptr, target, tlen))
            {
            tptr += tlen;
            while( tptr && *tptr)
               {
               if( *tptr == ',')
                  tptr++;
               while( *tptr == ' ')
                  tptr++;
               tptr2 = strchr( tptr, ',');
               tlen = (tptr2 ? tptr2 - tptr : strlen( tptr));
               if( tlen)
                  fprintf( ofile, "        <name>%.*s</name>\n", (int)tlen, tptr);
               tptr = tptr2;
               }
            }
}

/* Column 14 ('note 1') of punched-card astrometry contains either an
alphabetical note (see https://minorplanetcenter.net/iau/info/ObsNote.html
for their meanings) or a program code (was on the MPC site,  not
currently available).         */

static char program_code( const OBSERVE *obs)
{
   if( strchr( "ABGHOP", obs->note1) && strstr( "249 C49 C50", obs->mpc_code))
      return( obs->note1);    /* SOHO and STEREO program codes */
   return( isalpha( obs->note1) ? ' ' : obs->note1);
}

/* Outputs _only_ those observations from the station and program code
specified by the first observation.      */

static void create_ades_file_for_one_code( FILE *ofile,
                   const OBSERVE FAR *obs, int n_obs)
{
   char buff[200];
   const char *code = obs->mpc_code;
   const char progcode = program_code( obs);

   fprintf( ofile, "  <obsBlock>\n");
   fprintf( ofile, "    <obsContext>\n");
   fprintf( ofile, "      <observatory>\n");
   fprintf( ofile, "        <mpcCode>%s</mpcCode>\n", obs->mpc_code);
   fprintf( ofile, "      </observatory>\n");
   fprintf( ofile, "      <submitter>\n");
   output_names( ofile, obs, "CON ");
   fprintf( ofile, "      </submitter>\n");
   fprintf( ofile, "      <observers>\n");
   output_names( ofile, obs, "OBS ");
   fprintf( ofile, "      </observers>\n");
   fprintf( ofile, "      <measurers>\n");
   output_names( ofile, obs, "MEA ");
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
      if( !strcmp( obs->mpc_code, code) && program_code( obs) == progcode)
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
//       if( progcode != ' ')
//          fprintf( ofile, "        <prog>%c</prog>\n", progcode);
         fprintf( ofile, "        <obsTime>%s</obsTime>\n",
                  iso_time( buff, obs->jd, 3));         /* elem_out.cpp */
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
            {
            fprintf( ofile, "        <mag>%.*f</mag>\n",
                                    obs->mag_precision, obs->obs_mag);
            fprintf( ofile, "        <band>%c</band>\n", obs->mag_band);
            }
         fprintf( ofile, "      </optical>\n");
         }
      obs++;
      }
   fprintf( ofile, "    </obsData>\n");
   fprintf( ofile, "  </obsBlock>\n");
}

/* We dig through the entire set of observations looking for new codes.
When we find one,  we call the above function to create or add ADES data
for just that code.  We also add that code to the 'codes' string so we
do all this only once per code.   */

#define OBSCODE_BUFF_SIZE 20000

void create_ades_file( const char *filename, const OBSERVE FAR *obs,
                                           int n_obs)
{
   int i, j;
   char *codes = (char *)malloc( OBSCODE_BUFF_SIZE);
   char buff[200];
   FILE *ofile = fopen_ext( get_file_name( buff, filename), "tfcwb");

   *codes = '\0';
   setbuf( ofile, NULL);
   fprintf( ofile, "<?xml version=\"1.0\" ?>\n");
   fprintf( ofile, "<ades version=\"2017\">\n");
   for( i = 0; i < n_obs; i++)
      {
      char new_search[5];

      strcpy( new_search, obs[i].mpc_code);
      new_search[3] = program_code( obs + i);
      new_search[4] = '\0';
      j = 0;
      while( codes[j] && memcmp( codes + j, new_search, 4))
         j += 4;
      if( !codes[j])
         {
         assert( strlen( codes) + 1 < OBSCODE_BUFF_SIZE);
         strcat( codes, new_search);
         create_ades_file_for_one_code( ofile, obs + i, n_obs - i);
         }
      }
   free( codes);
   fprintf( ofile, "</ades>\n\n");
   fclose( ofile);
}
