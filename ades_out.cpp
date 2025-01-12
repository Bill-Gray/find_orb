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
#include <stringex.h>
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
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
double utc_from_td( const double jdt, double *delta_t);     /* ephem0.cpp */
char *iso_time( char *buff, const double jd, const int precision);   /* elem_out.c */
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
double original_observed_ra( const OBSERVE *obs);     /* ephem0.cpp */
double original_observed_dec( const OBSERVE *obs);    /* ephem0.cpp */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

static const char *_skip_number( const char *tptr)
{
   while( isdigit( *tptr))
      tptr++;
   if( *tptr == '.')
      {
      tptr++;
      while( isdigit( *tptr))
         tptr++;
      }
   return( tptr);
}

/* Column 14 ('note 1') of punched-card astrometry contains either an
alphabetical note (see https://minorplanetcenter.net/iau/info/ObsNote.html
for their meanings) or a program code :

https://minorplanetcenter.net//iau/lists/ProgramCodes.txt
https://minorplanetcenter.net/static/downloadable-files/program_codes.json

We do have the problem that many program codes are now letters.  If column
14 contains a letter,  and that letter is in use as a program code,  we
can't tell if the letter means "program code" or "note".  For my current
problem -- sungrazing astrometry -- column 14 will always indicate the
instrument used (see bottom of 'ObsNotes.html'), but that will not work
reliably with other obscodes.  */

static char program_code( const OBSERVE *obs)
{
   if( strchr( "ABGHOP", obs->note1) && strstr( "249 C49 C50", obs->mpc_code))
      return( obs->note1);    /* SOHO and STEREO program codes */
   return( isalpha( obs->note1) ? ' ' : obs->note1);
}

/* Given a line such as 'OBS E. E. Barnard, M. Wolf',  this will write out

        <name>E. E. Barnard</name>
        <name>M. Wolf</name>

   i.e.,  it 'ADES-izes' old-style observer details. */

static int dump_one_line_of_names( FILE *ofile, const char *line)
{
   int n_found = 0;

   if( !memcmp( line, "TEL ", 4))         /* telescope lines are 'special' */
      {
      const char *tptr = _skip_number( line + 4);

      if( tptr != line + 4 && *tptr == '-' && tptr[1] == 'm' && tptr[2] == ' ')
         {
         char aperture[9], f_ratio[9];
         const size_t aperture_len = tptr - (line + 4);

         assert( aperture_len < 7);
         strlcpy( aperture, line + 4, aperture_len + 1);
         tptr += 3;
         *f_ratio = '\0';
         if( tptr[0] == 'f' && tptr[1] == '/')
            {
            if( tptr[2] == '?')
               tptr += 4;
            else
               {
               const char *f_ratio_ptr = tptr + 2;
               tptr = _skip_number( f_ratio_ptr);
               assert( tptr > f_ratio_ptr);
               assert( tptr - f_ratio_ptr < 7);
               strlcpy( f_ratio, f_ratio_ptr, tptr - f_ratio_ptr + 1);
               tptr++;
               }
            }
         if( ofile)
            {
            int len = 0;
            const char *detector = "UNK";   /* default assumption */

            while( tptr[len] && tptr[len + 1] != '+')
               len++;
            fprintf( ofile, "        <aperture>%s</aperture>\n", aperture);
            if( *f_ratio)
               fprintf( ofile, "        <fRatio>%s</fRatio>\n", f_ratio);
            if( tptr[len] && tptr[len + 1] == '+')
               {
               if( strstr( tptr + len, "CCD"))
                  detector = "CCD";
               else if( strstr( tptr + len, "CMO"))
                  detector = "CMO";
               }
            fprintf( ofile, "        <detector>%s</detector>\n", detector);
            if( !memcmp( tptr, "coronagraph", 11))
               {
               fprintf( ofile, "        <name>%s</name>\n", tptr + 12);
               len = 11;
               }
            fprintf( ofile, "        <design>%.*s</design>\n", len, tptr);
            }
         n_found = 1;
         }
      return( n_found);
      }

   while( *line > ' ')
      line++;
   while( *line >= ' ')
      {
      size_t len = 0;

      while( *line == ' ')
         line++;
      while( line[len] != ',' && line[len] >= ' ')
         len++;
      if( len > 3 && line[1] == '.' && line[2] == ' ')
         {           /* very incomplete verification that it's really a name */
         if( ofile)
            fprintf( ofile, "        <name>%.*s</name>\n", (int)len, line);
         n_found++;
         }
      line += len;
      if( *line == ',')
         line++;
      }
   return( n_found);
}

static int output_names( FILE *ofile, const OBSERVE FAR *obs, const char *target)
{
   int i, n_found = 0;
   const size_t tlen = strlen( target);
   const char *tptr;
   const char **obs_details = obs->obs_details;
   const bool only_one_line = (!strcmp( target, "TEL ") || !ofile);

   if( obs_details)
      for( i = 0; (tptr = obs_details[i]) != NULL; i++)
         if( !strncmp( tptr, target, tlen))
            {
            n_found += dump_one_line_of_names( ofile, tptr);
            if( only_one_line)
               break;
            }

   if( !n_found)
      {
      FILE *ifile = fopen_ext( "details.txt", "fcrb");
      char buff[200];
      int got_it = 0;

      while( !got_it && fgets( buff, sizeof( buff), ifile))
         got_it = (!memcmp( buff, "COD ", 4) && !memcmp( buff + 4, obs->mpc_code, 3)
                  && (buff[7] != ' ' || buff[8] == program_code( obs)));
      if( got_it)
         while( fgets_trimmed( buff, sizeof( buff), ifile) && memcmp( buff, "COD ", 4))
            if( !strncmp( buff, target, tlen))
               {
               n_found += dump_one_line_of_names( ofile, buff);
               if( only_one_line)
                  break;
               }

      if( !n_found)     /* look for 'default' details */
         {
         while( fgets( buff, sizeof( buff), ifile))
            if( !memcmp( buff, "COD Def", 7))
               {
               while( fgets_trimmed( buff, sizeof( buff), ifile))
                  if( !strncmp( buff, target, tlen))
                     n_found += dump_one_line_of_names( ofile, buff);
               }
         }
      fclose( ifile);
      }
   return( n_found);
}

/* Outputs _only_ those observations from the station and program code
specified by the first observation.      */

static void create_ades_file_for_one_code( FILE *ofile,
                   const OBSERVE FAR *obs, int n_obs)
{
   char buff[200];
   const char *code = obs->mpc_code;
   const char progcode = program_code( obs);
   const bool is_sungrazer = (strstr( "249 C49 C50", code) != NULL);

   fprintf( ofile, "  <obsBlock>\n");
   fprintf( ofile, "    <obsContext>\n");
   fprintf( ofile, "      <observatory>\n");
   fprintf( ofile, "        <mpcCode>%s</mpcCode>\n", obs->mpc_code);
   fprintf( ofile, "      </observatory>\n");
   if( output_names( NULL, obs, "CON ") || is_sungrazer)
      {
      fprintf( ofile, "      <submitter>\n");
               /* Workaround for sungrazers without contact info */
      if( !output_names( ofile, obs, "CON ") && is_sungrazer)
         fprintf( ofile, "        <name>W. Gray</name>\n");
      fprintf( ofile, "      </submitter>\n");
      }
   if( output_names( NULL, obs, "OBS "))
      {
      fprintf( ofile, "      <observers>\n");
      output_names( ofile, obs, "OBS ");
      fprintf( ofile, "      </observers>\n");
      }
   if( output_names( NULL, obs, "MEA "))
      {
      fprintf( ofile, "      <measurers>\n");
      output_names( ofile, obs, "MEA ");
      fprintf( ofile, "      </measurers>\n");
      }
   if( output_names( NULL, obs, "TEL "))
      {
      fprintf( ofile, "      <telescope>\n");
      output_names( ofile, obs, "TEL ");
      fprintf( ofile, "      </telescope>\n");
      }
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
                        i + 1, (is_au ? 9 : 4), posn[i] * (is_au ? 1. : AU_IN_KM), i + 1);
            if( !(obs->flags & OBS_NO_VELOCITY))
               {
               double vel[3];

               compute_observer_vel( obs->jd, 3, 0., 0., 0., vel);
               for( i = 0; i < 3; i++)
                  vel[i] = obs->obs_vel[i] - vel[i]; /* make a geocentric velocity vector */
               ecliptic_to_equatorial( vel);
               for( i = 0; i < 3; i++)
                  fprintf( ofile, "        <vel%d>%.8f</vel%d>\n", i + 1,
                        vel[i] * (is_au ? 1. : AU_IN_KM / seconds_per_day), i + 1);
               }
            }
         else if( obs->note2 == 'V')
            {
            int i;

            fprintf( ofile, "        <sys>WGS84</sys>\n");
            fprintf( ofile, "        <ctr>399</ctr>\n");
            for( i = 0; i < 3; i++)
               {
               memcpy( buff, obs->second_line + i * 11 + 34, 10);
               buff[i == 2 ? 5 : 10] = '\0';
               text_search_and_replace( buff, " ", "");
               fprintf( ofile, "        <pos%d>%s</pos%d>\n",
                        i + 1, buff, i + 1);
               }
            }
//       if( progcode != ' ')
//          fprintf( ofile, "        <prog>%c</prog>\n", progcode);
         fprintf( ofile, "        <obsTime>%s</obsTime>\n",
                  iso_time( buff, utc_from_td( obs->jd, NULL), 3));
         fprintf( ofile, "        <ra>%.9f</ra>\n",
                                    original_observed_ra( obs) * 180. / PI);
         fprintf( ofile, "        <dec>%.9f</dec>\n",
                                    original_observed_dec( obs) * 180. / PI);
         fprintf( ofile, "        <rmsRA>%.3f</rmsRA>\n", obs->posn_sigma_1);
         fprintf( ofile, "        <rmsDec>%.3f</rmsDec>\n", obs->posn_sigma_2);
         fprintf( ofile, "        <rmsCorr>%.4f</rmsCorr>\n", correlation);
         if( !catalogue)
            catalogue = "UNK";
         strcpy( buff, catalogue);
         text_search_and_replace( buff, "-", "");
         text_search_and_replace( buff, " ", "");
         fprintf( ofile, "        <astCat>%s</astCat>\n", buff);
         if( obs->obs_mag != BLANK_MAG)
            {
            char mag_band = obs->mag_band;

            fprintf( ofile, "        <mag>%.*f</mag>\n",
                                    obs->mag_precision, obs->obs_mag);
            if( is_sungrazer && mag_band == ' ')
               mag_band = 'V';
            fprintf( ofile, "        <band>%c</band>\n", mag_band);
            }
         if( is_sungrazer)   /* for these,  the 'note' is always a 'program code' */
            fprintf( ofile, "        <notes>%c</notes>\n", progcode);
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
   fprintf( ofile, "<ades version=\"2022\">\n");
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
