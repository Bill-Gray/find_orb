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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include "watdefs.h"
#include "comets.h"
#include "afuncs.h"
#include "constant.h"

/* BC-405 gives orbital elements for 300 large asteroids at 40-day intervals,
running from JD 2378495.0 = 1799 Dec 30.5 to JD 2524615.0 = 2200 Jan 22.5.
That's data for 3654 epochs.  For each,  elements are given...

a = semimajor axis (in AU)
e = eccentricity
i = inclination (all angles in radians)
omega = arg perihelion
Omega = ascending node
M = mean anomaly

   ...for a total of 6 * 300 * 3654 = 6577200 lines.  The data are given in
ASCII,  in non-fixed-length format so you can't just hop into the file and
find the data for asteroid 141 at the 1878th interval.  So our first move
is to convert the data to binary,  8-byte double-precision floats,  resulting
in a file of 6 * 300 * 3654 = 52 617 600 bytes.  The input file is called
'asteroid_ephemeris.txt';  we create 'bc405.dat'.    */

const char *get_find_orb_text( const int index);      /* elem_out.cpp */
int detect_perturbers( const double jd, const double * __restrict xyz,
                       double *accel);          /* bc405.cpp */
double *get_asteroid_mass( const int astnum);   /* bc405.cpp */
int generic_message_box( const char *message, const char *box_type);
int asteroid_position_raw( const int astnum, const double jd,
                              double *posn, double *vel);      /* bc405.cpp */
int planet_posn( const int planet_no, const double jd, double *vect_2000);
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */

#define BC405_INVALID_CHUNK            (-1)
#define NO_BC405_FILE                  (-2)
#define BC405_COULDNT_OPEN_BINARY_FILE (-3)
#define MAX_BC405_N_ASTEROIDS    300

static int n_bc405_chunks = 3654;
static int bc405_n_asteroids = MAX_BC405_N_ASTEROIDS;
static double bc405_start_jd = 2378495.;
static double bc405_chunk_time = 40.;

static FILE *open_bc405_file( const bool shutting_down)
{
   const char *data_file_name = "bc405.dat";
   static FILE *ifile;
   static int failure_detected = 0;

   if( shutting_down)
      {
      if( ifile)
         fclose( ifile);
      ifile = NULL;
      return( NULL);
      }
   if( failure_detected)
      return( NULL);
   if( !ifile)
      ifile = fopen_ext( data_file_name, "crb");
   if( !ifile)
      {
      ifile = fopen_ext( "asteroid_ephemeris.txt", "crb");
      if( !ifile)       /* file name may be specified in environ.dat */
         {
         const char *ast_ephem_filename
                        = get_environment_ptr( "BC405_FILENAME");

         if( *ast_ephem_filename)
            ifile = fopen( ast_ephem_filename, "rb");
         }
      if( !ifile)       /* maybe the sub-ephemeris is available? */
         {
         ifile = fopen_ext( "bc405sub.dat", "crb");
         if( ifile)
            {
            int32_t temp;
            size_t count;

            fseek( ifile, -16, SEEK_END);
            count = fread( &bc405_start_jd, sizeof( double), 1, ifile);
            assert( count);
            count = fread( &temp, sizeof( int32_t), 1, ifile);
            assert( count);
            bc405_n_asteroids = (int)temp;
            count = fread( &temp, sizeof( int32_t), 1, ifile);
            assert( count);
            n_bc405_chunks = (int)temp;
            return( ifile);
            }
         }
      else
         {
         char buff[100];
         FILE *ofile = fopen_ext( data_file_name, "cwb");

         if( !ofile)
            return( NULL);
         while( fgets( buff, sizeof( buff), ifile))
            {
            const double z = atof( buff);

            fwrite( &z, 1, sizeof( double), ofile);
            }
         fclose( ifile);
         fclose( ofile);
         ifile = fopen_ext( data_file_name, "crb");
         assert( ifile != NULL);
         }
      }
   if( !ifile)             /* no asteroid ephems;  show err msg */
      {                    /* (see efindorb.txt)                */
      failure_detected = 1;
      generic_message_box( get_find_orb_text( 2021), "o");
      }
   return( ifile);
}

static void grab_elems( ELEMENTS *elems, FILE *fp, const int chunk_number)
{
   double array[6];

   memset( elems, 0, sizeof( ELEMENTS));
   if( fread( array, sizeof( double), 6, fp) == 6)
      {
      elems->major_axis = array[0];
      elems->ecc = array[1];
      elems->q = array[0] * (1. - array[1]);
      elems->incl = array[2];
      elems->arg_per = array[3];
      elems->asc_node = array[4];
      elems->mean_anomaly = array[5];
      elems->epoch = bc405_start_jd + (double)chunk_number * bc405_chunk_time;
      derive_quantities( elems, SOLAR_GM);
      elems->perih_time = elems->epoch - elems->mean_anomaly * elems->t0;
      }
   assert( array[0] > 1. && array[0] < 6.);
   assert( array[1] > 0. && array[1] < 1.);
}

#define N_CACHED_ELEMS 30

static void grab_cached_elems( ELEMENTS *elems, const int chunk_number,
                                 const int asteroid_number)
{
   static int astnums[N_CACHED_ELEMS], chunk_num[N_CACHED_ELEMS];
   static ELEMENTS *cache = NULL;
   FILE *fp;
   int i;

   if( !elems)
      {
      if( cache)
         free( cache);
      cache = NULL;
      for( i = 0; i < N_CACHED_ELEMS; i++)
         astnums[i] = chunk_num[i] = -1;
      return;
      }
   fp = open_bc405_file( false);
   assert( fp);
   if( !cache)
      cache = (ELEMENTS *)calloc( N_CACHED_ELEMS, sizeof( ELEMENTS));
   assert( cache);
   if( !cache)
      return;
   for( i = 0; i < N_CACHED_ELEMS && (asteroid_number != astnums[i]
               || chunk_number != chunk_num[i]); i++)
      ;
   if( i == N_CACHED_ELEMS)      /* looked in entire cache,  didn't find it */
      {
      i--;
      fseek( fp, (chunk_number * bc405_n_asteroids + asteroid_number)
                            * 6 * sizeof( double), SEEK_SET);
      grab_elems( &cache[i], fp, chunk_number);
      }
   *elems = cache[i];
   while( i > 0)
      {
      cache[i] = cache[i - 1];
      astnums[i] = astnums[i - 1];
      chunk_num[i] = chunk_num[i - 1];
      i--;
      }
   cache[0] = *elems;
   astnums[0] = asteroid_number;
   chunk_num[0] = chunk_number;
}

/* The idea of using 300 asteroid perturbers routinely is computationally
daunting.  After some testing and research,  I ended up with the following
scheme,  which appears to do a decent job of ensuring that asteroids
with significant perturbing effects are included,  whilst others don't
end up slowing integration by a noticeable amount.

   We only compute asteroid perturbations if the asteroid in question is
within a certain distance of the object being integrated.  The asteroid
orbital elements in BC405 are given at 40-day intervals,  so it makes
sense to pre-compute,  for each asteroid,  its range in x, y,  and z (in
heliocentric ecliptical coordinates) over a 40-day interval.

   Suppose that a given asteroid,  over that 40-day interval,  ranges in
x from -1.34 to -0.78 AU.  And suppose we've determined,  based on the
asteroid's mass,  that its perturbations should be included whenever the
object is within 0.27 AU.  In that case,  if the object we are integrating
has x > -0.51 or x < -1.61 AU,  we can ignore that asteroid.  Even if it
is within the correct range for x,  we still may find cause to ignore it
for y or z.  Furthermore,  if x, y, and z are scaled suitably,  we can
do all of these tests with integer arithmetic.  This will usually result
in the number of potentially significant perturbers being dropped from 300
to zero,  or sometimes a few.

   A few wrinkles to be considered : the scaled-integer xyz ranges are
stored in the file 'bc405pre.dat'.  This file is created if it does not
exist,  and filled with zeroes for all the 40-day ranges for all 300
asteroids.  Then,  as ranges are needed,  they are computed and written
back out to the file.  This wasn't a big deal to do,  and spared the
need to distribute a rather large data file.

   Also:  the "effective range" for perturbations,  0.27 in the above
example,  is scaled by the inverse of the _square root_ of the mass of
the asteroid.  Initially,  it did seem as if scaling by the inverse of
the mass of the asteroid would be correct.  However,  in a resonant
situation,  (GM)^-0.5 is actually right.  Currently,  it's set up so
that (2) Pallas has an effective range of 10 AU,  and all others are
scaled to this.

   And a final note : if you look at the comments for ASTEROID_PERT_LIST
in 'environ.def',  you'll see that one can override Find_Orb's judgment
and say "just use asteroids A, B, C, ...".  This helps if you're trying
to replicate results from a particular integrator.  (JPL,  for example,
has a "Big-16" list of asteroid perturbers they use frequently.)
Normally,  Find_Orb's method will be more accurate,  but it'll be
different,  making comparisons difficult.   */

static FILE *get_precomputed_data_fp( void)
{
   const char *data_file_name = "bc405pre.dat";
   FILE *ifile = fopen_ext( data_file_name, "cr+b");

   if( !ifile)             /* Create output binary file if not found, */
      {                    /* seeding with zeroes :                   */
      FILE *ofile = fopen_ext( data_file_name, "cw+b");
      int16_t posns[MAX_BC405_N_ASTEROIDS];
      int i;

      assert( ofile);
      for( i = 0; i < bc405_n_asteroids; i++)
         posns[i] = 0;
      for( i = 0; i < 3 * n_bc405_chunks; i++)
         {
         size_t count = fwrite( posns, sizeof( int16_t), bc405_n_asteroids, ofile);

         assert( count == (size_t)bc405_n_asteroids);
         }
      fclose( ofile);
      ifile = fopen_ext( data_file_name, "cr+b");
      assert( ifile);
      }
   return( ifile);
}

      /* Our pre-computed data for asteroid location limits is stored */
      /* in integer units,  with 1000 such units equalling one AU:    */
const double integer_scale = 1000.;

static int find_and_set_precomputed_data( FILE *precomputed_fp,
                         const int bc_chunk, int16_t *__restrict posns)
{
   const int chunk_size = bc405_n_asteroids * 3;
   FILE *bc405_elem_file = open_bc405_file( false);
   int i;

   assert( bc405_elem_file);
   if( bc_chunk < 0 || bc_chunk >= n_bc405_chunks)
      printf( "%d %d\n", bc_chunk, n_bc405_chunks);
   assert( bc_chunk >= 0 && bc_chunk < n_bc405_chunks);
   fseek( precomputed_fp, (size_t)( bc_chunk * chunk_size) * sizeof( int16_t), SEEK_SET);
   if( !fread( posns, chunk_size, sizeof( int16_t), precomputed_fp))
      return( -1);
   if( !posns[0] && !posns[1] && !posns[2])     /* still all zeroes; gotta compute */
      {
      int16_t * __restrict pptr = posns;

      fseek( bc405_elem_file, (size_t)bc_chunk * bc405_n_asteroids
                    * 6 * sizeof( double), SEEK_SET);
      for( i = 0; i < bc405_n_asteroids; i++)
         {
         ELEMENTS elems;
         double posn[4];

         grab_elems( &elems, bc405_elem_file, bc_chunk);
         comet_posn( &elems, elems.epoch - bc405_chunk_time / 2., posn);
         *pptr++ = (int16_t)( posn[0] * integer_scale);
         *pptr++ = (int16_t)( posn[1] * integer_scale);
         *pptr++ = (int16_t)( posn[2] * integer_scale);
         }
      fseek( precomputed_fp, -chunk_size * sizeof( int16_t), SEEK_CUR);
      fwrite( posns, chunk_size, sizeof( int16_t), precomputed_fp);
      }
   return( 0);
}

static int asteroid_numbers[MAX_BC405_N_ASTEROIDS];
int excluded_asteroid_number = 0;
static double *masses;

static double *load_asteroid_masses( void)
{
   FILE *ifile = fopen_ext( "mu1.txt", "fcrb");
   double *rval = NULL;

   assert( ifile);
   if( ifile)
      {
      int i = 0, ast_number;
      char buff[180];
      double mass;

      rval = (double *)calloc( bc405_n_asteroids, sizeof( double));
      assert( rval);
      while( i < bc405_n_asteroids && fgets( buff, sizeof( buff), ifile))
         if( *buff != ';'
               && sscanf( buff, "%d%lf", &ast_number, &mass) == 2)
            {
            asteroid_numbers[i] = ast_number;
            if( rval)
               rval[i++] = mass;
            }
      fclose( ifile);
/*    rval[0] = 4.747105847157599504245142421641e-10;  (solar masses) */
      }  /* Ceres: special 63 km^3/s^2 value:  not now in use */
         /* since found to be 62.68 km^3/s^2 by Dawn = 4.7229935634e-10 */
   assert( rval);
   return( rval);
}

int asteroid_position_raw( const int astnum, const double jd,
                              double *posn, double *vel)
{
   ELEMENTS elem;
   int chunk;

   open_bc405_file( false);
   chunk = (int)( (jd - bc405_start_jd) / bc405_chunk_time + .5);
   if( chunk < 0)
      chunk = 0;
   else if( chunk >= n_bc405_chunks - 1)
      chunk = n_bc405_chunks - 1;
   grab_cached_elems( &elem, chunk, astnum);
   comet_posn_and_vel( &elem, jd, posn, vel);
   return( 0);
}

double *get_asteroid_mass( const int astnum)
{
   int i;
   double *rval = NULL;

   if( !masses)
      masses = load_asteroid_masses( );
   for( i = 0; !rval && i < bc405_n_asteroids; i++)
      if( asteroid_numbers[i] == astnum)
         rval = &masses[i];
   return( rval);
}

int detect_perturbers( const double jd, const double * __restrict xyz,
                       double *accel)
{
   static int curr_chunk = -1;
   static int16_t posns0[MAX_BC405_N_ASTEROIDS * 3];
   static int16_t posns1[MAX_BC405_N_ASTEROIDS * 3];
   int16_t ixyz[3];
   static FILE *precomputed_fp;
   static bool bc405_available = true;
   static int n_asteroids_to_use = 0;
   double thresh;
   int i, load_posn0 = 0, load_posn1 = 0, chunk, n_fixed = 0;
   const char *fixed_perturber_list;
   int fixed_perturbers[MAX_BC405_N_ASTEROIDS];

   if( !bc405_available)
      return( NO_BC405_FILE);
   if( !xyz)              /* freeing memory, closing cached file pointers */
      {
      if( precomputed_fp)
         fclose( precomputed_fp);
      if( masses)
         free( masses);
      precomputed_fp = NULL;
      masses = NULL;
      open_bc405_file( true);
      grab_cached_elems( NULL, 0, 0);
      return( 0);
      }
   thresh = atof( get_environment_ptr( "ASTEROID_THRESH"));
   fixed_perturber_list = get_environment_ptr( "ASTEROID_PERT_LIST");
   if( !masses)
      masses = load_asteroid_masses( );
   if( !masses)
      {
      bc405_available = false;
      return( NO_BC405_FILE);
      }
   if( !open_bc405_file( false))
      {
      bc405_available = false;
      return( NO_BC405_FILE);
      }

   if( !precomputed_fp)
      precomputed_fp = get_precomputed_data_fp( );
   if( !precomputed_fp)
      {
      bc405_available = false;
      return( NO_BC405_FILE);
      }

   while( *fixed_perturber_list)  /* see ASTEROID_PERT_LIST comments in */
      {                           /* 'environ.def' for info on this */
      fixed_perturbers[n_fixed++] = atoi( fixed_perturber_list);
      while( *fixed_perturber_list && *fixed_perturber_list != ',')
         fixed_perturber_list++;
      if( *fixed_perturber_list == ',')
         fixed_perturber_list++;
      }

   chunk = (int)( (jd - bc405_start_jd) / bc405_chunk_time + .5);
   if( chunk < 0)
      chunk = 0;
                  /* Since we use both the 'current' _and_ the 'next' */
                  /* chunks,  we can't go past n_bc405_chunks - 2 :   */
   if( chunk > n_bc405_chunks - 2)
      chunk = n_bc405_chunks - 2;
   if( chunk == curr_chunk - 1)
      {
      memcpy( posns1, posns0, bc405_n_asteroids * 3 * sizeof( int16_t));
      load_posn0 = 1;
      }
   else if( chunk == curr_chunk + 1)
      {
      memcpy( posns0, posns1, bc405_n_asteroids * 3 * sizeof( int16_t));
      load_posn1 = 1;
      }
   else if( chunk != curr_chunk)
      load_posn0 = load_posn1 = 1;
   if( load_posn0)
      find_and_set_precomputed_data( precomputed_fp, chunk, posns0);
   if( load_posn1)
      find_and_set_precomputed_data( precomputed_fp, chunk + 1, posns1);
   curr_chunk = chunk;
   if( !n_asteroids_to_use)
      n_asteroids_to_use = atoi( get_environment_ptr( "BC405_ASTEROIDS"));
   if( !n_asteroids_to_use)
      n_asteroids_to_use = bc405_n_asteroids;
   if( !thresh)
      thresh = 10.;                              /* Pallas extends 10 AU;  all others */
   thresh *= integer_scale / sqrt( masses[1]);   /* scaled by sqrt of their masses    */
   for( i = 0; i < 3; i++)
      ixyz[i] = (int16_t)( integer_scale * xyz[i]);
   for( i = 0; i < n_asteroids_to_use; i++)
      {
      const int16_t *p0 = posns0 + i * 3;
      const int16_t *p1 = posns1 + i * 3;
      int j, possible_perturber = 1, fixed_perturber = 0;
      const double dthresh = thresh * sqrt( masses[i]) + .1 * integer_scale;
      const int16_t ithresh = (int16_t)( dthresh > 27000. ? 27000 : dthresh);

      for( j = 0; j < n_fixed; j++)
         if( asteroid_numbers[i] == fixed_perturbers[j])
            fixed_perturber = 1;
      if( *p0 > *p1)
         possible_perturber = (ixyz[0] + ithresh > *p1 && ixyz[0] - ithresh < *p0);
      else
         possible_perturber = (ixyz[0] + ithresh > *p0 && ixyz[0] - ithresh < *p1);
      if( asteroid_numbers[i] == excluded_asteroid_number)
         possible_perturber = fixed_perturber = 0;  /* don't let an asteroid perturb itself! */

      if( possible_perturber || fixed_perturber)
         {
         p0++;
         p1++;
         if( *p0 > *p1)
            possible_perturber = (ixyz[1] + ithresh > *p1 && ixyz[1] - ithresh < *p0);
         else
            possible_perturber = (ixyz[1] + ithresh > *p0 && ixyz[1] - ithresh < *p1);
         if( possible_perturber || fixed_perturber)
            {
            p0++;
            p1++;
            if( *p0 > *p1)
               possible_perturber = (ixyz[2] + ithresh > *p1 && ixyz[2] - ithresh < *p0);
            else
               possible_perturber = (ixyz[2] + ithresh > *p0 && ixyz[2] - ithresh < *p1);
            if( n_fixed)          /* fixed perturbers set;  only consider them */
               possible_perturber = 0;
            if( possible_perturber || fixed_perturber)
               {
               double asteroid_loc[4], dist2 = 0., delta[3], sun_dist2 = 0.;
               double factor1, factor2;

               planet_posn( i + 100, jd, asteroid_loc);
               for( j = 0; j < 3; j++)
                  {
                  delta[j] = asteroid_loc[j] - xyz[j];
                  sun_dist2 += asteroid_loc[j] * asteroid_loc[j];
                  dist2 += delta[j] * delta[j];
                  }
               factor1 = SOLAR_GM * masses[i] / (sun_dist2 * sqrt( sun_dist2));
               factor2 = SOLAR_GM * masses[i] / (dist2 * sqrt( dist2));
               for( j = 0; j < 3; j++)
                  {
                  accel[j + 3] += factor2 * delta[j];
                  accel[j + 3] -= factor1 * asteroid_loc[j];
                  }
#ifdef DEBUGGING_CODE
               if( dist2 < .05 * 0.05)
                  {
                  FILE *debug_file = fopen( "astpert.txt", "ab");

                  fprintf( debug_file, "%.5f: %3d, %f: mass %g\n",
                           JD_TO_YEAR( jd),
                           asteroid_numbers[i], sqrt( dist2) * AU_IN_KM, masses[i]);
                  fclose( debug_file);
                  }
#endif
               }
            }
         }
      }
   return( 0);
}

#ifdef TEST_CODE

/* When I first came up with the idea of filtering possible
asteroid perturbers by distance,  I wanted to be able to test
different points at different times to see how many objects would
be considered "close enough" to be possible perturbers.  I expected,
of course,  to see more perturbers for points inside the main
asteroid belt.  The following test code checked to see how many
potential perturbers would lie within 0.15 AU of a given point
at a given time.  */

static int inside_box( const int16_t * __restrict point,
                       const int16_t * __restrict p1, const int16_t * __restrict p2)
{
   int count = 3, rval = 1;
   const int16_t thresh = 150;        /* .15 AU */

   while( rval && count--)
      {
      if( *p2 > *p1)
         rval = (*point > *p1 - thresh && *point < *p2 + thresh);
      else
         rval = (*point > *p2 - thresh && *point < *p1 + thresh);
      p2++;
      p1++;
      point++;
      }
   return( rval);
}

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int main( const int argc, char **argv)
{
   const int asteroid_number = atoi( argv[1]);
   const double jd = atof( argv[2]);
   const int chunk_number = (int)( (jd - bc405_start_jd) / bc405_chunk_time + .5);
   FILE *fp = open_bc405_file( false);
   ELEMENTS elems;
   double posn[4];
   int16_t precomputed[bc405_n_asteroids * 3];

   fseek( fp, (chunk_number * bc405_n_asteroids + asteroid_number)
                       * 6 * sizeof( double), SEEK_SET);
   grab_elems( &elems, fp, chunk_number);
   printf( "Chunk %d; epoch %.1f\n", chunk_number, elems.epoch);
   printf( "Semimajor: %f\n", elems.major_axis);
   printf( "Ecc: %f\n", elems.ecc);
   printf( "Incl: %f\n", elems.incl * 180. / PI);
   printf( "Asc node: %f\n", elems.asc_node * 180. / PI);
   printf( "Arg per: %f\n", elems.arg_per * 180. / PI);
   printf( "Mean anom: %f\n", elems.mean_anomaly * 180. / PI);
   printf( "Tp = %f\n", elems.perih_time);

   comet_posn( &elems, jd, posn);
// ecliptic_to_equatorial( posn);
   printf( "x = %f; y = %f; z = %f\n", posn[0], posn[1], posn[2]);
   find_and_set_precomputed_data( chunk_number, precomputed);
   printf( "Precomputed: %d %d %d\n",
            (int)precomputed[asteroid_number * 3 + 0],
            (int)precomputed[asteroid_number * 3 + 1],
            (int)precomputed[asteroid_number * 3 + 2]);

   if( argc == 6)
      {
      int16_t prec2[bc405_n_asteroids * 3];
      int16_t iposn[3], *pptr1 = precomputed, *pptr2 = prec2;
      int n_found = 0, i;

      find_and_set_precomputed_data( chunk_number + 1, prec2);
      for( i = 0; i < 3; i++)
         iposn[i] = (int16_t)atoi( argv[i + 3]);
      for( i = 0; i < bc405_n_asteroids; i++, pptr1 += 3, pptr2 += 3)
         if( inside_box( iposn, pptr1, pptr2))
            n_found++;
      printf( "%d in box\n", n_found);
      }
   return( 0);
}
#endif
