/* pl_cache.cpp: computes and caches planetary positions

SEE PL_CACHE.TXT FOR A DISCUSSION OF WHAT THIS DOES.  It probably
won't make much sense to you if you don't.

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

/* #define TIMING_ON  */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#ifdef TIMING_ON
    #include <time.h>
#endif
#include "pl_cache.h"
#include "watdefs.h"
#include "stringex.h"
#include "lunar.h"
#include "afuncs.h"
#include "jpleph.h"

const char *get_find_orb_text( const int index);          /* elem_out.cpp */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
extern int debug_level;

int64_t planet_ns;
static void *jpl_eph = NULL;

#define J2000 2451545.0
#define J0 (J2000 - 2000. * 365.25)
#define JD_TO_YEAR( jd)  (((jd)-J2000) / 365.25 + 2000.)

FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
char *make_config_dir_name( char *oname, const char *iname);  /* miscell.cpp */
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
int generic_message_box( const char *message, const char *box_type);

int compute_rough_planet_loc( const double t_cen, const int planet_idx,
                                          double *vect);    /* sm_vsop.cpp */
int asteroid_position_raw( const int astnum, const double jd,
                              double *posn, double *vel);      /* bc405.cpp */
int64_t nanoseconds_since_1970( void);                      /* mpc_obs.c */
int format_jpl_ephemeris_info( char *buff);                 /* pl_cache.c */

static int planet_posn_raw( int planet_no, const double jd,
                            double *vect_2000)
{
   const int jpl_center = 11;         /* default to heliocentric */
   int rval = 0;
   static const char *jpl_filename = NULL;
   const int bc405_start = 100;
   const int calc_vel = (planet_no > PLANET_POSN_VELOCITY_OFFSET - 2);

   if( calc_vel)
      planet_no -= PLANET_POSN_VELOCITY_OFFSET;
   if( !planet_no)            /* the sun */
      {
      vect_2000[0] = vect_2000[1] = vect_2000[2] = 0.;
      if( !jd && jpl_eph)       /* return version data: */
         {
         vect_2000[0] = (double)jpl_get_long( jpl_eph, JPL_EPHEM_EPHEMERIS_VERSION);
         vect_2000[1] = jpl_get_double( jpl_eph, JPL_EPHEM_START_JD);
         vect_2000[2] = jpl_get_double( jpl_eph, JPL_EPHEM_END_JD);
         }
      return( 0);
      }

   if( planet_no >= bc405_start && planet_no < bc405_start + 300)
      {
      double temp_loc[4];

      rval = asteroid_position_raw( planet_no - bc405_start, jd,
               (calc_vel ? NULL : temp_loc),
               (calc_vel ? temp_loc : NULL));
      if( debug_level > 8)
         debug_printf( "JD %f, minor planet %d: (%f %f %f)\n",
                     jd, planet_no, temp_loc[0], temp_loc[1], temp_loc[2]);
      memcpy( vect_2000, temp_loc, 3 * sizeof( double));
      return( rval);
      }

   if( !jpl_filename)
      {
      FILE *ifile;

#if defined (_WIN32) || defined( __WATCOMC__)
      jpl_filename = get_environment_ptr( "JPL_FILENAME");
#else
      jpl_filename = get_environment_ptr( "LINUX_JPL_FILENAME");
#endif
      if( *jpl_filename)
         jpl_eph = jpl_init_ephemeris( jpl_filename, NULL, NULL);
      if( !jpl_eph)
         if( (ifile = fopen_ext( "jpl_eph.txt", "fcrb")) != NULL)
            {
            char buff[100];

            while( !jpl_eph && fgets_trimmed( buff, sizeof( buff), ifile))
               if( *buff && *buff != ';')
                  {
                  jpl_eph = jpl_init_ephemeris( buff, NULL, NULL);
                  if( !jpl_eph)
                     {
                     char tname[255];

                     make_config_dir_name( tname, buff);
                     jpl_eph = jpl_init_ephemeris( tname, NULL, NULL);
                     }
                  }
            if( debug_level)
               debug_printf( "Ephemeris file %s\n", buff);
            fclose( ifile);
            }
      if( debug_level && jpl_eph)
         {
         debug_printf( "\nEphemeris time span years %.3f to %.3f\n",
               (jpl_get_double( jpl_eph, JPL_EPHEM_START_JD) - J0) / 365.25,
               (jpl_get_double( jpl_eph, JPL_EPHEM_END_JD)   - J0) / 365.25);
         debug_printf( "Ephemeris version %ld\n", jpl_get_long( jpl_eph, JPL_EPHEM_EPHEMERIS_VERSION));
         debug_printf( "Kernel size %ld, record size %ld, swap_bytes %ld\n",
               jpl_get_long( jpl_eph, JPL_EPHEM_KERNEL_SIZE),
               jpl_get_long( jpl_eph, JPL_EPHEM_KERNEL_RECORD_SIZE),
               jpl_get_long( jpl_eph, JPL_EPHEM_KERNEL_SWAP_BYTES));
         debug_printf( "ncon = %ld AU=%f emrat = %f\n",
               jpl_get_long( jpl_eph, JPL_EPHEM_N_CONSTANTS),
               jpl_get_double( jpl_eph, JPL_EPHEM_AU_IN_KM),
               jpl_get_double( jpl_eph, JPL_EPHEM_EARTH_MOON_RATIO));
         }
      }

   if( jpl_eph)
      {
      double state[6];            /* DE gives both posn & velocity */
      int failure_code;

      if( planet_no < 0)          /* flag to unload everything */
         {
         jpl_close_ephemeris( jpl_eph);
         jpl_eph = NULL;
         jpl_filename = NULL;
         return( 0);
         }
      else if( planet_no == 10)
         failure_code = jpl_pleph( jpl_eph, jd, 10, 3, state, calc_vel);
      else
         failure_code = jpl_pleph( jpl_eph, jd,
              (planet_no == 3) ? 13 : planet_no, jpl_center, state, calc_vel);
      if( !failure_code)         /* we're done */
         {
         if( debug_level > 8)
            debug_printf( "JD %f, planet %d: (%f %f %f)\n",
                     jd, planet_no, state[0], state[1], state[2]);
         memcpy( vect_2000, state + calc_vel * 3, 3 * sizeof( double));
         equatorial_to_ecliptic( vect_2000);
         return( 0);
         }
      else
         if( debug_level)
            debug_printf( "Failed: JD %f, planet %d, code %d\n",
                           jd, planet_no, failure_code);
      }

   if( planet_no < 0)          /* flag to unload everything */
      rval = 0;                /* Currently,  nothing to unload */
   else
      {
      static bool error_message_shown = false;

      if( !error_message_shown)
         {
         error_message_shown = true;
         if( jpl_eph)
            {
            double jd_start, jd_end;
            char *buff = (char *)malloc( 10000);

            assert( buff);
            debug_printf( "Failure on planet %d, JD %f\n", planet_no, jd);
            get_jpl_ephemeris_info( NULL, &jd_start, &jd_end);
            snprintf( buff, 10000, get_find_orb_text( 2030), jpl_filename,
                        JD_TO_YEAR( jd_start),
                        JD_TO_YEAR( jd_end));
            generic_message_box( buff, "o");
            free( buff);
            }
         else
            generic_message_box( get_find_orb_text( 2029), "o");
         }
      if( planet_no > 0 && planet_no <= 10)
         {
         compute_rough_planet_loc( (jd - J2000) / 36525., planet_no, vect_2000);
         if( calc_vel)
            {
            const double delta_t = 10. / 1440.;    /* ten minute delta */
            double v[3];
            size_t i;

            compute_rough_planet_loc( (jd - J2000 + delta_t) / 36525., planet_no, v);
            for( i = 0; i < 3; i++)
               vect_2000[i] = (v[i] - vect_2000[i]) / delta_t;
            }
         rval = 0;
         }
      }

   return( rval);
}

#define POSN_CACHE struct posn_cache

POSN_CACHE
   {
   double jd;
   double vect[3];
   int planet_no;
   };

#define POSN_NODE struct posn_node

POSN_NODE
   {
   double min_jd;
   int used;
   POSN_CACHE *data;
   };

#define node_size                     1659
   /* When a node is 90% full,  it's time to split it. */
#define splitting_size    (node_size - (node_size / 10))
   /* If a node reaches splitting_size,  "spill over" to an adjacent  */
   /* node if it's less than half full : */
#define spillover_size    (node_size / 2)
int n_posns_cached = 0;

/* Hash the JD and planet number.  It seems a fair bit of time is
spent in this function,  so I spent a good bit of time trying to make
it as simple/fast as possible while still giving good distribution
(i.e.,  no more or not many more table collisions than would be
expected with a "perfectly randomizing" hash function.)  It helps
that node_size is not a power of two.

   Note that it works well _for the planet positions being hashed
here_.  Don't rely on it as a general-purpose hashing function!

   Updated 2015 Feb 12:  MSVC objected to simply setting

   const int32_t *dword_ptr = (int32_t *)&jd;

   and accessing the four-byte halves of jd directly,  so we're now
doing a totally pointless memcpy.         */

static inline int hash_function( const int planet_no, const double jd)
{
   int32_t dword_ptr[2];
   int rval;

   memcpy( dword_ptr, &jd, sizeof( double));
   rval = dword_ptr[0] ^ dword_ptr[1] ^ (planet_no << 8);

   rval &= 0x7fffffff;
   rval %= node_size;
   assert( rval >= 0);
   assert( rval < node_size);
   return( rval);
}

/* When a node gets to 'splitting_size' values,  we compact them
(removing unused entries),  then partition them in half.  We do
that partitioning,  at least for the nonce,  the lazy way : we
sort the entire array.  */

void shellsort_r( void *base, const size_t n_elements, const size_t elem_size,
         int (*compare)(const void *, const void *, void *), void *context);

static int compare_cached_posns( const void *a, const void *b, void *ignored_context)
{
   const POSN_CACHE *pa = (const POSN_CACHE *)a;
   const POSN_CACHE *pb = (const POSN_CACHE *)b;
   int rval;

   INTENTIONALLY_UNUSED_PARAMETER( ignored_context);
   if( pa->jd > pb->jd)
      rval = 1;
   else if( pa->jd < pb->jd)
      rval = -1;
   else
      rval = 0;
   return( rval);
}

static void collapse_and_partition( POSN_CACHE *ovals, const POSN_CACHE *ivals)
{
   int i, array_size;

   for( i = array_size = 0; i < node_size; i++)
      if( ivals[i].planet_no)
         ovals[array_size++] = ivals[i];
   assert( array_size == splitting_size);
   shellsort_r( ovals, array_size, sizeof( POSN_CACHE), compare_cached_posns, NULL);
}

/*   The following three long ints keep track of the number of searches
and probes done,  and the "worst-case" maximum number of probes required,
just to check that the hash function is truly random enough.      */

/* #define TEST_PLANET_CACHING_HASH_FUNCTION */

#ifdef TEST_PLANET_CACHING_HASH_FUNCTION
long total_n_searches = 0, total_n_probes = 0, max_probes_required = 0;
#endif

static int find_within_node( const int planet_no, const double jd, const POSN_CACHE *cache)
{
   int loc = hash_function( planet_no, jd);
   int n_probes = 1;

   while( cache[loc].planet_no)
      {
      if( cache[loc].planet_no == planet_no && cache[loc].jd == jd)
         break;
      n_probes++;
      loc = (loc + n_probes) % node_size;
      }
#if 0
   for( int i = 0; i < node_size; i++)
      if( i != loc)
         assert( cache[i].planet_no != planet_no || cache[i].jd != jd);
#endif
#ifdef TEST_PLANET_CACHING_HASH_FUNCTION
   if( max_probes_required < n_probes)
      max_probes_required = n_probes;
   total_n_searches++;
   total_n_probes += n_probes;
#endif
   return( loc);
}

/* At one time or another,  I've had concerns that the scheme for caching
planetary positions (essentially a single-level hashed B-tree;  see
'pl_cache.txt' for details) was buggy.  The check_integrity() function does
some basic consistency tests (are all the nodes within the time range
they're supposed to be?  is the number of used nodes what it's supposed to
be?)  If such concerns recur,  we can #define CHECK_CACHING_INTEGRITY
and see if it finds anything we need to worry about.

#define CHECK_CACHING_INTEGRITY
*/

#ifdef CHECK_CACHING_INTEGRITY
static void check_integrity( const POSN_NODE *nodes, const int n_nodes)
{
   int i;

   for( i = 0; i < n_nodes; i++, nodes++)
      {
      int j, n_used = 0;
      const double max_jd = (i < n_nodes - 1 ? nodes[1].min_jd : 1e+20);
      const POSN_CACHE *cptr = nodes->data;

      for( j = 0; j < node_size; j++, cptr++)
         if( cptr->planet_no)
            {
            assert( cptr->jd >= nodes->min_jd);
            assert( cptr->jd <= max_jd);
            n_used++;
            }
      assert( n_used == nodes->used);
      }
}
#endif         /* #ifdef CHECK_CACHING_INTEGRITY */

/* Computing planetary positions is somewhat expensive if we're using
JPL ephemerides,  and _very_ expensive if we aren't (the PS-1996 method
is used,  which involves lots of trig series).  And frequently,  we'll be
requesting the same data over and over (for example,  if we're integrating
over a particular time span repeatedly).  So it makes sense to cache the
planetary positions.

   Below,  this is done with a hash table in a very standard sort of way.
The planet_no and jd are hashed, we look in the 'cache' table,  we do a
quadratic search if there's an hash collision.  If we find the data,  we
return it.  If we don't find it,  we call the planet_posn_raw( ) (uncached)
function, and add the result to the cache, and _then_ return it.

   When the table is more than 80% full,  we double the table size,  dump
everything computed to date,  and start from scratch.  This is admittedly
mildly wasteful,  but I don't think the performance benefit of expanding
the cache and adding everything we've got back in would be worthwhile.
(It wouldn't be hard to do,  though.)

   If the cache gets above some limit (currently set to a million cached
positions,  or about 36 MBytes),  we stop growing the cache.  So if you
had a _really_ long integration,  the cache gets dumped,  rebuilt to the
same size,  dumped again,  built to the same size,  etc.

   Previously,  the data was stored using a balanced tree.  I don't know
what possessed me to do something that dumb.  (At the very least,  had
I keyed the tree using the above hash function,  entries to the tree would
have been nearly random,  and a plain old unbalanced tree would have worked
Just Fine.)    */

#define MAX_N_NODES 10000

int planet_posn( const int planet_no, const double jd, double *vect_2000)
{
   static POSN_NODE *nodes = NULL;
   static int n_nodes = 0, n_nodes_alloced = 0, curr_node = 0;
   int loc, rval = 0;
#ifdef TIMING_ON
   int64_t t_start;
#endif
   POSN_CACHE *cache;

   assert( fabs( jd) < 1e+9);
   if( !planet_no)            /* the sun */
      {
      vect_2000[0] = vect_2000[1] = vect_2000[2] = 0.;
      return( 0);
      }

   if( planet_no < 0 || n_nodes >= MAX_N_NODES)
      {                                  /* flag to unload everything */
      int i;

      for( i = 0; i < n_nodes; i++)
         if( nodes[i].data)
            free( nodes[i].data);
      if( nodes)
         free( nodes);
      nodes = NULL;
      n_posns_cached = 0;
      n_nodes = n_nodes_alloced = curr_node = 0;
      }

   if( planet_no < 0)
      {
      planet_posn_raw( -1, 0., NULL);
      return( 0);
      }

   if( (planet_no % PLANET_POSN_VELOCITY_OFFSET) == PLANET_POSN_EARTH
             || (planet_no % PLANET_POSN_VELOCITY_OFFSET) == PLANET_POSN_MOON)
      {
      double moon_loc[3];
      const int vel_offset = (planet_no > PLANET_POSN_VELOCITY_OFFSET ?
                  PLANET_POSN_VELOCITY_OFFSET : 0);

      rval = planet_posn( 3 + vel_offset, jd, vect_2000);   /* first,  get Earth-Moon */
      if( !rval)                               /* barycenter posn,  then */
         rval = planet_posn( 10 + vel_offset, jd, moon_loc);    /* lunar offset vect  */
      if( !rval)
         {
         size_t i;
         const double EARTH_MOON_BARYCENTER_FACTOR = 82.300679;
         const double factor = (planet_no % PLANET_POSN_VELOCITY_OFFSET == PLANET_POSN_EARTH ?
                     -1. / EARTH_MOON_BARYCENTER_FACTOR :
                 1. - 1. / EARTH_MOON_BARYCENTER_FACTOR);

         for( i = 0; i < 3; i++)
            vect_2000[i] += moon_loc[i] * factor;
         }
      return( rval);
      }

   if( !nodes || n_nodes == n_nodes_alloced - 1)
      {
      const unsigned new_n_alloced = 100 + 3 * n_nodes_alloced / 2;

      nodes = (POSN_NODE *)realloc( nodes, new_n_alloced * sizeof( POSN_NODE));
      assert( nodes);
      if( !n_nodes_alloced)      /* set up first node : */
         {
         n_nodes = 1;
         nodes[0].min_jd = -1e+10;
         nodes[0].used = 0;
         nodes[0].data = (POSN_CACHE *)calloc( node_size, sizeof( POSN_CACHE));
         }
      n_nodes_alloced = new_n_alloced;
      }

             /* Now,  find the right node in which to find/store this posn: */
   while( curr_node + 1 < n_nodes && nodes[curr_node + 1].min_jd <= jd)
      curr_node++;
   while( curr_node && jd < nodes[curr_node].min_jd)
      curr_node--;
   assert( jd >= nodes[curr_node].min_jd);
   assert( curr_node == n_nodes - 1 || jd < nodes[curr_node + 1].min_jd);

   cache = nodes[curr_node].data;
   loc = find_within_node( planet_no, jd, cache);

#ifdef TIMING_ON
   t_start = nanoseconds_since_1970( );
#endif

   if( !cache[loc].planet_no)
      {
      cache[loc].planet_no = planet_no;
      cache[loc].jd = jd;
      nodes[curr_node].used++;
      rval = planet_posn_raw( planet_no, jd, cache[loc].vect);
      n_posns_cached++;
      }
   else
      {
      assert( cache[loc].planet_no == planet_no);
      assert( cache[loc].jd == jd);
      memcpy( vect_2000, cache[loc].vect, 3 * sizeof( double));
      return( rval);
      }
#ifdef TIMING_ON
   planet_ns += nanoseconds_since_1970( ) - t_start;
#endif
   memcpy( vect_2000, cache[loc].vect, 3 * sizeof( double));
   assert( nodes[curr_node].used <= splitting_size);
#ifdef CHECK_CACHING_INTEGRITY
   if( n_posns_cached % 10000 == 0)
      check_integrity( nodes, n_nodes);
#endif
   if( nodes[curr_node].used == splitting_size)
      {
      POSN_CACHE *tcache = (POSN_CACHE *)calloc( node_size, sizeof( POSN_CACHE));
      int i, size1 = splitting_size / 2;
      const int size_left = (curr_node ? nodes[curr_node - 1].used : node_size);
      const int size_right = (curr_node < n_nodes - 1 ? nodes[curr_node + 1].used : node_size);

#ifdef CHECK_CACHING_INTEGRITY
      check_integrity( nodes, n_nodes);
#endif
      collapse_and_partition( tcache, cache);
      memset( cache, 0, node_size * sizeof( POSN_CACHE));
      nodes[curr_node].used = 0;
      if( size_left < size_right && size_left < spillover_size)
         {        /* "spill over" to left */
         curr_node--;
         size1 -= size_left / 2;
         if( debug_level > 5)
            debug_printf( "Spilling %d to left: %d\n", curr_node, size1);
         }
      else if( size_right <= size_left && size_right < spillover_size)
         {        /* "spill over" to right */
         size1 += size_right / 2;
         if( debug_level > 5)
            debug_printf( "Spilling %d to right: %d\n", curr_node, size1);
         }
      else        /* create new node */
         {
         memmove( nodes + curr_node + 2, nodes + curr_node + 1,
                  (n_nodes - curr_node - 1) * sizeof( POSN_NODE));
         nodes[curr_node + 1].data = (POSN_CACHE *)calloc( node_size, sizeof( POSN_CACHE));
         nodes[curr_node + 1].used = 0;
         if( debug_level > 5)
            debug_printf( "Splitting node %d\n", curr_node);
         n_nodes++;
         }
      while( tcache[size1].jd == tcache[size1 - 1].jd)
         size1--;
      for( i = 0; i < splitting_size; i++)
         {
         const int n = curr_node + (i < size1 ? 0 : 1);
         POSN_CACHE *cptr = nodes[n].data;
         const int new_loc = find_within_node( tcache[i].planet_no,
                           tcache[i].jd, cptr);

         cptr[new_loc] = tcache[i];
         nodes[n].used++;
         }
      nodes[curr_node + 1].min_jd = tcache[size1].jd;
#ifdef CHECK_CACHING_INTEGRITY
      check_integrity( nodes, n_nodes);
#endif
      free( tcache);
      }
   return( rval);
}

      /* In the following,  we get the earth's position for a particular    */
      /* instant,  just to ensure that JPL ephemerides (if any) are loaded. */
      /* Then we call with planet = JD = 0,  which causes the info about    */
      /* the JPL ephemerides to be put into the 'state vector'.             */
int get_jpl_ephemeris_info( int *de_version, double *jd_start, double *jd_end)
{
   double vect_2000[3];

   planet_posn_raw( 3, J2000, vect_2000);
   planet_posn_raw( 0, 0., vect_2000);
   if( de_version)
      *de_version = (int)vect_2000[0];
   if( jd_start)
      *jd_start = (int)vect_2000[1];
   if( jd_end)
      *jd_end = (int)vect_2000[2];
   return( 0);
}

int format_jpl_ephemeris_info( char *buff)
{
   int de_version;
   double jd_start, jd_end;
   const size_t buff_size = 250;

   get_jpl_ephemeris_info( &de_version, &jd_start, &jd_end);
   if( !de_version && !jd_start && !jd_end)
      strlcpy_err( buff, get_find_orb_text( 2056), buff_size);
   else
      snprintf_err( buff, buff_size,
            "\nUsing %s; covers years %.1f to %.1f\n",
            jpl_get_ephem_name( jpl_eph),
            JD_TO_YEAR( jd_start), JD_TO_YEAR( jd_end));
   return( de_version);
}
