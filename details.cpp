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

/* details.cpp : code to store/access 'observational details' (header)
data for 80-column MPC-formatted astrometry,  as documented at

https://www.minorplanetcenter.net/iau/info/ObsDetails.html */

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "stackall.h"
#include "details.h"

typedef struct
{
   char **lines;
   size_t n_lines;
   bool observations_found;
} mpc_code_details_t;

typedef struct
{
   void *stack;
   mpc_code_details_t *code_details;
   int n_code_details, n_curr;
   int curr_idx[10];
} observation_details_t;

#define min_n_code_details 16
#define min_code_lines 16

int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;

void *init_observation_details( void)
{
   void *stack = create_stack( 2000);
   observation_details_t *rval = (observation_details_t *)
                     stack_calloc( stack, sizeof( observation_details_t));

   rval->stack = stack;
   rval->n_curr = 0;    /* i.e.,  "no COD line yet" */
   rval->code_details = (mpc_code_details_t *)stack_calloc( stack,
                              min_n_code_details * sizeof( mpc_code_details_t));
   return( rval);
}

      /* A useful piece of bit-twiddling. X & (X-1) clears the lowest */
      /* bit of X.  If X is a power of two,  you've only got one bit  */
      /* to clear,  and the result is zero. */
#define is_power_of_two( X)   (!((X) & ((X) - 1)))

/* Binary-search to find the desired MPC code in the 'code_details' array.
Or where that code ought to be put in,  if we're adding a new code (in
which case *cmp will be nonzero.) */

static int code_cmp( const char *a, const char *b)
{
   int rval;

   while( *a == *b && *a)
      {
      a++;
      b++;
      }
   if( *a <= ' ' && *b <= ' ')
      rval = 0;
   else
      rval = *a - *b;
   return( rval);
}

static int find_code_details( const observation_details_t *det,
                 const char *mpc_code, int *cmp)
{
   int loc = -1, stepsize = (int)0x8000, compare = 1, loc1;

   while( compare && (stepsize >>= 1))
      if( (loc1 = loc + stepsize) < det->n_code_details)
         {
         assert( det->code_details[loc1].lines);
         compare = code_cmp( mpc_code, det->code_details[loc1].lines[0] + 4);
         if( compare >= 0)
            loc = loc1;
         }
   *cmp = compare;
   return( loc);
}

const char **get_code_details( const void *obs_details, const char *mpc_code)
{
   const observation_details_t *det = (const observation_details_t *)obs_details;
   int idx, compare;

#ifdef DEBUG_PRINTFS
   printf( "We've got %d details :", det->n_code_details);
   for( idx = 0; idx < det->n_code_details; idx++)
      printf( " '%s'", det->code_details[idx].lines[0]);
   printf( " And that's it\n");
#endif
   idx = find_code_details( det, mpc_code, &compare);
   if( !compare)
      return( (const char **)det->code_details[idx].lines);
   else
      return( NULL);
}

static int reset_mpc_code( observation_details_t *det, const char *mpc_code)
{
   int compare, idx = find_code_details( det, mpc_code, &compare);

#ifdef DEBUG_PRINTFS
   printf( "Resetting MPC code: idx %d, compare %d, %d details\n",
               idx, compare, det->n_code_details);
#endif
   if( compare)   /* not an MPC code we've already encountered */
      {
      idx++;
      det->n_code_details++;
      if( is_power_of_two( det->n_code_details) && det->n_code_details >= min_n_code_details)
         {
         mpc_code_details_t *new_array = (mpc_code_details_t *)
                  stack_alloc( det->stack, 2 * det->n_code_details * sizeof( mpc_code_details_t));

         memcpy( new_array, det->code_details, det->n_code_details * sizeof( mpc_code_details_t));
         det->code_details = new_array;
         }
#ifdef DEBUG_PRINTFS
      printf( "Moving memory\n");
#endif
      memmove( det->code_details + idx + 1, det->code_details + idx,
                  (det->n_code_details - idx) * sizeof( mpc_code_details_t));
#ifdef DEBUG_PRINTFS
      printf( "Moved memory\n");
#endif
      }
   det->code_details[idx].observations_found = false;
   det->code_details[idx].n_lines = 0;
   det->code_details[idx].lines = (char **)stack_calloc( det->stack,
                     (min_code_lines + 1) * sizeof( char *));
   det->curr_idx[det->n_curr++] = idx;
#ifdef DEBUG_PRINTFS
   printf( "Done\n");
#endif
   return( 0);
}

static int probable_mpc_record( const char *buff)
{
   int i, rval = 1;     /* assume it's a valid record */

   if( buff[12] != ' ' && buff[12] != '*' && buff[12] != '-')
      rval = 0;
   for( i = 77; i < 80; i++)
      if( buff[i] <= ' ' || buff[i] > 'z')
         rval = 0;
   return( rval);
}

int add_line_to_observation_details( void *obs_details, const char *iline)
{
   observation_details_t *det = (observation_details_t *)obs_details;
   size_t len = strlen( iline);
   const char *valid_lines = "COD CON OBS MEA TEL NET BND COM NUM ACK AC2 ";
   int i, compare = 1, rval;

   if( !memcmp( iline, "COM Sigmas", 10))    /* skip ADES comment lines */
      return( OBS_DETAILS_IRRELEVANT_LINE);
   if( !memcmp( iline, "COM vel (km/s)", 14))    /* skip spacecraft velocity lines */
      return( OBS_DETAILS_IRRELEVANT_LINE);
   while( len && (iline[len - 1] == 10 || iline[len - 1] == 13))
      len--;               /* drop trailing carriage returns/line feeds */
   if( !memcmp( iline, "COD ", 4))
      {
      char tbuff[4];

      det->n_curr = 0;
      if( memcmp( iline + 4, "Multi ", 6))
         {
         memcpy( tbuff, iline + 4, 3);
         tbuff[3] = '\0';
         reset_mpc_code( det, tbuff);
         }
      else       /* details for multiple codes;  see 'details.txt' */
         for( i = 10; i < (int)len - 2; i += 5)
            {
            memcpy( tbuff, iline + i, 3);
            tbuff[3] = '\0';
            reset_mpc_code( det, tbuff);
            }
      }
   if( !det->n_curr)        /* no COD line seen yet */
      return( OBS_DETAILS_IRRELEVANT_LINE);
   if( len == 80 && probable_mpc_record( iline))
      {
      const int idx = find_code_details( det, iline + 77, &compare);

      if( !compare)
         {
         mpc_code_details_t *mptr = det->code_details + idx;

#ifdef DEBUG_PRINTFS
         if( !mptr->observations_found)
            printf( "Got an observation for '%s'\n", mptr->lines[0]);
#endif
         mptr->observations_found = true;
         }
      return( OBS_DETAILS_MPC_80_COLUMN_LINE);
      }
   if( len < 4 || iline[3] != ' ')
      return( OBS_DETAILS_IRRELEVANT_LINE);
   for( i = 0; compare && valid_lines[i]; i += 4)
      compare = memcmp( iline, valid_lines + i, 4);

   if( !compare)    /* yup,  it's a valid header line */
      for( i = 0; i < det->n_curr; i++)
         {
         mpc_code_details_t *mptr = det->code_details + det->curr_idx[i];
         const bool reallocation_needed = (mptr->observations_found
                || (is_power_of_two( mptr->n_lines)
                    && mptr->n_lines >= min_code_lines));

         if( reallocation_needed)
            {
            size_t new_size = min_code_lines;
            char **tptr;

            while( new_size < mptr->n_lines)
               new_size <<= 1;
            tptr = (char **)stack_calloc( det->stack,
                              (2 * new_size + 1) * sizeof( char *));
            memcpy( tptr, mptr->lines, mptr->n_lines * sizeof( char *));
            mptr->lines = tptr;
            }
         mptr->lines[mptr->n_lines] = (char *)stack_calloc( det->stack, len + 1);
         memcpy( mptr->lines[mptr->n_lines], iline, len);
         mptr->n_lines++;
         rval = OBS_DETAILS_HEADER_LINE;
         }
   else
      rval = OBS_DETAILS_IRRELEVANT_LINE;
   return( rval);
}


void free_observation_details( void *obs_details)
{
   observation_details_t *det = (observation_details_t *)obs_details;

   destroy_stack( det->stack);
}

/* init_observation_details allocates the stack (see 'stackall.h/cpp'),
allocates its own return value from that stack,  sets up an initial
'code_details' array,  and nothing else.

Each line read from the file is then fed through get_code_details().  If the
line starts with COD,  we look for it in the existing code_details() array
(using get_code_details()) and either find it there or insert a new entry in
the array.  (For speed, 'code_details' is sorted so we can binary-search.)

   When we see COD,  we make sure that 'lines' is allocated to (say) eight
lines,  n_lines is zeroed,  and 'observations_found' is set to false.

   If it looks like an 80-column MPC observation that matches the current
MPC code,  we set 'observations_found' for that code to be true.

   If it's a CON, OBS, MEA, TEL,  etc. line,  and 'observations_found' is
false,  we can just add the new line to the mpc_code_details_t structure
for that MPC code.  (CON/OBS/MEA lines accumulate;  the others replace as
we see them.)

   If 'observations_found' is true,  we have to allocate a new
mpc_code_details_t structure and copy over the old one and change the
'code_details' pointer accordingly,  and _then_ add the new line. (After
which,  'observations_found' is reset to false.)

   The reason for this is that if we have an observation file such as...

COD XYZ
OBS B. Frank
NET GSC-1.1
(six observations from XYZ)
NET Gaia-DR2
COM Later observations
(three more observations from XYZ)

   ...we'll set up an mpc_code_details_t structure for XYZ when we see
the COD line,  and add the OBS and NET lines.  The first six observations
will be associated with that structure.

   When the second NET line is read,  though,  observations_found == true
tells us we need to allocate a new mpc_code_details_t,  copy over the
existing text,  replace the NET line,  and set (in the new struct)
observations_found = false.  Then we add in the COM line,  and the next
three observations are associated with the new structure.

   The result is that (in the above case) we only allocate two 'details'
structures,  and most of the lines they point to overlap.

   get_code_details() just does a binary search within the code_details
array.

   free_observation_details() should be a single line,  calling stack_free
for the 'stack' variable.  (Everything allocated here ought to be on the
stack-based allocator.)       */
