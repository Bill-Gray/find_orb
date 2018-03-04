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

#ifndef RESTRICT
   #if defined (__cplusplus) && !defined( _MSC_VER) && !defined( __WATCOMC__)
      #define RESTRICT __restrict
   #else
      #define RESTRICT
   #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

void init_mt64( const uint64_t seed, uint64_t *mt);
void init_mt64_by_array( const uint64_t init_key[],
          const uint64_t key_length, uint64_t *mt);
uint64_t mt64( uint64_t * RESTRICT mt);
uint64_t mt64_get_bits( uint64_t *mt, const int n_bits);
double mt64_double( uint64_t * RESTRICT mt);

#ifdef __cplusplus
}
#endif  /* #ifdef __cplusplus */

  /* An MT-64 state contains 312 64-bit values for the actual PRNG; */
  /*  plus an index value telling us where we are within the        */
  /* "cycle",  ranging from 1 to 312;  plus a value telling us how  */
  /* many leftover bits are stored for use in the mt64_get_bits()   */
  /* function;  plus the (up to) 64 bits actually stored.  That     */
  /* means we need an array of 312 + 1 + 1 + 1 = 315 int64_ts.      */

#define MT_STATE_SIZE    315
