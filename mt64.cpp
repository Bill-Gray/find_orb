/*
   A C-program for MT19937-64
   Coded by Takuji Nishimura and Makoto Matsumoto.  (With a lot of
   changes by Bill J. Gray;  see below.)

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_mt64(seed, state_array)
   or init_mt64_by_array( init_key, key_length, state_array).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
   (BEFORE PROVIDING FEEDBACK,  BE AWARE THAT THE FOLLOWING IS NOT THE
   ORIGINAL CODE!  It's been very heavily revised,  to the point where
   it may not be particularly familiar to the original authors.  You may
   want to send feedback to me,  Bill Gray,  pluto(at)projectpluto(dot)com.)

   2011 Apr 2:  BJG:  removed static variables to make code
   re-entrant (and therefore thread-safe,  I hope).  Changed
   'unsigned long long's to 'uint64_t'.  Used a UCONST( ) macro
   to allow the code to compile in Microsoft VC,  which has a somewhat
   different way (of course) for defining 64-bit constants.

      There is now a function to get N bits at a time,  where N <= 64.
   Certain ways of getting double-precision random numbers are given.
   I moved the 'test' routines into a separate mt64test.c module,
   added a timing test,  and created a header file.
*/

#include <stdint.h>
#include "mt64.h"

#ifndef UINT64_C
   #ifdef _MSC_VER
      #define UINT64_C( a) (a##ui64)
   #else
      #define UINT64_C( a) (a##ULL)
   #endif
#endif

#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
   #define DIAGNOSTIC_PRAGMAS_IN_FUNCTIONS_OK
#endif

#define MM 156
#define MATRIX_A UINT64_C( 0xB5026F5AA96619E9)
#define UM       UINT64_C( 0xFFFFFFFF80000000)   /* Most significant 33 bits */
#define LM       UINT64_C( 0x7FFFFFFF)           /* Least significant 31 bits */

#define MT_SIZE 312

#define LOC             mt[MT_SIZE]
#define N_BITS_LEFT     mt[MT_SIZE + 1]
#define LEFTOVER_BITS   mt[MT_SIZE + 2]

              /* Initializes mt[MT_SIZE] with a seed,  using a  */
              /* modified linear congruential generator.        */
void init_mt64( const uint64_t seed, uint64_t *mt)
{
   int i;

   mt[0] = seed;
   for( i = 1; i < MT_SIZE; i++)
        mt[i] =  (UINT64_C( 6364136223846793005) * (mt[i-1] ^ (mt[i-1] >> 62)) + i);
   LOC = MT_SIZE;
   N_BITS_LEFT = (uint32_t)0;
   LEFTOVER_BITS = (uint64_t)0;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_mt64_by_array( const uint64_t init_key[],
          const uint64_t key_length, uint64_t *mt)
{
    uint64_t i, j, k;

    init_mt64( UINT64_C( 19650218), mt);
    i=1; j=0;
    k = (MT_SIZE>key_length ? MT_SIZE : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * UINT64_C( 3935559000370003845)))
          + init_key[j] + j; /* non linear */
        i++; j++;
        if (i>=MT_SIZE) { mt[0] = mt[MT_SIZE-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=MT_SIZE-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * UINT64_C( 2862933555777941757)))
          - i; /* non linear */
        i++;
        if (i>=MT_SIZE) { mt[0] = mt[MT_SIZE-1]; i=1; }
    }

    mt[0] = UINT64_C( 1) << 63; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0, 2^64-1]-interval */
uint64_t mt64( uint64_t * RESTRICT mt)
{
    uint64_t x;

    if (LOC >= MT_SIZE) { /* generate MT_SIZE words at one time */
        static const uint64_t mag01[2]={ UINT64_C( 0), MATRIX_A};
        uint64_t * RESTRICT endptr = mt + MT_SIZE - MM;

        while( mt < endptr)
            {
            x = (mt[0] & UM) | (mt[1] & LM);
            *mt = mt[MM] ^ (x >> 1) ^ mag01[(unsigned)x & 1u];
            mt++;
            }
        endptr += MM - 1;
        while( mt < endptr)
            {
            x = (mt[0] & UM) | (mt[1] & LM);
            *mt = mt[MM - MT_SIZE] ^ (x >> 1) ^ mag01[(unsigned)x & 1u];
            mt++;
            }
        mt -= MT_SIZE - 1;
        x = (mt[MT_SIZE-1] & UM) | (mt[0] & LM);
        mt[MT_SIZE-1] = mt[MM-1] ^ (x >> 1) ^ mag01[(unsigned)x & 1u];

        LOC = 0;
    }

    x = mt[LOC++];

    x ^= (x >> 29) & UINT64_C( 0x5555555555555555);
    x ^= (x << 17) & UINT64_C( 0x71D67FFFEDA60000);
    x ^= (x << 37) & UINT64_C( 0xFFF7EEE000000000);
    x ^= (x >> 43);

    return x;
}

/* The following (currently unused) function allows one to
get n_bits of random bits at a time.  Thus,  for example,
mt64_get_bits( mt, 2) will return a random value from 0 to
3.  The benefit of this is that the full mt64() function
would then be called only once every 32 times,  and one
has somewhat faster code.

   Another example:  if you called mt64_get_bits( mt, 25)
three times in a row,  mt64( ) would be called the first time,
generating 64 new random bits.  25 would be returned to the
user,  with 39 saved in LEFTOVER_BITS (and therefore,
N_BITS_LEFT == 39).  On the second call,  25 more bits would
be extracted from LEFTOVER_BITS,  with N_BITS_LEFT reset
to 39-25=14.  On the third call,  we would see that we had
only fourteen random bits;  another call to mt64() would be
made,  with 25 of the resulting bits returned to the user
and the remaining 39 put back into LEFTOVER_BITS.  That
would reset N_BITS_LEFT to 14+39 = 53 random bits.

   The potential pitfall is that it assumes that mt64()
really returns random bits,  and that no unexpected
correlations will appear when examining them in,  say,
seven-bit increments. With a truly "random" source,  the
order in which bits are sent is obviously unimportant. But I
don't _know_ for a certainty that MT-64 is sufficiently
"random",  and most tests run on it have presumably looked
for correlations over 64-bit intervals rather than,  say,
seven-bit ones. As was said of the infamous RANDU generator
(in _Numerical Recipes in C_): "We guarantee that each number
is random individually,  but we don't guarantee that more
than one of them is random."  MT64 is immensely more "random"
than that generator,  and I believe it's safe to get bits in
this manner (and if it isn't,  the solution should be to use
a still better PRNG, rather than ignoring the problem;  the
"randomness" of a PRNG should not depend on the order in
which bits are examined.) But,  "caveat user". */

#ifdef CURRENTLY_UNUSED
uint64_t mt64_get_bits( uint64_t *mt, const int n_bits)
{
   uint64_t rval;

   if( (unsigned)n_bits <= N_BITS_LEFT)
      {                             /* yes, can just use bits     */
      rval = LEFTOVER_BITS;         /* computed on previous calls */
      LEFTOVER_BITS >>= n_bits;
      N_BITS_LEFT -= n_bits;
      }
   else                          /* no,  gotta generate new bits */
      {
      rval = mt64( mt);
      LEFTOVER_BITS ^= rval >> (n_bits - N_BITS_LEFT);
      N_BITS_LEFT += 64 - n_bits;
      }
            /* Mask off and return only the lowest 'n_bits' bits: */
   return( rval & (((uint64_t)1 << n_bits) - (uint64_t)1));
}
#endif

/* The following generates a double on the [0,1) interval.  I
tried three methods:  (1) the original one in the MT64 code,
which takes a pseudo-random number,  shifts it down,  casts it
to double, and multiplies it by 2^-53;  (2) a modification of a
method described in _Numerical Recipes_ that assumes 64-bit
doubles in IEEE format,  but which gets a type-punning warning
in gcc; and (3) almost the same method,  except by doing a
memcpy(), the type-punning warning is avoided.  (Any of the
three methods can be selected with suitable #defines.)

   All three methods work Just Fine.  The last two could have
trouble if doubles aren't 64 bits or if they aren't in IEEE
format.  And it turns out (after running timing tests) that
the performance of all three is essentially identical.  About
80% of the time is spent in mt64() anyway;  the efforts to
convert the result to a double take up the remaining 20%,
and that did not vary much between the three methods.

   Incidentally,  I also tried the commented-out version:

// const uint64_t rval = jflone | mt64_get_bits( mt, 52);

   My thought was that this would cut down on calls to mt64()
by a factor of 52/64 = .8125,  since we'd be only using 52
random bits at a time.  However,  mt64_get_bits() has some
overhead;  it actually was a little slower this way.  It
_could_ be a useful trick if one's PRNG was slower,  though.

   Another completely commented out version attempts to be
faster and provide somewhat better results.  It's actually
a little slower.  The "better results" refers to the fact that
the "usual" mt64_double has a granularity of 2^-53,  even
though doubles can represent numbers from .25 to .5 down to
2^-54,  those from .125 to .25 down to 2^-55,  and so on.
It's hard to persuade me that this is really a problem.
But if it were,  something resembling the commented-out
version would address it. */

#ifdef AVOID_TYPE_PUNNING
#include <string.h>     /* for memcpy( ) prototype */
#endif

double mt64_double( uint64_t * RESTRICT mt)
{
#ifdef NON_IEEE_WAY
    const double two_to_the_53rd_power = 9007199254740992.0;

    return (int64_t)(mt64( mt) >> 11) * (1.0 / two_to_the_53rd_power);
#else
   const uint64_t jflone = UINT64_C( 0x3ff0000000000000);
   const uint64_t jflmsk = UINT64_C( 0x000fffffffffffff);
   const uint64_t rval = jflone | (jflmsk & mt64( mt));
// const uint64_t rval = jflone | mt64_get_bits( mt, 52);

#ifndef AVOID_TYPE_PUNNING
   #ifdef DIAGNOSTIC_PRAGMAS_IN_FUNCTIONS_OK
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wstrict-aliasing"
   #endif
   return( *((const double * RESTRICT)&rval) - 1.);
   #ifdef DIAGNOSTIC_PRAGMAS_IN_FUNCTIONS_OK
      #pragma GCC diagnostic pop
   #endif
#else
   double d_rval;

   memcpy( &d_rval, &rval, sizeof( double));
   d_rval -= 1.;
   return( d_rval);
#endif
#endif
}

// double mt64_double( uint64_t * RESTRICT mt)
// {
//    const uint64_t jflmsk = UINT64_C( 0x000fffffffffffff);
//    uint64_t rval = mt64( mt);
//
//    if( rval & (UINT64_C( 1) << 63))    /* .5 <= returned value < 1 */
//       rval = UINT64_C( 0x3fe0000000000000) | (jflmsk & rval);
//    else if( rval & (UINT64_C( 1) << 62))    /* .25 <= returned value < .5 */
//       rval = UINT64_C( 0x3fd0000000000000) | (jflmsk & rval);
//    else if( rval & (UINT64_C( 1) << 61))    /* .125 <= returned value < .25 */
//       rval = UINT64_C( 0x3fc0000000000000) | (jflmsk & rval);
//    else
//       {
//       const double two_to_the_64th_power = 65536. * 65536. * 65536. * 65536.;
//       return( (double)( (int64_t)rval) * (1. / two_to_the_64th_power));
//       }
// #ifndef AVOID_TYPE_PUNNING
//    return( *((const double * RESTRICT)&rval));
// #else
//    double d_rval;
//
//    memcpy( &d_rval, &rval, sizeof( double));
//    return( d_rval);
// #endif
// }
