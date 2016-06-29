#ifndef RESTRICT
   #if defined (__cplusplus) && !defined( _MSC_VER)
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
