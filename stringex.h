/* Some extensions to the usual string/snprintf functions.  First,
the BSD strlcxx functions.  See comments in 'miscell.cpp'. */

size_t strlcpy( char *dst, const char *src, size_t dsize);   /* miscell.c */
size_t strlcat( char *dst, const char *src, size_t dsize);   /* miscell.c */

/* The above truncate silently.  On occasion,  one wants that.  On other
occasions,  you want the program to abort to let you know that a buffer
would have overflowed absent the checking;  in such cases,  use these : */

size_t strlcpy_err( char *dst, const char *src, size_t dsize); /* miscell.c */
size_t strlcat_err( char *dst, const char *src, size_t dsize); /* miscell.c */

/* Frequently,  'dsize' is simply sizeof( dst).  The following macros
can result in slightly easier to read code in such cases.      */

#define strlcat_error( a, b) strlcat_err( a, b, sizeof( a))
#define strlcpy_error( a, b) strlcpy_err( a, b, sizeof( a))

   /* Older MSVC/C++ lacks snprintf.  See 'ephem0.cpp' for details. */
#if defined(_MSC_VER) && _MSC_VER < 1900
int snprintf( char *string, const size_t max_len, const char *format, ...);
#endif

int snprintf_append( char *string, const size_t max_len,      /* ephem0.cpp */
                                   const char *format, ...)
#ifdef __GNUC__
         __attribute__ (( format( printf, 3, 4)))
#endif
;

int snprintf_err( char *string, const size_t max_len,      /* miscell.cpp */
                                   const char *format, ...)
#ifdef __GNUC__
         __attribute__ (( format( printf, 3, 4)))
#endif
;
