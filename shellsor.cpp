#include <string.h>

#ifdef _WIN32
#include <malloc.h>     /* for alloca() prototype */
#endif

/* The re-entrant version of qsort(),  qsort_r(),  is implemented in
subtly but infuriatingly different ways on different platforms (including
as 'qsort_s' on MS Windows).  You can "sort of" work around these,  as

https://github.com/noporpoise/sort_r

   shows.  But I found it easier to write the following.

https://sourceware.org/ml/libc-alpha/2008-12/msg00007.html

   makes an interesting argument that the Linux qsort_r,  which has the
"context" pointer in the comparison function as the third argument,  is
the "right" way of doing things.  (*BSD and Windows have the context
pointer as the first argument.)  The reason is that it means your
non-recursive sort can basically just call the recursive sort,  with
an ignored NULL context pointer.  (Which,  you'll see,  I've done below.)

   Note that Shellsort is decently fast,  but not nearly as fast as Quicksort.
If sorting speed matters in your application,  either ignore this or modify
the Shellsort to actually be qsort or something similarly fast... I'm tempted
to do so because Sort Algorithms Are Fun,  but I've got too much else to do. */

#if defined _GNU_SOURCE
   #define HAVE_REENTRANT_QSORT
   #include <stdlib.h>
#endif

void shellsort_r( void *base, const size_t n_elements, const size_t elem_size,
         int (*compare)(const void *, const void *, void *), void *context);
void shellsort( void *base, const size_t n_elements, const size_t elem_size,
         int (*compare)(const void *, const void *));    /* shellsor.cpp */

void shellsort_r( void *base, const size_t n_elements, const size_t elem_size,
         int (*compare)(const void *, const void *, void *), void *context)
{
/* printf( "Sort begins: %f, %u elements, size %u\n",
                    (double)clock( ) / (double)CLOCKS_PER_SEC,
                    (unsigned)n_elements, (unsigned)elem_size);   */
#ifdef HAVE_REENTRANT_QSORT
   qsort_r( base, n_elements, elem_size, compare, context);
#else
   size_t gap = 1;
   char *data = (char *)base;
   char *pivot = (char *)alloca( elem_size);

   while( gap < n_elements)
      gap = gap * 3 + 1;
   while( gap /= 3)
      {
      size_t j;
      const size_t spacing = elem_size * gap;

      for( j = gap; j < n_elements; j++)
         {
         char *tptr = data + j * elem_size;
         char *tptr2 = tptr - spacing;

         if( (compare)( tptr2, tptr, context) > 0)
            {
            memcpy( pivot, tptr, elem_size);
            memcpy( tptr, tptr2, elem_size);
            tptr = tptr2;
            tptr2 -= spacing;
            while( tptr2 >= base && (compare)( tptr2, pivot, context) > 0)
               {
               memcpy( tptr, tptr2, elem_size);
               tptr = tptr2;
               tptr2 -= spacing;
               }
            memcpy( tptr, pivot, elem_size);
            }
         }
      }
#endif
/* printf( "Sort ends: %f\n", (double)clock( ) / (double)CLOCKS_PER_SEC); */
}

void shellsort( void *base, const size_t n_elements, const size_t elem_size,
         int (*compare)(const void *, const void *))
{
   shellsort_r( base, n_elements, elem_size,
               (int (*)(const void *, const void *, void *))compare, NULL);
}
