#include <string.h>

#if defined( _WIN32) || defined( __WATCOMC__)
#include <malloc.h>     /* for alloca() prototype */
#else
#include <stdlib.h>
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
to do so because Sort Algorithms Are Fun,  but I've got too much else to do.

   A note about the gap sequence:  it used to be gap = gap * 3 + 1,  for
a sequence of 1, 4, 13, 40, 121,... I accepted this commonly-used sequence
uncritically for some years,  until testing found it to be _much_ slower
than gap = gap * 8 / 3 + 1 (i.e.,  a smaller gap exponent).  I did some
more testing to search for gaps that would be pairwise relatively prime,
and ideally relatively prime for nearby but non-adjacent gaps.  Hence
the value of gap0 = 250104703.  (For smaller arrays,  this will swiftly
be brought down to a "reasonable" size.)  The following explicitly lists
the gap sequence and their factorizations.

31011856950330734 = 2 193 80341598316919
11629446356374025 = 5 5 31 15005737234031
4361042383640259 = 3 1453680794546753
1635390893865097 = 1635390893865097
613271585199411 = 3 3 3 29 1229 15661 40693
229976844449779 = 7 7 4693404988771
86241316668667 = 103 1069 1709 458309
32340493750750 = 2 5 5 5 457 283067779
12127685156531 = 7 109891 15765863
4547881933699 = 13 31 1531 1721 4283
1705455725137 = 23 83 893376493
639545896926 = 2 3 3 17 61 281 121931
239829711347 = 97 3067 806153
89936141755 = 5 19 53 17862193
33726053158 = 2 7 2409003797
12647269934 = 2 13 19 61 419701
4742726225 = 5 5 189709049
1778522334 = 2 3 23 23 560341
666945875 = 5 5 5 5335567
250104703 = 97 103 25033
93789263 = 19 41 120397
35170973 = 35170973
13189114 = 2 6594557
4945917 = 3 223 7393
1854718 = 2 83 11173
695519 = 11 53 1193
260819 = 13 20063
97807 = 47 2081
36677 = 36677
13753 = 17 809
5157 = 3 3 3 191
1933 = 1933
724 = 2 2 181
271 = 271
101 = 101
37 = 37
13 = 13
4 = 2 2
1 = 1                    */

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
   size_t gap = 250104703;
   char *data = (char *)base;
   char *pivot = (char *)alloca( elem_size);

   while( gap < n_elements)
      gap = gap * 8 / 3 + 1;
   while( (gap = gap * 3 / 8) != 0)
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
