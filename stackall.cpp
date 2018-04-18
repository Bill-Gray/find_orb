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

#include <stdlib.h>
#include <string.h>
#include "stackall.h"

/* 'stack allocation' allows for simple,  fast memory allocation
in situations where you want to allocate many small bits of memory,
then free them all at once.  (In theory,  it could be a truly
'stack' sort of allocation,  in which you could deallocate memory
in the reverse order in which it was allocated.  There is actually
a function to do this,  not yet in use and untested.)

The idea is that you call create_stack( ) with a default size that
is reasonably larger than the average chunk of memory you expect to
allocate from the stack.  You then call stack_alloc( ) to allocate
all of those chunks.  Say you set the stack size to 2048 bytes,  and
you were allocating space for lines from a text file that averaged
about sixty bytes.  Every thirty or forty allocations,  a new block
of 2048 bytes would be malloc()ed and the next batch of stack allocations
would be drawn from that.  The blocks in question are maintained in a
linked list.

Once done with the allocations in question,  one would call stack_destroy(),
which would go through the linked list and free all the blocks.
*/

#define STACK struct stack

STACK
   {
   STACK *next;
   size_t size, used;
   };

void *create_stack( const size_t stack_size)
{
   STACK *rval = (STACK *)malloc( stack_size + sizeof( STACK));

   rval->size = stack_size;
   rval->used = 0;
   rval->next = NULL;
   return( rval);
}

void *stack_alloc( void *stack, const size_t nbytes)
{
   STACK *sptr = (STACK *)stack;
   const size_t size0 = sptr->size;
   char *rval;

   if( sptr->used + nbytes > size0)
      sptr = sptr->next;
   if( !sptr || sptr->used + nbytes > size0)
      {
      STACK *tptr = (STACK *)create_stack( (size0 > nbytes) ? size0 : nbytes);

      sptr = (STACK *)stack;
      tptr->next = sptr->next;
      sptr->next = tptr;
      sptr = tptr;
      }
   rval = (char *)( sptr + 1) + sptr->used;
   sptr->used += nbytes;
   return( rval);
}

void *stack_calloc( void *stack, const size_t nbytes)
{
   void *rval = stack_alloc( stack, nbytes);

   if( rval)
      memset( rval, 0, nbytes);
   return( rval);
}

void destroy_stack( void *stack)
{
   while( stack)
      {
      STACK *nextptr = ((STACK *)stack)->next;

      free( stack);
      stack = nextptr;
      }
}

#ifdef NOT_USED_YET
int stack_free( void *stack, void *ptr)
{
   STACK *sptr = (STACK *)stack;
   int rval = -1;

   while( rval && sptr)
      {
      size_t diff = (char *)ptr - (char *)( sptr + 1);

      if( diff < sptr->used)
         if( (char *)( sptr + 1) + diff == ptr)
            {
            STACK *tptr = sptr->next;

            sptr->used = diff;
            sptr->next = NULL;
            while( tptr)
               {
               STACK *next_ptr = tptr->next;

               free( tptr);
               tptr = next_ptr;
               }
            rval = 0;         /* got it */
            }
      sptr = sptr->next;
      }
   return( rval);
}
#endif
