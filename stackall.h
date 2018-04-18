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

void *create_stack( const size_t stack_size);
void *stack_alloc( void *stack, const size_t nbytes);
void *stack_calloc( void *stack, const size_t nbytes);
               /* same as stack_alloc( ),  but zeroes the memory */
void destroy_stack( void *stack);
#ifdef NOT_USED_YET
int stack_free( void *stack, void *ptr);
#endif
