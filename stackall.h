void *create_stack( const size_t stack_size);
void *stack_alloc( void *stack, const size_t nbytes);
void *stack_calloc( void *stack, const size_t nbytes);
               /* same as stack_alloc( ),  but zeroes the memory */
void destroy_stack( void *stack);
#ifdef NOT_USED_YET
int stack_free( void *stack, void *ptr);
#endif
