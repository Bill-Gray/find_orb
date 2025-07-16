/* #including this file,  instead of the 'usual' <assert.h>,  will
cause assert()s to always be evaluated,  even in non-debug mode.
This is useful if

   -- you really do always want to evaluate the assertions;  or
   -- you have variables that are only used in the assert()s,  and
which result in "unused variable" warnings in NDEBUG mode.  (Or,  with
the usual -Werror flag,  errors.)

   Such warnings/errors can be ugly to code around.  If the assert() in
question is inexpensive (not evaluated inside a performance-intensive
loop),  it may be easier to force the assert() to be evaluated,  causing
the warning or error to go away.

   The ideal thing would be to have an assert_even_in_debug() function or
macro.  This would be useful if,  for example,  the assert() is evaluated
frequently enough to be a performance issue.  At present,  I don't have that
problem;  if I eventually do,  I may re-visit this bit of code.      */

#ifdef NDEBUG          /* assert() would not normally be compiled */
   #define TEMP_NDEBUG NDEBUG
   #undef NDEBUG
   #include <assert.h>
   #define NDEBUG TEMP_NDEBUG
#else
   #include <assert.h>
#endif
