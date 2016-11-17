#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define MAX_POLY 10

int find_real_polynomial_roots( const double *poly, int poly_degree,
                                double *real_roots);        /* roots.cpp */

/* Compile with

g++ -Wall -O3 -Wextra -pedantic -o roottest roottest.cpp roots.cpp

This tests the code in 'roots.cpp' for finding all real roots of a
polynomial with all real coefficients.  Run it with without command
line arguments to get a description.   */

static int extract_polynomial( const char *buff, double *poly)
{
   unsigned i;

   for( i = 0; i < MAX_POLY; i++)
      poly[i] = 0.;
   while( *buff)
      {
      double coeff = atof( buff);

      if( coeff == 0.)
         coeff = (*buff == '-' ? -1. : 1.);
      while( *buff && *buff != 'x')
         buff++;
      if( !*buff)    /* constant term */
         i = 0;
      else if( buff[1] == '^')
         {
         buff += 2;
         i = (unsigned)atoi( buff);
         if( i <= 1)
            return( -1);
         while( *buff >= '0' && *buff <= '9')
            buff++;
         }
      else                    /* linear term */
         {
         buff++;
         i = 1;
         }
      poly[i] = coeff;
      }
   return( 0);
}

static int extract_full_polynomial( const char *buff, double *poly)
{
   unsigned i, j;
   bool first_time = true;

   if( *buff != '(')       /* no factors;  plain old polynomial */
      return( extract_polynomial( buff, poly));
   while( *buff)
      {
      char tbuff[100];

      while( *buff && *buff != '(')
         buff++;
      buff++;
      for( i = 0; buff[i] && buff[i] != ')'; i++)
         tbuff[i] = buff[i];
         ;
      if( !buff[i])     /* didn't find closing paren */
         return( -3);
      tbuff[i] = '\0';
      buff += i + 1;
      if( first_time)
         {
         first_time = false;
         extract_polynomial( tbuff, poly);
         }
      else        /* multiplying by a second,  or third,  or... term */
         {
         double tpoly[MAX_POLY], product[MAX_POLY];

         extract_polynomial( tbuff, tpoly);
         for( i = 0; i < MAX_POLY; i++)
            product[i] = 0.;
         for( i = 0; i < MAX_POLY; i++)
            for( j = 0; j + i < MAX_POLY; j++)
               product[i + j] += poly[i] * tpoly[j];
         for( i = 0; i < MAX_POLY; i++)
            poly[i] = product[i];
         }
      }
   return( 0);
}

static void remove_spaces( char *buff)
{
   char *optr = buff;

   while( *buff)
      {
      if( *buff != ' ')
         *optr++ = *buff;
      buff++;
      }
   *optr = '\0';
}

int main( const int argc, const char **argv)
{
   char poly_text[1000];
   double poly[MAX_POLY], roots[MAX_POLY];
   int err_code = -99, i, order = 0, n_roots;

   *poly_text = '\0';
   for( i = 1; i < argc; i++)
      strcat( poly_text, argv[i]);
   remove_spaces( poly_text);
   if( *poly_text)
      err_code = extract_full_polynomial( poly_text, poly);
   if( err_code)
      {
      printf( "Error in parsing polynomial: %d\n", err_code);
      printf( "'roottest' will find all _real_ roots of a given polynomial with\n");
      printf( "real coefficients,  guaranteed,  but no complex ones. For example :\n\n");
      printf( "./roottest \"(12 x^2 + x - 1)*(x-2)\"\n");
      printf( "./roottest +12 * x ^ 3 -23 * x ^ 2 -3 * x+2\n\n");
      printf( "would both result in finding roots at x=-1/3, x=1/4,  and x=2.\n");
      printf( "Method is described in 'roots.cpp'.\n");
      return( -1);
      }
   for( i = MAX_POLY - 1; i >= 1; i--)
      if( poly[i] != 0.)
         {
         if( poly[i] == 1.)
            printf( " + x");
         else if( poly[i] == -1.)
            printf( " - x");
         else
            printf( " %+g * x", poly[i]);
         if( i > 1)
            printf( " ^ %d", i);
         }
   if( poly[0] != 0.)
      printf( "%+g", poly[0]);
   printf( "\n");

   for( i = 0; i < MAX_POLY; i++)
      if( poly[i] != 0.)
         order = i;
   n_roots = find_real_polynomial_roots( poly, order, roots);
   for( i = 0; i < n_roots; i++)
      printf( "%f is a root\n", roots[i]);
   return( 0);
}

