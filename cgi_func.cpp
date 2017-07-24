#include <stdio.h>
#include <string.h>

/* Functions used for server-based versions of Find_Orb,  sat_id,  and
astcheck (see 'fo_serve.cpp',  'sat_id2.cpp',  cgicheck.cpp' respectively),
and which will probably be used for future server-based code.  For all
of these,  I wanted an avoid_runaway_process() function to ensure that if
something went wrong and the code got hung up,  it would abort after a
decent length of time and give an error message.  (Which really should be
different for each program;  at present,  you get the Find_Orb message.)

   get_multipart_form_data() reads and parses CGI form data,  including
getting uploaded file data.  See the aforementioned three files for
examples of usage.
*/

void avoid_runaway_process( const int max_time_to_run);   /* cgi_func.c */
int get_urlencoded_form_data( const char **idata,       /* cgi_func.c */
                              char *field, const size_t max_field,
                              char *buff, const size_t max_buff);
int get_multipart_form_data( const char *boundary, char *field,
                char *buff, char *filename, const size_t max_len);

#if defined( __linux) || defined( __unix__) || defined( __APPLE__)
#include <sys/time.h>         /* these allow resource limiting */
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>

static const char *email_address_mangled_for_no_spam =
   "p&#x202e;&ocirc;&#xe7;.&ouml;tulp&#x165;c&eacute;j&ocirc;&#x159;p&#x40;ot&uacute;l&#x202c;m";
         /* pluto at projectpluto address with diacritical marks and bidi  */
         /* text... should baffle all but the most devoted spammers.  Put */
         /* your mangled address here.   */

static void sighandler( int signum, siginfo_t *info, void *unused_ucontext)
                                 __attribute__ ((noreturn));

static void sighandler( int signum, siginfo_t *info, void *unused_ucontext)
{
   if( signum == SIGXCPU)
      {
      printf( "\n\n<h1> Ran out of time </h1>\n");
      printf( "<p> Find_Orb was unable to determine an orbit for these\n");
      printf( "observations.  That may just be because the problem exceeded\n");
      printf( "the capabilities of this server;  to avoid runaway processes,\n");
      printf( "there's an intentional fifteen-second cap on getting a solution.\n");
      printf( "<p> You should probably try again with a shorter\n");
      printf( "arc,  or else use the 'native' Windows or Linux or OS/X\n");
      printf( "flavor of Find_Orb.\n");
      }
   else
      {
      printf( "\n\n<h1> Unknown signal </h1>\n");
      printf( " <p> Got signal %d from process %ld.\n", signum,
                      (unsigned long)info->si_pid);
      printf( "I'm still learning a bit about these signals.  Please\n");
      printf( "copy/paste this message and e-mail it to\n");
      printf( "%s\n", email_address_mangled_for_no_spam);
      }
   exit( 0);
}

void avoid_runaway_process( const int max_time_to_run)
{
   struct rlimit r;
   struct sigaction act;

   r.rlim_cur = max_time_to_run;
   r.rlim_max = max_time_to_run + 5;
   setrlimit( RLIMIT_CPU, &r);

   memset(&act, 0, sizeof(act));

   act.sa_sigaction = sighandler;
   act.sa_flags = SA_SIGINFO;
   sigaction(SIGXCPU, &act, NULL);
}
#else
   /* In Win32,  you can't do this sort of process limiting     */
   /* (I think),  so we just have a dummy function.             */

void avoid_runaway_process( const int max_time_to_run)
{
}
#endif         /* _WIN32 */

static int get_urlencoded_piece( const char **idata,
              char *buff, size_t max_buff, const char end_char)
{
   int c;
   const char *tptr = *idata;

   max_buff--;
   while( (c = *tptr++) > 13 && max_buff-- && c != end_char)
      {
      if( c == '+')
         c = ' ';
      else if( c == '%')
         {
         int i, c1;

         c = 0;
         for( i = 0; i < 2; i++)
            {
            c *= 16;
            c1 = *tptr++;
            if( c1 >= '0' && c1 <= '9')
               c += c1 - '0';
            else if( c1 >= 'A' && c1 <= 'F')
               c += c1 - 'A' + 10;
            else                    /* wasn't an hex digit: */
               return( -1);         /* not supposed to happen */
            }
         }
      *buff++ = (char)c;
      }
   *buff =  '\0';
   *idata = tptr;
   return( c);
}

int get_urlencoded_form_data( const char **idata,
                              char *field, const size_t max_field,
                              char *buff, const size_t max_buff)
{
   int c;

   if( get_urlencoded_piece( idata, field, max_field, '=') != '=')
      return( -1);
   c = get_urlencoded_piece( idata, buff, max_buff, '&');
   return( c != '&' && c > 13);
}

int get_multipart_form_data( const char *boundary, char *field,
                char *buff, char *filename, const size_t max_len)
{
   char *tptr, *endptr;
   size_t bytes_read = 0, blen = 0;

   if( filename)
      *filename = '\0';
   while( boundary[blen] >= ' ')
      blen++;
   if( fgets( buff, (int)max_len, stdin)
                  && (tptr = strstr( buff, "name=\"")) != NULL
                  && (endptr = strchr( tptr + 6, '"')) != NULL)
      {
      char *filename_ptr = strstr( tptr, "filename=\"");

      *endptr = '\0';
      strcpy( field, tptr + 6);
      if( filename && filename_ptr
             && (endptr = strchr( filename_ptr + 10, '"')) != NULL)
         {
         *endptr = '\0';
         strcpy( filename, filename_ptr + 10);
         }
      if( fgets( buff, (int)max_len, stdin))
         {
         while( bytes_read < max_len - 1 &&
                 fgets( buff + bytes_read, (int)( max_len - bytes_read), stdin)
                 && memcmp( buff + bytes_read, boundary, blen))
            bytes_read += strlen( buff + bytes_read);
         }
      }
   else
      return( -1);
   while( bytes_read && (buff[bytes_read - 1] == 10
                             || buff[bytes_read - 1] == 13))
      bytes_read--;
   buff[bytes_read] = '\0';
   return( (int)bytes_read);
}
