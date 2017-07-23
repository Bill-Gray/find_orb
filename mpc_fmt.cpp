#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "afuncs.h"
#include "date.h"

bool is_valid_mpc_code( const char *mpc_code);        /* mpc_fmt.cpp */
double extract_date_from_mpc_report( const char *buff, unsigned *format);
int get_ra_dec_from_mpc_report( const char *ibuff,    /* mpc_fmt.cpp */
                       int *ra_format, double *ra, double *ra_precision,
                       int *dec_format, double *dec, double *dec_precision);

/* MPC has,  at least thus far,  only assigned MPC codes that are an uppercase
letter followed by two digits.  Look at 'rovers.txt',  and you'll see that
user-added "MPC codes" are not so limited;  e.g.,  "2.2" and "0.6" are codes
used for the 2.2-m and 0.6-m telescopes at Mauna Kea.  */

bool is_valid_mpc_code( const char *mpc_code)
{
   int i;

   for( i = 0; i < 3; i++)
      if( mpc_code[i] <= ' ' || mpc_code[i] > 'z')
         return( false);
   return( true);
}


static inline int get_two_digits( const char *iptr)
{
   return( (int)iptr[0] * 10 + (int)iptr[1] - (int)'0' * 11);
}

double minimum_observation_year = -1e+9;   /* set in console Find_Orb's */
double maximum_observation_year =  1e+9;   /* command line      */

#define J2000 2451545.
#define YEAR_TO_JD( year) (J2000 + (year - 2000.) * 365.25)
#define MINIMUM_OBSERVATION_JD YEAR_TO_JD( minimum_observation_year)
#define MAXIMUM_OBSERVATION_JD YEAR_TO_JD( maximum_observation_year)

/* The date/time in an 80-column MPC report line is stored in columns 16-32.
MPC expects this to be in the form YYYY MM DD.dddddd,  with the date usually
given to five digits = .864 second precision;  the sixth digit should only
be given if you're really sure the time is good to that level.  In some
cases (mostly older,  poorly timed observations),  fewer digits are given.
A few ancient observations are just given to the nearest day.

Find_Orb also supports some NON-STANDARD,  FIND_ORB ONLY (i.e.,  MPC will
reject anything you do in these formats).  These formats support greater
precision and spare you the need to convert dates from HH:MM:SS.sss format
to decimal-day format,  and/or the need to convert from JD or MJD format
to YYYY MM SS format.  To demonstrate,  the following lines give the same
date/time in the various formats,  give or take a little rounding error :

            2013 02 13.141593     (MPC's expected format)
            2456336.641592653     (Julian Day format)
            M056336.141592653     (MJD = Modified Julian Day)
            K130213.141592653     (CYYMMDD.dddddd)
            K130213:032353605     (CYYMMDD HH:MM:SS.sss)

The last two use the MPC "mutant hex" convention of K for 21st century
years and J for twentieth-century years,  followed by a two-digit year.
After the two-digit month and two-digit day of month,  one can have
a decimal point and decimal day _or_ a colon and HHMMSS,  and optionally,
up to millisecond precision.  (Though you can and should supply fewer
digits if -- as will almost always be the case -- your timing is not
really that exact.)  The last one is a little confusing,  given the
placement of the ':':  it actually means "2013 02 13 03:23:53.605".

Note that these non-standard formats allow precision up to 10^-9 day
(86.4 microseconds) or,  for the last format,  one millisecond.  This
precision should be good enough for anybody.  Indeed,  that last digit
is just past the precision limit of a 64-bit float.  (We could get
around this by having dates relative to J2000 instead of using MJD
dates,  and/or by using 80-bit "long doubles".  See my get_time.cpp
code in the 'lunar' library for details.  But I don't see it as a big
issue,  at least not yet.)

For each of these formats,  a very little bit of format checking is
done (make sure digits are in certain key places,  and that the full
line is exactly 80 bytes).  Some malformed records _can_ slip past!

Note also that MPC-formatted times are on the UTC scale,  but don't
provide a way to express observations made on the sheer lunacy that is
a leap second.  The HH:MM:SS formats below _should_ let you do this,
but (as yet) do not.  (The other formats flat out can't express that
extra 86401st second in a day.)


      49    M056336.641592653     (MJD, 10^-9 day)
      48    M056336.64159265      (MJD, 10^-8 day)
      47    M056336.6415926       (MJD, 10^-7 day)
      46    M056336.641592        (MJD, 10^-6 day)
      45    M056336.64159         (MJD, 10^-5 day)
      44    M056336.6415          (MJD, 10^-4 day)
      43    M056336.641           (MJD,  .001 day)
      42    M056336.64            (MJD,   .01 day)
      41    M056336.6             (MJD,    .1 day... not supported)
      40    M056336.              (MJD,     1 day... not supported)

      39    K130213.141592653     (High-prec, 10^-9 day)
      38    K130213.14159265      (High-prec, 10^-8 day)
      37    K130213.1415926       (High-prec, 10^-7 day)
      36    K130213.141592        (High-prec, 10^-6 day)
      35    K130213.14159         (High-prec, 10^-5 day)
      34    K130213.1415          (High-prec, 10^-4 day)
      33    K130213.141           (High-prec,  .001 day)
      32    K130213.14            (High-prec,   .01 day)
      31    K130213.1             (High-prec,    .1 day... not supported)
      30    K130213.              (High-prec,     1 day... not supported)

      23    K130213:032353605     (CYYMMDD HH:MM:SS.sss)
      22    K130213:03235360      (CYYMMDD HH:MM:SS.ss)
      21    K130213:0323536       (CYYMMDD HH:MM:SS.s)
      20    K130213:032353        (CYYMMDD HH:MM:SS)

      19    2456336.641592653     (Julian Day, 10^-9 day)
      18    2456336.64159265      (Julian Day, 10^-8 day)
      17    2456336.6415926       (Julian Day, 10^-7 day)
      16    2456336.641592        (Julian Day, 10^-6 day)
      15    2456336.64159         (Julian Day, 10^-5 day)
      14    2456336.6415          (Julian Day, 10^-4 day)
      13    2456336.641           (Julian Day, 10^-3 day)
      12    2456336.64            (Julian Day, 10^-2 day... not supported)
      11    2456336.6             (Julian Day, 10^-1 day... not supported)

       6    2013 02 13.141593     (MPC's expected format, 10^-6 day)
       5    2013 02 13.14159      (MPC's expected format, 10^-5 day)
       4    2013 02 13.1415       (MPC's expected format, 10^-4 day)
       3    2013 02 13.141        (MPC's expected format, 10^-3 day)
       2    2013 02 13.14         (MPC's expected format, 10^-2 day)
       1    2013 02 13.1          (MPC's expected format, 10^-1 day)
       0    2013 02 13.           (MPC's expected format, 10^-0 day) */

double extract_date_from_mpc_report( const char *buff, unsigned *format)
{
   double rval = 0.;
   int year = 0, month = 0;
   size_t start_of_decimals = 0;
   unsigned format_found = 0;
   char tbuff[18];
   const size_t len = strlen( buff);
   unsigned i, bit, digits_mask = 0;

   if( len != 80)             /* check for correct length */
      return( 0.);
   if( buff[12] != ' ' && buff[12] != '*' && buff[12] != '-')
      return( 0.);
   if( !is_valid_mpc_code( buff + 77))
      return( 0.);
   memcpy( tbuff, buff + 15, 17);
   for( i = 0, bit = 1; i < 17; i++, bit <<= 1)
      if( isdigit( tbuff[i]))
         digits_mask |= bit;
   tbuff[17] = '\0';
   if( tbuff[4] == ' ')
      {                       /* standard YYYY MM DD.dddddd format */
      if( (digits_mask & 0x3ff) == 0x36f    /* i.e.,  'dddd dd dd' */
                            && tbuff[7] == ' ' && tbuff[10] == '.')
         {
         int divisor = 1;

         year = atoi( tbuff);
         month = atoi( tbuff + 5);
//       rval = atof( tbuff + 8);
                     /* atof( ) is a little slow,  so we use a little more */
         for( i = 11; i < 17 && tbuff[i] != ' '; i++)  /* code in exchange */
            divisor *= 10;                             /* for better speed */
         rval = (double)atoi( tbuff + 8) +
                           (double)atoi( tbuff + 11) / (double)divisor;
         format_found = 0;
         start_of_decimals = 11;
         }
      }
   else if( *tbuff >= 'H' && *tbuff <= 'K')  /* 18th through 21st century */
      {                                          /* CYYMMDD format */
      if( (tbuff[7] == '.' || tbuff[7] == ':')
               && (digits_mask & 0x3ff) == 0x37e)  /* i.e, 'Zdddddd.dd' */
         {
         year = (*tbuff - 'J') * 100 + 1900 +
                    get_two_digits( tbuff + 1);
         month = get_two_digits( tbuff + 3);
         rval = atof( tbuff + 5);
         if( tbuff[7] == ':')
            {
            rval += (double)get_two_digits( tbuff + 8) / hours_per_day
               + (double)get_two_digits( tbuff + 10) / minutes_per_day
               + (double)get_two_digits( tbuff + 12) / seconds_per_day;
            tbuff[13] = '.';
            rval += atof( tbuff + 13) / seconds_per_day;
            format_found = 20;      /* formats 20-23;  see above */
            start_of_decimals = 14;
            }
         else     /* decimal formats 32-40 */
            {
            format_found = 30;
            start_of_decimals = 8;
            }
         }
      }
   else if( tbuff[7] == '.')        /* MJD or JD format */
      {
      if( (digits_mask & 0x3fe) == 0x37e)   /* i.e., 'zdddddd.dd' */
         {
         if( *tbuff == 'M')    /* MJD */
            {
            format_found = 40;
            rval = 2400000.5 + atof( tbuff + 1);
            }
         else
            {
            format_found = 10;
            rval = atof( tbuff);      /* plain ol' JD */
            }
         start_of_decimals = 8;
         }
      }
   if( format)
      {
      if( start_of_decimals)
         while( isdigit( tbuff[start_of_decimals++]))
            format_found++;
      *format = format_found;
      }

   if( month >= 1 && month <= 12 && rval > 0. && rval < 99.)
      rval += (double)dmy_to_day( 0, month, year,
                                    CALENDAR_JULIAN_GREGORIAN) - .5;

   if( rval < MINIMUM_OBSERVATION_JD || rval > MAXIMUM_OBSERVATION_JD)
      rval = 0.;
             /* Radar obs are always given to the nearest UTC second. So  */
             /* some rounding is usually required with MPC microday data. */
   if( rval && (buff[14] == 'R' || buff[14] == 'r'))
      {
      const double time_of_day = rval - floor( rval);
      const double resolution = 1. / seconds_per_day;
      const double half = .5 / seconds_per_day;

      rval += half - fmod( time_of_day + half, resolution);
      }
   return( rval);
}

/* get_ra_dec() looks at an RA or dec from an MPC report and returns
its precision.  It interprets the formats used by MPC,  plus a lot of
"extended" formats that can be useful if your input data is in other
formats and/or has extra digits.  The return value has the following
meanings ('z' = 'hours or degrees',  i.e.,  this format can apply to
both RAs and decs.)
   hh mm ss.sss       3    (MPC uses this for 'precise' RAs)
   zz mm ss.ss        2    (MPC uses for most RAs & 'precise' decs)
   zz mm ss.s         1    (MPC uses for most decs & low-precision RA)
   zz mm ss           0    (Used by MPC, only rarely)
   zz mm             -1    (Used _very_ rarely by MPC)
   zz mm.m           -2    (Maybe used by MPC once or twice)
   zz mm.mm          -3       The following are for Find_Orb only
   zz mm.mmm         -4
   zz mm.mmmm        -5
   zz mm.mmmmm       -6
   zz mm.mmmmmm      -7
   zz.               100
   zz.z              101
   zz.zz             102
   zz.zzz            103
   zz.zzzz           104... can go up to nine places = 109 in RA,
                           or to 108 = eight places in dec
   ddd.              200  (used for RA only)
   ddd.d             201
   ddd.dd            202... can go up to eight places = 208
   HHMMSSs           307 (RA to .1 second)
   ddmmSSs           307 (dec to .1 arcsecond)
   HHMMSSss          308 (RA to .01 second)
   ddmmSSss          308 (dec to .01 arcsecond)
           ... and so forth until...
   ddmmSSsssss       311 (dec to 10 microarcseconds)
   HHMMSSssssss      312 (RA to one microsecond)
   Undetermined      -99

   Please note that (at least thus far) I've only seen the first six
cases used by MPC,  and they will probably balk at anything sent in
anything but the first four formats.

   The remaining formats have been quite useful in situations where I
find data in a non-MPC format;  I can leave it in decimal degrees or
whatnot,  and can accommodate extra digits for super-high-accuracy
data.  That's why formats 307-312 were added;  they accommodate some
highly precise VLBA astrometry that really _is_ good to the tens of
microarcseconds level.  See

http://iau-comm4.jpl.nasa.gov/plan-eph-data/vlbaobs.html

   And Gaia,  at least in theory,  will be of a similar level of
accuracy;  we need to Be Prepared for that.
   The precision is stored and used to "recreate" the RA/dec in the
original form.  (It could,  and probably should,  also be used to
weight the observations.)                                             */

#define BAD_RA_DEC_FMT           -99

static double get_ra_dec( const char *ibuff, int *format, double *precision)
{
   char buff[13];
   double rval;
   unsigned i = 0, n_digits = 0;
   const bool is_dec = (*ibuff == '-' || *ibuff == '+');
   const bool is_negative = (*ibuff == '-');

   *precision = 1.;   /* in arcseconds */
   if( is_dec)
      ibuff++;
   memcpy( buff, ibuff, 12);
   buff[12] = '\0';
   rval = atof( buff);
   while( isdigit( buff[i]))
      i++;
   if( i > 7)        /* "packed" highly precise RA/decs described above */
      {
      unsigned tval;
      double factor = 1. / 3600.;

      *format = 300 + i;
      n_digits = i - 6;
      buff[6] = '\0';
      tval = atoi( buff);
      rval = (double)( tval / 10000)
           + (double)( (tval / 100) % 100) / 60.
           + (double)( tval % 100) / 3600.;
      for( i = 6; i < 12 && isdigit( ibuff[i]); i++)
         {
         factor *= .1;
         rval += (double)( ibuff[i] - '0') * factor;
         }
      buff[6] = ibuff[6];
//    debug_printf( "Extended: '%s' = %.8f\n", ibuff, rval);
      }
   else if( buff[2] == '.')        /* decimal degrees or hours */
      {
      *precision = 3600.;
      for( i = 3; isdigit( buff[i]) && i < 12; i++)
         n_digits++;
      *format = 100 + n_digits;
      }
   else if( buff[3] == '.')        /* decimal degrees for RA,  ddd.ddd.... */
      {
      *precision = 3600. / 15.;
      rval /= 15.;
      for( i = 4; isdigit( buff[i]) && i < 12; i++)
         n_digits++;
      *format = 200 + n_digits;
      }
   else if( buff[2] == ' ')            /* zz mm(.mmm...) or zz mm ss(.sss...)  */
      {
      rval += atof( buff + 3) / 60.;
      if( buff[5] == ' ' && isdigit( buff[7]))  /* i.e., seconds are given */
         {
         rval += atof( buff + 6) / 3600.;
         for( i = 9; i < 12 && isdigit( buff[i]); i++)
            n_digits++;
         *format = n_digits;
         }
      else           /* minutes: */
         {
         for( i = 6; i < 12 && isdigit( buff[i]); i++)
            n_digits++;
         *format = -((int)n_digits + 1);
         *precision = 60.;
         }
      }
   else
      *format = BAD_RA_DEC_FMT;
   while( n_digits--)
      *precision *= .1;
   if( is_negative)
      rval = -rval;
   return( rval);
}

/* Extracts both the RA and the dec from an MPC-formatted report,  subject
to all the weirdnesses described above.  Return value is zero if both
RA and dec are successfully found;  -1 if the RA is bad;  -2 if the dec
is bad;  -3 if neither was read. */

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

int get_ra_dec_from_mpc_report( const char *ibuff,
                       int *ra_format, double *ra, double *ra_precision,
                       int *dec_format, double *dec, double *dec_precision)
{
   int rval = 0;

   *ra  = get_ra_dec( ibuff + 32, ra_format, ra_precision) * (PI / 12.);
   *ra_precision *= 15.;     /* cvt minutes/seconds to arcmin/arcsec */
   if( *ra_format == BAD_RA_DEC_FMT)
      rval = -1;
   *dec =  get_ra_dec( ibuff + 44, dec_format, dec_precision) * (PI / 180.);
   if( *dec_format == BAD_RA_DEC_FMT)
      rval -= 2;
   return( rval);
}
