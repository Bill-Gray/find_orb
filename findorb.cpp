/* findorb.cpp: main driver for console Find_Orb

Copyright (C) 2010, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */

#define _XOPEN_SOURCE_EXTENDED   1
#define PDC_NCMOUSE
#define PDC_FORCE_UTF8
#define MOUSE_MOVEMENT_EVENTS_ENABLED

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#if defined( _WIN32)
   #ifdef MOUSE_MOVED
      #undef MOUSE_MOVED
   #endif
   #ifdef COLOR_MENU
      #undef COLOR_MENU
   #endif
#endif
   #include "curses.h"

      /* The 'usual' Curses library provided with Linux lacks a few things */
      /* that PDCurses and MyCurses have, such as definitions for ALT_A    */
      /* and such.  'curs_lin.h' fills in these gaps.   */
#ifndef ALT_A
   #include "curs_lin.h"
#endif

#ifndef BUTTON_CTRL
   #define BUTTON_CTRL BUTTON_CONTROL
#endif

   /* On version 1 of ncurses,  the mouse mask is constrained to 32 bits and
   there's no way to express the 'fifth button', which would otherwise give
   you a mouse-wheel-up event.  Instead,  the button is returned as zero.
   Note that that's also what you get with a 'mouse move' event,  and we can't
   tell them apart.  So when 'mouse move' events are enabled -- currently done
   only when there's a popup window on-screen -- we can't detect mouse wheel
   up events.  At least,  not on version 1 of ncurses. */

#ifdef BUTTON5_PRESSED
   #define button5_pressed (button & BUTTON5_PRESSED)
#else
   #define button5_pressed (!button)
   #define BUTTON5_PRESSED 0
#endif

#ifdef __PDCURSES__
   #define BUTTON_MODIFIERS  (BUTTON_MODIFIER_SHIFT | BUTTON_MODIFIER_CONTROL | BUTTON_MODIFIER_ALT)
#else       /* ncurses lacks these */
   #define BUTTON_MODIFIERS 0
#endif

#if !defined( MOUSE_WHEEL_SCROLL)
   #define MOUSE_WHEEL_SCROLL 0
#endif

#define CSI "\x1b["
#define OSC "\x1b]"

#define default_mouse_events (BUTTON1_CLICKED | BUTTON1_DOUBLE_CLICKED \
                            | BUTTON2_CLICKED | BUTTON2_DOUBLE_CLICKED \
                            | BUTTON3_CLICKED | BUTTON3_DOUBLE_CLICKED \
                            | MOUSE_WHEEL_SCROLL | BUTTON_MODIFIERS    \
                            | BUTTON4_PRESSED | BUTTON5_PRESSED)

/* Some terminals don't actually report mouse movements.  'Hint' text
relies on such reports;  we show no hints until at least one movement
report is received. */

static bool _mouse_movements_are_reported = false;

#include <wchar.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>
#include <sys/stat.h>
#include <locale.h>
#include "watdefs.h"
#include "sigma.h"
#include "afuncs.h"
#include "comets.h"
#include "mpc_obs.h"
#include "date.h"
#include "monte0.h"
#include "stringex.h"
#include "constant.h"

int debug_level = 0;

extern unsigned perturbers;

#define AUTO_REPEATING                31002
#define KEY_ADD_MENU_LINE             31004
#define KEY_REMOVE_MENU_LINE          31005
#define KEY_CYCLE_RESID_DISPLAY_UP    31006
#define KEY_CYCLE_RESID_DISPLAY_DN    31007
#define KEY_OBSCODE_CLICKED           31008

/* You can cycle between showing only the station data for the currently
selected observation;  or the "normal" having,  at most,  a third of
the residual area devoted to station data;  or having most of that area
devoted to station data.   */

#define SHOW_MPC_CODES_ONLY_ONE        0
#define SHOW_MPC_CODES_NORMAL          1
#define SHOW_MPC_CODES_MANY            2

#define button1_events (BUTTON1_CLICKED | BUTTON1_DOUBLE_CLICKED \
      | BUTTON1_TRIPLE_CLICKED | BUTTON1_PRESSED | BUTTON1_RELEASED)

#define BUTTON_PRESS_EVENT (BUTTON1_PRESSED | BUTTON2_PRESSED | BUTTON3_PRESSED)

      /* Not all Curses allow the following attributes. */

#ifndef A_ITALIC
   #define A_ITALIC   0
#endif
#ifndef A_LEFTLINE
   #define A_LEFTLINE 0
#endif
#ifndef A_RIGHTLINE
   #define A_RIGHTLINE 0
#endif
#ifndef A_UNDERLINE
   #define A_UNDERLINE 0
#endif
#ifndef A_OVERLINE
   #define A_OVERLINE 0
#endif

#ifndef max
#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)
#endif

#define CTRL(c) ((c) & 0x1f)

void ensure_config_directory_exists(); /* miscell.c */
static int user_select_file( char *filename, const char *title, const int flags);
double get_planet_mass( const int planet_idx);                /* orb_func.c */
int simplex_method( OBSERVE FAR *obs, int n_obs, double *orbit,
               const double r1, const double r2, const char *constraints);
int superplex_method( OBSERVE FAR *obs, int n_obs, double *orbit, const char *constraints);
static void show_a_file( const char *filename, const int flags);
static void put_colored_text( const char *text, const int line_no,
               const int column, const int n_bytes, const int color);
int find_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
             const double r1, const double angle_param);   /* orb_func.cpp */
int search_for_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
              const double r1, double *angle_param);  /* orb_func.cpp */
void create_ades_file( const char *filename, const OBSERVE FAR *obs, int n_obs);
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile);
int generic_message_box( const char *message, const char *box_type);
int write_excluded_observations_file( const OBSERVE *obs, int n_obs);
int debug_printf( const char *format, ...)                 /* mpc_obs.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
int fetch_astrometry_from_mpc( FILE *ofile, const char *desig);
static void get_mouse_data( int *mouse_x, int *mouse_y, int *mouse_z, mmask_t *button);
int make_pseudo_mpec( const char *mpec_filename, const char *obj_name);
                                              /* ephem0.cpp */
int store_defaults( const ephem_option_t ephemeris_output_options,
         const int element_format, const int element_precision,
         const double max_residual_for_filtering,
         const double noise_in_sigmas);           /* elem_out.cpp */
int get_defaults( ephem_option_t *ephemeris_output_options, int *element_format,
         int *element_precision, double *max_residual_for_filtering,
         double *noise_in_sigmas);                /* elem_out.cpp */
int text_search_and_replace( char FAR *str, const char *oldstr,
                                     const char *newstr);   /* ephem0.cpp */
int sort_obs_by_date_and_remove_duplicates( OBSERVE *obs, const int n_obs);
double utc_from_td( const double jdt, double *delta_t);     /* ephem0.cpp */
void fix_home_dir( char *filename);                /* ephem0.cpp */
int write_environment_pointers( void);             /* mpc_obs.cpp */
int add_ephemeris_details( FILE *ofile, const double start_jd,  /* ephem0.c */
                                               const double end_jd);
void set_distance( OBSERVE FAR *obs, double r);             /* orb_func.c */
void set_statistical_ranging( const int new_using_sr);      /* elem_out.cpp */
int link_arcs( OBSERVE *obs, int n_obs, const double r1, const double r2);
int find_circular_orbits( OBSERVE FAR *obs1, OBSERVE FAR *obs2,
               double *orbit, const int desired_soln);   /* orb_fun2.cpp */
void set_up_observation( OBSERVE FAR *obs);               /* mpc_obs.cpp */
char *get_file_name( char *filename, const char *template_file_name);
double euler_function( const OBSERVE FAR *obs1, const OBSERVE FAR *obs2);
int find_relative_orbit( const double jd, const double *ivect,
               ELEMENTS *elements, const int ref_planet);     /* runge.cpp */
int find_parabolic_orbit( OBSERVE FAR *obs, const int n_obs,
            double *orbit, const int direction);         /* orb_func.cpp */
int format_jpl_ephemeris_info( char *buff);
double improve_along_lov( double *orbit, const double epoch, const double *lov,
          const unsigned n_params, unsigned n_obs, OBSERVE *obs);
bool is_topocentric_mpc_code( const char *mpc_code);
int64_t nanoseconds_since_1970( void);                      /* mpc_obs.c */
int metropolis_search( OBSERVE *obs, const int n_obs, double *orbit,
               const double epoch, int n_iterations, double scale);
const char *get_find_orb_text( const int index);
int set_tholen_style_sigmas( OBSERVE *obs, const char *buff);  /* mpc_obs.c */
FILE *fopen_ext( const char *filename, const char *permits);   /* miscell.cpp */
int find_vaisala_orbit( double *orbit, const OBSERVE *obs1,   /* orb_func.c */
                     const OBSERVE *obs2, const double solar_r);
int extended_orbit_fit( double *orbit, OBSERVE *obs, int n_obs,
                  const unsigned fit_type, double epoch);     /* orb_func.c */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
int load_environment_file( const char *filename);          /* mpc_obs.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);
int orbital_monte_carlo( const double *orbit, OBSERVE *obs, const int n_obs,
         const double curr_epoch, const double epoch_shown);   /* orb_func.cpp */
char *make_config_dir_name( char *oname, const char *iname);    /* miscell.cpp */
int reset_astrometry_filename( int *argc, const char **argv);
int set_language( const int language);                      /* elem_out.cpp */
static void show_splash_screen( void);
void shellsort_r( void *base, const size_t n_elements, const size_t esize,
         int (*compare)(const void *, const void *, void *), void *context);
static int count_wide_chars_in_utf8_string( const char *iptr, const char *endptr);
char **load_file_into_memory( const char *filename, size_t *n_lines,
                        const bool fail_if_not_found);      /* mpc_obs.cpp */
void make_observatory_info_text( char *text, const size_t textlen,
             const OBSERVE *obs, int n_obs, const char *mpc_code);
void size_from_h_text( const double abs_mag, char *obuff,
                                 const int obuff_size);  /* ephem0.c */
int find_fcct_biases( const double ra, const double dec, const char catalog,
                 const double jd, double *bias_ra, double *bias_dec);
int select_tracklet( OBSERVE *obs, const int n_obs, const int idx);
int get_orbit_from_mpcorb_sof( const char *object_name, double *orbit,
             ELEMENTS *elems, const double full_arc_len, double *max_resid);
int improve_sr_orbits( sr_orbit_t *orbits, OBSERVE FAR *obs,
               const unsigned n_obs, const unsigned n_orbits,  /* orb_func.c */
               const double noise_in_sigmas, const int writing_sr_elems);
int save_ephemeris_settings( const ephem_option_t ephemeris_output_options,
      const int n_steps, const char *obscode, const char *step_size,
      const char *ephem_start, const char *config);      /* elem_out.cpp */
int load_ephemeris_settings( ephem_option_t *ephemeris_output_options,
      int *n_steps, char *obscode, char *step_size, char *ephem_start,
      const char *config);                               /* elem_out.cpp */
void compute_effective_solar_multiplier( const char *constraints);  /* runge.c */

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

int getnstr_ex( char *str, int *loc, int maxlen, const int size);  /* getstrex.c */

#ifdef __cplusplus
}
#endif  /* #ifdef __cplusplus */

double comet_g_func( const long double r);                   /* runge.cpp */

extern int n_orbit_params;
extern double maximum_jd, minimum_jd;        /* orb_func.cpp */

#define COLOR_BACKGROUND            1
#define COLOR_ORBITAL_ELEMENTS      2
#define COLOR_FINAL_LINE            3
#define COLOR_SELECTED_OBS          4
#define COLOR_HIGHLIT_BUTTON        5
#define COLOR_EXCLUDED_OBS          6
#define COLOR_OBS_INFO              7
#define COLOR_MESSAGE_TO_USER       8
#define COLOR_RESIDUAL_LEGEND       9
#define COLOR_MENU                 10
#define COLOR_SCROLL_BAR           11
/* #define COLOR_DEFAULT_INQUIRY      12   Defined in 'mpc_obs.h' */
/* #define COLOR_ATTENTION            13   Defined in 'mpc_obs.h' */
#define COLOR_MPC_CODES            14
/* colors 15, 16, 17, and 18 also used for MPC codes;  see 'command.txt' */

#define SHOW_FILE_IS_EPHEM          1
#define SHOW_FILE_IS_CALENDAR       2

int curses_kbhit( )
{
   int c;

   nodelay( stdscr, TRUE);
   c = getch( );
   nodelay( stdscr, FALSE);
   if( c != ERR)     /* no key waiting */
      ungetch( c);
   return( c);
}

int curses_kbhit_without_mouse( )
{
   const int rval = curses_kbhit( );

   return( rval != KEY_MOUSE ? rval : 0);
}

static int extended_getch( void)
{
#ifdef _WIN32
   int rval = getch( );

   if( !rval)
      rval = 256 + getch( );
#else
   int rval = getch( );

   if( rval == 27)
      {
      nodelay( stdscr, TRUE);
      rval = getch( );
      nodelay( stdscr, FALSE);
      if( rval == ERR)    /* no second key found */
         rval = 27;       /* just a plain ol' Escape */
      else
         rval += (ALT_A - 'a');
      }
#endif
   return( rval);
}

static int *store_curr_screen( void)
{
   const int xsize = getmaxx( stdscr), ysize = getmaxy( stdscr);
   int x, y;
   int *rval = (int *)malloc( xsize * ysize * sizeof( chtype)
                    + 2 * sizeof( int));
   chtype *cptr = (chtype *)( rval + 2);

   rval[0] = xsize;
   rval[1] = ysize;
   for( y = 0; y < ysize; y++)
      for( x = 0; x < xsize; x++)
         *cptr++ = mvinch( y, x);
   return( rval);
}

#ifndef __PDCURSES__
/* Some ncurses platforms are squirrelly about how they handle
mouse movements.  Console commands have to be issued to turn
the mouse on or off.  This is still getting some testing,  and
is commented out by default. */

#ifdef MOUSE_MOVEMENT_EVENTS_ENABLED
   #define VT_IGNORE_ALL_MOUSE      CSI "?1003l\n"
   #define VT_RECEIVE_ALL_MOUSE     CSI "?1003h\n"
#endif
#endif

#ifdef __WATCOMC__
#undef endwin
extern "C" {
PDCEX  int     endwin_u64_4302(void);
}
#define endwin endwin_u64_4302
#endif

static int full_endwin( void)
{
#ifdef VT_IGNORE_ALL_MOUSE
   printf( VT_IGNORE_ALL_MOUSE);
#endif
   return( endwin( ));
}

#ifndef _WIN32
static int restart_curses( void)
{
#ifdef VT_RECEIVE_ALL_MOUSE
   printf( VT_RECEIVE_ALL_MOUSE);
#endif
   return( refresh( ));
}
#endif

static void restore_screen( const int *screen)
{
   const int xsize = getmaxx( stdscr), ysize = getmaxy( stdscr);
   int y;
   const int n_out = (xsize > screen[0] ? screen[0] : xsize);
   const chtype *cptr = (const chtype *)( screen + 2);

   clear( );
   for( y = 0; y < ysize && y < screen[1]; y++, cptr += screen[0])
      mvaddchnstr( y, 0, cptr, n_out);
}

int clipboard_to_file( const char *filename, const int append,
                           const bool use_selection); /* clipfunc.cpp */
int copy_file_to_clipboard( const char *filename);    /* clipfunc.cpp */

static bool curses_running = false;

static const char *help_file_name = NULL;
static int mpc_code_select = 0;

#define HINT_TEXT -1

static int full_inquire( const char *prompt, char *buff, const int max_len,
                     const int color, const int line0, const int col0)
{
   int rval = -1;

   if( !curses_running)    /* for error messages either before initscr() */
      {                    /* or after endwin( )                         */
      printf( "%s", prompt);
      getchar( );
      return( 0);
      }

   while( rval < 0)
      {
      int i, j, n_lines = 1, box_size = (buff ? 14 : 0);
      const int side_borders = 1;   /* leave a blank on either side */
      int real_width, line = line0, col = col0;
      char tbuff[200];
      int *buffered_screen;
      int x, y, z;
      mmask_t button;

      i = 0;
      while( prompt[i])
         {
         int new_size;

         j = i;
         while( prompt[j] && prompt[j] != '\n')
            j++;
         new_size = count_wide_chars_in_utf8_string(
                        prompt + i, prompt + j);
         if( help_file_name && !i)
            new_size += 4;       /* ensure room on top line for [?] */
         if( box_size < new_size)
            box_size = new_size;
         i = j;
         if( prompt[i])   /* skip trailing '\n's */
            {
            n_lines++;
            i++;
            }
         }
      if( box_size > getmaxx( stdscr) - 2)
         box_size = getmaxx( stdscr) - 2;

      real_width = side_borders * 2 + box_size;
      if( line == -1)         /* just center the box */
         {
         line = (getmaxy( stdscr) - n_lines) / 2;
         col = (getmaxx( stdscr) - box_size) / 2;
         }
      else              /* pop-up; may move above or below */
         {
         if( line + n_lines >= getmaxy( stdscr))
            line -= n_lines;
         else
            line++;
         if( col + real_width >= getmaxx( stdscr))
            {
            col -= real_width;
            if( col < 1)
               col = 1;
            }
         else
            col++;
         }
      assert( n_lines > 0 && n_lines < 200);
      assert( real_width > 0 && real_width < (int)sizeof( tbuff));
      tbuff[real_width] = '\0';
      buffered_screen = store_curr_screen( );
      col -= side_borders;
      for( i = 0; prompt[i]; )
         {
         int n_spaces, color_to_use = color;
         int n_wchars;

         for( j = i; prompt[j] && prompt[j] != '\n'; j++)
            ;
         memset( tbuff, ' ', side_borders);
         memcpy( tbuff + side_borders, prompt + i, j - i);
         n_wchars = count_wide_chars_in_utf8_string(
                        prompt + i, prompt + j);
         n_spaces = box_size + side_borders - n_wchars;
         if( n_spaces > 0)
            memset( tbuff + side_borders + j - i, ' ', n_spaces);
         else
            n_spaces = 0;
         if( !i)
            color_to_use |= A_OVERLINE;
         if( !prompt[j] || !prompt[j + 1])
            color_to_use |= A_UNDERLINE;
         put_colored_text( tbuff, line, col,
                side_borders + j - i + n_spaces, color_to_use);
         if( !i && help_file_name)
            put_colored_text( "[?]", line, col + real_width - 4,
                             3, color_to_use | A_REVERSE);
         i = j;
         if( prompt[i] == '\n')
            i++;
         line++;
         }
      if( buff)         /* we're asking for text from the user */
         {
         int loc = 0;

         memset( tbuff, ' ', real_width);
         *buff = '\0';
         do
            {
            put_colored_text( tbuff, line + 1, col, real_width, color);
            put_colored_text( tbuff, line + 2, col, real_width, color);
            put_colored_text( "[OK]",
                 line + 2, col + side_borders, 4, color | A_REVERSE);
            put_colored_text( "[Cancel]",
                 line + 2, col + real_width - side_borders - 8, 8, color | A_REVERSE);
            put_colored_text( tbuff, line, col, real_width, color);
            move( line, col + side_borders);
            attrset( COLOR_PAIR( COLOR_BACKGROUND));
            rval = getnstr_ex( buff, &loc, max_len, box_size);
            if( rval == KEY_MOUSE)
               {
               get_mouse_data( &x, &y, &z, &button);
               if( button & BUTTON1_CLICKED)
                  {
                  x -= col;
                  if( y == line - n_lines && x >= real_width - 4 && x < real_width - 1)
                     rval = KEY_F( 1);
                  if( y == line + 2 && x >= side_borders && x < real_width - 1)
                     {
                     if( x < side_borders + 4)
                        rval = 0;         /* OK clicked */
                     else if( x >= real_width - 8 - side_borders)
                        rval = 27;        /* Cancel clicked */
                     }
                  }
               }
            if( rval == KEY_F( 1) && help_file_name)
               {
               int *buffered_screen_2 = store_curr_screen( );

               show_a_file( help_file_name, 0);
               restore_screen( buffered_screen_2);
               free( buffered_screen_2);
               }
            }
            while( rval > 0 && rval != 27);
         }
      else        /* we just want the user to pick a line */
         {
         int highlit_line = -1;     /* initially,  no line is highlit */
         int highlit_x = 0;
         bool show_help = false;

         curs_set( 0);        /* turn cursor off */
         if( max_len == HINT_TEXT)  /* with hint text,  we quit as soon as anything happens */
            {
            while( curses_kbhit( ) == ERR)
               napms( 50);
            rval = '?';
            }
         else do
            {
            rval = extended_getch( );
            if( rval == KEY_RESIZE)
               {
               rval = -1;
               resize_term( 0, 0);
               }
                    /* If you click within the inquiry box borders, */
                    /* you get KEY_F(1) on the top line, KEY_F(2)   */
                    /* on the second line,  etc.                    */
            if( rval == KEY_MOUSE)
               {
               int curr_line, pass;

               get_mouse_data( &x, &y, &z, &button);
               curr_line = y;
               x -= col;
               y -= line - n_lines;
               if( button & BUTTON4_PRESSED)   /* actually 'wheel up' */
                  rval = KEY_UP;
               else if( button5_pressed)       /* actually 'wheel down' */
                  rval = KEY_DOWN;
               else if( y >= 0 && y < n_lines && x >= 0 && x < real_width)
                  {
                  if( mpc_code_select)
                     {
                     int selection;

                     if( x == real_width - 1)
                        x--;
                     selection = x / 4 + y * (box_size / 4 + 1);
                     if( y == n_lines - 1)
                        rval = 27;
                     else if( selection < mpc_code_select)
                        rval = KEY_F( 1 + selection);
                     else
                        curr_line = -1;
                     }
                  else
                     rval = KEY_F( y + 1);
                  }
               else
                  curr_line = -1;
               if( button & (REPORT_MOUSE_POSITION | BUTTON_PRESS_EVENT))
                  rval = KEY_MOUSE;          /* ignore mouse moves */
               if( curr_line != highlit_line || /* move the highlight */
                           (mpc_code_select && x / 4 != highlit_x / 4))
                  for( pass = 0; pass < 2; pass++)
                     {
                     if( highlit_line != -1)
                        {
                        const attr_t attr = (pass ? A_REVERSE : A_NORMAL);

                        if( mpc_code_select && highlit_line < line - 1)
                           mvchgat( highlit_line, col + highlit_x - highlit_x % 4,
                                               5, attr, color, NULL);
                        else
                           mvchgat( highlit_line, col, real_width, attr, color,
                                                     NULL);
                        if( highlit_line == line - n_lines && help_file_name)
                           mvchgat( highlit_line, col + real_width - 4, 3,
                                 attr ^ (A_REVERSE | A_NORMAL), color, NULL);
                        }
                     highlit_line = curr_line;
                     highlit_x = x;
                     }
               if( !y && x >= real_width - 4 && x < real_width - 1
                      && help_file_name && rval == KEY_F( 1))
                  show_help = true;
               }
            }
            while( rval == KEY_MOUSE);
         curs_set( 1);        /* turn cursor back on */
         if( rval == '?' && help_file_name)
            show_help = true;
         if( show_help)
            {
            show_a_file( help_file_name, 0);
            rval = -1;
            }
         }
      restore_screen( buffered_screen);
      free( buffered_screen);
      if( max_len != HINT_TEXT)
         flushinp( );
      refresh( );
      }
   help_file_name = NULL;
   return( rval);
}

int inquire( const char *prompt, char *buff, const int max_len,
                     const int color)
{
   extern char *mpec_error_message;

   assert( prompt);
   if( !mpec_error_message && color == COLOR_ATTENTION)
      {
      const size_t len = strlen( prompt) + 1;

      mpec_error_message = (char *)malloc( len);
      strlcpy_err( mpec_error_message, prompt, len);
      }
   return( full_inquire( prompt, buff, max_len, color, -1, -1));
}

static int select_mpc_code( const OBSERVE *obs, const int n_obs, int curr_obs)
{
   const int max_n_codes = 500;
   const int buffsize = max_n_codes * 4 + 70;
   char *buff = (char *)malloc( buffsize);
   int n_codes = 0, i, nx = 0, c;

   *buff = '\0';
   for( i = 0; i < n_obs && n_codes < max_n_codes; i++)
      if( !strstr( buff, obs[i].mpc_code))
         {
         int j = 0;
         char *tptr = buff;

         while( j < n_codes && memcmp( tptr, obs[i].mpc_code, 4) < 0)
            {
            j++;
            tptr += 4;
            }
         memmove( tptr + 4, tptr, (n_codes - j) * 4 + 1);
         memcpy( tptr, obs[i].mpc_code, 3);
         tptr[3] = ' ';
         n_codes++;
         if( nx * nx < n_codes && (nx + 1) * 4 < COLS - 3)
            nx++;
         }
   for( i = nx; i < n_codes; i += nx)
      buff[i * 4 - 1] = '\n';
   buff[n_codes * 4 - 1] = '\0';
   strlcat_err( buff, "\nCancel", buffsize);
   mpc_code_select = n_codes;
   c = inquire( buff, NULL, 0, COLOR_DEFAULT_INQUIRY);
   mpc_code_select = 0;
   c -= KEY_F( 1);
   if( c >= 0 && c < n_codes)
      {
      char *tptr = buff + c * 4;

      for( i = 0; i < n_obs; i++)
         {
         curr_obs = (curr_obs + 1) % n_obs;
         if( !memcmp( obs[curr_obs].mpc_code, tptr, 3))
            break;
         }
      }
   free( buff);
   return( curr_obs);
}

/* In the (interactive) console Find_Orb,  these allow some functions
in orb_func.cpp to show info as orbits are being computed.  In the
non-interactive 'fo' code,  they're mapped to do nothing (see 'fo.cpp'). */

void refresh_console( void)
{
   refresh( );
}

void move_add_nstr( const int col, const int row, const char *msg, const int n_bytes)
{
   attrset( COLOR_PAIR( COLOR_FINAL_LINE));
   mvaddnstr( col, row, msg, n_bytes);
}

double current_jd( void);                       /* elem_out.cpp */

static int extract_date( const char *buff, double *jd)
{
   int rval = 0, is_ut;

               /* If the date seems spurious,  use 'now' as our zero point: */
   if( *jd < minimum_jd || *jd > maximum_jd || *jd == -.5 || *jd == 0.)
      *jd = current_jd( );
   *jd = get_time_from_string( *jd, buff,
                     FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, &is_ut);
   rval = 2;
   if( *jd == 0. || is_ut < 0)
      rval = -1;
   return( rval);
}

void compute_variant_orbit( double *variant, const double *ref_orbit,
                     const double n_sigmas);       /* orb_func.cpp */

static double *set_up_alt_orbits( const double *orbit, unsigned *n_orbits)
{
   extern int available_sigmas;
   extern double *sr_orbits;
   extern unsigned n_sr_orbits;

   switch( available_sigmas)
      {
      case COVARIANCE_AVAILABLE:
         memcpy( sr_orbits, orbit, MAX_N_PARAMS * sizeof( double));
         compute_variant_orbit( sr_orbits + n_orbit_params, sr_orbits, 1.);
         *n_orbits = 2;
         break;
      case SR_SIGMAS_AVAILABLE:
         *n_orbits = n_sr_orbits;
         break;
      default:
         memcpy( sr_orbits, orbit, MAX_N_PARAMS * sizeof( double));
         *n_orbits = 1;
         break;
      }
   return( sr_orbits);
}

/* Outputs residuals in the 'short',  three-columns-wide format
used in MPECs and pseudo-MPECs.   */

static void create_resid_file( const OBSERVE *obs, const int n_obs,
         const char *input_filename, int residual_format)
{
   char buff[260];
   extern const char *residual_filename;

   residual_format &= (RESIDUAL_FORMAT_TIME_RESIDS | RESIDUAL_FORMAT_MAG_RESIDS
               | RESIDUAL_FORMAT_NORMALIZED
               | RESIDUAL_FORMAT_PRECISE | RESIDUAL_FORMAT_OVERPRECISE);
   residual_format |= RESIDUAL_FORMAT_SHORT;
   write_residuals_to_file( get_file_name( buff, residual_filename),
                    input_filename, n_obs, obs, residual_format);
}

static void set_ra_dec_format( void)
{
   char buff[1000];
   const char *lines[60];
   const char *iline = get_find_orb_text( 2054), *tptr = iline;
   const char *ra_dec_fmt = "RA_DEC_FORMAT";
   const char *curr_format = get_environment_ptr( ra_dec_fmt);
   size_t len = 0, n_lines = 0;
   int c;

   while( *tptr)
      {
      const char *tptr2 = strchr( tptr, ':');

      assert( tptr2);
      snprintf( buff + len, sizeof( buff) - len, "%d ( ) ", (int)n_lines + 1);
      assert( n_lines < sizeof( lines) / sizeof( lines[0]));
      lines[n_lines++] = tptr;
      if( strlen( curr_format) == (size_t)( tptr2 - tptr) &&
                  !memcmp( curr_format, tptr, tptr2 - tptr))
         buff[len + 3] = '*';
      len += 6;
      tptr2++;
      while( *tptr2 >= ' ')
         buff[len++] = *tptr2++;
      if( *tptr2)
         buff[len++] = *tptr2++;
      tptr = tptr2;
      }
   buff[len] = '\0';
   help_file_name = "radecfmt.txt";
   c = inquire( buff, NULL, 0, COLOR_DEFAULT_INQUIRY);
   if( c >= '0' && c < '1' + (int)n_lines)
      c += KEY_F( 1) - '1';
   if( c >= KEY_F( 1) && c <= (int)KEY_F( n_lines))
      {
      const char *tptr2;

      tptr = lines[c - KEY_F( 1)];
      tptr2 = strchr( tptr, ':');
      assert( tptr2);
      memcpy( buff, tptr, tptr2 - tptr);
      buff[tptr2 - tptr] = '\0';
      set_environment_ptr( ra_dec_fmt, buff);
      }
}

static void select_angular_motion_units( void)
{
   char buff[1000], *tptr = buff, curr_units[4];
   int c;

   strlcpy( curr_units, get_environment_ptr( "MOTION_UNITS"),
                     sizeof( curr_units));
   if( !*curr_units)
      strlcpy( curr_units, "'/h", sizeof( curr_units));
   strlcpy_error( buff, get_find_orb_text( 2078));
   if( (tptr = strstr( buff, curr_units)) != NULL)
      {                    /* mark currently selected units */
      while( *tptr != '(')
         tptr--;
      tptr[1] = 'o';
      }
   c = inquire( buff, NULL, 0, COLOR_DEFAULT_INQUIRY);
   if( c >= '0' && c < '9')
      c += KEY_F( 1) - '1';
   c -= KEY_F( 0);
   if( c >= 0)
      {
      tptr = buff;
      while( *tptr && c)
         {
         while( *tptr && *tptr != '\n')
            tptr++;
         if( *tptr == '\n')
            tptr++;
         c--;
         }
      if( *tptr)
         {
         size_t i = 0;

         tptr += 6;
         while( tptr[i] > ' ')
            i++;
         tptr[i] = '\0';
         set_environment_ptr( "MOTION_UNITS", tptr);
         }
      }
}

/* For 'radio button' dialogs such as the comet model,  resid type,  and
ephemeris type (see 'efindorb.txt'),  this sets the appropriate pseudo-radio
button. */

static char *_set_radio_button( char *text, const int option_num)
{
   char tbuff[10], *line_ptr;

   snprintf_err( tbuff, sizeof( tbuff), "%d ( )", option_num);
   line_ptr = strstr( text, tbuff);
   assert( line_ptr);
   if( line_ptr)
      line_ptr[3] = '*';
   return( line_ptr);
}

static void show_calendar( void)
{
   FILE *ifile = fopen_ext( "calend.txt", "clrb");

   if( ifile)
      fclose( ifile);
   show_a_file( (ifile ? "calend.txt" : "calendar.txt"), SHOW_FILE_IS_CALENDAR);
}

/* Here's a simplified example of the use of the 'ephemeris_in_a_file'
   function... nothing fancy,  but it shows how it's used.  */

static char mpc_code[80];
static char ephemeris_start[80], ephemeris_step_size[300];
static int n_ephemeris_steps;
static ephem_option_t ephemeris_output_options;

static void create_ephemeris( const double *orbit, const double epoch_jd,
         OBSERVE *obs, const int n_obs, const char *obj_name,
         const char *input_filename, const int residual_format)
{
   int c = 1;
   char buff[2000];
   double jd_start = 0., jd_end = 0., step = 0.;
   bool show_advanced_options = false;
   int offset_line = 0;
   char message_to_user[180];

   *message_to_user = '\0';
   while( c > 0)
      {
      int format_start;
      unsigned i, n_lines;
      int vect_frame = -1, planet_idx;
      double vect_dist_units = 0., vect_time_units = 0.;
      const int ephem_type = (int)(ephemeris_output_options & 7);
      bool reset_vect_units = false;
      extern double ephemeris_mag_limit;
      const char *tptr, *err_msg = NULL;
      char *end_of_location_text;
      char vect_epoch[9];
      const bool is_topocentric =
               is_topocentric_mpc_code( mpc_code);

      put_colored_text( message_to_user, getmaxy( stdscr) - 1, 0, -1,
                        COLOR_MESSAGE_TO_USER);
      *message_to_user = '\0';
      jd_start = 0.;
      format_start = extract_date( ephemeris_start, &jd_start);
      step = get_step_size( ephemeris_step_size, NULL, NULL);
      if( format_start == 1 || format_start == 2)
         {
         if( step && format_start == 1)  /* time was relative to 'right now' */
            jd_start = floor( (jd_start - .5) / step) * step + .5;
         snprintf( buff, sizeof( buff), " (Ephem start: JD %.5f = ", jd_start);
         full_ctime( buff + strlen( buff), jd_start,
                    FULL_CTIME_DAY_OF_WEEK_FIRST | CALENDAR_JULIAN_GREGORIAN);
         strcat( buff, ")\n");
         jd_end = jd_start + step * (double)n_ephemeris_steps;
         snprintf_append( buff, sizeof( buff), " (Ephem end:   JD %.5f = ", jd_end);
         full_ctime( buff + strlen( buff), jd_end,
                    FULL_CTIME_DAY_OF_WEEK_FIRST | CALENDAR_JULIAN_GREGORIAN);
         strcat( buff, ")\n");
         }
      else
         {
         strlcpy_error( buff, "(Ephemeris starting time isn't valid)\n");
         jd_start = jd_end = 0.;
         }

      snprintf_append( buff, sizeof( buff), "T  Ephem start: %s\n", ephemeris_start);
      snprintf_append( buff, sizeof( buff), "N  Number steps: %d\n",
                                        n_ephemeris_steps);
      snprintf_append( buff, sizeof( buff) , "S  Step size: %.60s\n", ephemeris_step_size);
      snprintf_append( buff, sizeof( buff), "L  Location: (%s) ", mpc_code);
      planet_idx = put_observer_data_in_text( mpc_code, buff + strlen( buff));
      strcat( buff, "\n");
      end_of_location_text = buff + strlen( buff);
      if( ephem_type == OPTION_STATE_VECTOR_OUTPUT
                   || ephem_type == OPTION_POSITION_OUTPUT)
         {
         const char *vect_opts = get_environment_ptr( "VECTOR_OPTS");
         const char *otext;
         int end_loc;

         sscanf( vect_opts, "%d,%lf,%lf%n", &vect_frame, &vect_dist_units,
                        &vect_time_units, &end_loc);
         snprintf_append( buff, sizeof( buff), "F  %s\n",
                     vect_frame ? "Ecliptic frame" : "Equatorial frame");
         if( vect_opts[end_loc] == ',')
            strlcpy_error( vect_epoch, vect_opts + end_loc + 1);
         else
            strlcpy_error( vect_epoch, "2000");
         if( *vect_epoch == 'm')
            strlcat_error( buff, "P  Mean coords of date\n");
         else if( *vect_epoch == 't')
            strlcat_error( buff, "P  True coords of date\n");
         else
            {
            const double epoch_year = atof( vect_epoch);

            assert( epoch_year > 1700. && epoch_year < 2300.);
            snprintf_append( buff, sizeof( buff), "P  Coord epoch %s\n",
                                       vect_epoch);
            }
         if( vect_dist_units == 1.)
            otext = "AU";
         else if( fabs( vect_dist_units - AU_IN_KM) < 1.)  /* allow for roundoff */
            {
            otext = "km";
            vect_dist_units = AU_IN_KM;
            }
         else
            {
            otext = "Earth radii";
            vect_dist_units = AU_IN_KM / EARTH_RADIUS_IN_KM;
            }
         snprintf_append( buff, sizeof( buff), "D  Distances in %s\n", otext);
         if( ephem_type == OPTION_STATE_VECTOR_OUTPUT)
            {
            if( vect_time_units == 1.)
               otext = "days";
            else if( vect_time_units == 24.)
               otext = "hours";
            else if( vect_time_units == 1440.)
               otext = "minutes";
            else
               {
               vect_time_units = 86400.;
               otext = "seconds";
               }
            snprintf_append( buff, sizeof( buff), "U  Time units %s\n", otext);
            }
         }
      if( ephem_type == OPTION_OBSERVABLES)    /* for other tables,        */
         {                          /* these options are irrelevant:       */
         snprintf_append( buff, sizeof( buff), "Z [%c] Total motion & position angle\n",
                  (ephemeris_output_options & OPTION_MOTION_OUTPUT) ? '*' : ' ');
         snprintf_append( buff, sizeof( buff), "O [%c] Motions in RA and dec\n",
                  (ephemeris_output_options & OPTION_SEPARATE_MOTIONS) ? '*' : ' ');
         if( is_topocentric)
            snprintf_append( buff, sizeof( buff), "A [%c] Alt/az info\n",
                  (ephemeris_output_options & OPTION_ALT_AZ_OUTPUT) ? '*' : ' ');
         snprintf_append( buff, sizeof( buff), "R [%c] Radial velocity\n",
                  (ephemeris_output_options & OPTION_RADIAL_VEL_OUTPUT) ? '*' : ' ');
         snprintf_append( buff, sizeof( buff), "P [%c] Phase angle\n",
                  (ephemeris_output_options & OPTION_PHASE_ANGLE_OUTPUT) ? '*' : ' ');
         snprintf_append( buff, sizeof( buff), "V [%c] Visibility indicator\n",
               (ephemeris_output_options & OPTION_VISIBILITY) ? '*' : ' ');
         if( is_topocentric)
            snprintf_append( buff, sizeof( buff), "U [%c] Suppress unobservables\n",
                  (ephemeris_output_options & OPTION_SUPPRESS_UNOBSERVABLE) ? '*' : ' ');
         snprintf_append( buff, sizeof( buff), "F Suppress when fainter than mag: %.1f\n",
                  ephemeris_mag_limit);
         snprintf_append( buff, sizeof( buff), "D [%c] Positional sigmas\n",
                  (ephemeris_output_options & OPTION_SHOW_SIGMAS) ? '*' : ' ');
         snprintf_append( buff, sizeof( buff), "0 [%c] Show advanced options\n",
                  show_advanced_options ? '*' : ' ');
         if( show_advanced_options)
            {
            snprintf_append( buff, sizeof( buff), "! [%c] Show explanations at end of ephems\n",
                  (ephemeris_output_options & OPTION_EXPLANATIONS) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "B [%c] Phase angle bisector\n",
                  (ephemeris_output_options & OPTION_PHASE_ANGLE_BISECTOR) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "H [%c] Heliocentric ecliptic\n",
                  (ephemeris_output_options & OPTION_HELIO_ECLIPTIC) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "X [%c] Topocentric ecliptic\n",
                  (ephemeris_output_options & OPTION_TOPO_ECLIPTIC) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "G [%c] Ground track\n",
                  (ephemeris_output_options & OPTION_GROUND_TRACK) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "J [%c] Lunar elongation\n",
                  (ephemeris_output_options & OPTION_LUNAR_ELONGATION) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "Y [%c] Computer-friendly output\n",
                  (ephemeris_output_options & OPTION_COMPUTER_FRIENDLY) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "W [%c] Round to nearest step\n",
                  (ephemeris_output_options & OPTION_ROUND_TO_NEAREST_STEP) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "I [%c] Space velocity\n",
                  (ephemeris_output_options & OPTION_SPACE_VEL_OUTPUT) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "1 [%c] RA/decs\n",
                  (ephemeris_output_options & OPTION_SUPPRESS_RA_DEC) ? ' ' : '*');
            snprintf_append( buff, sizeof( buff), "2 [%c] delta (dist from observer)\n",
                  (ephemeris_output_options & OPTION_SUPPRESS_DELTA) ? ' ' : '*');
            snprintf_append( buff, sizeof( buff), "3 [%c] r (dist from sun)\n",
                  (ephemeris_output_options & OPTION_SUPPRESS_SOLAR_R) ? ' ' : '*');
            snprintf_append( buff, sizeof( buff), "4 [%c] elong\n",
                  (ephemeris_output_options & OPTION_SUPPRESS_ELONG) ? ' ' : '*');
            snprintf_append( buff, sizeof( buff), "9 [%c] Sky brightness\n",
                  (ephemeris_output_options & OPTION_SKY_BRIGHTNESS) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "K [%c] Galactic coordinates\n",
                  (ephemeris_output_options & OPTION_GALACTIC_COORDS) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "%% [%c] Galactic confusion\n",
                  (ephemeris_output_options & OPTION_GALACTIC_CONFUSION) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "@ [%c] Comet options\n",
                  (ephemeris_output_options & OPTION_SUN_TARGET_PA) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "5 [%c] Sun altitude\n",
                  (ephemeris_output_options & OPTION_SUN_ALT) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "6 [%c] Sun azimuth\n",
                  (ephemeris_output_options & OPTION_SUN_AZ) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "7 [%c] Moon altitude\n",
                  (ephemeris_output_options & OPTION_MOON_ALT) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "8 [%c] Moon azimuth\n",
                  (ephemeris_output_options & OPTION_MOON_AZ) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "{ [%c] SNR\n",
                  (ephemeris_output_options & OPTION_SNR) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "} [%c] Exposure time\n",
                  (ephemeris_output_options & OPTION_EXPOSURE_TIME) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "= [%c] Constellation\n",
                  (ephemeris_output_options & OPTION_CONSTELLATION) ? '*' : ' ');
            snprintf_append( buff, sizeof( buff), "$ [%c] RVel, delta sigmas\n",
                  (ephemeris_output_options & OPTION_RV_AND_DELTA_SIGMAS) ? '*' : ' ');
            }
         }
      for( i = n_lines = 0; buff[i]; i++)
         if( buff[i] == '\n')
            n_lines++;
      if( n_lines > (unsigned)LINES - 3)
         {
         char *tptr1, *tptr = end_of_location_text;

         if( offset_line < 0)
            offset_line = 0;
         if( offset_line > (int)n_lines - (LINES - 4))
            offset_line = (int)n_lines - ( LINES - 4);

         for( i = 0; i < (unsigned)offset_line; i++)
            {
            tptr = strchr( tptr, '\n');
            assert( tptr);
            tptr++;
            }
         tptr1 = tptr;
         for( i = 0; i < (unsigned)( LINES - 10); i++)
            {
            tptr1 = strchr( tptr1, '\n');
            if( !tptr1)
               fprintf( stderr, "\ni %u offset_line %u n_lines %u LINES %d\n",
                           i, offset_line, n_lines, LINES);
            assert( tptr1);
            tptr1++;
            if( !*tptr1)
               break;
            }
         *tptr1 = '\0';
         memmove( end_of_location_text, tptr, tptr1 - tptr + 1);
         strlcat_error( buff, "Up/down arrows or mouse wheel for other options\n");
         }
      tptr = get_find_orb_text( 2064);
      i = ephem_type;
      while( i && *tptr)
         if( *tptr++ == '\n')
            i--;
      if( *tptr == '0' + ephem_type)
         {
         tptr += 6;
         i = 0;
         while( tptr[i] >= ' ')
            i++;
         }
      else
         {
         tptr = (char *)"Unknown";
         i = 7;
         }
      snprintf_append( buff, sizeof( buff), "C  Ephemeris type : %.*s\n", i, tptr);
      snprintf_append( buff, sizeof( buff), "M  Make ephemeris\n");
      snprintf_append( buff, sizeof( buff), "Q  Return to main display");
      n_lines += 4;
      help_file_name = "dosephem.txt";
      c = inquire( buff, NULL, 0, COLOR_DEFAULT_INQUIRY);
                     /* Convert mouse clicks inside the 'dialog box'     */
                     /* to the corresponding first letter on that line   */
                     /* (except for the first two lines,  which wouldn't */
                     /* work;  they start with spaces) :                 */
      if( c >= KEY_F( 3) && c <= (int)KEY_F( n_lines))
         {
         unsigned n = c - KEY_F( 1);

         for( i = 0; buff[i] && n; i++)
            if( buff[i] == '\n')
               n--;
         c = tolower( buff[i]);
         }
      if( c >= ALT_0 && c <= ALT_7)
         {
         ephemeris_output_options &= ~7;
         ephemeris_output_options |= (c - ALT_0);
         }
      else if( c == ALT_A)
         {                    /* loading configuration */
         c = inquire( "Hit letter of configuration to load",
                         NULL, 30, COLOR_MESSAGE_TO_USER);
         if( c >= 'a' && c <= 'z')
            {
            char tstr[2];
            int err_code;

            tstr[0] = (char)( c + 'A' - 'a');
            tstr[1] = '\0';
            snprintf( buff, sizeof( buff), "STORED_EPHEM_%c", c);
            err_code = load_ephemeris_settings( &ephemeris_output_options, &n_ephemeris_steps,
                           mpc_code, ephemeris_step_size, ephemeris_start, tstr);
            snprintf( message_to_user, sizeof( message_to_user),
                  (err_code ? "Preset %c not found" : "Preset %c loaded"), c);
            }
         }
      else if( c >= CTRL( 'A') && c <= CTRL( 'Z'))
         {                    /* saving configuration */
         char tstr[2];

         tstr[0] = (char)( c + 'A' - CTRL( 'A'));
         tstr[1] = '\0';
         save_ephemeris_settings( ephemeris_output_options, n_ephemeris_steps,
                        mpc_code, ephemeris_step_size, ephemeris_start, tstr);
         snprintf( message_to_user, sizeof( message_to_user),
                  "Preset %c set to current ephemeris options", c);
         }
      else switch( c)
         {
         case '0':
            show_advanced_options = !show_advanced_options;
            break;
         case '1':
            if( ephemeris_output_options & OPTION_SUPPRESS_RA_DEC)
               set_ra_dec_format( );
            ephemeris_output_options ^= OPTION_SUPPRESS_RA_DEC;
            break;
         case '2':
            ephemeris_output_options ^= OPTION_SUPPRESS_DELTA;
            break;
         case '3':
            ephemeris_output_options ^= OPTION_SUPPRESS_SOLAR_R;
            break;
         case '4':
            ephemeris_output_options ^= OPTION_SUPPRESS_ELONG;
            break;
         case '5':
            ephemeris_output_options ^= OPTION_SUN_ALT;
            break;
         case '6':
            ephemeris_output_options ^= OPTION_SUN_AZ;
            break;
         case '7':
            ephemeris_output_options ^= OPTION_MOON_ALT;
            break;
         case '8':
            ephemeris_output_options ^= OPTION_MOON_AZ;
            break;
         case '9':
            ephemeris_output_options ^= OPTION_SKY_BRIGHTNESS;
            break;
         case '{':
            ephemeris_output_options ^= OPTION_SNR;
            break;
         case '}':
            ephemeris_output_options ^= OPTION_EXPOSURE_TIME;
            break;
         case '!':
            ephemeris_output_options ^= OPTION_EXPLANATIONS;
            break;
         case 'a': case 'A':
            ephemeris_output_options ^= OPTION_ALT_AZ_OUTPUT;
            break;
         case 'b': case 'B':
            ephemeris_output_options ^= OPTION_PHASE_ANGLE_BISECTOR;
            break;
         case 'c':
            {
            strlcpy_error( buff, get_find_orb_text( 2064));
            _set_radio_button( buff, ephem_type);
            help_file_name = "eph_type.txt";
            c = inquire( buff, NULL, 0, COLOR_DEFAULT_INQUIRY);
            if( c >= KEY_F( 1))
               i = c - KEY_F( 1);
            else if( c >= '0' && c < '8')
               i = c - '0';
            if( i < 7)
               {
               ephemeris_output_options &= ~7;
               ephemeris_output_options |= i;
               }
            }
            break;
         case 'C':
            show_calendar( );
            break;
         case 'd': case 'D':
            if( vect_frame  > -1)
               {
               if( vect_dist_units == 1.)
                  vect_dist_units = AU_IN_KM;
               else if( vect_dist_units == AU_IN_KM)
                  vect_dist_units = AU_IN_KM / EARTH_RADIUS_IN_KM;
               else
                  vect_dist_units = 1.;
               reset_vect_units = true;
               }
            else
               ephemeris_output_options ^= OPTION_SHOW_SIGMAS;
            break;
         case 'e': case 'E': case KEY_F( 2):
            help_file_name = "timehelp.txt";
            if( !inquire( "Enter end of ephemeris (YYYY MM DD, or JD, or 'now'):",
                     buff, sizeof( buff), COLOR_MESSAGE_TO_USER))
               {
               format_start = extract_date( buff, &jd_end);
               if( (format_start == 1 || format_start == 2)
                               && jd_end > minimum_jd && jd_end < maximum_jd)
                  {
                  n_ephemeris_steps = (int)ceil( (jd_end - jd_start) / step);
                  if( n_ephemeris_steps < 0)
                     {
                     n_ephemeris_steps = 1 - n_ephemeris_steps;
                     if( *ephemeris_step_size == '-')   /* eliminate neg sign */
                        memmove( ephemeris_step_size, ephemeris_step_size + 1,
                                             strlen( ephemeris_step_size));
                     else
                        {                            /* insert neg sign */
                        memmove( ephemeris_step_size + 1, ephemeris_step_size,
                                             strlen( ephemeris_step_size) + 1);
                        *ephemeris_step_size = '-';
                        }
                     }
                  }
               else
                  err_msg = "Couldn't parse the time you entered!";
               }
            break;
         case 'f': case 'F':
            if( vect_frame > -1)
               {
               vect_frame ^= 1;
               reset_vect_units = true;
               }
            else if( !inquire( "Mag limit for ephemerides: ", buff, sizeof( buff),
                     COLOR_MESSAGE_TO_USER) && atof( buff))
               ephemeris_mag_limit = atof( buff);
            break;
         case 'g': case 'G':
            ephemeris_output_options ^= OPTION_GROUND_TRACK;
            break;
         case 'h': case 'H':
            ephemeris_output_options ^= OPTION_HELIO_ECLIPTIC;
            break;
         case 'i': case 'I':
            ephemeris_output_options ^= OPTION_SPACE_VEL_OUTPUT;
            break;
         case 'j': case 'J':
            ephemeris_output_options ^= OPTION_LUNAR_ELONGATION;
            break;
         case 'k': case 'K':
            ephemeris_output_options ^= OPTION_GALACTIC_COORDS;
            break;
         case 'l': case 'L':
            if( !inquire( "Enter MPC code: ", buff, sizeof( buff), COLOR_MESSAGE_TO_USER))
               {
               if( strlen( buff) < 3)
                  err_msg = "MPC codes must be at least three characters long";
               else if( strlen( buff) < 5 || !memcmp( buff, "Ast", 3)
                           || !get_lat_lon_info( NULL, buff))
                  strlcpy_error( mpc_code, buff);
               else if( strlen( buff) > 4)
                  {
                  if( !get_observer_data( buff, buff, NULL))
                     {
                     buff[4] = '\0';
                     if( buff[3] == ' ')
                        buff[3] = '\0';
                     strlcpy_error( mpc_code, buff);
                     }
                  else
                     err_msg = "Didn't understand that observatory code";
                  }
               }
            break;
         case '@':        /* toggle comet-related options */
            ephemeris_output_options ^= OPTION_SUN_TARGET_PA |
                           OPTION_SUN_HELIO_VEL_PA | OPTION_ORBIT_PLANE_ANGLE;
            break;
         case ALT_D:
            strlcpy_error( ephemeris_start, "+0");
            ephemeris_output_options &= ~7;
                     /* FALLTHRU */
         case ALT_G:
            strlcpy_error( mpc_code, "500");
            break;
         case '%':
            ephemeris_output_options ^= OPTION_GALACTIC_CONFUSION;
            break;
         case ALT_L:
            strlcpy_error( mpc_code, "Lun");
            break;
         case ALT_O:                         /* set 'observables' */
            ephemeris_output_options &= ~7;
            break;
         case 'm': case 'M':
            {
            if( jd_start < minimum_jd || jd_start > maximum_jd)
               err_msg = "You need to set a valid starting date!";
            else if( !n_ephemeris_steps)
               err_msg = "You need to set the number of ephemeris steps!";
            else if( !step && *ephemeris_step_size != 'a')
               err_msg = "You need to set a valid step size!";
            else if( planet_idx == -1)
               err_msg = "You need to set a valid observer location!";
            else                 /* yes,  we can make an ephemeris */
               c = -2;
            }
            break;
         case 'n': case 'N': case KEY_F( 4):
            if( !inquire( "Number of steps:", buff, sizeof( buff), COLOR_MESSAGE_TO_USER)
                             && atoi( buff) > 0)
               n_ephemeris_steps = atoi( buff);
            break;
         case 'o': case 'O':
            ephemeris_output_options ^= OPTION_SEPARATE_MOTIONS;
            if( ephemeris_output_options & OPTION_SEPARATE_MOTIONS)
               select_angular_motion_units( );
            break;
         case 'p': case 'P':
            if( vect_frame > -1)
               {
               if( !inquire( "Epoch of vectors (t(rue), m(ean),  or a year) :",
                              buff, sizeof( vect_epoch), COLOR_MESSAGE_TO_USER))
                  {
                  if( (!buff[1] && (buff[0] == 't' || buff[0] == 'm'))
                        || (atof( buff) > 1700 && atof( buff) < 2300.))
                     {
                     strlcpy_error( vect_epoch, buff);
                     reset_vect_units = true;
                     }
                  else
                     err_msg = "Vector epoch must be t(rue), m(ean),  or a year!";
                  }
               }
            else
               ephemeris_output_options ^= OPTION_PHASE_ANGLE_OUTPUT;
            break;
         case 'r': case 'R':
            ephemeris_output_options ^= OPTION_RADIAL_VEL_OUTPUT;
            break;
         case 's': case 'S': case KEY_F( 5):
            strlcpy_error( buff, ephemeris_step_size);
            if( !inquire( "Enter step size in days: ",
                                      buff, sizeof( buff), COLOR_MESSAGE_TO_USER))
               {
               if( !strcmp( buff, "t"))
                  {
                  user_select_file( buff + 2, "Select ephemeris times file", 0);
                  if( !buff[2])        /* cancelled file selection */
                     *buff = '\0';
                  else
                     buff[1] = ' ';
                  }
               if( *buff)
                  strlcpy_error( ephemeris_step_size, buff);
               }
            break;
         case '-':         /* reverse ephemeris direction */
            {
            const size_t loc = (ephemeris_step_size[0] == '-' ? 1 : 0);

            memmove( ephemeris_step_size + 1 - loc, ephemeris_step_size + loc,
                           strlen( ephemeris_step_size) + 1);
            if( !loc)
               *ephemeris_step_size = '-';
            }
            break;
         case 't': case 'T': case KEY_F( 1): case KEY_F( 3):
            help_file_name = "timehelp.txt";
            strlcpy_err( buff, ephemeris_start, sizeof( ephemeris_start));
            if( !inquire( "Enter start of ephemeris (YYYY MM DD, or JD, or 'now'):",
                         buff, sizeof( ephemeris_start), COLOR_MESSAGE_TO_USER))
               {
               double jd = 0.;
               const int format = extract_date( buff, &jd);

               if( (format == 1 || format == 2)
                               && jd > minimum_jd && jd < maximum_jd)
                  strlcpy_err( ephemeris_start, buff, sizeof( ephemeris_start));
               else
                  err_msg = "Couldn't parse the time you entered!";
               }
            break;
         case ALT_T:
            strlcpy_error( ephemeris_start, "+0");
            break;
         case 'u': case 'U':
            if( vect_time_units)
               {
               if( vect_time_units == 1.)
                  vect_time_units = 24.;
               else if( vect_time_units == 24.)
                  vect_time_units = 24. * 60.;
               else if( vect_time_units == 24. * 60.)
                  vect_time_units = 24. * 60. * 60.;
               else
                  vect_time_units = 1;
               reset_vect_units = true;
               }
            else
               ephemeris_output_options ^= OPTION_SUPPRESS_UNOBSERVABLE;
            break;
         case 'v': case 'V':
            ephemeris_output_options ^= OPTION_VISIBILITY;
            break;
         case 'w': case 'W':
            ephemeris_output_options ^= OPTION_ROUND_TO_NEAREST_STEP;
            break;
         case 'x': case 'X':
            ephemeris_output_options ^= OPTION_TOPO_ECLIPTIC;
            break;
         case 'y': case 'Y':
            ephemeris_output_options ^= OPTION_COMPUTER_FRIENDLY;
            break;
         case 'z': case 'Z':
            ephemeris_output_options ^= OPTION_MOTION_OUTPUT;
            if( ephemeris_output_options & OPTION_MOTION_OUTPUT)
               select_angular_motion_units( );
            break;
         case '#':
            ephemeris_output_options ^= OPTION_MOIDS;
            break;
         case '=':
            ephemeris_output_options ^= OPTION_CONSTELLATION;
            break;
         case '$':
            ephemeris_output_options ^= OPTION_RV_AND_DELTA_SIGMAS;
            break;
         case 'q': case 'Q':
         case 27:
            c = -1;
            break;
#ifdef KEY_EXIT
         case KEY_EXIT:
            exit( 0);
            break;
#endif
         case KEY_UP:
            offset_line--;
            break;
         case KEY_DOWN:
            offset_line++;
            break;
         default:
         case '\'': case '"': case '<':
         case ',': case '.': case '>':
         case '^': case '&': case '*':
         case '(': case ')': case '[': case ']':
            show_a_file( "dosephem.txt", 0);
            break;
         }
      if( err_msg)
         inquire( err_msg, NULL, 0, COLOR_FINAL_LINE);
      if( reset_vect_units)
         {
         snprintf( buff, sizeof( buff), "%d,%f,%f,%s", vect_frame,
                     vect_dist_units, vect_time_units, vect_epoch);
         set_environment_ptr( "VECTOR_OPTS", buff);
         }
      }

   if( c == -2)         /* yes,  we're making an ephemeris */
      {
      const double *orbits_to_use = orbit;
      extern const char *ephemeris_filename;
      unsigned n_orbits = 1;
      const bool is_observables =
               ((ephemeris_output_options & 7) == OPTION_OBSERVABLES);

      if( is_observables && (ephemeris_output_options & OPTION_SHOW_SIGMAS))
         orbits_to_use = set_up_alt_orbits( orbit, &n_orbits);
      if( ephemeris_in_a_file_from_mpc_code(
            get_file_name( buff, ephemeris_filename),
            orbits_to_use, obs, n_obs,
            epoch_jd, jd_start, ephemeris_step_size,
            n_ephemeris_steps, mpc_code,
            ephemeris_output_options, n_orbits))
         inquire( "Ephemeris generation failed!  Hit any key:", NULL, 0,
                           COLOR_MESSAGE_TO_USER);
      else
         {
         show_a_file( get_file_name( buff, ephemeris_filename),
                        (is_observables? SHOW_FILE_IS_EPHEM : 0));
         create_resid_file( obs, n_obs, input_filename, residual_format);
         make_pseudo_mpec( get_file_name( buff, "mpec.htm"), obj_name);
         if( ephemeris_output_options
                  & (OPTION_STATE_VECTOR_OUTPUT | OPTION_POSITION_OUTPUT))
            {
            FILE *ofile = fopen_ext( ephemeris_filename, "tfca");

            add_ephemeris_details( ofile, jd_start, jd_end);
            fclose( ofile);
            }
         }
      }
}

static void element_frame_dialog_text( char *buff)
{
   const int curr_frame = atoi( get_environment_ptr( "ELEMENTS_FRAME"));

   strcpy( buff, get_find_orb_text( 2034));
   _set_radio_button( buff, curr_frame);
}

static void select_element_frame( void)
{
   int c;

   do
      {
      char buff[300];

      element_frame_dialog_text( buff);
      help_file_name = "frame_he.txt";
      c = inquire( buff, NULL, 0, COLOR_DEFAULT_INQUIRY);
      if( c >= KEY_F( 1) && c <= KEY_F( 5))
         c += '0' - KEY_F( 1);
      if( c >= '0' && c <= '3')
         {
         char obuff[2];

         obuff[0] = (char)c;
         obuff[1] = '\0';
         set_environment_ptr( "ELEMENTS_FRAME", obuff);
         }
      }
      while( !c);
}

static void object_comment_text( char *buff, const OBJECT_INFO *id)
{
   int i;

   snprintf_err( buff, 12, "%d obs; ", id->n_obs);
   make_date_range_text( buff + strlen( buff), id->jd_start, id->jd_end);
   if( id->jd_updated != id->jd_end)
      {
      char time_buff[80];

      full_ctime( time_buff, id->jd_updated, FULL_CTIME_FORMAT_HH_MM
                  | FULL_CTIME_NO_YEAR | CALENDAR_JULIAN_GREGORIAN);
      snprintf_append( buff, 90, " Updated %s", time_buff);
      }
   for( i = 0; id->mpc_codes[i]; i += 3)
      snprintf_append( buff, 90, " %.3s", id->mpc_codes + i);
}

/* select_object_in_file( ) uses the find_objects_in_file( ) function to
   get a list of all the objects listed in a file of observations.  It
   prints out their IDs and asks the user to hit a key to select one.
   The name of that object is then copied into the 'obj_name' buffer.

   At present,  it's used only in main( ) when you first start findorb,
   so you can select the object for which you want an orbit.         */

extern bool force_bogus_orbit;

int select_object_in_file( OBJECT_INFO *ids, const int n_ids)
{
   static int choice = 0, show_packed = 0;
   int rval = -1;
   int sort_order = OBJECT_INFO_COMPARE_PACKED;
   char search_text[20];

   *search_text = '\0';
   if( ids && n_ids)
      {
      int i, curr_page = 0, err_message = 0;
      static int force_full_width_display = 0;

      clear( );
      while( rval == -1)
         {
         const int xmax = getmaxx( stdscr);
         const int n_lines = getmaxy( stdscr) - 3;
         int column_width = (force_full_width_display ? 40 : 16);
         int c, n_cols = xmax / column_width;
         char *buff = (char *)malloc( xmax > 89 ? xmax + 1 : 91);
         const int x0 = xmax - 25;  /* column where buttons start */

         if( choice < 0)
            choice = 0;
         if( choice > n_ids - 1)
            choice = n_ids - 1;
         while( curr_page > choice)
            curr_page -= n_lines;
         while( curr_page + n_lines * n_cols <= choice)
            curr_page += n_lines;
         if( curr_page < 0)      /* ensure that we wrap around: */
            curr_page = 0;
         if( n_ids < n_lines * n_cols)
            n_cols = (n_ids - 1) / n_lines + 1;
         column_width = xmax / n_cols;
//       if( column_width > 80)
//          column_width = 80;
         for( i = 0; i < n_lines * n_cols; i++)
            {
            char desig[181];
            int color = COLOR_BACKGROUND;
            const bool in_last_column = (i >= n_lines * (n_cols - 1));
            int this_column_width = (in_last_column ?
                    xmax - (n_cols - 1) * column_width : column_width);

            if( this_column_width > (int)sizeof( desig) - 1)
               this_column_width = (int)sizeof( desig) - 1;
            if( i + curr_page < n_ids)
               {
               if( show_packed)
                  snprintf_err( desig, sizeof( desig), "'%s'",
                             ids[i + curr_page].packed_desig);
               else
                  strlcpy_error( desig, ids[i + curr_page].obj_name);
               if( i + curr_page == choice)
                  {
                  snprintf_err( buff, xmax, "Object %d of %d: %s",
                              choice + 1, n_ids, desig);
                  buff[x0] = '\0';
                  put_colored_text( buff, n_lines + 1, 0, -1,
                                                COLOR_SELECTED_OBS);
                  object_comment_text( buff, ids + choice);
                  buff[x0] = '\0';
                  put_colored_text( buff, n_lines + 2, 0, -1,
                                                COLOR_SELECTED_OBS);
                  color = COLOR_HIGHLIT_BUTTON;
                  }
               else
                  if( ids[i + curr_page].solution_exists)
                     color = COLOR_OBS_INFO;
               if( force_full_width_display || this_column_width > 40)
                  {
                  strlcat_error( desig, " ");
                  object_comment_text( desig + strlen( desig), ids + i + curr_page);
                  }
               desig[this_column_width - 1] = ' ';
               desig[this_column_width] = '\0';    /* trim to fit available space */
               }
            else                        /* just want to erase space: */
               *desig = '\0';
            snprintf_err( buff, xmax + 1, "%-*s", this_column_width, desig);
            put_colored_text( buff, i % n_lines,
                   (i / n_lines) * column_width, this_column_width, color);
            }

              /* "Use arrow keys to select object", etc. */
         put_colored_text( get_find_orb_text( 2033),
                                 n_lines, 0, -1, COLOR_SELECTED_OBS);
         if( err_message)
            put_colored_text( "Not a valid choice",
                                 n_lines + 2, 0, -1, COLOR_FINAL_LINE);
         put_colored_text( "Quit", n_lines + 2, x0 + 20, 4, COLOR_HIGHLIT_BUTTON);
         if( curr_page + i < n_ids)
            {
            put_colored_text( "Next", n_lines + 2, x0 + 15, 4, COLOR_HIGHLIT_BUTTON);
            put_colored_text( "End", n_lines + 2, x0 + 6, 3, COLOR_HIGHLIT_BUTTON);
            }
         if( curr_page)
            {
            put_colored_text( "Prev", n_lines + 2, x0 + 10, 4, COLOR_HIGHLIT_BUTTON);
            put_colored_text( "Start", n_lines + 2, x0, 5, COLOR_HIGHLIT_BUTTON);
            }
         put_colored_text( "[?]", 0, xmax - 4, 3,
                              A_REVERSE | COLOR_BACKGROUND);
         if( *search_text)
            put_colored_text( search_text, n_lines + 1, x0,
                      (int)strlen( search_text), COLOR_FINAL_LINE);
         put_colored_text( "Open...", n_lines + 1, x0 + 12, 7, COLOR_HIGHLIT_BUTTON);
         put_colored_text( "HELP", n_lines + 1, x0 + 20, 4, COLOR_HIGHLIT_BUTTON);
         flushinp( );
         do
            {
            c = extended_getch( );
            err_message = 0;
            if( c == KEY_MOUSE)
               {
               int x, y, z;
               mmask_t button;

               get_mouse_data( &x, &y, &z, &button);
               if( button & REPORT_MOUSE_POSITION)
                  c = 0;
               else if( button & BUTTON4_PRESSED)   /* actually 'wheel up' */
                  c = KEY_UP;
               else if( button5_pressed)   /* actually 'wheel down' */
                  c = KEY_DOWN;
               else if( !y && x >= xmax - 4 && x < xmax - 1)
                  c = '?';                 /* clicked on [?] at upper right */
               else if( y < n_lines)
                  choice = curr_page + y + (x / column_width) * n_lines;
               else if( y == n_lines + 1 && x >= x0 - 12 && x < x0 + 19)
                  c = ALT_F;           /* clicked on 'Open...' */
               else if( y == n_lines + 1 || y == n_lines)
                  c = '?';
               else if( y == n_lines + 2 && x >= x0)
                  {
                  const int dx = x - x0;

                  if( dx >= 20)
                     c = 27;          /* quit */
                  else if( dx >= 15)
                     c = KEY_NPAGE;   /* 'next page' */
                  else if( dx >= 10)
                     c = KEY_PPAGE;   /* 'prev page' */
                  else if( dx >= 6)
                     c = KEY_END;     /* end of list */
                  else
                     c = KEY_HOME;    /* start of list */
                  }
               if( button & BUTTON1_DOUBLE_CLICKED)
                  rval = choice;
               }
            } while( !c);
         if( c == 127 || c == KEY_BACKSPACE)
            c = 8;
                     /* if a letter/number is hit,  look for an obj that */
                     /* starts with that letter/number: */
         if( (c >= ' ' && c <= 'z' && c != '?') || c == 8)
            {
            size_t len = strlen( search_text);

            if( c != 8)
               search_text[len++] = (char)c;
            else if( len)
               len--;
            search_text[len] = '\0';
            for( i = 0; i < n_ids; i++)
               if( !memcmp( search_text, ids[i].obj_name, len))
                  {
                  choice = i;
                  c = 0;
                  break;
                  }
            if( c)         /* you can also enter the number of the object */
               {           /* within the on-screen list to select it */
               i = 0;
               while( isdigit( search_text[i]))
                  i++;
               if( !search_text[i])
                  {
                  i = atoi( search_text);
                  if( i > 0 && i <= n_ids)
                     {
                     c = 0;
                     choice = i - 1;
                     }
                  }
               }
            }
         else
            *search_text = '\0';
#ifdef ALT_0
         if( c >= ALT_0 && c <= ALT_9)
            choice = (n_ids - 1) * (c - ALT_0 + 1) / 11;
#endif
         switch( c)
            {
            case '?':
               show_a_file( "obj_help.txt", 0);
               break;
            case 9:
               force_full_width_display ^= 1;
               break;
            case ';':
               sort_order = (sort_order + 1) % 3;
               sort_object_info( ids, n_ids, sort_order);
               break;
            case '!':
               force_bogus_orbit = true;
                     /* FALLTHRU */
            case ' ':
            case 10:
            case 13:
               rval = choice;
               break;
#ifdef KEY_C2
            case KEY_C2:
#endif
            case KEY_DOWN:
               choice++;
               break;
#ifdef KEY_A2
            case KEY_A2:
#endif
            case KEY_UP:
               choice--;
               break;
#ifdef KEY_B1
            case KEY_B1:
#endif
            case KEY_LEFT:
               choice -= n_lines;
               break;
#ifdef KEY_B3
            case KEY_B3:
#endif
            case KEY_RIGHT:
               choice += n_lines;
            break;
            case KEY_C3:         /* "PgDn" = lower right key in keypad */
            case KEY_NPAGE:
            case 'n':
               choice += n_lines * n_cols;
               break;
            case KEY_A3:         /* "PgUp" = upper right key in keypad */
            case KEY_PPAGE:
            case 'p':
               choice -= n_lines * n_cols;
               break;
            case KEY_C1:
            case KEY_END:
            case 'e':
               choice = n_ids - 1;
               break;
            case KEY_A1:
            case KEY_HOME:
            case 's':
               choice = 0;
               break;
            case ',':
               show_packed ^= 1;
               break;
            case ALT_F:
            case 'o': case 'O':
               rval = -3;
               break;
#ifndef _WIN32
            case KEY_F(8):     /* show original screens */
               full_endwin( );
               extended_getch( );
               restart_curses( );
               break;
#endif
            case KEY_MOUSE:      /* already handled above */
               break;
#ifdef KEY_RESIZE
            case KEY_RESIZE:
               resize_term( 0, 0);
               break;
#endif
            case 'q': case 'Q': case 27:
#ifdef KEY_EXIT
            case KEY_EXIT:
#endif
               rval = -2;
               break;
            }
         free( buff);
         }
      if( debug_level > 3)
         debug_printf( "rval = %d; leaving select_object_in_file\n", rval);
      }
   return( rval);
}

void format_dist_in_buff( char *buff, const double dist_in_au); /* ephem0.c */

#define is_power_of_two( X)   (!((X) & ((X) - 1)))

typedef struct
   {
   int key;
   unsigned line, col1, col2;
   } command_area_t;

static command_area_t *command_areas;
static unsigned n_command_areas;

static void add_cmd_area( const unsigned key,
              const unsigned line, const unsigned col1, const unsigned len)
{
   command_area_t *tptr;

   n_command_areas++;
   if( is_power_of_two( n_command_areas))
      command_areas = (command_area_t *)realloc( command_areas,
                         n_command_areas * 2 * sizeof( command_area_t));
   tptr = command_areas + n_command_areas - 1;
   tptr->key = key;
   tptr->line = line;
   tptr->col1 = col1;
   tptr->col2 = col1 + len;
   tptr[1].key = 0;        /* move the null terminator ahead */
}

static int get_character_code( const char *buff)
{
   int rval;

   if( *buff <= ' ')
      rval = 0;
   else if( !memcmp( buff, "Lt", 2))
      rval = KEY_LEFT;
   else if( !memcmp( buff, "Rt", 2))
      rval = KEY_RIGHT;
   else if( !memcmp( buff, "Up", 2))
      rval = KEY_UP;
   else if( !memcmp( buff, "Dn", 2))
      rval = KEY_DOWN;
   else if( !memcmp( buff, "Alt-", 4))
      {
      if( !memcmp( buff + 4, "Up", 2))
         rval = ALT_UP;
      else if( !memcmp( buff + 4, "Dn", 2))
         rval = ALT_DOWN;
      else if( buff[4] >= 'A' && buff[4] <= 'Z')
         rval = ALT_A + buff[4] - 'A';
      else
         rval = ALT_0 + buff[4] - '0';
      }
   else if( *buff == 'F')
      rval = KEY_F( atoi( buff + 1));
   else if( *buff == 'U' && buff[1] == '+')
      rval = atoi( buff + 2);
   else if( !memcmp( buff, "Ctrl-", 5))
      {
      if( !memcmp( buff + 4, "Up", 2))
         rval = CTL_UP;
      else if( !memcmp( buff + 4, "Dn", 2))
         rval = CTL_DN;
      else
         rval = buff[5] - 64;
      }
   else if( !memcmp( buff, "Shift-", 6))
      {
      if( !memcmp( buff + 6, "Up", 2))
         rval = KEY_SR;
      else if( !memcmp( buff + 6, "Dn", 2))
         rval = KEY_SF;
      else if( !memcmp( buff + 6, "Lf", 2))
         rval = KEY_SLEFT;
      else if( !memcmp( buff + 6, "Rt", 2))
         rval = KEY_SRIGHT;
      else if( buff[6] == 'F')
         rval = KEY_F( 12 + atoi( buff + 7));
      else        /* shouldn't happen */
         {
         rval = 0;
         fprintf( stderr, "Buff '%s'\n", buff);
         assert( 0);
         }
      }
   else
      rval = *buff;
   return( rval);
}

static unsigned show_basic_info( const OBSERVE FAR *obs, const int n_obs,
                                          const unsigned max_lines_to_show)
{
   char buff[81];
   double r1, r2;
   unsigned line = 1, column = 24;
   unsigned max_column = (unsigned)getmaxx( stdscr);
   FILE *ifile = fopen_ext( "command.txt", "fcrb");

   get_r1_and_r2( n_obs, obs, &r1, &r2);    /* orb_func.cpp */
   strlcpy_error( buff, "R1:");
   format_dist_in_buff( buff + 3, r1);
   put_colored_text( buff, 0, 0, 15, COLOR_BACKGROUND);

   strlcpy_error( buff, "  R2:");
   format_dist_in_buff( buff + 5, r2);
   put_colored_text( buff, 0, 10, -1, COLOR_BACKGROUND);

   n_command_areas = 0;
   add_cmd_area( 'r', 0, 0, 23);
   add_cmd_area( '?', 0, max_column - 4, 3);
   put_colored_text( "[?]", 0, max_column - 4, 3, COLOR_MENU);
   if( ifile)
      {
      while( fgets_trimmed( buff, sizeof( buff), ifile)
                     && memcmp( buff, "End", 3))
         if( *buff != ' ')
            {
            const unsigned len = (unsigned)strlen( buff + 15);
            unsigned max_column_this_line = max_column, key;

            if( line == 1)
               max_column_this_line -= (max_lines_to_show > 1 ? 12 : 8);
            if( column + len >= max_column_this_line)
               {
               if( line == max_lines_to_show)
                  {
                  fclose( ifile);
                  add_cmd_area( KEY_ADD_MENU_LINE, 0, max_column - 8, 3);
                  put_colored_text( "[+]", 0, max_column - 8, 3, COLOR_MENU);
                  return( line);
                  }
               if( line == 1)        /* we could subtract lines */
                  {
                  add_cmd_area( KEY_REMOVE_MENU_LINE, 0, max_column - 12, 3);
                  put_colored_text( "[-]", 0, max_column - 12, 3, COLOR_MENU);
                  }
               column = 0;
               line++;
               put_colored_text( "", line - 1, column, -1, COLOR_BACKGROUND);
               }
            key = get_character_code( buff);
            add_cmd_area( key, line - 1, column, len);
            put_colored_text( buff + 15, line - 1, column, len, COLOR_MENU);
            column += len + 1;
            }
      fclose( ifile);
      }
   if( line != 1)
      max_column -= 4;
   return( line);
}

static int select_central_object( char *buff, const size_t buffsize, const bool dialog_text_only)
{
   extern int forced_central_body;
   int i, c;
   const char *hotkeys = "AB0123456789L", *tptr;

   *buff = '\0';
   for( i = -2; i <= 10; i++)
      {
      bool is_selected = (i == forced_central_body);

      snprintf_append( buff, buffsize, "%c (%c) ",
               hotkeys[i + 2], is_selected ? '*' : ' ');
      switch( i)
         {
         case -2:
            strlcat( buff, "Automatic\n", buffsize);
            break;
         case -1:
            strlcat( buff, "Barycentric\n", buffsize);
            break;
         case 0:
            strlcat( buff, "Heliocentric\n", buffsize);
            break;
         default:
            {
            strlcat( buff, get_find_orb_text( 99108 + i - 1), buffsize);
            strlcat( buff, "\n", buffsize);
            }
         }
      }
   if( dialog_text_only)
      return( 0);
   strlcat( buff, get_find_orb_text( 96002), buffsize);    /* 'Cancel' */
   c = inquire( buff, NULL, 0, COLOR_DEFAULT_INQUIRY);
   if( c && (tptr = strchr( hotkeys, toupper( c))) != NULL)
      c = KEY_F( 1) + (int)( tptr - hotkeys);
   if( c >= KEY_F( 1) && c <= KEY_F( 13))
      forced_central_body = c - KEY_F( 3);
   return( forced_central_body);
}

const char *find_nth_utf8_char( const char *itext, size_t n);

void show_perturbers( const unsigned line)
{
   int i;

   for( i = 0; i < 11; i++)
      {
      int color = COLOR_BACKGROUND;
      const int shift_amt = (i == 10 ? 20 : i + 1);
      char buff[20];
      char perturber_letter = (i == 10 ? 'a' : '0' + (i + 1) % 10);

      strlcpy_error( buff, "(o)");
      if( (perturbers >> shift_amt) & 1)
         color = COLOR_HIGHLIT_BUTTON;
      else
         buff[1] = perturber_letter;
      put_colored_text( buff, line, i * 7, 3, color);
      add_cmd_area( perturber_letter, line, i * 7, 6);
      strlcpy_error( buff, get_find_orb_text( 99108 + i));
      strcpy( (char *)find_nth_utf8_char( buff, 3), " ");
      put_colored_text( buff, line, i * 7 + 3, (int)strlen( buff),
                                       COLOR_BACKGROUND);
      }
               /* Clear to end of line: */
   put_colored_text( " ", line, i * 7, -1, COLOR_BACKGROUND);
}

int max_mpc_color_codes = 5;
static MPC_STATION *mpc_color_codes = NULL;

/* Show the text for a given observation... which will be in
   'default_color' unless it's excluded.  If it is,  the residuals
   are shown in COLOR_EXCLUDED_OBS. */

#if defined( __APPLE__) || defined( __WATCOMC__)
   static int make_unicode_substitutions = 0;
#else
   static int make_unicode_substitutions = 1;
#endif

static const char *legend =
"   YYYY MM DD.DDDDD   RA (J2000)   dec      sigmas   mag     ref Obs     Xres  Yres   delta  R";

static void show_residual_text( char *buff, const int line_no,
           const int column, const int default_color,
           const int is_included)
{
   const int resid_column = 73;  /* see above legend */
   char *tptr = buff + resid_column - 2;
   const int residual_field_size = 13;
   int resid_color = default_color;
   char tbuff[40];

#ifdef DOESNT_WORK_QUITE_YET
   text_search_and_replace( buff, "\xb5", "\xc2\xb5");
#endif
   put_colored_text( buff, line_no, column, (int)strlen( buff), default_color);

   memcpy( tbuff, tptr, residual_field_size);
   if( !is_included)
      {
      resid_color = COLOR_EXCLUDED_OBS;
      tbuff[0] = '(';                       /* put ()s around excluded obs */
      tbuff[residual_field_size - 1] = ')';
      }
   tbuff[residual_field_size] = '\0';
            /* cvt 'u' = 'micro' = 'mu' to the Unicode U+00B5,  in UTF-8: */
   if( make_unicode_substitutions)
      text_search_and_replace( tbuff, "u", "\xc2\xb5");
   put_colored_text( tbuff, line_no, column + resid_column - 2,
            (int)strlen( tbuff), resid_color);
}

static void show_mpc_code_in_color( const char *mpc_code,
               const int y, const int x)
{
   put_colored_text( mpc_code, y, x, 3, 512 + COLOR_MPC_CODES
                      + find_mpc_color( mpc_color_codes, mpc_code));
}

#define SORT_BY_SCORE 0
#define SORT_BY_NAME  1

static void sort_mpc_codes( const int n_to_sort, const int sort_flag)
{
   int i, do_swap;

   for( i = 0; i < n_to_sort - 1; i++)
      {
      if( sort_flag == SORT_BY_SCORE)
         do_swap = (mpc_color_codes[i].score < mpc_color_codes[i + 1].score);
      else
         do_swap = (strcmp( mpc_color_codes[i].code,
                            mpc_color_codes[i + 1].code) > 0);
      if( do_swap)
         {
         MPC_STATION temp = mpc_color_codes[i];

         mpc_color_codes[i] = mpc_color_codes[i + 1];
         mpc_color_codes[i + 1] = temp;
         if( i)
            i -= 2;
         }
      }
}

static void add_to_mpc_color( const char *code, const int score_increment)
{
   int i;

   for( i = 0; mpc_color_codes[i].code[0]; i++)
      if( !strcmp( mpc_color_codes[i].code, code))
         mpc_color_codes[i].score += score_increment;
}

static void show_right_hand_scroll_bar( const int line_start,
      const int lines_to_show, const int first_line,
      const int n_lines)
{
   if( lines_to_show < n_lines)
      {
      int i;
      const int scroll0 =        first_line            * (lines_to_show - 2) / n_lines + 1;
      const int scroll1 = (first_line + lines_to_show) * (lines_to_show - 2) / n_lines + 1;

      for( i = 0; i < lines_to_show; i++)
         {
         int color = COLOR_SCROLL_BAR;
         const char *text = "|";
         const char *end_text = "-";

         if( !i)
            {
            color = COLOR_OBS_INFO;
            text = (first_line ? "^" : end_text);
            }
         else if( i == lines_to_show - 1)
            {
            color = COLOR_OBS_INFO;
            text = (first_line + lines_to_show == n_lines ? end_text : "v");
            }
         else
            if( i >= scroll0 && i <= scroll1)
               color = COLOR_HIGHLIT_BUTTON;
         if( text == end_text)
            color = COLOR_MENU;
         put_colored_text( text, i + line_start, getmaxx( stdscr) - 1, -1,
                              color);
         }
      }
}

static int show_station_info( const OBSERVE FAR *obs, const int n_obs,
              const int top_line_residual_area,
              const int curr_obs,
              const int list_codes)
{
   int line_no, i;
   int n_obs_shown = getmaxy( stdscr) - top_line_residual_area;
   const int n_mpc_codes = find_mpc_color( mpc_color_codes, NULL);
   int n_stations_shown = (list_codes == SHOW_MPC_CODES_MANY ?
                              n_obs_shown - 3 : n_obs_shown / 3);

   if( n_obs_shown < 0)
      return( 0);
   if( n_stations_shown > n_mpc_codes)
      n_stations_shown = n_mpc_codes;
   if( list_codes == SHOW_MPC_CODES_ONLY_ONE || n_stations_shown < 1)
      n_stations_shown = 1;
   if( !obs)                     /* just getting a count of # stns shown */
      return( n_stations_shown);

                  /* set 'scores' for each code to equal # of times */
                  /* that code was used: */
   for( i = 0; i < n_obs; i++)
      add_to_mpc_color( obs[i].mpc_code, 1);

   sort_mpc_codes( n_mpc_codes, SORT_BY_SCORE);
                /* ...then sort the ones we're going to display by name: */
   sort_mpc_codes( n_stations_shown, SORT_BY_NAME);
   line_no = getmaxy( stdscr) - n_stations_shown;
   for( i = 0; i < n_stations_shown; i++, line_no++)
      {
      char buff[200];
      const int is_curr_code = !strcmp( mpc_color_codes[i].code,
                                        obs[curr_obs].mpc_code);

//    put_colored_text( "", line_no, 0, -1, COLOR_BACKGROUND);
      snprintf_err( buff, sizeof( buff), "(%s)", mpc_color_codes[i].code);
      put_colored_text( "(   ) ", line_no, 0, 6, COLOR_BACKGROUND);
      show_mpc_code_in_color( mpc_color_codes[i].code, line_no, 1);
      put_observer_data_in_text( mpc_color_codes[i].code, buff);
      put_colored_text( buff, line_no, 6, -1,
             (is_curr_code ? COLOR_FINAL_LINE : COLOR_BACKGROUND));
      }
   return( n_stations_shown);
}


static void show_one_observation( OBSERVE obs, const int line,
                        const int residual_format, bool is_underlined)
{
   char buff[200];
   const int mpc_column = (int)( strstr( legend, "Obs") - legend);
   int color = COLOR_BACKGROUND;       /* show in 80-column MPC */
   char resid_data[70];      /* format, w/added data if it fits */
   const int dropped_start = 12;     /* ...but omit designation */

   put_colored_text( "", line, 0, -1, COLOR_BACKGROUND);
   add_to_mpc_color( obs.mpc_code, 1000);
   format_observation( &obs, buff,
                        (residual_format & ~(3 | RESIDUAL_FORMAT_TIME))
                        | RESIDUAL_FORMAT_FOUR_DIGIT_YEARS);
   strlcpy_error( resid_data, buff + 49);
   *resid_data = ' ';
   switch( GET_RESID_TIME_FORMAT( residual_format))
      {
      case 1:     /* force all to MPC80-style decimal day format */
         if( obs.time_precision == 20)
            obs.time_precision = 5;
         else if( obs.time_precision > 6)
            obs.time_precision = 6;
         break;
      case 2:     /* force all to HHMMSS format */
         if( obs.time_precision >= 4 && obs.time_precision <= 6)
            obs.time_precision += 16;
         else if( obs.time_precision < 20 || obs.time_precision > 23)
            obs.time_precision = 23;
         break;
      case 3:     /* force all to MJD */
         if( obs.time_precision <= 6)   /* any MPC80 decimal day format */
            obs.time_precision += 40;      /* corresponding MJD format */
         else if( obs.time_precision >= 20 && obs.time_precision < 24)
            obs.time_precision += 25;
         else
            obs.time_precision = 49;      /* max precision MJD format */
         break;
      }
                     /* show corresponding 1s or 0.1s HHMMSS fmt */
   recreate_observation_line( buff, &obs, residual_format);
   assert( 80 == strlen( buff));
   memmove( buff, buff + dropped_start, strlen( buff + dropped_start) + 1);
   strcat( buff, resid_data);
   if( obs.flags & OBS_IS_SELECTED)
      color = COLOR_SELECTED_OBS;
   if( is_underlined)
      color |= A_UNDERLINE;
   show_residual_text( buff, line, 0, color, obs.is_included);
   show_mpc_code_in_color( buff + mpc_column, line, mpc_column);

   if( residual_format & RESIDUAL_FORMAT_SHOW_DESIGS)
      {
      char *tptr = obs.packed_id;

      if( *tptr && *tptr == ' ')
         while( tptr[1] == ' ')
            tptr++;
      put_colored_text( tptr, line, 44, 9, color + 1);
      }
}

static void show_observations( const OBSERVE *obs, const int first_obs_idx,
                int line_no, const int residual_format, const int n_obs_shown,
                const int n_obs)
{
   int i;

   for( i = first_obs_idx; i < first_obs_idx + n_obs_shown; i++)
      {
      add_cmd_area( KEY_OBSCODE_CLICKED, line_no, 65, 3);
      show_one_observation( obs[i], line_no++, residual_format,
            i == n_obs - 1);
      }
}

static void show_final_line( const int n_obs,
                                      const int curr_obs, const int color)
{
   char buff[90];
   int len;

   snprintf( buff, sizeof( buff), " %d/%d", curr_obs + 1, n_obs);
   len = (int)strlen( buff);
   put_colored_text( buff, getmaxy( stdscr) - 1,
                              getmaxx( stdscr) - len, len, color);
}

static void show_residual_legend( const int line_no, const int residual_format)
{
   char buff[290];

   strlcpy_error( buff, legend);
   if( residual_format & (RESIDUAL_FORMAT_TIME_RESIDS | RESIDUAL_FORMAT_NORMALIZED))
      {         /* residuals in time & cross-time, not RA and dec */
      char *text_loc;

      while( (text_loc = strstr( buff, "Xres")) != NULL)
         if( residual_format & RESIDUAL_FORMAT_TIME_RESIDS)
            *text_loc = 'T';
         else
            memcpy( text_loc, "Res1", 4);
      while( (text_loc = strstr( buff, "Yres")) != NULL)
         if( residual_format & RESIDUAL_FORMAT_TIME_RESIDS)
            *text_loc = 'C';
         else
            memcpy( text_loc, "Res2", 4);
      }
   if( residual_format & RESIDUAL_FORMAT_MAG_RESIDS)
      text_search_and_replace( buff, "delta  R", "delta mresid");

   switch( GET_RESID_TIME_FORMAT( residual_format))
      {
      case 2:     /* HMS form */
         text_search_and_replace( buff, "YYYY MM DD.DDDDD ",
                                     "CYYMMDD:HHMMSSsss");
         break;
      case 3:     /* MJD form */
         text_search_and_replace( buff, "YYYY MM DD.DDDDD ",
                                     "MJDDDDD.ddddddddd");
         break;
      }

   if( residual_format & RESIDUAL_FORMAT_EXTRA)
      strcat( buff, (residual_format & RESIDUAL_FORMAT_TIME_RESIDS)
                  ? "      Xres  Yres  total Mres"
                  : "      Tres  Cres  total Mres");

   if( residual_format & RESIDUAL_FORMAT_SHOW_DESIGS)
      text_search_and_replace( buff, "sigmas", "desigs");
   else
      {
      extern int sigmas_in_columns_57_to_65;

      if( !sigmas_in_columns_57_to_65)
         text_search_and_replace( buff, "sigmas", "------");
      }
   put_colored_text( buff, line_no, 0, -1, COLOR_RESIDUAL_LEGEND);
   add_cmd_area( 't', line_no, 73, 10);  /* click on 'Xres  Yres' resets resid format */
   add_cmd_area( '=', line_no, 85, 7);   /* click on ' delta ' toggles mag resids */
   add_cmd_area( ALT_K, line_no, 44, 8);  /* 'sigmas' toggles sigma display */
   add_cmd_area( '&', line_no, 21, 17);  /* 'ra dec' toggles forced punch card fmt */
   add_cmd_area( ALT_X, line_no, 65, 3);  /* Show list of MPC codes */
   add_cmd_area( 'b', line_no, 0, 20);    /* Reset obs time format _or_ scan obs for time */
}

static int find_rgb( const int irgb);
static void set_color_table( void);

static void drop_starting_columns( char *buff, const int start_column)
{
   char *tptr = (char *)find_nth_utf8_char( buff, (size_t)start_column);

   memmove( buff, tptr, strlen( tptr) + 1);
}

/* Ephemeris files may contain a six-character hexadecimal color,
computed using sky background brightness data,  prefaced by
'$'.  This should be removed before the line is shown to a user.
The RGB value is returned,  and may (or may not) be made use of. */

static int remove_rgb_code( char *buff, int *offset)
{
   unsigned rval = (unsigned)-1;
   size_t i;
   char *loc = buff;

   while( (loc = strchr( loc, '$')) != NULL)
      {
      i = 1;
      while( i < 7 && isxdigit( loc[i]))
         i++;
      if( i == 7)    /* yes,  it's a $ followed by six hex digits */
         {
         sscanf( loc + 1, "%06x", &rval);
         if( offset)
            *offset = (int)( loc - buff);
         memmove( loc, loc + 7, strlen( loc + 6));
         return( (int)rval);
         }
      loc++;      /* loop to look for further '$'s   */
      }
   return( -1);
}

static void show_a_file( const char *filename, const int flags)
{
   FILE *ifile = fopen_ext( filename, "tclrb");
   char buff[560], err_text[100];
   int line_no = 0, keep_going = 1;
   int n_lines = 0, msg_num = 0;
   bool search_text_found = true;
   int n_lines_alloced = 0, search_text_length = 0;
   int *index = NULL, find_text = 0, *backup_screen;
   char search_text[100];
   int calendar_line = -99, calendar_col = -1;
   int calendar_cell_width, calendar_cell_height;
   int start_column = 0;
   size_t max_column_shown = 0;

   if( !ifile)
      ifile = fopen_ext( filename, "clrb");
   if( !ifile)
      {
      snprintf_err( buff, sizeof( buff), "Couldn't open '%s'", filename);
      generic_message_box( buff, "o");
      return;
      }
   *search_text = '\0';
   while( fgets( buff, sizeof( buff), ifile))
      {
      if( n_lines == n_lines_alloced)
         {
         n_lines_alloced = 2 * n_lines_alloced + 50;
         index = (int *)realloc( index, (n_lines_alloced + 1) * sizeof( int));
         }
      n_lines++;
      index[n_lines] = ftell( ifile);
      }
   if( !n_lines)
      {
      fclose( ifile);
      return;
      }
   index[0] = 0;
   *err_text = '\0';
   backup_screen = store_curr_screen( );
   if( flags & SHOW_FILE_IS_CALENDAR)
      {
      const long jd = (long)( current_jd( ) + .5);
      long year, calendar_year0;
      int day, month, day0, n_read;

      fseek( ifile, 0L, SEEK_SET);
      if( !fgets( buff, sizeof( buff), ifile))
         exit( -3);
      n_read = sscanf( buff, "%ld %d %d", &calendar_year0,
                  &calendar_cell_width, &calendar_cell_height);
      assert( 3 == n_read);
      day_to_dmy( jd, &day, &month, &year, CALENDAR_GREGORIAN);
      day0 = (int)( jd - (long)day + 2) % 7;
      day0 += day - 1;
      calendar_line = ((year - calendar_year0) * 12 + month - 1);
      calendar_line *= calendar_cell_height * 5 + 5;
      calendar_line += ((day0 / 7) % 5) * calendar_cell_height + 6;
      calendar_col = (day0 % 7) * calendar_cell_width;
      line_no = calendar_line;
      }
   while( keep_going)
      {
      int i, c;
      const char *msgs[] = { "1Cursor keys to move",
                             "2Already at end of file!",
                             "2Already at top of file!" };
      const int is_ephem = (flags & SHOW_FILE_IS_EPHEM);
      const int top_possible_line = (is_ephem ? 3 : 0);
      const int n_lines_to_show = getmaxy( stdscr) - 1;
      int top_line;
      int color_start = 18;      /* see 'command.txt' */
      short color_pair_idx = (short)color_start;
      const bool color_visibility =
                  n_lines_to_show + color_start < COLORS;

      clear( );
      if( line_no < top_possible_line)
         line_no = top_possible_line;
      if( line_no > n_lines - 1)
         line_no = n_lines - 1;
      top_line = line_no - n_lines_to_show / 2;
      if( top_line >= n_lines - n_lines_to_show)
         top_line = n_lines - n_lines_to_show;
      if( top_line < 0)
         top_line = 0;
      if( start_column > (int)max_column_shown - COLS + 2)
         start_column = (int)max_column_shown - COLS + 2;
      if( start_column < 0)
         start_column = 0;
      max_column_shown = 0;
      if( is_ephem)
         {
         fseek( ifile, 0L, SEEK_SET);
         for( i = 0; i < 3; i++)
            {
            fgets_trimmed( buff, sizeof( buff), ifile);
            if( max_column_shown < strlen( buff))
               max_column_shown = strlen( buff);
            drop_starting_columns( buff, start_column);
            put_colored_text( buff, i, 0, -1, COLOR_BACKGROUND);
            }
         }
      fseek( ifile, index[top_line], SEEK_SET);
      for( i = 0; i < n_lines_to_show
                        && fgets_trimmed( buff, sizeof( buff), ifile); i++)
         {
         const int curr_line = top_line + i;
         int color_col[20], rgb[20], n_colored = 0, j;

         while( (rgb[n_colored] = remove_rgb_code( buff, color_col + n_colored)) >= 0)
            n_colored++;
         if( i < 3)     /* don't try to 'colorize' top lines */
            n_colored = 0;
         if( max_column_shown < strlen( buff))
            max_column_shown = strlen( buff);
         drop_starting_columns( buff, start_column);
         if( i >= 3 || !is_ephem)
            put_colored_text( buff, i, 0, -1,
               (line_no == curr_line ? COLOR_ORBITAL_ELEMENTS : COLOR_BACKGROUND));
         assert( n_colored < 2);
         if( color_visibility)
            for( j = 0; j < n_colored; j++)
               if( color_col[j] >= start_column)
                  {
                  const int blue =  ((rgb[j] >> 3) & 0x1f);
                  const int green = ((rgb[j] >> 11) & 0x1f);
                  const int red =   ((rgb[j] >> 19) & 0x1f);
                  const int text_color = find_rgb(
                                    blue + red + green > 48 ? 0 : 0xffffff);

#ifdef __PDCURSESMOD__
                  init_extended_pair( color_pair_idx, text_color, find_rgb( rgb[j]));
#else
                  init_pair( color_pair_idx, text_color, find_rgb( rgb[j]));
#endif
                  mvchgat( i, color_col[j] - start_column, 2, A_NORMAL, color_pair_idx, NULL);
                  color_pair_idx++;
                  }
         if( calendar_col >= 0 && curr_line > calendar_line - 2
                   && curr_line < calendar_line + calendar_cell_height - 2)
            mvchgat( i, calendar_col + 1, calendar_cell_width - 1, 0, COLOR_OBS_INFO, NULL);
         }
               /* show "scroll bar" to right of text: */
      show_right_hand_scroll_bar( 0, n_lines_to_show, top_line, n_lines);
      snprintf_err( buff, sizeof( buff), "   Line %d of %d", line_no, n_lines);
      put_colored_text( buff, i, getmaxx( stdscr) - (int)strlen( buff) - 1,
                                 (int)strlen( buff), COLOR_FINAL_LINE);
      put_colored_text( "Back", i, 25, 4, COLOR_FINAL_LINE);
      if( line_no < n_lines - 1)
         {
         put_colored_text( "pgDown", i, 30, 6, COLOR_FINAL_LINE);
         put_colored_text( "End", i, 37, 3, COLOR_FINAL_LINE);
         }
      if( line_no)
         {
         put_colored_text( "pgUp", i, 41, 4, COLOR_FINAL_LINE);
         put_colored_text( "Top", i, 46, 3, COLOR_FINAL_LINE);
         }
      put_colored_text( "Save", i, 50, 4, COLOR_FINAL_LINE);
      put_colored_text( "Copy", i, 55, 4, COLOR_FINAL_LINE);

      strlcpy_error( buff, msgs[msg_num] + 1);
      if( *err_text)
         strlcpy_error( buff, err_text);
      else if( search_text[0])
         {
         strlcat_error( buff, "  Find: ");
         strlcat_error( buff, search_text);
         if( !search_text_found)
            strlcat_error( buff, "   Search text not found");
         }
      search_text_found = true;
      put_colored_text( buff, i, 0, (int)strlen( buff), msgs[msg_num][0] - '0');
      *err_text = '\0';
      msg_num = 0;
      flushinp( );
      do
         {
         c = extended_getch( );
         if( c == KEY_MOUSE)
            {
            int x, y, z;
            mmask_t button;

            get_mouse_data( &x, &y, &z, &button);
            if( button & BUTTON4_PRESSED)   /* actually 'wheel up' */
               c = KEY_UP;
            else if( button5_pressed)       /* actually 'wheel down' */
               c = KEY_DOWN;
            else if( button & BUTTON1_CLICKED)
               {
               if( y == i)
                  {
                  if( x >= 25 && x <= 28)       /* "Quit" */
                     c = 27;
                  else if( x >= 30 && x <= 35)  /* "pgDown" */
                     c = KEY_NPAGE;
                  else if( x >= 37 && x <= 39)  /* "End" */
                     c = KEY_END;
                  else if( x >= 41 && x <= 44)  /* "pgUp" */
                     c = KEY_PPAGE;
                  else if( x >= 46 && x <= 48)  /* "Top" */
                     c = KEY_HOME;
                  else if( x >= 50 && x <= 53)  /* "Save" */
                     c = ALT_S;
                  else if( x >= 55 && x <= 58)  /* "Copy" */
                     c = ALT_C;
                  }
               else if( x == getmaxx( stdscr) - 1)  /* clicked scroll bar */
                  line_no = y * n_lines / n_lines_to_show;
               else
                  line_no = y + top_line;
               }
            else
               c = 0;         /* ignorable mouse event */
            }
         } while( !c);
#ifdef ALT_0
      if( c >= ALT_0 && c <= ALT_9)
         line_no = (n_lines - 1) * (c - ALT_0 + 1) / 11;
#endif
      if( !find_text && c >= 'a' && c <= 'z')
         c += ALT_A - 'a';
      switch( c)
         {
         case KEY_C1:
         case KEY_END:
            line_no = n_lines - 1;
            break;
#ifdef KEY_B3
         case KEY_B3:
#endif
         case KEY_RIGHT:
         case 9:     /* tab */
            start_column += 8;
            break;
#ifdef KEY_B1
         case KEY_B1:
#endif
         case KEY_LEFT:
         case KEY_BTAB:
            if( start_column)
               start_column = (start_column - 1) & ~0x7;
            break;
         case KEY_A1:
         case KEY_HOME:
            line_no = top_possible_line;
            break;
         case KEY_UP:
#ifdef KEY_A2
         case KEY_A2:
#endif
            if( line_no > top_possible_line)
               line_no--;
            else
               msg_num = 2;
            break;
         case KEY_DOWN:
#ifdef KEY_C2
         case KEY_C2:
#endif
            if( line_no >= n_lines - 1)
               msg_num = 1;
            else
               line_no++;
            break;
         case KEY_C3:         /* "PgDn" = lower right key in keypad */
         case KEY_NPAGE:
            if( line_no >= n_lines - 1)
               msg_num = 1;
            else
               line_no += n_lines_to_show - 1;
            break;
         case KEY_A3:         /* "PgUp" = upper right key in keypad */
         case KEY_PPAGE:
            if( line_no > top_possible_line)
               line_no -= n_lines_to_show - 1;
            else
               msg_num = 2;
            break;
#ifdef KEY_RESIZE
         case KEY_RESIZE:
            resize_term( 0, 0);
            break;
#endif
         case 27:
#ifdef KEY_EXIT
         case KEY_EXIT:
#endif
            keep_going = 0;
            break;
         case ALT_C:
         case ALT_S:
            {
            if( c == ALT_C)
               make_config_dir_name( buff, "tfile");
            else
               user_select_file( buff, "Save file", 1);
            if( *buff)
               {
               FILE *ofile = fopen( buff, "wb");

               if( ofile)
                  {
                  char tbuff[500];

                  fseek( ifile, 0L, SEEK_SET);
                  while( fgets( tbuff, sizeof( tbuff), ifile))
                     {
                     remove_rgb_code( tbuff, NULL);
                     fputs( tbuff, ofile);
                     }
                  fclose( ofile);
                  }
               }
            if( c == ALT_C && copy_file_to_clipboard( buff))
               inquire( get_find_orb_text( 2039), NULL, 0,
                                    COLOR_MESSAGE_TO_USER);
            }
            break;
         case ALT_Q:
            keep_going = 0;
            break;
         case 127:           /* backspace takes on different values */
         case KEY_BACKSPACE: /* on PDCurses,  ncurses,  etc.        */
         case 8:
            if( search_text_length)
               search_text_length--;
            else
               find_text = 0;
            break;
         case CTRL( 'F'):
            line_no++;
            find_text = 2;
            break;
         case '/':
            find_text = 1;
            search_text_length = 0;
            break;
         case 13:
         case 10:
            find_text = 0;
            search_text_length = 0;
            break;
         case KEY_MOUSE:
            break;
         default:
            if( find_text)
               {
               search_text[search_text_length++] = (char)c;
               find_text = 2;
               }
            break;
         }
      search_text[search_text_length] = '\0';
      if( find_text == 2)
         {
         fseek( ifile, index[line_no], SEEK_SET);
         i = line_no;
         while( fgets_trimmed( buff, sizeof( buff), ifile)
                     && !strstr( buff, search_text))
            i++;
         if( i != n_lines)       /* didn't reach end of file */
            line_no = i;
         else
            search_text_found = false;
         find_text = 1;
         }
      }
   restore_screen( backup_screen);
   set_color_table( );
   free( backup_screen);
   fclose( ifile);
   free( index);
}

static int get_epoch_range_of_included_obs( const OBSERVE FAR *obs,
                  const int n_obs, double *start_jd, double *end_jd)
{
   int idx1, idx2, rval;

   rval = get_idx1_and_idx2( n_obs, obs, &idx1, &idx2);
   if( start_jd)
      *start_jd = obs[idx1].jd;
   if( end_jd)
      *end_jd   = obs[idx2].jd;
   return( rval);
}

static void get_mouse_data( int *mouse_x, int *mouse_y,
                            int *mouse_z, mmask_t *button)
{
   MEVENT mouse_event;
   int i;

   getmouse( &mouse_event);
   *mouse_x = mouse_event.x;
   *mouse_y = mouse_event.y;
   *mouse_z = mouse_event.z;
   *button  = mouse_event.bstate;
   if( *button & REPORT_MOUSE_POSITION)
      _mouse_movements_are_reported = true;
         /* Mouse events are sometimes 'lost' in some terminals.  If we get
         a button release and,  based on past history,  thought the button
         should have already been released,  we turn it into a click.  This
         is an imperfect workaround. */
   for( i = 0; i < 3; i++)
      {
      static int old_buttons_pressed = 0;
      static const mmask_t releases[3] = { BUTTON1_RELEASED, BUTTON2_RELEASED, BUTTON3_RELEASED };
      static const mmask_t presses[3] = { BUTTON1_PRESSED, BUTTON2_PRESSED, BUTTON3_PRESSED };
      static const mmask_t clicks[3] = { BUTTON1_CLICKED, BUTTON2_CLICKED, BUTTON3_CLICKED };
      static const mmask_t button_up[3] = {
                  BUTTON1_RELEASED | BUTTON1_CLICKED | BUTTON1_DOUBLE_CLICKED,
                  BUTTON2_RELEASED | BUTTON2_CLICKED | BUTTON2_DOUBLE_CLICKED,
                  BUTTON3_RELEASED | BUTTON3_CLICKED | BUTTON3_DOUBLE_CLICKED };
      const int mask = 1 << i;

      if( (*button & releases[i]) && !(old_buttons_pressed & mask))
         *button ^= releases[i] | clicks[i];    /* turn a release into a click */
      if( *button & button_up[i])
         old_buttons_pressed &= ~mask;
      else if( *button & presses[i])
         old_buttons_pressed |= mask;
      }
}

static void put_colored_text( const char *text, const int line_no,
               const int column, const int n_bytes, const int color)
{
   const attr_t attrib_mask = A_BOLD | A_BLINK | A_LEFTLINE | A_REVERSE
            | A_RIGHTLINE | A_ITALIC | A_UNDERLINE | A_OVERLINE;

   attrset( COLOR_PAIR( color & 255));
   attron( color & attrib_mask);
   if( n_bytes > 0)
      {
      const int len = getmaxx( stdscr) - column;

      mvaddnstr( line_no, column, text, (n_bytes < len) ? n_bytes : len);
      }
   else              /* clear to end of line */
      {
      const int len = (int)mbstowcs( NULL, text, 0);
      int remains = getmaxx( stdscr) - len - column;

      if( len == -1)
         debug_printf( "Bad text, col %d, line %d: '%s'\n",
                     column, line_no, text);
      assert( len != -1);
#ifdef UNNECESSARY_CODE
      if( remains < 0)
         {
         len += remains;
         remains = 0;
         }
#endif
      mvaddstr( line_no, column, text);
      while( remains > 0)
         {
         char buff[30];
         const int n_bytes = (remains > 29 ? 29 : remains);

         memset( buff, ' ', n_bytes);
         buff[n_bytes] = '\0';
         addnstr( buff, n_bytes);
         remains -= n_bytes;
         }
      }
   attroff( color & attrib_mask);
}

OBSERVE *add_observations( FILE *ifile, OBSERVE *obs,
                  const OBJECT_INFO *ids, int *n_obs)
{
   OBSERVE *obs2 = load_observations( ifile, ids->packed_desig, ids->n_obs);
   extern int n_obs_actually_loaded;

   if( debug_level)
      printf( "Got %d new obs\n", n_obs_actually_loaded);
   fclose( ifile);
   obs = (OBSERVE *)realloc( obs,
                 (*n_obs + n_obs_actually_loaded) * sizeof( OBSERVE));
   memcpy( obs + *n_obs, obs2, n_obs_actually_loaded * sizeof( OBSERVE));
   *n_obs += n_obs_actually_loaded;
   free( obs2);
   *n_obs = sort_obs_by_date_and_remove_duplicates( obs, *n_obs);
   return( obs);
}

int find_first_and_last_obs_idx( const OBSERVE *obs, const int n_obs,
         int *last);       /* elem_out.cpp */
double mid_epoch_of_arc( const OBSERVE *obs, const int n_obs);


/* I really should use getopt() or a portable variant.  However,  this has
been sufficiently effective thus far... */

static const char *get_arg( const int argc, const char **argv, const int idx)
{
   const char *rval;

   argv += idx;
   if( argv[0][2] || idx == argc - 1)
      rval = argv[0] + 2;
   else
      {
      if( argv[1][0] == '-')
         rval = "";
      else
         rval = argv[1];
      }
   return( rval);
}

static bool force_eight_color_mode = false;

      /* Either locate a palette color close to the desired RGB,  or
         (on platforms where one can do so) allocate a new palette entry
         with the desired RGB value.      */

static int find_rgb( const int irgb)
{
   int best_match = 0, best_dist = 4000, i;
   static int n_colors = 0;
   short rgb0[3];

   if( irgb < 0)        /* resetting colors */
      {
      n_colors = 0;
      return( 0);
      }
   rgb0[0] = (int)( irgb >> 16);
   rgb0[1] = (int)( (irgb >> 8) & 0xff);
   rgb0[2] = (int)( irgb & 0xff);
#ifdef __PDCURSESMOD__
   if( COLORS == 256 + (1 << 24))
      return( 256 + rgb0[0] + (rgb0[1] << 8) + (rgb0[2] << 16));
#endif
   if( COLORS <= 8 || force_eight_color_mode)
      return( (rgb0[0] > 127 ? 1 : 0) | (rgb0[1] > 127 ? 2 : 0)
                                      | (rgb0[2] > 127 ? 4 : 0));
   if( COLORS >= 256)
      {
      if( rgb0[0] == rgb0[1] && rgb0[0] == rgb0[2]    /* gray-scale */
                                 && rgb0[0] > 4 && rgb0[0] < 247)
         return( 232 + (rgb0[0] - 4) / 10);
      else
         {
         const int r_idx = (rgb0[0] + 10) / 50;
         const int g_idx = (rgb0[1] + 10) / 50;
         const int b_idx = (rgb0[2] + 10) / 50;

         return( 16 + b_idx + 6 * g_idx + 36 * r_idx);
         }
      }
   for( i = 0; i < 3; i++)
      rgb0[i] = rgb0[i] * 200 / 51;
   for( i = 0; i < n_colors && i < 256; i++)
      {
      short rgb[3];
      int j, dist = 0;

      color_content( (short)i, rgb, rgb + 1, rgb + 2);
      for( j = 0; j < 3; j++)
         dist += abs( rgb[j] - rgb0[j]);
      if( best_dist > dist)
         {
         best_dist = dist;
         best_match = i;
         }
      }
            /* We may have found a suitable color already allocated.  Or
               we may have no choice but to accept the best fit we found. */
   if( best_dist < 15 || n_colors == COLORS)
      return( best_match);
   assert( !force_eight_color_mode);
   init_color( n_colors, rgb0[0], rgb0[1], rgb0[2]);
   n_colors++;
   return( n_colors - 1);
}

#ifndef __PDCURSES__
      /* In ncurses,  use the xterm command for setting title and */
      /* hope it works (it usually will)                          */
static void PDC_set_title( const char *title)
{
   printf( OSC "\x32;%s\a", title);
   fflush( stdout);
}
#endif

static void set_color_table( void)
{
   FILE *ifile = fopen_ext( "command.txt", "fcrb");
   char buff[90];
   const char char_to_search_for = (COLORS > 8 ? 'c' : '8');

   find_rgb( -1);
   while( fgets( buff, sizeof( buff), ifile) && memcmp( buff, "End c", 5))
      if( *buff == char_to_search_for)
         {
         int idx;
         unsigned fore_rgb, back_rgb;

         if( 3 == sscanf( buff + 1, "%d %x %x", &idx, &fore_rgb, &back_rgb))
#ifdef __PDCURSESMOD__
            init_extended_pair( idx, find_rgb( fore_rgb), find_rgb( back_rgb));
#else
            init_pair( idx, find_rgb( fore_rgb), find_rgb( back_rgb));
#endif
         }
   fclose( ifile);
}

static int *get_key_remap_table( void)
{
   FILE *ifile = fopen_ext( "command.txt", "fcrb");
   char buff[90];
   int *rval = (int *)calloc( 200, sizeof( int)), from, to;
   size_t i = 0;

   while( fgets( buff, sizeof( buff), ifile) && memcmp( buff, "Start re", 8))
      ;
   while( fgets( buff, sizeof( buff), ifile))
      if( (from = get_character_code( buff)) > 0
                  && (to = get_character_code( buff + 15)) > 0)
         {
         rval[i++] = from;
         rval[i++] = to;
         }
   rval = (int *)realloc( rval, (i + 2) * sizeof( int));
   fclose( ifile);
   return( rval);
}

static SCREEN *screen_ptr;

static inline int initialize_curses( const int argc, const char **argv)
{
#ifdef __PDCURSES__
   ttytype[0] = 20;    /* Window must have at least 20 lines in Win32a */
   ttytype[1] = 55;    /* Window can have a max of 55 lines in Win32a */
   ttytype[2] = 70;    /* Window must have at least 70 columns in Win32a */
   ttytype[3] = (char)200; /* Window can have a max of 200 columns in Win32a */
#else
   PDC_set_title( get_find_orb_text( 18));
#endif

#ifdef XCURSES
   resize_term( 50, 98);
   Xinitscr( argc, (char **)argv);
   screen_ptr = SP;
#else
   INTENTIONALLY_UNUSED_PARAMETER( argc);
   INTENTIONALLY_UNUSED_PARAMETER( argv);

   char xterm_256color_name[20];

   strlcpy_error( xterm_256color_name, "xterm-256color");
   if( force_eight_color_mode || !(screen_ptr = newterm( xterm_256color_name, stdout, stdin)))
      screen_ptr = newterm( NULL, stdout, stdin);
#endif
   if( debug_level > 2)
      debug_printf( "Curses initialised, ");
   cbreak( );
   if( debug_level > 7)
      debug_printf( "cbreak, ");
   noecho( );
   if( debug_level > 7)
      debug_printf( "noecho, ");
   clear( );
   if( debug_level > 7)
      debug_printf( "clear, ");
   curses_running = true;
   if( debug_level > 2)
      debug_printf( "(2), ");
   start_color( );
#ifdef __PDCURSES__
   PDC_set_blink( TRUE);
   PDC_set_title( get_find_orb_text( 18));
                              /* "Find_Orb -- Orbit Determination" */
#endif
#ifdef VT_RECEIVE_ALL_MOUSE
   printf( VT_RECEIVE_ALL_MOUSE);
#endif
   if( debug_level > 2)
      debug_printf( "(3)\n");
   keypad( stdscr, 1);
#ifdef MOUSE_MOVEMENT_EVENTS_ENABLED
   mousemask( default_mouse_events | REPORT_MOUSE_POSITION, NULL);
#else
   mousemask( default_mouse_events, NULL);
#endif
   set_color_table( );
   return( 0);
}

#ifdef _WIN32
static int user_select_file( char *filename, const char *title, const int flags)
{
   const bool is_save_dlg = (flags & 1);
   OPENFILENAME ofns;
   wchar_t old_path[_MAX_DIR];

   _wgetcwd( old_path, _MAX_DIR);
   memset( &ofns, 0, sizeof( ofns));
   ofns.lStructSize = sizeof( ofns );
   ofns.lpstrFile = filename;
   ofns.nMaxFile = 256;
   ofns.lpstrTitle = title;
   ofns.Flags =  OFN_EXPLORER | OFN_HIDEREADONLY;
   if( is_save_dlg)
      GetSaveFileName( &ofns );
   else
      GetOpenFileName( &ofns );
   _wchdir( old_path);
   return 0;
}
#else

/* In non-Windows situations,  file selection is delegated to the 'zenity'
program. If that's unavailable,  we try 'yad' (fork of zenity with
essentially the same options),  then 'kdialog' (used on KDE),  then
'Xdialog'.  If all else fails, we go to the "traditional" curses
'dialog' program.  (I may add other possibilities as I find them.
The Curses 'dialog' is pretty bad.)

Note that a Windows version of Zenity is available,  which may come in
handy at some point (not used yet) :

https://github.com/maravento/winzenity */

#ifndef __WATCOMC__

static int try_a_file_dialog_program( char *filename, const char *command)
{
   FILE *f = popen( command, "r");

   assert( f);
   if( !fgets_trimmed( filename, 256, f))
      *filename = '\0';
   return( (pclose( f) & 0x4000) ? -1 : 0);
}
#endif

static int user_select_file( char *filename, const char *title, const int flags)
{
#ifndef __WATCOMC__
   const bool is_save_dlg = (flags & 1);
   char cmd[256];
   int rval;
   const bool x_is_running = (NULL != getenv( "DISPLAY"));

   if( !x_is_running)
      {
      inquire( "Enter filename: ", filename, 80, COLOR_DEFAULT_INQUIRY);
      fix_home_dir( filename);
      return( 0);
      }
   strlcpy_error( cmd, "zenity --file-selection");
   if( is_save_dlg)
      strlcat_error( cmd, " --save --confirm-overwrite");
   snprintf_append( cmd, sizeof( cmd), " --title \"%s\"", title);
   strlcat_error( cmd, " 2>/dev/null");
   if( !try_a_file_dialog_program( filename, cmd))
      return( 0);

   memcpy( cmd, "yad   ", 6);
   if( !try_a_file_dialog_program( filename, cmd))
      return( 0);

   snprintf( cmd, sizeof( cmd),
            "kdialog --get%sfilename :find_orb --title \"%s\"",
            (is_save_dlg ? "save" : "open"), title);
   strlcat_error( cmd, " 2>/dev/null");
   if( !try_a_file_dialog_program( filename, cmd))
      return( 0);

   snprintf( cmd, sizeof( cmd), "Xdialog --stdout --title \"%s\"", title);
   strlcat_error( cmd, " --fselect ~ 0 0");
   rval = try_a_file_dialog_program( filename, cmd);
   if( !rval)
      return( 0);

         /* dialog and Xdialog take the same options : */
   full_endwin( );
   snprintf_err( strchr( cmd, '~'), 12, "~ %d %d",
                          getmaxy( stdscr) - 15, getmaxx( stdscr) - 3);
   rval = try_a_file_dialog_program( filename, cmd + 1);
   restart_curses( );
   if( !rval)
      return( 0);
#endif

            /* We have none of the 'usual' file dialogues available
            to us.  So we'll just ask outright for the file name : */
   if( !inquire( "Enter file name :", filename, 100, COLOR_DEFAULT_INQUIRY)
                     && *filename)
      {
      if( *filename == '~')
         text_search_and_replace( filename, "~", getenv( "HOME"));
      return( 0);
      }
   else
      return( -1);
}
#endif

extern const char *elements_filename;

#define DISPLAY_OBSERVATION_DETAILS  2
#define DISPLAY_ORBITAL_ELEMENTS     4

static int count_wide_chars_in_utf8_string( const char *iptr, const char *endptr)
{
   int rval = 0;
   char tbuff[200];

   assert( endptr);
   assert( iptr);
   assert( endptr >= iptr);
   assert( endptr < iptr + sizeof( tbuff));
   memcpy( tbuff, iptr, endptr - iptr);
   tbuff[endptr - iptr] = '\0';
   while( iptr < endptr)
      {
      switch( ((unsigned char)*iptr) >> 4)
         {
         case 0xf:          /* four-byte token;  U+10000 to U+1FFFFF */
            iptr += 4;
            break;
         case 0xe:          /* three-byte token; U+0800 to U+FFFF */
            iptr += 3;
            break;
         case 0xc:          /* two-byte token: U+0080 to U+03FF */
         case 0xd:          /* two-byte token: U+0400 to U+07FF */
            iptr += 2;
            break;
         default:          /* "ordinary" ASCII (U+0 to U+7F) */
            iptr++;        /* single-byte token              */
            break;
         }
      rval++;
      }
   assert( endptr == iptr);
   return( rval);
}

static int toggle_selected_observations( OBSERVE *obs, const unsigned n_obs,
                                 unsigned *n_found)
{
   unsigned n_on = 0, n_off = 0, i;
   int rval;

   for( i = 0; i < n_obs; i++)
      if( obs[i].flags & OBS_IS_SELECTED)
         {
         if( obs[i].is_included)
            n_on++;
         else
            n_off++;
         }
   rval = (n_off > n_on);
   for( i = 0; i < n_obs; i++)
      if( obs[i].flags & OBS_IS_SELECTED)
         obs[i].is_included = rval;
   if( n_found)
      *n_found = n_on + n_off;
   return( rval);
}


/* In Windows,  a canonical file path will contain :\.  Under everything
else I know about,  it won't.  */

static bool filename_fits_current_os( const char *filename)
{
#ifdef _WIN32
   return( strstr( filename, ":\\"));
#else
   return( !strstr( filename, ":\\"));
#endif
}

static OBJECT_INFO *load_file( char *ifilename, int *n_ids, char *err_buff,
                    const size_t err_buff_size,
                    const bool drop_single_obs, const bool already_got_obs)
{
   OBJECT_INFO *ids;
   extern char *temp_obs_filename;     /* miscell.cpp */
   size_t n_lines, i;
   const char *prev_fn = "previous.txt";
   char **prev_files = load_file_into_memory( prev_fn, &n_lines, false);
   bool is_temp_file = false;

   if( !prev_files)
      prev_files = load_file_into_memory( "previous.def", &n_lines, true);
   assert( prev_files);
   *err_buff = '\0';
   if( *ifilename == '@')     /* can use '@' followed by end of filename */
      {                      /* if that file is a recently accessed one */
      size_t len, len1 = strlen( ifilename) - 1;
      struct stat file_info;

      for( i = n_lines - 1; i; i--)
         if( (len = strlen( prev_files[i])) >= len1
                      && !strcmp( prev_files[i] + len - len1, ifilename + 1)
                      && !stat( prev_files[i], &file_info))
            {
            strcpy( ifilename, prev_files[i]);
            break;
            }
      }
   if( !*ifilename)
      {
      const size_t buffsize = 8000;
      char *buff = (char *)malloc( buffsize);
      int c, base_key = (int)KEY_F( 5);
      size_t n_prev = 0, prev_idx[60];
      struct stat file_info;
      FILE *ifile;
      const char *hotkeys = "0123456789abdeghijklmoprtuvwxyz;',./=-[]()*&^%$#@!:<>", *tptr;

      help_file_name = "openfile.txt";
      clear( );
      ifile = fopen_ext( help_file_name, "fclrb");
      i = 0;
      while( fgets_trimmed( buff, buffsize, ifile) && *buff != '$')
         put_colored_text( buff, (int)i++, 0, -1, COLOR_BACKGROUND);
      fclose( ifile);

      strlcpy_err( buff, get_find_orb_text( 2031), buffsize);
#ifdef _WIN32
         {
         char *sel_ptr = strstr( buff, "S ");

         assert( sel_ptr);          /* remove 'selection' option; */
         sel_ptr[-1] = '\0';        /* it doesn't work in Windows */
         base_key--;
         }
#endif
      for( i = n_lines - 1; i && hotkeys[n_prev]
                                     && n_prev + 12 < (size_t)getmaxy( stdscr); i--)
         if( prev_files[i][0] != '#' && filename_fits_current_os( prev_files[i]))
            {
            if( !stat( prev_files[i], &file_info))
               {
               snprintf_append( buff, buffsize, "\n%c ", hotkeys[n_prev]);
               strcat( buff, prev_files[i]);
#ifndef _WIN32
               text_search_and_replace( buff, getenv( "HOME"), "~");
#endif
               prev_idx[n_prev++] = i;
               }
            else      /* file doesn't exist anymore;  remove from list */
               prev_files[i] = NULL;
            }
      strlcat_err( buff, (already_got_obs ? "\nQ Cancel" : "\nQ Quit"), buffsize);

      while( 'C' == (c = inquire( buff, NULL, 30, COLOR_DEFAULT_INQUIRY)))
         show_calendar( );
      free( buff);
      if( c >= ' ' && c < 127)
         {
         tptr = strchr( hotkeys, c);
         if( tptr)
            c = base_key + (int)( tptr - hotkeys);
         }
      i = c - base_key;
      if( i < n_prev)
         strcpy( ifilename, prev_files[prev_idx[i]]);
      else if( i == n_prev)
         c = 'q';
      switch( c)
         {
         case 'F': case 'f': case KEY_F( 1):
            user_select_file( ifilename, "Open astrometry file", 0);
            break;
         case 'N': case 'n': case KEY_F( 2):
            {
            char object_name[83];

            help_file_name = "obj_name.txt";
            if( !inquire( "Enter object name :", object_name, 40, COLOR_DEFAULT_INQUIRY)
                     && *object_name)
               {
               FILE *ofile = fopen_ext( temp_obs_filename, "fwb");

               assert( ofile);
               fetch_astrometry_from_mpc( ofile, object_name);
               fclose( ofile);
               strcpy( ifilename, temp_obs_filename);
               }
            }
            break;
         case 'c': case KEY_F( 3):
            strcpy( ifilename, "c");
            break;
#ifndef _WIN32
         case 's': case KEY_F( 4):
            strcpy( ifilename, "s");
            break;
#endif
         case 'q': case 'Q': case 27:
            free( prev_files);
            *err_buff = '\0';    /* signals 'cancel' */
            return( NULL);
            break;
         }
      }
   if( !*ifilename)
      {
      free( prev_files);
      strlcpy_err( err_buff,
               "'findorb' needs the name of an input file of MPC-formatted\n"
               "astrometry as a command-line argument.\n", err_buff_size);
      return( NULL);
      }

   if( !strcmp( ifilename, "c") || !strcmp( ifilename, "c+")
          || !strcmp( ifilename, "s") || !strcmp( ifilename, "s+"))
      {
      const int err_code = clipboard_to_file( temp_obs_filename, ifilename[1] == '+',
                              *ifilename == 's');

      if( err_code)
         inquire( get_find_orb_text(
#ifdef _WIN32
                        2076  /* "Clipboard is empty" */
#else
                        2039  /* "You need xclip to get at the clipboard" */
#endif
                                    ), NULL, 0, COLOR_MESSAGE_TO_USER);
      strcpy( ifilename, temp_obs_filename);
      is_temp_file = true;
      }

   show_splash_screen( );
   ids = find_objects_in_file( ifilename, n_ids, NULL);
   if( *n_ids > 0 && drop_single_obs)
      {
      int i, j;

      for( i = j = 0; i < *n_ids; i++)
         if( ids[i].n_obs - ids[i].solution_exists > 1)
            ids[j++] = ids[i];
      *n_ids = j;
      }
   if( debug_level > 2)
      debug_printf( "%d objects in file\n", *n_ids);
   if( *n_ids <= 0)
      {        /* no objects found,  or file not found */
      const char *err_msg;

      if( *n_ids == -1)
         err_msg = "Couldn't locate the file '%s'\n";
      else
         err_msg = "No objects found in file '%s'\n";
      snprintf_err( err_buff, err_buff_size, err_msg, ifilename);
      if( ids)
         free( ids);
      ids = NULL;
      }
   else if( !is_temp_file)
      {
      FILE *ofile = fopen_ext( prev_fn, "fcw");
#if defined( _WIN32)
      char canonical_path[MAX_PATH];
      DWORD len =
                 GetFullPathNameA( ifilename, MAX_PATH, canonical_path, NULL);

      assert( len && len < MAX_PATH);
#else
   #if defined( __WATCOMC__)
      char *canonical_path = ifilename;
   #else
      char *canonical_path = realpath( ifilename, NULL);

      assert( canonical_path);
   #endif
#endif
      set_solutions_found( ids, *n_ids);
      for( i = 0; i < n_lines; i++)
         if( prev_files[i] && strcmp( prev_files[i], canonical_path))
            fprintf( ofile, "%s\n", prev_files[i]);
      if( strcmp( canonical_path, temp_obs_filename))
         fprintf( ofile, "%s\n", canonical_path);
      fclose( ofile);
#ifndef _WIN32
      free( canonical_path);
#endif
      }
   free( prev_files);
   return( ids);
}

static int non_grav_menu( char *message_to_user)
{
   char buff[400];
   int c;
   size_t i = 0;
   extern int force_model;
   static int models[] = { FORCE_MODEL_NO_NONGRAVS, FORCE_MODEL_SRP,
            FORCE_MODEL_SRP_TWO_PARAM, FORCE_MODEL_SRP_THREE_PARAM,
            FORCE_MODEL_COMET_TWO_PARAM, FORCE_MODEL_COMET_THREE_PARAM,
            FORCE_MODEL_COMET_FOUR_PARAM, FORCE_MODEL_YARKO_A2,
            FORCE_MODEL_DELTA_V };
   const size_t n_models = sizeof( models) / sizeof( models[0]);

   while( i < n_models && force_model != models[i])
      i++;
   assert( i < n_models);
   strlcpy_error( buff, get_find_orb_text( 2060));
   _set_radio_button( buff, (int)i);
   help_file_name = "nongravs.txt";
   c = full_inquire( buff, NULL, 0, COLOR_MENU, -1, -1);
   if( c >= KEY_F(1) && c <= KEY_F(8))
      c -= KEY_F( 1);
   else
      c -= '0';
   if( c >= 0 && c < (int)n_models)
      {
      char *tptr;

      force_model = models[c];
      n_orbit_params = 6 + (force_model & 0xf);
      strlcpy_error( buff, get_find_orb_text( 2060));
      tptr = _set_radio_button( buff, c) + 6;
      for( i = 0; tptr[i] >= ' '; i++)
         message_to_user[i] = tptr[i];
      message_to_user[i] = '\0';
      }
   else
      *message_to_user = '\0';
   return( force_model);
}

static int _obs_format_menu( int residual_format, const bool is_time_fmt)
{
   char tbuff[300];
   const int curr_format = (is_time_fmt ?
                        GET_RESID_TIME_FORMAT( residual_format) :
                        GET_RESID_RA_DEC_FORMAT( residual_format));
   const int n_formats = 4;    /* at least for now,  for both */
   int c;

   strlcpy_error( tbuff, get_find_orb_text( is_time_fmt ? 2070 : 2068));
   assert( curr_format >= 0);
   assert( curr_format < n_formats);
   _set_radio_button( tbuff, curr_format);
// help_file_name = (is_time_fmt ? "time_fmt.txt" : "radecfmt.txt");
   c = full_inquire( tbuff, NULL, 0, COLOR_MENU, -1, -1);
   if( c >= KEY_F(1) && c <= KEY_F(14))
      c -= KEY_F( 1);
   else
      c -= '0';
   if( c >= 0 && c < n_formats)
      {
      residual_format &= (is_time_fmt ? ~RESIDUAL_FORMAT_TIME :
                                        ~RESIDUAL_FORMAT_RA_DEC);
      residual_format |= c << (is_time_fmt ? 11 : 14);
      }
   return( residual_format);
}

static int resid_format_menu( char *message_to_user, int resid_format)
{
   char buff[300];
   int c, option_selected = 0;
   size_t i;

   strcpy( buff, get_find_orb_text( 2062));
   if( resid_format & RESIDUAL_FORMAT_NORMALIZED)
      option_selected = 2;
   else if( resid_format & RESIDUAL_FORMAT_TIME_RESIDS)
      option_selected = 1;
   _set_radio_button( buff, option_selected);
   help_file_name = "residfmt.txt";
   c = full_inquire( buff, NULL, 0, COLOR_MENU, -1, -1);
   if( c >= KEY_F(1) && c <= KEY_F(3))
      c += '0' - KEY_F( 1);
   if( c >= '0' && c <= '2')
      {
      char *line_ptr;

      resid_format &= ~(RESIDUAL_FORMAT_TIME_RESIDS | RESIDUAL_FORMAT_NORMALIZED);
      switch( c)
         {
         case '0':
            break;
         case '1':
            resid_format |= RESIDUAL_FORMAT_TIME_RESIDS;
            break;
         case '2':
            resid_format |= RESIDUAL_FORMAT_NORMALIZED;
            break;
         }
      strlcpy_error( buff, get_find_orb_text( 2062));
      line_ptr = _set_radio_button( buff, c - '0') + 6;
      for( i = 0; line_ptr[i] >= ' '; i++)
         message_to_user[i] = line_ptr[i];
      message_to_user[i] = '\0';
      }
   else
      *message_to_user = '\0';
   return( resid_format);
}

static void setup_elements_dialog( char *buff, const char *constraints)
{
   int pass;
   char tbuff[300], *tptr;
   FILE *ifile = fopen_ext( get_file_name( tbuff, elements_filename), "tfcrb");

   strcpy( buff, get_find_orb_text( 2037));
   text_search_and_replace( buff, "$REF",
                  get_environment_ptr( "REFERENCE"));
   if( !*constraints)
      constraints = "(none)";
   text_search_and_replace( buff, "$CON", constraints);

   for( pass = 0; pass < 2; pass++)
      {
      if( pass)
         select_central_object( tbuff, sizeof( tbuff), true);
      else
         element_frame_dialog_text( tbuff);
      tptr = strchr( tbuff, '*');
      if( *tptr)     /* found the element frame */
         {
         size_t i;

         tptr += 3;
         for( i = 0; tptr[i] >= ' '; i++)
            tbuff[i] = tptr[i];
         tbuff[i] = '\0';
         }
      else
         strlcpy_error( tbuff, "(?!)");
      text_search_and_replace( buff, (pass ? "$CEN" : "$FRA"), tbuff);
      }
   while( fgets( tbuff, sizeof( tbuff), ifile))
      if( !memcmp( tbuff, "Epoch ", 6)
                     && (tptr = strstr( tbuff + 6, " TT")) != NULL)
         {
         *tptr = '\0';
         text_search_and_replace( buff, "$EPO", tbuff + 6);
         }
   fclose( ifile);
}

static void eop_info_text( char *buff, const size_t buffsize)
{
   int eop_range[3];

   load_earth_orientation_params( NULL, eop_range);
   if( eop_range[0])
      {
      char date_buff[3][50];
      size_t i;

      for( i = 0; i < 3; i++)
         full_ctime( date_buff[i], 2400000.5 + (double)eop_range[i],
               FULL_CTIME_DATE_ONLY | FULL_CTIME_YMD);
      snprintf_err( buff, (int)buffsize,
               "EOPs run from %s to %s\n(%s with extrapolation)\n",
                        date_buff[0], date_buff[2], date_buff[1]);
      }
   else
      strlcpy_err( buff, "No EOPs available\n", buffsize);
}

static void debias_info_text( char *buff, const size_t buffsize)
{
   const int debias_version = find_fcct_biases( 1., 1., 0,
                           0., NULL, NULL);

   if( debias_version > 2000)
      {
      const char *ver = "?unknown?\n";

      strlcpy_err( buff, "Astrometric debiasing version ", buffsize);
      if( debias_version == 2018)
         ver = "EFCC18\n";
      else if( debias_version == 2014)
         ver = "FCCT14\n";
      strlcat_err( buff, ver, buffsize);
      }
   else
      strlcpy_err( buff, "No astrometric debiasing applied\n", buffsize);
}

static void show_splash_screen( void)
{
   FILE *ifile = fopen_ext( "splash.txt", "crb");

   if( ifile)
      {
      char buff[200], eop_line_1[200], *eop_line_2;
      char debias_text[100], jpl_ephem_text[250];
      bool show_it = false;
      int lines, cols, i;

      debias_info_text( debias_text, sizeof( debias_text));
      eop_info_text( eop_line_1, sizeof( eop_line_1));
      eop_line_2 = strchr( eop_line_1, '\n');
      if( eop_line_2)
         *eop_line_2++ = '\0';
      if( !format_jpl_ephemeris_info( jpl_ephem_text))
         strlcpy_error( jpl_ephem_text, " No JPL DE ephemeris set up");
      clear( );
      while( !show_it && fgets( buff, sizeof( buff), ifile))
         if( 2 == sscanf( buff, "%d %d", &lines, &cols))
            {
            show_it = (lines < LINES && cols < COLS);
            for( i = 0; i < lines && fgets_trimmed( buff, sizeof( buff), ifile); i++)
               if( show_it)
                  {
                  text_search_and_replace( buff, "$v", find_orb_version_jd( NULL));
                  text_search_and_replace( buff, "$d", debias_text);
                  text_search_and_replace( buff, "$e1", eop_line_1);
                  text_search_and_replace( buff, "$e2", (eop_line_2 ? eop_line_2 : ""));
                  text_search_and_replace( buff, "$j", jpl_ephem_text + 1);
                  put_colored_text( buff, (LINES - lines) / 2 + i, 0, -1, COLOR_BACKGROUND);
                  }
            }
      doupdate( );
      fclose( ifile);
      }
}

static int find_command_area( const unsigned mouse_x, const unsigned mouse_y,
                              size_t *index)
{
   size_t i;
   int rval = KEY_MOUSE;

   if( index)
      *index = (size_t)( -1);
   for( i = 0; command_areas[i].key; i++)
      if( mouse_y == command_areas[i].line &&
                    mouse_x >= command_areas[i].col1 &&
                    mouse_x < command_areas[i].col2)
         {
         rval = command_areas[i].key;
         if( index)
            *index = i;
         }
   return( rval);
}

static void show_hint_text( const unsigned mouse_x, const unsigned mouse_y,
               unsigned line, unsigned col, unsigned len, const char *text)
{
   chtype prev_scr[100];
   unsigned i;
   int pass;

   assert( len < 100);
   mvinchnstr( line, col, prev_scr, len);
   for( pass = 0; pass < 2; pass++)
      {
      for( i = 0; i < len; i++)
         prev_scr[i] ^= A_REVERSE;
      mvaddchnstr( line, col, prev_scr, len);
      if( !pass)
         full_inquire( text, NULL, HINT_TEXT, COLOR_MENU, mouse_y, mouse_x);
      }
}

static void show_hint_from_command_area( const unsigned mouse_x,
               const unsigned mouse_y, const size_t index, const char *hint_text)
{
   show_hint_text( mouse_x, mouse_y, command_areas[index].line,
                     command_areas[index].col1,
                     command_areas[index].col2 - command_areas[index].col1,
                     hint_text);
}

static int show_hint( const unsigned mouse_x, const unsigned mouse_y, size_t *index)
{
   const int cmd = find_command_area( mouse_x, mouse_y, index);

   if( cmd != KEY_MOUSE && cmd > 0)
      {
      char buff[100];
      FILE *ifile = fopen_ext( "command.txt", "fcrb");
      int got_it = 0;

      while( fgets_trimmed( buff, sizeof( buff), ifile)
                  && memcmp( buff, "Start h", 7))
         ;   /* just skip lines until we get to the section we want */
      while( !got_it && fgets_trimmed( buff, sizeof( buff), ifile)
                         && memcmp( buff, "Start re", 8))
         if( *buff > ' ')
            got_it = (get_character_code( buff) == cmd);
      fclose( ifile);
      if( got_it)
         {
         char hotkey[15];

         memcpy( hotkey, buff, 15);
         *strchr( hotkey, ' ') = '\0';
         memmove( buff, buff + 15, strlen( buff + 14));
         if( hotkey[1] != '+')      /* skip the U+nnnn cases */
            snprintf_append( buff, sizeof( buff), "\nHotkey for this is %s", hotkey);
         }
      else
         *buff = '\0';
      if( *buff)      /* yes,  we have a hint to show */
         show_hint_from_command_area( mouse_x, mouse_y, *index, buff);
      }
   return( cmd);
}

   /* On any platform with ASLR,  the address of 'zval' will be
   cyptographically selected at startup time.  (Or at least,  some
   bits of it will be.)  It should be more than random enough for
   our not-very-security-crucial needs.  Gets a compile failure on
   Windows,  though,  so we'll just use the clock there. */

static unsigned get_random_seed( void)
{

#ifdef _WIN32
   int64_t zval = nanoseconds_since_1970( );
#else
   unsigned long zval;

   zval = (unsigned long)&zval;
#endif
   return( (unsigned)( zval ^ (zval >> 31)));
}

int sanity_test_observations( const char *filename);

/* main( ) begins by using the select_object_in_file( ) function (see above)
   to determine which object is going to be analyzed in the current 'run'.
   It loads up the observations for the object to be analyzed,
   and computes an initial guess for the orbit.

   It then goes into a loop,  setting up and displaying a screenful of
   data,  then asking the user to hit a key.  Based on which key is hit,
   an action such as taking an Herget step or a "full step",  or resetting
   the epoch,  or quitting,  is taken. */

int main( int argc, const char **argv)
{
   char obj_name[84], tbuff[500], orbit_constraints[90];
   char ifilename[256];
   unsigned n_command_lines = 1;
   int c = 1, element_precision,  add_off_on = -1;
   bool get_new_object = true, get_new_file = true;
   unsigned top_line_orbital_elements = 0;
   unsigned top_line_residuals;
   bool is_monte_orbit = false;
   unsigned list_codes = SHOW_MPC_CODES_NORMAL;
   int i, quit = 0, n_obs = 0, clock_line = 0;
   OBSERVE FAR *obs = NULL;
   int curr_obs = 0;
   double epoch_shown, curr_epoch, orbit[MAX_N_PARAMS];
   double r1 = 1., r2 = 1.;
   char message_to_user[180];
   int update_element_display = 1, gauss_soln = 0, residual_line = 0;
   int residual_format = RESIDUAL_FORMAT_80_COL, bad_elements = 0;
   int element_format = 0, debug_mouse_messages = 0, prev_getch = 0;
   int auto_repeat_full_improvement = 0, n_ids = 0, planet_orbiting = 0;
   OBJECT_INFO *ids = NULL;
   double noise_in_sigmas = 1.;
   double monte_data[MONTE_DATA_SIZE];
   extern int monte_carlo_object_count;
   extern char default_comet_magnitude_type;
   extern double max_monte_rms;
   extern int use_config_directory;          /* miscell.c */
   double max_residual_for_filtering = 2.5;
   bool show_commented_elements = false;
   bool drop_single_obs = true;
   int sort_obs_by_code = 0, *key_remaps;
   int n_stations_shown = 0, top_obs_shown = 0, n_obs_shown = 0;
   bool single_obs_selected = false;
   extern unsigned random_seed;
   unsigned mouse_x = (unsigned)-1, mouse_y = 0, mouse_z = 0;
   mmask_t button = 0;

   random_seed = get_random_seed( );
   if( !strcmp( argv[0], "find_orb"))
      use_config_directory = true;
   else
      use_config_directory = false;
   ensure_config_directory_exists();
   if( !setlocale( LC_ALL, "C.UTF-8"))
      setlocale( LC_ALL, "en_US.utf8");
   *ifilename = '\0';

   if( reset_astrometry_filename( &argc, argv))
      drop_single_obs = false;

   for( i = 1; i < argc; i++)       /* check to see if we're debugging: */
      if( argv[i][0] == '-')
         {
         const char *arg = get_arg( argc, argv, i);

         assert( arg);
         switch( argv[i][1])
            {
            case '1':
               drop_single_obs = false;
               break;
            case '8':
               force_eight_color_mode = true;
               break;
            case 'a':
               {
               extern int separate_periodic_comet_apparitions;

               separate_periodic_comet_apparitions ^= 1;
               }
               break;
            case 'c':
               {
               extern const char *combine_all_observations;

               combine_all_observations = arg;
               if( strlen( arg) > 12)
                  {           /* warn that packed desig is overlong */
                  snprintf_err( tbuff, sizeof( tbuff), get_find_orb_text( 2066),
                        arg, (int)strlen( arg));
                  inquire( tbuff, NULL, 30, COLOR_DEFAULT_INQUIRY);
                  return( -1);
                  }
               }
               break;
            case 'd':
               debug_level = atoi( arg);
               if( !debug_level)
                  debug_level = 1;
               debug_printf( "findorb: debug_level = %d; %s %s\n",
                           debug_level, __DATE__, __TIME__);
               break;
            case 'D':
               if( load_environment_file( arg))
                  {
                  fprintf( stdout, "Couldn't read environment file '%s'\n", arg);
                  return( -1);
                  }
               break;
            case 'f':            /* obj designation;  fall through, */
               break;            /* handle below */
            case 'i':
               {
               extern int ignore_prev_solns;

               ignore_prev_solns = 1;
               }
               break;
            case 'j':
               {
               extern bool force_final_full_improvement;

               force_final_full_improvement = true;
               }
               break;
            case 'l':
               {
               extern char findorb_language;

               findorb_language = *arg;
               }
               break;
            case 'L':
               if( !setlocale( LC_ALL, arg))
                  debug_printf( "Couldn't set locale '%s'\n", arg);
               break;
            case 'm':
               {
               extern int integration_method;

               integration_method = atoi( arg);
               }
               break;
            case 'n':
               max_mpc_color_codes = atoi( arg);
               break;
            case 'N':
               {
               extern const char *fullname_pattern;

               fullname_pattern = arg;
               }
               break;
            case 'O':          /* write output files to specified dir */
               {
               extern const char *output_directory;

               output_directory = arg;
               }
               break;
            case 'o':            /* obj designation / ephemeris from orbital */
               break;            /* elems:  fall through, handle below */
            case 'p':
               {
               extern int process_count;

               process_count = atoi( arg);
               }
               break;
            case 'P':
               {
               extern const char *desig_pattern;

               desig_pattern = arg;
               }
               break;
            case 'q':
               {
               extern bool take_first_soln;

               take_first_soln = true;
               }
               break;
            case 'r':
               {
               const char *comma = strchr( arg, ',');

               if( comma)
                  {
                  extern double minimum_observation_jd;  /* default is 1     */
                  extern double maximum_observation_jd;  /* default is +1e+9. */
                  const size_t len = comma - arg;

                  assert( len < sizeof( tbuff));
                  memcpy( tbuff, arg, len);
                  tbuff[len] = '\0';
                  if( len)
                     minimum_observation_jd =
                          get_time_from_string( 0., tbuff,
                          FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, NULL);
                  if( comma[1])
                     maximum_observation_jd =
                          get_time_from_string( 0., comma + 1,
                          FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, NULL);
                  else
                     maximum_observation_jd = 1e+9;
                  }
               }
               break;
            case 's':
               sanity_test_observations( ifilename);
               printf( "Sanity check complete\n");
               exit( 0);
//             break;
            case 'S':
               {
               extern int sanity_check_observations;

               sanity_check_observations = 0;
               }
               break;
            case 'u':
               make_unicode_substitutions ^= 1;
               break;
            case 'v':
               {
               extern const char *state_vect_text;

               state_vect_text = arg;
               }
               break;
            case 'x':
               {
               extern bool saving_elements_for_reuse;

               saving_elements_for_reuse = true;
               }
               break;
            case 'y':
               {
               extern int n_extra_full_steps;

               n_extra_full_steps = atoi( arg);
               }
               break;
            case 'z':
               {
               extern const char *alt_config_directory;

               use_config_directory = true;
               alt_config_directory = arg;
               }
               break;
            default:
               printf( "Unknown command-line option '%s'\n", argv[i]);
               return( -1);
            }
         }
   sscanf( get_environment_ptr( "CONSOLE_OPTS"), "%s %*d %d %u",
               mpc_code, &residual_format, &list_codes);
   text_search_and_replace( mpc_code, "_", " ");

   residual_format |= RESIDUAL_FORMAT_80_COL;      /* force 80-column mode */

   for( i = 1; i < argc; i++)
      if( argv[i][0] != '-')
         {
         const char *tptr = strchr( argv[i], '=');

         if( tptr)
            {
            const size_t len = tptr - argv[i];

            memcpy( tbuff, argv[i], len);
            tbuff[len] = '\0';
            set_environment_ptr( tbuff, argv[i] + len + 1);
            }
         else if( !*ifilename)
            strlcpy_error( ifilename, argv[i]);
         }

   get_defaults( &ephemeris_output_options, &element_format,
         &element_precision, &max_residual_for_filtering,
         &noise_in_sigmas);

   strlcpy_err( ephemeris_start, get_environment_ptr( "EPHEM_START"),
                  sizeof( ephemeris_start));
   strlcpy_err( ephemeris_step_size, get_environment_ptr( "EPHEM_STEP_SIZE"),
                  sizeof( ephemeris_step_size));
   sscanf( get_environment_ptr( "EPHEM_STEPS"), "%d %79s",
               &n_ephemeris_steps, ephemeris_step_size);
   if( debug_level)
      debug_printf( "Options read\n");

   i = load_up_sigma_records( "sigma.txt");

   if( debug_level)
      debug_printf( "%d sigma recs read\n", i);
   key_remaps = get_key_remap_table( );

   initialize_curses( argc, argv);

   *message_to_user = '\0';
   while( !quit)
      {
      int line_no = 0, mouse_wheel_motion;
      extern bool saving_elements_for_reuse;

      prev_getch = c;
      while( get_new_file || get_new_object)
         {
         while( get_new_file)
            {
            OBJECT_INFO *new_ids;
            int n_new_ids;

            new_ids = load_file( ifilename, &n_new_ids, tbuff, 200, drop_single_obs,
                                     (ids ? true : false));
            if( !new_ids && !*tbuff && !ids)   /* at startup,  and hit Quit */
               goto Shutdown_program;
            if( !new_ids && *tbuff)
               {
               *ifilename = '\0';
               inquire( tbuff, NULL, 30, COLOR_DEFAULT_INQUIRY);
               }
            if( new_ids)
               {
               if( ids)
                  free( ids);
               ids = new_ids;
               n_ids = n_new_ids;
               get_new_object = true;
               get_new_file = false;
               }
            else if( !*tbuff)    /* hit Cancel,  going back to curr obj */
               {
               get_new_file = false;
               get_new_object = !obs;
               }
            }
         if( debug_level > 3)
            debug_printf( "get_new_object = %d\n", (int)get_new_object);
         if( get_new_object)
            {
            int id_number = 0;

            assert( ids);
            assert( n_ids);
            strlcpy_err( tbuff, get_find_orb_text( 18), sizeof( tbuff));
            strlcat_err( tbuff, " | ", sizeof( tbuff));
            i = (int)strlen( ifilename);
            while( i && ifilename[i - 1] != '/' && ifilename[i - 1] != '\\')
               i--;
            strlcat_err( tbuff, ifilename + i, sizeof( tbuff));
            PDC_set_title( tbuff);
            if( n_ids > 1)
               id_number = select_object_in_file( ids, n_ids);
            if( debug_level > 3 && id_number >= 0)
               debug_printf( "id_number = %d; '%s'\n", id_number,
                                       ids[id_number].obj_name);
            get_new_object = false;
            *orbit_constraints = '\0';
            if( id_number == -3)          /* 'Open...' clicked */
               {
               *ifilename = '\0';
               get_new_object = get_new_file = true;
               }
            else if( id_number < 0)
               goto Shutdown_program;
            else
               {
               FILE *ifile;
               long file_offset;

               strlcpy_error( obj_name, ids[id_number].obj_name);
               snprintf_err( tbuff, sizeof( tbuff), "Loading '%s'...", obj_name);
               put_colored_text( tbuff, getmaxy( stdscr) - 3,
                                    0, -1, COLOR_FINAL_LINE);
               if( debug_level)
                  debug_printf( "%s: ", tbuff);
               refresh( );
               monte_carlo_object_count = 0;

               ifile = fopen( ifilename, "rb");
                   /* Start quite a bit ahead of the actual data,  just in case */
                   /* there's a #Sigma: or something in the observation header */
                   /* to which we should pay attention:                        */
               file_offset = ids[id_number].file_offset - 40000L;
               if( file_offset < 0L)
                  file_offset = 0L;
               fseek( ifile, file_offset, SEEK_SET);
               if( file_offset)        /* read and discard partial line */
                  if( !fgets( tbuff, sizeof( tbuff), ifile))
                     return( -4);
               if( obs)
                  unload_observations( obs, n_obs);

               strlcpy_error( tbuff, get_find_orb_text( 18));
               strlcat_error( tbuff, " | ");
               strlcat_error( tbuff, ids[id_number].obj_name);
               PDC_set_title( tbuff);
               show_splash_screen( );
               obs = load_object( ifile, ids + id_number, &curr_epoch,
                                                     &epoch_shown, orbit);
               assert( obs);
               fclose( ifile);
               n_obs = ids[id_number].n_obs;
               if( !curr_epoch || !epoch_shown || !obs || n_obs < 2)
                  debug_printf( "Curr epoch %f; shown %f; obs %p; %d obs\n",
                                 curr_epoch, epoch_shown, (void *)obs, n_obs);
               if( debug_level)
                  debug_printf( "got obs; ");
               if( mpc_color_codes)
                  free( mpc_color_codes);
               if( max_mpc_color_codes)
                  mpc_color_codes = find_mpc_color_codes( n_obs, obs,
                             max_mpc_color_codes);
               if( debug_level)
                  debug_printf( "got color codes; ");
               get_r1_and_r2( n_obs, obs, &r1, &r2);    /* orb_func.cpp */
               if( debug_level)
                  debug_printf( "R1 = %f; R2 = %f\n", r1, r2);
               for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
                  ;
               curr_obs = top_obs_shown = i;
               update_element_display = 1;
               clear( );
               }
            force_bogus_orbit = false;
            }
         if( !n_obs)
            {
            get_new_file = true;
            *ifilename = '\0';
            }
         }
      if( curr_obs > n_obs - 1)
         curr_obs = n_obs - 1;
      if( curr_obs < 0)
         curr_obs = 0;
      if( single_obs_selected)
         {
         for( i = 0; i < n_obs; i++)
            obs[i].flags &= ~OBS_IS_SELECTED;
         obs[curr_obs].flags |= OBS_IS_SELECTED;
         }
      if( debug_level > 2)
         debug_printf( "update_element_display = %d\n", update_element_display);
      if( residual_format & (RESIDUAL_FORMAT_PRECISE | RESIDUAL_FORMAT_OVERPRECISE))
         element_format |= ELEM_OUT_PRECISE_MEAN_RESIDS;
      else
         element_format &= ~ELEM_OUT_PRECISE_MEAN_RESIDS;
      if( residual_format & RESIDUAL_FORMAT_NORMALIZED)
         element_format |= ELEM_OUT_NORMALIZED_MEAN_RESID;
      else
         element_format &= ~ELEM_OUT_NORMALIZED_MEAN_RESID;
      if( update_element_display)
         bad_elements = write_out_elements_to_file( orbit, curr_epoch, epoch_shown,
             obs, n_obs, orbit_constraints, element_precision,
             is_monte_orbit, element_format);
      saving_elements_for_reuse = false;
      is_monte_orbit = false;
      if( debug_level > 2)
         debug_printf( "elements written\n");
      update_element_display = 0;
      n_command_lines = show_basic_info( obs, n_obs, n_command_lines);
      show_perturbers( n_command_lines);
      if( debug_level)
         refresh( );
      line_no = n_command_lines + 1;
         {
         char *tptr = tbuff;
         FILE *ofile = fopen_ext( "comments.txt", "tfcw");

         clock_line = 0;
         if( sort_obs_by_code)
            shellsort_r( obs, n_obs, sizeof( OBSERVE), compare_observations,
                                                   &sort_obs_by_code);
         generate_obs_text( obs, n_obs, tbuff, sizeof( tbuff));
         if( make_unicode_substitutions)
            {
                     /* cvt +/- to the Unicode U+00B1,  in UTF-8: */
               text_search_and_replace( tbuff, " +/- ", " \xc2\xb1 ");
                     /* cvt OEM degree symbol (0xf8) to U+00B0,  in UTF-8: */
               text_search_and_replace( tbuff, "\xf8", "\xc2\xb0 ");
            }
         if( sort_obs_by_code)
            shellsort_r( obs, n_obs, sizeof( OBSERVE), compare_observations, NULL);
         while( *tptr)
            {
            size_t i;
            char *tptr2;

            for( i = 0; (unsigned char)tptr[i] >= ' '; i++)
               ;
            tptr[i] = '\0';
            tptr2 = strstr( tptr, "Sigma");
            if( tptr2)
               add_cmd_area( '%', line_no, (unsigned)( tptr2 - tptr), 5);
            put_colored_text( tptr, line_no++, 0, -1,
                     (tptr[i + 1] ? COLOR_OBS_INFO : COLOR_OBS_INFO | A_UNDERLINE));
            fprintf( ofile, "%s\n", tptr);
            tptr += i + 1;
            while( *tptr == 10 || *tptr == 13)
               tptr++;
            if( i < 72 && !clock_line)
               clock_line = line_no - 1;
            }
         fclose( ofile);
         if( debug_level)
            refresh( );
         top_line_orbital_elements = line_no;
         }
         {
         FILE *ifile;

         ifile = fopen_ext( get_file_name( tbuff, elements_filename), "tfcrb");
         if( ifile)
            {
            unsigned iline = 0;
            unsigned right_side_line = 0;
            unsigned right_side_col = 0;
            const unsigned spacing = 4;  /* allow four columns between */
                        /* 'standard' and 'extended' (right-hand) text */
            int elem_color = (bad_elements ?
                              A_BLINK + COLOR_OBS_INFO : COLOR_ORBITAL_ELEMENTS);

            while( iline < 20 && fgets_trimmed( tbuff, sizeof( tbuff), ifile))
               {
               if( show_commented_elements && *tbuff == '#')
                  {
                  if( tbuff[3] != '$' && memcmp( tbuff, "# Find", 6)
                                      && memcmp( tbuff, "# Scor", 6))
                     {
                     put_colored_text( tbuff + 2, line_no + iline, 0, -1, elem_color);
                     iline++;
                     }
                  }
               else if( !show_commented_elements && *tbuff != '#')
                  {
                  char *tptr;
                  int n_chars_to_highlight = 0;
                  const char *reference = get_environment_ptr( "REFERENCE");

                  if( !memcmp( tbuff, "IMPACT", 6))
                     elem_color = COLOR_ATTENTION + A_BLINK;
                  if( !memcmp( tbuff, "Epoch", 5))
                     add_cmd_area( 'e', line_no + iline, 0, 5);
                  if( !memcmp( tbuff + 3, "Peri", 4))
                     add_cmd_area( '+', line_no + iline, 3, 13);
                  tptr = strstr( tbuff, " mean residual");
                  if( tptr)
                     add_cmd_area( 't', line_no + iline,
                                          (int)( tptr - tbuff), 20);
                  tptr = strstr( tbuff, "(J2000 ");
                  if( tptr)
                     add_cmd_area( ALT_N, line_no + iline,
                                          (int)( tptr - tbuff) - 3, 20);
                  tptr = strstr( tbuff, reference);
                  if( tptr)
                     add_cmd_area( ALT_R, line_no + iline,
                                          (int)( tptr - tbuff) - 3, (int)strlen( reference) + 6);
                  tptr = strstr( tbuff, "G ");
                  if( tptr)
                     {
                     char *tptr2 = strstr( tbuff, "H ");

                     if( tptr && tptr2)
                        {
                        add_cmd_area( CTRL( 'G'), line_no + iline,
                                          (int)( tptr - tbuff) - 1, 9);
                        add_cmd_area( -1,         line_no + iline,
                                          (int)( tptr2 - tbuff), 9);
                        }
                     }
                  if( make_unicode_substitutions)
                     {
                     text_search_and_replace( tbuff, " +/- ", " \xc2\xb1 ");
                     text_search_and_replace( tbuff, "^2", "\xc2\xb2");
                     }
                  put_colored_text( tbuff, line_no + iline, 0, -1, elem_color);
                  if( right_side_col < (unsigned)strlen( tbuff) + spacing)
                     right_side_col = (unsigned)strlen( tbuff) + spacing;
                  tptr = strstr( tbuff, "Earth MOID:");
                  if( tptr && atof( tptr + 11) < .01)
                     n_chars_to_highlight = 20;
                  if( !tptr)
                     {
                     tptr = strstr( tbuff, "AMR ");
                     if( tptr)
                        {
                        const double amr = atof( tptr + 4);

                        if( amr < 0. || amr > 90.)    /* unrealistic values */
                           n_chars_to_highlight = 20;
                        }
                     }
                  if( n_chars_to_highlight)
                     put_colored_text( tptr, line_no + iline,
                              count_wide_chars_in_utf8_string( tbuff, tptr),
                              n_chars_to_highlight, COLOR_ATTENTION);
                  iline++;
                  }
               else
                  if( *tbuff == '#' && (unsigned)getmaxx( stdscr)
                                 >= strlen( tbuff + 2) + right_side_col
                              && right_side_line < iline)
                     if( !memcmp( tbuff + 2, "Find", 4)
                       || !memcmp( tbuff + 2, "Eart", 4)
                       ||  (!memcmp( tbuff + 2, "Tiss", 4) && atof( tbuff + 32) < 5.)
                       || !memcmp( tbuff + 2, "Barbee", 6)
                       || !memcmp( tbuff + 2, "Perihe", 6)
                       || !memcmp( tbuff + 2, "Score", 5)
                       || !memcmp( tbuff + 2, "Diame", 5)
                       || !memcmp( tbuff + 2, "MOID", 4))
                        {
                        put_colored_text( tbuff + 2, line_no + right_side_line,
                                right_side_col, -1, elem_color);
                        right_side_line++;
                        }
               }
            line_no += iline;
            fclose( ifile);
            }
         if( debug_level)
            refresh( );
         residual_line = line_no;
         }
      line_no = residual_line;
      show_residual_legend( line_no, residual_format);
      line_no++;
      if( debug_level)
         refresh( );
      if( debug_level > 2)
         debug_printf( "resid legend shown\n");

      if( sort_obs_by_code)
         shellsort_r( obs, n_obs, sizeof( OBSERVE), compare_observations,
                                                     &sort_obs_by_code);

      top_line_residuals = line_no;
/*    if( c != KEY_TIMER) */
         {
         const int n_mpc_codes = find_mpc_color( mpc_color_codes, NULL);
         int lines_available;

         n_stations_shown = show_station_info( NULL, n_obs,
                     line_no, curr_obs, list_codes);

         lines_available = getmaxy( stdscr) - n_stations_shown - line_no;
         n_obs_shown = (lines_available > n_obs ? n_obs : lines_available);
         if( single_obs_selected)
            {
            if( top_obs_shown > curr_obs)
               top_obs_shown = curr_obs;
            if( top_obs_shown < curr_obs - n_obs_shown + 1)
               top_obs_shown = curr_obs - n_obs_shown + 1;
            }

         if( top_obs_shown > n_obs - lines_available)
            top_obs_shown = n_obs - lines_available;
         if( top_obs_shown < 0)
            top_obs_shown = 0;
         for( i = 0; i < n_mpc_codes; i++)
            mpc_color_codes[i].score = 0;
         show_observations( obs, top_obs_shown, top_line_residuals,
                                 residual_format, n_obs_shown, n_obs);
         show_station_info( obs, n_obs,
                     line_no, curr_obs, list_codes);

         show_right_hand_scroll_bar( line_no,
                            n_obs_shown, top_obs_shown, n_obs);
            /* At present,  any unused space between the observation data and
            the station info is left blank.  But as the commented-out lines
            suggest,  it could be filled with... something,  TBD.  */
         *tbuff = '\0';
         for( i = 0; i < lines_available - n_obs_shown; i++)
            {
//          snprintf( tbuff, sizeof( tbuff), "Line %d of %d",
//                      i + 1, lines_available - n_obs_shown);
            put_colored_text( tbuff, top_line_residuals + n_obs_shown + i,
                          0, -1, COLOR_BACKGROUND);
            }
         }
      single_obs_selected = false;
      if( debug_level > 2)
         debug_printf( "resids shown\n");
      if( debug_level)
         refresh( );
      show_final_line( n_obs, curr_obs, COLOR_FINAL_LINE);
      if( sort_obs_by_code)
         shellsort_r( obs, n_obs, sizeof( OBSERVE), compare_observations, NULL);
      if( debug_level)
         refresh( );
      if( *message_to_user)
         {
         int xloc;

         if( add_off_on >= 0)
            strlcat_error( message_to_user, (add_off_on ? " on" : " off"));
         xloc = getmaxx( stdscr) - (int)strlen( message_to_user) - 1;
         put_colored_text( message_to_user, getmaxy( stdscr) - 1,
                           (xloc < 0 ? 0 : xloc), -1, COLOR_MESSAGE_TO_USER);
         }
      add_off_on = -1;
      move( getmaxy( stdscr) - 1, 0);
      if( c == AUTO_REPEATING)
         if( curses_kbhit( ) != ERR)
            {
            extended_getch( );
            c = 0;
            }
      if( c != AUTO_REPEATING)
         {
         int n_ticks_elapsed = 0, n_ticks_mouse_stationary = 0;

         c = 0;
         while( !c)
            {
            if( curses_kbhit( ) != ERR)
               {
               c = extended_getch( );
               if( c == KEY_MOUSE)
                  {
                  get_mouse_data( (int *)&mouse_x, (int *)&mouse_y, (int *)&mouse_z, &button);
                  if( (button & BUTTON4_PRESSED) || button5_pressed)  /* 'wheel up'/'dn' */
                     if( curses_kbhit( ) == KEY_MOUSE)
                        c = 0;      /* avoid duplicated wheel mouse events */
                  }
               }
            else
               {
               if( (n_ticks_elapsed % 20 == 0) && clock_line)
                  {                          /* update clock once a second */
                  const time_t t0 = time( NULL);

                  strlcpy_error( tbuff, ctime( &t0) + 11);
                  tbuff[9] = '\0';
                  put_colored_text( tbuff, clock_line,
                                 getmaxx( stdscr) - 9, 9, COLOR_OBS_INFO);
                  }
               napms( 50);      /* a 'tick' is 50 milliseconds long */
               n_ticks_elapsed++;
               n_ticks_mouse_stationary++;
               if( n_ticks_mouse_stationary == 20 /* = 1 second w/no motion */
                              && _mouse_movements_are_reported)
                  {
                  size_t index;
                  const int cmd = show_hint( mouse_x, mouse_y, &index);
                  int col1 = 0, col2 = 0;

                  if( cmd == -1)
                     {
                     const double abs_mag = calc_absolute_magnitude( obs, n_obs);

                     size_from_h_text( abs_mag, tbuff, 80);
                     show_hint_from_command_area( mouse_x, mouse_y, index, tbuff);
                     }
                  if( cmd == KEY_MOUSE)
                     {
                     *tbuff = '\0';
                     i = mouse_y - top_line_residuals;
                     if( i >= n_obs_shown)
                        {
                        char buff[5];

                        for( i = 0; i < 5; i++)
                           buff[i] = (char)mvinch( mouse_y, i);
                        if( buff[0] == '(' && buff[4] == ')')
                           {
                           buff[4] = '\0';
                           make_observatory_info_text( tbuff, sizeof( tbuff),
                                                       obs, n_obs, buff + 1);
                           col2 = 5;
                           }
                        }
                     else if( i >= 0)
                        {
                        i += top_obs_shown;
                        if( mouse_x < 3)
                           *tbuff = '\0';       /* nothing to do */
                        else if( mouse_x < 20)
                           {
                           double dt, utc;

                           utc = utc_from_td( obs[i].jd, NULL);
                           dt = current_jd( ) - utc;

                           full_ctime( tbuff, utc, CALENDAR_JULIAN_GREGORIAN
                                      | FULL_CTIME_YMD | FULL_CTIME_LEADING_ZEROES
                                      | FULL_CTIME_MICRODAYS);
                           strcat( tbuff, " = ");
                           full_ctime( tbuff + strlen( tbuff), utc, FULL_CTIME_TIME_ONLY);
                           if( dt < 1.)
                              snprintf_append( tbuff, sizeof( tbuff), "\n%.2f hours ago", dt * 24.);
                           else if( dt < 7.)
                              snprintf_append( tbuff, sizeof( tbuff), "\n%.2f days ago", dt);
                           if( utc > 2400001.0)
                              snprintf_append( tbuff, sizeof( tbuff), "\nMJD %.6f", utc - 2400000.5);
                           else
                              snprintf_append( tbuff, sizeof( tbuff), "\nJD %.6f", utc);
                           col1 = 3;
                           col2 = 20;
                           }
                        else if( mouse_x < 42)
                           {
                           char formatted[82];

                           recreate_observation_line( formatted, obs + i, 3 << 14);   /* HH MM SS.s */
                           formatted[56] = '\0';
                           strlcpy_error( tbuff, formatted + 32);
                           strlcat_error( tbuff, "\n");
                           recreate_observation_line( formatted, obs + i, 1 << 14);   /* decimal degrees */
                           formatted[56] = '\0';
                           strlcat_error( tbuff, formatted + 32);
                           col1 = 20;
                           col2 = 44;
                           }
                        if( mouse_x >= 53 && mouse_x < 59 && obs[i].computed_mag != BLANK_MAG)
                           {
                           snprintf_err( tbuff, sizeof( tbuff), "computed mag %.2f",
                                                            obs[i].computed_mag);
                           col1 = 53;
                           col2 = 59;
                           }
                        }
                     if( *tbuff)
                        show_hint_text( mouse_x, mouse_y, mouse_y, col1, col2 - col1, tbuff);
                     }
                  }
               }
            if( c == KEY_MOUSE && (button & REPORT_MOUSE_POSITION))
               c = n_ticks_mouse_stationary = 0;      /* suppress mouse moves */
            }
         auto_repeat_full_improvement = 0;
         }
      for( i = 0; key_remaps[i]; i += 2)
         if( c == key_remaps[i])
            c = key_remaps[i + 1];
      if( (button & BUTTON4_PRESSED) || button5_pressed)  /* 'wheel up'/'dn' */
         {
         mouse_wheel_motion = ((button & BUTTON_CTRL) ? 5 : 1);
         if( !button5_pressed)         /* 'wheel down' */
            mouse_wheel_motion = -mouse_wheel_motion;
         }
      else
         mouse_wheel_motion = 0;

      *message_to_user = '\0';
      if( c == KEY_MOUSE && !(button & REPORT_MOUSE_POSITION))
         {
         c = find_command_area( mouse_x, mouse_y, NULL);
         if( c == KEY_OBSCODE_CLICKED)
            {
            curr_obs = top_obs_shown + (mouse_y - top_line_residuals);
            c = ( button & (BUTTON1_RELEASED | BUTTON1_CLICKED) ? 'T' : 'X');
            }
         if( c == ALT_X && mouse_wheel_motion)
            c = (mouse_wheel_motion < 0 ? KEY_F( 4) : KEY_F( 5));
         if( c == 't' && mouse_wheel_motion)
            c = (mouse_wheel_motion < 0 ? KEY_CYCLE_RESID_DISPLAY_UP
                                        : KEY_CYCLE_RESID_DISPLAY_DN);
         if( c < 0)      /* informational hint text; no command */
            c = KEY_MOUSE;
         }

      if( c == KEY_MOUSE && !(button & REPORT_MOUSE_POSITION))
         {
         const unsigned station_start_line = getmaxy( stdscr) - n_stations_shown;

         if( debug_mouse_messages)
            snprintf_err( message_to_user, sizeof( message_to_user),
                              "x=%d y=%d z=%d button=%lx",
                              mouse_x, mouse_y, mouse_z, (unsigned long)button);
         if( mouse_wheel_motion && mouse_y < station_start_line)
            {
            int n_selected = 0;

            top_obs_shown += mouse_wheel_motion;
            for( i = 0; i < n_obs; i++)
               if( obs[i].flags & OBS_IS_SELECTED)
                  n_selected++;
            if( n_selected == 1)
               {
               curr_obs += mouse_wheel_motion;
               single_obs_selected = true;
               }
            }
         else if( mouse_y >= station_start_line)
            {
            int c1;
            const char *search_code =
                    mpc_color_codes[mouse_y - station_start_line].code;
            char *tptr;

            strlcpy_error( tbuff, get_find_orb_text( 2050));
            tptr = tbuff + 1;
            for( i = 0; i < (int)(list_codes & 3); i++)
                tptr = strstr( tptr, "( )") + 1;
            *tptr = '*';
            text_search_and_replace( tbuff, "$", search_code);
            if( mouse_wheel_motion)
               c1 = (mouse_wheel_motion > 0 ? KEY_F( 4) : KEY_F( 5));
            else
               {
               help_file_name = "mpc_area.txt";
               c1 = full_inquire( tbuff, NULL, 0, COLOR_MENU, mouse_y, mouse_x);
               }
            if( c1 >= KEY_F(1) && c1 <= KEY_F(3))
               list_codes = c1 - KEY_F(1);
            else if( c1 == KEY_F(4) || c1 == KEY_F(5))
               {                              /* find next or prev */
               for( i = 0; i < n_obs; i++)    /* obs from this code */
                  {
                  curr_obs += (c1 == KEY_F( 5) ? n_obs - 1 : 1);
                  curr_obs %= n_obs;
                  if( !strcmp( obs[curr_obs].mpc_code, search_code))
                     break;
                  }
               single_obs_selected = true;
               }
            else if( c1 == KEY_F(6))
               c = ALT_X;
            }
         else if( mouse_y >= top_line_residuals)
            {
            const unsigned max_x = getmaxx( stdscr);

            if( mouse_x == max_x - 1)    /* clicked on 'scroll bar' right of obs */
               {
               int dir = 0;

               if( mouse_y == top_line_residuals)
                  dir = -1;        /* similar to mouse wheel up */
               else if( mouse_y == top_line_residuals + n_obs_shown - 1)
                  dir =  1;        /* similar to mouse wheel down */
               else
                  {              /* clicked on scroll bar */
                  top_obs_shown =
                          (mouse_y - top_line_residuals) * n_obs / n_obs_shown;
                  top_obs_shown -= n_obs_shown / 2;  /* Center display */
                  }
               if( button & BUTTON1_DOUBLE_CLICKED)
                  dir += dir;
               if( button & BUTTON1_TRIPLE_CLICKED)
                  dir *= 3;
               top_obs_shown += dir;
               }
            else
               {              /* clicked among the observations */
               int new_curr = top_obs_shown + (mouse_y - top_line_residuals);

               if( new_curr < n_obs)  /* "normal" click in the observations area */
                  {
                  if( button & BUTTON1_DOUBLE_CLICKED)
                     obs[new_curr].is_included ^= 1;
                  else if( button & BUTTON_CTRL)
                     obs[new_curr].flags ^= OBS_IS_SELECTED;
                  else        /* "ordinary",  unshifted or ctrled click */
                     {
                     int idx1 = new_curr, idx2 = new_curr;

                     if( button & (BUTTON1_RELEASED | BUTTON2_RELEASED | BUTTON3_RELEASED
                                          | BUTTON_SHIFT))
                        {                          /* selected a range of obs */
                        idx1 = min( curr_obs, new_curr);
                        idx2 = max( curr_obs, new_curr);
                        }
                     if( button & BUTTON_CTRL)
                        obs[new_curr].flags ^= OBS_IS_SELECTED;
                     else if( button & BUTTON_SHIFT)
                        {
                        int n_selected = 0, n_unselected;

                        for( i = idx1; i <= idx2; i++)
                           if( obs[i].flags & OBS_IS_SELECTED)
                              n_selected++;
                        n_unselected = (idx2 - idx1 + 1) - n_selected;
                        for( i = idx1; i <= idx2; i++)
                           if( n_selected > n_unselected)
                              obs[i].flags &= ~OBS_IS_SELECTED;
                           else
                              obs[i].flags |= OBS_IS_SELECTED;
                        }
                     else if( button & (BUTTON1_RELEASED | BUTTON1_CLICKED))
                        {
                        for( i = 0; i < n_obs; i++)
                           if( i >= idx1 && i <= idx2)
                              obs[i].flags |= OBS_IS_SELECTED;
                           else
                              obs[i].flags &= ~OBS_IS_SELECTED;
                        }
                     if( button & (BUTTON2_RELEASED | BUTTON2_CLICKED
                                 | BUTTON3_RELEASED | BUTTON3_CLICKED))
                        {                 /* right or middle button click/release */
                        char *search_code = obs[new_curr].mpc_code;

                        strlcpy_error( tbuff, get_find_orb_text( 2022));
                        text_search_and_replace( tbuff, "$", search_code);
                        i = full_inquire( tbuff, NULL, 0, COLOR_MENU, mouse_y, mouse_x);
                        switch( i)
                           {
                           case KEY_F( 1) :       /* toggle obs */
                              c = 'x';
                              break;
                           case KEY_F( 2) :       /* set uncertainty */
                              c = '%';
                              break;
                           case KEY_F( 3) :       /* select all this obscode */
                              {
                              int n_selected = 0, n_deselected = 0;

                              for( i = 0; i < n_obs; i++)
                                 if( !FSTRCMP( obs[i].mpc_code, search_code))
                                    {
                                    if( obs[i].flags & OBS_IS_SELECTED)
                                       n_selected++;
                                    else
                                       n_deselected++;
                                    }
                              for( i = 0; i < n_obs; i++)
                                 if( !FSTRCMP( obs[i].mpc_code, search_code))
                                    {
                                    if( n_selected < n_deselected)
                                       obs[i].flags |= OBS_IS_SELECTED;
                                    else
                                       obs[i].flags &= ~OBS_IS_SELECTED;
                                    }
                              snprintf_err( message_to_user, sizeof( message_to_user),
                                      "%d obs from (%s) %sselected\n",
                                      n_selected + n_deselected,
                                      search_code, (n_selected < n_deselected) ? "" : "de");
                              }
                              break;
                           case KEY_F( 4) :       /* prev this obscode */
                           case KEY_F( 5) :       /* next this obscode */
                              c = i;         /* works out perfectly,  by luck */
                              break;
                           }
                        }
                     }
                  curr_obs = new_curr;
                  }
               }
            }
         else if( mouse_y > n_command_lines
               && mouse_y < top_line_orbital_elements)   /* in obs details area: */
            c = ALT_Q;         /* toggle display header/'traditional' data */
         else if( c == KEY_MOUSE
                  && mouse_y >= top_line_orbital_elements)
            {
            if( button & (BUTTON2_RELEASED | BUTTON2_CLICKED
                        | BUTTON3_RELEASED | BUTTON3_CLICKED))
               {                 /* right or middle button click/release */
               setup_elements_dialog( tbuff, orbit_constraints);
               help_file_name = "elem_pop.txt";
               i = full_inquire( tbuff, NULL, 0,
                               COLOR_MENU, mouse_y, mouse_x);
               c = KEY_MOUSE;
               switch( i)
                  {
                  case KEY_F( 1) :     /* elems->clipboard */
                     c = ALT_I;
                     break;
                  case KEY_F( 2) :     /* save elems       */
                     c = 's';
                     break;
                  case KEY_F( 3) :     /* reference        */
                     c = ALT_R;
                     break;
                  case KEY_F( 4) :     /* epoch            */
                     c = 'e';
                     break;
                  case KEY_F( 5) :     /* frame            */
                     c = ALT_N;
                     break;
                  case KEY_F( 6) :     /* constraints      */
                     c = 'l';
                     break;
                  case KEY_F( 7) :     /* add digit prec   */
                     c = 'p';
                     break;
                  case KEY_F( 8) :     /* sub digit prec   */
                     c = 'P';
                     break;
                  case KEY_F( 9) :     /* elem centre      */
                     c = '+';
                     break;
                  default:
                     break;
                  }
               }
            else
               c = CTRL( 'B');      /* toggle commented elems */
            }
         }           /* end of mouse-handling code */

      if( c >= '1' && c <= '9')
         perturbers ^= (1 << (c - '0'));
#ifdef ALT_0
      else if( c >= ALT_0 && c <= ALT_9)
         top_obs_shown = (n_obs - 1) * (c - ALT_0) / 10;
#endif
      else switch( c)
         {
         case '0':
            perturbers ^= 1024;
            break;
         case KEY_C1:
         case KEY_END:
            curr_obs = n_obs - 1;
            single_obs_selected = true;
            break;
         case KEY_A1:
         case KEY_HOME:
            curr_obs = 0;
            single_obs_selected = true;
            break;
         case KEY_UP:
#ifdef KEY_A2
         case KEY_A2:
#endif
         case KEY_LEFT:
#ifdef KEY_B1
         case KEY_B1:
#endif
            curr_obs--;
            single_obs_selected = true;
            break;
         case KEY_DOWN:
#ifdef KEY_C2
         case KEY_C2:
#endif
         case KEY_RIGHT:
#ifdef KEY_B3
         case KEY_B3:
#endif
            curr_obs++;
            single_obs_selected = true;
            break;
         case KEY_C3:         /* "PgDn" = lower right key in keypad */
         case KEY_NPAGE:
            curr_obs += n_obs_shown;
            single_obs_selected = true;
            break;
         case KEY_A3:         /* "PgUp" = upper right key in keypad */
         case KEY_PPAGE:
            curr_obs -= n_obs_shown;
            single_obs_selected = true;
            break;
         case ALT_W:
            {
            extern bool use_sigmas;

            use_sigmas = !use_sigmas;
            strlcpy_error( message_to_user, get_find_orb_text( 19));
                         /* "Weighting of posns/mags/times is" */
            add_off_on = use_sigmas;
            }
            break;
#ifdef KEY_SUP
         case KEY_SUP:
#endif
         case KEY_SR:
            obs[curr_obs].is_included ^= 1;
            for( i = 0; i < curr_obs; i++)
               obs[i].is_included = obs[curr_obs].is_included;
            strlcpy_error( message_to_user, get_find_orb_text( 20));
            add_off_on = obs[curr_obs].is_included;
            break;
         case KEY_F(16):     /* Shift-F4:  select another object to add in */
            if( (i = select_object_in_file( ids, n_ids)) >= 0)
               {
               FILE *ifile = fopen( ifilename, "rb");

               obs = add_observations( ifile, obs, ids + i, &n_obs);
               if( debug_level)
                  printf( "Now have %d obs\n", n_obs);
               if( mpc_color_codes)
                  free( mpc_color_codes);
               if( max_mpc_color_codes)
                  mpc_color_codes = find_mpc_color_codes( n_obs, obs,
                             max_mpc_color_codes);
               if( debug_level)
                  debug_printf( "got color codes; ");
               for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
                  ;
               curr_obs = i;
               set_locs( orbit, curr_epoch, obs, n_obs);
               update_element_display = 1;
               clear( );
               }
            break;
#ifdef KEY_SDOWN
         case KEY_SDOWN:
#endif
         case KEY_SF:            /* turn on/off all obs after curr one */
            obs[curr_obs].is_included ^= 1;
            for( i = curr_obs; i < n_obs; i++)
               obs[i].is_included = obs[curr_obs].is_included;
            strlcpy_error( message_to_user, "All subsequent observations toggled");
            add_off_on = obs[curr_obs].is_included;
            break;
         case KEY_F(3):          /* turn on/off all obs w/same observatory ID */
            obs[curr_obs].is_included ^= 1;
            for( i = 0; i < n_obs; i++)
               if( !FSTRCMP( obs[i].mpc_code, obs[curr_obs].mpc_code))
                  obs[i].is_included = obs[curr_obs].is_included;
            strlcpy_error( message_to_user, "All observations from xxx toggled");
            FMEMCPY( message_to_user + 22, obs[curr_obs].mpc_code, 3);
            add_off_on = obs[curr_obs].is_included;
            break;
         case KEY_F(4):          /* find prev obs from this observatory */
         case KEY_F(5):          /* find next obs from this observatory */
            {
            const int dir = (c == KEY_F(4) ? n_obs - 1 : 1);
            int new_obs = (curr_obs + dir) % n_obs;

            while( new_obs != curr_obs &&
                        FSTRCMP( obs[new_obs].mpc_code, obs[curr_obs].mpc_code))
               new_obs = (new_obs + dir) % n_obs;
            curr_obs = new_obs;
            single_obs_selected = true;
            }
            break;
         case KEY_F(6):          /* find prev excluded obs */
         case KEY_F(7):          /* find next excluded obs */
            {
            const int dir = (c == KEY_F(6) ? n_obs - 1 : 1);
            int new_obs = (curr_obs + dir) % n_obs;

            while( new_obs != curr_obs && obs[new_obs].is_included)
               new_obs = (new_obs + dir) % n_obs;
            curr_obs = new_obs;
            single_obs_selected = true;
            }
            break;
         case '*':
         case '^':
            non_grav_menu( message_to_user);
            if( *message_to_user)      /* new force model selected */
               {
               extern int force_model;

               for( i = 6; i < n_orbit_params; i++)
                  orbit[i] = 0.;
               if( force_model == FORCE_MODEL_DELTA_V)
                  {        /* we need a starting estimate of the maneuver time */
                  if( !inquire( get_find_orb_text( 2099), tbuff, sizeof( tbuff),
                            COLOR_DEFAULT_INQUIRY) && *tbuff)
                     orbit[9] = get_time_from_string( obs->jd, tbuff,
                             FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, NULL);
                  if( orbit[9] < obs->jd || orbit[9] > obs[n_obs - 1].jd)
                     force_model = FORCE_MODEL_NO_NONGRAVS;
                  }
               }
            break;
#ifndef _WIN32
         case KEY_F(8):     /* show original screens */
            full_endwin( );
            extended_getch( );
            restart_curses( );
            break;
#endif
         case 'a': case 'A':
            perturbers ^= (7 << 20);
            strlcpy_error( message_to_user, "Asteroids toggled");
            add_off_on = (perturbers >> 20) & 1;
            break;
         case '!':
            perturbers = ((perturbers == 0x3fe) ? 0 : 0x3fe);
            break;
         case KEY_F(9):           /* find start of included arc */
            for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
               ;
            curr_obs = i;
            single_obs_selected = true;
            break;
         case KEY_F(10):          /* find end of included arc */
            for( i = n_obs - 1; i > 0 && !obs[i].is_included; i--)
               ;
            curr_obs = i;
            single_obs_selected = true;
            break;
         case CTRL( 'G'):
            if( !inquire( "Slope parameter( G): ",
                                tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY)
                        && *tbuff)
               {
               extern double asteroid_magnitude_slope_param;

               asteroid_magnitude_slope_param = atof( tbuff);
               }
            break;
         case KEY_F(11):
            auto_repeat_full_improvement ^= 1;
            strlcpy_error( message_to_user, "Automatic full improvement repeat is");
            add_off_on = auto_repeat_full_improvement;
            break;
         case 'D':
            {
            extern const char *observe_filename;

            create_obs_file( obs, n_obs, 0, residual_format);
            show_a_file( get_file_name( tbuff, observe_filename), 0);
            }
            break;
         case 'd':
            {
            extern const char *residual_filename;

            create_resid_file( obs, n_obs, ifilename, residual_format);
            show_a_file( get_file_name( tbuff, residual_filename), 0);
            }
            break;
         case 'e': case'E':
            {
            double new_jd = epoch_shown;

            help_file_name = "timehelp.txt";
            if( inquire( "Enter new epoch,  as YYYY MM DD, or JD,  or 'now':",
                             tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
               *tbuff = '\0';       /* cancelled entry */
            else if( !tbuff[1] && strchr( "sme", tbuff[0]))
               {
               int first, last;

               first = find_first_and_last_obs_idx( obs, n_obs, &last);
               if( *tbuff == 's')
                  epoch_shown = obs[first].jd;
               if( *tbuff == 'e')
                  epoch_shown = obs[last].jd;
               if( *tbuff == 'm')
                  epoch_shown = (obs[first].jd + obs[last].jd) / 2.;
               epoch_shown = floor( epoch_shown) + .5;
               }
            else if( *tbuff == 'p')
               {
               ELEMENTS elem;
               double orbit2[MAX_N_PARAMS];

               memcpy( orbit2, orbit, n_orbit_params * sizeof( double));
               integrate_orbit( orbit2, curr_epoch, epoch_shown);
               find_relative_orbit( epoch_shown, orbit2, &elem,
                        (tbuff[1] == 'g' ? 3 : 0));
               epoch_shown = floor( elem.perih_time * 10. + 0.5) / 10.;
               }
            else if( extract_date( tbuff, &new_jd) >= 0)
               if( new_jd > minimum_jd && new_jd < maximum_jd)
                  epoch_shown = new_jd;
            }
            update_element_display = 1;
            break;
         case 'f': case 'F':        /* do a "full improvement" step */
         case '|':                  /* Monte Carlo step */
         case CTRL( 'A'):           /* statistical ranging */
         case AUTO_REPEATING:
            {
            double *stored_ra_decs = NULL;
            int err = 0;
            extern int using_sr;
            const clock_t t0 = clock( );

            if( c == '|' || c == CTRL( 'A'))
               {
               const double rms = compute_rms( obs, n_obs);

               set_statistical_ranging( c == CTRL( 'A'));
               c = '|';
               if( !inquire( "Gaussian noise level (sigmas): ",
                             tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
                  {
                  noise_in_sigmas = atof( tbuff);
                  if( noise_in_sigmas)
                     {
                     max_monte_rms =
                              sqrt( noise_in_sigmas * noise_in_sigmas
                                           + rms * rms);
                     c = AUTO_REPEATING;
                     }
                  }
               }
            if( c == AUTO_REPEATING)
               stored_ra_decs =
                   add_gaussian_noise_to_obs( n_obs, obs, noise_in_sigmas);
            push_orbit( curr_epoch, orbit);
            if( c == AUTO_REPEATING && using_sr)
               {
               sr_orbit_t sr;

               find_nth_sr_orbit( &sr, obs, n_obs, monte_carlo_object_count);
               memcpy( orbit, sr.orbit, 6 * sizeof( double));
               adjust_herget_results( obs, n_obs, orbit);
                        /* epoch is that of first included observation: */
               get_epoch_range_of_included_obs( obs, n_obs, &curr_epoch, NULL);
               }
            else
               {
               double saved_orbit[MAX_N_PARAMS];
               const double mid_epoch = mid_epoch_of_arc( obs, n_obs);
//             const double mid_epoch = curr_epoch;
//             const double mid_epoch = epoch_shown;

               memcpy( saved_orbit, orbit, n_orbit_params * sizeof( double));
               integrate_orbit( orbit, curr_epoch, mid_epoch);
                        /* Only request sigmas for i=1 (last pass... */
                        /* _only_ pass for a real 'full improvement') */
               for( i = (c == AUTO_REPEATING ? 5 : 1); i && !err; i--)
                  {
                  int sigma_type = NO_ORBIT_SIGMAS_REQUESTED;

                  if( i == 1 && c != AUTO_REPEATING &&
                              (element_format & ELEM_OUT_ALTERNATIVE_FORMAT))
                     sigma_type = ORBIT_SIGMAS_REQUESTED;
                  err = full_improvement( obs, n_obs, orbit, mid_epoch,
                              orbit_constraints, sigma_type, epoch_shown);
                  }
               if( err)    /* Full Monte Carlo isn't working.  Let's try SR, */
                  {        /* & recover from error by using the saved orbit */
                  debug_printf( "Full improvement fail %d\n", err);
                  set_statistical_ranging( 1);
                  memcpy( orbit, saved_orbit, n_orbit_params * sizeof( double));
                  }
               else
                  curr_epoch = mid_epoch;
               }
            get_r1_and_r2( n_obs, obs, &r1, &r2);
            if( c == AUTO_REPEATING)
               {
               restore_ra_decs_mags_times( n_obs, obs, stored_ra_decs);
               free( stored_ra_decs);
               }
            update_element_display = (err ? 0 : 1);
            if( c == AUTO_REPEATING && !err)
               {
               ELEMENTS elem;
               double rel_orbit[MAX_N_PARAMS], orbit2[MAX_N_PARAMS];
               int curr_planet_orbiting;
               extern int n_clones_accepted;

               is_monte_orbit = true;
               memcpy( orbit2, orbit, n_orbit_params * sizeof( double));
               integrate_orbit( orbit2, curr_epoch, epoch_shown);
               snprintf_err( message_to_user, sizeof( message_to_user),
                           (using_sr ? "Stat Ranging %d/%d" : "Monte Carlo %d/%d"),
                           n_clones_accepted,
                           monte_carlo_object_count);
               curr_planet_orbiting = find_best_fit_planet( epoch_shown,
                                  orbit2, rel_orbit);
               if( !monte_carlo_object_count)
                  planet_orbiting = curr_planet_orbiting;

               if( planet_orbiting == curr_planet_orbiting)
                  {
                  elem.gm = get_planet_mass( planet_orbiting);
                  calc_classical_elements( &elem, rel_orbit, epoch_shown, 1);
                  add_monte_orbit( monte_data, &elem, monte_carlo_object_count);
                  }
               if( monte_carlo_object_count > 3)
                  {
                  double sigmas[MONTE_N_ENTRIES];
                  FILE *monte_file = fopen_ext( get_file_name( tbuff, "monte.txt"), "tfcwb");

                  fprintf( monte_file,
                          "Computed from %d orbits around object %d\n",
                          monte_carlo_object_count, planet_orbiting);
                  compute_monte_sigmas( sigmas, monte_data,
                                             monte_carlo_object_count);
                  dump_monte_data_to_file( monte_file, sigmas,
                        elem.major_axis, elem.ecc, planet_orbiting);
                  fclose( monte_file);
                  }
               }
            else
               {
               strlcpy_error( message_to_user,
                               (err ? "Full step FAILED" : "Full step taken"));
               if( curses_kbhit_without_mouse( ) > 0)
                  {
                  extended_getch( );
                  strlcpy_error( message_to_user, "User interrupted step");
                  }
               snprintf_append( message_to_user, sizeof( message_to_user), "(%.5f s)",
                      (double)( clock( ) - t0) / (double)CLOCKS_PER_SEC);
               }
            }
            break;
         case 'g': case 'G':        /* do a method of Gauss soln */
            {
            double new_epoch;

            perturbers = 0;
            push_orbit( curr_epoch, orbit);
            if( c == 'G' && ! inquire( "Initial Gauss rho: ", tbuff,
                       sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
               {
               orbit[0] = atof( tbuff);
               gauss_soln = -1;
               }
            new_epoch = convenient_gauss( obs, n_obs, orbit, 1., gauss_soln);
            if( !new_epoch && gauss_soln > 0)
               {
               gauss_soln = 0;
               new_epoch = convenient_gauss( obs, n_obs, orbit, 1., gauss_soln);
               }
            gauss_soln++;
            if( new_epoch)
               {
               curr_epoch = new_epoch;
               set_locs( orbit, curr_epoch, obs, n_obs);
               get_r1_and_r2( n_obs, obs, &r1, &r2);
               update_element_display = 1;
               strlcpy_error( message_to_user, "Gauss solution found");
               }
            else
               strlcpy_error( message_to_user, "Gauss method failed!");
            }
            break;
         case '#':
         case '/':
            push_orbit( curr_epoch, orbit);
                     /* epoch is that of first valid observation: */
            for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
                ;
            if( c == '#')
               {
               simplex_method( obs + i, n_obs - i, orbit, r1, r2, orbit_constraints);
               strlcpy_error( message_to_user, "Simplex method used");
               }
            else
               {
               integrate_orbit( orbit, curr_epoch, obs[i].jd);
               superplex_method( obs + i, n_obs - i, orbit, orbit_constraints);
               strlcpy_error( message_to_user, "Superplex method used");
               }
//          integrate_orbit( orbit, obs[i].jd, curr_epoch);
            curr_epoch = obs[i].jd;
            get_r1_and_r2( n_obs, obs, &r1, &r2);    /* orb_func.cpp */
            update_element_display = 1;
            break;
         case 'h':               /* do a method of Herget step */
         case 'H':               /* same, but don't linearize */
         case ':':               /* just linearize */
            {
            int err = 0;

            push_orbit( curr_epoch, orbit);
            if( c != ':')
               {
               double d_r1, d_r2;

               err = herget_method( obs, n_obs, r1, r2, orbit, &d_r1, &d_r2,
                                          orbit_constraints);
               if( !err)
                  {
                  r1 += d_r1;
                  r2 += d_r2;
                  herget_method( obs, n_obs, r1, r2, orbit, NULL, NULL, NULL);
                  }
               }
            if( c != 'H')
               if( !err || (c == ':' && n_obs == 2))
                  err = adjust_herget_results( obs, n_obs, orbit);
                     /* epoch is that of first valid observation: */
            if( !err)
               {
               for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
                  ;
               curr_epoch = obs[i].jd;
               update_element_display = 1;
               strlcpy_error( message_to_user, (c == ':') ? "Orbit linearized" :
                                                  "Herget step taken");
               }
            else
               strlcpy_error( message_to_user, (c == ':') ? "Linearizing FAILED" :
                                                  "Herget step FAILED");
            }
            break;
         case '<':
            push_orbit( curr_epoch, orbit);
            herget_method( obs, n_obs, r1, r2, orbit, NULL, NULL, NULL);
            for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
               ;
            curr_epoch = obs[i].jd;
            update_element_display = 1;
            snprintf_err( message_to_user, sizeof( message_to_user),
                                     "Radii set: %f %f", r1, r2);
            break;
         case 'l': case 'L':
            if( !*orbit_constraints)
               {
               if( inquire(
   "Enter limits on the orbit (e.g.,  'e=0' or 'q=2.3' or 'q=.7,P=1.4').\n"
   "Constraints can be placed on e, q, Q, P, a, n, O, o, or i.",
                     orbit_constraints, sizeof( orbit_constraints),
                     COLOR_DEFAULT_INQUIRY))
                  *orbit_constraints = '\0';    /* cancelled constraints */
               }
            else
               {
               *orbit_constraints = '\0';
               strlcpy_error( message_to_user, "Orbit is now unconstrained");
               }
            if( !strcmp( orbit_constraints, "K"))
               strlcpy_error( orbit_constraints, "e=1,b=282.81,l=35.22");
            if( !strcmp( orbit_constraints, "Me"))
               strlcpy_error( orbit_constraints, "e=1,i=72,O=72");
            if( !strcmp( orbit_constraints, "Ma"))
               strlcpy_error( orbit_constraints, "e=1,i=26,O=81"); /* q=.049? */
            compute_effective_solar_multiplier( orbit_constraints);
            break;
         case 'm': case 'M':
            create_obs_file( obs, n_obs, 0, residual_format);
            create_ephemeris( orbit, curr_epoch, obs, n_obs, obj_name,
                           ifilename, residual_format);
            break;
         case 'o':             /* select a new file */
            get_new_file = get_new_object = true;
            *ifilename = '\0';
            break;
         case 'n': case 'N':   /* select a new object from the input file */
            get_new_object = true;
            update_element_display = 1;
            pop_all_orbits( );
            break;
         case 'P':         /* show one less digit of precision in elements */
         case 'p':         /* show one more digit of precision in elements */
            if( c == 'P' && element_precision > 1)
               element_precision--;
            else if( c == 'p' && element_precision < 15)
               element_precision++;
            update_element_display = 1;
            snprintf( message_to_user, sizeof( message_to_user),
                                    "%d digits\n", element_precision);
            break;
         case CTRL( 'P'):
            if( !inquire( "Blunder probability: ", tbuff, sizeof( tbuff),
                            COLOR_DEFAULT_INQUIRY) && *tbuff)
               {
               extern double probability_of_blunder;

               probability_of_blunder = atof( tbuff);
               }
            break;
         case CTRL( 'F'):
            {
            extern double **eigenvects;
            static double total_sigma = 0.;

            if( eigenvects && !inquire( "Enter sigma: ", tbuff,
                        sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
               {
               if( *tbuff == '0')
                  total_sigma = 0.;
               else
                  {
                  const double sigma = atof( tbuff);

                  compute_variant_orbit( orbit, orbit, sigma);
                  set_locs( orbit, curr_epoch, obs, n_obs);
                  total_sigma += sigma;
                  update_element_display = 1;
                  snprintf( message_to_user,  sizeof( message_to_user),
                           "Epoch = %f (%f); sigma %f",
                           curr_epoch, epoch_shown, total_sigma);
                  }
               }
            }
            break;
         case ALT_F:
            {
            extern double **eigenvects;

            if( eigenvects && eigenvects[0])
               {
               const double n_sigmas = improve_along_lov( orbit, curr_epoch,
                     eigenvects[0], n_orbit_params, n_obs, obs);

               snprintf( message_to_user, sizeof( message_to_user),
                                       "Adjusted by %f sigmas", n_sigmas);
               update_element_display = 1;
               }
            else
               strlcpy_error( message_to_user, "No LOV along which to improve");
            }
            break;
         case CTRL( 'D'):
            if( !inquire( "Number SR orbits: ", tbuff, sizeof( tbuff),
                            COLOR_DEFAULT_INQUIRY))
               {
               unsigned max_orbits, n_improvements = 0;

               if( sscanf( tbuff, "%u,%u", &max_orbits, &n_improvements))
                  {
                  sr_orbit_t *orbits = (sr_orbit_t *)calloc( max_orbits, sizeof( sr_orbit_t));
                  int n_found;

                  for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
                     ;
                  n_found = get_sr_orbits( orbits, obs + i, n_obs - i, 0, max_orbits,
                           86400., 1., 1);
                  snprintf_err( message_to_user, sizeof( message_to_user),
                                   "%d orbits computed: best score=%.3f\n",
                                   n_found, orbits[0].score);
                  debug_printf( "%d found, will do %u improvements\n",
                              n_found, n_improvements);
                  if( n_found > 0)
                     {
                     while( n_improvements--)
                        improve_sr_orbits( orbits, obs + i,
                                    n_obs - i, n_found, 1., 0);
                     curr_epoch = obs[i].jd;
                     for( i = 0; i < n_found; i++)
                        if( orbits[i].score < orbits[0].score)
                           {
                           const sr_orbit_t torbit = orbits[0];

                           orbits[0] = orbits[i];
                           orbits[i] = torbit;
                           }
                     push_orbit( curr_epoch, orbit);
                     memcpy( orbit, orbits[0].orbit, 6 * sizeof( double));
                     update_element_display = 1;
                     set_locs( orbit, curr_epoch, obs, n_obs);
                     show_a_file( "sr_elems.txt", 0);
                     }
                  free( orbits);
                  }
               }
            break;
         case CTRL( 'R'):
            {
            extern double sr_min_r, sr_max_r;
            extern int using_sr;

            if( !inquire( "Enter SR R1, R2: ", tbuff, sizeof( tbuff),
                            COLOR_DEFAULT_INQUIRY) &&
                   3 == sscanf( tbuff, "%lf %lf %lf", &sr_min_r, &sr_max_r, &max_monte_rms))
               using_sr = 1;
            }
            break;
         case 'r': case 'R':
            if( !inquire( "Enter new R1, R2: ", tbuff, sizeof( tbuff),
                            COLOR_DEFAULT_INQUIRY)
                   && sscanf( tbuff, "%lf%n", &r1, &i) == 1)
               {
               int j;

               while( tbuff[i] == ' ')
                  i++;
               if( tolower( tbuff[i]) == 'k')
                  {
                  r1 /= AU_IN_KM;
                  i++;
                  if( tolower( tbuff[i]) == 'm')
                     i++;
                  }
               while( tbuff[i] == ',' || tbuff[i] == ' ')
                  i++;
               if( !tbuff[i])    /* only one distance entered */
                  r2 = r1;
               else if( sscanf( tbuff + i, "%lf%n", &r2, &j) == 1)
                  {
                  i += j;
                  while( tbuff[i] == ' ')
                     i++;
                  if( tolower( tbuff[i]) == 'k')
                     r2 /= AU_IN_KM;
                  }
               snprintf_err( message_to_user, sizeof( message_to_user),
                                                "R1 = %f; R2 = %f", r1, r2);
               }
            else if( *tbuff == 'g')
               {
               extern double general_relativity_factor;

               general_relativity_factor = atof( tbuff + 1);
               }
            else if( *tbuff == 'l')
               {
               extern double levenberg_marquardt_lambda;
               const double orig_lambda = levenberg_marquardt_lambda;

               if( tbuff[1] == '*')
                  levenberg_marquardt_lambda *= atof( tbuff + 2);
               else if( tbuff[1] == '-')
                  levenberg_marquardt_lambda = pow( 0.1, atof( tbuff + 2));
               else
                  levenberg_marquardt_lambda = atof( tbuff + 1);
               snprintf( message_to_user, sizeof( message_to_user),
                     "Lambda changed from %.3e to %.3e", orig_lambda,
                                             levenberg_marquardt_lambda);
               }
            update_element_display = 1;
            break;
         case 's': case 'S':     /* save orbital elements to a file */
            {
            if( !inquire( "Enter filename for saving elements: ", tbuff,
                              sizeof( tbuff), COLOR_DEFAULT_INQUIRY) && *tbuff)
               {
               FILE *ofile;
               char filename[85];
               double orbit2[MAX_N_PARAMS];

#ifndef _WIN32
               fix_home_dir( tbuff);
#endif
               sscanf( tbuff, "%79s", filename);
               ofile = fopen( filename, "wb");
               if( ofile)
                  {
                  FILE *ifile = fopen_ext( get_file_name( tbuff, elements_filename), "tfcrb");

                  while( fgets( tbuff, sizeof( tbuff), ifile))
                     fputs( tbuff, ofile);
                  fclose( ifile);
                  fclose( ofile);
                  }
               memcpy( orbit2, orbit, n_orbit_params * sizeof( double));
//             debug_printf( "Before integration (epoch %f):\n%20.14f %20.14f %20.14f\n",
//                               curr_epoch,
//                               orbit2[0], orbit2[1], orbit2[2]);
//             debug_printf( "%20.14f %20.14f %20.14f\n",
//                               orbit2[3], orbit2[4], orbit2[5]);
               integrate_orbit( orbit2, curr_epoch, epoch_shown);
//             debug_printf( "After  integration (%f):\n%20.14f %20.14f %20.14f\n",
//                               epoch_shown,
//                               orbit2[0], orbit2[1], orbit2[2]);
//             debug_printf( "%20.14f %20.14f %20.14f\n",
//                               orbit2[3], orbit2[4], orbit2[5]);
               }
            }
            break;
         case 't':
            residual_format = resid_format_menu( message_to_user, residual_format);
            if( *message_to_user)
               update_element_display = 1;
            break;
         case 'T':
            snprintf_err( message_to_user, sizeof( message_to_user),
                                       "%u observation(s) in tracklet",
                      select_tracklet( obs, n_obs, curr_obs));  /* orb_func.cpp */
            break;
         case KEY_CYCLE_RESID_DISPLAY_UP:
            if( residual_format & RESIDUAL_FORMAT_NORMALIZED)
                residual_format ^= RESIDUAL_FORMAT_NORMALIZED;
            else if( residual_format & RESIDUAL_FORMAT_TIME_RESIDS)
                residual_format ^= (RESIDUAL_FORMAT_NORMALIZED | RESIDUAL_FORMAT_TIME_RESIDS);
            else
                residual_format ^= RESIDUAL_FORMAT_TIME_RESIDS;
            update_element_display = 1;
            break;
         case KEY_CYCLE_RESID_DISPLAY_DN:
            if( residual_format & RESIDUAL_FORMAT_NORMALIZED)
                residual_format ^= (RESIDUAL_FORMAT_NORMALIZED | RESIDUAL_FORMAT_TIME_RESIDS);
            else if( residual_format & RESIDUAL_FORMAT_TIME_RESIDS)
                residual_format ^= RESIDUAL_FORMAT_TIME_RESIDS;
            else
                residual_format ^= RESIDUAL_FORMAT_NORMALIZED;
            update_element_display = 1;
            break;
         case 127:           /* backspace takes on different values */
         case KEY_BACKSPACE: /* on PDCurses,  ncurses,  etc.        */
         case 8:
            if( !pop_orbit( &curr_epoch, orbit))
               {
               strlcpy_error( message_to_user, "Last orbit operation undone");
               update_element_display = 1;
               set_locs( orbit, curr_epoch, obs, n_obs);
               }
            else
               strlcpy_error( message_to_user, "No more orbits to undo!");
            break;
         case 'V':         /* apply Vaisala method, without linearizing */
            if( !inquire( "Enter peri/apohelion distance: ", tbuff, sizeof( tbuff),
                        COLOR_DEFAULT_INQUIRY))
               {
               double vaisala_dist;
               extern int vaisala_center_object;
               const char *centers = "mveMjsunpl", *tptr = strchr( centers, *tbuff);

               if( tptr)
                  vaisala_center_object = (int)( tptr - centers) + 1;
               vaisala_dist = atof( tbuff + (vaisala_center_object != 0));
               if( vaisala_dist > 300.)   /* larger distances are assumed */
                  vaisala_dist /= AU_IN_KM;        /* to be in kilometers */
               if( vaisala_dist)
                  {
                  curr_epoch = obs->jd;
                  find_vaisala_orbit( orbit, obs, obs + n_obs - 1, vaisala_dist);
                  update_element_display = 1;
                  set_locs( orbit, curr_epoch, obs, n_obs);
                  }
               vaisala_center_object = 0;
               }
            break;
         case 'v':         /* apply Vaisala method */
            extended_orbit_fit( orbit, obs, n_obs, FIT_VAISALA_FULL, curr_epoch);
            update_element_display = 1;
            break;
         case ALT_H:
            curr_epoch = obs->jd;
            extended_orbit_fit( orbit, obs, n_obs,
                     atoi( get_environment_ptr( "H")), curr_epoch);
            update_element_display = 1;
            break;
#ifdef NOT_IN_USE
         case 'v':         /* apply Vaisala method */
            {
            double vaisala_dist, angle_param;
            int n_fields, success = 0;

            inquire( "Enter peri/apohelion distance: ", tbuff, sizeof( tbuff),
                                    COLOR_DEFAULT_INQUIRY);
            n_fields = sscanf( tbuff, "%lf,%lf", &vaisala_dist, &angle_param);
            push_orbit( curr_epoch, orbit);
            if( n_fields == 1)      /* simple Vaisala */
               {
               herget_method( obs, n_obs, -vaisala_dist, 0., orbit, NULL, NULL,
                                                   NULL);
               if( c != 'V')
                  adjust_herget_results( obs, n_obs, orbit);
               success = 1;
               }
            else if( n_fields == 2)
               {
               if( vaisala_dist > 6400.)     /* assume input in km, not AU */
                  vaisala_dist /= AU_IN_KM;
               if( angle_param == -1 || angle_param == 1)
                  {
                  set_distance( obs, vaisala_dist);
                  if( !find_parabolic_orbit( obs, n_obs, orbit,
                             angle_param > 0.))
                     success = 1;
                  }
               else if( angle_param == 2)
                  {
                  const int retval = search_for_trial_orbit( orbit, obs, n_obs,
                                              vaisala_dist, &angle_param);

                  if( retval)
                     snprintf( message_to_user, sizeof( message_to_user),
                                                "Trial orbit error %d\n", retval);
                  else
                     {
                     snprintf( message_to_user, sizeof( message_to_user),
                                                "Minimum at %f\n", angle_param);
                     success = 1;
                     }
                  }
               else                          /* "trial orbit" method */
                  {
                  const int retval = find_trial_orbit( orbit, obs, n_obs,
                                              vaisala_dist, angle_param);

                  if( retval)
                     snprintf_err( message_to_user, sizeof( message_to_user),
                                                "Trial orbit error %d\n", retval);
                  else
                     {
                     if( c != 'V')
                        adjust_herget_results( obs, n_obs, orbit);
                     success = 1;
                     }
                  }
               }
            if( success == 1)
               {
                        /* epoch is that of first included observation: */
               get_epoch_range_of_included_obs( obs, n_obs, &curr_epoch, NULL);
               get_r1_and_r2( n_obs, obs, &r1, &r2);
               update_element_display = 1;
               }
            }
            break;
#endif
         case 'w': case 'W':
            {
            double worst_rms = 0., rms, top_rms = 1e+20;

            if( prev_getch == 'w')
               {
               if( residual_format & RESIDUAL_FORMAT_NORMALIZED)
                  top_rms = compute_weighted_rms( obs + curr_obs, 1, NULL);
               else
                  top_rms = compute_rms( obs + curr_obs, 1);
               }

            for( i = 0; i < n_obs; i++)
                if( obs[i].is_included)
                   {
                   if( residual_format & RESIDUAL_FORMAT_NORMALIZED)
                      rms = compute_weighted_rms( obs + i, 1, NULL);
                   else
                      rms = compute_rms( obs + i, 1);
                   if( rms > worst_rms && rms < top_rms)
                      {
                      worst_rms = rms;
                      curr_obs = i;
                      }
                   }
            single_obs_selected = true;
            strlcpy_error( message_to_user, "Worst observation found");
            }
            break;
         case 'X':         /* toggle tracklet */
            for( i = 0; i < n_obs; i++)
               if( obs[i].flags & OBS_IS_SELECTED)
                  obs[i].flags |= OBS_TEMP_USE_FLAG;
               else
                  obs[i].flags &= ~OBS_TEMP_USE_FLAG;
            select_tracklet( obs, n_obs, curr_obs);  /* orb_func.cpp */
                     /* FALLTHRU */
         case 'x':
            {
            unsigned n_found;

            add_off_on = toggle_selected_observations( obs, n_obs, &n_found);
            snprintf_err( message_to_user, sizeof( message_to_user),
                                       "%u observation(s) toggled", n_found);
            if( c == 'X')     /* restore original selection */
               for( i = 0; i < n_obs; i++)
                  {
                  if( obs[i].flags & OBS_TEMP_USE_FLAG)
                     obs[i].flags |= OBS_IS_SELECTED;
                  else
                     obs[i].flags &= ~OBS_IS_SELECTED;
                  }
            }
            break;
         case 'y': case 'Y':
            show_a_file( "gauss.out", 0);
            break;
         case 'z': case 'Z':
            {
            double state2[MAX_N_PARAMS], delta_squared = 0;
            int64_t t0;
            const int64_t one_billion = (int64_t)1000000000;

            inquire( "Time span: ", tbuff, sizeof( tbuff),
                              COLOR_DEFAULT_INQUIRY);
            memcpy( state2, orbit, n_orbit_params * sizeof( double));
            t0 = nanoseconds_since_1970( );
            integrate_orbit( state2, curr_epoch, curr_epoch + atof( tbuff));
            integrate_orbit( state2, curr_epoch + atof( tbuff), curr_epoch);
            for( i = 0; i < 3; i++)
               {
               state2[i] -= orbit[i];
               delta_squared += state2[i] * state2[i];
               }
            t0 = nanoseconds_since_1970( ) - t0;
            snprintf_err( message_to_user, sizeof( message_to_user),
                              "Change = %.3e AU = %.3e km; %lu.%lu seconds",
                              sqrt( delta_squared),
                              sqrt( delta_squared) * AU_IN_KM,
                              (unsigned long)( t0 / one_billion),
                              (unsigned long)( t0 % one_billion));
            snprintf_append( message_to_user, sizeof( message_to_user),
                     " JD %.4f to %.4f", curr_epoch, curr_epoch + atof( tbuff));
            }
            break;
         case ALT_D:
            if( !inquire( "Debug level: ",
                                tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
               {
               char *equals_ptr = strchr( tbuff, '=');

               if( !equals_ptr)
                  debug_level = atoi( tbuff);
               else
                  {
                  *equals_ptr = '\0';
                  set_environment_ptr( tbuff, equals_ptr + 1);
                  }
               }
            break;
         case ALT_J:
            if( !inquire( "J2 multiplier: ",
                                tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
               {
               if( *tbuff == 'r')
                  {
                  unsigned j;

                  for( i = j = 0; i < n_obs; i++)
                     if( obs[i].note2 == 'R')
                        obs[j++] = obs[i];
                  n_obs = j;
                  update_element_display = 1;
                  }
               else
                  {
                  extern double j2_multiplier;

                  j2_multiplier = atof( tbuff);
                  }
               }
            break;
         case ALT_S:
            {
            extern bool use_symmetric_derivatives;

            use_symmetric_derivatives = !use_symmetric_derivatives;
            strlcpy_error( message_to_user, "Symmetric derivatives are");
            add_off_on = use_symmetric_derivatives;
            }
            break;
         case '$':
            {
            extern double integration_tolerance;

            snprintf_err( tbuff, sizeof( tbuff), "Integration tolerance: (currently %.4g)",
                                 integration_tolerance);
            if( !inquire( tbuff, tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY)
                                 && atof( tbuff) > 0.)
               integration_tolerance = atof( tbuff);
            }
            break;
         case '%':
            if( !inquire( "Uncertainty of selected observation(s), in arcsec: ",
                                tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
               {
               int n_selected = 0;
               double new_sig;
               bool is_mul = false;

               if( *tbuff == 't' || *tbuff == 'm')
                  new_sig = atof( tbuff + 1);
               else if( *tbuff == '*')
                  {
                  is_mul = true;
                  new_sig = atof( tbuff + 1);
                  }
               else
                  new_sig = atof( tbuff);
               for( i = 0; new_sig > 0. && i < n_obs; i++)
                  if( obs[i].flags & OBS_IS_SELECTED)
                     {
                     n_selected++;
                     switch( *tbuff)
                        {
                        case 'm':
                           obs[i].mag_sigma = new_sig;
                           snprintf_err( message_to_user, sizeof( message_to_user),
                                  "Magnitude sigma reset to %.3e", new_sig);
                           break;
                        case 't':
                           obs[i].time_sigma = new_sig / seconds_per_day;
                           snprintf_err( message_to_user, sizeof( message_to_user),
                                  "Time sigma reset to %.3e seconds", new_sig);
                           break;
                        default:
                           if( is_mul)
                              {
                              obs[i].posn_sigma_1 *= new_sig;
                              obs[i].posn_sigma_2 *= new_sig;
                              }
                           else
                              set_tholen_style_sigmas( obs + i, tbuff);
                           strlcpy_error( message_to_user, "Positional uncertainty reset");
                           break;
                        }
                     }
               snprintf_append( message_to_user, sizeof( message_to_user),
                        " for %d observations", n_selected);
               }
            break;
         case '"':
            debug_mouse_messages ^= 1;
            strlcpy_error( message_to_user, "Mouse debugging");
            add_off_on = debug_mouse_messages;
            break;
         case '@':
            {
            extern int setting_outside_of_arc;

            setting_outside_of_arc ^= 1;
            strlcpy_error( message_to_user, "Setting outside of arc turned");
            add_off_on = setting_outside_of_arc;
            }
            break;
         case '(':
            set_locs( orbit, curr_epoch, obs, n_obs);
            strlcpy_error( message_to_user, "Full arc set");
            break;
         case ')':
            if( !inquire( "Enter name of file to be displayed: ",
                               tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
               show_a_file( tbuff, 0);
            break;
         case '`':
            {
            default_comet_magnitude_type =
                        'N' + 'T' - default_comet_magnitude_type;
            if( default_comet_magnitude_type == 'N')
               strlcpy_error( message_to_user, "Using nuclear mags for comets");
            else
               strlcpy_error( message_to_user, "Using total mags for comets");
            }
         case KEY_MOUSE:   /* already handled above */
            break;
         case 'b':
            if( mouse_wheel_motion && mouse_x >= 5 && mouse_x < 17)
               {       /* clicked on YY MM DD.DDD */
               static const double time_diffs[] = { 3650., 365., 180., 90., 30.,
                              0., 10., 1., .3, 0.1, 0.01, 0.001 };
               const double step = time_diffs[mouse_x - 5];
               const double target_jd = obs[curr_obs].jd
                               + (double)mouse_wheel_motion * step;

               if( mouse_wheel_motion > 0)

                  while( curr_obs < n_obs - 1 && obs[curr_obs].jd < target_jd)
                     curr_obs++;
               else
                  while( curr_obs > 0 && obs[curr_obs].jd > target_jd)
                     curr_obs--;
               single_obs_selected = true;
               strlcpy_error( message_to_user, "Stepped through observations");
               }
            else
               {
               residual_format = _obs_format_menu( residual_format, true);
               strlcpy_error( message_to_user, "Obs time format reset");
               }
            break;
         case 27:
         case 'q': case 'Q':
#ifdef KEY_EXIT
         case KEY_EXIT:
#endif
            quit = 1;
            break;
         case '=':
            residual_format ^= RESIDUAL_FORMAT_MAG_RESIDS;
            strlcpy_error( message_to_user, "Magnitude residual display turned");
            add_off_on = (residual_format & RESIDUAL_FORMAT_MAG_RESIDS);
            break;
         case '+':
            select_central_object( tbuff, sizeof( tbuff), false);
            update_element_display = 1;
            break;
         case '[':
            show_a_file( get_file_name( tbuff, "covar.txt"), 0);
            break;
         case '-':
            {
            static const char *messages[3] = {
                           "One MPC code listed",
                           "Normal MPC code listing",
                           "Many MPC codes listed" };

            list_codes = (list_codes + 1) % 3;
            strlcpy_error( message_to_user, messages[list_codes]);
            }
            break;
         case ',':
            show_a_file( "debug.txt", 0);
            break;
         case KEY_F( 15):        /* shift-f3 : show Curses info */
            {
            snprintf_err( tbuff, sizeof( tbuff), "%s\n%s\n%s\n",
                                 longname( ), termname( ), curses_version( ));
            snprintf_append( tbuff, sizeof( tbuff),
                        "%d pairs of %d colors\n", COLOR_PAIRS, COLORS);
            if( can_change_color( ))
               strlcat_error( tbuff, "Colors are changeable\n");
            inquire( tbuff, NULL, 0, COLOR_DEFAULT_INQUIRY);
            }
            break;
#ifdef KEY_RESIZE
         case KEY_RESIZE:
            resize_term( 0, 0);
            snprintf_err( message_to_user, sizeof( message_to_user),
                        "KEY_RESIZE: %d x %d",
                        getmaxx( stdscr), getmaxy( stdscr));
            break;
#endif
         case '{':
            {
            const int using_normalized_resids =
                     (residual_format & RESIDUAL_FORMAT_NORMALIZED);
            const char *prompt = (using_normalized_resids ?
                     "Enter cutoff in sigmas :" : "Enter cutoff in arcseconds :");

            if( !inquire( prompt, tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY)
                    && atof( tbuff))
               {
               const double cutoff = atof( tbuff);

               if( cutoff > 0.)
                  {
                  const int n_changed = filter_obs( obs, n_obs, cutoff,
                                    !using_normalized_resids);

                  if( !n_changed)
                     strlcpy_error( message_to_user, "No filtering done");
                  else if( n_changed == -1)
                     strlcpy_error( message_to_user, "Filtering FAILED");
                  else
                     snprintf_err( message_to_user, sizeof( message_to_user),
                                "Rejections at %.3f %s; %d changes", cutoff,
                                 (using_normalized_resids ? "sigmas" : "arcsecs"),
                                 n_changed);
                  }
               }
            }
            break;
         case '}':
            {
            if( residual_format & RESIDUAL_FORMAT_OVERPRECISE)
               {
               residual_format ^= RESIDUAL_FORMAT_OVERPRECISE;
               strlcpy_error( message_to_user, "Normal resids");
               }
            else if( residual_format & RESIDUAL_FORMAT_PRECISE)
               {
               residual_format ^=
                      (RESIDUAL_FORMAT_OVERPRECISE ^ RESIDUAL_FORMAT_PRECISE);
               strlcpy_error( message_to_user, "Super-precise resids");
               }
            else
               {
               residual_format ^= RESIDUAL_FORMAT_PRECISE;
               strlcpy_error( message_to_user, "Precise resids");
               }
            add_off_on = 1;
            }
            break;
         case '_':
            link_arcs( obs, n_obs, r1, r2);
            show_a_file( "gauss.out", 0);
            break;
         case '>':
            perturbers = 1;
            break;
         case ALT_B:  /* set observed RA/dec = computed for selected obs */
         case ALT_C:  /* set observed RA/dec = computed for all obs */
            for( i = 0; i < n_obs; i++)
               if( c == ALT_C || (obs[i].flags & OBS_IS_SELECTED))
                  {
                  obs[i].ra  = obs[i].computed_ra;
                  obs[i].dec = obs[i].computed_dec;
                  }
            snprintf_err( message_to_user, sizeof( message_to_user),
                      "RA/decs for %s obs set to computed values",
                      (c == ALT_C ? "all" : "selected"));
            break;
#ifdef ALT_MINUS
         case ALT_MINUS:
            snprintf_err( message_to_user, sizeof( message_to_user), "Extended %d obs",
                   extend_orbit_solution( obs, n_obs, 100., 365.25 * 20.));
            break;
#endif
         case ALT_T:     /* remove timing errors for all obs */
         case ALT_V:     /* remove timing errors for selected obs */
            for( i = 0; i < n_obs; i++)
               if( c == ALT_T || (obs[i].flags & OBS_IS_SELECTED))
                  {
                  MOTION_DETAILS m;

                  compute_observation_motion_details( obs + i, &m);
                  obs[i].jd += m.time_residual / seconds_per_day;
                  set_up_observation( obs + i);         /* mpc_obs.cpp */
                  snprintf( message_to_user, sizeof( message_to_user),
                           "Timing errors zeroed for %s obs",
                           (c == ALT_T ? "all" : "selected"));
                  }
            break;
         case KEY_F(23):    /* shift-f11: ephemeride-less pseudo-MPEC */
            {
            extern const char *ephemeris_filename;
            extern const char *residual_filename;
            const char *path = get_environment_ptr( "SOHO_MPEC");
            double orbit2[MAX_N_PARAMS];

            memcpy( orbit2, orbit, n_orbit_params * sizeof( double));
            integrate_orbit( orbit2, curr_epoch, epoch_shown);
            create_obs_file( obs, n_obs, 0, residual_format);
#ifdef _WIN32                /* MS is different. */
            _unlink( get_file_name( tbuff, ephemeris_filename));
#else
            unlink( get_file_name( tbuff, ephemeris_filename));
#endif
            write_residuals_to_file( get_file_name( tbuff, residual_filename),
                             ifilename, n_obs, obs, RESIDUAL_FORMAT_SHORT);
            snprintf_err( tbuff, sizeof( tbuff), "%s%s.htm", path, obs->packed_id);
            debug_printf( "Creating '%s'\n", tbuff);
            text_search_and_replace( tbuff, " ", "");
            make_pseudo_mpec( tbuff, obj_name);
            strlcpy_error( message_to_user, "Ephemeride-less pseudo-MPEC made");
            }
            break;
         case KEY_F(22):    /* shift-f10 */
            snprintf_err( message_to_user, sizeof( message_to_user), "Euler = %f",
                           euler_function( obs, obs + n_obs - 1));
            break;
         case KEY_F(20):    /* shift-f8 */
            if( !inquire( "MPC code for which to search: ",
                               tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
               {
               if( strlen( tbuff) == 3)
                  {
                  curr_obs = (curr_obs + 1) % n_obs;
                  for( i = 0; i < n_obs && strcmp( tbuff,
                            obs[curr_obs].mpc_code); i++)
                     curr_obs = (curr_obs + 1) % n_obs;
                  }
               else if( !strcmp( tbuff, "r"))
                  {
                  curr_obs = (curr_obs + 1) % n_obs;
                  for( i = 0; i < n_obs && obs[curr_obs].note2 != 'R'; i++)
                     curr_obs = (curr_obs + 1) % n_obs;
                  }
               else
                  {
                  const double jd = get_time_from_string( obs[curr_obs].jd,
                                 tbuff, CALENDAR_JULIAN_GREGORIAN, NULL);

                  if( jd)
                     for( curr_obs = 0; curr_obs < n_obs &&
                            obs[curr_obs].jd < jd; curr_obs++)
                        ;
                  }
               single_obs_selected = true;
               }
            break;
         case KEY_F(19):    /* shift-f7 */
            element_format ^= ELEM_OUT_ALTERNATIVE_FORMAT;
            strlcpy_error( message_to_user, "Alternative element format");
            add_off_on = (element_format & ELEM_OUT_ALTERNATIVE_FORMAT);
            update_element_display = 1;
            break;
         case KEY_F(12):
            {
            int last, first;
            OBSERVE tobs1, tobs2;
            static int desired_soln = 0;

            first = find_first_and_last_obs_idx( obs, n_obs, &last);
            tobs1 = obs[first];
            tobs2 = obs[last];
            if( find_circular_orbits( &tobs1, &tobs2, orbit, desired_soln++))
               strlcpy_error( message_to_user, "No circular orbits found");
            else
               {
               curr_epoch = tobs1.jd - tobs1.r / AU_PER_DAY;
               set_locs( orbit, curr_epoch, obs, n_obs);
               get_r1_and_r2( n_obs, obs, &r1, &r2);
               update_element_display = 1;
               }
            }
            break;
         case KEY_F(21):    /* shift-f9 */
            {
            extern bool all_reasonable;

            all_reasonable = !all_reasonable;
            strlcpy_error( message_to_user,  "'All reasonable'");
            add_off_on = all_reasonable;
            }
            break;
         case 'u': case 'U':
            full_improvement( NULL, 0, NULL, 0., NULL, 0, 0.);
            break;
         case ALT_E:
            full_ctime( tbuff, curr_epoch, FULL_CTIME_YMD | FULL_CTIME_MILLIDAYS);
            snprintf_err( message_to_user, sizeof( message_to_user),
                        "Curr epoch %.40s = %.3f ", tbuff, curr_epoch);
            break;
         case ALT_M:
            if( !inquire( "Number Metropolis steps: ", tbuff, sizeof( tbuff),
                        COLOR_DEFAULT_INQUIRY) && (i = atoi( tbuff)) > 0)
               {
               metropolis_search( obs, n_obs, orbit, curr_epoch, i, 1.);
               update_element_display = 1;
               }
            break;
         case CTRL( 'B'):
            show_commented_elements = !show_commented_elements;
            break;
         case ALT_L:
            {
            const char *language_letters = "efirsd";

            c = inquire( get_find_orb_text( 2035), NULL, 0, COLOR_DEFAULT_INQUIRY);
            if( c >= KEY_F( 1) && c <= KEY_F( 6))
               c = language_letters[c - KEY_F( 1)];
            if( strchr( language_letters, c))
               {
               set_language( c);
               update_element_display = 1;
               }
            }
            break;
         case CTRL( 'K'):
            {
            extern int apply_debiasing;

            apply_debiasing ^= 1;
            strlcpy_error( message_to_user, "FCCT14 debiasing is");
            add_off_on = apply_debiasing;
            }
            break;
         case KEY_F(18):    /* shift-f6 */
            residual_format ^= RESIDUAL_FORMAT_COMPUTER_FRIENDLY;
            break;
         case '~':
            residual_format ^= RESIDUAL_FORMAT_EXTRA;
            break;
         case ALT_O:
            for( i = 0; i < n_obs; i++)
               if( !strcmp( obs[i].mpc_code, obs[curr_obs].mpc_code))
                  obs[i].flags |= OBS_IS_SELECTED;
               else
                  obs[i].flags &= ~OBS_IS_SELECTED;
            break;
         case ALT_N:
            select_element_frame( );
            update_element_display = 1;
            break;
         case ALT_G:
            orbital_monte_carlo( orbit, obs, n_obs, curr_epoch, epoch_shown);
            update_element_display = 1;
            strlcpy_error( message_to_user, "Orbital MC generated");
            break;
         case ALT_I:
            {
            make_config_dir_name( tbuff, elements_filename);
            i = copy_file_to_clipboard( tbuff);
            if( i)
               inquire( get_find_orb_text( 2039), NULL, 0,
                                    COLOR_MESSAGE_TO_USER);
            else
               strlcpy_error( message_to_user, "Elements copied to clipboard");
            }
            break;
         case ALT_K:
            {
            extern int sigmas_in_columns_57_to_65;

            strlcpy_error( message_to_user, "Sigma display");
            if( residual_format & RESIDUAL_FORMAT_SHOW_DESIGS)
               {
               residual_format ^= RESIDUAL_FORMAT_SHOW_DESIGS;
               sigmas_in_columns_57_to_65 = 0;
               }
            else
               {
               sigmas_in_columns_57_to_65 ^= 1;
               if( !sigmas_in_columns_57_to_65)
                  {
                  residual_format |= RESIDUAL_FORMAT_SHOW_DESIGS;
                  strlcpy_error( message_to_user, "Showing packed desigs in obs area");
                  }
               }
            if( message_to_user[1] == 'i')
               add_off_on = sigmas_in_columns_57_to_65;
            }
            break;
         case ALT_Q:
            {
            extern int show_observational_details;

            show_observational_details ^= 1;
            }
            break;
         case '&':
            residual_format = _obs_format_menu( residual_format, false);
            strlcpy_error( message_to_user, "Obs RA/dec format reset");
            break;
         case 9:
            sort_obs_by_code = !sort_obs_by_code;
            strlcpy_error( message_to_user, sort_obs_by_code ?
                     "Obs sorted by MPC code" : "Obs sorted by date");
            break;
         case ALT_R:
            if( !inquire( get_find_orb_text( 2038),
                               tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY))
               {
               set_environment_ptr( "REFERENCE", tbuff);
               update_element_display = 1;
               }
            break;
         case 'k':
            write_environment_pointers( );
            show_a_file( "env.txt", 0);
            break;
         case ALT_Z:
            {
            extern int integration_method;

            integration_method ^= 1;
            strlcpy_error( message_to_user, integration_method ?
                           "Using PD89" : "Using RKF");
            }
            break;
         case ALT_X:
            curr_obs = select_mpc_code( obs, n_obs, curr_obs);
            single_obs_selected = true;
            break;
         case KEY_ADD_MENU_LINE:
            n_command_lines++;
            strlcpy_error( message_to_user, "Adding a menu line");
            break;
         case KEY_REMOVE_MENU_LINE:
            n_command_lines--;
            strlcpy_error( message_to_user, "Removing a menu line");
            break;
         case ALT_P:    /* see 'environ.def' comments for COMET_CONSTRAINTS */
            {
            char prompt[200];

            snprintf( prompt, sizeof( prompt),
                        "Enter comet r0, m, n, k (current values are %s):",
                        get_environment_ptr( "COMET_CONSTANTS"));
            if( !inquire( prompt, tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY) && *tbuff)
               {
               set_environment_ptr( "COMET_CONSTANTS", tbuff);
               comet_g_func( 0.);
               }
            }
            break;
         case CTRL( 'X'):
            {
            saving_elements_for_reuse = true;
            strlcpy_error( message_to_user, "Elements saved for re-use");
            update_element_display = 1;
            }
            break;
         case CTRL( 'Y'):
            snprintf_err( message_to_user, sizeof( message_to_user),
                     "Exclusion file of %d obs written",
                                     write_excluded_observations_file( obs, n_obs));
            break;
         case CTL_UP:
            if( curr_obs)
               {
               curr_obs--;
               while( curr_obs && !opposition_break( obs + curr_obs - 1))
                  curr_obs--;
               single_obs_selected = true;
               }
            break;
#ifdef CTL_DOWN        /* PDCurses uses this #define */
         case CTL_DOWN:
#else                  /* ncurses uses this #define */
         case CTL_DN:
#endif
            if( curr_obs < n_obs - 1)
               {
               curr_obs++;
               while( curr_obs < n_obs - 1 && !opposition_break( obs + curr_obs))
                  curr_obs++;
               single_obs_selected = true;
               }
            break;
         case '\'':
            if( !inquire( "Enter override absolute mag (H) value:",
                       tbuff, sizeof( tbuff), COLOR_DEFAULT_INQUIRY) && *tbuff)
            {
            extern double override_abs_mag;

            override_abs_mag = atof( tbuff);
            }
            break;
         case 'i': case 'I':
            {
            extern int show_alt_info;

            show_alt_info ^= 1;
            strlcpy_error( message_to_user, "Display of 'alternative' observation info is");
            add_off_on = show_alt_info;
            }
            break;
         case 'c': case 'C':
            show_calendar( );
            break;
         case '.':
            show_splash_screen( );     /* just to test */
            do
               {
               c = extended_getch( );
               if( c == KEY_MOUSE)
                  {
                  get_mouse_data( (int *)&mouse_x, (int *)&mouse_y, (int *)&mouse_z, &button);
                  if( button & REPORT_MOUSE_POSITION)
                     c = 0;         /* ignore mouse moves */
                  }
               } while( !c);
            break;
         case ALT_Y:
            {
            char object_name[80];
            double unused_max_resid;
            ELEMENTS elems;

            get_object_name( object_name, obs->packed_id);
            if( get_orbit_from_mpcorb_sof( object_name, orbit,
                               &elems, 0.1, &unused_max_resid))
               {
               curr_epoch = epoch_shown = elems.epoch;
               update_element_display = 1;
               set_locs( orbit, curr_epoch, obs, n_obs);
               strlcpy_error( message_to_user, "Elements copied from 'mpcorb.sof'");
               }
            else
               snprintf_err( message_to_user, sizeof( message_to_user),
                       "Didn't find elems for '%s' in 'mpcorb.sof'", object_name);
            }
            break;
         case 'j': case 'J':
            {
            extern int object_type;
            const char *text;

            if( object_type == OBJECT_TYPE_ASTEROID)
               {
               object_type = OBJECT_TYPE_COMET;
               default_comet_magnitude_type = 'T';
               text = "comet total";
               }
            else if( default_comet_magnitude_type == 'T')
               {
               default_comet_magnitude_type = 'N';
               text = "comet nuclear";
               }
            else
               {
               object_type = OBJECT_TYPE_ASTEROID;
               text = "asteroid";
               }
            calc_absolute_magnitude( obs, n_obs);
            snprintf_err( message_to_user, sizeof( message_to_user),
                       "Using %s mag model", text);
            update_element_display = 1;
            }
            break;
         case '\\':
         case 'O':
         case ';': case ']':
         case CTRL( 'E'): case CTRL( 'J'): case CTRL( 'L'):
         case CTRL( 'N'): case CTRL( 'O'): case CTRL( 'Q'):
         case CTRL( 'S'): case CTRL( 'T'): case CTRL( 'U'):
         case CTRL( 'V'): case CTRL( 'W'): case CTRL( 'Z'):
         case CTRL( '_'): case CTRL( ']'):
         case ALT_A:
         case CTL_LEFT: case CTL_RIGHT:
         case KEY_F( 1):
         case KEY_F( 2):
         case KEY_F( 13):        /* shift-f1 */
         case KEY_F( 14):        /* shift-f2 */
         case KEY_F( 24):        /* shift-f12 */
         case KEY_DC:            /* delete key */
         case CTL_DEL:
         case ALT_DEL:
         case KEY_B2:            /* central key on numeric keypad */
         case KEY_ENTER:         /* on numeric keypad */
         case ALT_UP:
         case ALT_LEFT: case ALT_RIGHT:
         case KEY_SRIGHT: case KEY_SLEFT:
#ifdef __PDCURSES__
         case PADPLUS: case PADMINUS: case PADSLASH:
#endif
         case ALT_DOWN:
         case KEY_F(17):    /* shift-f5 */
         default:
            debug_printf( "Key %d hit\n", c);
            show_a_file( "dos_help.txt", 0);
            snprintf_err( message_to_user, sizeof( message_to_user),
                               "Key %d ($%x, o%o) '%s' hit",
                               c, c, c, keyname( c));
            break;
         }
      }
   attrset( COLOR_PAIR( COLOR_BACKGROUND));
   show_final_line( n_obs, curr_obs, COLOR_BACKGROUND);
Shutdown_program:
   full_endwin( );                 /* terminals to end mouse movement */
   delscreen( screen_ptr);  /* Not really needed,  but ensures Valgrind  */
   curses_running = false;               /* says all memory was freed */
   if( obs && n_obs)
      {
      create_obs_file( obs, n_obs, 0, residual_format);
      create_ades_file( "observe.xml",  obs, n_obs);
      }
   unload_observations( obs, n_obs);

   text_search_and_replace( mpc_code, " ", "_");
   snprintf_err( tbuff, sizeof( tbuff), "%s 0 %d %d",
                              mpc_code, residual_format, list_codes);
   set_environment_ptr( "CONSOLE_OPTS", tbuff);
   store_defaults( ephemeris_output_options, element_format,
         element_precision, max_residual_for_filtering,
         noise_in_sigmas);
   set_environment_ptr( "EPHEM_START", ephemeris_start);
   snprintf_err( tbuff, sizeof( tbuff), "%d", n_ephemeris_steps);
   set_environment_ptr( "EPHEM_STEPS", tbuff);
   set_environment_ptr( "EPHEM_STEP_SIZE", ephemeris_step_size);
   if( ids)
      free( ids);
   if( mpc_color_codes)
      free( mpc_color_codes);
   free( command_areas);
   free( key_remaps);
   clean_up_find_orb_memory( );
   return( 0);
}
