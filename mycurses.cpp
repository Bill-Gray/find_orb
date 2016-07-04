/* mycurses.cpp: really basic DOS Curses

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

#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <i86.h>
#define STDSCR_DEFINED    1
#include "mycurses.h"

static int scr_xsize, scr_ysize, echo_on = 0;
static int initial_cursor_x, initial_cursor_y;
static char *backup_screen;
WINDOW *stdscr;

/* At present,  I don't see any actual use for curscr in this setup.
Everything goes directly to the screen.  If,  at some point,  I attempt to
make some more full-fledged curses clone,  I may have to reconsider.  */
#ifdef DEFINE_CURSCR
WINDOW *curscr;
#endif

static int get_cursor_loc( int *row, int *col);
static int set_cursor_loc( const int row, const int col);
static int set_cursor_mode( const int startrow, const int endrow );

#define DISPMEM ((unsigned char *)0xb8000)

WINDOW *resize_window( WINDOW *old, const int ysize, const int xsize)
{
   WINDOW *rval = (WINDOW *)calloc( 1,
                      sizeof( WINDOW) + xsize * ysize * sizeof( attr_t));

   if( rval)
      {
      if( old)
         memcpy( rval, old, sizeof( WINDOW));
      rval->xsize = xsize;
      rval->ysize = ysize;
      rval->curr_x = rval->curr_y = 0;
      }
   if( old)
      free( old);
   return( rval);
}

/* Once upon a time,  I had a graphics card that supported _nine_ different
text modes.  My current,  theoretically more advanced,  card supports only
80-wide with 25, 28,  or 50 lines.  They call this 'progress'... anyway,
I've included all the data for the nine possible modes. */

#define TOTAL_N_TEXT_MODES       9

static int default_mode_index = 0;
static const unsigned char text_mode_xsizes[TOTAL_N_TEXT_MODES] =
                 { 80, 80, 80, 100, 80, 132, 132, 132, 132 };
static const unsigned char text_mode_ysizes[TOTAL_N_TEXT_MODES] =
                 { 25, 28, 50, 37, 60, 44, 50, 25, 28 };


static void init_stdscr( void)
{
   const char *default_mode_text = getenv( "MYCURSES_MODE");

   if( default_mode_text)
      default_mode_index = (int)( *default_mode_text - '0');
   resize_term( text_mode_ysizes[default_mode_index],
                text_mode_xsizes[default_mode_index]);
   backup_screen = (char *)malloc( scr_xsize * scr_ysize * 2);
   memcpy( backup_screen, DISPMEM, scr_xsize * scr_ysize * 2);
   stdscr->wrap_on = 1;
}

WINDOW *initscr( void )
{
   get_cursor_loc( &initial_cursor_y, &initial_cursor_x);
   init_stdscr( );
   return( stdscr);
}

int wmove( WINDOW *w, const int y, const int x)
{
   int rval = ERR;

   if( x >= 0 && x < w->xsize && y >= 0 && y < w->ysize)
      {
      w->curr_x = x;
      w->curr_y = y;
      rval = 0;
      }
   return( rval);
}

int waddnstr( WINDOW* w, const char *str, int len)
{
   attr_t *dataptr = w->data + (w->curr_x + w->curr_y * w->xsize);
   int max_out;

   if( w->curr_y >= w->ysize)
      return( -1);
   max_out = w->xsize - w->curr_x;
   if( w->wrap_on)
      max_out += (w->ysize - w->curr_y - 1) * w->xsize;
   if( len < 0)
      len = strlen( str);
   if( len > max_out)
      len = max_out;
   w->curr_x += len;
   if( w->wrap_on)
      while( w->curr_x > w->xsize)
         {
         w->curr_x -= w->xsize;
         w->curr_y++;
         }
   while( *str && len--)
      *dataptr++ = (attr_t)((unsigned char)*str++) | w->attr;
   return( 0);
}

int waddch( WINDOW *w, const char c)
{
   return( waddnstr( w, &c, 1));
}

int start_color( void )
{
   return( 0);
}

int wclear( WINDOW *w)
{
   int i = w->xsize * w->ysize;
   attr_t *tptr = w->data;

   while( i--)
      *tptr++ = 0xf20;
   return( 0);
}

attr_t curses_attrs[COLOR_PAIRS];

int init_pair( const short idx, const short foreground, const short background)
{
   if( idx >= 0 && idx < COLOR_PAIRS)
      curses_attrs[idx] = ((attr_t)foreground << 8)
                        | ((attr_t)background << 12);
   return( 0);
}

int wattrset( WINDOW *w, const attr_t attr)
{
   w->attr = attr;
   return( 0);
}

int wattroff( WINDOW *w, const attr_t attr)
{
   w->attr = (attr_t)( w->attr & ~attr);
   return( 0);
}

int wattron( WINDOW *w, const attr_t attr)
{
   w->attr = (attr_t)( w->attr | attr);
   return( 0);
}

int wrefresh_rect( const WINDOW *w, int xmin, int ymin, int xmax, int ymax)
{
   int y;

   if( w->yoffset < 0)
      ymin = -w->yoffset;
   if( ymax > scr_ysize - w->yoffset)
      ymax = scr_ysize - w->yoffset;
   if( w->xoffset < 0)
      xmin = -w->xoffset;
   if( xmax > scr_xsize - w->xoffset)
      xmax = scr_xsize - w->xoffset;

   if( xmin < xmax)
      for( y = ymin; y < ymax; y++)
         {
         const attr_t *tptr = w->data + y * w->xsize + xmin;
         unsigned char *dispmem = DISPMEM
                   + ((y + w->yoffset) * scr_xsize + (xmin + w->xoffset)) * 2;
         const attr_t *endptr = tptr + xmax - xmin;

         while( tptr != endptr)
            {
            *dispmem++ = (unsigned char)*tptr;
            *dispmem++ = (unsigned char)(*tptr++ >> 8);
            }
         }
   return( 0);
}

int wset_rect_attr( WINDOW *w, int xmin, int ymin, int xmax, int ymax)
{
   int y;

   if( w->yoffset < 0)
      ymin = -w->yoffset;
   if( ymax > scr_ysize - w->yoffset)
      ymax = scr_ysize - w->yoffset;
   if( w->xoffset < 0)
      xmin = -w->xoffset;
   if( xmax > scr_xsize - w->xoffset)
      xmax = scr_xsize - w->xoffset;

   if( xmin < xmax)
      for( y = ymin; y < ymax; y++)
         {
         attr_t *tptr = w->data + y * w->xsize + xmin;
         attr_t *endptr = tptr + (xmax - xmin);

         while( tptr != endptr)
            {
            *tptr = (w->attr | (*tptr & 0xff));
            tptr++;
            }
         }
   return( 0);
}

int wrefresh( const WINDOW *w)
{
   return( wrefresh_rect( w, 0, 0, w->xsize, w->ysize));
}

int refresh( void)
{
   if( !stdscr)
      init_stdscr( );
   set_cursor_loc( stdscr->curr_y + stdscr->yoffset,
                   stdscr->curr_x + stdscr->xoffset);
   return( wrefresh( stdscr));
}

int noecho( void )
{
   echo_on = 0;
   return( 0);
}

int echo( void )
{
   echo_on = 1;
   return( 1);
}

static int set_cursor_mode( const int startrow, const int endrow )
{
   union _REGS regset;

   regset.h.ah = 0x01;
   regset.h.ch = (unsigned char) startrow;
   regset.h.cl = (unsigned char) endrow;
   int386(0x10, &regset, &regset);
   return( 0 );
}

static int set_cursor_loc( const int row, const int col)
{
   union _REGS regset;

   regset.h.ah = 0x02;
   regset.h.bh = 0;
// regset.h.bh = SP->video_page;
   regset.h.dh = (unsigned char) row;
   regset.h.dl = (unsigned char) col;
   int386( 0x10, &regset, &regset);
   return( 0);
}

static int get_cursor_loc( int *row, int *col)
{
   union _REGS regset;

   regset.h.ah = 0x03;
   regset.h.bh = 0;
// regset.h.bh = SP->video_page;
   int386(0x10, &regset, &regset);
   *row = regset.h.dh;
   *col = regset.h.dl;
   return( 0);
}

int wgetnstr( WINDOW *w, char *str, int max_len)
{
   int c = 0;
   int len = 0, insert_mode = 0;

   set_cursor_mode( 14, 16);
   while( len < max_len && c != 13)
      {
      int x = w->curr_x, y = w->curr_y;

      set_cursor_loc( y + w->yoffset, x + w->xoffset);
      c = getch( );
      if( !c)
         c = getch( ) + 256;
      if( c == ( 82 + 256))         /* Insert key */
         {
         insert_mode ^= 1;
         set_cursor_mode( 14 - insert_mode * 8, 16);
         }
      else if( c != 13)
         {
         if( c != 8)
            {
            str[len++] = (char)c;
            waddch( w, (char)c);
            }
         else if( len > 0)
            {
            len--;
            str[len] = '\0';
            x--;
            mvwaddch( w, y, x, ' ');
            wmove( w, y, x);
            }
         if( echo_on)
            wrefresh_rect( w, x, y, x + 1, y + 1);
         }
      }
   str[len] = '\0';
   if( insert_mode)
      set_cursor_mode( 14, 16);
   return( 0);
}

int flushinp( void )
{
   while( kbhit( ))
      getch( );
   return( 0);
}

int resize_term( const int ysize, const int xsize)
{
   scr_xsize = xsize;
   scr_ysize = ysize;
   stdscr = resize_window( stdscr, ysize, xsize);
#ifdef DEFINE_CURSCR
   curscr = resize_window( curscr, ysize, xsize);
#endif
   return( 0);
}

int PDC_set_scrn_mode( const int new_mode)
{
   return( 0);
}

int endwin( void )
{
// set_text_mode( TEXT_MODE_80x25);
   set_text_mode( default_mode_index);
   set_cursor_loc( initial_cursor_y, initial_cursor_x);
   if( stdscr)
      {
      free( stdscr);
      stdscr = NULL;
      }
#ifdef DEFINE_CURSCR
   if( curscr)
      {
      free( curscr);
      curscr = NULL;
      }
#endif
   if( backup_screen)
      {
      memcpy( DISPMEM, backup_screen, scr_xsize * scr_ysize * 2);
      free( backup_screen);
      }
   return( 0);
}

int set_text_mode( const unsigned mode)
{
   static const unsigned char modes[TOTAL_N_TEXT_MODES] =
                 { 3, 3, 3, 42, 38, 34, 34, 34, 34 };
            /* for a 16-high font, use 0x14 */
            /* for a 14-high font, use 0x11 */
            /* for an 8-high font, use 0x12 */
   static const unsigned char font_sizes[TOTAL_N_TEXT_MODES] =
                         { 0x14, 0x11, 0x12, 0, 0, 0, 0x12, 0x14, 0x11 };
   union _REGS regset;

   regset.w.ax = modes[mode];
   int386( 0x10, &regset, &regset);
   if( font_sizes[mode])
      {
      regset.h.ah = 0x11;
      regset.h.al = font_sizes[mode];
      regset.h.bl = 0x00;
      int386( 0x10, &regset, &regset);
      }
   resize_term( text_mode_ysizes[mode], text_mode_xsizes[mode]);
   return( 0);
}

int get_text_mode_sizes( const int mode, int *ysize, int *xsize)
{
   int rval = 0;

   if( mode < 0 || mode >= N_TEXT_MODES)
      rval = -1;
   else
      {
      if( xsize)
         *xsize = (int)text_mode_xsizes[mode];
      if( ysize)
         *ysize = (int)text_mode_ysizes[mode];
      }
   return( rval);
}

int putwin(WINDOW *win, FILE *filep)
{
   fwrite( win, sizeof( WINDOW), 1, filep);
   fwrite( win->data, sizeof( attr_t), win->xsize * win->ysize, filep);
   return( OK);
}

static int read_window( WINDOW *win, FILE *filep)
{
   fread( win, sizeof( WINDOW), 1, filep);
   fread( win->data, sizeof( attr_t), win->xsize * win->ysize, filep);
   return( OK);
}

WINDOW *getwin(FILE *filep)
{
   WINDOW temp, *rval;

   fread( &temp, sizeof( WINDOW), 1, filep);
   rval = (WINDOW *)calloc(
            sizeof( WINDOW) + temp.xsize * temp.ysize * sizeof( attr_t), 1);
   memcpy( rval, &temp, sizeof( WINDOW));
   fread( rval->data, sizeof( attr_t), temp.xsize * temp.ysize, filep);
   return( rval);
}

int scr_dump(const char *filename)
{
   int rval = ERR;

   if( stdscr)
      {
      FILE *ofile = fopen( filename, "wb");

      if( ofile)
         {
         rval = putwin( stdscr, ofile);
         fclose( ofile);
         }
      }
   return( rval);
}

int scr_restore(const char *filename)
{
   int rval = ERR;

   if( stdscr)
      {
      FILE *ifile = fopen( filename, "rb");

      if( ifile)
         {
         rval = read_window( stdscr, ifile);
         fclose( ifile);
         }
      }
   return( rval);
}

char *longname( void)
{
   return( "MyCurses");
}

            /* Cannibalized,  with little modification,  from PDCurses. */

int curs_set( const int visibility)
{
    union _REGS regs;
    int start, end;

    switch (visibility)
    {
        case 0:  /* invisible */
            start = 32;
            end = 0;  /* was 32 */
            break;
        case 2:  /* highly visible */
            start = 0;   /* full-height block */
            end = 7;
            break;
        default:  /*' normal' visibility */
            start = 0;
            end = 31;
            break;
    }

    /* if scrnmode is not set, some BIOSes hang */

    regs.h.ah = 0x01;
    regs.h.al = (unsigned char)3;
    regs.h.ch = (unsigned char)start;
    regs.h.cl = (unsigned char)end;
    int386(0x10, &regs, &regs);

    return 1;
}

      /* *inch functions basically copied from PDCurses: */

chtype winch(WINDOW *win)
{
    if( !win)
        return (chtype)ERR;

    return win->data[win->curr_x + win->curr_y * win->xsize];
}

chtype inch(void)
{
    return winch(stdscr);
}

chtype mvinch(int y, int x)
{
    if (move(y, x) == ERR)
        return (chtype)ERR;

    return( winch( stdscr));
}

chtype mvwinch(WINDOW *win, int y, int x)
{
    if (wmove(win, y, x) == ERR)
        return (chtype)ERR;

    return( winch( win));
}

int waddchnstr(WINDOW *win, const chtype *ch, int n)
{
    attr_t *tptr = win->data + win->curr_x + win->curr_y * win->xsize;

    if (n == -1 || n > win->xsize - win->curr_x)
        n = win->xsize - win->curr_x;
    if( n < 0)
        return( ERR);
    while( n--)
        *tptr++ = (attr_t)*ch++;
    return( OK);
}

#ifdef __WATCOMC__
#define ZKEY_CLOCK  (*(long *)((long)0x46c))
#else
#define ZKEY_CLOCK  (*(long far *)((long)0x46c))
#endif

int napms( int n_milliseconds)
{
   const long t_end = ZKEY_CLOCK + (long)n_milliseconds / 55;
   long curr_t = 0;

   while( curr_t < t_end)
      curr_t = ZKEY_CLOCK;
   return( 0);
}

int init_color(short color, short red, short green, short blue)
{
   return( 0);
}

mmask_t mousemask(mmask_t ignored1, mmask_t *ignored2)
{
   return( (mmask_t)0);
}
