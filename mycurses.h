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

#include <stdio.h>      /* Required for scr_dump() & scr_restore() below   */

#define attr_t       unsigned long

typedef struct _win      /* definition of a window.      */
{
   int xsize, ysize, xoffset, yoffset;
   int curr_x, curr_y, use_keypad, wrap_on;
   attr_t attr;
   attr_t data[1];
} WINDOW;

#ifndef STDSCR_DEFINED
extern WINDOW *stdscr, *curscr;
#endif

/* Shamelessly lifted from (PD)Curses.h: */
#ifndef FALSE
# define FALSE 0
#endif
#ifndef TRUE
# define TRUE 1
#endif
#ifndef NULL
# define NULL (void *)0
#endif
#ifndef ERR
# define ERR (-1)
#endif
#ifndef OK
# define OK 0
#endif

typedef unsigned short chtype;

#ifdef __cplusplus
extern "C" {
#endif
int wmove( WINDOW *w, const int y, const int x);
int waddch( WINDOW*w, const char c);
int wrefresh( const WINDOW *w);
int refresh( void);
int start_color( void );
int init_pair( const short idx, const short foreground, const short background);
int wattrset( WINDOW *w, const attr_t attr);
int wattroff( WINDOW*, const attr_t attr);
int wattron( WINDOW*, const attr_t attr);
int endwin( void );
WINDOW *initscr( void );
int noecho( void );
int echo( void );
int wgetnstr( WINDOW *w, char *str, int max_len);
int wclear( WINDOW *w);
int flushinp( void );
int waddnstr( WINDOW* w, const char *str, int len);
int resize_term( const int ysize, const int xsize);
int PDC_set_scrn_mode( const int new_mode);
WINDOW *resize_window( WINDOW *old, const int ysize, const int xsize);
int set_text_mode( const unsigned mode);
int get_text_mode_sizes( const int mode, int *ysize, int *xsize);
int wset_rect_attr( WINDOW *w, int xmin, int ymin, int xmax, int ymax);
int putwin( WINDOW *win, FILE *filep);
WINDOW *getwin( FILE *filep);
int scr_dump( const char *filename);
int scr_restore( const char *filename);
char *longname( void);
int curs_set( const int visibility);
chtype inch(void);
chtype mvinch(int, int);
chtype mvwinch(WINDOW *, int, int);
chtype winch(WINDOW *);
int waddchnstr(WINDOW *, const chtype *, int);
bool can_change_color( void);

#define keypad(w,flag)          (w->use_keypad = flag)
#define attroff(attr)           wattroff( stdscr, attr )
#define attron(attr)            wattron( stdscr, attr )
#define attrset(attr)           wattrset( stdscr, attr )
#define clear()                 wclear( stdscr )
#define getnstr(str,num)        wgetnstr( stdscr, str, num )
#define getmaxx(w)              (w)->xsize
#define getmaxy(w)              (w)->ysize
#define getmaxyx(w,y,x)         ( y = (w)->ysize, x = (w)->xsize )
#define move(y,x)               wmove( stdscr, y, x )
#define addch( c )              waddch( stdscr, c )
#define mvaddch(y,x,c)          (move( y, x )==ERR?ERR:addch( c ))
#define mvwaddch(w,y,x,c)       (wmove( w, y, x )==ERR?ERR:waddch( w, c ))
#define mvaddchstr(y,x,c)       (move( y, x )==ERR?ERR:addchnstr( c, -1 ))
#define mvaddchnstr(y,x,c,n)    (move( y, x )==ERR?ERR:addchnstr( c, n ))
#define mvaddnstr(y,x,str,n)    (move( y, x )==ERR?ERR:addnstr( str, n ))
#define addnstr(str, n)         waddnstr( stdscr, str, n )
#define addstr(str)             addnstr( str, -1)
#define mvaddstr(y,x,str)       mvaddnstr(y,x,str,-1)
#define waddstr(w, str)         waddnstr( w, str, -1)
#define cbreak()
#define addchnstr( str, n)      waddchnstr( stdscr, str, n)
int     napms(int);

#ifndef STDSCR_DEFINED
extern attr_t curses_attrs[];
#endif

#define COLOR_PAIRS        256
#define COLOR_PAIR( X)          (curses_attrs[X])
#define A_BLINK          0x8000
#define A_STANDOUT       0x0800
#define A_BOLD           0x0800
#define A_CHARTEXT       0x00ff

#define KEY_DOWN       ( 80 + 256)
#define KEY_UP         ( 72 + 256)
#define KEY_LEFT       ( 75 + 256)
#define KEY_RIGHT      ( 77 + 256)
#define KEY_C3         ( 81 + 512)
#define KEY_NPAGE      ( 81 + 256)
#define KEY_A3         ( 73 + 512)
#define KEY_PPAGE      ( 73 + 256)
#define KEY_C1         ( 79 + 512)
#define KEY_END        ( 79 + 256)
#define KEY_A1         ( 71 + 512)
#define KEY_HOME       ( 71 + 256)
#define KEY_IC         ( 82 + 256)
      /* KEY_F( 1...12) = F1...12    */
      /* KEY_F(13...24) = Shift-F1...12    */
      /* KEY_F(25...36) = Ctrl-F1...12    */
      /* KEY_F(37...48) = Alt-F1...12    */

#define KEY_F1_TO_10(n)  (( 58 + 256)+(n)-(((n)-1)/12)*2+((n)>12 ? 15 : 0))
#define KEY_F11_OR_12(n) ((122 + 256)+(n)-(((n)-1)/12)*10)
#define KEY_F(n)     (((n)-1) % 12 >= 10 ? KEY_F11_OR_12(n) : KEY_F1_TO_10(n))
#define ALT_A     ( 30 + 256)
#define ALT_B     ( 48 + 256)
#define ALT_C     ( 46 + 256)
#define ALT_D     ( 32 + 256)
#define ALT_E     ( 18 + 256)
#define ALT_F     ( 33 + 256)
#define ALT_G     ( 34 + 256)
#define ALT_H     ( 35 + 256)
#define ALT_I     ( 23 + 256)
#define ALT_J     ( 36 + 256)
#define ALT_K     ( 37 + 256)
#define ALT_L     ( 38 + 256)
#define ALT_M     ( 50 + 256)
#define ALT_N     ( 49 + 256)
#define ALT_O     ( 24 + 256)
#define ALT_P     ( 25 + 256)
#define ALT_Q     ( 16 + 256)
#define ALT_R     ( 19 + 256)
#define ALT_S     ( 31 + 256)
#define ALT_T     ( 20 + 256)
#define ALT_U     ( 22 + 256)
#define ALT_V     ( 47 + 256)
#define ALT_W     ( 17 + 256)
#define ALT_X     ( 45 + 256)
#define ALT_Y     ( 21 + 256)
#define ALT_Z     ( 44 + 256)
#define ALT_1     (120 + 256)
#define ALT_2     (121 + 256)
#define ALT_3     (122 + 256)
#define ALT_4     (123 + 256)
#define ALT_5     (124 + 256)
#define ALT_6     (125 + 256)
#define ALT_7     (126 + 256)
#define ALT_8     (127 + 256)
#define ALT_9     (128 + 256)
#define ALT_0     (129 + 256)

#define CTL_PAD0        (146 + 256)          /* ctl-keypad 0 */
#define CTL_PAD1        (117 + 256)          /* ctl-keypad 1 */
#define CTL_PAD2        (145 + 256)          /* ctl-keypad 2 */
#define CTL_PAD3        (118 + 256)          /* ctl-keypad 3 */
#define CTL_PAD4        (115 + 256)          /* ctl-keypad 4 */
#define CTL_PAD5        (143 + 256)          /* ctl-keypad 5 */
#define CTL_PAD6        (116 + 256)          /* ctl-keypad 6 */
#define CTL_PAD7        (119 + 256)          /* ctl-keypad 7 */
#define CTL_PAD8        (141 + 256)          /* ctl-keypad 8 */
#define CTL_PAD9        (132 + 256)          /* ctl-keypad 9 */

#define CTL_LEFT        (115 + 256)   /* Control-Left-Arrow   PC only  */
#define CTL_RIGHT       (116 + 256)   /* Control-Right-Arrow  PC only  */
#define CTL_UP          401     /* ctl-up arrow                  */
#define CTL_DOWN        397     /* ctl-down arrow                */
#define CTL_INS         402     /* ctl-insert                    */
#define ALT_DEL         419     /* alt-delete                    */
#define ALT_INS         418     /* alt-insert                    */
#define CTL_PGUP        (132 + 256)   /* Control-PgUp         PC only  */
#define CTL_PGDN        (118 + 256)   /* Control-PgDn         PC only  */
#define CTL_HOME        (119 + 256)   /* Control-Home         PC only  */
#define CTL_END         (117 + 256)   /* Control-End          PC only  */
#define CTL_PADCENTER   399     /* ctl-enter on keypad           */
#define CTL_PADPLUS     400     /* ctl-plus on keypad            */
#define CTL_PADMINUS    398     /* ctl-minus on keypad           */
#define CTL_PADSLASH    26      /* ctl-slash on keypad... but also Ctrl-Z */
#define CTL_PADSTAR     406     /* ctl-star on keypad            */
#define ALT_PADPLUS     334     /* alt-plus on keypad            */
#define ALT_PADMINUS    330     /* alt-minus on keypad           */
#define ALT_PADSLASH    300     /* alt-slash on keypad           */
#define ALT_PADENTER    422     /* alt-enter on keypad           */
#define ALT_LBRACKET    282     /* alt-left bracket              */
#define ALT_RBRACKET    283     /* alt-right bracket             */
#define ALT_MINUS       386     /* alt-minus                     */
#define ALT_EQUAL       387     /* alt-equal                     */
#define ALT_HOME        407     /* alt-home                      */
#define ALT_PGUP        409     /* alt-pgup                      */
#define ALT_PGDN        417     /* alt-pgdn                      */
#define ALT_END         415     /* alt-end                       */
#define ALT_UP          408     /* alt-up arrow                  */
#define ALT_DOWN        416     /* alt-down arrow                */
#define ALT_RIGHT       411     /* alt-right arrow               */
#define ALT_LEFT        413     /* alt-left arrow                */

#define ZKEY_CTRL_A     1
#define ZKEY_CTRL_B     2
#define ZKEY_CTRL_C     3
#define ZKEY_CTRL_D     4
#define ZKEY_CTRL_E     5
#define ZKEY_CTRL_F     6
#define ZKEY_CTRL_G     7
#define ZKEY_CTRL_H     8
#define ZKEY_CTRL_I     9
#define ZKEY_CTRL_J    10
#define ZKEY_CTRL_K    11
#define ZKEY_CTRL_L    12
#define ZKEY_CTRL_M    13
#define ZKEY_CTRL_N    14
#define ZKEY_CTRL_O    15
#define ZKEY_CTRL_P    16
#define ZKEY_CTRL_Q    17
#define ZKEY_CTRL_R    18
#define ZKEY_CTRL_S    19
#define ZKEY_CTRL_T    20
#define ZKEY_CTRL_U    21
#define ZKEY_CTRL_V    22
#define ZKEY_CTRL_W    23
#define ZKEY_CTRL_X    24
#define ZKEY_CTRL_Y    25
#define ZKEY_CTRL_Z    26
#define KEY_BTAB        ( 15 + 256)
#define ALT_COMMA       ( 51 + 256)
#define ALT_STOP        ( 52 + 256)
#define CTL_DEL         (147 + 256)
#define KEY_DC          ( 83 + 256)
#define PADENTER    13

#define KEY_MOUSE       (255 + 256)

#define KEY_RESIZE      (256 + 256)

#define COLOR_BLACK      0
#define COLOR_BLUE       1
#define COLOR_GREEN      2
#define COLOR_CYAN       3
#define COLOR_RED        4
#define COLOR_MAGENTA    5
#define COLOR_YELLOW     6
#define COLOR_WHITE      7

#define COLORS 16

#define TEXT_MODE_80x25          0
#define TEXT_MODE_80x28          1
#define TEXT_MODE_80x50          2
#define TEXT_MODE_100x37         3
#define TEXT_MODE_80x60          4
#define TEXT_MODE_132x44         5
#define TEXT_MODE_132x50         6
#define TEXT_MODE_132x25         7
#define TEXT_MODE_132x28         8

/* with my current video card,  only the first three are supported: */
#define N_TEXT_MODES             3


   /* None of the following work;  they're all dummy functions */
   /* and macros supplied to avoid compiler warnings/errors.   */
typedef unsigned long mmask_t;

int init_color(short, short, short, short);
mmask_t mousemask(mmask_t, mmask_t *);

#define BUTTON_CONTROL      0x0010  /* PDCurses */
#define ALL_MOUSE_EVENTS        0x1fffffffL

#ifdef __cplusplus
}
#endif
