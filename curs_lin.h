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

#define ALT_BASE 0x1000

#define ALT_A  (ALT_BASE + 'a')
#define ALT_B  (ALT_BASE + 'b')
#define ALT_C  (ALT_BASE + 'c')
#define ALT_D  (ALT_BASE + 'd')
#define ALT_E  (ALT_BASE + 'e')
#define ALT_F  (ALT_BASE + 'f')
#define ALT_G  (ALT_BASE + 'g')
#define ALT_H  (ALT_BASE + 'h')
#define ALT_I  (ALT_BASE + 'i')
#define ALT_J  (ALT_BASE + 'j')
#define ALT_K  (ALT_BASE + 'k')
#define ALT_L  (ALT_BASE + 'l')
#define ALT_M  (ALT_BASE + 'm')
#define ALT_N  (ALT_BASE + 'n')
#define ALT_O  (ALT_BASE + 'o')
#define ALT_P  (ALT_BASE + 'p')
#define ALT_Q  (ALT_BASE + 'q')
#define ALT_R  (ALT_BASE + 'r')
#define ALT_S  (ALT_BASE + 's')
#define ALT_T  (ALT_BASE + 't')
#define ALT_U  (ALT_BASE + 'u')
#define ALT_V  (ALT_BASE + 'v')
#define ALT_W  (ALT_BASE + 'w')
#define ALT_X  (ALT_BASE + 'x')
#define ALT_Y  (ALT_BASE + 'y')
#define ALT_Z  (ALT_BASE + 'z')
#define ALT_0  (ALT_BASE + '0')
#define ALT_1  (ALT_BASE + '1')
#define ALT_2  (ALT_BASE + '2')
#define ALT_3  (ALT_BASE + '3')
#define ALT_4  (ALT_BASE + '4')
#define ALT_5  (ALT_BASE + '5')
#define ALT_6  (ALT_BASE + '6')
#define ALT_7  (ALT_BASE + '7')
#define ALT_8  (ALT_BASE + '8')
#define ALT_9  (ALT_BASE + '9')

#define ALT_COMMA     (ALT_BASE + ',')
#define ALT_STOP      (ALT_BASE + '.')

#ifdef OLD_VALUES_SEE_COMMENT_BELOW
#define CTL_PGUP      0x22a
#define CTL_PGDN      0x225
#define CTL_RIGHT     0x22f
#define CTL_LEFT      0x220
#define CTL_UP        0x235
#define CTL_DN        0x20c
#define CTL_DEL       0x206
#define ALT_PGUP      0x228
#define ALT_PGDN      0x223
#define ALT_RIGHT     0x22d
#define ALT_LEFT      0x21e
#define ALT_UP        0x233
#define ALT_DOWN      0x20a
// #define CTL_LEFT      (ALT_BASE + 257)
// #define CTL_RIGHT     (ALT_BASE + 258)
// #define CTL_PGUP      (ALT_BASE + 259)
// #define CTL_PGDN      (ALT_BASE + 260)
// #define CTL_DEL       (ALT_BASE + 261)

#define CTL_PAD2     0x020c
#define CTL_PAD3     0x0225
#define CTL_PAD4     0x0220
#define CTL_PAD6     0x022f
#define CTL_PAD8     0x0235
#define CTL_PAD9     0x022a
      /*  CTL_PAD1, 5, 7 all map to 0x105B and can't be distinguished */
#endif

/* 2016 Jul 31:  after upgrading from Xubuntu 14.04 to 16.04,  some key
code changed.  CTL_PGUP and CTL_PGDN no longer are recognized by ncurses.
The following values got bumped up by one.  I don't have a good solution
to this at present;  it appears,  at least,  as if ncurses does not have a
good way of describing which of these keys has been hit,  and one has to
have version-specific code.  My regard for the developers of ncurses is
such that I doubt they did something that dumb;  this is probably
something I've just not figured out yet.      */

#define CTL_RIGHT     0x230
#define CTL_LEFT      0x221
#define CTL_UP        0x236
#define CTL_DN        0x20d
#define CTL_DEL       0x207
#define ALT_PGUP      0x229
#define ALT_PGDN      0x224
#define ALT_RIGHT     0x22e
#define ALT_LEFT      0x21f
#define ALT_UP        0x234
#define ALT_DOWN      0x20b

#define ALT_DEL       0x205
#define ALT_INS       0x21a

#define CTL_PAD2     0x020d
#define CTL_PAD3     0x0226
#define CTL_PAD4     0x0221
#define CTL_PAD6     0x0230
#define CTL_PAD8     0x0236
#define CTL_PAD9     0x022b
