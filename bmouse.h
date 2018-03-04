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

#define BMOUSE struct bmouse

#ifdef __cplusplus
extern "C" {            /* Assume C declarations for C++ */
#endif

BMOUSE
   {
   unsigned xmin, xmax, ymin, ymax, x, y, sensitivity, prev_x, prev_y;
   unsigned pressed, released, accum_released, button;
   int dx, dy, diff_mode, mouse_installed;
   unsigned pressed_x, pressed_y;
   unsigned private_unsigned[2];
   int char_code;
   };

#define NO_MOUSE_FOUND -1

int init_mouse( BMOUSE *b);
int mouse_read( BMOUSE *b);
int mouse_set_position( BMOUSE *b, unsigned x, unsigned y);
#ifdef __cplusplus
}                       /* End of extern "C"  */
#endif   /* __cplusplus */
