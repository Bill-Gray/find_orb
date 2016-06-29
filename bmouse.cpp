#include <stdio.h>
#include <stdlib.h>
#include <dos.h>
#include "bmouse.h"

#define MOUSEPORT 51
#define leftshift 0

static int mouse_call( int command, int ms2, int ms3)
{
   union REGS in,out;

#ifdef __WATCOMC__
   in.w.ax = (unsigned short)command;
   in.w.bx = (unsigned)0;
   in.w.cx = (unsigned short)ms2;
   in.w.dx = (unsigned short)ms3;
   in.w.cflag = 0;

#ifdef __386__
   int386( MOUSEPORT, &in, &out);
#else
   int86( MOUSEPORT, &in, &out);
#endif
   return( out.w.ax);
#else
         /* non-Watcom */
   in.x.ax = (unsigned)command;
   in.x.bx = (unsigned)0;
   in.x.cx = (unsigned)ms2;
   in.x.dx = (unsigned)ms3;
   in.x.cflag = 0;

   int86( MOUSEPORT, &in, &out);
   return( out.x.ax);
#endif
}

int init_mouse( BMOUSE *b)
{
   if( b->sensitivity == 0)
      b->sensitivity = 3;

   b->pressed = b->released = b->button = b->accum_released = 0;
   b->dx = b->dy = b->private_unsigned[0] = b->private_unsigned[1] = 0;
   b->pressed_x = b->pressed_y = 0;
   b->prev_x = b->x;
   b->prev_y = b->y;
   if( b->diff_mode)
      {
      b->xmin = b->ymin = 0;
      b->xmax = b->ymax = 1000;
      b->x = b->y = 500;
      }

   if( b->mouse_installed == 2)
      return( 0);
   if( !mouse_call( 0, 0, 0))      /* install & reset the mouse */
      return( NO_MOUSE_FOUND);

   b->mouse_installed = 1;
   mouse_call( 2, 0, 0);   /* get rid of the cursor */
// mouse_call( 1, 0, 0);   /* show mouse cursor */

   mouse_call( 7, b->xmin << leftshift, b->xmax << leftshift);     /* set the x-range */
   mouse_call( 8, b->ymin << leftshift, b->ymax << leftshift);     /* set the y-range */
   mouse_call( 15, b->sensitivity, b->sensitivity);
   mouse_call( 4, b->x << leftshift, b->y << leftshift);           /* set starting loc */
   return( 0);
}

void mouse_update_buttons( BMOUSE *b, int button)
{
   b->pressed = button & (~b->button);
   b->released = (~button) & b->button;
   if( b->button && !button)     /* all buttons have been released */
      {
      b->accum_released = b->private_unsigned[0];
      b->private_unsigned[0] = 0;
      }
   else
      {
      b->accum_released = 0;
      if( button && !b->private_unsigned[0])
         {
         b->pressed_x = b->x;
         b->pressed_y = b->y;
         }
      b->private_unsigned[0] |= button;
      }
   b->button = button;
}

int mouse_read( BMOUSE *b)
{
   int button;
   unsigned x, y;
   union REGS in,out;

   if( !b->mouse_installed)
      return( -1);
#ifdef __WATCOMC__
   in.w.ax = 3;   /* get mouse location */
   in.w.cflag=0;
#ifdef __386__
   int386( MOUSEPORT, &in, &out);
#else
   int86( MOUSEPORT, &in, &out);
#endif
   button = (int)out.w.bx;
   x = out.w.cx >> leftshift;
   y = out.w.dx >> leftshift;
#else
         /* non-Watcom */
   in.x.ax = 3;   /* get mouse location */
   in.x.cflag=0;
   int86( MOUSEPORT, &in, &out);
   button = (int)out.x.bx;
   x = out.x.cx >> leftshift;
   y = out.x.dx >> leftshift;
#endif
   b->dx = (int)( x - b->x);
   b->dy = (int)( y - b->y);
   b->prev_x = b->x;
   b->prev_y = b->y;
   if( b->diff_mode && (b->dx || b->dy))
      mouse_set_position( b, 500, 500);
   else
      {
      b->x = x;
      b->y = y;
      }
   mouse_update_buttons( b, button);
   return( 0);
}

int mouse_set_position( BMOUSE *b, unsigned x, unsigned y)
{
   if( x > b->xmax) x = b->xmax;
   if( y > b->ymax) y = b->ymax;
   if( x < b->xmin) x = b->xmin;
   if( y < b->ymin) y = b->ymin;
   b->prev_x = b->x;
   b->prev_y = b->y;
   b->x = x;
   b->y = y;
   if( !b->mouse_installed)
      return( -1);
   if( b->mouse_installed == 1)
      mouse_call( 4, x << leftshift, y << leftshift);
   return( 0);
}
