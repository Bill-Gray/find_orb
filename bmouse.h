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
