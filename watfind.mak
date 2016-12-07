# Make file for find_orb with the Watcom C/C++ compiler.
# Version for my machine, using 'mycurses'.  For pdcurses,  switch to
# the other CFLAGS definition,  and link to pdcurses.lib instead of
# mycurses.obj.

!ifdef PDCURSES
LINK_OPTS=l pdcurses.lib f clipfunc.obj
CURSES_OBJ=clipfunc.obj
!else
LINK_OPTS=f mycurses.obj, bmouse.obj
CURSES_OBJ=mycurses.obj bmouse.obj
!endif

all: find_orb.exe

LINKOPTS=option stub=dos32a option map=find_orb.map option stack=20000

OBJS=findorb.obj b32_eph.obj bc405.obj collide.obj conv_ele.obj &
  eigen.obj elem2tle.obj elem_out.obj ephem0.obj gauss.obj &
  jpleph.obj lsquare.obj moid4.obj monte0.obj mpc_obs.obj &
  mt64.obj orb_func.obj orb_fun2.obj pl_cache.obj  &
  roots.obj runge.obj sm_vsop.obj sr.obj tle_out.obj sigma.obj $(CURSES_OBJ)

find_orb.exe: $(OBJS)
   wlink N find_orb.exe @wat_find.lnk l wafuncs.lib $(LINK_OPTS)

clean:
   rm $(OBJS)
   rm #(CURSES_OBJ) findorb.exe

CFLAGS=/Ox /W3 /4r /s /j /zq

.cpp.obj:
   wcc386 $(CFLAGS) $<

