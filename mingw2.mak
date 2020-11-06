# Make file for console Find_Orb,  using regular Curses
# Use 'xlinmake' for version with XCurses

all: find_orb fo

CC=g++

CFLAGS=-c -O3 -Wall

OBJS=b32_eph.o bc405.o clipfunc.o collide.o conv_ele.o eigen.o \
	elem2tle.o elem_out.o ephem0.o gauss.o \
	jpleph.o lsquare.o moid4.o monte0.o mpc_obs.o \
	orb_func.o orb_fun2.o pl_cache.o roots.o runge.o \
	sm_vsop.o sr.o tle_out.o weight.o

find_orb:          findorb.o $(OBJS) lunar.a
	$(CC) -o find_orb findorb.o $(OBJS) lunar.a c:\pdcurses\win32a\pdcurses.a

fo:          fo.o $(OBJS) lunar.a
	$(CC) -o fo fo.o $(OBJS) lunar.a

.cpp.o:
	$(CC) $(CFLAGS) $<
