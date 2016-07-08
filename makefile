# Make file for console Find_Orb,  using regular Curses
# Use 'bsdmake' for BSD
# GNU MAKE Makefile for Find_Orb
#
# Usage: make -f [path/]linmake [CLANG=Y] [XCOMPILE=Y] [MSWIN=Y] [X=Y] [tgt]
#
#	where tgt can be any of:
# [all|find_orb|fo|fo_serve|clean|clean_temp]
#
#	'XCOMPILE' = cross-compile for Windows,  using MinGW,  on a Linux box
#	'MSWIN' = compile for Windows,  using MinGW and PDCurses,  on a Windows machine
#	'CLANG' = use clang instead of GCC;  Linux only
# 'X' = use PDCurses instead of ncurses
# None of these: compile using g++ on Linux,  for Linux
#
# 'clean_temp' removes various temporary files made by Find_Orb and friends:
# orbital elements,  covariance matrices,  etc.  'clean' removes these _and_
# the more usual object and executable files.

CURSES_LIB=-lncursesw
CC=g++
LIBSADDED=-lm
EXE=
OBJSADDED=
RM=rm -f

ifdef CLANG
	CC=clang
endif

ifdef MSWIN
	LIBSADDED=
	EXE=.exe
	OBJSADDED=clipfunc.o
	CURSES_LIB=pdcurses.a -static-libgcc
endif

ifdef X
	ADDED_CFLAGS=-DXCURSES -DPDC_WIDE -I../PDCurses
	CURSES_LIB=-lXCurses -lXaw -lXmu -lXt -lX11 -lSM -lICE -lXext -lXpm
endif

ifdef XCOMPILE
	CC=x86_64-w64-mingw32-g++
	ADDED_CFLAGS=-DUTF8 -DPDC_WIDE -I/usr/local/include
 OBJSADDED=clipfunc.o
	EXE=.exe
	LIBSADDED=
	CURSES_LIB=pdcurses.a -static-libgcc
endif

all: fo$(EXE) find_orb$(EXE) fo_serve.cgi

CFLAGS=-c -O3 -Wall -pedantic -Wextra -Wno-unused-parameter

OBJS=b32_eph.o bc405.o bias.o collide.o conv_ele.o eigen.o \
	elem2tle.o elem_out.o ephem0.o gauss.o geo_pot.o healpix.o \
	lsquare.o moid4.o monte0.o mpc_obs.o mt64.o \
	orb_func.o orb_fun2.o pl_cache.o roots.o  \
	runge.o sigma.o sm_vsop.o sr.o $(OBJSADDED)

LIBS=-llunar -ljpl -lsatell

find_orb$(EXE):          findorb.o $(OBJS)
	$(CC) -o find_orb$(EXE) findorb.o $(OBJS) $(CURSES_LIB) $(LIBSADDED) $(LIBS)

fo$(EXE):          fo.o $(OBJS)
	$(CC) -o fo$(EXE) fo.o $(OBJS) $(LIBSADDED) $(LIBS)

fo_serve.cgi:          fo_serve.o cgi_func.o $(OBJS)
	$(CC) -o fo_serve.cgi fo_serve.o cgi_func.o $(OBJS) $(LIBSADDED) $(LIBS)

clean:
	$(RM) $(OBJS) fo.o findorb.o fo_serve.o find_orb$(EXE) fo$(EXE)
	$(RM) fo_serve.cgi cgi_func.o
	$(RM) covar.txt covar?.txt debug.txt eleme?.txt elements.txt
	$(RM) ephemeri.txt gauss.out guide.txt guide?.txt monte.txt monte?.txt
	$(RM) mpc_f?.txt mpc_fmt.txt mpc_s?.txt mpec.htm obser?.txt observe.txt
	$(RM) residual.txt sr_el?.txt state.txt state?.txt virtu?.txt virtual.txt
	$(RM) sr_elems.txt mpcorb.dat

clean_temp:
	$(RM) covar.txt covar?.txt debug.txt eleme?.txt elements.txt
	$(RM) ephemeri.txt gauss.out guide.txt guide?.txt monte.txt monte?.txt
	$(RM) mpc_f?.txt mpc_fmt.txt mpc_s?.txt mpec.htm obser?.txt observe.txt
	$(RM) residual.txt sr_el?.txt state.txt state?.txt virtu?.txt virtual.txt
	$(RM) sr_elems.txt mpcorb.dat

.cpp.o:
	$(CC) $(CFLAGS) $(ADDED_CFLAGS) $<
