# Make file for console Find_Orb,  using regular Curses
# Use 'bsdmake' for BSD
# GNU MAKE Makefile for Find_Orb
#
# Usage: make -f [path/]linmake [CLANG=Y] [XCOMPILE=Y] [MSWIN=Y] [X=Y] [tgt]
#
#	where tgt can be any of:
# [all|find_orb|fo|fo_serve|clean|clean_temp|eph2tle|cssfield]
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
LIBSADDED=-L $(INSTALL_DIR)/lib -lm
EXE=
OBJSADDED=
RM=rm -f

ifdef CLANG
	CC=clang
endif

# I'm using 'mkdir -p' to avoid error messages if the directory exists.
# It may fail on very old systems,  and will probably fail on non-POSIX
# systems.  If so,  change to '-mkdir' and ignore errors.

ifdef MSWIN
	LIBSADDED=-static-libgcc
	EXE=.exe
	OBJSADDED=clipfunc.o
	CURSES_LIB=-lpdcurses
	MKDIR=-mkdir
else
	MKDIR=mkdir -p
endif

# You can have your include files in ~/include and libraries in
# ~/lib,  in which case only the current user can use them;  or
# (with root privileges) you can install them to /usr/local/include
# and /usr/local/lib for all to enjoy.

ifdef GLOBAL
	INSTALL_DIR=/usr/local
else
	INSTALL_DIR=~
endif

ifdef X
	ADDED_CFLAGS=-DXCURSES -DPDC_WIDE -I../PDCurses
	CURSES_LIB=-lXCurses -lXaw -lXmu -lXt -lX11 -lSM -lICE -lXext -lXpm
endif

ifdef XCOMPILE
	CC=x86_64-w64-mingw32-g++
	ADDED_CFLAGS=-DUTF8 -DPDC_WIDE -I $(INSTALL_DIR)/include -I../PDCurses
 OBJSADDED=clipfunc.o
	EXE=.exe
	CURSES_LIB=-lpdcurses -static-libgcc
	LIBSADDED=-L $(INSTALL_DIR)/win_lib -lm -lgdi32 -luser32 -mwindows
endif

all: fo$(EXE) find_orb$(EXE) fo_serve.cgi eph2tle$(EXE)

CFLAGS=-c -O3 -Wall -pedantic -Wextra -Wno-unused-parameter -I $(INSTALL_DIR)/include

OBJS=b32_eph.o bc405.o bias.o collide.o conv_ele.o details.o eigen.o \
	elem2tle.o elem_out.o elem_ou2.o ephem0.o errors.o gauss.o   \
	geo_pot.o healpix.o lsquare.o miscell.o moid4.o monte0.o \
	mpc_obs.o mt64.o nanosecs.o orb_func.o orb_fun2.o pl_cache.o roots.o  \
	runge.o shellsor.o sigma.o simplex.o sm_vsop.o sr.o stackall.o $(OBJSADDED)

LIBS=$(LIBSADDED) -llunar -ljpl -lsatell

find_orb$(EXE):          findorb.o $(OBJS)
	$(CC) -o find_orb$(EXE) findorb.o $(OBJS)  $(LIBS) $(CURSES_LIB)

fo$(EXE):          fo.o $(OBJS)
	$(CC) -o fo$(EXE) fo.o $(OBJS) $(LIBS)

eph2tle$(EXE):          eph2tle.o conv_ele.o elem2tle.o simplex.o lsquare.o
	$(CC) -o eph2tle$(EXE) eph2tle.o conv_ele.o elem2tle.o simplex.o lsquare.o $(LIBS)

cssfield$(EXE):          cssfield.o
	$(CC) -o cssfield$(EXE) cssfield.o $(LIBS)

fo_serve.cgi:          fo_serve.o $(OBJS)
	$(CC) -o fo_serve.cgi fo_serve.o $(OBJS) $(LIBS)

clean:
	$(RM) $(OBJS) fo.o findorb.o fo_serve.o find_orb$(EXE) fo$(EXE)
	$(RM) fo_serve.cgi eph2tle.o eph2tle$(EXE) cssfield$(EXE)
	$(RM) cssfield.o

install:
	-cp find_orb $(HOME)/bin
	cp fo       $(HOME)/bin
	$(MKDIR) $(IDIR)
	cp command.txt details.txt dosephem.txt dos_help.txt ?findorb.txt           $(IDIR)
	cp environ.def geo_rect.txt header.htm jpl_eph.txt mpcorb.hdr               $(IDIR)
	cp mu1.txt observer.txt obslinks.htm ObsCodes.htm ObsCodesF.html            $(IDIR)
	cp  odd_name.txt rovers.txt scopes.txt sigma.txt xdesig.txt                 $(IDIR)

uninstall:
	rm -f $(HOME)/bin/find_orb
	rm -f $(HOME)/bin/fo
	rm -f $(IDIR)/*
	rmdir $(IDIR)

.cpp.o:
	$(CC) $(CFLAGS) $(ADDED_CFLAGS) $<
