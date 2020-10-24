# Make file for console Find_Orb,  using regular Curses
# Use 'bsdmake' for BSD
# GNU MAKE Makefile for Find_Orb
#
# Usage: make -f [path/]linmake [CLANG=Y] [W32=Y] [W64=Y] [MSWIN=Y] [X=Y] [VT=Y] [tgt]
#
#	where tgt can be any of:
# [all|find_orb|fo|fo_serve|clean|clean_temp|eph2tle|cssfield|neat_xvt]
#
#	'W32'/'W64' = cross-compile for 32- or 64-bit Windows,  using MinGW-w64,
#      on a Linux box
#	'MSWIN' = compile for Windows,  using MinGW and PDCurses,  on a Windows machine
#	'CLANG' = use clang instead of GCC;  Linux only
# 'X' = use PDCurses instead of ncurses
# 'VT' = use PDCurses with VT platform (see github.com/Bill-Gray/PDCurses/vt)
# 'CC=g++-4.8' = use that version of g++;  helpful when testing older compilers
# None of these: compile using g++ on Linux,  for Linux
#
# 'clean_temp' removes various temporary files made by Find_Orb and friends:
# orbital elements,  covariance matrices,  etc.  'clean' removes these _and_
# the more usual object and executable files.

CURSES_LIB=-lncursesw
CC=g++
LIBSADDED=-L $(INSTALL_DIR)/lib -lm
EXE=
RM=rm -f

ifeq ($(OS),Windows_NT)
    detected_OS := Windows
else
    detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

ifeq ($(detected_OS),Darwin)  # Mac OS X uses 'ncurses', not 'ncursesw'
	CURSES_LIB=-lncurses
endif

ifdef CLANG
	CC=clang
endif

# I'm using 'mkdir -p' to avoid error messages if the directory exists.
# It may fail on very old systems,  and will probably fail on non-POSIX
# systems.  If so,  change to '-mkdir' and ignore errors.

ifdef MSWIN
	LIBSADDED=-static-libgcc
	EXE=.exe
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
	CURSES_FLAGS=-DXCURSES -I../PDCurses
	CURSES_LIB=-lXCurses -lXaw -lXmu -lXt -lX11 -lSM -lICE -lXext -lXpm
endif

ifdef VT
	CURSES_FLAGS=-DVT -I$(HOME)/PDCurses
	CURSES_LIB=-lpdcurses
endif

LIB_DIR=$(INSTALL_DIR)/lib

ifdef W64
	CC=x86_64-w64-mingw32-g++
	CURSES_FLAGS=-I $(INSTALL_DIR)/include -I../PDCurses
	EXE=.exe
	CURSES_LIB=-lpdcurses
	LIB_DIR=$(INSTALL_DIR)/win_lib
	LIBSADDED=-L $(LIB_DIR) -lm -lgdi32 -luser32 -mwindows -static-libgcc
	FO_EXE=fo64.exe
	FIND_ORB_EXE=find_o64.exe
endif

ifdef W32
	CC=i686-w64-mingw32-g++
	CURSES_FLAGS=-I $(INSTALL_DIR)/include -I../PDCurses
	EXE=.exe
	CURSES_LIB=-lpdcurses
	LIB_DIR=$(INSTALL_DIR)/win_lib32
	LIBSADDED=-L $(LIB_DIR) -lm -lgdi32 -luser32 -mwindows -static-libgcc
	FO_EXE=fo32.exe
	FIND_ORB_EXE=find_o32.exe
endif

ifndef FO_EXE
	FO_EXE=fo
	FIND_ORB_EXE=find_orb
endif

all: $(FO_EXE) $(FIND_ORB_EXE) fo_serve.cgi eph2tle$(EXE)

CFLAGS=-c -Wall -pedantic -Wextra -I $(INSTALL_DIR)/include

ifdef DEBUG
	CFLAGS += -g -O0
else
	CFLAGS += -O3
endif

OBJS=ades_out.o b32_eph.o bc405.o bias.o collide.o conv_ele.o details.o eigen.o \
	elem2tle.o elem_out.o elem_ou2.o ephem0.o errors.o expcalc.o gauss.o \
	geo_pot.o healpix.o lsquare.o miscell.o         monte0.o \
	mpc_obs.o nanosecs.o orb_func.o orb_fun2.o pl_cache.o roots.o  \
	runge.o shellsor.o sigma.o simplex.o sm_vsop.o sr.o stackall.o

LIBS=$(LIBSADDED) -llunar -ljpl -lsatell

$(FIND_ORB_EXE):          findorb.o clipfunc.o $(OBJS)
	$(CC) -o $(FIND_ORB_EXE) findorb.o clipfunc.o $(OBJS) $(LIBS) $(CURSES_LIB)

findorb.o:         findorb.cpp
	$(CC) $(CFLAGS) $(CURSES_FLAGS) $<

clipfunc.o:        clipfunc.cpp
	$(CC) $(CFLAGS) $(CURSES_FLAGS) $<

$(FO_EXE):          fo.o $(OBJS)
	$(CC) -o $(FO_EXE) fo.o $(OBJS) $(LIBS)

eph2tle$(EXE):          eph2tle.o conv_ele.o elem2tle.o simplex.o lsquare.o
	$(CC) -o eph2tle$(EXE) eph2tle.o conv_ele.o elem2tle.o simplex.o lsquare.o $(LIBS)

cssfield$(EXE):          cssfield.o
	$(CC) -o cssfield$(EXE) cssfield.o $(LIBS)

expcalc$(EXE):          expcalc.cpp
	$(CC) -o expcalc$(EXE) -Wall -Wextra -pedantic -DTEST_CODE expcalc.cpp

roottest$(EXE):          roottest.o
	$(CC) -o roottest$(EXE) roottest.o roots.o

neat_xvt$(EXE):          neat_xvt.o
	$(CC) -o neat_xvt$(EXE) neat_xvt.o

fo_serve.cgi:          fo_serve.o $(OBJS)
	$(CC) -o fo_serve.cgi fo_serve.o $(OBJS) $(LIBS)

IDIR=$(HOME)/.find_orb

clean:
	$(RM) $(OBJS) fo.o findorb.o fo_serve.o $(FIND_ORB_EXE) $(FO_EXE)
	$(RM) fo_serve.cgi eph2tle.o eph2tle$(EXE) cssfield$(EXE)
	$(RM) clipfunc.o cssfield.o neat_xvt.o neat_xvt$(EXE)

clean_temp:
	$(RM) $(IDIR)/bc405pre.txt
	$(RM) $(IDIR)/cmt_sof.txt
	$(RM) $(IDIR)/covar.txt
	$(RM) $(IDIR)/covar?.txt
	$(RM) $(IDIR)/debug.txt
	$(RM) $(IDIR)/eleme?.txt
	$(RM) $(IDIR)/elements.txt
	$(RM) $(IDIR)/ephemeri.txt
	$(RM) $(IDIR)/gauss.out
	$(RM) $(IDIR)/guide.txt
	$(RM) $(IDIR)/guide?.txt
	$(RM) $(IDIR)/monte.txt
	$(RM) $(IDIR)/monte?.txt
	$(RM) $(IDIR)/mpcorb.dat
	$(RM) $(IDIR)/mpc_f?.txt
	$(RM) $(IDIR)/mpc_fmt.txt
	$(RM) $(IDIR)/mpc_s?.txt
	$(RM) $(IDIR)/mpec.htm
	$(RM) $(IDIR)/obser?.txt
	$(RM) $(IDIR)/observe.txt
	$(RM) $(IDIR)/residual.txt
	$(RM) $(IDIR)/sof.txt
	$(RM) $(IDIR)/sof?.txt
	$(RM) $(IDIR)/sofv?.txt
	$(RM) $(IDIR)/sr_el?.txt
	$(RM) $(IDIR)/sr_elems.txt
	$(RM) $(IDIR)/state.txt
	$(RM) $(IDIR)/state?.txt
	$(RM) $(IDIR)/vectors.txt
	$(RM) $(IDIR)/virtu?.txt
	$(RM) $(IDIR)/virtual.txt

install:
	$(MKDIR) $(IDIR)
ifdef EXE
	cp $(FIND_ORB_EXE) $(IDIR)
	cp $(FO_EXE)       $(IDIR)
else
	cp $(FIND_ORB_EXE) $(HOME)/bin
	cp $(FO_EXE)       $(HOME)/bin
endif
	cp command.txt details.txt dosephem.txt dos_help.txt ?findorb.txt           $(IDIR)
	cp environ.def frame_he.txt geo_rect.txt header.htm jpl_eph.txt             $(IDIR)
	cp link_def.json mpc_area.txt mpcorb.hdr                                    $(IDIR)
	cp mu1.txt observer.txt obslinks.htm ObsCodes.htm ObsCodesF.html            $(IDIR)
	cp obj_help.txt odd_name.txt openfile.txt previous.def rovers.txt           $(IDIR)
	cp sat_xref.txt scope.json scopes.txt sigma.txt                             $(IDIR)
	cp xdesig.txt bright2.pgm elem_pop.txt                                      $(IDIR)

uninstall:
ifdef EXE
	rm -f $(IDIR)/$(FIND_ORB_EXE)
	rm -f $(IDIR)/$(FO_EXE)
else
	rm -f $(HOME)/bin/$(FIND_ORB_EXE)
	rm -f $(HOME)/bin/$(FO_EXE)
endif
	rm -f $(IDIR)/*
	rmdir $(IDIR)

.cpp.o:
	$(CC) $(CFLAGS) $<
