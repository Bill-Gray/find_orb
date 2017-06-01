# Make file for eph2tle with g++
# (eph2tle reads an ephemeris of state vectors for a geocentric object and
# computes TLEs fitted to them.  See 'eph2tle.cpp' for details)

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

ifdef XCOMPILE
	CC=x86_64-w64-mingw32-g++
	EXE=.exe
	LIBSADDED=-static-libgcc -L $(INSTALL_DIR)/win_lib -lm -mwindows
endif

all: eph2tle$(EXE)

OBJS=eph2tle.o elem2tle.o lsquare.o conv_ele.o

eph2tle$(EXE): $(OBJS)
	$(CC) -o eph2tle$(EXE) $(OBJS) -llunar -lsatell $(LIBSADDED)

CFLAGS=-O3 -Wall -pedantic -Wextra -c -I $(INSTALL_DIR)/include

.cpp.o:
	$(CC) $(CFLAGS) $<

clean:
	rm -f eph2tle$(EXE) $(OBJS)
