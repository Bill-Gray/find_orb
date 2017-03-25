# Make file for vec2tle with g++
# (vec2tle reads an ephemeris of state vectors for a geocentric object and
# computes TLEs fitted to them.  See 'vec2tle.cpp' for details)

ifdef CLANG
CC=clang
ADDED_LIBS=-lm
else
CC=g++
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

all: vec2tle

OBJS=vec2tle.o elem2tle.o lsquare.o conv_ele.o

vec2tle: $(OBJS)
	$(CC) -o vec2tle $(OBJS) -llunar -lsatell $(ADDED_LIBS) -L $(INSTALL_DIR)/lib

CFLAGS=-O3 -Wall -c -I $(INSTALL_DIR)/include $<

.cpp.o:
	$(CC) $(CFLAGS)

clean:
	rm -f vec2tle $(OBJS)
