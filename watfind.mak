# Make file for find_orb with the Watcom C/C++ compiler.
# Version using PDCurses.

PDC_DIR=..\PDCurs~1
PDC_LIB=$(PDC_DIR)\dosvga\PDCURSES.LIB
LIB_DIR=..\watlib
LINK_OPTS=l find_orb.lib l $(LIB_DIR)\WJPLEPH.LIB l $(LIB_DIR)\WAFUNCS.LIB l $(PDC_LIB) l $(LIB_DIR)\WSATLIB.LIB
# LINK_OPTS=l wafuncs.lib l wjpleph.lib l wsatlib.lib l find_orb.lib l pdcurses.lib f clipfunc.obj
CURSES_OBJ=clipfunc.obj getstrex.obj
CURSES_OBJ_COMMAS=clipfunc.obj, getstrex.obj

all: find_orb.exe fo.exe fo_serve.exe

LINKOPTS=option stub=dos32a option map=find_orb.map option stack=20000 f $(PDC_LIB)

OBJS=ades_out.obj b32_eph.obj bc405.obj bias.obj collide.obj conv_ele.obj &
  details.obj eigen.obj elem2tle.obj elem_out.obj elem_ou2.obj ephem0.obj &
  errors.obj expcalc.obj gauss.obj geo_pot.obj healpix.obj lsquare.obj &
  miscell.obj monte0.obj mpc_obs.obj &
!ifeq DOS Y
  $(CURSES_OBJ) &
!endif
  orb_func.obj orb_fun2.obj pl_cache.obj roots.obj &
  runge.obj simplex.obj sm_vsop.obj sr.obj shellsor.obj sigma.obj stackall.obj

find_orb.lib: $(OBJS)
   %write wfind.lrf $(OBJS)
   wlib -q -n -t $@ @wfind.lrf
   -rm wfind.lrf

find_orb.exe: findorb.obj find_orb.lib $(CURSES_OBJ)
   %write wfind.lrf $(LINK_OPTS)
   wlink name find_orb.exe f findorb.obj f $(CURSES_OBJ_COMMAS) @wfind.lrf
   -rm wfind.lrf

fo.exe: fo.obj find_orb.lib
   wlink N fo.exe f fo.obj $(LINK_OPTS)

fo_serve.exe: fo_serve.obj find_orb.lib cgi_func.obj
   wlink N fo_serve.exe f fo.obj f cgi_func.obj $(LINK_OPTS)

clean:
   -rm $(OBJS)
   -rm $(CURSES_OBJ)
   -rm cgi_func.obj
   -rm covar.txt
   -rm debug.txt
   -rm elements.txt
   -rm findorb.obj
   -rm find_orb.exe
   -rm find_orb.lib
   -rm find_orb.map
   -rm -v fo.exe
   -rm fo.obj
   -rm fo_serve.exe
   -rm fo_serve.obj
   -rm -v gauss.out
   -rm -v mpc_fmt.txt
   -rm -v ephemeri.txt guide.txt monte.txt mpec.htm residual.txt
   -rm -v mpcorb.dat
   -rm -v observe.txt
   -rm state.txt
   -rm wfind.lrf

CFLAGS=/Ox /W3 /4r /s /j /zq -i=..\include

.cpp.obj:
   wpp386 $(CFLAGS) $<

findorb.obj : findorb.cpp
   wpp386 $(CFLAGS) -i=$(PDC_DIR) $<

getstrex.obj : getstrex.cpp
   wpp386 $(CFLAGS) -i=$(PDC_DIR) $<

