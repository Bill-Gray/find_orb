# Make file for find_orb with the Watcom C/C++ compiler.
# Version using PDCurses.

LINK_OPTS=l pdcurses.lib f clipfunc.obj
CURSES_OBJ=clipfunc.obj

all: find_orb.exe fo.exe fo_serve.exe

LINKOPTS=option stub=dos32a option map=find_orb.map option stack=20000

OBJS=b32_eph.obj bc405.obj bias.obj collide.obj conv_ele.obj &
  eigen.obj elem2tle.obj elem_out.obj elem_ou2.obj ephem0.obj &
  gauss.obj geo_pot.obj healpix.obj jpleph.obj lsquare.obj &
  miscell.obj moid4.obj monte0.obj mpc_obs.obj mt64.obj &
  orb_func.obj orb_fun2.obj pl_cache.obj roots.obj &
  runge.obj sm_vsop.obj sr.obj tle_out.obj sigma.obj

find_orb.lib: $(OBJS)
   %write wfind.lrf $(OBJS)
   wlib -q -n -t $@ @wfind.lrf
   -del wfind.lrf

find_orb.exe: findorb.obj find_orb.lib $(CURSES_OBJ)
   wlink N find_orb.exe f findorb.obj l find_orb.lib l wafuncs.lib $(LINK_OPTS)

fo.exe: fo.obj find_orb.lib
   wlink N fo.exe f fo.obj l find_orb.lib l wafuncs.lib

fo_serve.exe: fo_serve.obj find_orb.lib cgi_func.obj
   wlink N fo_serve.exe f fo.obj l find_orb.lib l wafuncs.lib f cgi_func.obj

clean:
   -rm $(OBJS)
   -rm $(CURSES_OBJ)
   -rm cgi_func.obj
   -rm fo.exe
   -rm fo.obj
   -rm fo_serve.exe
   -rm fo_serve.obj
   -rm find_orb.exe
   -rm find_orb.lib
   -rm find_orb.map
   -rm debug.txt
   -rm observe.txt
   -rm covar.txt
   -rm guide.txt
   -rm elements.txt
   -rm mpc_fmt.txt
   -rm gauss.out
   -rm sof.txt
   -rm wfind.lrf
   -rm findorb.obj
   -rm monte.txt

CFLAGS=/Ox /W3 /4r /s /j /zq

.cpp.obj:
   wpp386 $(CFLAGS) $<

