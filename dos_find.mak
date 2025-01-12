# Makefile for Microsoft C/C++ to build an executable linked to PDCursesMod.
# Add BITS_32=Y for 32-bit code.  Assumes various libraries ('pdcurses',
# 'lunar' or 'lunar64', 'sat_code', 'jpl_eph') will be found.  Note that
# I'm using the PDCursesMod from my GitHub repository.
#    Should work reasonably well in the future;  I have automated CI
# builds running for it on GitHub(R) now.

all: find_orb.exe fo.exe

OBJS=ades_out.obj bc405.obj bias.obj collide.obj   \
  conv_ele.obj details.obj eigen.obj elem2tle.obj elem_out.obj  \
  elem_ou2.obj ephem0.obj errors.obj expcalc.obj gauss.obj  \
  geo_pot.obj healpix.obj lsquare.obj miscell.obj  \
  monte0.obj mpc_obs.obj  \
  orb_func.obj orb_fun2.obj pl_cache.obj roots.obj runge.obj \
  shellsor.obj sigma.obj simplex.obj sm_vsop.obj sr.obj stackall.obj

CCLIBS      = user32.lib gdi32.lib advapi32.lib shell32.lib comdlg32.lib
!ifdef BITS_32
CFLAGS=-Ot -W3 -nologo -MT -I../PDCursesMod
ADD_LIBS    = sat_code32.lib jpleph32.lib lunar.lib
RM=rm
!else
CFLAGS=-Ot -W3 -nologo -MT -I../PDCursesMod -D_CRT_SECURE_NO_WARNINGS
ADD_LIBS    = sat_code64.lib jpleph64.lib lunar64.lib
RM=del
!endif

eph2tle.exe: eph2tle.obj conv_ele.obj elem2tle.obj simplex.obj lsquare.obj
   link /out:eph2tle.exe eph2tle.obj conv_ele.obj elem2tle.obj simplex.obj \
                            lsquare.obj $(ADD_LIBS)

cssfield.exe: cssfield.cpp
   cl $(CFLAGS) cssfield.cpp $(ADD_LIBS)

find_orb.exe:               findorb.obj $(OBJS) clipfunc.obj getstrex.obj
     link /out:find_orb.exe findorb.obj $(OBJS) clipfunc.obj getstrex.obj \
                       pdcurses.lib user32.lib winmm.lib $(CCLIBS) $(ADD_LIBS)

fo.exe:                     fo.obj $(OBJS)
     link /out:fo.exe       fo.obj $(OBJS) $(ADD_LIBS)

fo_serve.exe:               fo_serve.obj $(OBJS)
     link /out:fo_serve.exe fo_serve.obj $(OBJS) $(ADD_LIBS)

roottest.exe: roottest.obj roots.obj
   link /out:roottest.exe roottest.obj roots.obj

.cpp.obj:
   cl -c $(CFLAGS) $<

clean:
   $(RM) $(OBJS)
   $(RM) clipfunc.obj cssfield.obj eph2tle.obj roottest.obj
   $(RM) fo.obj findorb.obj fo_serve.obj getstrex.obj
   $(RM) covar.txt covar?.txt debug.txt eleme?.txt elements.txt
   $(RM) ephemeri.txt find_orb.exe fo.exe gauss.out monte.txt monte?.txt
   $(RM) mpc_f?.txt mpc_fmt.txt mpec.htm obser?.txt observe.txt
   $(RM) residual.txt state.txt state?.txt virtu?.txt virtual.txt
   $(RM) sr_elems.txt mpcorb.dat fo_serve.cgi find_orb.res
   $(RM) elements.txt covar.txt gauss.out
   $(RM) find_orb.exp vc*.pdb obs_temp.txt guide.txt
   $(RM) find_orb.map find_orb.pdb find_orb.lib vc*.idb
   $(RM) cssfield.exe eph2tle.exe roottest.exe
   $(RM) find_o32.exe find_o64.exe

clean_temp:
   $(RM) bc405pre.txt
   $(RM) cmt_sof.txt
   $(RM) combined.json
   $(RM) covar.txt
   $(RM) covar.json
   $(RM) debug.txt
   $(RM) dummy.txt
   $(RM) elements.txt
   $(RM) elements.json
   $(RM) elem_short.json
   $(RM) ephemeri.txt
   $(RM) ephemeri.json
   $(RM) gauss.out
   $(RM) guide.txt
   $(RM) linkage.json
   $(RM) lock.txt
   $(RM) monte.txt
   $(RM) mpcorb.dat
   $(RM) mpc_fmt.txt
   $(RM) mpec.htm
   $(RM) observe.txt
   $(RM) observe.xml
   $(RM) obs_temp.txt
   $(RM) residual.txt
   $(RM) sof.txt
   $(RM) sr_elems.txt
   $(RM) state.txt
   $(RM) total.json
   $(RM) vectors.txt
   $(RM) virtual.txt
