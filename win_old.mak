# "Old" make file,  for Windows ME and MSVC 5.0.  Still works,
# as of early 2015,  to create a 32-bit Windows Find_Orb.

all: find_o32.exe

RM=rm

OBJS=about.obj clipfunc.obj ephem.obj generic.obj monte.obj  \
  orbitdlg.obj settings.obj stdafx.obj \
  b32_eph.obj bc405.obj collide.obj conv_ele.obj eigen.obj \
  elem2tle.obj elem_out.obj ephem0.obj gauss.obj \
  jpleph.obj lsquare.obj moid4.obj monte0.obj mpc_obs.obj \
  orb_func.obj orb_fun2.obj pl_cache.obj roots.obj runge.obj \
  sigma.obj sm_vsop.obj sr.obj tle_out.obj

!ifdef DEBUG
CFLAGS=/nologo /MDd /W3 /Gm /GX /Zi /Od /D "_DEBUG" /D "_AFXDLL"\
 /Fp"find_orb.pch" /FD /c
LINK32_FLAGS=lunar.lib /nologo /stack:0x8800 /subsystem:windows\
 /profile /map:"find_orb.map" /debug /machine:IX86\
 /def:"find_orb.def" /out:"find_orb.exe"
!else
CFLAGS=-nologo -MT -W3 -GX -O1 -D "NDEBUG" -FD -c -D_CRT_SECURE_NO_WARNINGS
LINK32_FLAGS=lunar.lib /nologo /stack:0x8800 /subsystem:windows\
 /incremental:no /machine:IX86 /def:"find_orb.def" /out:"find_orb.exe"
!endif

find_o32.exe:              find_orb.obj $(OBJS) find_orb.res
     link $(LINK32_FLAGS)  find_orb.obj $(OBJS) find_orb.res
     del find_o32.exe
     ren find_orb.exe find_o32.exe
     rm find_orb.exp
     rm find_orb.lib

.cpp.obj:
   cl $(CFLAGS) $<

clean:
   $(RM) $(OBJS)
   $(RM) fo.obj findorb.obj fo_serve.obj
   $(RM) covar.txt covar?.txt debug.txt eleme?.txt elements.txt
   $(RM) ephemeri.txt find_orb$(EXE) fo$(EXE) gauss.out monte.txt monte?.txt
   $(RM) mpc_f?.txt mpc_fmt.txt mpec.htm obser?.txt observe.txt
   $(RM) residual.txt state.txt state?.txt virtu?.txt virtual.txt
   $(RM) sr_elems.txt mpcorb.dat fo_serve.cgi find_orb.res
   $(RM) mpc_fmt.txt elements.txt covar.txt gauss.out *.idb

