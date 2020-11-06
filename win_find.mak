# Makefile for Microsoft C/C++ to build a Windows executable with MFC.
# Add BITS_32=Y for 32-bit code.  Assumes various libraries
# ('lunar' or 'lunar64', 'sat_code', 'jpl_eph') will be found.

OBJS=about.obj b32_eph.obj bc405.obj bias.obj clipfunc.obj \
  collide.obj conv_ele.obj details.obj eigen.obj elem2tle.obj elem_ou2.obj \
  elem_out.obj ephem0.obj ephem.obj expcalc.obj generic.obj gauss.obj \
  geo_pot.obj healpix.obj lsquare.obj miscell.obj \
  monte0.obj errors.obj monte.obj  \
  mpc_obs.obj nanosecs.obj orbitdlg.obj \
  orb_func.obj orb_fun2.obj pl_cache.obj roots.obj  \
  runge.obj settings.obj shellsor.obj stackall.obj \
  sigma.obj simplex.obj sm_vsop.obj sr.obj stdafx.obj

COMMON_FLAGS=-nologo -W3 -EHsc -c -FD -D_CRT_SECURE_NO_WARNINGS

!ifdef BITS_32
BITS=32
# COMMON_LINK=lunar64.lib jpleph.lib sat_code.lib /nologo /stack:0x8800 /subsystem:windows /MACHINE:IX86
!else
BITS=64
!endif

RM=del
EXE_NAME=find_o$(BITS).exe
COMMON_LINK=lunar$(BITS).lib jpleph$(BITS).lib sat_code$(BITS).lib /nologo /stack:0x8800 /subsystem:windows -ENTRY:"wWinMainCRTStartup"

!ifdef DEBUG
CFLAGS=-MDd -Gm -Zi -Od -D "_DEBUG" -Fp"find_orb.pch" $(COMMON_FLAGS) -D "_AFXDLL"
LINK_FLAGS=$(COMMON_LINK) -profile -map:"find_orb.map" -debug
!else
CFLAGS=-MT -O1 -D "NDEBUG" $(COMMON_FLAGS)
LINK_FLAGS=$(COMMON_LINK) -incremental:no
!endif

all: $(EXE_NAME)

$(EXE_NAME):              find_orb.obj $(OBJS) find_orb.res
     link $(LINK_FLAGS) find_orb.obj $(OBJS) find_orb.res
     del $(EXE_NAME)
     ren find_orb.exe $(EXE_NAME)

.cpp.obj:
   cl $(CFLAGS) $<

clean:
   $(RM) $(OBJS)
   $(RM) fo.obj find_orb.obj fo_serve.obj
   $(RM) covar.txt covar?.txt debug.txt eleme?.txt elements.txt
   $(RM) ephemeri.txt find_orb$(EXE) fo$(EXE) gauss.out monte.txt monte?.txt
   $(RM) mpc_f?.txt mpc_fmt.txt mpec.htm obser?.txt observe.txt
   $(RM) residual.txt state.txt state?.txt virtu?.txt virtual.txt
   $(RM) sr_elems.txt mpcorb.dat fo_serve.cgi find_orb.res
   $(RM) elements.txt covar.txt gauss.out
   $(RM) find_orb.exp vc50.pdb obs_temp.txt guide.txt
   $(RM) find_orb.map find_orb.pdb find_orb.lib vc*.idb

clean_temp:
   $(RM) covar.txt covar?.txt debug.txt eleme?.txt elements.txt
   $(RM) ephemeri.txt                         gauss.out monte.txt monte?.txt
   $(RM) mpc_f?.txt mpc_fmt.txt mpec.htm obser?.txt observe.txt
   $(RM) residual.txt state.txt state?.txt virtu?.txt virtual.txt
   $(RM) sr_elems.txt mpcorb.dat              find_orb.res
   $(RM) elements.txt covar.txt gauss.out
   $(RM) find_orb.exp vc50.pdb obs_temp.txt guide.txt
   $(RM) find_orb.map find_orb.pdb find_orb.lib vc*.idb

