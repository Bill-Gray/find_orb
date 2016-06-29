# Microsoft Developer Studio Project File - Name="Find_orb" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=Find_orb - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Find_orb.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Find_orb.mak" CFG="Find_orb - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Find_orb - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "Find_orb - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Find_orb - Win32 Release"

# PROP BASE Use_MFC 1
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 1
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MT /W4 /GX /O1 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "JSAT_OBS" /FR /Yu"STDAFX.H" /FD /c
# ADD CPP /nologo /MT /W3 /GX /O1 /I "\myincl" /D "NDEBUG" /FR /Yu"STDAFX.H" /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /o NUL /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /o NUL /win32
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 oldnames.lib \lunar\winlunar.lib /nologo /stack:0x2800 /subsystem:windows /machine:IX86
# ADD LINK32 \lunar\posttest\lunar.lib /nologo /stack:0x8800 /subsystem:windows /machine:IX86
# SUBTRACT LINK32 /pdb:none
# Begin Special Build Tool
SOURCE=$(InputPath)
PostBuild_Cmds=copy release\find_orb.exe find_o32.exe
# End Special Build Tool

!ELSEIF  "$(CFG)" == "Find_orb - Win32 Debug"

# PROP BASE Use_MFC 1
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 2
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /FR /Yu"STDAFX.H" /FD /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /Zi /Od /I "\myincl" /D "_DEBUG" /D "_AFXDLL" /FR /Yu"STDAFX.H" /FD /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /o NUL /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /o NUL /win32
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 oldnames.lib \lunar\winlunar.lib /nologo /stack:0x2800 /subsystem:windows /debug /machine:IX86 /pdbtype:sept
# ADD LINK32 \lunar\posttest\lunar.lib /nologo /stack:0x8800 /subsystem:windows /profile /map /debug /machine:IX86
# SUBTRACT LINK32 /nodefaultlib

!ENDIF 

# Begin Target

# Name "Find_orb - Win32 Release"
# Name "Find_orb - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;hpj;bat"
# Begin Source File

SOURCE=.\About.cpp

!IF  "$(CFG)" == "Find_orb - Win32 Release"

# ADD CPP /Zp8

!ELSEIF  "$(CFG)" == "Find_orb - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\b32_eph.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\bc405.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\clipfunc.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\collide.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Conv_ele.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\eigen.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\elem2tle.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\elem_out.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Ephem.cpp

!IF  "$(CFG)" == "Find_orb - Win32 Release"

# SUBTRACT CPP /YX /Yc /Yu

!ELSEIF  "$(CFG)" == "Find_orb - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Ephem0.cpp

!IF  "$(CFG)" == "Find_orb - Win32 Release"

# ADD CPP /O1
# SUBTRACT CPP /YX /Yc /Yu

!ELSEIF  "$(CFG)" == "Find_orb - Win32 Debug"

# ADD CPP /O2
# SUBTRACT CPP /YX /Yc /Yu

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Find_orb.cpp
# End Source File
# Begin Source File

SOURCE=.\Find_orb.def
# End Source File
# Begin Source File

SOURCE=.\Find_orb.rc
# End Source File
# Begin Source File

SOURCE=.\Gauss.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Generic.cpp
# End Source File
# Begin Source File

SOURCE=.\Get_pert.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\jpl\pl\Jpleph.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=..\Image\Lsquare.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\moid4.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Monte.cpp
# End Source File
# Begin Source File

SOURCE=.\monte0.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Mpc_obs.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=..\prng\mt64.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\orb_fun2.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Orb_func.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Orbitdlg.cpp
# End Source File
# Begin Source File

SOURCE=.\pl_cache.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Roots.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Runge.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Settings.cpp
# End Source File
# Begin Source File

SOURCE=.\sigma.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\sm_vsop.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Sr.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\Stdafx.cpp
# ADD BASE CPP /Yc"STDAFX.H"
# ADD CPP /Yc"STDAFX.H"
# End Source File
# Begin Source File

SOURCE=..\sattest\neoklis\tle_out.cpp
# SUBTRACT CPP /YX /Yc /Yu
# End Source File
# Begin Source File

SOURCE=.\weight.cpp
# PROP Exclude_From_Build 1
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\de_cache.h
# End Source File
# Begin Source File

SOURCE=.\Generic.h
# End Source File
# Begin Source File

SOURCE=.\Monte.h
# End Source File
# Begin Source File

SOURCE=.\monte0.h
# End Source File
# Begin Source File

SOURCE=.\Settings.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;cnt;rtf;gif;jpg;jpeg;jpe"
# Begin Source File

SOURCE=.\RES\FIND_ORB.ICO
# End Source File
# Begin Source File

SOURCE=.\res\find_orb.rc2
# End Source File
# End Group
# End Target
# End Project
