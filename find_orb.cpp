/* find_orb.cpp : Defines the class behaviors for the application.

Copyright (C) 2010, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.    */


#include "stdafx.h"
#include <stdint.h>
#include "mpc_obs.h"
#include "find_orb.h"
#include "orbitdlg.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CFind_orbApp

BEGIN_MESSAGE_MAP(CFind_orbApp, CWinApp)
   //{{AFX_MSG_MAP(CFind_orbApp)
   //}}AFX_MSG_MAP
   // Standard file based document commands
   ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
   ON_COMMAND(ID_FILE_OPEN, CWinApp::OnFileOpen)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFind_orbApp construction

CFind_orbApp::CFind_orbApp()
{
   // TODO: add construction code here,
   // Place all significant initialization in InitInstance
}

/////////////////////////////////////////////////////////////////////////////
// The one and only CFind_orbApp object

CFind_orbApp NEAR theApp;

/////////////////////////////////////////////////////////////////////////////
// CFind_orbApp initialization

int debug_level = 0;

BOOL CFind_orbApp::InitInstance()
{
   // Standard initialization
   // If you are not using these features and wish to reduce the size
   //  of your final executable, you should remove from the following
   //  the specific initialization routines you do not need.

   LoadStdProfileSettings();  // Load standard INI file options (including MRU)

   if (m_lpCmdLine[0] != '\0')
   {
      // TODO: add command line processing here
      const char *language = strstr( (const char *)m_lpCmdLine, "-l");
      extern char findorb_language;       /* defaults to 'e' for English */

      if( language && language[2])
         findorb_language = language[2];
   }

   COrbitDlg dlg;

   dlg.m_r1 = dlg.m_r2 = "1";
   dlg.m_step_size = 3.;
   dlg.DoModal( );

   return TRUE;
}
