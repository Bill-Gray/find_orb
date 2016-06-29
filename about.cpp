/* about.cpp: 'about' dialog for Windows Find_Orb

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

// about.cpp : implementation file
//

#include "stdafx.h"
#include "find_orb.h"
#include "about.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

int reset_dialog_language( CDialog *dlg, const int dlg_number);  /* elem_out.cpp */

/////////////////////////////////////////////////////////////////////////////
// CAbout dialog

CAbout::CAbout(CWnd* pParent /*=NULL*/)
   : CDialog(CAbout::IDD, pParent)
{
   //{{AFX_DATA_INIT(CAbout)
      // NOTE: the ClassWizard will add member initialization here
   //}}AFX_DATA_INIT
}

void CAbout::DoDataExchange(CDataExchange* pDX)
{
   CDialog::DoDataExchange(pDX);
   //{{AFX_DATA_MAP(CAbout)
      // NOTE: the ClassWizard will add DDX and DDV calls here
   //}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAbout, CDialog)
   //{{AFX_MSG_MAP(CAbout)
   //}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CAbout message handlers

BOOL CAbout::OnInitDialog()
{
   CDialog::OnInitDialog();
   // TODO: Add extra initialization here
   reset_dialog_language( this, 98000);

   return TRUE;  // return TRUE  unless you set the focus to a control
}
