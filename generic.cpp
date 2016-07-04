/* generic.cpp: dialog to give a prompt and get a string from user

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

// Generic.cpp : implementation file
//

#include "stdafx.h"
#include "find_orb.h"
#include "Generic.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CGenericEntry dialog


CGenericEntry::CGenericEntry(CWnd* pParent /*=NULL*/)
   : CDialog(CGenericEntry::IDD, pParent)
{
   //{{AFX_DATA_INIT(CGenericEntry)
   m_text = _T("");
   //}}AFX_DATA_INIT
}


void CGenericEntry::DoDataExchange(CDataExchange* pDX)
{
   CDialog::DoDataExchange(pDX);
   //{{AFX_DATA_MAP(CGenericEntry)
   DDX_Text(pDX, IDC_EDIT1, m_text);
   //}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CGenericEntry, CDialog)
   //{{AFX_MSG_MAP(CGenericEntry)
   //}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CGenericEntry message handlers

BOOL CGenericEntry::OnInitDialog()
{
   CDialog::OnInitDialog();

   // TODO: Add extra initialization here

   SetDlgItemText( IDC_STATIC1, m_caption);
   return TRUE;  // return TRUE unless you set the focus to a control
                 // EXCEPTION: OCX Property Pages should return FALSE
}
