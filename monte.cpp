/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA. */

// Monte.cpp : implementation file
//

#include "stdafx.h"
#include "find_orb.h"
#include "Monte.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CMonteCarlo dialog


CMonteCarlo::CMonteCarlo(CWnd* pParent /*=NULL*/)
   : CDialog(CMonteCarlo::IDD, pParent)
{
   //{{AFX_DATA_INIT(CMonteCarlo)
   m_gaussian_noise = _T("");
   m_max_ecc = 0.0;
   m_max_incl = 0.0;
   m_min_range = _T("");
   m_statistical = FALSE;
   m_max_range = _T("");
   //}}AFX_DATA_INIT
}


void CMonteCarlo::DoDataExchange(CDataExchange* pDX)
{
   CDialog::DoDataExchange(pDX);
   //{{AFX_DATA_MAP(CMonteCarlo)
   DDX_Text(pDX, IDC_GAUSSIAN, m_gaussian_noise);
   DDX_Text(pDX, IDC_MAX_ECC, m_max_ecc);
   DDX_Text(pDX, IDC_MAX_INCL, m_max_incl);
   DDX_Text(pDX, IDC_MIN_RANGE, m_min_range);
   DDX_Check(pDX, IDC_STATISTICAL, m_statistical);
   DDX_Text(pDX, IDC_MAX_RANGE, m_max_range);
   //}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CMonteCarlo, CDialog)
   //{{AFX_MSG_MAP(CMonteCarlo)
      // NOTE: the ClassWizard will add message map macros here
   //}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMonteCarlo message handlers
