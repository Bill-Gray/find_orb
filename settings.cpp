/* settings.cpp: 'settings' dialog for Windows Find_Orb

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

// Settings.cpp : implementation file
//

#include "stdafx.h"
#include "find_orb.h"
#include "Settings.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CSettings dialog

int reset_dialog_language( CDialog *dlg, const int dlg_number);  /* elem_out.cpp */

CSettings::CSettings(CWnd* pParent /*=NULL*/)
   : CDialog(CSettings::IDD, pParent)
{
   //{{AFX_DATA_INIT(CSettings)
   m_element_precision = 0;
   m_constraints = _T("");
   m_max_residual = _T("");
   m_monte_noise = 0.0;
   m_reference = _T("");
   m_srp = FALSE;
   m_comet_mags_total = -1;
   m_physical_model = -1;
   m_precise_residuals = FALSE;
   m_use_blunder_method = -1;
   m_overobserving_ceiling = 0;
   m_overobserving_time_span = 0.0;
   m_probability_of_blunder = 0.0;
   m_alternative_elements = FALSE;
   m_use_weights = FALSE;
   m_element_center = -1;
   m_debiasing = FALSE;
   //}}AFX_DATA_INIT
}


void CSettings::DoDataExchange(CDataExchange* pDX)
{
   CDialog::DoDataExchange(pDX);
   //{{AFX_DATA_MAP(CSettings)
   DDX_Text(pDX, IDC_ELEMENT_PRECISION, m_element_precision);
   DDX_Text(pDX, IDC_CONSTRAINTS, m_constraints);
   DDX_Text(pDX, IDC_MAX_RESIDUAL, m_max_residual);
   DDX_Text(pDX, IDC_MONTE_NOISE, m_monte_noise);
   DDX_Text(pDX, IDC_REFERENCE, m_reference);
   DDX_Radio(pDX, IDC_RADIO1, m_comet_mags_total);
   DDX_Radio(pDX, IDC_RADIO3, m_physical_model);
   DDX_Check(pDX, IDC_PRECISE_RESIDUALS, m_precise_residuals);
   DDX_Radio(pDX, IDC_RADIO6, m_use_blunder_method);
   DDX_Text(pDX, IDC_OVEROBSERVING_CEILING, m_overobserving_ceiling);
   DDX_Text(pDX, IDC_OVEROBSERVING_TIME_SPAN, m_overobserving_time_span);
   DDX_Text(pDX, IDC_PROBABILITY_OF_BLUNDER, m_probability_of_blunder);
   DDX_Check(pDX, IDC_ALT_ELEM, m_alternative_elements);
   DDX_Check(pDX, IDC_USE_WEIGHTS, m_use_weights);
   DDX_CBIndex(pDX, IDC_COMBO1, m_element_center);
   DDX_Check(pDX, IDC_DEBIASING, m_debiasing);
   //}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CSettings, CDialog)
   //{{AFX_MSG_MAP(CSettings)
   //}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSettings message handlers

BOOL CSettings::OnInitDialog()
{
   CDialog::OnInitDialog();

   // TODO: Add extra initialization here
   reset_dialog_language( this, 96000);

   return TRUE;  // return TRUE unless you set the focus to a control
                 // EXCEPTION: OCX Property Pages should return FALSE
}

