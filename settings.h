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

#if !defined(AFX_SETTINGS_H__7104F85B_624C_4D56_983E_7A9DA275F443__INCLUDED_)
#define AFX_SETTINGS_H__7104F85B_624C_4D56_983E_7A9DA275F443__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// Settings.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CSettings dialog

class CSettings : public CDialog
{
// Construction
public:
   CSettings(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
   //{{AFX_DATA(CSettings)
   enum { IDD = IDD_SETTINGS };
   BOOL   m_heliocentric;
   int      m_element_precision;
   CString   m_constraints;
   CString   m_max_residual;
   double   m_monte_noise;
   CString   m_reference;
   BOOL   m_srp;
   int      m_comet_mags_total;
   int      m_physical_model;
   BOOL   m_precise_residuals;
   int      m_use_blunder_method;
   int      m_overobserving_ceiling;
   double   m_overobserving_time_span;
   double   m_probability_of_blunder;
   BOOL   m_alternative_elements;
   BOOL   m_use_weights;
   int      m_element_center;
   BOOL   m_debiasing;
   int      m_language;
   //}}AFX_DATA


// Overrides
   // ClassWizard generated virtual function overrides
   //{{AFX_VIRTUAL(CSettings)
   protected:
   virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
   //}}AFX_VIRTUAL

// Implementation
protected:

   // Generated message map functions
   //{{AFX_MSG(CSettings)
   virtual BOOL OnInitDialog();
   //}}AFX_MSG
   DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SETTINGS_H__7104F85B_624C_4D56_983E_7A9DA275F443__INCLUDED_)
