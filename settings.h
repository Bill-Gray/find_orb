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
