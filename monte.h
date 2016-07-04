#if !defined(AFX_MONTE_H__817675B2_4DB0_474A_BB96_CF20EB5651AB__INCLUDED_)
#define AFX_MONTE_H__817675B2_4DB0_474A_BB96_CF20EB5651AB__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// Monte.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CMonteCarlo dialog

class CMonteCarlo : public CDialog
{
// Construction
public:
   CMonteCarlo(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
   //{{AFX_DATA(CMonteCarlo)
   enum { IDD = IDD_MONTE_CARLO };
   CString   m_gaussian_noise;
   double   m_max_ecc;
   double   m_max_incl;
   CString   m_min_range;
   BOOL   m_statistical;
   CString   m_max_range;
   //}}AFX_DATA


// Overrides
   // ClassWizard generated virtual function overrides
   //{{AFX_VIRTUAL(CMonteCarlo)
   protected:
   virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
   //}}AFX_VIRTUAL

// Implementation
protected:

   // Generated message map functions
   //{{AFX_MSG(CMonteCarlo)
      // NOTE: the ClassWizard will add member functions here
   //}}AFX_MSG
   DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MONTE_H__817675B2_4DB0_474A_BB96_CF20EB5651AB__INCLUDED_)
