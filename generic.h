#if !defined(AFX_GENERIC_H__17872CB7_EBB5_49E6_8086_CA2FB7B215C1__INCLUDED_)
#define AFX_GENERIC_H__17872CB7_EBB5_49E6_8086_CA2FB7B215C1__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
// Generic.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CGenericEntry dialog

class CGenericEntry : public CDialog
{
// Construction
public:
   CGenericEntry(CWnd* pParent = NULL);   // standard constructor
   CString m_caption;

// Dialog Data
   //{{AFX_DATA(CGenericEntry)
   enum { IDD = IDD_GENERIC_ENTRY };
   CString   m_text;
   //}}AFX_DATA


// Overrides
   // ClassWizard generated virtual function overrides
   //{{AFX_VIRTUAL(CGenericEntry)
   protected:
   virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
   //}}AFX_VIRTUAL

// Implementation
protected:

   // Generated message map functions
   //{{AFX_MSG(CGenericEntry)
   virtual BOOL OnInitDialog();
   //}}AFX_MSG
   DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_GENERIC_H__17872CB7_EBB5_49E6_8086_CA2FB7B215C1__INCLUDED_)
