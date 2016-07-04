// about.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CAbout dialog

class CAbout : public CDialog
{
// Construction
public:
   CAbout(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
   //{{AFX_DATA(CAbout)
   enum { IDD = IDD_ABOUT };
      // NOTE: the ClassWizard will add data members here
   //}}AFX_DATA

// Implementation
protected:
   virtual void DoDataExchange(CDataExchange* pDX);   // DDX/DDV support

   // Generated message map functions
   //{{AFX_MSG(CAbout)
   virtual BOOL OnInitDialog();
   //}}AFX_MSG
   DECLARE_MESSAGE_MAP()
};
