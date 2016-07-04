// find_orb.h : main header file for the FIND_ORB application
//

#ifndef __AFXWIN_H__
   #error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"      // main symbols

/////////////////////////////////////////////////////////////////////////////
// CFind_orbApp:
// See find_orb.cpp for the implementation of this class
//

class CFind_orbApp : public CWinApp
{
public:
   CFind_orbApp();

// Overrides
   virtual BOOL InitInstance();

private:
// Implementation

   //{{AFX_MSG(CFind_orbApp)
   //}}AFX_MSG
   DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////
