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
