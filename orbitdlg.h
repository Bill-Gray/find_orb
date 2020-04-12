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

// orbitdlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// COrbitDlg dialog

         /* Somewhere along the line -- not exactly sure where -- */
         /* OnTimer went from taking an unsigned int to taking a  */
         /* pointer to an unsigned int.                           */
#if _MSC_VER <= 1100
   #define TIMER_TYPE unsigned
#else
   #define TIMER_TYPE UINT_PTR
#endif

class COrbitDlg : public CDialog
{
// Construction
public:
   COrbitDlg(CWnd* pParent = NULL);   // standard constructor

// Dialog Data
   //{{AFX_DATA(COrbitDlg)
   enum { IDD = IDD_FIND_ORB };
   double   m_step_size;
   CString   m_epoch;
   CString   m_r1;
   CString   m_r2;
   //}}AFX_DATA

// Implementation
protected:
   int RunMonteCarlo( void);
   int precise_residuals;
   int monte_carlo;
   int show_commented_elements;
   int compute_covariance;
   virtual void DoDataExchange(CDataExchange* pDX);   // DDX/DDV support
   void ResetPerturbers();
   void Reset_r1_and_r2( );
   int ImproveOrbitSolution( int full_step, int n_repeats);
   void UpdateElementDisplay( int update_orbit);
   void UpdateResidualDisplay();
   void LoadAnObject( const int obj_idx);
   void LoadAFile( const char *filename);
   int GetPerturberMask();
   CString curr_file_name;
   CString curr_object_name;
   CString constraints;
// CFont elements_font;
   time_t curr_file_time;
   double max_residual_for_filtering;
   double orbit[6], orbit_epoch, monte_noise;
   int n_obs, element_format, element_precision, n_objects;
   ephem_option_t ephemeris_output_options;
   int bad_elements;
   double blunder_probability;
   void FAR *obs_data;
   OBJECT_INFO *obj_info;
   CRect OriginalDlgRect;
   void OnGetMinMaxInfo(MINMAXINFO FAR* lpMMI);
   void AdjustControls( void);

   // Generated message map functions
   //{{AFX_MSG(COrbitDlg)
   afx_msg void OnClickedFullStep();
   afx_msg void OnClickedHerget();
   afx_msg void OnClickedOpen();
   virtual BOOL OnInitDialog();
   afx_msg void OnDblclkObject();
   afx_msg void OnClickedSave();
   afx_msg void OnSelchangeResiduals();
   afx_msg void OnClickedMakeEphemeris();
   afx_msg void OnClickedSaveResids();
   afx_msg void OnDblclkResiduals();
 afx_msg void OnClickedAbout();
   afx_msg void OnClickedVaisala();
   afx_msg void OnClickedAutoSolve();
// afx_msg void OnChar(UINT nChar, UINT nRepCnt, UINT nFlags);
   afx_msg void OnMonteCarlo();
   afx_msg void OnTimer( TIMER_TYPE pfnIDEvent);
   afx_msg void OnGauss();
   afx_msg void OnWorst();
   afx_msg void OnDoubleclickedWorst();
   afx_msg void OnFilterObs();
   afx_msg void OnSelcancelListAsteroids();
   afx_msg void OnAsteroids();
   afx_msg void OnSettings();
   afx_msg void OnSelchangeListAsteroids();
   afx_msg void OnSetWeight();
   afx_msg void OnToggleObs();
   afx_msg void OnOrbitalElements();
   afx_msg void OnDestroy();
   afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
   afx_msg void OnSize(UINT nType, int cx, int cy);
   afx_msg void OnTogglePerturbers();
   afx_msg void OnSimplex();
   afx_msg void OnStatisticalRanging();
   afx_msg HBRUSH OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor);
   afx_msg void OnKillfocusEpoch();
   afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
   //}}AFX_MSG
   DECLARE_MESSAGE_MAP()
};
