// ephem.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CEphem dialog

class CEphem : public CDialog
{
// Construction
public:
   CEphem(COrbitDlg* pParent = NULL);   // standard constructor
   double orbit[6], epoch;
   double jd, ht_in_meters;
   OBSERVE *obs;
   int n_obs;
   const char *obj_name;
   int GetEphemerisBitmask( );

// Dialog Data
   //{{AFX_DATA(CEphem)
   enum { IDD = IDD_MAKE_EPHEMERIS };
   CString   m_day;
   int      m_number_steps;
   CString   m_lat;
   CString   m_lon;
   CString   m_ephem_step;
   BOOL   m_alt_az;
   BOOL   m_motion;
   int      m_ephem_type;
   BOOL   m_phase_angle;
   BOOL   m_radial_velocity;
   BOOL   m_separate_motions;
   BOOL   m_round_step;
   BOOL   m_topo_ecliptic;
   BOOL   m_helio_ecliptic;
   BOOL   m_phase_angle_bisector;
   BOOL   m_suppress_unobservable;
   BOOL   m_visibility;
   CString   m_mpc_code;
   int      m_use_mpc_code;
   double   m_mag_limit;
   BOOL   m_show_sigmas;
   BOOL   m_human_readable;
   BOOL   m_ground_track;
   BOOL   m_speed;
   //}}AFX_DATA


// Implementation
protected:
   int ephemeris_and_pseudo_mpec_made;
   CFont list_box_font;
   virtual void DoDataExchange(CDataExchange* pDX);   // DDX/DDV support
   void CreateB32Ephemeris( const char *filename);
   void set_jd_from_xtrols( char *err_msg);
   CRect OriginalDlgRect;
   void OnGetMinMaxInfo(MINMAXINFO FAR* lpMMI);

   // Generated message map functions
   //{{AFX_MSG(CEphem)
   afx_msg void OnClickedSave();
   afx_msg void OnClickedGo();
   virtual BOOL OnInitDialog();
   afx_msg void OnPseudoMpec();
   afx_msg void OnSize(UINT nType, int cx, int cy);
   afx_msg void OnCopy();
   afx_msg void OnDestroy();
   afx_msg void OnNow();
   //}}AFX_MSG
   DECLARE_MESSAGE_MAP()
};
