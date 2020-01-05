/* ephem.cpp: code to show ephemerides in Windows Find_Orb

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
// ephem.cpp : implementation file
//

#include "stdafx.h"
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "find_orb.h"
#include "mpc_obs.h"
#include "orbitdlg.h"
#include "ephem.h"
#include "watdefs.h"
#include "lunar.h"
#include "afuncs.h"
#include "date.h"
#include "comets.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

#ifndef CA2T
   #define CT2A(x) ((const char *)x)
   #define CA2T(x) ((const char *)x)
#endif


#define PI 3.141592653589793238462643383279502884197169399375105

const char *get_find_orb_text( const int index);      /* elem_out.cpp */
double current_jd( void);                       /* elem_out.cpp */
void remove_trailing_cr_lf( char *buff);      /* ephem0.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);
int make_pseudo_mpec( const char *mpec_filename, const char *obj_name);
                                              /* ephem0.cpp */
int reset_dialog_language( CDialog *dlg, const int dlg_number);  /* elem_out.cpp */
int lat_alt_to_parallax( const double lat, const double ht_in_meters,
            double *rho_cos_phi, double *rho_sin_phi, const int planet_idx);
int debug_printf( const char *format, ...)                 /* runge.cpp */
#ifdef __GNUC__
         __attribute__ (( format( printf, 1, 2)))
#endif
;
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile);
int add_ephemeris_details( FILE *ofile, const double start_jd,  /* b32_eph.c */
                                               const double end_jd);
int create_b32_ephemeris( const char *filename, const double epoch,
                const double *orbit, const int n_steps,         /* b32_eph.c */
                const double ephem_step, const double jd_start);
void set_window_placement( HWND hwnd, const char *ibuff);      /* orbitdlg.c */
char *get_placement_text( HWND hwnd);                          /* orbitdlg.c */
char *mpc_station_name( char *station_data);       /* mpc_obs.cpp */
int remove_rgb_code( char *buff);                              /* ephem.cpp */
void compute_variant_orbit( double *variant, const double *ref_orbit,
                     const double n_sigmas);                /* orb_func.cpp */

/////////////////////////////////////////////////////////////////////////////
// CEphem dialog

CEphem::CEphem(COrbitDlg* pParent /*=NULL*/)
   : CDialog(CEphem::IDD, pParent)
{
   //{{AFX_DATA_INIT(CEphem)
   m_day = "";
   m_number_steps = 0;
   m_lat = "";
   m_lon = "";
   m_ephem_step = _T("");
   m_ephem_type = -1;
   m_phase_angle = FALSE;
   m_radial_velocity = FALSE;
   m_separate_motions = FALSE;
   m_round_step = FALSE;
   m_topo_ecliptic = FALSE;
   m_helio_ecliptic = FALSE;
   m_phase_angle_bisector = FALSE;
   m_suppress_unobservable = FALSE;
   m_visibility = FALSE;
   m_mpc_code = _T("");
   m_use_mpc_code = -1;
   m_mag_limit = 0.0;
   m_show_sigmas = FALSE;
   m_human_readable = FALSE;
   m_ground_track = FALSE;
   m_speed = FALSE;
   //}}AFX_DATA_INIT
   OriginalDlgRect.top = OriginalDlgRect.bottom = 0;
   ephemeris_and_pseudo_mpec_made = 0;
}

void CEphem::DoDataExchange(CDataExchange* pDX)
{
   CDialog::DoDataExchange(pDX);
   //{{AFX_DATA_MAP(CEphem)
   DDX_Text(pDX, IDC_EPHEM_DAY, m_day);
   DDX_Text(pDX, IDC_EPHEM_NUM_STEPS, m_number_steps);
   DDX_Text(pDX, IDC_LAT, m_lat);
   DDX_Text(pDX, IDC_LON, m_lon);
   DDX_Text(pDX, IDC_EPHEM_STEP, m_ephem_step);
   DDX_Check(pDX, IDC_ALT_AZ, m_alt_az);
   DDX_Check(pDX, IDC_MOTION, m_motion);
   DDX_Radio(pDX, IDC_OBSERVABLES, m_ephem_type);
   DDX_Check(pDX, IDC_PHASE_ANGLE, m_phase_angle);
   DDX_Check(pDX, IDC_RADIAL_VELOCITY, m_radial_velocity);
   DDX_Check(pDX, IDC_SEPARATE_MOTIONS, m_separate_motions);
   DDX_Check(pDX, IDC_ROUND_STEP, m_round_step);
   DDX_Check(pDX, IDC_TOPO_ECLIPTIC, m_topo_ecliptic);
   DDX_Check(pDX, IDC_HELIO_ECLIPTIC, m_helio_ecliptic);
   DDX_Check(pDX, IDC_PHASE_ANGLE_BISECTOR, m_phase_angle_bisector);
   DDX_Check(pDX, IDC_SUPPRESS_UNOBSERVABLE, m_suppress_unobservable);
   DDX_Check(pDX, IDC_VISIBILITY, m_visibility);
   DDX_Text(pDX, IDC_MPC_CODE, m_mpc_code);
   DDX_Radio(pDX, IDC_LOCATION1, m_use_mpc_code);
   DDX_Text(pDX, IDC_MAG_LIMIT, m_mag_limit);
   DDX_Check(pDX, IDC_SIGMAS, m_show_sigmas);
   DDX_Check(pDX, IDC_HUMAN_READABLE, m_human_readable);
   DDX_Check(pDX, IDC_GROUND_TRACK, m_ground_track);
   DDX_Check(pDX, IDC_SPEED, m_speed);
   //}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CEphem, CDialog)
   //{{AFX_MSG_MAP(CEphem)
   ON_BN_CLICKED(IDC_SAVE, OnClickedSave)
   ON_BN_CLICKED(IDC_GO, OnClickedGo)
   ON_BN_CLICKED(IDC_MPEC, OnPseudoMpec)
   ON_WM_SIZE()
   ON_BN_CLICKED(IDC_COPY, OnCopy)
   ON_WM_DESTROY()
   ON_WM_GETMINMAXINFO()
   ON_WM_CHAR()
   ON_BN_CLICKED(IDC_NOW, OnNow)
   //}}AFX_MSG_MAP
END_MESSAGE_MAP()

void get_file_from_dialog( int is_open, const char *default_ext,
                           const char *filter, char *buff, const char *path);

/////////////////////////////////////////////////////////////////////////////
// CEphem message handlers

extern const char *ephemeris_filename;         /* "ephemeri.txt" */
int save_ephemeris_file( const char *filename);   /* ephem0.cpp */

void CEphem::OnClickedSave()
{
   // TODO: Add your control notification handler code here
   char filename[_MAX_DIR];

   get_file_from_dialog( FALSE, "", "*.*", filename, NULL);
   if( *filename)
      if( *m_lat == 'b')
         CreateB32Ephemeris( filename);
      else
         save_ephemeris_file( filename);
}

void CEphem::set_jd_from_xtrols( char *err_msg)
{
   static const double jan_1970 = 2440587.5;
   static const double curr_time = current_jd( );

   if( err_msg)
      *err_msg = '\0';
   UpdateData( TRUE);     /* Get changes made in edit boxes */
#ifdef _UNICODE
   CT2A ascii( m_day);
   jd = get_time_from_string( curr_time, ascii.m_psz,
                       FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, NULL);
#else
   jd = get_time_from_string( curr_time, m_day,
                       FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, NULL);
#endif
}

double atof_cstring( const CString str)
{
   return(  _tcstod( (const TCHAR *)str, NULL));
}


const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */

void CEphem::CreateB32Ephemeris( const char *filename)
{
   int rval;

   set_jd_from_xtrols( NULL);
   rval = create_b32_ephemeris( filename, epoch, orbit, m_number_steps,
               atof_cstring( m_ephem_step), jd);
}

static const char *mpec_filename = "mpec.htm";

int CEphem::GetEphemerisBitmask( )
{
   int rval = m_ephem_type;

   if( m_alt_az)
      rval |= OPTION_ALT_AZ_OUTPUT;
   if( m_motion)
      rval |= OPTION_MOTION_OUTPUT;
   if( m_separate_motions)
      rval |= OPTION_SEPARATE_MOTIONS;
   if( m_round_step)
      rval |= OPTION_ROUND_TO_NEAREST_STEP;
   if( m_radial_velocity)
      rval |= OPTION_RADIAL_VEL_OUTPUT;
   if( m_phase_angle)
      rval |= OPTION_PHASE_ANGLE_OUTPUT;
   if( m_phase_angle_bisector)
      rval |= OPTION_PHASE_ANGLE_BISECTOR;
   if( m_helio_ecliptic)
      rval |= OPTION_HELIO_ECLIPTIC;
   if( m_topo_ecliptic)
      rval |= OPTION_TOPO_ECLIPTIC;
   if( m_suppress_unobservable)
      rval |= OPTION_SUPPRESS_UNOBSERVABLE;
   if( m_visibility)
      rval |= OPTION_VISIBILITY;
   if( m_show_sigmas)
      rval |= OPTION_SHOW_SIGMAS;
   if( !m_human_readable)
      rval |= OPTION_COMPUTER_FRIENDLY;
   if( m_ground_track)
      rval |= OPTION_GROUND_TRACK;
   if( m_speed)
      rval |= OPTION_SPACE_VEL_OUTPUT;
   return( rval);
}

void CEphem::OnClickedGo()
{
   char buff[200], *err_msg = NULL;
   // TODO: Add your control notification handler code here
   set_jd_from_xtrols( buff);
   extern double minimum_jd, maximum_jd;

   const double step_size = get_step_size( CT2A( m_ephem_step), NULL, NULL);
// CT2A mephem_step( m_ephem_step);
// const double step_size = get_step_size( mephem_step.m_psz, NULL, NULL);

   if( *buff)
      err_msg = buff;
   else if( !step_size)
      {
      strcpy( buff, get_find_orb_text( 7));   /* "No step size specified!" */
      err_msg = buff;
      }
   else if( !m_number_steps)
      {
      strcpy( buff, get_find_orb_text( 6));
      err_msg = buff;   /* "Ephemeris must contain at least one entry!" */
      }
   else if( jd < minimum_jd || jd > maximum_jd)
      {
      strcpy( buff, get_find_orb_text( 14));
      err_msg = buff;   /* "Invalid date!" */
      }

   if( !err_msg)
      {
      double rho_sin_phi, rho_cos_phi, lon = 0.;
      char note_text[80];
      int options = GetEphemerisBitmask( );
      int planet_no;

      if( m_use_mpc_code)
         {
#ifdef _UNICODE
         CT2A mpc_code( m_mpc_code);
         planet_no = get_observer_data( mpc_code.m_psz, buff, &lon,
                                   &rho_cos_phi, &rho_sin_phi);
#else
         planet_no = get_observer_data( m_mpc_code, buff, &lon,
                                   &rho_cos_phi, &rho_sin_phi);
#endif
         if( planet_no < 0)
            MessageBox( CA2T( err_msg),
                    _T( "MPC code unknown;  using lat/lon instead\n"), MB_OK);
         }
      else
         planet_no = -1;
      *note_text = '\0';
      if( planet_no < 0)
         {
         planet_no = 3;
         double lat;
#ifdef _UNICODE
         CT2A mlat( m_lat);
         CT2A mlon( m_lon);
#endif

         if( *m_lat == 'n' || *m_lat == 'N')
            lat = _tcstod( (const TCHAR *)m_lat + 1, NULL);
//          lat = atof( (const char *)m_lat + 1);
         else if( *m_lat == 's' || *m_lat == 'S')
            lat = -_tcstod( (const TCHAR *)m_lat + 1, NULL);
         else
            {
            strcpy( buff, get_find_orb_text( 12));
            err_msg = buff;    /* Latitude must start with an 'N' or 'S'! */
            lat = 0.;
            }

         if( !err_msg)
            if( *m_lon == 'e' || *m_lon == 'E')
               lon = _tcstod( (const TCHAR *)m_lon + 1, NULL);
            else if( *m_lon == 'w' || *m_lon == 'W')
               lon = -_tcstod( (const TCHAR *)m_lon + 1, NULL);
            else
               {
               strcpy( buff, get_find_orb_text( 13));
               err_msg = buff; /* Longitude must start with an 'E' or 'W'! */
               }
         lat *= PI / 180.;
         lon *= PI / 180.;
         lat_alt_to_parallax( lat, 0., &rho_cos_phi, &rho_sin_phi, 3);
#ifdef _UNICODE
         sprintf( note_text, "For %s, %s", mlon.m_psz, mlat.m_psz);
#else
         sprintf( note_text, "For %s, %s", m_lon, m_lat);
#endif
         }
      else
         strcpy( note_text, mpc_station_name( buff)); /* copy in observer loc */
      if( !err_msg)
         {
         double temp_orbit[12];
         unsigned n_orbits = 1;
         extern int available_sigmas;
         extern double ephemeris_mag_limit;

         memcpy( temp_orbit, orbit, 6 * sizeof( double));
         if( available_sigmas == COVARIANCE_AVAILABLE
                     && (options & 7) == OPTION_OBSERVABLES)
            {
            compute_variant_orbit( temp_orbit + 6, temp_orbit, 1.);
            n_orbits = 2;
            }
         ephemeris_mag_limit = m_mag_limit;
         ephemeris_in_a_file( ephemeris_filename, temp_orbit, obs, n_obs,
                     planet_no, epoch, jd,
                     CT2A( m_ephem_step), lon, rho_cos_phi, rho_sin_phi,
                     m_number_steps, note_text,
                     (ephem_option_t)options, n_orbits);

         CListBox* pListBox = (CListBox*)GetDlgItem( IDC_LIST1);
         FILE *ifile = fopen( ephemeris_filename, "r");

         pListBox->ResetContent( );
         while( fgets_trimmed( buff, sizeof( buff), ifile))
            if( *buff != '#' && *buff)
               {
               remove_rgb_code( buff);
               pListBox->AddString( CA2T( buff));
               }
         fclose( ifile);
         make_pseudo_mpec( mpec_filename, obj_name);      /* ephem0.cpp */
         ephemeris_and_pseudo_mpec_made = 1;
         if( options & (OPTION_STATE_VECTOR_OUTPUT | OPTION_POSITION_OUTPUT))
            {
            FILE *ofile = fopen( ephemeris_filename, "a");

            add_ephemeris_details( ofile, jd,
                          jd + step_size * (double)m_number_steps);
            fclose( ofile);
            }
         }
      }
   if( err_msg)
      MessageBox( CA2T( err_msg), _T( "Ephemeris"), MB_OK);
}

void load_font_from_text( LOGFONT *lf, const char *font_specifier)
{
   memset( lf, 0, sizeof(LOGFONT));        // Clear out structure.
   if( font_specifier && *font_specifier)
      {
      int bytes, fields[8];

      sscanf( font_specifier, "%ld %ld %ld %ld %ld %d %d %d %d %d %d %d %d %n",
            &lf->lfHeight, &lf->lfWidth,
            &lf->lfEscapement, &lf->lfOrientation, &lf->lfWeight,
            fields, fields + 1, fields + 2, fields + 3, fields + 4,
            fields + 5, fields + 6, fields + 7, &bytes);
      strcpy( (char *)lf->lfFaceName, font_specifier + bytes);
      lf->lfItalic          = (BYTE)fields[0];
      lf->lfUnderline       = (BYTE)fields[1];
      lf->lfStrikeOut       = (BYTE)fields[2];
      lf->lfCharSet         = (BYTE)fields[3];
      lf->lfOutPrecision    = (BYTE)fields[4];
      lf->lfClipPrecision   = (BYTE)fields[5];
      lf->lfQuality         = (BYTE)fields[6];
      lf->lfPitchAndFamily  = (BYTE)fields[7];
      debug_printf( "Face name '%s'\n", lf->lfFaceName);
      }
   else
      {
      lf->lfHeight = -12;                      // Request a 12-pixel-high font
#ifdef _UNICODE
      wcscpy(lf->lfFaceName, L"Courier New");       // Request font
#else
      strcpy(lf->lfFaceName, "Courier New");       // Request font
#endif
      lf->lfPitchAndFamily = FIXED_PITCH | FF_MODERN;
      lf->lfWeight = FW_NORMAL;
      lf->lfCharSet = DEFAULT_CHARSET;
      }
}

BOOL CEphem::OnInitDialog()
{
   CDialog::OnInitDialog();

   // TODO: Add extra initialization here
   reset_dialog_language( this, 97000);
   GetWindowRect( &OriginalDlgRect);
   ScreenToClient( &OriginalDlgRect);
   set_window_placement( this->m_hWnd,
                  get_environment_ptr( "EPHEM_WINDOW"));

   LOGFONT lf;                             // Used to create the CFont.
   const char *font_str = get_environment_ptr( "EPHEM_FONT");
   CListBox* pListBox = (CListBox*)GetDlgItem( IDC_LIST1);

   load_font_from_text( &lf, font_str);
   list_box_font.CreateFontIndirect(&lf);    // Create the font.
   pListBox->SetFont( &list_box_font);     // set the font
   pListBox->SetItemHeight( 0, 1 - lf.lfHeight);
   return TRUE;  // return TRUE  unless you set the focus to a control
}

void CEphem::OnPseudoMpec()
{
   // TODO: Add your control notification handler code here

   if( !ephemeris_and_pseudo_mpec_made)
      MessageBox( _T( "You must make an ephemeris before\nmaking a pseudo-MPEC."),
                        _T( "Find_Orb"), MB_OK);
   else
      if( (int)ShellExecuteA( NULL, NULL, mpec_filename, NULL,
                                "", SW_SHOWNORMAL) <= 32)
         ShellExecuteA( NULL, "explore", mpec_filename, NULL, "", SW_SHOWNORMAL);
#ifdef OBSOLETE_VERSION
      {
      char filename[_MAX_DIR];

      get_file_from_dialog( FALSE, "", "*.*", filename, NULL);
      if( *filename)
         {
         unlink( filename);
         rename( mpec_filename, filename);
         }
      }
#endif
}

int copy_file_to_clipboard( const char *filename);    /* ephem0.cpp */
int clipboard_to_file( const char *filename, const int append);      /* mpc_obs.cpp */

void CEphem::OnCopy()
{
   // TODO: Add your control notification handler code here
   if( !ephemeris_and_pseudo_mpec_made)
      MessageBox( _T( "You must make an ephemeris before\ncopying to the clipboard."),
                  _T( "Find_Orb"), MB_OK);
   else
      {
      int rval;
      const char *temp_name = "ephem.tmp";

      save_ephemeris_file( temp_name);
      rval = copy_file_to_clipboard(  temp_name);
      if( rval)
         {
         char buff[80];

         sprintf( buff, "rval %d", rval);
         MessageBox( CA2T( buff), _T( "Find_Orb"), MB_OK);
         }
      _unlink( temp_name);
      }
}

void CEphem::OnSize(UINT nType, int cx, int cy)
{
   CDialog::OnSize(nType, cx, cy);

   // TODO: Add your message handler code here
   CListBox* pListBox = (CListBox*)GetDlgItem( IDC_LIST1);

   if( pListBox)
      {
      CRect rect;

      pListBox->GetWindowRect( &rect);
      ScreenToClient( &rect);
      pListBox->SetWindowPos( this, 0, 0,
                  cx,
                  cy - rect.top,
                  SWP_NOMOVE | SWP_NOZORDER);
      }

   // TODO: Add your message handler code here
}

/* Solution found from http://www.codeguru.com/forum/showthread.php?t=318933 */

void CEphem::OnGetMinMaxInfo(MINMAXINFO FAR* lpMMI)
{
  // set the minimum tracking width
  // and the minimum tracking height of the window
  if( OriginalDlgRect.Height( ))
     {
     lpMMI->ptMinTrackSize.x = OriginalDlgRect.Width( );
//   lpMMI->ptMaxTrackSize.x = OriginalDlgRect.Width( );
     lpMMI->ptMinTrackSize.y = OriginalDlgRect.Height( );
     }
}

void CEphem::OnDestroy()
{
   char *place_text = get_placement_text( this->m_hWnd);
   set_environment_ptr( "EPHEM_WINDOW", place_text);
   free( place_text);

   CDialog::OnDestroy();
   // TODO: Add your message handler code here
}

void CEphem::OnNow()
{
   // TODO: Add your control notification handler code here
   m_day = "+0";
   UpdateData( FALSE);     /* Move change to edit box */
}
