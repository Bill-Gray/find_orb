/* orbitdlg.cpp: main dialog for Windows Find_Orb

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

// orbitdlg.cpp : implementation file
//

#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <direct.h>
#include <assert.h>
#include <dos.h>
#include <sys/stat.h>
#include "watdefs.h"
#include "find_orb.h"
#include "mpc_obs.h"
#include "orbitdlg.h"
#include "afuncs.h"
#include "comets.h"
#include "ephem.h"
#include "about.h"
#include "date.h"
#include "generic.h"
#include "settings.h"
#include "sigma.h"
#include "monte0.h"

#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

#ifndef CA2T
   #define CT2A(x) ((const char *)x)
   #define CA2T(x, encoding) ((const char *)x)
#endif

extern double solar_pressure[];
extern int n_extra_params;
extern int prev_shifted_residual = -1;
extern unsigned perturbers;
static const TCHAR *program_name = _T( "Find_Orb");

#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078
#define AU_IN_LIGHT_YEAR ((365.25 * seconds_per_day * SPEED_OF_LIGHT) / AU_IN_KM)

#define SRP1AU 2.3e-7
             /* "Solar radiation pressure at 1 AU",  in             */
             /* kg*AU^3 / (m^2*d^2),  from a private communication  */
             /* from Steve Chesley; see orb_func.cpp for details    */

#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_MINOR_AXIS 6356755.
#define EARTH_AXIS_RATIO (EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS)

#ifndef _UNICODE
void utf8_to_win1252( char *text);                    /* elem_out.cpp */
#endif
int format_jpl_ephemeris_info( char *buff);           /* pl_cache.cpp */
int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit);
int find_trial_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
                 const double r1, const double angle_param);
void improve_parabolic( OBSERVE FAR *obs, int n_obs, double *orbit, double epoch);
int set_locs( const double *orbit, double t0, OBSERVE FAR *obs, int n_obs);
double initial_orbit( OBSERVE FAR *obs, int n_obs, double *orbit);
int get_r1_and_r2( const int n_obs, const OBSERVE FAR *obs,
                             double *r1, double *r2);    /* orb_func.cpp */
int write_residuals_to_file( const char *filename, const char *ast_filename,
          const int n_obs, const OBSERVE FAR *obs_data, const int short_form);
                                                /* ephem0.cpp */
const char *get_environment_ptr( const char *env_ptr);     /* mpc_obs.cpp */
void set_environment_ptr( const char *env_ptr, const char *new_value);
double convenient_gauss( const OBSERVE FAR *obs, int n_obs, double *orbit,
                  const double mu, const int desired_soln); /* gauss.cpp */
int debug_printf( const char *format, ...);                /* runge.cpp */
char *fgets_trimmed( char *buff, size_t max_bytes, FILE *ifile); /*elem_out.c*/
int get_jpl_ephemeris_info( int *de_version, double *jd_start, double *jd_end);
void set_statistical_ranging( const int new_using_sr);      /* elem_out.cpp */
int find_nth_sr_orbit( double *orbit, OBSERVE FAR *obs, int n_obs,
                            const int orbit_number);       /* orb_func.cpp */
void set_window_placement( HWND hwnd, const char *ibuff);  /* orbitdlg.cpp */
char *get_placement_text( HWND hwnd);                      /* orbitdlg.cpp */
int store_defaults( const int ephemeris_output_options,
         const int element_format, const int element_precision,
         const double max_residual_for_filtering,
         const double noise_in_arcseconds);           /* elem_out.cpp */
int get_defaults( int *ephemeris_output_options, int *element_format,
         int *element_precision, double *max_residual_for_filtering,
         double *noise_in_arcseconds);                /* elem_out.cpp */
const char *get_find_orb_text( const int index);      /* elem_out.cpp */
void get_find_orb_text_filename( char *filename);     /* elem_out.cpp */
int reset_dialog_language( CDialog *dlg, const int dlg_number);  /* elem_out.cpp */

/////////////////////////////////////////////////////////////////////////////
// COrbitDlg dialog

COrbitDlg::COrbitDlg(CWnd* pParent /*=NULL*/)
   : CDialog(COrbitDlg::IDD, pParent)
{
   extern char default_comet_magnitude_type;

   //{{AFX_DATA_INIT(COrbitDlg)
   m_step_size = 0;
   m_epoch = "";
   m_r1 = "1.";
   m_r2 = "1.";
   //}}AFX_DATA_INIT
   n_obs = monte_carlo = 0;
   compute_covariance = 1;
   element_format = 0;
   element_precision = 5;
   show_commented_elements = 0;
   obs_data = NULL;
   max_residual_for_filtering = 1.;
   OriginalDlgRect.top = OriginalDlgRect.bottom = 0;
   constraints = "";
   monte_noise = .5;
   n_objects = 0;
   precise_residuals = 0;
   obj_info = NULL;
   ephemeris_output_options = 0;
   get_defaults( &ephemeris_output_options, &element_format,
         &element_precision, &max_residual_for_filtering,
         &monte_noise);                               /* elem_out.cpp */

}

void COrbitDlg::DoDataExchange(CDataExchange* pDX)
{
   CDialog::DoDataExchange(pDX);
   //{{AFX_DATA_MAP(COrbitDlg)
   DDX_Text(pDX, IDC_EPOCH, m_epoch);
   DDV_MaxChars(pDX, m_epoch, 17);
   DDX_Text(pDX, IDC_R1, m_r1);
   DDX_Text(pDX, IDC_R2, m_r2);
   //}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(COrbitDlg, CDialog)
   //{{AFX_MSG_MAP(COrbitDlg)
   ON_BN_CLICKED(IDC_FULL_STEP, OnClickedFullStep)
   ON_BN_CLICKED(IDC_HERGET, OnClickedHerget)
   ON_BN_CLICKED(IDC_OPEN, OnClickedOpen)
   ON_LBN_DBLCLK(IDC_LIST_ASTEROIDS, OnDblclkObject)
   ON_BN_CLICKED(IDC_SAVE, OnClickedSave)
   ON_LBN_SELCHANGE(IDC_RESIDUALS, OnSelchangeResiduals)
   ON_BN_CLICKED(IDC_MAKE_EPHEMERIS, OnClickedMakeEphemeris)
   ON_BN_CLICKED(IDC_SAVE_RESIDS, OnClickedSaveResids)
   ON_LBN_DBLCLK(IDC_RESIDUALS, OnDblclkResiduals)
   ON_BN_CLICKED(IDC_ABOUT, OnClickedAbout)
   ON_BN_CLICKED(IDC_VAISALA, OnClickedVaisala)
   ON_BN_CLICKED(IDC_AUTO_SOLVE, OnClickedAutoSolve)
   ON_WM_CHAR()
   ON_BN_CLICKED(IDC_MONTE_CARLO, OnMonteCarlo)
   ON_WM_TIMER()
   ON_BN_CLICKED(IDC_GAUSS, OnGauss)
   ON_BN_CLICKED(IDC_WORST, OnWorst)
   ON_BN_DOUBLECLICKED(IDC_WORST, OnDoubleclickedWorst)
   ON_BN_CLICKED(IDC_FILTER_OBS, OnFilterObs)
   ON_LBN_SELCANCEL(IDC_LIST_ASTEROIDS, OnSelcancelListAsteroids)
   ON_BN_CLICKED(IDC_ASTEROIDS, OnAsteroids)
   ON_BN_CLICKED(IDC_SETTINGS, OnSettings)
   ON_LBN_SELCHANGE(IDC_LIST_ASTEROIDS, OnSelchangeListAsteroids)
   ON_BN_CLICKED(IDC_SET_WEIGHT, OnSetWeight)
   ON_BN_CLICKED(IDC_TOGGLE_OBS, OnToggleObs)
   ON_BN_CLICKED(IDC_ORBITAL_ELEMENTS, OnOrbitalElements)
   ON_WM_DESTROY()
   ON_WM_RBUTTONUP()
   ON_WM_SIZE()
   ON_BN_CLICKED(IDC_TOGGLE_PERTURBERS, OnTogglePerturbers)
   ON_BN_CLICKED(IDC_SIMPLEX, OnSimplex)
   ON_BN_CLICKED(IDC_STATISTICAL_RANGING, OnStatisticalRanging)
   ON_WM_CTLCOLOR()
   ON_EN_KILLFOCUS(IDC_EPOCH, OnKillfocusEpoch)
   ON_BN_DOUBLECLICKED(IDC_WORST, OnDoubleclickedWorst)
   ON_WM_GETMINMAXINFO()
   ON_WM_LBUTTONUP()
   //}}AFX_MSG_MAP
END_MESSAGE_MAP()

void get_file_from_dialog( int is_open, const char *default_ext,
                           const char *filter, char *buff, const char *path)
{
   char old_path[_MAX_DIR];
#ifndef _WIN32
   unsigned int n_drives;
#endif

   *buff = '\0';
   _getcwd( old_path, _MAX_DIR);
   if( path && *path)
      {
      size_t i;
      char path2[_MAX_DIR];

      for( i = strlen( path); i && path[i - 1] != '\\'; i--)
         ;
      strcpy( path2, path);
      path2[i] = '\0';
      _chdir( path2);
      }
   CFileDialog dlg( is_open, CA2T( default_ext, CP_UTF8), CA2T( filter, CP_UTF8));
   if( dlg.DoModal( ) == IDOK)
      strcpy( buff, CT2A( dlg.GetPathName( )));

   _chdir( old_path);
#ifndef _WIN32
   _dos_setdrive( *old_path - 'A' + 1, &n_drives);
#endif
}

static double extract_epoch( const char *epoch_text)
{
   const double jan_1970 = 2440587.5;
   const double initial_jd =
            jan_1970 + (double)(time( NULL) / seconds_per_day) + .5;
   const double rval = get_time_from_string( initial_jd,
                   epoch_text, FULL_CTIME_YMD | CALENDAR_JULIAN_GREGORIAN, NULL);

   return( rval > 1. ? rval : initial_jd);
}

static void format_distance( char *buff, const double ival)
{
   if( ival < .001)       /* about 150,000 km */
      sprintf( buff, "%dkm", (long)( ival * AU_IN_KM));
   else if( ival < 9999.)
      sprintf( buff, "%.4lf", ival);
   else if( ival < 9999 * AU_IN_LIGHT_YEAR)
      sprintf( buff, "%.4lf LY", ival / AU_IN_LIGHT_YEAR);
   else
      strcpy( buff, "<HUGE>");
}

void COrbitDlg::Reset_r1_and_r2( void)
{
   double r1, r2;
   char buff[80];

   get_r1_and_r2( n_obs, (OBSERVE FAR *)obs_data, &r1, &r2);    /* orb_func.cpp */
   format_distance( buff, r1);
   m_r1 = buff;
   if( m_r2[0] < 'A' || m_r2[0] > 'z')
      {
      format_distance( buff, r2);
      m_r2 = buff;
      }
   UpdateData( FALSE);  /* 'False' indicates 'move data to xtrols' */
}

static double parse_distance_text( const char *buff)
{
   double rval = atof( buff);
   const char end_char = buff[strlen( buff) - 1];

   if( end_char == 'm')
      rval /= AU_IN_KM;
   else if( end_char == 'Y')
      rval *= AU_IN_LIGHT_YEAR;
   return( rval);
}

int get_idx1_and_idx2( const int n_obs, const OBSERVE FAR *obs,
                                int *idx1, int *idx2);      /* elem_out.c */
TCHAR *tget_find_orb_text( const int index)      /* elem_out.cpp */
{
   static TCHAR buff[100];

   _tcscpy( buff, CA2T( get_find_orb_text( index), CP_UTF8));
   return( buff);
}

int COrbitDlg::ImproveOrbitSolution( int full_step, int n_repeats)
{
   int is_vaisala = (full_step & 2);
   int is_linearizing = (full_step & 4);
   int rval = 0;

   full_step &= 1;
   UpdateData( TRUE);     /* Get changes made in edit boxes */
   if( n_obs && obs_data)
      {
      OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;
      const char *limited_orbit = NULL;
      char tbuff[80];
      int i, idx1, idx2;

      strcpy( tbuff, CT2A( m_r1));
      for( i = 0; tbuff[i] && tbuff[i] != '='; i++)
         ;
      if( tbuff[i] == '=')
         {
         tbuff[i] = '\0';
         set_environment_ptr( tbuff, CT2A( m_r1) + i + 1);
         }

      if( m_r2[0] == 's')
         n_extra_params = atoi( CT2A( m_r2) + 1);
      else if( m_r2[0] == 'k')
         {
         extern int integration_method;

         integration_method = atoi( CT2A( m_r2) + 1);
         }
      else if( m_r2[0] == 't')
         {
         extern double integration_tolerance;

         integration_tolerance = atof( CT2A( m_r2) + 1);
         }
#ifdef OBSOLETE_DEBUGGING_CODE
      else if( m_r2[0] == 'r')
         {
         extern double relativistic_factor;

         relativistic_factor = atof( CT2A( m_r2) + 1);
         }
#endif
      else if( m_r2[0] == 'j')
         {
         extern double j2_multiplier;

         j2_multiplier = atof( CT2A( m_r2) + 1);
         }
      else if( m_r2[0] >= 'A' && m_r2[0] <= 'z')
         limited_orbit = CT2A( m_r2);
      SetCursor( AfxGetApp( )->LoadStandardCursor( IDC_WAIT));
      GetPerturberMask( );
      get_idx1_and_idx2( n_obs, obs, &idx1, &idx2);

      if( !full_step)     /* Vaisala or Herget orbit */
         {
         if( is_linearizing)
            adjust_herget_results( obs, n_obs, orbit);
         else if( m_r1[0] != '<' && m_r2[0] != '<')
            {
            double r1 = parse_distance_text( CT2A( m_r1));
            double r2 = parse_distance_text( CT2A( m_r2));
            double d_r1, d_r2;

            if( limited_orbit)
               {
               for( i = n_obs - 1; i && !obs[i].is_included; i--)
                  ;
               r2 = obs[i].r;
               }
            if( is_vaisala)
               {
               double angle_param;

               if( _stscanf( m_r1, _T( "%lf,%lf"), &r1, &angle_param) == 2)
                  {
                  if( !find_trial_orbit( orbit, obs, n_obs, r1, angle_param))
                     {
                     r1 = obs->r;
                     r2 = obs[n_obs - 1].r;
                     }
                  }
               else
                  {
                  herget_method( obs, n_obs, -r1, r2, orbit, &d_r1, &d_r2,
                                            CT2A( constraints));
                  r1 = d_r1;
                  r2 = d_r2;
                  }
               }
            else
               {
               herget_method( obs, n_obs, r1, r2, orbit, &d_r1, &d_r2,
                                            CT2A( constraints));
               herget_method( obs, n_obs, r1 + d_r1, r2 + d_r2, orbit,
                                          NULL, NULL, NULL);
               }
            }
         orbit_epoch = obs[idx1].jd;
         }
      else
         {
         while( !rval && (GetAsyncKeyState( VK_SHIFT) & 0x8001) == 0x8001)
            rval = full_improvement( obs, n_obs, orbit, orbit_epoch,
                                    CT2A( constraints),
                                    NO_ORBIT_SIGMAS_REQUESTED, orbit_epoch);
         while( !rval && n_repeats--)
            {
            int sigma_request = NO_ORBIT_SIGMAS_REQUESTED;
            const double epoch_in_edit_box = extract_epoch( CT2A( m_epoch));

                          /* on last run,  request correct sigmas */
            if( !n_repeats && compute_covariance)
               sigma_request = ((element_format & ELEM_OUT_HELIOCENTRIC_ONLY) ?
                        HELIOCENTRIC_SIGMAS_ONLY : ORBIT_SIGMAS_REQUESTED);
            rval = full_improvement( obs, n_obs, orbit, orbit_epoch,
                                    CT2A( constraints),
                                    sigma_request, epoch_in_edit_box);
            }
         }

      Reset_r1_and_r2( );
//    if( !monte_carlo)
         {
         UpdateElementDisplay( 1);
         UpdateResidualDisplay( );
         }
      SetCursor( AfxGetApp( )->LoadStandardCursor( IDC_ARROW));
      }
   else                    /* "No orbit to improve!" */
      MessageBox( tget_find_orb_text( 1), _T( "Find_Orb"), MB_OK);
   return( rval);
}

/////////////////////////////////////////////////////////////////////////////
// COrbitDlg message handlers

void COrbitDlg::OnClickedFullStep()
{
   // TODO: Add your control notification handler code here

   ImproveOrbitSolution( 1, 1);
   if( *get_environment_ptr( "LINEAR_FULL"))
      ImproveOrbitSolution( 4, 1);     /* ...then linearize */
}

void COrbitDlg::OnClickedHerget()
{
   // TODO: Add your control notification handler code here

   ImproveOrbitSolution( 0, 1);
   if( *get_environment_ptr( "LINEAR_HERGET"))
      ImproveOrbitSolution( 4, 1);     /* ...then linearize */
}

void COrbitDlg::ResetPerturbers()
{
   int i;

   for( i = 0; i < 11; i++)
      {
      int is_on = ((perturbers >> (i + 1)) & 1);
      CButton *pButton = (CButton *)GetDlgItem( i + IDC_PLANET1);

      if( i == 10)            /* asteroids are a special case */
         is_on = (perturbers >> 20) & 7;
      pButton->SetCheck( is_on);           /* turn off all the perturbers */
      }
   SetDlgItemText( IDC_TOGGLE_PERTURBERS,
                           tget_find_orb_text( perturbers > 1 ? 22 : 21));
}

int COrbitDlg::GetPerturberMask()
{
   int i, rval = 0;

   for( i = IDC_PLANET1; i <= IDC_PLANET10; i++)
      if( ((CButton *)GetDlgItem( i))->GetCheck( ))
         rval |= (2 << (i - IDC_PLANET1));
   if( ((CButton *)GetDlgItem( IDC_ASTEROIDS))->GetCheck( ))
      rval |= (7 << 20);
   perturbers = rval;
   SetDlgItemText( IDC_TOGGLE_PERTURBERS,
                           tget_find_orb_text( perturbers > 1 ? 22 : 21));
   return( rval);
}

void COrbitDlg::LoadAFile( const char *filename)
{
   CListBox* pListBox = (CListBox*)GetDlgItem( IDC_LIST_ASTEROIDS);

   if( obj_info)
      free( obj_info);
   obj_info = find_objects_in_file( filename, &n_objects, NULL);
   if( obj_info && n_objects && pListBox)
      {
      int i;
      struct stat s;

      if( !stat( filename, &s))
         curr_file_time = s.st_mtime;
      pListBox->ResetContent( );
      curr_file_name = filename;
      curr_object_name = "";

      for( i = 0; i < n_objects; i++)
         pListBox->AddString( CA2T( obj_info[i].obj_name, CP_UTF8));
      pListBox->SetSel( 1, TRUE);
      if( n_objects == 1)
         LoadAnObject( 0);
      }
   else
      MessageBox( tget_find_orb_text( obj_info ? 10 : 11), program_name, MB_OK);
}

void COrbitDlg::OnClickedOpen()
{
   // TODO: Add your control notification handler code here
   char filename[_MAX_DIR];
   static char path[_MAX_DIR];

   get_file_from_dialog( TRUE, "", "*.*", filename, path);
   if( *filename)
      {
      LoadAFile( filename);
      strcpy( path, filename);
      }
}

void add_version_and_de_text( char *buff);

int reset_dialog_language( CDialog *dlg, const int dlg_number)
{
   char buff[400];
   int in_dialog = 0;
   FILE *ifile;

   get_find_orb_text_filename( buff);             /* elem_out.cpp */
   ifile = fopen( buff, "rb");
   assert( ifile);
   while( fgets_trimmed( buff, sizeof( buff), ifile))
      {
      const int number = atoi( buff) - dlg_number;

      if( number >= -1 && number < 999)
         {
#ifndef _UNICODE
         utf8_to_win1252( buff);
#endif
         if( !number)
            dlg->SetWindowText( CA2T( buff + 8, CP_UTF8));
         else
            dlg->SetDlgItemText( number, CA2T( buff + 8, CP_UTF8));
         if( number == IDC_VERSION_INFO)
            {                                /* add JPL DE info */
            int dlg_id;

            add_version_and_de_text( buff + 8);
            if( dlg_number == 99000)     /* main dialog */
               dlg_id = IDC_STATION_INFO;
            else                         /* "about" box */
               dlg_id = IDC_VERSION_INFO;
            dlg->SetDlgItemText( dlg_id, CA2T( buff + 8, CP_UTF8));
            }
         }
      }
   fclose( ifile);
   return( 0);
}

void set_window_placement( HWND hwnd, const char *ibuff)
{
   if( ibuff && *ibuff)
      {
      WINDOWPLACEMENT place;
      char *tptr = (char *)&place;
      int tval, i;

      for( i = 0; i < sizeof( WINDOWPLACEMENT); i++)
         {
         sscanf( ibuff + i * 3, "%x", &tval);
         *tptr++ = (char)tval;
         }
      SetWindowPlacement( hwnd, &place);
      }
}

void load_font_from_text( LOGFONT *lf, const char *font_specifier);
         /* see 'ephem.cpp' */

BOOL COrbitDlg::OnInitDialog()
{
   CDialog::OnInitDialog();

   // TODO: Add extra initialization here
   FILE *startup;
   const char FAR *cmd_line = CT2A( AfxGetApp( )->m_lpCmdLine);
   const char FAR *cmd_language = NULL;
   int i;
   extern char findorb_language;       /* defaults to 'e' for English */

   GetWindowRect( &OriginalDlgRect);
   ScreenToClient( &OriginalDlgRect);
   set_window_placement( this->m_hWnd,
                  get_environment_ptr( "WINDOW_PLACEMENT"));

   for( i = 0; cmd_line[i] && !cmd_language; i++)
      if( cmd_line[i] == '-' && cmd_line[i + 1] == 'l')
         cmd_language = cmd_line + i;
   if( cmd_language)
      findorb_language = cmd_language[2];
   else if( startup = fopen( "startup.mar", "rb"))
      {
      char buff[140];

      while( fgets( buff, 140, startup))
         if( !memcmp( buff, "51 language", 11))
            findorb_language = buff[12];
      fclose( startup);
      }
   reset_dialog_language( this, 99000);
   perturbers = 0;
   m_step_size = 3.;
   srand( (unsigned)time( NULL));         /* for Monte Carlo code */
   load_up_sigma_records( "sigma.txt");
   m_r1 = "1";
   m_r2 = "1";
   UpdateData( FALSE);     /* 'False' indicates 'move data to edit boxes' */
   ResetPerturbers( );
   SetTimer( 1, 500, NULL);
   SetDlgItemText( IDC_ORBIT1, CA2T(
               "Orbital elements will appear here after you open a file,\n"
               "select an object from it,  and an orbit is computed.\n\n",
               CP_UTF8));
   if( *cmd_line != '-' && *cmd_line)
      LoadAFile( cmd_line);
   AdjustControls( );
   return TRUE;  // return TRUE  unless you set the focus to a control
}

extern const char *elements_filename;

void COrbitDlg::UpdateElementDisplay( int update_orbit)
{
   int i, index, obs_format, n_selected, *selections;
   char *obuff = (char *)calloc( 12, 80);
   OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;
   CListBox* pListBox = (CListBox*)GetDlgItem( IDC_RESIDUALS);
// int residual_stops[14] = {  0,  10,  20, 53, 60,  80,  90, 102,
//                           130, 145, 155, 180, 210, 270 };
/*                            YY   MM   DD. X  mpc   HH   MM   SS.SSS   */
   int residual_stops[14] = {  0,  10,  20, 58, 65,  85,  95, 107,
                             140, 155, 165, 195, 225, 285 };
/*                           ddd   mm  ss   dx   dy                     */

   if( n_obs >= 2 && update_orbit)
      {
      double epoch_in_edit_box;
      FILE *ifile;
      char *orbit_buff = (char *)malloc( 2000), *tptr = orbit_buff;

      assert( orbit_buff);
      if( !orbit_buff)
         return;
      UpdateData( TRUE);     /* Get changes made in edit boxes */
      epoch_in_edit_box = extract_epoch( CT2A( m_epoch));
      bad_elements = write_out_elements_to_file( orbit,
                     orbit_epoch, epoch_in_edit_box,
                     (OBSERVE FAR *)obs, n_obs,
                     CT2A( constraints), element_precision, monte_carlo,
                     element_format);
      ifile = fopen( elements_filename, "rb");
      fgets_trimmed( tptr, 200, ifile);           /* "Orbital elements:" */
      SetDlgItemText( IDC_ORBITAL_ELEMENTS, CA2T( tptr, CP_UTF8));
      i = 0;
      while( i < 11 && fgets( tptr, 200, ifile))
         {
         bool show_this_line;

         if( show_commented_elements)
            {
            show_this_line = (*tptr == '#' && tptr[3] != '$'
                   && memcmp( tptr, "# Find", 6)
                   && memcmp( tptr, "# Scor", 6));
            if( show_this_line)
               memmove( tptr, tptr + 2, strlen( tptr + 1));
            }
         else
            show_this_line = (*tptr != '#');
         if( show_this_line)
            {
            tptr += strlen( tptr);
            i++;
            }
         }
      *tptr = '\0';

      SetDlgItemText( IDC_ORBIT1, CA2T( orbit_buff, CP_UTF8));

      fclose( ifile);
      if( monte_carlo)
         {
         extern int monte_carlo_object_count;
         extern int using_sr;
         char buff[10];

         sprintf( buff, "%d", monte_carlo_object_count);
         SetDlgItemText(
               (using_sr ? IDC_STATISTICAL_RANGING : IDC_MONTE_CARLO),
                                        CA2T( buff, CP_UTF8));
         }
      }

   index = pListBox->GetTopIndex( );
   n_selected = pListBox->GetSelCount( );
   selections = (int *)calloc( n_selected, sizeof( int));
   assert( selections);    /* store existing selections,  rebuild list box, */
   if( !selections)        /* & restore the selected items  */
      return;
   pListBox->GetSelItems( n_selected, selections);
   pListBox->ResetContent( );
   pListBox->SetTabStops( 14, (LPINT)residual_stops);
   obs_format = RESIDUAL_FORMAT_FULL_WITH_TABS;
   if( precise_residuals)
      obs_format |= RESIDUAL_FORMAT_PRECISE;

   for( i = 0; i < n_obs; i++)
      {
//    if( !obs_format)
//       {
//       recreate_observation_line( obuff, obs + i);
//       memmove( obuff, obuff + 12, strlen( obuff + 11));
//       }
//    else
         format_observation( obs + i, obuff, obs_format);
      pListBox->AddString( CA2T( obuff, CP_UTF8));
      }

   for( i = 0; i < n_selected; i++)
      pListBox->SetSel( selections[i], TRUE);
   free( selections);
   if( index >= 0 && index < n_obs)
      pListBox->SetTopIndex( index);
   free( obuff);
}

static void put_epoch_text( char *buff, const double jd)
{
   size_t i;

   full_ctime( buff, jd, FULL_CTIME_YMD | FULL_CTIME_DATE_ONLY | 0x40
                         | CALENDAR_JULIAN_GREGORIAN);
   for( i = strlen( buff); buff[i - 1] == '0'; i--)
      ;                    /* trim trailing zeroes */
   if( buff[i - 1] == '.')
      i--;
   buff[i] = '\0';
}

void COrbitDlg::LoadAnObject( const int obj_idx)
{
   OBSERVE FAR *obs;
   char buff[90];
   long file_offset;
   FILE *ifile = fopen( CT2A( curr_file_name), "rb");
   extern int n_obs_actually_loaded;
   const char *override_epoch_text = get_environment_ptr( "EPOCH");
   double epoch_shown;
   double override_epoch = (*override_epoch_text ?
                     extract_epoch( override_epoch_text) : 0.);

   curr_object_name = obj_info[obj_idx].obj_name;
   SetCursor( AfxGetApp( )->LoadStandardCursor( IDC_WAIT));
            /* Start a bit ahead of the actual data,  just in case */
            /* there's a #Sigma or similar command in there: */
   file_offset = obj_info[obj_idx].file_offset - 40L;
   if( file_offset < 0L)
      file_offset = 0L;
   fseek( ifile, file_offset, SEEK_SET);

   constraints = "";
   obs = load_object( ifile, obj_info + obj_idx, &orbit_epoch,
                                           &epoch_shown, orbit);
   n_obs = obj_info[obj_idx].n_obs = n_obs_actually_loaded;
   fclose( ifile);
   m_step_size = 3.;
   put_epoch_text( buff, epoch_shown);
   set_locs( orbit, orbit_epoch, obs, n_obs);
   m_epoch = buff;
   SetCursor( AfxGetApp( )->LoadStandardCursor( IDC_ARROW));
   if( !n_obs)
      return;
   if( obs_data)
      FFREE( obs_data);
   SetCursor( AfxGetApp( )->LoadStandardCursor( IDC_WAIT));
   obs_data = obs;

   ((CListBox *)GetDlgItem( IDC_RESIDUALS))->ResetContent( );
   m_r1 = "1";
   m_r2 = "1";
   UpdateData( FALSE);     /* 'False' indicates 'move data to edit boxes' */
   ResetPerturbers( );
   Reset_r1_and_r2( );

   UpdateElementDisplay( 1);
   UpdateResidualDisplay( );
   SetCursor( AfxGetApp( )->LoadStandardCursor( IDC_ARROW));
}

void COrbitDlg::OnDblclkObject()
{
   CListBox* pListBox = (CListBox*)GetDlgItem( IDC_LIST_ASTEROIDS);
   int selected = pListBox->GetCurSel( );
   struct stat s;

                  /* If the file was modified,  we should reload everything: */
   if( !stat( CT2A( curr_file_name), &s))
      if( curr_file_time != s.st_mtime)
         {
         char obj_name[80];
         int i;

         strcpy( obj_name, obj_info[selected].obj_name);
         debug_printf( "Modification time changed: reloading file\n");
         LoadAFile( CT2A( curr_file_name));
         selected = 0;
         for( i = 0; i < n_objects; i++)
            if( !strcmp( obj_name, obj_info[i].obj_name))
               selected = i;
         }
   LoadAnObject( selected);
}

void COrbitDlg::OnClickedSave()
{
   // TODO: Add your control notification handler code here
   char filename[_MAX_DIR];

   get_file_from_dialog( FALSE, "", "*.*", filename, NULL);

   if( *filename)
      {
      double orbit2[6];
      double epoch_in_edit_box;
      FILE *ofile;

      UpdateData( TRUE);     /* Get changes made in edit boxes */
      epoch_in_edit_box = extract_epoch( CT2A( m_epoch));
      memcpy( orbit2, orbit, 6 * sizeof( double));
      integrate_orbit( orbit2, orbit_epoch, epoch_in_edit_box);
      if( ofile = fopen( filename, "w"))
         {
         char buff[200];
         FILE *ifile = fopen( elements_filename, "rb");

         if( ifile)
            {
            while( fgets( buff, sizeof( buff), ifile))
               fputs( buff, ofile);
            fclose( ifile);
            }
         fclose( ofile);
         }
      else
         MessageBox( CA2T( filename, CP_UTF8), _T( "Not opened"), MB_OK);
      store_solution( (const OBSERVE FAR *)obs_data, n_obs, orbit2,
                      epoch_in_edit_box, perturbers);
      }
}

void COrbitDlg::UpdateResidualDisplay()
{
   CListBox *pResidualBox = (CListBox *)GetDlgItem( IDC_RESIDUALS);
   const int n_selected = pResidualBox->GetSelCount( );
   char buff[640];
   OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;
   int i;

   // TODO: Add your control notification handler code here
// if( n_selected == 1 && (GetAsyncKeyState( VK_SHIFT) & 0x8001) == 0x8001)
//    {
//    int selected;
//
//    pResidualBox->GetSelItems( 1, &selected);
//    if( prev_shifted_residual == -1)
//       prev_shifted_residual = selected;
//    else
//       {
//       int temp, n1 = prev_shifted_residual, n2 = selected;
//
//       if( n2 < n1)
//          {
//          temp = n1;
//          n1 = n2;
//          n2 = temp;
//          }
//       for( i = n1; i <= n2; i++)
//          obs[i].is_included ^= 1;
//       UpdateElementDisplay( 0);
//
//       prev_shifted_residual = -1;
//       }
//    }
   for( i = 0; i < n_obs; i++)
      if( pResidualBox->GetSel( i))
         obs[i].flags |= OBS_IS_SELECTED;
      else
         obs[i].flags &= ~OBS_IS_SELECTED;
   generate_obs_text( obs, n_obs, buff);

   SetDlgItemText( IDC_STATION_INFO, CA2T( buff, CP_UTF8));
   AdjustControls( );
   GetDlgItem( IDC_STATION_INFO)->Invalidate();
}

void COrbitDlg::OnSelchangeResiduals()
{
   // TODO: Add your control notification handler code here
   UpdateResidualDisplay( );
}

#define TRUE_OR_FALSE( X)  ((X) ? TRUE : FALSE)

void COrbitDlg::OnClickedMakeEphemeris()
{
   // TODO: Add your control notification handler code here
   if( obs_data && n_obs)
      {
      CEphem dlg;
      int lines_read = 0;
      char tstr[90], object_name[80];
      const char *envptr;
      FILE *ifile = fopen( "startup.mar", "rb");
      extern const char *residual_filename;
      extern double ephemeris_mag_limit;

                                 /* For purposes of making a pseudo-MPEC: */
      write_residuals_to_file( residual_filename, CT2A( curr_file_name), n_obs,
                        (OBSERVE FAR *)obs_data, RESIDUAL_FORMAT_SHORT);
      create_obs_file( (OBSERVE FAR *)obs_data, n_obs, 0);
      dlg.epoch = orbit_epoch;
      dlg.m_number_steps = 10;
      dlg.m_ephem_step = _T( "1");
      dlg.m_ephem_type =(ephemeris_output_options & 7);
      dlg.m_alt_az = ((ephemeris_output_options & OPTION_ALT_AZ_OUTPUT) ? TRUE : FALSE);
      dlg.m_motion = ((ephemeris_output_options & OPTION_MOTION_OUTPUT) ? TRUE : FALSE);
      dlg.m_separate_motions = ((ephemeris_output_options & OPTION_SEPARATE_MOTIONS) ? TRUE : FALSE);
      dlg.m_round_step = ((ephemeris_output_options & OPTION_ROUND_TO_NEAREST_STEP) ? TRUE : FALSE);
      dlg.m_radial_velocity = ((ephemeris_output_options & OPTION_RADIAL_VEL_OUTPUT) ? TRUE : FALSE);
      dlg.m_phase_angle = ((ephemeris_output_options & OPTION_PHASE_ANGLE_OUTPUT) ? TRUE : FALSE);
      dlg.m_phase_angle_bisector =
                  TRUE_OR_FALSE( ephemeris_output_options & OPTION_PHASE_ANGLE_BISECTOR);
      dlg.m_helio_ecliptic =
                  TRUE_OR_FALSE( ephemeris_output_options & OPTION_HELIO_ECLIPTIC);
      dlg.m_show_sigmas =
                  TRUE_OR_FALSE( ephemeris_output_options & OPTION_SHOW_SIGMAS);
      dlg.m_topo_ecliptic =
                  TRUE_OR_FALSE( ephemeris_output_options & OPTION_TOPO_ECLIPTIC);
      dlg.m_human_readable =
                 !TRUE_OR_FALSE( ephemeris_output_options & OPTION_COMPUTER_FRIENDLY);
      dlg.m_ground_track =
                 TRUE_OR_FALSE( ephemeris_output_options & OPTION_GROUND_TRACK);
      dlg.m_speed =
                 TRUE_OR_FALSE( ephemeris_output_options & OPTION_SPACE_VEL_OUTPUT);
      get_object_name( object_name, ((OBSERVE FAR *)obs_data)->packed_id);
      dlg.obj_name = object_name;
      dlg.m_mag_limit = ephemeris_mag_limit;
      envptr = get_environment_ptr( "EPHEM_START");
      if( envptr && *envptr)
         strcpy( tstr, envptr);
      else
         full_ctime( tstr, dlg.epoch,
            FULL_CTIME_YMD | FULL_CTIME_DATE_ONLY | CALENDAR_JULIAN_GREGORIAN);
      dlg.m_day = tstr;
      envptr = get_environment_ptr( "EPHEM_STEPS");
      if( envptr && *envptr)
         {
         sscanf( envptr, "%d %s", &dlg.m_number_steps, tstr);
         dlg.m_ephem_step = CA2T( tstr, CP_UTF8);
         }
      dlg.m_lon = dlg.m_lat = _T( "");
      dlg.obs = (OBSERVE *)obs_data;
      dlg.n_obs = n_obs;
      if( ifile)
         {
         while( fgets( tstr, 90, ifile))
            {
            if( !memcmp( tstr, "11 lat/lon", 10))
               {
               double lat, lon;

               sscanf( tstr + 12, "%lf%lf", &lon, &lat);
               sprintf( tstr, "%c %.4lf", (lon < 0. ? 'W' : 'E'),
                           fabs( lon));
               dlg.m_lon = tstr;
               sprintf( tstr, "%c %.4lf", (lat < 0. ? 'S' : 'N'),
                           fabs( lat));
               dlg.m_lat = tstr;
               }
            lines_read++;
            }
         fclose( ifile);
         }

      envptr = get_environment_ptr( "EPHEM_LAT");
      if( envptr && *envptr)
         dlg.m_lat = envptr;
      envptr = get_environment_ptr( "EPHEM_LON");
      if( envptr && *envptr)
         dlg.m_lon = envptr;
      envptr = get_environment_ptr( "EPHEM_MPC_CODE");
      if( envptr && *envptr)
         {
         dlg.m_use_mpc_code = (int)( *envptr - '0');
         dlg.m_mpc_code = CA2T( envptr + 2, CP_UTF8);
         }

      memcpy( dlg.orbit, orbit, 6 * sizeof( double));
      GetPerturberMask( );
      if( dlg.DoModal( ) == IDOK)
         {
//       sprintf( tstr, "%d %s",
//                      dlg.m_number_steps,  CT2A( dlg.m_ephem_step));
         sprintf( tstr, "%d ", dlg.m_number_steps);
         strcat( tstr, CT2A( dlg.m_ephem_step));
         set_environment_ptr( "EPHEM_STEPS", tstr);
         set_environment_ptr( "EPHEM_START", CT2A( dlg.m_day));
         set_environment_ptr( "EPHEM_LON", CT2A( dlg.m_lon));
         set_environment_ptr( "EPHEM_LAT", CT2A( dlg.m_lat));
         sprintf( tstr, "%d ", dlg.m_use_mpc_code);
         strcat( tstr, CT2A( dlg.m_mpc_code));
         debug_printf( "New MPC code text: '%s'\n", tstr);
         set_environment_ptr( "EPHEM_MPC_CODE", tstr);
         ephemeris_output_options = dlg.GetEphemerisBitmask();
         }
      }
   else
      {
      char buff[80];

      strcpy( buff, get_find_orb_text( 3));  /* "No orbit to make an ephemeris!" */
      MessageBox( CA2T( buff, CP_UTF8), _T( "Find_Orb"), MB_OK);
      }
}

void COrbitDlg::OnClickedSaveResids()
{
   // TODO: Add your control notification handler code here
   char filename[_MAX_DIR];

   if( !n_obs)
      {
      strcpy( filename, get_find_orb_text( 2));  /* "No residuals to save!" */
      MessageBox( CA2T( filename, CP_UTF8), _T( "Save residuals"), MB_OK);
      return;
      }
   get_file_from_dialog( FALSE, "", "*.*", filename, NULL);
   if( *filename)
      {
      int residual_format;

#ifdef _WIN32        /* MS is different. */
      if( _stricmp( filename + strlen( filename) - 4, ".res"))
#else
      if( stricmp( filename + strlen( filename) - 4, ".res"))
#endif
         residual_format = RESIDUAL_FORMAT_FULL_NO_TABS;
      else
         residual_format = RESIDUAL_FORMAT_SHORT;
      write_residuals_to_file( filename, CT2A( curr_file_name),
                    n_obs, (OBSERVE FAR *)obs_data, residual_format);
      }
}

void COrbitDlg::OnDblclkResiduals()
{
   // TODO: Add your control notification handler code here
   OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;
   CListBox* pListBox = (CListBox*)GetDlgItem( IDC_RESIDUALS);
   int selected = pListBox->GetCurSel( );

   obs[selected].is_included ^= 1;
   UpdateElementDisplay( 0);
}

void COrbitDlg::OnClickedAbout()
{
   // TODO: Add your control notification handler code here
   CAbout dlg;

   dlg.DoModal( );
}

void COrbitDlg::OnClickedVaisala()
{
   // TODO: Add your control notification handler code here
   CGenericEntry dlg;

   dlg.m_caption = "Enter peri/aphelion distance:";
   dlg.m_text = "2.3";
   if( dlg.DoModal( ) == IDOK && atof( CT2A( dlg.m_text)))
      {
      m_r1 = dlg.m_text;
      UpdateData( FALSE);  /* 'False' indicates 'move data to xtrols' */
      ImproveOrbitSolution( 2, 1);     /* initial Vaisala... */
      ImproveOrbitSolution( 4, 1);     /* ...then linearize */
      }
}

#define EARTH_CLOSE_APPROACH .01

void COrbitDlg::OnClickedAutoSolve()
{
   // TODO: Add your control notification handler code here
   if( !n_obs || !obs_data)
      {
      char buff[80];

      strcpy( buff, get_find_orb_text( 1));    /* "No orbit to improve!" */
      MessageBox( CA2T( buff, CP_UTF8), _T( "Find_Orb"), MB_OK);
      }
   else
      {
      int pass, bug_out = 0, perturbers_used = 0;

      for( pass = 0; !bug_out && pass < 2; pass++)
         {
         double rms[20];
         int i, iter, done = 0;
         OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;

         rms[0] = compute_rms( obs, n_obs);
         for( iter = 1; !bug_out && iter < 20 && !done; iter++)
            {
            ImproveOrbitSolution( pass, 1);
            rms[iter] = compute_rms( (OBSERVE FAR *)obs_data, n_obs);
            if( iter > 1 && rms[iter] < 10.)
               if( rms[iter] > rms[iter - 2] * .9) /* not getting much improvement */
                  done = 1;
            if( rms[iter] > 100000.)         /* wups!  explosion occurred */
               bug_out = 1;
            }
         if( pass)       /* may need to repeat */
            {
            int n1, n2, perturbers_to_use = 0;
            double curr_arc;

            for( i = 0; i < n_obs; i++)
               if( obs[i].is_included && obs[i].r < EARTH_CLOSE_APPROACH)
                  perturbers_to_use = 4 | (1 << 10);    /* earth & moon */
            for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
               ;
            n1 = i;
            for( i = n_obs - 1; i && !obs[i].is_included; i--)
               ;
            n2 = i;

            curr_arc = obs[n2].jd - obs[n1].jd;
            if( curr_arc > 60.)
               perturbers_to_use |= (1 << 4);       /* Jupiter */
            if( curr_arc > 250.)
               perturbers_to_use = 255;            /* Mercury...Neptune */
            if( perturbers_to_use == perturbers_used)
               if( n1 > 0 || n2 < n_obs - 1)   /* gotta extend the arc */
                  {
                  int m1, m2;

                  for( i = n1; i >= 0 && obs[n2].jd - obs[i].jd < 4. * curr_arc; i--)
                     obs[i].is_included = 1;
                  m1 = i + 1;
                  for( i = n2; i < n_obs && obs[i].jd - obs[n1].jd < 4. * curr_arc; i++)
                     obs[i].is_included = 1;
                  m2 = i - 1;
                  if( m1 == n1 && m2 == n2)     /* gotta expand arc by at least */
                     {                          /* one observation: */
                     if( !m1)
                        m2++;
                     else if( m2 == n_obs - 1)
                        m1--;
                     else
                        {
                        if( obs[m1].jd - obs[m1 - 1].jd >
                                                obs[m2 + 1].jd - obs[m2].jd)
                           m2++;
                        else
                           m1--;
                        }
                     obs[m1].is_included = obs[m2].is_included = 1;
                     }
                           /* reconsider the arc: */
                  curr_arc = obs[m2].jd - obs[m1].jd;
                  if( curr_arc > 60.)
                     perturbers_to_use |= (1 << 4);       /* Jupiter */
                  if( curr_arc > 250.)
                     perturbers_to_use = 255;          /* Mercury...Neptune */
                  pass = 0;         /* ensure another iteration */
                  }

            if( perturbers_to_use != perturbers_used)
               {
               for( i = 0; i < 9; i++)
                  if( (perturbers_to_use >> i) & 1)
                     ((CButton *)GetDlgItem( IDC_PLANET1 + i))->SetCheck( 1);
               pass = 0;         /* ensure another iteration */
               perturbers_used = perturbers_to_use;
               }
            }
         }
      }
}

#ifdef REMOVED_CODE
void COrbitDlg::OnChar(UINT nChar, UINT nRepCnt, UINT nFlags)
{
   // TODO: Add your message handler code here and/or call default
   TCHAR text[90];

   _stprintf( text, _T( "Got %02lx (%c)\n"), nChar, nChar);
   MessageBox( text, program_name, MB_OK);
   CDialog::OnChar(nChar, nRepCnt, nFlags);
}
#endif

void COrbitDlg::OnMonteCarlo()
{
   // TODO: Add your control notification handler code here
   extern int using_sr;

// if( using_sr != 0)
      set_statistical_ranging( 0);
   monte_carlo ^= 1;
   if( !monte_carlo)
      SetDlgItemText( IDC_MONTE_CARLO, _T( "Monte Carlo"));
}

int COrbitDlg::RunMonteCarlo( void)
{
   if( n_obs && obs_data)
      {
      OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;
      double *stored_ra_decs;
      extern int using_sr;
      extern int monte_carlo_object_count;
      static double nominal_orbit[6];

      if( !monte_carlo_object_count)
         {
         const double rms = compute_rms( (OBSERVE FAR *)obs_data, n_obs);
         extern double max_monte_rms;

         max_monte_rms = sqrt( rms * rms + monte_noise * monte_noise);
         memcpy( nominal_orbit, orbit, 6 * sizeof( double));
         }
      stored_ra_decs = add_gaussian_noise_to_obs( n_obs, obs, monte_noise);
      if( !using_sr)
         {
         memcpy( orbit, nominal_orbit, 6 * sizeof( double));
         compute_covariance = 0;
         ImproveOrbitSolution( 1, 2);
         compute_covariance = 1;
         }
      else
         {
         int i;

         find_nth_sr_orbit( orbit, obs, n_obs, monte_carlo_object_count);
         adjust_herget_results( obs, n_obs, orbit);
                        /* epoch is that of first included observation: */
         for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
            ;
         integrate_orbit( orbit, obs[i].jd, orbit_epoch);
         set_locs( orbit, orbit_epoch, obs, n_obs);
         Reset_r1_and_r2( );
         UpdateElementDisplay( 1);
         UpdateResidualDisplay( );
         }
      restore_ra_decs_mags_times( n_obs, obs, stored_ra_decs);
      free( stored_ra_decs);
      }
   return( 0);
}

/* We get a timer tick every 500 milliseconds (see above SetTimer call). */
/* The following code does Monte Carlo orbits until a quarter second has */
/* passed.  The idea is that in this way,  the program will consume about */
/* half the CPU time,  rather than either hogging too much time or not    */
/* making full use of what's available.   */

void COrbitDlg::OnTimer( TIMER_TYPE pnIDEvent)
{
   // TODO: Add your message handler code here and/or call default
   extern char *runtime_message;

   CDialog::OnTimer(pnIDEvent);
   if( monte_carlo)
      {
      clock_t t0 = clock( );

      while( clock( ) < t0 + CLOCKS_PER_SEC / 4)
         RunMonteCarlo( );
      }
   if( runtime_message)
      {
      SetDlgItemText( IDC_STATION_INFO, CA2T( runtime_message, CP_UTF8));
      AdjustControls( );
      }
}

void COrbitDlg::OnGauss()
{
   OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;
   static int soln_number = 0;
   double new_epoch;

   // TODO: Add your control notification handler code here
   constraints = "";
   new_epoch = convenient_gauss( obs, n_obs, orbit, 1., soln_number++);
   if( !new_epoch && soln_number > 1)
      {                       /* wrap around to the first Gauss soln */
      soln_number = 0;
      new_epoch = convenient_gauss( obs, n_obs, orbit, 1., soln_number++);
      }
   if( !new_epoch)
      {
      char buff[80];

      strcpy( buff, get_find_orb_text( 4));   /* "Method of Gauss failed!" */
      MessageBox( CA2T( buff, CP_UTF8), program_name, MB_OK);
      }
   else
      {
      orbit_epoch = new_epoch;
      set_locs( orbit, orbit_epoch, obs, n_obs);
      Reset_r1_and_r2( );
      UpdateElementDisplay( 1);
      }
}

int find_worst_observation( const OBSERVE FAR *obs, const int n_obs);

void COrbitDlg::OnWorst()
{
   // TODO: Add your control notification handler code here

   const int worst = find_worst_observation( (OBSERVE FAR *)obs_data, n_obs);

   if( worst >= 0)
      {
      CListBox* pListBox = (CListBox*)GetDlgItem( IDC_RESIDUALS);

//    pListBox->SetCurSel( worst);     Valid for single-selection box only
      pListBox->SelItemRange( FALSE, 0, pListBox->GetCount( ) - 1);
      pListBox->SetSel( worst, TRUE);
      UpdateResidualDisplay( );
      }
}

void COrbitDlg::OnDoubleclickedWorst()
{
   // TODO: Add your control notification handler code here

}


void COrbitDlg::OnFilterObs()
{
   // TODO: Add your control notification handler code here
   extern int use_blunder_method;

   if( use_blunder_method)
      {
      use_blunder_method = 2;
      ImproveOrbitSolution( 1, 1);
      use_blunder_method = 1;
      }
   else        /* "traditional" filtering */
      {
      const int rval = filter_obs( (OBSERVE FAR *)obs_data, n_obs,
                              max_residual_for_filtering);

      if( rval != FILTERING_CHANGES_MADE)
         {
         char buff[80];

         strcpy( buff, get_find_orb_text( 5));      /* "No changes made!" */
         MessageBox( CA2T( buff, CP_UTF8), program_name, MB_OK);
         }
      else
         ImproveOrbitSolution( 1, 1);
      }
}

void COrbitDlg::OnSelcancelListAsteroids()
{
   // TODO: Add your control notification handler code here

}

void COrbitDlg::OnAsteroids()
{
   // TODO: Add your control notification handler code here

}

void COrbitDlg::OnSettings()
{
   // TODO: Add your control notification handler code here
   CSettings dlg;
   extern char default_comet_magnitude_type;
   extern double probability_of_blunder;
   char tstr[30];
   extern double overobserving_time_span;
   extern unsigned overobserving_ceiling;
   extern int use_blunder_method;
   extern bool use_sigmas;
   extern int forced_central_body;
   extern int apply_debiasing;

   if( element_format & ELEM_OUT_HELIOCENTRIC_ONLY)
      dlg.m_element_center = forced_central_body + 2;
   else
      dlg.m_element_center = 0;
   dlg.m_element_precision = element_precision;
   dlg.m_constraints = constraints;
   sprintf( tstr, "%.2lf", max_residual_for_filtering);
   dlg.m_max_residual = tstr;
   dlg.m_monte_noise = monte_noise;
   dlg.m_reference = get_environment_ptr( "REFERENCE");
   dlg.m_physical_model = n_extra_params;
   dlg.m_comet_mags_total = (default_comet_magnitude_type == 'T');
   dlg.m_precise_residuals = precise_residuals;
   dlg.m_probability_of_blunder = probability_of_blunder * 100.;  /* cvt to % */
   dlg.m_overobserving_time_span = overobserving_time_span;
   dlg.m_overobserving_ceiling = overobserving_ceiling;
   dlg.m_use_blunder_method = use_blunder_method;
   dlg.m_alternative_elements = ((element_format & ELEM_OUT_ALTERNATIVE_FORMAT) != 0);
   dlg.m_use_weights = use_sigmas;
   dlg.m_debiasing = apply_debiasing;
   if( dlg.DoModal( ) == IDOK)
      {
      element_precision = dlg.m_element_precision;
      if( !dlg.m_element_center)
         element_format = 0;
      else
         {
         element_format = ELEM_OUT_HELIOCENTRIC_ONLY;
         forced_central_body = dlg.m_element_center - 2;
         }
      if( dlg.m_alternative_elements)
         element_format |= ELEM_OUT_ALTERNATIVE_FORMAT;
      max_residual_for_filtering = atof( CT2A( dlg.m_max_residual));
      monte_noise = dlg.m_monte_noise;
      constraints = CT2A( dlg.m_constraints);
      if( n_extra_params != dlg.m_physical_model)
         solar_pressure[0] = solar_pressure[1] = solar_pressure[2] = 0.;
      n_extra_params = dlg.m_physical_model;
      set_environment_ptr( "REFERENCE", CT2A( dlg.m_reference));
      default_comet_magnitude_type =
                  (dlg.m_comet_mags_total ? 'T' : 'N');
      precise_residuals = dlg.m_precise_residuals;
      probability_of_blunder = dlg.m_probability_of_blunder / 100.;
      overobserving_time_span = dlg.m_overobserving_time_span;
      overobserving_ceiling = dlg.m_overobserving_ceiling;
      use_blunder_method = dlg.m_use_blunder_method;
      use_sigmas = (dlg.m_use_weights ? true : false);
      apply_debiasing = dlg.m_debiasing;
      UpdateElementDisplay( 1);
      }
}

void COrbitDlg::OnSelchangeListAsteroids()
{
   // TODO: Add your control notification handler code here
   const int selected = ((CListBox*)GetDlgItem( IDC_LIST_ASTEROIDS))->GetCurSel( );

   if( selected >= 0 && selected < n_objects)
      {
      char buff[440];
      OBJECT_INFO *ids = obj_info + selected;

      sprintf( buff, "Object %d of %d: %s\n", selected + 1, n_objects,
                              ids->obj_name);
      sprintf( buff + strlen( buff), "%d observations; ", ids->n_obs);
      make_date_range_text( buff + strlen( buff),
                                      (double)ids->jd_start,
                                      (double)ids->jd_end);
      strcat( buff, "\n");
      add_version_and_de_text( buff + strlen( buff));
      SetDlgItemText( IDC_STATION_INFO, CA2T( buff, CP_UTF8));
      GetDlgItem( IDC_STATION_INFO)->Invalidate();
      AdjustControls( );
      }
}

void COrbitDlg::OnSetWeight()
{
   // TODO: Add your control notification handler code here
   CListBox* pListBox = (CListBox*)GetDlgItem( IDC_RESIDUALS);
   int n_selected;

   n_selected = pListBox->GetSelCount( );
   if( !n_selected)
      MessageBox(
          _T( "No observations selected!  Select one or more\n")
          _T( "observations,  click on this button,  and you'll\n")
          _T( "be prompted to enter uncertainties for them.\n"),
              program_name, MB_OK);
   else
      {
      CGenericEntry dlg;

      dlg.m_caption = _T( "Enter uncertainty in arcseconds:");
      dlg.m_text = _T( "1");
      if( dlg.DoModal( ) == IDOK)
         {
         double tval = atof( CT2A( dlg.m_text));
         int *selections;

         if( !tval)        /* 'm' (mag) or 't' (time) uncertainty */
            tval = atof( CT2A( dlg.m_text) + 1);
         if( tval &&
                (selections = (int *)calloc( n_selected, sizeof( int))))
            {
            int i;
            OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;

            pListBox->GetSelItems( n_selected, selections);
            for( i = 0; i < n_selected; i++)
               switch( dlg.m_text[0])
                  {
                  case 'm': case 'M':
                     obs[selections[i]].mag_sigma = tval;
                     break;
                  case 't': case 'T':
                     obs[selections[i]].time_sigma = tval / seconds_per_day;
                     break;
                  default:
                     obs[selections[i]].posn_sigma = tval;
                     break;
                  }
            free( selections);
            UpdateElementDisplay( 0);
            UpdateResidualDisplay( );
            }
         }
      }
}

void COrbitDlg::OnToggleObs()
{
   // TODO: Add your control notification handler code here
   CListBox* pListBox = (CListBox*)GetDlgItem( IDC_RESIDUALS);
   int n_selected;

   n_selected = pListBox->GetSelCount( );
   if( !n_selected)
      MessageBox( \
           _T( "No observations selected!  Select one or more\n")
           _T( "observations,  click on this button,  and\n")
           _T( "they will be toggled.\n"), program_name, MB_OK);
   else
      {
      int *selections, n_included = 0, i;
      OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;

      selections = (int *)calloc( n_selected, sizeof( int));
      assert( selections);
      if( selections)
         {
         pListBox->GetSelItems( n_selected, selections);
         for( i = 0; i < n_selected; i++)
            if( obs[selections[i]].is_included)
               n_included++;
         for( i = 0; i < n_selected; i++)
            obs[selections[i]].is_included = (n_included <= n_selected / 2);
         free( selections);
         }
      UpdateElementDisplay( 0);
      UpdateResidualDisplay( );
      }
}

void COrbitDlg::OnOrbitalElements()
{
   // TODO: Add your control notification handler code here

   MessageBox( _T( "Orbital elements clicked"), program_name, MB_OK);
}

char *get_placement_text( HWND hwnd)
{
   char *buff = (char *)malloc( 3 * sizeof( WINDOWPLACEMENT) + 3);
   int i;
   WINDOWPLACEMENT place;
   unsigned char *tptr = (unsigned char *)&place;

   if( buff)
      {
      place.length = sizeof( WINDOWPLACEMENT);
      GetWindowPlacement( hwnd, &place);
      for( i = 0; i < sizeof( WINDOWPLACEMENT); i++)
         sprintf( buff + i * 3, "%02x ", (unsigned)tptr[i]);
      }
   return( buff);
}

void COrbitDlg::OnDestroy()
{
   char *place_text;
   extern char default_comet_magnitude_type;

   free_sigma_recs( );
#ifdef OBSOLETE
   sprintf( buff, "%c,%d,%d,%d,%.2lf,%.2lf,%d",
               default_comet_magnitude_type,
               element_format, element_precision,
               ephemeris_output_options,
               max_residual_for_filtering, monte_noise,
               precise_residuals);
   set_environment_ptr( "SETTINGS", buff),
#endif
   store_defaults( ephemeris_output_options, element_format,
         element_precision, max_residual_for_filtering,
         monte_noise);
   place_text = get_placement_text( this->m_hWnd);
   set_environment_ptr( "WINDOW_PLACEMENT", place_text);
   free( place_text);

   CDialog::OnDestroy();

   // TODO: Add your message handler code here

}

void COrbitDlg::OnLButtonUp(UINT nFlags, CPoint point)
{
   // TODO: Add your message handler code here and/or call default
   CWnd *child = this->ChildWindowFromPoint( point);

   if( child == GetDlgItem( IDC_ORBITAL_ELEMENTS))
      {
      show_commented_elements = !show_commented_elements;
//    MessageBox( _T( "Left button in orbital elements"),
//                _T( "Find_Orb"), MB_OK);
      UpdateElementDisplay( 1);
      AdjustControls( );
      }

   CDialog::OnLButtonUp(nFlags, point);
}

/* A right-click over the observation/station info,  or over the orbital
elements,  leads to a small popup wherein one can choose to save the data
in those areas to a file or copy said text to the clipboard: */

int copy_buffer_to_clipboard( const char *contents, const long length);
int clipboard_to_file( const char *filename, const int append);      /* ephem0.cpp */

void COrbitDlg::OnRButtonUp(UINT nFlags, CPoint point)
{
   // TODO: Add your message handler code here and/or call default

   CWnd *child = this->ChildWindowFromPoint( point);
   int curr_control = 0;;

   if( child == GetDlgItem( IDC_STATION_INFO))
      curr_control = IDC_STATION_INFO;
   if( child == GetDlgItem( IDC_ORBITAL_ELEMENTS))
      curr_control = IDC_ORBIT1;

   if( curr_control)
      {
      CMenu mnuPopup, *mnuPopupMenu;
      CPoint tpoint = point;
      DWORD selection;

      mnuPopup.LoadMenu( IDR_POPUP1);
                    // Get a pointer to the first item of the menu
      mnuPopupMenu = mnuPopup.GetSubMenu(0);
      ClientToScreen( &tpoint);
      selection = mnuPopupMenu->TrackPopupMenu(
            TPM_LEFTALIGN | TPM_RIGHTBUTTON | TPM_NONOTIFY | TPM_RETURNCMD,
                              tpoint.x, tpoint.y, this);
      if( selection == IDM_READ_CLIPBOARD)
         {
         const char *filename = "obs_temp.txt";

         if( !clipboard_to_file( filename, 0))
            LoadAFile( filename);
         }
      else if( selection == IDM_SAVE_TO_FILE &&
                     curr_control == IDC_ORBIT1)
         OnClickedSave( );
      else if( curr_control)
         if( selection == IDM_SAVE_TO_FILE || selection == IDM_COPY_TO_CLIPBOARD)
            {
            CString str;

            this->GetDlgItemText( curr_control, str);
            if( selection == IDM_COPY_TO_CLIPBOARD)
               {
               int rval = copy_buffer_to_clipboard( CT2A( str),
                                (long)strlen( CT2A( str)));

               if( rval)
                  {
                  char tbuff[80];

                  sprintf( tbuff, "Clipboard rval %d\n", rval);
                  AfxMessageBox( CA2T( tbuff, CP_UTF8));
                  }
               }
            else
               {
               char filename[_MAX_DIR];

               get_file_from_dialog( FALSE, "", "*.*", filename, NULL);

               if( *filename)
                  {
                  FILE *ofile = fopen( filename, "w");

                  fwrite( str, 1, strlen( CT2A( str)), ofile);
                  fclose( ofile);
                  }
               }
            }
      }
#ifdef FIND_OUT_WHAT_RIGHT_CLICKS
   else
      {
      TCHAR buff[80];

      _sntprintf( buff, sizeof( buff) / sizeof( buff[0]),
                                    _T( "Right-click %p\n"), child);
      MessageBox( buff, program_name, MB_OK);
      }
#endif
   CDialog::OnRButtonUp(nFlags, point);
}

#ifdef NO_LONGER_NEEDED_I_HOPE
BOOL COrbitDlg::OnCommand( UINT wParam, LONG lParam)
{
   CWnd::OnCommand( wParam, lParam);

   return TRUE;
}
#endif

int debug_printf( const char *format, ...);                /* runge.cpp */

static int count_lines_in_buffer( const TCHAR *buff)
{
   int i, rval;

   for( i = rval = 0; buff[i]; i++)
      if( buff[i] == 10)
         rval++;
   if( i && buff[i - 1] != 10)      /* last line lacks a LF */
      rval++;
   return( rval);
}


/* GetWindowRect gives a rectangle in screeen coords relative to the
upper left corner of the screen.
   SetWindowPos moves and/or resizes.  The position is in client coords;
width & height are in pixels.
   MoveWindow:  For a top-level window, the position and dimensions are
relative to the upper-left corner of the screen. For a child window,
they are relative to the upper-left corner of the parent window's client
area.  */

void COrbitDlg::AdjustControls( void)
{
   CListBox* pListBox = (CListBox*)GetDlgItem( IDC_RESIDUALS);
   CWnd *station_info = GetDlgItem( IDC_STATION_INFO);
   CWnd *orbit_elems = GetDlgItem( IDC_ORBIT1);
   CWnd *orbit_elems_group = GetDlgItem( IDC_ORBITAL_ELEMENTS);

   // TODO: Add your message handler code here
   if( pListBox && station_info && orbit_elems && orbit_elems_group)
      {
      CRect stn_rect, list_box_rect;
      CRect orbit_rect, orbit_group_rect;
      static int font_height;
      CRect dlg_rect;
      int cy, new_height, new_elems_ht;
      int y0;
      TCHAR buff[1000];

      station_info->GetWindowRect( &stn_rect);
      ScreenToClient( &stn_rect);
      pListBox->GetWindowRect( &list_box_rect);
      ScreenToClient( &list_box_rect);
      orbit_elems->GetWindowRect( &orbit_rect);
      ScreenToClient( &orbit_rect);
      orbit_elems_group->GetWindowRect( &orbit_group_rect);
      ScreenToClient( &orbit_group_rect);
      GetWindowRect( &dlg_rect);
      ScreenToClient( &dlg_rect);

      if( !font_height)    /* we start out with a six-line rect */
         font_height = stn_rect.Height( ) / 7;
      GetDlgItemText( IDC_STATION_INFO, buff, sizeof( buff));
      new_height = (count_lines_in_buffer( buff) + 2) * font_height;
               /* need an extra half line in Win64 version : */
      new_height += font_height / 2;
      GetDlgItemText( IDC_ORBIT1, buff, sizeof( buff));
      new_elems_ht = (count_lines_in_buffer( buff) + 0) * font_height;
      cy = dlg_rect.Height( ) - new_height;
      station_info->SetWindowPos( NULL, stn_rect.left, cy,
                                stn_rect.Width( ), new_height, SWP_NOZORDER);

//    pListBox->SetWindowPos( this, 0, 0,
//                list_box_rect.Width( ),             /* don't change width */
//                cy - list_box_rect.top,
//                SWP_NOMOVE | SWP_NOZORDER);
      y0 = orbit_rect.top + new_elems_ht;
      pListBox->SetWindowPos( this, list_box_rect.left,     /* don't move "sideways" */
                  y0,
                  list_box_rect.Width( ),             /* don't change width */
                  cy - y0,
                  SWP_NOZORDER);
      orbit_elems->SetWindowPos( this, 0, 0,
                  orbit_rect.Width( ), new_elems_ht, SWP_NOMOVE | SWP_NOZORDER);
      orbit_elems_group->SetWindowPos( this, 0, 0,
                  orbit_group_rect.Width( ),
                  new_elems_ht + orbit_group_rect.Height( ) - orbit_rect.Height( ),
                  SWP_NOMOVE | SWP_NOZORDER);
      }
}

void COrbitDlg::OnSize(UINT nType, int cx, int cy)
{
   CDialog::OnSize(nType, cx, cy);

   AdjustControls( );
}

/* Solution found from http://www.codeguru.com/forum/showthread.php?t=318933 */
/* Modified using http://www.flounder.com/getminmaxinfo.htm ,  which has
   lots of useful info on resizing windows & keeping controls to match */

void COrbitDlg::OnGetMinMaxInfo(MINMAXINFO FAR* lpMMI)
{
  // set the minimum tracking width
  // and the minimum tracking height of the window
  if( OriginalDlgRect.Height( ))
     {
     lpMMI->ptMinTrackSize.x = OriginalDlgRect.Width( );
     lpMMI->ptMaxTrackSize.x = OriginalDlgRect.Width( );
     lpMMI->ptMinTrackSize.y = OriginalDlgRect.Height( );
     }
}

void COrbitDlg::OnTogglePerturbers()
{
   // TODO: Add your control notification handler code here
   int default_perturbers = 0x7fe;   /* Merc-Pluto including the moon */
                       /* One _could_ use defaults of 0x3fe (same, but EMB) */
                       /* or 0x5fe (Merc-Nept and including the moon) */
                       /* 0x7007fe adds in Ceres,  Pallas,  and Vesta */

   sscanf( get_environment_ptr( "DEFAULT_PERTURBERS"), "%x",
                                            &default_perturbers);
   perturbers = (perturbers > 1 ? 0: default_perturbers);
   ResetPerturbers( );
}

int simplex_method( OBSERVE FAR *obs, int n_obs, double *orbit,  /* orb_func.cpp */
               const double r1, const double r2, const char *constraints);

void COrbitDlg::OnSimplex()
{
   // TODO: Add your control notification handler code here
   OBSERVE FAR *obs = (OBSERVE FAR *)obs_data;
   int i;

   UpdateData( TRUE);     /* Get changes made in edit boxes */
   SetCursor( AfxGetApp( )->LoadStandardCursor( IDC_WAIT));
   GetPerturberMask( );
   simplex_method( obs, n_obs, orbit,
                        parse_distance_text( CT2A( m_r1)),
                        parse_distance_text( CT2A( m_r2)),
                        CT2A( constraints));
   for( i = 0; i < n_obs - 1 && !obs[i].is_included; i++)
      ;
   integrate_orbit( orbit, obs[i].jd, orbit_epoch);
   set_locs( orbit, orbit_epoch, obs, n_obs);
   Reset_r1_and_r2( );
   UpdateElementDisplay( 1);
   UpdateResidualDisplay( );
   SetCursor( AfxGetApp( )->LoadStandardCursor( IDC_ARROW));
}

void COrbitDlg::OnStatisticalRanging()
{
   // TODO: Add your control notification handler code here
   extern int using_sr;

// if( using_sr != 1)
      set_statistical_ranging( 1);
   monte_carlo ^= 1;
   if( !monte_carlo)
      SetDlgItemText( IDC_STATISTICAL_RANGING, _T( "Stat Rang&ing"));
}

HBRUSH COrbitDlg::OnCtlColor(CDC* pDC, CWnd* pWnd, UINT nCtlColor)
{
   HBRUSH hbr = CDialog::OnCtlColor(pDC, pWnd, nCtlColor);

   // TODO: Change any attributes of the DC here
   if( pWnd->GetDlgCtrlID() == IDC_ORBIT1)
      {
      pDC->SetTextColor( bad_elements ? RGB( 255, 0, 0) : RGB( 0, 0, 0));
      pDC->SetBkMode( TRANSPARENT);
      }

   return hbr;
}

void COrbitDlg::OnKillfocusEpoch()
{
   UpdateElementDisplay( 1);
#ifdef DONT_THINK_THIS_IS_NEEDED
   UpdateData( TRUE);     /* Get changes made in edit boxes */

   // TODO: Add your control notification handler code here
   const double epoch_in_edit_box = extract_epoch( m_epoch);

   integrate_orbit( orbit, orbit_epoch, epoch_in_edit_box);
   orbit_epoch = epoch_in_edit_box;
#endif
}

