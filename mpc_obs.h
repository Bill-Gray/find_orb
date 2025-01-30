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

#define OBSERVE struct observe

OBSERVE
   {
   double jd, obs_posn[3], obs_vel[3], vect[3], ra, dec, obs_mag;
   double r,  obj_posn[3], obj_vel[3], solar_r, computed_ra, computed_dec;
   double posn_sigma_1, posn_sigma_2;         /* in arcseconds */
            /* Usually, posn_sigma_1 = RA sigma, posn_sigma_2 = dec sigma. */
   double posn_sigma_theta;   /* tilt angle of uncertainty ellipse */
   double mag_sigma;
   double time_sigma;         /* in days */
   double unc_time;           /* uncTime (semi-systematic timing err) in days */
   double computed_mag;
   double ra_bias, dec_bias;     /* in arcseconds */
   char *second_line;
   int flags, is_included;
   int ref_center;       /* 399 = geocenter, 10 = heliocenter, 0 = SSB... */
   int time_precision, ra_precision, dec_precision, mag_precision;
   char mpc_code[4], packed_id[13], reference[6];
   char columns_57_to_65[10];
   char mag_band, astrometric_net_code, discovery_asterisk, note1, note2, satellite_obs;
   const char **obs_details;
   char *ades_ids;
   };

#define OBJECT_INFO struct object_info

OBJECT_INFO
   {
   double jd_start, jd_end, jd_updated;
   long file_offset;
   char *obj_name;
   int n_obs;
   char packed_desig[13], solution_exists;
   char mpc_codes[16];        /* first five three-char codes,  plus a null terminator */
   };

#define OBJECT_INFO_COMPARE_PACKED           0
#define OBJECT_INFO_COMPARE_LAST_OBSERVED    1
#define OBJECT_INFO_COMPARE_NAME             2
#define OBJECT_INFO_COMPARE_LAST_UPDATED     3

#define MOTION_DETAILS struct motion_details

MOTION_DETAILS
   {
   double xresid, yresid;
   double ra_motion, dec_motion, total_motion;     /* in arcmin/hour */
   double position_angle_of_motion;                /* in degrees,  0-360 */
   double radial_vel;                              /* in km/s */
   double time_residual;                           /* in seconds */
   double cross_residual;                          /* in arcseconds */
   };

#define RADAR_INFO struct radar_info

RADAR_INFO
   {
   double doppler_comp;    /* all in Hz */
   double doppler_obs;
   double doppler_sigma;
   double freq_hz;
   double rtt_comp;        /* all in seconds */
   double rtt_obs;
   double rtt_sigma;
   };

#ifndef MPC_FUNC_H_INCLUDED
   #include "mpc_func.h"
#endif

int compute_radar_info( const OBSERVE *obs, RADAR_INFO *rinfo);

/* So far,  there can be zero,  one,  two,  or three nongravitational
parameters in Find_Orb.  I have some force models in mind that will
require four.  Six will,  I think,  be more than enough (though easy
enough to increase if I'm wrong).  Six nongravs, plus a 'traditional'
state vector,  means a possible total of twelve parameters for an
orbit;  i.e.,  we could have a 12x12 covariance matrix. */

#define MAX_N_NONGRAV_PARAMS 6
#define MAX_N_PARAMS 12
#define ORBIT_CENTER_AUTO  -2

typedef uint64_t ephem_option_t;

/* Bitfield options for ephemeris_in_a_file( ): */
/* Bottom three bits define an ephemeris type.  "Observables" are the */
/* usual RA/dec,  radial velocity,  etc. type output.  "State vector  */
/* output" results in the time (as a JD) and position (in AU,  relative */
/* to the observer,  in Cartesian coordinates) being output.  "Position */
/* output" is the same thing,  minus the velocity components.  "MPCORB  */
/* output" means that the orbital elements will be written out for each */
/* ephemeris step,  on a single line.  "8-line output" is almost the    */
/* same,  except that the elements are written out in the MPC's usual   */
/* eight-line form.  "Close approaches" will result in the range minima */
/* (times and distances) being output.                                  */

#define OPTION_OBSERVABLES             0
#define OPTION_STATE_VECTOR_OUTPUT     1
#define OPTION_POSITION_OUTPUT         2
#define OPTION_MPCORB_OUTPUT           3
#define OPTION_8_LINE_OUTPUT           4
#define OPTION_CLOSE_APPROACHES        5
#define OPTION_FAKE_ASTROMETRY         6

#define EPHEM_OPTION_BIT( N)  (((ephem_option_t)1) << (N))

#define OPTION_ALT_AZ_OUTPUT            EPHEM_OPTION_BIT( 3)
#define OPTION_RADIAL_VEL_OUTPUT        EPHEM_OPTION_BIT( 4)
#define OPTION_MOTION_OUTPUT            EPHEM_OPTION_BIT( 5)
#define OPTION_PHASE_ANGLE_OUTPUT       EPHEM_OPTION_BIT( 6)
#define OPTION_GROUND_TRACK             EPHEM_OPTION_BIT( 8)
#define OPTION_SEPARATE_MOTIONS         EPHEM_OPTION_BIT( 9)

#define OPTION_ROUND_TO_NEAREST_STEP    EPHEM_OPTION_BIT( 10)
#define OPTION_PHASE_ANGLE_BISECTOR     EPHEM_OPTION_BIT( 11)
#define OPTION_HELIO_ECLIPTIC           EPHEM_OPTION_BIT( 12)
#define OPTION_TOPO_ECLIPTIC            EPHEM_OPTION_BIT( 13)

#define OPTION_VISIBILITY               EPHEM_OPTION_BIT( 14)
#define OPTION_SUPPRESS_UNOBSERVABLE    EPHEM_OPTION_BIT( 15)
#define OPTION_SHOW_SIGMAS              EPHEM_OPTION_BIT( 16)
#define OPTION_COMPUTER_FRIENDLY        EPHEM_OPTION_BIT( 17)
      /* Above option means 'ephems are written in format easy for  */
      /* software to read,  instead of in a human-readable format'. */

      /* Added 2015 May 4 at suggestion of Denis Denisenko          */
#define OPTION_MOIDS                    EPHEM_OPTION_BIT( 18)
#define OPTION_SPACE_VEL_OUTPUT         EPHEM_OPTION_BIT( 19)
#define OPTION_LUNAR_ELONGATION         EPHEM_OPTION_BIT( 20)

#define OPTION_SUPPRESS_RA_DEC          EPHEM_OPTION_BIT( 21)
#define OPTION_SUPPRESS_DELTA           EPHEM_OPTION_BIT( 22)
#define OPTION_SUPPRESS_SOLAR_R         EPHEM_OPTION_BIT( 23)
#define OPTION_SUPPRESS_ELONG           EPHEM_OPTION_BIT( 24)

#define OPTION_SUN_ALT                  EPHEM_OPTION_BIT( 25)
#define OPTION_SUN_AZ                   EPHEM_OPTION_BIT( 26)
#define OPTION_MOON_ALT                 EPHEM_OPTION_BIT( 27)
#define OPTION_MOON_AZ                  EPHEM_OPTION_BIT( 28)
#define OPTION_SKY_BRIGHTNESS           EPHEM_OPTION_BIT( 29)


#define OPTION_SUN_TARGET_PA            EPHEM_OPTION_BIT( 30)
#define OPTION_SUN_HELIO_VEL_PA         EPHEM_OPTION_BIT( 31)
#define OPTION_ORBIT_PLANE_ANGLE        EPHEM_OPTION_BIT( 32)
#define OPTION_GALACTIC_COORDS          EPHEM_OPTION_BIT( 33)
#define OPTION_GALACTIC_CONFUSION       EPHEM_OPTION_BIT( 34)
#define OPTION_SNR                      EPHEM_OPTION_BIT( 35)
#define OPTION_EXPOSURE_TIME            EPHEM_OPTION_BIT( 36)
#define OPTION_EXPLANATIONS             EPHEM_OPTION_BIT( 37)
#define OPTION_CONSTELLATION            EPHEM_OPTION_BIT( 38)
#define OPTION_RV_AND_DELTA_SIGMAS      EPHEM_OPTION_BIT( 39)

#define ORBIT_SIGMAS_REQUESTED         1
#define NO_ORBIT_SIGMAS_REQUESTED    (-1)

/* Bitfields for 'flags' parameter of OBSERVE */
/*   The OBS_IS_COMET flag is now obsolete */
/* #define OBS_IS_COMET       1   */
#define OBS_DONT_USE       2
   /* Following used for obs from spacecraft that lack offset data */
#define OBS_NO_OFFSET      4
#define OBS_IS_SELECTED    8
#define OBS_NO_VELOCITY    0x10

   /* Following is used for newer NEODyS/AstDyS data for which */
   /* FCCT14 or VFCC17 over-observing correction has already been */
   /* applied;  we shouldn't 'correct' a second time */
#define OBS_ALREADY_CORRECTED_FOR_OVEROBSERVING  0x10


   /* Following flag used temporarily within some functions to note, */
   /* e.g.,  that a particular observation has been processed        */
#define OBS_TEMP_USE_FLAG                        0x20

extern int object_type;

      /* The above should be one of these three values.  Natural */
      /* satellites are folded into the 'asteroid' umbrella :    */
#define OBJECT_TYPE_ASTEROID     0
#define OBJECT_TYPE_COMET        1
#define OBJECT_TYPE_ARTSAT       2

int compare_observations( const void *a, const void *b, void *context);

/* 'Context' can be NULL,  or a pointer to the following integers,  with
other sort orders (perhaps reversed,  or by residuals) possible later. */

#define SORT_OBS_BY_DATE                0
#define SORT_OBS_BY_CODE_THEN_DATE      1
#define SORT_OBS_RADAR_LAST             2

#ifdef SEEK_CUR
OBSERVE FAR *load_observations( FILE *ifile, const char *packed_desig,
                        const int n_obs);
#endif
int unload_observations( OBSERVE FAR *obs, const int n_obs);
OBJECT_INFO *find_objects_in_file( const char *filename,
                                         int *n_found, const char *station);
void sort_object_info( OBJECT_INFO *ids, const int n_ids,
                                          int compare_by_last_obs_time);
int get_object_name( char *obuff, const char *packed_desig);
int get_observer_data( const char FAR *mpc_code, char *buff, mpc_code_t *cinfo);
void recreate_observation_line( char *obuff, const OBSERVE FAR *obs,
                           const int residual_format);   /* ephem0.cpp */
int put_observer_data_in_text( const char FAR *mpc_code, char *buff);

void create_obs_file( const OBSERVE FAR *obs, int n_obs, const int append,
                  const int resid_format);            /* ephem0.cpp */
void create_obs_file_with_computed_values( const OBSERVE FAR *obs,
                  int n_obs, const int append,
                  const int resid_format);            /* ephem0.cpp */
int find_worst_observation( const OBSERVE FAR *obs, const int n_obs);
double calc_absolute_magnitude( OBSERVE FAR *obs, int n_obs);
double compute_rms( const OBSERVE FAR *obs, const int n_obs);
double compute_weighted_rms( const OBSERVE FAR *obs, const int n_obs, int *n_resids);
bool opposition_break( const OBSERVE *obs);              /* elem_out.cpp */
int herget_method( OBSERVE FAR *obs, int n_obs, double r1, double r2,
         double *orbit, double *d_r1, double *d_r2, const char *limited_orbit);
int adjust_herget_results( OBSERVE FAR *obs, int n_obs, double *orbit);
void improve_parabolic( OBSERVE FAR *obs, int n_obs, double *orbit, double epoch);
int full_improvement( OBSERVE FAR *obs, int n_obs, double *orbit,
                 const double epoch, const char *limited_orbit,
                 const int sigmas_requested, const double epoch2);
int set_locs( const double *orbit, const double t0, OBSERVE FAR *obs,
                                   const int n_obs);
void make_date_range_text( char *obuff, const double jd1, const double jd2);
                                                        /* orb_func.cpp */
int get_r1_and_r2( const int n_obs, const OBSERVE FAR *obs,
                             double *r1, double *r2);    /* elem_out.cpp */
int get_idx1_and_idx2( const int n_obs, const OBSERVE FAR *obs,
                                  int *idx1, int *idx2);  /* elem_out.cpp */
double initial_orbit( OBSERVE FAR *obs, int n_obs, double *orbit);
double get_step_size( const char *stepsize, char *step_units,
                                 int *step_digits);          /* ephem0.cpp */
int ephemeris_in_a_file_from_mpc_code( const char *filename,
         const double *orbit,
         OBSERVE *obs, const int n_obs,
         const double epoch_jd, const double jd_start, const char *stepsize,
         const int n_steps, const char *mpc_code,
         ephem_option_t options, const unsigned n_objects);
int find_best_fit_planet( const double jd, const double *ivect,
                     double *rel_vect);     /* runge.cpp */
int integrate_orbit( double *orbit, const double t0, const double t1);
int generate_obs_text( const OBSERVE FAR *obs, const int n_obs, char *buff,
                                          const size_t buffsize);
double convenient_gauss( const OBSERVE FAR *obs, int n_obs, double *orbit,
                  const double mu, const int desired_soln); /* gauss.cpp */
void set_solutions_found( OBJECT_INFO *ids, const int n_ids);
OBSERVE FAR *load_object( FILE *ifile, OBJECT_INFO *id,
                       double *curr_epoch, double *epoch_shown, double *orbit);
int store_solution( const OBSERVE FAR *obs, const int n_obs, const double *orbit,
       const double orbit_epoch, const int perturbers);
int compute_observation_motion_details( const OBSERVE FAR *obs,
               MOTION_DETAILS *m);                    /* mpc_obs.cpp */
int compute_observer_loc( const double jde, const int planet_no,
               const double rho_cos_phi,                    /* mpc_obs.cpp */
               const double rho_sin_phi, const double lon, double FAR *offset);
int compute_observer_vel( const double jde, const int planet_no,
               const double rho_cos_phi,                    /* mpc_obs.cpp */
               const double rho_sin_phi, const double lon, double FAR *offset);
int get_findorb_text( char *buff, const int ival);    /* ephem.cpp */
int write_out_elements_to_file( const double *orbit,
            const double curr_epoch,
            const double epoch_shown,
            OBSERVE FAR *obs, const int n_obs, const char *constraints,
            const int precision, const int monte_carlo,
            const int options);    /* elem_out.cpp */
int extend_orbit_solution( OBSERVE FAR *obs, const int n_obs,
            const double limit, const double time_limit);
int clean_up_find_orb_memory( void);         /* orb_func.cpp */

#define ELEM_OUT_NORMALIZED_MEAN_RESID 16
#define ELEM_OUT_PRECISE_MEAN_RESIDS   8
#define ELEM_OUT_ALTERNATIVE_FORMAT    4
#define ELEM_OUT_NO_COMMENT_DATA       2
/* #define ELEM_OUT_HELIOCENTRIC_ONLY     1     Now obsolete */

#define NO_SIGMAS_AVAILABLE            0
#define COVARIANCE_AVAILABLE           1
#define MONTE_CARLO_SIGMAS_AVAILABLE   2
#define SR_SIGMAS_AVAILABLE            3

#define BLANK_MAG 99.99

/*
   Lowest two bits of the residual_format field:
      resid_format = 0 -> full-line output without tabs;
      resid_format = 1 -> full-line output with tabs;
      resid_format = 2 -> short MPC-like output,  CYY res res form
      resid_format = 3 -> standard 80-column format

   if( resid_format & 4),  four-digit year
   if( resid_format & 8),  residuals are expressed in time and cross-track
         instead of the default residuals in RA and dec
   if( resid_format & 16),  date is expressed as HH:MM:SS (instead of as a
         decimal fraction of a day)
   if( resid_format & 32),  magnitude residuals are shown instead of posn
   if( resid_format & 64),  resids less than an arcsec are shown to 0.001"
   if( resid_format & 128),  very low resids are shown in milli,  micro,
                  nano,  or pico-arcsec.  This is just for checking certain
                  roundoff errors;  you'd not normally do this!
   if( resid_format & 0x100),  residuals are shown with extra digits,  but
        in 'computer-friendly' (no suffixes,  etc.) form.
   if( resid_format & 0x200),  residuals are shown on-screen in RA/dec,
        time/cross-track,  and total angular/magnitude.
*/

#define RESIDUAL_FORMAT_FULL_NO_TABS          0
#define RESIDUAL_FORMAT_FULL_WITH_TABS        1
#define RESIDUAL_FORMAT_SHORT                 2
#define RESIDUAL_FORMAT_80_COL                3
#define RESIDUAL_FORMAT_FOUR_DIGIT_YEARS      4
#define RESIDUAL_FORMAT_TIME_RESIDS           8
         /* 0x10 reserved     */
#define RESIDUAL_FORMAT_MAG_RESIDS         0x20
#define RESIDUAL_FORMAT_PRECISE            0x40
#define RESIDUAL_FORMAT_OVERPRECISE        0x80
#define RESIDUAL_FORMAT_COMPUTER_FRIENDLY  0x100
#define RESIDUAL_FORMAT_EXTRA              0x200
#define RESIDUAL_FORMAT_NORMALIZED         0x400
      /* three bits reserved for time format.
            0 = times in reported format
            1 = times always in decimal days (MPC80 form)
            2 = times always in HHMMSS.sss (ADES form)
            3 = times in MJD
            4-7 = not currently used */
#define RESIDUAL_FORMAT_TIME               0x3800
#define GET_RESID_TIME_FORMAT( residual_format)   ((residual_format >> 11) & 7)

      /* ...and similarly,  three bits reserved for RA format.
            0 = RA/decs in reported format
            1 = RA/decs always in decimal degrees (ADES form)
            2 = RA/decs always in HH.hhhh, dd.dddd
            3 = RA/decs always in HH MM SS.s, dd mm ss.sss (MPC80 form)
            4-7 = not currently used */
#define RESIDUAL_FORMAT_RA_DEC            0x1c000
#define GET_RESID_RA_DEC_FORMAT( residual_format)   ((residual_format >> 14) & 7)
#define RESIDUAL_FORMAT_SHOW_DESIGS              0x20000

int write_residuals_to_file( const char *filename, const char *ast_filename,
        const int n_obs, const OBSERVE FAR *obs_data, const int resid_format);
void format_observation( const OBSERVE FAR *obs, char *text,
                                   const int resid_format);   /* ephem0.cpp */

#define MPC_STATION struct mpc_station

MPC_STATION
   {
   char code[4];
   int color;
   int score;
   };

int find_mpc_color( const MPC_STATION *sdata, const char *mpc_code);
MPC_STATION *find_mpc_color_codes( const int n_obs, const OBSERVE FAR *obs,
                   const int max_n_colors);           /* elem_out.cpp */

int filter_obs( OBSERVE FAR *obs, const int n_obs,           /* orb_fun2.cpp */
                  const double max_residual_in_sigmas, const int filter_type);
   /* Currently,  filter_type = 0 -> in sigmas; = 1 => in arcsec */


    /* Functions used to store and restore orbits for the 'undo' function. */
    /* Orbits are stored on a stack and can be retrieved from it.          */
void push_orbit( const double epoch, const double *orbit); /* orb_fun2.cpp */
int pop_orbit( double *epoch, double *orbit);              /* orb_fun2.cpp */
void pop_all_orbits( void);                                /* orb_fun2.cpp */

typedef struct
{
   double rparam, vparam, orbit[6], score;
} sr_orbit_t;

int find_nth_sr_orbit( sr_orbit_t *orbit, OBSERVE FAR *obs, int n_obs,
                            const int orbit_number);         /* orb_func.cpp */
int get_sr_orbits( sr_orbit_t *orbits, OBSERVE FAR *obs,     /* orb_func.cpp */
               const unsigned n_obs, const unsigned starting_orbit,
               const unsigned max_orbits, const double max_time,
               const double noise_in_sigmas, const int writing_sr_elems);

/* For a 'classical' Herget,  we fit two parameters: R1 and R2,  assuming that they
   have zero residuals.
   For an 'extended' Herget,  it's six parameters: R1, R2,  deltaRAs and delta_decs
   (the residuals for those two observations).
   For an 'adjustment',  four params:  just the deltaRAs and delta_decs.
   For a traditional Vaisala,  no parameters (q/Q is fed),  then 'adjustment' is
      done to get a minimum chi-squared.
   For a parabolic orbit at fixed distance,  four parameters (deltaRAs and delta_decs)
      are fitted,  with the second distance changed using Euler's relationship.  So it
      can be just a tweaked version of the 'adjustment'.
   For a best-fit parabolic orbit,  the distance at first observation is a fifth
      parameter.
   NOTE that not all of these possible fits are implemented yet!

   The last hex digit gives the number of parameters;  preceding digit
   gives the fit type.  For a "full" fit,  four parameters are added:  offsets
   in RA and dec for two observations.  Thus,  FIT_CLASSIC_HERGET fits just two
   parameters:  distances as of two observations.  FIT_HERGET_FULL fits those
   _and_ the four offset parameters (so it ends up as being a plain-vanilla
   unconstrained six-parameter fit.)   */

#define FIT_CLASSIC_HERGET          0x12
#define FIT_HERGET_FULL             0x16

/* You can call find_parameterized_orbit() with params = NULL and
   parameter_type == FIT_NOTHING.  In that case,  we fall through a lot
   of code and end up just getting the orbit connecting obs1 to obs2. */

#define FIT_NOTHING                 0x00

/* The following leaves the distances to the two observations unchanged.
   The result is almost identical to that from the 'linearizing' trick
   I've long applied to Herget orbits to distribute the residuals among
   all observations. */
#define FIT_FIXED_DISTANCES         0x24

/* A plain vanilla FIT_VAISALA has one parameter,  q/Q,  and the assumption
   that the two observations are "perfect".  FIT_VAISALA_FULL also adjusts
   those two observations,  adding four more parameters as described above,
   to distribute the residual error among all observations. */

#define FIT_VAISALA                 0x31
#define FIT_VAISALA_FULL            0x35

/* There are two possible parabolic orbit fits for a given distance,
one headed "toward" you and one "away" from you.  (I know this is always
the case when the orbital arc is far enough from the sun,  but haven't
been able to mathematically prove (yet) that it is true near the sun...
which may turn out to matter for sungrazers.  Also note that I'm only
considering orbits where the transfer orbit between the two observations
does not include the sun:  i.e.,  the change in true anomaly is between
-180 and +180 degrees.)  */

#define FIT_PARABOLIC1              0x41
#define FIT_PARABOLIC1_FULL         0x45
#define FIT_PARABOLIC2              0x51
#define FIT_PARABOLIC2_FULL         0x55

/* If you're fitting a parabolic orbit,  but with the initial distance
fixed,  you have one less parameter.  (In the "non-full" cases,  you
have no parameters to fit at all:  the distance to the first observation,
and the direction to the second,  completely specify the orbit.) */

#define FIT_PARABOLIC1_FIXED        0x60
#define FIT_PARABOLIC1_FIXED_FULL   0x65
#define FIT_PARABOLIC2_FIXED        0x70
#define FIT_PARABOLIC2_FIXED_FULL   0x74

/* Dunno if we'll actually do this,  but... we _could_ fit circular orbits.
To do so,  we'd apply the four RA/dec offset parameters,  then find the
circular orbit connecting the two resulting adjusted observations.  (Or
some other set of four parameters.)    */

// #define FIT_CIRCULAR_ORBIT          0x?4

/* For the following constants,  bits 0-3 give the number of 'additional'
non-grav parameters,  and bit 4 indicates if the force is an inverse
square one. Higher bits can distinguish,  for example,  SRP from
Yarkovsky with A2 (both inverse-square and both one added parameter). */

#define FORCE_MODEL_NO_NONGRAVS        0
#define FORCE_MODEL_SRP                0x11
#define FORCE_MODEL_SRP_TWO_PARAM      0x12
#define FORCE_MODEL_SRP_THREE_PARAM    0x13
#define FORCE_MODEL_COMET_TWO_PARAM    0x02
#define FORCE_MODEL_COMET_THREE_PARAM  0x03
#define FORCE_MODEL_COMET_FOUR_PARAM   0x04

/* If we think there was a delta-v,  due to spacecraft maneuver or impact
or something else,  it can be described in four parameters : three for
the amount of the delta-V,  and a fourth saying when it happened.  At
some point,  there will probably be a model with Yet Another Parameter
so that we can solve for the circumstances of the delta-V _and_ the
object's area/mass ratio at the same time.  This 'Delta-V plus SRP'
model has not yet been implemented.    */

#define FORCE_MODEL_DELTA_V            0x204
#define FORCE_MODEL_DELTA_V_SRP        0x205

/* For some rocks,  Yarkovsky can be modelled as an A2 (along-orbit)
inverse square force.  */

#define FORCE_MODEL_YARKO_A2           0x111

extern int force_model;

bool is_inverse_square_force_model( void);

const char *find_orb_version_jd( double *jd);

      /* In the console version of Find_Orb,  the following two functions */
      /* get remapped to Curses functions.  In the non-interactive one,   */
      /* they're mapped to 'do-nothings'.  See fo.cpp & find_orb.cpp.     */
void refresh_console( void);
void move_add_nstr( const int col, const int row, const char *msg, const int n_bytes);

#define COLOR_DEFAULT_INQUIRY       12
#define COLOR_ATTENTION             13

int inquire( const char *prompt, char *buff, const int max_len,
                     const int color);
