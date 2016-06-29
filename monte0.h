#define MONTE_TP           0
#define MONTE_ECC          1
#define MONTE_q            2
#define MONTE_Q            3
#define MONTE_INV_A        4
#define MONTE_INCL         5
#define MONTE_MEAN_ANOM    6
#define MONTE_ARG_PER      7
#define MONTE_ASC_NODE     8

#define MONTE_N_ENTRIES    9
#define MONTE_DATA_SIZE   (3 * MONTE_N_ENTRIES)

double *add_gaussian_noise_to_obs( int n_obs, OBSERVE *obs,
                 const double noise_in_sigmas);             /* monte0.cpp */
void add_monte_orbit( double *monte_data, const ELEMENTS *elem,
                  const int n_orbits);                      /* monte0.cpp */
void compute_monte_sigmas( double *sigmas, const double *monte_data,
                  const int n_orbits);                      /* monte0.cpp */
void restore_ra_decs_mags_times( unsigned n_obs, OBSERVE *obs,
                           const double *stored_ra_decs);
void put_orbital_elements_in_array_form( const ELEMENTS *elem,
                  double *output_array);                    /* monte0.cpp */
double dump_monte_data_to_file( FILE *ofile, const double *sigmas,
            const double semimajor_axis, const double ecc,
            const int planet_orbiting);                       /* monte0.cpp */
char * put_double_in_buff( char *buff, const double ival);    /* monte0.cpp */

