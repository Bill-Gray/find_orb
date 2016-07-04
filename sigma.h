int load_up_sigma_records( const char *filename);
void free_sigma_recs( void);
double get_observation_sigma( const double jd, const int mag_in_tenths,
                  const char *mpc_code, double *mag_sigma,
                  double *time_sigma, const char program_code);
