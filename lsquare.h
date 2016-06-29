void *lsquare_init( const int n_params);
int lsquare_add_observation( void *lsquare, const double residual,
                                    const double weight, const double *obs);
int lsquare_solve( const void *lsquare, double *result);
void lsquare_free( void *lsquare);
double *lsquare_covariance_matrix( const void *lsquare);
double *lsquare_wtw_matrix( const void *lsquare);
