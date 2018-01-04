#ifndef MPC_FMT_H_INCLUDED
#define MPC_FMT_H_INCLUDED

bool is_valid_mpc_code( const char *mpc_code);        /* mpc_fmt.cpp */
double extract_date_from_mpc_report( const char *buff, unsigned *format);
int get_ra_dec_from_mpc_report( const char *ibuff,    /* mpc_fmt.cpp */
                       int *ra_format, double *ra, double *ra_precision,
                       int *dec_format, double *dec, double *dec_precision);
#endif
