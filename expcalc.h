typedef struct
{
   const char *mpc_code;
   char filter;
   double primary_diam, obstruction_diam;   /* in cm */
   double aperture, fwhm;                   /* arcsec */
   double qe;                               /* unitless,  0-1 */
   double readnoise;                        /* counts per pixel */
   double pixel_size;                       /* arcsec */
   double sky_brightness;                   /* magnitudes/arcsec^2 */
   double airmass;                          /* ~ 1/sin(alt) */
} expcalc_config_t;

int find_expcalc_config_from_mpc_code( const char *mpc_code, expcalc_config_t *c);
double mag_from_snr_and_exposure( const expcalc_config_t *c,
                              const double snr, const double exposure);
double snr_from_mag_and_exposure( const expcalc_config_t *c,
                              const double mag, const double exposure);
double exposure_from_snr_and_mag( const expcalc_config_t *c,
                              const double snr, const double mag);
