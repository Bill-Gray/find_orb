/* details.h : code to store/access 'observational details' (header)
data for 80-column MPC-formatted astrometry,  as documented at

http://www.minorplanetcenter.net/iau/info/ObsDetails.html */

void *init_observation_details( void);
int add_line_to_observation_details( void *obs_details, const char *iline);
const char **get_code_details( const void *obs_details, const char *mpc_code);
void free_observation_details( void *obs_details);

/* Return codes for add_line_to_observation_details : */

#define OBS_DETAILS_IRRELEVANT_LINE             0
#define OBS_DETAILS_MPC_80_COLUMN_LINE          1
#define OBS_DETAILS_HEADER_LINE                 2
