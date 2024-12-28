/* Miscellaneous constants used throughout Find_Orb */

#define J2000 2451545.
#define JD_TO_YEAR( jd)  (2000. + ((jd) - J2000) / 365.25)
#define YEAR_TO_JD( year) (J2000 + (year - 2000.) * 365.25)
#define METERS_PER_KM      1000.

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078

#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
      /* GAUSS_K is a fixed constant. SOLAR_GM = 0.0002959122082855911025, */
      /* in AU^3/day^2; = 1.3271243994E+11 km3/s2 */

#define EARTH_RADIUS_IN_KM    6378.140
#define SUN_RADIUS_IN_KM      695700.
#define LUNAR_RADIUS_IN_KM    1737.4
#define EARTH_RADIUS_IN_AU    (EARTH_RADIUS_IN_KM / AU_IN_KM)
#define SUN_RADIUS_IN_AU      (SUN_RADIUS_IN_KM / AU_IN_KM)
#define LUNAR_RADIUS_IN_AU    (LUNAR_RADIUS_IN_KM / AU_IN_KM)

/* "Solar constant" : at one AU from the sun,  a square meter
collects 1361 watts (total solar irradiance at all wavelengths). */

#define SOLAR_CONSTANT     1361.       /* W/m^2 = kg/s^3 at one AU */

const double SRP1AU = SOLAR_CONSTANT * seconds_per_day * seconds_per_day
                            / (AU_IN_METERS * SPEED_OF_LIGHT * METERS_PER_KM);

/* "Solar radiation pressure at 1 AU",  in kg*AU^3 / (m^2*d^2),
from a private communication from Steve Chesley.  This means
that if you had a one-kilogram object showing one square
meter of surface area to the sun,  and it was one AU from the
sun,  and it absorbed all the solar radiation (and re-radiated
it isotropically;  i.e.,  the re-radiated photons didn't cause
any net force),  then that object would accelerate away from
the sun at about 2.3e-7 AU/day^2 (a.k.a. 4.6e-6 m/s^2).

   Reality will not be so precise;  real objects reflect and
re-radiate photons,  instead of politely absorbing them and
then re-radiating them isotropically.

    A final comment : non-gravs of the A1, A2, A3 form give the
acceleration the object would have at one AU from the sun,  in
units of AU/day.  If the non-gravs are of the 1/r^2 model (rock-like)
rather than,  say,  the Sekanina-Marsden model for an outgassing
comet,  then for an area/mass ratio of z m^2/kg,  A1 = SRP1AU * z.
Or,  alternatively,  A1 = SRP1AU * AMR.  (Roughly speaking,  anyway;
in such scenarios,  A2 and A3 are fitted parameters,  so the radial
component is slightly different.)  Thus,  for example,  one gets an
area/mass ratio for 1I/`Oumuamua of 1.15 m^2/kg,  but if you fit A1
and A2,  you get A1 = 2.39e-7 AU/day^2... somewhat lower than you'd
expect if you just multiplied 1.15 by SRP1AU,  but pretty close;  the
tangential components are rarely large.  */
