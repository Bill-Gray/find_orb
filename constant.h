/* Miscellaneous constants used throughout Find_Orb */

#define J2000 2451545.
#define JD_TO_YEAR( jd)  (2000. + ((jd) - J2000) / 365.25)
#define YEAR_TO_JD( year) (J2000 + (year - 2000.) * 365.25)

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078

#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
      /* GAUSS_K is a fixed constant. SOLAR_GM = 0.0002959122082855911025, */
      /* in AU^3/day^2; = 1.3271243994E+11 km3/s2 */

#define EARTH_RADIUS_IN_KM    6378.140
#define SUN_RADIUS_IN_KM      695700.
#define EARTH_RADIUS_IN_AU    (EARTH_RADIUS_IN_KM / AU_IN_KM)
#define SUN_RADIUS_IN_AU      (SUN_RADIUS_IN_KM / AU_IN_KM)

const double SRP1AU = 2.3e-7;

/* "Solar radiation pressure at 1 AU",  in kg*AU^3 / (m^2*d^2),
from a private communication from Steve Chesley.  This means
that if you had a one-kilogram object showing one square
meter of surface area to the sun,  and it was one AU from the
sun,  and it absorbed all the solar radiation (and re-radiated
it isotropically;  i.e.,  the re-radiated photons didn't cause
any net force),  then that object would accelerate away from
the sun at 2.3e-7 AU/day^2 (a.k.a. 4.562e-6 m/s^2).

  One can derive SRP1AU  from basic principles.   The 'solar
constant' is C = 1367.6 AU^2*W/m^2; i.e.,  if you set up a
one square meter solar panel with 100% efficiency one AU from
the Sun, pointed straight at the sun, it would generate
1367.6 watts. To get the above number, one uses

SRP1AU = C * d^2 / (meters_per_AU * c)

   C = 1367.6 AU^2*W/m^2 = 1367.6 AU^2*kg/s^3
   d = 86400 seconds/day
   meters_per_AU = 1495978707 m/AU
   c = 299792458 m/s

   The result indicates that if the solar panel in question
had a mass of one kilogram,  it would accelerate away from the
sun at 2.27636e-7 AU/day^2.  I think Steve rounded off with a
one-percent error because that's a decent match to the accuracy
you can expect with real-world objects that reflect and re-radiate
photons,  instead of politely absorbing them and then re-radiating
them isotropically.

  Interestingly,  this also means that an object with area/mass =
1286 m^2/kg would have SRP balancing the sun's gravity.  Which
would make it big and light.  Solar sails aren't easy.

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
