/* sm_vsop.cpp: _very_ rough ephems for planets and earth's moon

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

#include <math.h>
#include <stdlib.h>
#include "watdefs.h"
#include "afuncs.h"     /* for vector3_length() proto */

/* Following are the first few terms of the VSOP87A Cartesian coordinate
   series for each planet (Mercury through Neptune;  VSOP doesn't do
   Pluto).  Data can be found at https://cdsarc.cds.unistra.fr/viz-bin/cat/VI/81.

Planetary Theories in rectangular and spherical variables: VSOP87 solution.
     Bretagnon P., Francou G.
    <Astron. Astrophys. 202, 309 (1988)> =1988A&A...202..309B

Theorie du mouvement de l'ensemble des planetes (VSOP82).
     Bretagnon P.
    <Astron. Astrophys. 114, 278 (1982)> =1982A&A...114..278B

   The 'axis' is 1 for x, 2 for y, 3 for z.

   The position is computed as the sum of terms of the form

   amplitude * cos( phase + freq * t_cen)

   where t_cen = (JD - 2451545) / 36525 = centuries from J2000.  Except that
   the first term,  with axis = 0,  just includes constant offsets:  the
   'amplitude' = xoffset, 'phase' = yoffset, 'freq' = zoffset.

   Coordinates are J2000 ecliptic.

   All this is used to get _very_ rough planetary positions for use in
   deciding if any perturbers are needed for the initial orbit determination
   (basically,  "is this coordinate  within a few Hill spheres of a planet").
   We can usually use this to say:  "No planet is conceivably close enough
   to exert serious short-term perturbations right now",  saving us any
   really computationally intensive effort (most of the time),  even the
   effort of finding more accurate planetary positions.  The full series
   contain thousands of terms,  including some where the terms are multiplied
   by t_cen ^ n (n=1...5).

   Terms of the following magnitude or greater are included:
      Mercury     .006 AU
      Venus       .003
      Earth       .0001 *
      Mars        .01
      Jupiter     .02
      Saturn      .025
      Uranus      .02
      Neptune     .03

   * Present,  but commented out,  are terms for the Earth down to .00001 AU.
   But you'd only get that level of accuracy within a decade or so of J2000;
   beyond that,  some secular terms would be needed.

   Code is also included for approximate lunar and Plutonian ephemerides.  */

int compute_rough_planet_loc( const double t_cen, const int planet_idx,
                                          double *vect);    /* sm_vsop.c */
int check_for_perturbers( const double t_cen, const double *vect); /* sm_vsop*/

/* #define LONG_EARTH_SERIES  1  */

#define SMALL_VSOP_TERM struct small_vsop_term

SMALL_VSOP_TERM
   {
   int axis;
   double amplitude, phase, freq;
   };

#define N_MERC_TERMS 6
#define N_VENU_TERMS 4
#ifdef LONG_EARTH_SERIES
   #define N_EART_TERMS 21
#else
   #define N_EART_TERMS 7
#endif
#define N_MARS_TERMS 6
#define N_JUPI_TERMS 6
#define N_SATU_TERMS 10
#define N_URAN_TERMS 12
#define N_NEPT_TERMS 10

static const int n_terms[9] = { 0, N_MERC_TERMS, N_VENU_TERMS, N_EART_TERMS,
       N_MARS_TERMS, N_JUPI_TERMS, N_SATU_TERMS, N_URAN_TERMS, N_NEPT_TERMS };

static const SMALL_VSOP_TERM merc_terms[N_MERC_TERMS] = {
 { 0, -0.02625615963, -0.11626131831, -0.00708734365 },
 { 1,  0.37546291728, 4.39651506942,   2608.790314157420  },
 { 1,  0.03825746672, 1.16485604339,   5217.580628314840  },
 { 2,  0.37953642888, 2.83780617820,   2608.790314157420  },
 { 2,  0.03854668215, 5.88780608966,   5217.580628314840  },
 { 3,  0.04607665326, 1.99295081967,   2608.790314157420  } };

static const SMALL_VSOP_TERM venu_terms[N_VENU_TERMS] = {
 { 0,  0.00486448018, -0.00549506273, -0.00035588343 },
 { 1,  0.72211281391, 3.17575836361,   1021.328554621100  },
 { 2,  0.72324820731, 1.60573808356,   1021.328554621100  },
 { 3,  0.04282990302, 0.26703856476,   1021.328554621100  } };

#ifdef LONG_EARTH_SERIES
         /* Might come in handy someday,  for computing basic orbits without */
         /* PS-1996 or DE ephems.  At this level,  though,  we may have to   */
         /* include some of the currently ignored secular terms.             */
static const SMALL_VSOP_TERM eart_terms[N_EART_TERMS] = {
   { 0,  0.00561144206, -0.02442699036, -0.00000001086 },
   { 1,  0.99982928844, 1.75348568475,    628.307584999140  },
   { 1,  0.00835257300, 1.71034539450,   1256.615169998280  },
   { 1,  0.00010466628, 1.66722645223,   1884.922754997420  },
   { 1,  0.00003110838, 0.66875185215,   8399.684731811189  },
   { 1,  0.00002552498, 0.58310207301,     52.969096509460  },
   { 1,  0.00002137256, 1.09235189672,    157.734354244780  },
   { 1,  0.00001709103, 0.49540223397,    627.955273164240  },
   { 1,  0.00001707882, 6.15315547484,    628.659896834040  },
   { 1,  0.00001445242, 3.47272783760,    235.286615377180  },
   { 1,  0.00001091006, 3.68984782465,    522.369391980220  },
   { 2,  0.99989211030, 0.18265890456,    628.307584999140  },
   { 2,  0.00835292314, 0.13952878991,   1256.615169998280  },
   { 2,  0.00010466965, 0.09641690558,   1884.922754997420  },
   { 2,  0.00003110838, 5.38114091484,   8399.684731811189  },
   { 2,  0.00002570338, 5.30103973360,     52.969096509460  },
   { 2,  0.00002147473, 2.66253538905,    157.734354244780  },
   { 2,  0.00001709219, 5.20780401071,    627.955273164240  },
   { 2,  0.00001707987, 4.58232858766,    628.659896834040  },
   { 2,  0.00001440265, 1.90068164664,    235.286615377180  },
   { 2,  0.00001135092, 5.27313415220,    522.369391980220  } };
#else
static const SMALL_VSOP_TERM eart_terms[N_EART_TERMS] = {
 { 0,  0.00561144206, -0.02442699036, -0.00000001086 },
 { 1,  0.99982928844, 1.75348568475,    628.307584999140  },
 { 1,  0.00835257300, 1.71034539450,   1256.615169998280  },
 { 1,  0.00010466628, 1.66722645223,   1884.922754997420  },
 { 2,  0.99989211030, 0.18265890456,    628.307584999140  },
 { 2,  0.00835292314, 0.13952878991,   1256.615169998280  },
 { 2,  0.00010466965, 0.09641690558,   1884.922754997420  } };
#endif

#ifdef INCLUDE_EARTH_MOON_BARYCENTER
         /* no real use for these,  except for reference: */
static const SMALL_VSOP_TERM emb_terms[5] = {
 { 0,  0.00561144161, -0.02442698841, -0.00000001086 },
 { 1,  0.99982927460, 1.75348568475,    628.307584999140  },
 { 1,  0.00835257300, 1.71034539450,   1256.615169998280  },
 { 2,  0.99989209645, 0.18265890456,    628.307584999140  },
 { 2,  0.00835292314, 0.13952878991,   1256.615169998280  } };
#endif

static const SMALL_VSOP_TERM mars_terms[N_MARS_TERMS] = {
 { 0, -0.19502945246,  0.08655481102,  0.00660669541 },
 { 1,  1.51769936383, 6.20403346548,    334.061242669980  },
 { 1,  0.07070919655, 0.25870338558,    668.122485339960  },
 { 2,  1.51558976277, 4.63212206588,    334.061242669980  },
 { 2,  0.07064550239, 4.97051892902,    668.122485339960  },
 { 3,  0.04901207220, 3.76712324286,    334.061242669980  } };

static const SMALL_VSOP_TERM jupi_terms[N_JUPI_TERMS] = {
 { 0, -0.36662642320, -0.09363670616,  0.00859031952 },
 { 1,  5.19663470114, 0.59945082355,     52.969096509460  },
 { 1,  0.12593937922, 0.94911583701,    105.938193018920  },
 { 2,  5.19520046589, 5.31203162731,     52.969096509460  },
 { 2,  0.12592862602, 5.66160227728,    105.938193018920  },
 { 3,  0.11823100489, 3.55844646343,     52.969096509460  } };

static const SMALL_VSOP_TERM satu_terms[N_SATU_TERMS] = {
 { 0,  0.04244797817, -0.79387988806,  0.01214249867 },
 { 1,  9.51638335797, 0.87441380794,     21.329909543800  },
 { 1,  0.26412374238, 0.12390892620,     42.659819087600  },
 { 1,  0.06760430339, 4.16767145778,     20.618554843720  },
 { 1,  0.06624260115, 0.75094737780,     22.041264243880  },
 { 2,  9.52986882699, 5.58600556665,     21.329909543800  },
 { 2,  0.26441781302, 4.83528061849,     42.659819087600  },
 { 2,  0.06916653915, 2.55279408706,     20.618554843720  },
 { 2,  0.06633570703, 5.46258848288,     22.041264243880  },
 { 3,  0.41356950940, 3.60234142982,     21.329909543800  } };

static const SMALL_VSOP_TERM uran_terms[N_URAN_TERMS] = {
 { 0,  1.32272523872, -0.16256125476, -0.01774318778 },
 { 1, 19.17370730359, 5.48133416489,      7.478159856730  },
 { 1,  0.44402496796, 1.65967519586,     14.956319713460  },
 { 1,  0.14668209481, 3.42395862804,      7.329712585900  },
 { 1,  0.14130269479, 4.39572927934,      7.626607127560  },
 { 1,  0.06201106178, 5.14043574125,       .148447270830  },
 { 2, 19.16518231584, 3.91045677002,      7.478159856730  },
 { 2,  0.44390465203, 0.08884111329,     14.956319713460  },
 { 2,  0.14755940186, 1.85423280679,      7.329712585900  },
 { 2,  0.14123958128, 2.82486076549,      7.626607127560  },
 { 2,  0.06250078231, 3.56960243857,       .148447270830  },
 { 3,  0.25878127698, 2.61861272578,      7.478159856730  } };

static const SMALL_VSOP_TERM nept_terms[N_NEPT_TERMS] = {
 { 0, -0.27080164222, -0.30205857683,  0.01245978462 },
 { 1, 30.05890004476, 5.31211340029,      3.813303563780  },
 { 1,  0.13505661755, 3.50078975634,      7.626607127560  },
 { 1,  0.15726094556, 0.11319072675,      3.664856292950  },
 { 1,  0.14935120126, 1.08499403018,      3.961750834610  },
 { 2, 30.06056351665, 3.74086294714,      3.813303563780  },
 { 2,  0.13506391797, 1.92953034883,      7.626607127560  },
 { 2,  0.15706589373, 4.82539970129,      3.664856292950  },
 { 2,  0.14936165806, 5.79694900665,      3.961750834610  },
 { 3,  0.92866054405, 1.44103930278,      3.813303563780  } };

static const SMALL_VSOP_TERM *vsop_terms[] = { NULL,
         merc_terms, venu_terms, eart_terms, mars_terms,
         jupi_terms, satu_terms, uran_terms, nept_terms };

      /* See comments below on Meeus' ' _Astronomical Algorithms_. */
static const double lunar_lon_r_terms[] = {
 0.10975981, -0.0001397435,   8328.69142695, 2.35555590, /*  0  0  1  0 */
 0.02223597, -0.0000247270,   7214.06286667, 1.75819228, /*  2  0 -1  0 */
 0.01148975, -0.0000197594,  15542.75429363, 4.11374817, /*  2  0  0  0 */
 0.00372834, -0.0000038097,  16657.38285391, 4.71111180, /*  0  0  2  0 */
-0.00323088,  0.0000003268,    628.30195517, 6.24006013, /*  0  1  0  0 */
-0.00199547, -0.0000000210,  16866.93231626, 3.25581047, /*  0  0  0  2 */
#ifdef SMALL_TERMS_LISTED_FOR_REFERENCE_ONLY
         /* You _could_ remove the above #ifdef and following #endif and */
         /* include these smaller terms.  But it's probably unnecessary; */
         /* use JPL ephemerides to get still better accuracy with faster */
         /* computations to boot.  The extra terms are included just in  */
         /* case I come up with a good use for them someday.  This very  */
         /* truncated series is enough to answer the question:  "Is the  */
         /* Moon close enough to be a significant perturber?"            */
 0.00102613,  0.0000016455,  -1114.62856028, 5.68582168, /*  2  0 -2  0 */
 0.00099599, -0.0000010170,   6585.76091150, 1.80131746, /*  2 -1 -1  0 */
 0.00093064, -0.0000011413,  23871.44572058, 0.18611877, /*  2  0  1  0 */
 0.00079863, -0.0000013676,  14914.45233846, 4.15687335, /*  2 -1  0  0 */
-0.00071424, -0.0000008665,  -7700.38947179, 3.88450423, /*  0  1 -1  0 */
-0.00060598,  0.0000007269,   7771.37714681, 5.19846674, /*  1  0  0  0 */
-0.00053028,  0.0000007002,   8956.99338212, 2.31243072, /*  0  1  1  0 */
 0.00026751,  0.0000000690,  -1324.17802264, 0.85793771, /*  2  0  0 -2 */
-0.00021865,  0.0000000000,  25195.62374322, 5.61136636, /*  0  0  1  2 */
 0.00019164,  0.0000005325,  -8538.24088931, 5.38293074, /*  0  0  1 -2 */
 0.00018631, -0.0000002325,  22756.81716030, 5.87194045, /*  4  0 -1  0 */
 0.00017513, -0.0000001551,  24986.07428086, 0.78348239, /*  0  0  3  0 */
 0.00014919, -0.0000001446,  14428.12573334, 3.51638455, /*  4  0 -2  0 */
-0.00013767,  0.0000001618,   7842.36482184, 1.71506710, /*  2  1 -1  0 */
-0.00011809,  0.0000002060,  16171.05624879, 4.07062299, /*  2  1  0  0 */
-0.00009011, -0.0000000560,   -557.31428014, 2.84291084, /*  1  0 -1  0 */
 0.00008704, -0.0000001115,   8399.67910198, 5.15534156, /*  1  1  0  0 */
 0.00007044, -0.0000000858,  23243.14376541, 0.22924395, /*  2 -1  1  0 */
 0.00006971, -0.0000000698,  32200.13714754, 2.54167467, /*  2  0  2  0 */
 0.00006739, -0.0000000779,  31085.50858725, 1.94431104, /*  4  0  0  0 */
 0.00006397,  0.0000000963,  -9443.31998724, 3.33026579, /*  2  0 -3  0 */
-0.00004693, -0.0000000468, -16029.08089874, 1.52894833, /*  0  1 -2  0 */
-0.00004541,  0.0000000000,  24080.99518293, 5.01400274, /*  2  0 -1  2 */
 0.00004171,  0.0000000672,  -1742.93051545, 5.72894686, /*  2 -1 -2  0 */
-0.00004098,  0.0000000423,  16100.06857377, 1.27083733, /*  1  0  1  0 */
 0.00003903, -0.0000000661,  14286.15038329, 4.19999853, /*  2 -2  0  0 */
-0.00003700,  0.0000000384,  17285.68480908, 4.66798662, /*  0  1  2  0 */
-0.00003611,  0.0000000000,   1256.60391033, 6.19693495, /*  0  2  0  0 */
 0.00003574, -0.0000000331,   5957.45895634, 1.84444264, /*  2 -2 -1  0 */
-0.00003094,  0.0000000276,   7004.51340432, 3.21349361, /*  2  0  1 -2 */
-0.00002784,  0.0000000000,  32409.68660989, 1.08637333, /*  2  0  0  2 */
 0.00002121, -0.0000000265,  22128.51520513, 5.91506563, /*  4 -1 -1  0 */
-0.00001937,  0.0000000000,  33524.31517017, 1.68373696, /*  0  0  2  2 */
-0.00001557,  0.0000000218,  14985.44001348, 0.67347371, /*  3  0 -1  0 */
-0.00001414,  0.0000000175,  24499.74767575, 0.14299359, /*  2  1  1  0 */
 0.00001325, -0.0000000127,  13799.82377817, 3.55950973, /*  4 -1 -2  0 */
-0.00001244, -0.0000000142,  -7072.08751662, 3.84137905, /*  0  2 -1  0 */
-0.00001222,  0.0000000157,   8470.66677701, 1.67194192, /*  2  2 -1  0 */
 0.00001206,  0.0000000000,   -486.32660512, 5.64269650, /*  2  1 -2  0 */
 0.00001040,  0.0000000000,  -1952.47997780, 0.90106289, /*  2 -1  0 -2 */
 0.00000958, -0.0000000095,  39414.20001421, 4.29986694, /*  4  0  1  0 */
 0.00000937, -0.0000000075,  33314.76570782, 3.13903829, /*  0  0  4  0 */
 0.00000908, -0.0000000105,  30457.20663208, 1.98743622, /*  4 -1  0  0 */
-0.00000850, -0.0000000116,  -8886.00570710, 0.48735494, /*  1  0 -2  0 */
-0.00000696,  0.0000000000,   -695.87606747, 0.81481253, /*  2  1  0 -2 */
-0.00000665, -0.0000000296,   -209.54946235, 1.45530133, /*  0  0  2 -2 */
 0.00000613,  0.0000000000,  16728.37052893, 1.22771215, /*  1  1  1  0 */
-0.00000593,  0.0000000000,   6656.74858653, 4.60110312, /*  3  0 -2  0 */
 0.00000576,  0.0000000000,   6099.43430639, 1.16082865, /*  4  0 -3  0 */
 0.00000571,  0.0000000000,  31571.83519237, 2.58479984, /*  2 -1  2  0 */
-0.00000564,  0.0000000078,   9585.29533729, 2.26930554, /*  0  2  1  0 */
 0.00000522,  0.0000000000,     70.98767503, 2.79978566, /*  1  1 -1  0 */
 0.00000513,  0.0000000000,  40528.82857449, 4.89723056, /*  2  0  3  0 */
 0.00000000,  0.0000000585,  -9652.86944959, 4.78556712, /*  2  0 -1 -2 */
#endif /* SMALL_TERMS_LISTED_FOR_REFERENCE_ONLY */
0 };

static const double lunar_lat_terms[] = {
 0.08950261,   8433.46615813, 1.62790523, /*  0  0  0  1 */
 0.00489743,  16762.15758509, 3.98346113, /*  0  0  1  1 */
 0.00484666,   -104.77473118, 0.72765067, /*  0  0  1 -1 */
 0.00302356,   7109.28813550, 2.48584294, /*  2  0  0 -1 */
 0.00096714,  15647.52902480, 3.38609751, /*  2  0 -1  1 */
 0.00080758,  -1219.40329146, 0.13028704, /*  2  0 -1 -1 */
#ifdef SMALL_TERMS_LISTED_FOR_REFERENCE_ONLY
         /* See comments above. */
 0.00056851,  23976.22045176, 5.74165341, /*  2  0  0  1 */
 0.00030016,  25090.84901204, 0.05583172, /*  0  0  2  1 */
 0.00016172,  15437.97956245, 4.84139884, /*  2  0  1 -1 */
 0.00015397,   8223.91669578, 3.08320656, /*  0  0  2 -1 */
 0.00014340,   6480.98618033, 2.52896812, /*  2 -1  0 -1 */
 0.00007547,  -9548.09471841, 4.05791645, /*  2  0 -2 -1 */
 0.00007330,  32304.91187871, 1.81402400, /*  2  0  1  1 */
-0.00005863,   7737.59009066, 2.44271776, /*  2  1  0 -1 */
 0.00004299,  15019.22706963, 3.42922269, /*  2 -1 -1  1 */
 0.00003859,  23347.91849659, 5.78477859, /*  2 -1  0  1 */
 0.00003604,  -1847.70524663, 0.17341222, /*  2 -1 -1 -1 */
-0.00003264, -16133.85562992, 2.25659900, /*  0  1 -1 -1 */
 0.00003190,  14323.35100217, 4.24403522, /*  4  0 -1 -1 */
-0.00003131,   9061.76811330, 1.58478005, /*  0  1  0  1 */
-0.00003053,  25300.39847439, 4.88371570, /*  0  0  0  3 */
-0.00002731,    733.07668634, 5.51240946, /*  0  1 -1  1 */
-0.00002602,  16204.84330494, 0.54318667, /*  1  0  0  1 */
-0.00002574,  17390.45954025, 3.94033595, /*  0  1  1  1 */
-0.00002461,    523.52722399, 0.68452549, /*  0  1  1 -1 */
-0.00002346,  -7805.16420296, 4.61215489, /*  0  1  0 -1 */
-0.00002330,   -662.08901132, 3.57056151, /*  1  0  0 -1 */
 0.00001932,  33419.54043900, 2.41138762, /*  0  0  3  1 */
 0.00001782,  22652.04242912, 0.31640581, /*  4  0  0 -1 */
 0.00001454,  31190.28331843, 1.21666038, /*  4  0 -1  1 */
 0.00001356, -16971.70704744, 3.75502551, /*  0  0  1 -3 */
 0.00001171,  22861.59189147, 5.14428979, /*  4  0 -2  1 */
 0.00001059,  -9757.64418077, 5.51321778, /*  2  0  0 -3 */
 0.00001040,  23766.67098940, 0.91376943, /*  2  0  2 -1 */
 0.00000857,  14809.67760728, 4.88452402, /*  2 -1  1 -1 */
-0.00000787,   7318.83759785, 1.03054161, /*  2  0 -2  1 */
 0.00000766,  16552.60812273, 5.43876246, /*  0  0  3 -1 */
 0.00000737,  40633.60330567, 4.16957990, /*  2  0  2  1 */
 0.00000735, -17876.78614537, 1.70236055, /*  2  0 -3 -1 */
-0.00000639,  16275.83097997, 3.34297233, /*  2  1 -1  1 */
-0.00000613,  24604.52240692, 5.69852823, /*  2  1  0  1 */
 0.00000578,  39518.97474538, 3.57221628, /*  4  0  0  1 */
 0.00000550,  31676.60992354, 1.85714918, /*  2 -1  1  1 */
 0.00000527,   5852.68422516, 2.57209330, /*  2 -2  0 -1 */
-0.00000494,  33629.08990135, 0.95608629, /*  0  0  1  3 */
-0.00000400,  16066.28151762, 4.79827366, /*  2  1  1 -1 */
 0.00000389,    -33.78705615, 3.52743633, /*  1  1  0 -1 */
 0.00000389,  16833.14526011, 0.50006149, /*  1  1  0  1 */
-0.00000384, -24462.54705687, 6.18422840, /*  0  1 -2 -1 */
-0.00000384,  -8919.79276325, 4.01479127, /*  2  1 -2 -1 */
-0.00000323,  24533.53473190, 2.89874257, /*  1  0  1  1 */
 0.00000316, -10176.39667358, 4.10104163, /*  2 -1 -2 -1 */
-0.00000309,  25719.15096721, 0.01270654, /*  0  1  2  1 */
 0.00000307,   5994.65957521, 1.88847932, /*  4  0 -2 -1 */
 0.00000290,  13695.04904700, 4.28716040, /*  4 -1 -1 -1 */
-0.00000286,   7666.60241564, 5.92611741, /*  1  0  1 -1 */
 0.00000230,  30980.73385608, 2.67196171, /*  4  0  1 -1 */
-0.00000208,  -8990.78043827, 1.21500561, /*  1  0 -1 -1 */
 0.00000201,  22023.74047395, 0.35953099, /*  4 -1  0 -1 */
 0.00000187,  22719.61654142, 5.82790377, /*  2 -2  0  1 */
#endif /* SMALL_TERMS_LISTED_FOR_REFERENCE_ONLY */
0. };

      /* Uses Meeus' _Astronomical Algorithms_, chapter 45,  'Position of */
      /* the Moon'.  The first few terms give a position to within about */
      /* 500 km.  If you want more than that,  you should probably be using */
      /* ELP-82 or JPL ephems.  I've included a term to convert Meeus'   */
      /* ecliptic-of-date result to J2000,  roughly. */

static void polar3_to_cartesian_r( double *vect,
                 const double lon, const double lat, const double r)
{
   *vect++ = cos( lat) * cos( lon) * r;
   *vect++ = cos( lat) * sin( lon) * r;
   *vect++ = sin( lat) * r;
}

const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078;

static void compute_rough_lunar_loc( const double t_cen, double *vect)
{
   const double precession_per_cy = (4700.29 * pi / 180.) / 3600.;
   const double lp0 = 218.3164591 * pi / 180.;
   const double lp1 = 481267.88134236 * pi / 180. - precession_per_cy;
   double lon = lp0 + lp1 * t_cen, lat = 0., r = 385000.56 / AU_IN_KM;
   const double *tptr;

   for( tptr = lunar_lon_r_terms; *tptr != 0.; tptr += 4)
      {
      const double angle = tptr[3] + tptr[2] * t_cen;

      lon += tptr[0] * sin( angle);
      r += tptr[1] * cos( angle);
      }

   for( tptr = lunar_lat_terms; *tptr; tptr += 3)
      lat += tptr[0] * sin( tptr[2] + tptr[1] * t_cen);
   polar3_to_cartesian_r( vect, lon, lat, r);
}

   /* Pluto formula from Meeus, _Astronomical Algorithms_, p 264 */

static void compute_rough_pluto_loc( const double t_cen, double *vect)
{
   const double p_arg = 238.96 * pi / 180. + (144.96 * pi / 180.) * t_cen;
   const double sin_p = sin( p_arg);
   const double cos_p = cos( p_arg);
   const double sin_2p = 2. * cos_p * sin_p;
   const double cos_2p = 2. * cos_p * cos_p - 1.;
   const double sin_3p = sin_p * cos_2p + sin_2p * cos_p;
   const double cos_3p = cos_p * cos_2p - sin_p * sin_2p;
   const double lon = p_arg - .001884 * pi / 180.
            - 19.799805 * pi / 180. * sin_p
            + 19.850055 * pi / 180. * cos_p
            +  0.897144 * pi / 180. * sin_2p
            -  4.954829 * pi / 180. * cos_2p
            +  0.611149 * pi / 180. * sin_3p
            +  1.211027 * pi / 180. * cos_3p;
   const double lat = -3.908239 * pi / 180.
            -  5.452852 * pi / 180. * sin_p
            - 14.974862 * pi / 180. * cos_p
            +  3.527812 * pi / 180. * sin_2p
            +  1.672790 * pi / 180. * cos_2p
            -  1.050748 * pi / 180. * sin_3p
            +  0.327647 * pi / 180. * cos_3p;
   const double r = 40.7241346
            +  6.6865439 * sin_p
            +  6.8951812 * cos_p
            -  1.1827535 * sin_2p
            -  0.0332538 * cos_2p
            +  0.1593179 * sin_3p
            -  0.1438890 * cos_3p;

   polar3_to_cartesian_r( vect, lon, lat, r);
}

int compute_rough_planet_loc( const double t_cen, const int planet_idx,
                                          double *vect)
{
   if( planet_idx == 10)
      compute_rough_lunar_loc( t_cen, vect);
   else if( planet_idx == 9)
      compute_rough_pluto_loc( t_cen, vect);
   else
      {
      int i;
      const SMALL_VSOP_TERM *vptr = vsop_terms[planet_idx];

      vect[0] = vptr->amplitude;
      vect[1] = vptr->phase;
      vect[2] = vptr->freq;
      vptr++;
      for( i = n_terms[planet_idx] - 1; i; i--, vptr++)
         vect[vptr->axis - 1] +=
               vptr->amplitude * cos( vptr->phase + vptr->freq * t_cen);
      }
   return( 0);
}

#define MERCURY_PERIHELION          0.30749951
#define MERCURY_APHELION            0.46669835
#define VENUS_PERIHELION            0.71843270
#define VENUS_APHELION              0.72823128
#define EARTH_PERIHELION            0.97
#define EARTH_APHELION              1.03
#define MARS_PERIHELION             1.38133346
#define MARS_APHELION               1.66599116
#define JUPITER_PERIHELION          4.95155843
#define JUPITER_APHELION            5.45516759
#define SATURN_PERIHELION           9.02063224
#define SATURN_APHELION            10.0535084
#define URANUS_PERIHELION          18.2860560
#define URANUS_APHELION            20.0964719
#define NEPTUNE_PERIHELION         29.8107953
#define NEPTUNE_APHELION           30.3271317
#define PLUTO_PERIHELION           29.6583407
#define PLUTO_APHELION             49.3050329

/* Perturbations from a given planet ought to be included if its rough
position is within the following ranges.  These include the fact that
the positions are somewhat approximate,  plus a little room beyond the
Hill sphere,  but they're a tad arbitrary/subject to tweaking.  Set 'em
too high,  and perturbations will be included when not needed.  Too low,
and perturbations will be omitted when they _are_ needed.

   The following ranges,  no longer used,  are essentially twice the Hill
sphere radii (except for Mercury,  Mars,  and Pluto,  for which the
uncertainty in the position formula is the limiting factor).  The new
ranges are probably overkill.  However,  it's not a very expensive check.
*/

#ifdef OBSOLETE_RANGES
#define MERCURY_RANGE            .02
#define VENUS_RANGE              .02
#define EARTH_RANGE              .02
#define MARS_RANGE               .02
#define JUPITER_RANGE            .7
#define SATURN_RANGE             .8
#define URANUS_RANGE             .9
#define NEPTUNE_RANGE           1.4
#define PLUTO_RANGE             1
#endif

#define MERCURY_RANGE            .1
#define VENUS_RANGE              .15
#define EARTH_RANGE              .2
#define MARS_RANGE               .1
#define JUPITER_RANGE            2
#define SATURN_RANGE             2
#define URANUS_RANGE             3
#define NEPTUNE_RANGE            3
#define PLUTO_RANGE             1

#define MERCURY_INNER_LIMIT   (MERCURY_PERIHELION - MERCURY_RANGE)
#define MERCURY_OUTER_LIMIT   (MERCURY_APHELION   + MERCURY_RANGE)
#define VENUS_INNER_LIMIT     (VENUS_PERIHELION   - VENUS_RANGE  )
#define VENUS_OUTER_LIMIT     (VENUS_APHELION     + VENUS_RANGE  )
#define EARTH_INNER_LIMIT     (EARTH_PERIHELION   - EARTH_RANGE  )
#define EARTH_OUTER_LIMIT     (EARTH_APHELION     + EARTH_RANGE  )
#define MARS_INNER_LIMIT      (MARS_PERIHELION    - MARS_RANGE   )
#define MARS_OUTER_LIMIT      (MARS_APHELION      + MARS_RANGE   )
#define JUPITER_INNER_LIMIT   (JUPITER_PERIHELION - JUPITER_RANGE)
#define JUPITER_OUTER_LIMIT   (JUPITER_APHELION   + JUPITER_RANGE)
#define SATURN_INNER_LIMIT    (SATURN_PERIHELION  - SATURN_RANGE )
#define SATURN_OUTER_LIMIT    (SATURN_APHELION    + SATURN_RANGE )
#define URANUS_INNER_LIMIT    (URANUS_PERIHELION  - URANUS_RANGE )
#define URANUS_OUTER_LIMIT    (URANUS_APHELION    + URANUS_RANGE )
#define NEPTUNE_INNER_LIMIT   (NEPTUNE_PERIHELION - NEPTUNE_RANGE)
#define NEPTUNE_OUTER_LIMIT   (NEPTUNE_APHELION   + NEPTUNE_RANGE)
#define PLUTO_INNER_LIMIT     (PLUTO_PERIHELION - PLUTO_RANGE)
#define PLUTO_OUTER_LIMIT     (PLUTO_APHELION   + PLUTO_RANGE)

/* The following code does a passable job of determining, automatically,
if the perturbations from a given planet are "strong enough" to merit
including in an initial orbit computation.  The idea is to avoid
computing perturbations from all objects at all times,  even when those
perturbations aren't all that strong.

   In determining whether any planet may be close enough to exert
serious perturbations,  we first check the square of its distance from
the sun;  if that's between the limits we're using for a given planet,
we then have to compute the location of that planet.  (We check the
square of the distance against the square of the planet's q-range and
Q+range, just to avoid taking a totally unnecessary square root.)

   Frequently,  that's enough to tell us that perturbations can be
ignored at present.  For example,  the outer limit for Mars is its
aphelion distance plus MARS_RANGE,  or 1.68599 AU.  The inner limit
for worrying about Jupiter is its perihelion distance minus JUPITER_RANGE,
or 4.25155 AU.  So for most of the Main Belt,  the following function
would just look at the distance of the planet from the sun,  and decide
to skip perturbations.

   Even if we're in the right heliocentric range,  we don't necessarily
compute the entire planetary position.  If, after computing the
approximate z-coordinate of the planet, we see that it can't be within
the limits set above,  then there's no point in going on to compute the
x and y coordinates.

   Incidentally,  z is checked first,  because that coordinate requires
the fewest terms.  (For the earth,  it requires none,  since the z-coord
in the ecliptic frame of the earth is basically zero at this rough level.)
If z is okay,  then y is checked,  then x.

   Pluto is a bit of an odd,  special case,  since its heliocentric
radius range overlaps that of Neptune.  Thus the extra code for it at
the end. */

int check_for_perturbers( const double t_cen, const double *vect)
{
   const double r2 = vect[0] * vect[0] + vect[1] * vect[1] + vect[2] * vect[2];
   int planet_to_check = 0, axis;

   if( r2 > JUPITER_INNER_LIMIT * JUPITER_INNER_LIMIT)
      {                 /* consider gas giants only */
      if( r2 < JUPITER_OUTER_LIMIT * JUPITER_OUTER_LIMIT)
         planet_to_check = 5;
      else if( r2 > SATURN_INNER_LIMIT * SATURN_INNER_LIMIT)
         {              /* Saturn,  Uranus,  or Neptune? */
         if( r2 < SATURN_OUTER_LIMIT * SATURN_OUTER_LIMIT)
            planet_to_check = 6;
         else if( r2 > URANUS_INNER_LIMIT * URANUS_INNER_LIMIT)
            {
            if( r2 < URANUS_OUTER_LIMIT * URANUS_OUTER_LIMIT)
               planet_to_check = 7;
            else if( r2 > NEPTUNE_INNER_LIMIT * NEPTUNE_INNER_LIMIT &&
                        r2 < NEPTUNE_OUTER_LIMIT * NEPTUNE_OUTER_LIMIT)
               planet_to_check = 8;
            }
         }
      }
   else        /* check inner planets */
      if( r2 < EARTH_OUTER_LIMIT * EARTH_OUTER_LIMIT)
         {
         if( r2 > EARTH_INNER_LIMIT * EARTH_INNER_LIMIT)
            planet_to_check = 3;
         else if( r2 > VENUS_INNER_LIMIT * VENUS_INNER_LIMIT)
            {
            if( r2 < VENUS_OUTER_LIMIT * VENUS_OUTER_LIMIT)
               planet_to_check = 2;
            }
         else if( r2 < MERCURY_OUTER_LIMIT * MERCURY_OUTER_LIMIT &&
                  r2 > MERCURY_INNER_LIMIT * MERCURY_INNER_LIMIT)
            planet_to_check = 1;
         }
      else if( r2 < MARS_OUTER_LIMIT * MARS_OUTER_LIMIT &&
               r2 > MARS_INNER_LIMIT * MARS_INNER_LIMIT)
         planet_to_check = 4;

   for( axis = 3; planet_to_check && axis > 0; axis--)
      {
      double delta = -vect[axis - 1];
      int i;
      const SMALL_VSOP_TERM *vptr = vsop_terms[planet_to_check];
      const double ranges[] = {0., MERCURY_RANGE, VENUS_RANGE, EARTH_RANGE,
                  MARS_RANGE, JUPITER_RANGE, SATURN_RANGE,
                  URANUS_RANGE, NEPTUNE_RANGE };

      if( axis == 3)
         delta += vptr->freq;
      else if( axis == 2)
         delta += vptr->phase;
      else if( axis == 1)
         delta += vptr->amplitude;
      vptr++;
      for( i = n_terms[planet_to_check] - 1; i; i--, vptr++)
         if( vptr->axis == axis)
            delta += vptr->amplitude * cos( vptr->phase + vptr->freq * t_cen);
                     /* Is the difference delta within limits? */
      if( delta > ranges[planet_to_check] || delta < -ranges[planet_to_check])
         planet_to_check = 0;
      }
   if( !planet_to_check && r2 > PLUTO_INNER_LIMIT * PLUTO_INNER_LIMIT
                        && r2 < PLUTO_OUTER_LIMIT * PLUTO_OUTER_LIMIT)
      {
      double pluto_loc[3];

      compute_rough_pluto_loc( (2448908.5 - 2451545.0) / 36525., pluto_loc);

      planet_to_check = 9;
      compute_rough_pluto_loc( t_cen, pluto_loc);
      for( axis = 0; planet_to_check && axis < 3; axis++)
         if( pluto_loc[axis] - vect[axis] > PLUTO_RANGE ||
                   pluto_loc[axis] - vect[axis] < -PLUTO_RANGE)
            planet_to_check = 0;
      }
   return( planet_to_check);
}

#ifdef TEST_MAIN

/* A little test routine I used to verify that the above functions were
computing reasonably correct planetary positions. */

#include <stdio.h>

void main( int argc, char **argv)
{
   const double jd = atof( argv[1]);
   const double j2000 = 2451545.0;
   const double t_cen = (jd - j2000) / 36525.;
   double vect[3];
   int i;

   for( i = 1; i < 9; i++)
      {
      compute_rough_planet_loc( t_cen, i, vect);
      printf( "%d   %8.4lf %8.4lf %8.4lf\n", i, vect[0], vect[1], vect[2]);
      }
   if( argc >= 5)
      {
      for( i = 0; i < 3; i++)
         vect[i] = atof( argv[i + 2]);
      printf( "Perturber = %d\n", check_for_perturbers( t_cen, vect));
      }
}
#endif
