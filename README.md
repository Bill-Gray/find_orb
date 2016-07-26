# find_orb

### Orbit determination from observations

More about Find_Orb (what it does,  how to use it,  pre-built Windows
executables) at

http://www.projectpluto.com/find_orb.htm

This project includes code for building Linux,  Windows,  and BSD
(and possibly OS/X) versions of the interactive Find_Orb orbit
determination software, the non-interactive `fo` software,  and
the `fo_serve.cgi` program that underlies the
[on-line Find_Orb service](http://www.projectpluto.com/fo.htm).  Be
warned,  though, that (at least thus far) only the Linux and
BSD versions (OS/X is possible,  but untested) can be built with
what's currently posted.

This project depends on three of my other projects :

- [`jpl_eph`](https://github.com/Bill-Gray/jpl_eph) (code to read JPL ephemerides)
- [`sat_code`](https://github.com/Bill-Gray/sat_code) (code for Earth-orbiting satellite ephemerides)
- [`lunar`](https://github.com/Bill-Gray/lunar) (basic astronomical ephemeris/time functions)

I'd suggest getting all three of these,  either by cloning them or
just downloading the ZIPballs,  and running `make` and `sudo make
install` for each.  That will copy the relevant include files and
libraries to `/usr/local/include` and `/usr/local/lib`.  This appears
to work for Linux and BSD,  and may work with OS/X as well (I don't
have Apple products to try it out).  For BSD (and possibly OS/X),
use `gmake` instead of `make`.

Obviously, I'll have to change the procedure for Windows.

Once you have those three projects built and installed,  get this
project and run `make`,  and Find_Orb should be built,  as well as the
aforementioned `fo` and `fo_serve.cgi`.  Note that at present,  there
is no `make install` for Find_Orb yet;  that's on my to-do list.
