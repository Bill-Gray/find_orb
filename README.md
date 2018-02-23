# find_orb

### Orbit determination from observations

More about Find_Orb (what it does,  how to use it,  pre-built Windows
executables) at

[`https://www.projectpluto.com/find_orb.htm`](https://www.projectpluto.com/find_orb.htm)

This project includes code for building Linux,  Windows,  and BSD
(and possibly OS/X) versions of the interactive Find_Orb orbit
determination software, the non-interactive `fo` software,  and
the `fo_serve.cgi` program that underlies the
[on-line Find_Orb service](https://www.projectpluto.com/fo.htm).
(For a lot of purposes,  the on-line service or pre-built .EXEs
may be all you really need.)

Right now,  only the Linux and BSD versions can be built with
what's posted here.  (You can probably build for OS/X,  too,
but I've not heard any reports on that recently.) I've not
gotten around to documenting the Windows build process yet;  I
just provide the aforementioned pre-built EXEs.  Sadly,  the
Windows builds require the Microsoft compiler,  due to an
unfortunate early choice to use MFC.

You can [click here for directions on building Find_Orb from
the source in this repository.](https://projectpluto.com/find_sou.htm).

As is described at the above link,  this project depends on three
of my other projects :

- [`jpl_eph`](https://github.com/Bill-Gray/jpl_eph) (code to read JPL ephemerides)
- [`sat_code`](https://github.com/Bill-Gray/sat_code) (code for Earth-orbiting satellite ephemerides)
- [`lunar`](https://github.com/Bill-Gray/lunar) (basic astronomical ephemeris/time functions)

At some point,  I'll probably document the build procedure for
Windows,  but it does appear that most Windows users just want
pre-built executables.
