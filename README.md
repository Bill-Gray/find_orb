# find_orb
Orbit determination from observations

   More about Find_Orb (including pre-built Windows executables) at

http://www.projectpluto.com/find_orb.htm

   This project includes code for building Linux,  Windows,  and BSD
versions of the interactive Find_Orb orbit determination software,
the non-interactive 'fo' software,  and the 'fo_serve.cgi' program
that underlies the on-line Find_Orb service.  Be warned,  though,
that (at least thus far) only the Linux versions of all these can
be built with what's currently posted.

   This project depends on three of my other projects :

   -- jpl_eph (code to read JPL ephemerides)
   -- sat_code (code for Earth-orbiting satellite ephemerides)
   -- lunar (basic astronomical ephemeris/time functions)

   I'd suggest getting all of these,  and running 'make' and
'sudo make install' for each.  That will put the relevant include
files and libraries where they're supposed to be for a Linux box.
(I'm still puzzling out where they're supposed to go on BSD or OS/X.)

   Once you have those three projects built and installed,  get this
project and run 'make',  and you should be off to the races (there is
no 'make install' for Find_Orb yet;  that's on my to-do list.)
