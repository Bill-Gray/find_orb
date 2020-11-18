/* Copyright (C) 2018, Project Pluto

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.

   Microsoft Windows(R) lacks wget or cURL.  This is by no means
a full replacement for those programs.   However,  run

win_get (url) (filename)

   and the specified URL will be downloaded and saved under that
file name.  Add another command-line argument for 'verbose' output.
Compile with one of :

cl -W4 -Ox win_get.c urlmon.lib
x86_64-w64-mingw32-gcc -Wextra -Wall -O3 -pedantic -o win_get.exe win_get.c -lurlmon

   Let's start by suppressing MSVC's beloved nuisance warnings. */

#define _CRT_SECURE_NO_WARNINGS

#include "windows.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

int main( const int argc, const char **argv)
{
   const char *output_name = (argc < 2 ? "zzz1" : argv[2]);
   HRESULT rval = URLDownloadToFile( NULL, argv[1], output_name, 0, NULL);

   if( argc > 3)     /* extra parameter interpreted to mean 'verbose' */
      {
      FILE *ifile;

      printf( "Downloading '%s'\n", argv[1]);
      printf( "URLD %d\n", (int)rval);
      ifile = fopen( output_name, "rb");
      assert( ifile);
      fseek( ifile, 0L, SEEK_END);
      printf( "Size %ld\n", (long)ftell( ifile));
      fclose( ifile);
      }
   return( rval != S_OK);
}
