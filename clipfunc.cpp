/* clipfunc.cpp: functions for getting/saving data to/from the clipboard

Copyright (C) 2012, Project Pluto

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

#include <stdio.h>
#include "watdefs.h"

int copy_buffer_to_clipboard( const char *contents, const long length);
int copy_file_to_clipboard( const char *filename);    /* clipfunc.cpp */
int clipboard_to_file( const char *filename, const int append,
                           const bool use_selection); /* clipfunc.cpp */

#ifdef _WIN32
#include <windows.h>

      /* Following is copied shamelessly from pdcclip.c from PDCurses. */
      /* With one strange addition:  "solo" line feeds are expanded to */
      /* CR/LFs,  resulting in the output being 'n_added' bytes larger. */

int copy_buffer_to_clipboard( const char *contents, const long length)
{
   HGLOBAL ptr1;
   LPTSTR ptr2;
   int i, n_added;

   if (!OpenClipboard(NULL))
      return -3;

   for( i = n_added = 0; i < length; i++)
      if( contents[i] != 13 && contents[i + 1] == 10)
         n_added++;
   ptr1 = GlobalAlloc(GMEM_MOVEABLE|GMEM_DDESHARE,
      (length + 1 + n_added) * sizeof(TCHAR));

   if (!ptr1)
      return -4;

   ptr2 = (LPTSTR)GlobalLock(ptr1);
   if( !ptr2)
   {
      GlobalFree(ptr1);
      return -6;
   }

   for( i = n_added = 0; i < length; i++)
      {
      ptr2[i + n_added] = contents[i];
      if( contents[i] != 13 && contents[i + 1] == 10)
         {
         n_added++;
         ptr2[i + n_added] = 13;
         }
      }
   ptr2[length + n_added] = 0;
// memcpy((char *)ptr2, contents, length + 1);
   GlobalUnlock(ptr1);
   EmptyClipboard();

   if (!SetClipboardData( CF_TEXT, ptr1))
   {
      GlobalFree(ptr1);
      return -5;
   }

   CloseClipboard();

   return 0;
}

int copy_file_to_clipboard( const char *filename)
{
   FILE *ifile = fopen( filename, "rb");
   int rval;

   if( ifile)
      {
      size_t filesize;
      char *buff;

      fseek( ifile, 0L, SEEK_END);
      filesize = ftell( ifile);
      fseek( ifile, 0L, SEEK_SET);
      buff = (char *)malloc( filesize);
      if( buff)
         {
         fread( buff, 1, filesize, ifile);
         rval = copy_buffer_to_clipboard( buff, (long)filesize);
         free( buff);
         }
      else
         rval = -2;
      fclose( ifile);
      }
   else
      rval = -1;
   return( rval);
}

int clipboard_to_file( const char *filename, const int append,
                                          const bool use_selection)
{
    int rval;

            /* sadly,  we can't get at the selected text in MSWin */
    INTENTIONALLY_UNUSED_PARAMETER( use_selection);
    if( !OpenClipboard(NULL))
        rval =  -1;
    else
        {
        void *text;
        FILE *ofile;
        const char *permits = (append ? "ab" : "wb");
        HANDLE handle = GetClipboardData( CF_OEMTEXT);

        if( !handle)
            rval = -2;
        else if( (text = GlobalLock( handle)) == NULL)
            rval = -4;
        else if( (ofile = fopen( filename, permits)) != NULL)
            {
            rval = 0;
            fwrite( text, 1, strlen( (char *)text), ofile);
            fclose( ofile);
            }
        else
            rval = -3;
        CloseClipboard();
        }
    return rval;
}
#elif defined __PDCURSES__ && !defined VT

#include <stdlib.h>
#include "curses.h"

int clipboard_to_file( const char *filename, const int append,
                                          const bool use_selection)
{
   long size = -99;
   char *contents;
   int err_code;

   INTENTIONALLY_UNUSED_PARAMETER( use_selection);
   err_code = PDC_getclipboard( &contents, &size);
   if( err_code == PDC_CLIP_SUCCESS)
      {
      FILE *ofile = fopen( filename, "wb");

      if( ofile)
         {
         fwrite( contents, size, 1, ofile);
         fclose( ofile);
         }
      PDC_freeclipboard( contents);
      }
   return( err_code);
}

int copy_file_to_clipboard( const char *filename)
{
   FILE *ifile = fopen( filename, "rb");
   int err_code = -1;

   if( ifile)
      {
      size_t length, bytes_read;
      char *buff;

      fseek( ifile, 0L, SEEK_END);
      length = (size_t)ftell( ifile);
      fseek( ifile, 0L, SEEK_SET);
      buff = (char *)malloc( length + 1);
      assert( buff);
      buff[length] = '\0';
      bytes_read = fread( buff, 1, length, ifile);
      assert( bytes_read == length);
      fclose( ifile);
      err_code = PDC_setclipboard( buff, length);
      free( buff);
      }
   return( err_code);
}
#else    /* non-PDCurses, non-Windows:  use xclip */
#include <stdlib.h>

int clipboard_to_file( const char *filename, const int append,
                                          const bool use_selection)
{
   char cmd[80];

   snprintf( cmd, sizeof( cmd), "xclip -o %s >%c %s",
            (use_selection ? "" : "-se c"),
            (append ? '>' : ' '), filename);
   return( system( cmd));
}

int copy_file_to_clipboard( const char *filename)
{
   char cmd[80];

   snprintf( cmd, sizeof( cmd), "xclip -i %s", filename);
   return( system( cmd));
}
#endif
