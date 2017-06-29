/* clipfunc.cpp: functions for getting/saving data to/from the Windows clipboard

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
#include <windows.h>

int copy_buffer_to_clipboard( const char *contents, const long length);
int copy_file_to_clipboard( const char *filename);    /* clipfunc.cpp */
int clipboard_to_file( const char *filename, const int append); /* clipfunc.cpp */

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
   GlobalFree(ptr1);

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

int clipboard_to_file( const char *filename, const int append)
{
    int rval;

    if( !OpenClipboard(NULL))
        rval =  -1;
    else
        {
        char *text;
        FILE *ofile;
        const char *permits = (append ? "ab" : "wb");

        if( (text = (char *)GetClipboardData( CF_OEMTEXT)) == NULL)
            rval = -2;
        else if( (ofile = fopen( filename, permits)) != NULL)
            {
            rval = 0;
            fwrite( text, 1, strlen( text), ofile);
            fclose( ofile);
            }
        else
            rval = -3;
        CloseClipboard();
        }
    return rval;
}
