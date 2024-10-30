/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' charCp SOF: Start Of File
'   - Copies contents of string 1 to string using chars
'   - These functions are slower than strcpy and strlen,
'     but they allow deliminators to be used.
'   o header:
'     - included libraries
'   o fun01: cpLen_charCp
'    - copies cpStr into dupStr using characters
'   o fun02: cpDelim_charCp
'     - copies cpStr into dupStr until delimC is found
'   o fun03: lenStr_charCp
'     - finds the length of a string using characters
'   o fun04: lenStrNull_charCp
'     - finds the length of a string using characters, but
'       also stops at null '\0'
'   o fun05: eql_charCp
'     - checks to see if two strings are equal
'   o fun06: match_charCp
'     - checks to see if two strings are equal, but does
'       not check to see if there is anything past the
'       query's deliminator. This is to deal with strings
'       with differnt deliminators.
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "charCp.h"

/*-------------------------------------------------------\
| Fun01: cpLen_charCp
|   - copies cpStr into dupStr using characters
| Input:
|   - dupStr:
|     o pointer to string to copy cpStr into
|   - cpStr:
|     o pointer to string to copy
|   - lenUI:
|     o length of cpStr
| Output:
|   - Modifies:
|     o dupStr to hold lenUI characters from cpStr
\-------------------------------------------------------*/
void
cpLen_charCp(
   signed char *dupStr,
   signed char *cpStr,
   unsigned int lenUI
){
   dupStr[lenUI] = '\0';

   while(lenUI > 0)
   { /*Loop: copy strings*/
      --lenUI;
      dupStr[lenUI] = cpStr[lenUI];
   } /*Loop: copy strings*/
} /*cpLen_charCp*/

/*-------------------------------------------------------\
| Fun02: cpDelim_charCp
|   - copies cpStr into dupStr until delimC is found
| Input:
|   - dupStr:
|     o pointer to string to copy cpStr into
|   - cpStr:
|     o pointer to string to copy
|   - delimUL:
|     o deliminatro to end at (as long). Use makeULDelim
|       to build this deliminator
| Output:
|   - Modifies:
|     o dupStr to hold the characters from cpStr
|   - Returns:
|     o length of copied string
\-------------------------------------------------------*/
unsigned int
cpDelim_charCp(
   signed char *dupStr,
   signed char *cpStr,
   signed char delimSC
){
   signed char *dupTmpStr = dupStr;

   while(*cpStr != delimSC)
      *dupTmpStr++ = *cpStr++;

   *dupTmpStr = '\0';
   return (unsigned int) (dupTmpStr - dupStr);
} /*cpDelim_charCp*/

/*-------------------------------------------------------\
| Fun03: lenStr_charCp
|   - finds the length of a string using characters
| Input:
|   - inStr:
|     o c-string or string with deliminator to find length
|       of
|   - delimSC:
|     o deliminator marking the end of the string
| Output:
|   - returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
lenStr_charCp(
   signed char *inStr,
   signed char delimSC
){
   signed char *tmpStr = inStr;
   while(*tmpStr++ != delimSC) ;
   return (unsigned int) (tmpStr - inStr - 1);
} /*lenStr_charCp*/

/*-------------------------------------------------------\
| Fun04: lenStrNull_charCp
|   - finds the length of a string using characters, but
|     also stops at null '\0'
| Input:
|   - inStr:
|     o c-string or string with deliminator to find length
|       of
|   - delimSC:
|     o deliminator marking the end of the string
| Output:
|   - returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
lenStrNull_charCp(
   signed char *inStr,
   signed char delimSC
){
   signed char *tmpStr = inStr;
   while(*tmpStr++ != delimSC && *tmpStr != '\0') ;
   return (unsigned int) (tmpStr - inStr - 1);
} /*lenStrNull_charCp*/

/*-------------------------------------------------------\
| Fun05: eql_charCp
|   - checks to see if two strings are equal
| Input:
|   - qryStr:
|     o the query c-string to compare against the
|       reference
|   - refStr:
|     o the reference c-string to compare against the
|       query
| Output:
|   - Returns:
|     o < 0 for qry less than reference
|     o 0 for qry is same as the reference
|     o >0 for qry is greater than the reference
\-------------------------------------------------------*/
signed int
eql_charCp(
   signed char *qryStr,
   signed char *refStr,
   signed char delimSC
){
   signed int siChar = 0;

   while(qryStr[siChar] == refStr[siChar])
   { /*Looop*/
      if(qryStr[siChar] == delimSC)
         break;

      ++siChar;
   } /*Looop*/

   return (signed int) (qryStr[siChar] - refStr[siChar]);
} /*eql_charCp*/

/*-------------------------------------------------------\
| Fun06: match_charCp
|   - checks to see if two strings are equal, but does
|     not check to see if there is anything past the
|     query's deliminator. This is to deal with strings
|     with differnt deliminators.
| Input:
|   - qryStr:
|     o the query c-string to compare against the
|       reference
|   - refStr:
|     o the reference c-string to compare against the
|       query
|   - delimSC:
|     o deliminator to stop at
|   - endI:
|     o Will have the position before the deliminator.
|       This assumes that you know what you are doing,
|       so if you input endI > 0; this will assume you
|       do not want to compare at 0. This is handy for
|       looping through headers in tsv files.
| Output:
|   - Returns:
|     o < 0 for qry less than reference
|     o 0 for qry is same as the reference
|     o >0 for qry is greater than the reference
\-------------------------------------------------------*/
signed int
match_charCp(
   signed char *qryStr,
   signed char *refStr,
   signed char delimSC,
   signed int *endSI
){
   signed char *qryTmpStr = qryStr;

   while(*qryTmpStr++ == *refStr++)
   { /*Loop: compare query and reference*/
      if(*qryTmpStr == delimSC)
      { /*If: I found the  deliminator*/
         --qryTmpStr;
         --refStr;
         break;
      } /*If: I found the  deliminator*/
   } /*Loop: compare query and reference*/

   *endSI += qryTmpStr - qryStr;
   return *qryTmpStr - *refStr;
} /*match_charCp*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconvient / not possible, this code is under the
:   MIT license.
: 
: Public domain:
: 
: This is free and unencumbered software released into the
:   public domain.
: 
: Anyone is free to copy, modify, publish, use, compile,
:   sell, or distribute this software, either in source
:   code form or as a compiled binary, for any purpose,
:   commercial or non-commercial, and by any means.
: 
: In jurisdictions that recognize copyright laws, the
:   author or authors of this software dedicate any and
:   all copyright interest in the software to the public
:   domain. We make this dedication for the benefit of the
:   public at large and to the detriment of our heirs and
:   successors. We intend this dedication to be an overt
:   act of relinquishment in perpetuity of all present and
:   future rights to this software under copyright law.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO
:   EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM,
:   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
:   CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
:   IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
:   DEALINGS IN THE SOFTWARE.
: 
: For more information, please refer to
:   <https://unlicense.org>
: 
: MIT License:
: 
: Copyright (c) 2024 jeremyButtler
: 
: Permission is hereby granted, free of charge, to any
:   person obtaining a copy of this software and
:   associated documentation files (the "Software"), to
:   deal in the Software without restriction, including
:   without limitation the rights to use, copy, modify,
:   merge, publish, distribute, sublicense, and/or sell
:   copies of the Software, and to permit persons to whom
:   the Software is furnished to do so, subject to the
:   following conditions:
: 
: The above copyright notice and this permission notice
:   shall be included in all copies or substantial
:   portions of the Software.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
:   AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
:   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
