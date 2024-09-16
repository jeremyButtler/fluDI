/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' ulCp SOF: Start Of File
'   - uses longs to copy contents of a string
'   - These functions are slower than strcpy and strlen,
'     but they allow deliminators to be used.
'   o header:
'     - included libraries
'   o fun01: cpLen_ulCp
'     - copies cpStr into dupStr using unsigned longs
'   o fun02: mkDelim_ulCp
'     - makes an unsigned long delimintor from a character
'       deliminator for use in cpStrDelim
'   o .c note01: ifDelim_ulCp
'     - logic for detecting if the input deliminator is in
'       input unsigned long
'   o fun03: cpDelim_ulCp
'     - copies string until deliminator is found
'   o fun04: lenStr_ulCp
'     - finds the length of a string using unsigned longs
'   o fun05: lenStrNull_ulCp
'     - finds the length of a string using unsigned longs
'       but also stops at null '\0'
'   o fun06: endLine_ulCp
'     - finds the end of a c-string. This assumes that
'       the line ends in an '\0' or an '\n'
'   o fun07: cpLine_ulCp
'     - copies string until end of line
'   o fun08: cpWhite_ulCp
'     - copies string until white space
'   o fun09: rmWhite_ulCp
'     - removes white space from c-string
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "ulCp.h"

/*-------------------------------------------------------\
| Fun01: cpLen_ulCp
|   - copies cpStr into dupStr using unsigned longs
| Input:
|   - dupStr:
|     o pointer to string to copy cpStr into
|   - cpStr:
|     o pointer to string to copy
|   - lenUI:
|     o length of cpStr
| Output:
|   - Modifies:
|     o  dupStr to hold lenUI characters from cpStr
\-------------------------------------------------------*/
void
cpLen_ulCp(
   signed char *dupStr,
   signed char *cpStr,
   unsigned int lenUI
){
   unsigned long *cpUL = (unsigned long *) cpStr;
   unsigned long *dupUL = (unsigned long *) dupStr;
   unsigned int uiChar = 0;

   for(
      uiChar= 0;
      uiChar < (lenUI >> def_shiftULBy_ulCp);
      ++uiChar
   ) dupUL[uiChar] = cpUL[uiChar];

   for(
      uiChar = ( lenUI - (lenUI & def_modUL_ulCp) );
      uiChar < lenUI;
     ++uiChar
   ) dupStr[uiChar] = cpStr[uiChar];

   dupStr[lenUI] = '\0';
} /*cpLen_ulCp*/

/*-------------------------------------------------------\
| Fun02: mkDelim_ulCp
|   - makes an unsigned long delimintor from a character
|     deliminator for use in cpStrDelim
| Input:
|   - delimC:
|     o deliminator to use to end a string copy
| Output:
|   - Modifies:
|     o dupStr to hold the characters from cpStr
\-------------------------------------------------------*/
unsigned long
mkDelim_ulCp(
   signed char delimSC
){
   unsigned long uiRep = 0;
   unsigned long delimUL = (unsigned char) delimSC;

   for(
      uiRep = sizeof(unsigned long) << 2;
      uiRep >= def_bitsPerChar_ulCp;
      uiRep >>= 1
   ) delimUL |= (delimUL << uiRep);

   return delimUL;
} /*mkDelim_ulCp*/
   /*Logic:
   `   - macDelimUL  = (unsigned char) delimC:
   `     o the (unsigned char) is needed to avoid negative
   `       values being converted to their long
   `       equivulents
   `   - sizeof(unsigned long) << 2:
   `     o Gets me to the half way point in bits. So
   `       8 (64 bit) becomes 32, 4 (32 bit) becomes 16,
   `       an 2 (16 bit) becomes 8
   `   - iRepMac << 1:
   `     o Divides the shift range by 2 each round, so
   `       that each round fills in a different part of
   `       the long
   `   - macDelimUL |= (macDelimUL << iRepMac);
   `     o Moves the deliminator up (shifts) to an
   `       unfilled position.
   */

/*-------------------------------------------------------\
| Note01: ifDelim_ulCp
|   - logic for detecting if the input deliminator is in
|     input unsigned long
\-------------------------------------------------------*/
/*Logic:
`  * checkUL holds the comparision/is temporary variable
`  * strUL are the bytes being compared
`  * delimUL is the deliminator to check for (as ulong)
`  - checkUL = *strUL ^ delimUL:
`    o Converts all values matching delimUL to 0, and
`      everything else to > 0. Saves result to retUL
`  - checkUL = (checkUL & def_one_ulCp) + checkUL
`    o convets odd numbers to even numbers. this
`      means there is no 1;
`    o this is needed to ensure -1 will not cause an
`      overflow when 1 is next to an deliminator. The
`      & 1 ensures that -1 will alwas remove the sign
`      bit, if one character was one off.
`      - only use if you plan to find the deliminator
`        position using the long. Else, remove
`  - checkUL = ((checkUL & 1) + checkUL)
`    o subtracts one from every character. This results
`      in an overflow for only 0, while every other
`      value is >= 1, but can go to 0 if next number
`      overflowed. This sets the value 0 numbers to -1
`  - checkUL = checkUL & def_highBit_ulCp
`    o This keeps only the overflow bit (sign bit). This
`      means that every non-deliminator value becomes
`      0.
*/

/*-------------------------------------------------------\
| Fun03: cpDelim_ulCp
|   - copies string until deliminator is found
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
|   - delimUL:
|     o delminator to end at (as long). Use makeULDelim
|       to build this deliminator
|   - delimSC:
|     o delminator (as char) to stop copying at
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
unsigned int
cpDelim_ulCp(
   signed char *dupStr,
   signed char *cpStr,
   unsigned long delimUL,
   signed char delimSC
){
   unsigned long *cpUL = (unsigned long *) (cpStr);
   unsigned long *dupUL = (unsigned long *) (dupStr);

   signed char *dupTmpStr = 0;
   signed char *cpTmpStr = 0;

   unsigned int uiChar = 0;
   unsigned long checkUL = 0;

   /*see note01 ifDelim_ulCp for logic)*/
   checkUL = *cpUL ^ delimUL;
   checkUL -= def_one_ulCp;
   checkUL = checkUL & def_highBit_ulCp;

   while(! checkUL)
   { /*Loop: Copy cpStr to dupStr*/
      *dupUL++ = *cpUL++;

      checkUL = *cpUL ^ delimUL;
      checkUL -= def_one_ulCp;
      checkUL = checkUL & def_highBit_ulCp;
   } /*Loop: Copy cpStr to dupStr*/

   cpTmpStr = (signed char *) cpUL;
   dupTmpStr = (signed char *) dupUL;

   while(*cpTmpStr != delimSC)
      *dupTmpStr++ = *cpTmpStr++;

   uiChar = dupTmpStr - dupStr;
   dupStr[uiChar] = '\0';
   return uiChar; /*the number of characters copied*/
} /*cpDelim_ulCp*/

/*-------------------------------------------------------\
| Fun04: lenStr_ulCp
|   - finds the length of a string using unsigned longs
| Input:
|   - inStr:
|     o c-string or string with deliminator to find length
|       of
|   - delimUL:
|     o deliminator (as unsigned long (fun02)) at end of
|       string
|   - delimSC:
|     o deliminator (as char) at end of string
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
lenStr_ulCp(
   signed char *inStr,
   unsigned long delimUL,
   signed char delimSC
){
   unsigned long *ptrUL = (unsigned long *) inStr;
   unsigned long checkUL = 0;
   unsigned int uiLenStr = 0;

   checkUL = *ptrUL ^ delimUL;
   checkUL -= def_one_ulCp;
   checkUL = checkUL & def_highBit_ulCp;

   while(! checkUL)
   { /*Loop: find length*/
      uiLenStr += def_charInUL_ulCp;
      ++ptrUL;

      checkUL = *ptrUL ^ delimUL;
      checkUL -= def_one_ulCp;
      checkUL = checkUL & def_highBit_ulCp;
   } /*Loop: find length*/

   while(inStr[uiLenStr] != delimSC)
      ++uiLenStr;

   return uiLenStr;
} /*lenStr_ulCp*/

/*-------------------------------------------------------\
| Fun05: lenStrNull_ulCp
|   - finds the length of a string using unsigned longs
|     but also stops at null '\0'
| Input:
|   - inStr:
|     o c-string or string with deliminator to find length
|       of
|   - delimUL:
|     o deliminator (as unsigned long (fun02)) at end of
|       string
|   - delimSC:
|     o deliminator (as char) at end of string
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
lenStrNull_ulCp(
   signed char *inStr,
   unsigned long delimUL,
   signed char delimSC
){
   unsigned long *ptrUL = (unsigned long *) inStr;
   unsigned long checkUL = 0;
   unsigned int uiLenStr = 0;

   checkUL = *ptrUL | delimUL; /*so null goes to 0*/
   checkUL ^= delimUL;
   checkUL -= def_one_ulCp;
   checkUL = checkUL & def_highBit_ulCp;

   while(! checkUL)
   { /*Loop: find length*/
      uiLenStr += def_charInUL_ulCp;
      ++ptrUL;

      checkUL = *ptrUL | delimUL; /*so null goes to 0*/
      checkUL ^= delimUL;
      checkUL -= def_one_ulCp;
      checkUL = checkUL & def_highBit_ulCp;
   } /*Loop: find length*/

   while(inStr[uiLenStr] != delimSC)
      ++uiLenStr;

   return uiLenStr;
} /*lenStrNull_ulCp*/

/*-------------------------------------------------------\
| Fun06: endLine_ulCp
|   - finds the end of a c-string. This assumes that the
|     line ends in an '\0' or an '\n'
| Input:
|   - inStr:
|     o c-string or string to look for end in
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
endLine_ulCp(
   signed char *inStr
){
   unsigned long *ptrUL = (unsigned long *) inStr;
   unsigned int uiLenStr = 0;

   while( !
     (
          (
               ((*ptrUL) & (~ def_newline_ulCp))
             - def_one_ulCp
          ) & def_highBit_ulCp
     )
   ){ /*Loop: find end of line or null*/
      uiLenStr += def_charInUL_ulCp;
      ++ptrUL;
   } /*Loop: find end of line or null*/

   while(inStr[uiLenStr] & (~ '\n'))
      ++uiLenStr;

   return uiLenStr;
} /*endLine_ulCp*/

/*-------------------------------------------------------\
| Fun07: cpLine_ulCp
|   - copies string until end of line
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
unsigned int
cpLine_ulCp(
   signed char *dupStr,
   signed char *cpStr
){
   unsigned long *cpUL = (unsigned long *) (cpStr);
   unsigned long *dupUL = (unsigned long *) (dupStr);

   signed char *dupTmpStr = 0;
   signed char *cpTmpStr = 0;

   unsigned int uiChar = 0;
   unsigned long checkUL = 0;

   /*see note01 ifDelim_ulCp for logic)*/
   checkUL = *cpUL & (~ def_newline_ulCp);
   checkUL -= def_one_ulCp;
   checkUL &= def_highBit_ulCp;

   while(! checkUL)
   { /*Loop: copy line with longs*/
      *dupUL++ = *cpUL++;

      checkUL = *cpUL & (~ def_newline_ulCp);
      checkUL -= def_one_ulCp;
      checkUL &= def_highBit_ulCp;
   } /*Loop: copy line with longs*/

   cpTmpStr = (signed char *) cpUL;
   dupTmpStr = (signed char *) dupUL;

   while(*cpTmpStr & (~ '\n'))
      *dupTmpStr++ = *cpTmpStr++;

   uiChar = dupTmpStr - dupStr;
   dupStr[uiChar] = '\0';
   return uiChar; /*the number of characters copied*/
} /*cpLine_ulCp*/

/*-------------------------------------------------------\
| Fun08: cpWhite_ulCp
|   - copies string until white space
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
unsigned int
cpWhite_ulCp(
   signed char *dupStr,
   signed char *cpStr
){
   unsigned long *cpUL = (unsigned long *) (cpStr);
   unsigned long *dupUL = (unsigned long *) (dupStr);

   signed char *dupTmpStr = 0;
   signed char *cpTmpStr = 0;

   unsigned int uiChar = 0;
   unsigned long checkUL = 0;

   /*see note01 ifDelim_cpUL for logic)*/
   checkUL = *cpUL + def_31_ulCp;
   checkUL &= def_gt63_ulCp;
   checkUL |= (checkUL << 1);
   checkUL &= def_highBit_ulCp;
   checkUL ^= def_highBit_ulCp;

   /* Logic:
   `   - checkUL = cpUL + def_31_ulCp
   `     o converts 33 or above to 64 or >
   `       - shifts non-space one bit
   `     o 64 goes to 128 or greater (8th bit)
   `   - checkUL &= def_gt63_ulCp:
   `     o everything beneath 63 goes to 0
   `   - checkUL |= (checkUL << 1):
   `     o moves 7th bit (64) to 8th bit (128)
   `   - checkUL &= def_highBit_ulCp:
   `     o clears all bits but 8th bith (128)
   `     o at this point have only 0x80 or 0x00
   `   - checkUL ^= def_highBit_ulCp:
   `     o sets everything to zero
   */

   while(! checkUL)
   { /*Loop: copy line with longs*/
      *dupUL++ = *cpUL++;

      checkUL = *cpUL + def_31_ulCp;
      checkUL &= def_gt63_ulCp;
      checkUL |= (checkUL << 1);
      checkUL &= def_highBit_ulCp;
      checkUL ^= def_highBit_ulCp;
   } /*Loop: copy line with longs*/

   cpTmpStr = (signed char *) cpUL;
   dupTmpStr = (signed char *) dupUL;

   while(*cpTmpStr > 33)
      *dupTmpStr++ = *cpTmpStr++;

   uiChar = dupTmpStr - dupStr;
   dupStr[uiChar] = '\0';
   return uiChar; /*the number of characters copied*/
} /*cpWhite_ulCp*/


/*-------------------------------------------------------\
| Fun09: rmWhite_ulCp
|   - removes white space from c-string
| Input:
|   - inStr:
|     o c-string to remove whites space
| Output:
|   - Modifies:
|     o  inStr to have white space removed
|   - Returns:
|     o number of characters in inStr
\-------------------------------------------------------*/
unsigned int
rmWhite_ulCp(
   signed char *inStr
){
   unsigned long *cpUL = (unsigned long *) inStr;
   unsigned long *dupUL = (unsigned long *) inStr;

   signed char *dupStr = inStr;
   signed char *cpStr = inStr;

   unsigned int uiChar = 0;
   unsigned long checkUL = 0;

   while(*cpStr != '\0')
   { /*Loop: remove white space*/
      checkUL = *cpUL + def_31_ulCp;
      checkUL &= def_gt63_ulCp;
      checkUL |= (checkUL << 1);
      checkUL &= def_highBit_ulCp;
      checkUL ^= def_highBit_ulCp;
      /*see fun08 for logic*/

      if(! checkUL)
         *dupUL++ = *cpUL++;
      else
      { /*Else: have white space, manually copy*/
         cpStr = (signed char *) cpUL;
         dupStr = (signed char *) dupUL;
         
         for(
            uiChar = 0;
            uiChar < sizeof(unsigned long);
            ++uiChar
         ){ /*Loop: remove white space*/
            if(*cpStr == '\0')
               break; /*finished*/
            else if(*cpStr < 32)
               ++cpStr;
            else
               *cpStr++ = *dupStr++;
         } /*Loop: remove white space*/

         cpUL = (unsigned long *) cpStr;
         dupUL = (unsigned long *) dupStr;
      } /*Else: have white space, manually copy*/
   } /*Loop: remove white space*/

   *dupStr = '\0';
   return dupStr - inStr;
} /*rmWhite_ulCp*/

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
