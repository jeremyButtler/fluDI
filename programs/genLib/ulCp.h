/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' ulCp SOF: Start Of File
'   - uses longs to copy contents of a string
'   - These functions are slower than strcpy and strlen,
'     but they allow deliminators to be used.
'   o header:
'     - Included libraries, definitions and preprocessor
'       checks
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
'   o fun10: swapDelim_ulCp
'     - swaps two strings until deliminator is found
'   o fun11: cmpDelim_ulCp
'     - compares two strings until deliminator is found
'   o fun12: eqlNull_ulCp
'     - compares two strings until null is found
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Included libraries, definitions and preprocessor
|     checks
\-------------------------------------------------------*/

#ifndef UNSINGED_LONG_COPY_H
#define UNSINGED_LONG_COPY_H

/*unsigned long shortcuts*/
/*the compiler should truncate these as needed*/
#define def_null_ulCp 0
#define def_newline_ulCp (unsigned long)0x0a0a0a0a0a0a0a0a
#define def_one_ulCp (unsigned long) 0x0101010101010101
#define def_tab_ulCp (unsigned long) 0x0909090909090909
#define def_highBit_ulCp (unsigned long)0x8080808080808080

/*for white space checks*/
#define def_31_ulCp (unsigned long) 0x0f0f0f0f0f0f0f0f
#define def_gt63_ulCp (unsigned long) 0xc0c0c0c0c0c0c0c0
   /*0xc0 is 12*/


#define def_bitsPerChar_ulCp 8
#define def_charInUL_ulCp sizeof(unsigned long)
#define def_modUL_ulCp (sizeof(unsigned long ) -1)
   /*The mod values are 7 (longs), 3 (ints), and 1
   `  (shorts), so, 1 - number of bytes will give me the
   `  mod value
   */

/*This should be compliled in.
`  This will only work for up to 64 bit numbers, after
`  that it will fail. However, it is very unlikely a
`  long will by 128 bits.
*/
#define def_shiftULBy_ulCp ((sizeof(unsigned long) >> 1) - (sizeof(unsigned long) >> 3))

/*Logic:
`   - sizeof(unsigned long) >> 1:
`     o This gets the number of bytes to shift to get
`       an 32 bit or 16 bit number to one bytes. Basically
`       it is divide by 2
`   - (sizeof(unsigned long) >> 3)
`     o This adds in 1 if I have eight bytes. It is needed
`       as the 64 bit correction step
`   - (sizeof(long) >> 1) - (sizeof(long) >> 3)
`     o Gives me the number of bits I need to shift to
`       get the number of characters in an long
`   - (sizeof(long) >> 3) << (sizeof(long) >> 4)
`     o This is not used, but would correct for a 128 bit
`       number. Each time I increase the bit size I need
`       to add a correction in
*/

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

/*-------------------------------------------------------\
| Fun10: swapDelim_ulCp
|   - swaps two strings until deliminator is found
| Input:
|   - firstStr:
|     o Pointer to string to first string to swap
|   - secStr:
|     o Pointer to second string to swap
|   - delimUL:
|     o delminator to end at (as long). Use makeULDelim
|       to build this deliminator
|   - delimSC:
|     o delminator (as char) to stop copying at
| Output:
|   - Modifies:
|     o firstStr to have secStr string
|     o secStr to have firstStr string
| Note:
|   - This will likely not be very good at swapping short
|     strings.
\-------------------------------------------------------*/
void
swapDelim_ulCp(
   signed char *firstStr,
   signed char *secStr,
   unsigned long delimUL,
   signed char delimSC
);

/*-------------------------------------------------------\
| Fun11: cmpDelim_ulCp
|   - compares two strings until deliminator is found
| Input:
|   - qryStr:
|     o Pointer to query string
|   - refStr:
|     o Pointer to reference strin
|   - delimUL:
|     o delminator to end at (as long). Use makeULDelim
|       to build this deliminator
|   - delimSC:
|     o delminator (as char) to stop copying at
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query > reference
|     o < 0 if query < reference
| Note:
|   - This will likely not be very good at comparing
|     short strings.
\-------------------------------------------------------*/
signed long
eql_ulCp(
   signed char *qryStr,
   signed char *refStr,
   unsigned long delimUL,
   signed char delimSC
);

/*-------------------------------------------------------\
| Fun12: eqlNull_ulCp
|   - compares two strings until null is found
| Input:
|   - qryStr:
|     o Pointer to query string
|   - refStr:
|     o Pointer to reference strin
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query > reference
|     o < 0 if query < reference
| Note:
|   - this will likely not be very good at comparing
|     short strings.
\-------------------------------------------------------*/
signed long
eqlNull_ulCp(
   signed char *qryStr,
   signed char *refStr
);

#endif

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
