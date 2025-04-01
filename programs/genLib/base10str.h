/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' base10str SOF: Start Of Header
'   - functions to converts c-strings to numbers
'   o header:
'     - has includes and the defined size of each number
'   o fun01: strToUL_base10str
'     - converts a base 10 c-string to an unsigned long
'   o fun02: strToUI_base10str
'     - converts a base 10 c-string to an unsigned int
'   o fun03: strToUS_base10str
'     - converts a base 10 c-string to an unsigned short
'   o fun04: strToUC_base10str
'     - converts a base 10 c-string to an unsigned char
'   o fun05: strToSL_base10str
'     - converts base 10 c-string to a signed long
'   o fun06: strToSI_base10str
'     - converts base 10 c-string to a signed int
'   o fun07: strToSS_base10str
'     - converts base 10 c-string to a signed short
'   o fun08: strToSC_base10str
'     - converts base 10 c-string to a signed char
'   o fun09: strToF_base10str
'     - converts base 10 c-string to float
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - Has includes and the defined size of each number
|  o header sec-01:
|    - guards
|  o header sec-02:
|    - maximum digits and values for unsigned data types
|  o header sec-03:
|    - maximum digits and values for signed data types
\-------------------------------------------------------*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec-01:
^  - guards
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifndef BASE_10_STRING_TO_NUMBER_H
#define BASE_10_STRING_TO_NUMBER_H

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec-02:
^  - maximum digits and values for unsigned data types
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define def_floatPercision_base10str 6
   /*think 6 decimal digits of percision*/
#define def_maxF_base10str (24 / 3.3333)
   /*floats are always 24 bit for integer, but 8 bit
   `   for mentesa
   */

/*Some tricks I use that will be optimized out at
` compiler time:
*/
#define def_maxUL_base10str (unsigned long) -1
#define def_ulBits_base10str (sizeof(unsigned long) << 3)
#define def_ulDig_base10str (def_ulBits_base10str /3.3333)

#define def_maxUI_base10str (unsigned int) -1
#define def_uiBits_base10str (sizeof(unsigned int) << 3)
#define def_uiDig_base10str (def_uiBits_base10str /3.3333)

#define def_maxUS_base10str (unsigned short) -1
#define def_usBits_base10str (sizeof(unsigned short) << 3)
#define def_usDig_base10str (def_usBits_base10str /3.3333)

#define def_maxUC_base10str (unsigned char) -1
#define def_ucBits_base10str (sizeof(unsigned char) << 3)
#define def_ucDig_base10str (def_ucBits_base10str /3.3333)

/* Logic:
`  - def_maxType_base10str:
`    o (unsigned type) -1:
`      - gives largest number possible for the datatype
`  - def_typeBits_base10str:
`    o sizeof(type) << 3:
`      - gives number of bits (bytes * 8) in data type
`  - def_ulDig_base10str:
`    o def_typeBits_base10str / 3.3333
`      - gives a rough idea of the maximum number of base
`        10 digits for a datatype (ie: 10^x)
*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec-03:
^  - Maximum digits and values for signed data types
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define def_maxSL_base10str (signed long) (def_maxUL_base10str >> 1)
#define def_slBits_base10str ((sizeof(signed long) <<3)-1)
#define def_slDig_base10str (def_slBits_base10str /3.3333)


#define def_maxSI_base10str (signed int) (def_maxUI_base10str >> 1)
#define def_siBits_base10str ((sizeof(signed int) <<3) -1)
#define def_siDig_base10str (def_siBits_base10str /3.3333)


#define def_maxSS_base10str (signed short) (def_maxUS_base10str >> 1)
#define def_ssBits_base10str ((sizeof(signed short)<<3)-1)
#define def_ssDig_base10str (def_ssBits_base10str /3.3333)


#define def_maxSC_base10str (signed char) (def_maxUC_base10str >> 1)
#define def_scBits_base10str ((sizeof(signed char) <<3)-1)
#define def_scDig_base10str (def_scBits_base10str /3.3333)

/* Logic:
`  - def_maxType_base10str:
`    o (signed type) ((unsigned type) -1 >> 1):
`      - sets all bits for type to 1 and then removes the
`        negative bit, giving the largest number possible
`  - def_typeBits_base10str:
`    o (sizeof(type) << 3) - 1:
`      - gives number of bits (bytes * 8) in data type and
`        removes the negative bit (-1)
`  - def_ulDig_base10str:
`    o def_typeBits_base10str / 3.3333
`      - gives a rough idea of the maximum number of base
`        10 digits for a datatype (ie: 10^x)
*/

/*-------------------------------------------------------\
| Fun01: strToUL_base10str
|   - converts a base 10 c-string to an unsigned long
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retUL:
|     o pointer to unisnged long to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retUL to hold the unsigned long
| Note:
|   - only converts till retUL goes out off bounds
\-------------------------------------------------------*/
unsigned int
strToUL_base10str(
   signed char *inStr,
   unsigned long *retUL
);

/*-------------------------------------------------------\
| Fun02: strToUI_base10str
|   - converts a base 10 c-string to an unsigned int
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retUI:
|     o pointer to unisnged int to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retUI to hold the unsigned int
| Note:
|   - only converts till retUI goes out off bounds
\-------------------------------------------------------*/
unsigned int
strToUI_base10str(
   signed char *inStr,
   unsigned int *retUI
);

/*-------------------------------------------------------\
| Fun03: strToUS_base10str
|   - converts a base 10 c-string to an unsigned short
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retUS:
|     o pointer to unisnged short to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retUS to hold the unsigned short
| Note:
|   - only converts till retUS goes out off bounds
\-------------------------------------------------------*/
unsigned int
strToUS_base10str(
   signed char *inStr,
   unsigned short *retUS
);

/*-------------------------------------------------------\
| Fun04: strToUC_base10str
|   - converts a base 10 c-string to an unsigned char
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retUC:
|     o pointer to unisnged char to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retUC to hold the unsigned char
| Note:
|   - only converts till retUC goes out off bounds
\-------------------------------------------------------*/
unsigned int
strToUC_base10str(
   signed char *inStr,
   unsigned char *retUC
);

/*-------------------------------------------------------\
| Fun05: strToSL_base10str
|   - converts base 10 c-string to a signed long
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retSL:
|     o pointer to signed long to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retSL to hold the signed long
| Note:
|   - only converts until retSL overflows
\-------------------------------------------------------*/
unsigned int
strToSL_base10str(
   signed char *inStr,
   signed long *retSL
);

/*-------------------------------------------------------\
| Fun06: strToSI_base10str
|   - converts base 10 c-string to a signed int
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retSI:
|     o pointer to signed int to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retSI to hold the signed int
| Note:
|   - only converts until retSI overflows
\-------------------------------------------------------*/
unsigned int
strToSI_base10str(
   signed char *inStr,
   signed int *retSI
);

/*-------------------------------------------------------\
| Fun07: strToSS_base10str
|   - converts base 10 c-string to a signed short
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retSS:
|     o pointer to signed short to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retSS to hold the signed short
| Note:
|   - only converts until retSS overflows
\-------------------------------------------------------*/
unsigned short
strToSS_base10str(
   signed char *inStr,
   signed short *retSS
);

/*-------------------------------------------------------\
| Fun08: strToSC_base10str
|   - converts base 10 c-string to a signed char
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retSC:
|     o pointer to signed char to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retSC to hold the signed char
| Note:
|   - only converts until retSC overflows
\-------------------------------------------------------*/
unsigned char
strToSC_base10str(
   signed char *inStr,
   signed char *retSC
);

/*-------------------------------------------------------\
| Fun09: strToF_base10str
|   - converts base 10 c-string to float
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retF:
|     o pointer to float to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retF to hold the float
| Note:
|   - only converts until retF overflows
\-------------------------------------------------------*/
unsigned char
strToF_base10str(
   signed char *inStr,
   float *retF
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
