/*########################################################
# Name: genMath.h
#   - Has my math functions I use. These are often
#     branchless (noBrach) variations.
# Note:
#   Most of the branchless code is from or modified from:
#   - https://graphics.stanford.edu/~seander/bithacks.html
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o fun01: hibit2_genMath
'     - Does a branchles method to find the higest bit in
'       the element (integer log base 2)
'   o fun02: ab_genMath
'     - Does a branchless absolute value
'   o .h fun03: max_genMath
'     - Find the maximum value (branchless)
'   o .h fun04: ifmax_genMath
'     - Set a value (ret) to a value based on which value
'       is greater.
'   o .h fun05: min_genMath
'     - Find the Minimum value (branchless)
'   o .h fun06: ifMin_genMath
'     - Set a value (ret) to a value based on which value
'       is less than.
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

#include "genMath.h"

/*-------------------------------------------------------\
| Fun01: hibit2_genMath
|   - Does a branchles method to find the higest bit in
|     the element (integer log base 2)
| Input:
|   - numI:
|     o Integer (non-long) number to get the highest bit
|       for
| Output:
|   - Returns:
|     o the position of the highest bit
\-------------------------------------------------------*/
unsigned int
hibit2_genMath(
   unsigned int numUI
){
   unsigned int retUI = 0;
   unsigned int shiftUI = 0;

   retUI = (numUI > 0xFFFF) << 4;
   numUI >>= retUI;

   shiftUI = (numUI > 0xFF) << 3;
   numUI >>= shiftUI;
   retUI |= shiftUI;

   shiftUI = (numUI > 0xF) << 2;
   numUI >>= shiftUI;
   retUI |= shiftUI;

   shiftUI = (numUI > 0x3) << 1;
   numUI >>= shiftUI;
   retUI |= shiftUI;

   return retUI | (numUI >> 1);
} /*hibit2_genMath*/

/*-------------------------------------------------------\
| Fun02: ab_genMath
|   - does a branchless absolute value
| Input:
|   - numI:
|     o integer (non-long) to get the absolute value to
| Output:
|   - Returns:
|     o absolute value of numI
\-------------------------------------------------------*/
signed int
ab_genMath(
   signed int numSI
){
   signed int maskSI =
      numSI >> ( (sizeof(signed int) * 8) - 1 );

   return (numSI + maskSI) ^ maskSI;
} /*ab_genMath*/

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
