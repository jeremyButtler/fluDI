/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' cigToEqx SOF: Start Of File
'   - functions to convert an regular sam cigar to an eqx
'     cigar
'   o header:
'     - guards and forward declerations
'   o fun01: samCigToEqxCig
'     - convert an sam entry cigar to an eqx cigar
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards and forward declerations
\-------------------------------------------------------*/

#ifndef CIGAR_TO_EQX_CIGAR_H
#define CIGAR_TO_EQX_CIGAR_H

struct samEntry;
struct seqST;

/*-------------------------------------------------------\
| Fun01: samCigToEqxCig 
|   - convert an sam entry cigar to an eqx cigar
| Input:
|   - samSTPtr:
|     o pointer to an samEntry struct to convert cigar
|       into an eqx cigar entry
|   - refStr:
|     o pointer to an reference sequence to use in the
|       conversion
|   - strBuffPtr:
|     o pointer to an c-string buffer to hold the cigar
|       entry temporarly
|   - valArySI:
|     o pointer to an int array to hold the cigar
|       entry temporarly
|   - lenBuffSI:
|     o pointer to current length of strBuffPtr and
|       valArySI
|   - nAsSnpSC:
|     o 1: Treat anymous bases as SNPs always
|     o 0: Treat annoymous bases as matches, if they match
| Output:
|   - Modifies:
|     o cigTypeStr in samSTPtr to be an eqx entry
|     o cigValAryI in samSTPtr to have counts for the
|       eqx entry
|     o strBuffPtr to be larger if needed
|     o valArySI to be larger if needed
|     o lenBuffSI to hold new buffer size if needed
|   - Returns:
|     o 0 for success
|     o 1 for memory error
|       - does not free strBuffPtr (you need to free)
|       - does not free valArySI   (you need to free)
\-------------------------------------------------------*/
signed char
cigToEqx(
   struct samEntry *samSTPtr,
   signed char *refStr,
   signed char **strBuffPtr,
   signed int **valArySI,
   signed int *lenBuffSI,
   signed char nAsSnpSC
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
