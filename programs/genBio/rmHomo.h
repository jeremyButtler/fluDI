/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' rmHomo SOF: Start Of File
'   - removes indels in homopolymers
'   o header:
'     - guards and forward declarations
'   o fun01: indel_rmHomo
'     - remove indels from homopolymers
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards and forward declarations
\-------------------------------------------------------*/

#ifndef REMOVE_INDELS_IN_HOMOPOLYMERS_H
#define REMOVE_INDELS_IN_HOMOPOLYMERS_H

struct samEntry;

/*-------------------------------------------------------\
| Fun01: indel_rmHomo
|   - remove indels from homopolymers
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to sequence to remove
|       indels from homopolymers
|   - refStr:
|     o c-string with reference sequence in samSTPtr is
|       mapped to
|   - minHomoSI:
|     o minimum homopolymer length to remove indel from
|   - maxIndelSI:
|     o maximum indel size to remove
|   - maskSC:
|     o base to mask deletions with
|     o if 0 if(! maskSC), then uses reference base
|   - seqBuffStrPtr:
|     o c-string pointer to buffer to use in copy sequence
|     o size to avoid reallocation of buffer
|       * sequence length + number deletions + 8
|   - sizeSeqUIPtr:
|     o unsigned long pointer to with size of
|       seqBuffStrPtr
|   - qBuffStrPtr:
|     o c-string pointer to buffer to use in copy q-score
|       entry
|     o size to avoid reallocation of buffer
|       * sequence length + number deletions + 8
|   - sizeQUIPtr:
|     o unsigned long pointer to with size of qBuffStrPtr
| Output:
|   - Modifies:
|     o seqStr in samSTPtr to have no homopolymer indels
|     o seqBuffStrPtr and qBuffStrPtr to be resized (or
|       allocated) to sequence length + deletions
|     o sizeQUIPtr and sizeSeqUIPtr to have their buffer
|       size if seqBuffStrPtr or qBuffStrPtr are resized
|   - Returns:
|     o 0 for no errors
|     o 1 for memory errors
\-------------------------------------------------------*/
signed char
indel_rmHomo(
   struct samEntry *samSTPtr,/*sequence with indels*/
   signed char *refStr,      /*reference sequence*/
   signed int minHomoSI,     /*min homoplymer length*/
   signed int maxIndelSI,    /*indel to large to remove*/
   signed char maskSC,       /*base to mask dels with*/
   signed char **seqBuffStrPtr, /*sequence buffer*/
   unsigned int *sizeSeqUIPtr,  /*size of seqBuffStrPtr*/
   signed char **qBuffStrPtr,   /*q-score buffer*/
   unsigned int *sizeQUIPtr     /*size of qBuffStrPtr*/
);
   /*having users provide temporary buffers was choosen,
   `  because I need a temporary buffer to get sam
   `  entries. Often times this temporary buffer will
   `  be able to be easily split into two separate
   `  buffers of the correct size. So, no point adding
   `  more memory and waste time allocating memory.
   `  This does come at the cost of having to always
   `  copy the sequence though.
   */

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
: Copyright (c) 2025 jeremyButtler
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
