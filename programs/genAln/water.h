/*########################################################
# Name water
# Use:
#  o Holds functions to do a Smith Waterman pairwise
#    alignment
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'   o header:
'     - forward declerations and guards
'   o .c fun01: scoreIndel_water
'     - gets the indel score for a water alignment
'   o fun01 water:
'     - run a memory efficent Waterman Smith alignment on
'       input sequences
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - forward declerations and guards
\-------------------------------------------------------*/

#ifndef WATERMAN_SMITH_ALIGNMENT_H
#define WATERMAN_SMITH_ALIGNMENT_H

typedef struct seqST seqST;
typedef struct alnSet alnSet;
typedef struct dirMatrix dirMatrix;

#define def_memErr_water 1

/*-------------------------------------------------------\
| Fun01: water
|   - run a Waterman Smith alignment on input sequences
| Input:
|   - qrySTPtr:
|     o pointer to seqST with query sequence 
|       - qrySTPtr->offsetUL; first query base to align
|       - qrySTPtr->endAlnUL; last query base to align
|   - refSTPtr:
|     o pointer to seqST with reference sequence 
|       - refSTPtr->offsetUL; 1st reference base to align
|       - refSTPtr->endAlnUL; last reference base to align
|   - matrixSTPtr:
|     o pointer to dirMatrix to use for the alingment
|   - alnSet:
|     o pointer to alnSet with alignment settings
| Output:
|  - Modifies:
|    o errSC in matrixSTPtr to def_memErr_water for
|      memory errors
|    o allocates memory for dirMatrixSC and scoreAryUL
|      if they are to small
|    o updates lenMatrixUL and lenScoreUL if dirMatrixSC
|      or scoreAryUL are resized
|  - Returns:
|    o score for alignment
\-------------------------------------------------------*/
signed long
water(
   struct seqST *qrySTPtr, /*query sequence and data*/
   struct seqST *refSTPtr, /*ref sequence and data*/
   struct dirMatrix *matrixSTPtr, /*direction matrix*/
   struct alnSet *settings     /*settings for alignment*/
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
