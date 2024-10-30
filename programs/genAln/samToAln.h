/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' samToAln SOF: Start Of File
'   - has functions to convert a sam entry and the
'     reference sequence to an alignment
'   o header:
'     - forward declerations and guards
'   o fun01: psamPg_samToAln
'     - prints sam file program entries in the header
'   o fun02: phead_samToAln
'     - prints expanded cigar format header for alignment
'   o fun03: startLine_samToAln
'     - setup start of line in buffer
'   o fun04: paln_samToAln
'     - prints an alignment in expanded cigar format
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - forward declerations and guards
\-------------------------------------------------------*/

#ifndef SAM_FILE_TO_ALIGNMENT
#define SAM_FILE_TO_ALIGNMENT

typedef samEntry samEntry;
typedef seqST seqST;
typedef alnSet alnSet;

/*-------------------------------------------------------\
| Fun01: psamPg_samToAln
|   - prints sam file program entries in the header
| Input:
|   - samSTPtr:
|     o pointer to a sam entry struct with first header
|   - buffStrPtr:
|     o pointer to buffer to use in read sam file headers
|   - lenBuffULPtr:
|     o pointer to unsigned long with length of buffStrPtr
|   - samFILE:
|     o sam file to get header from
|   - outFILE:
|     o file to print header to
| Output:
|   - Modifies:
|     o samSTPtr to have first entry after header
|     o buffStrPtr to be larger (if needed)
|     o lenBuffULPtr to have current buffStrPtr size
|   - Prints:
|     o aln master header to outFILe
|   - Returns:
|     o 0 for no errors
|     o def_EOF_samEntry for eof
|     o def_memErr_samEntry for memory errors
\-------------------------------------------------------*/
signed char
psamPg_samToAln(
   struct samEntry *samSTPtr,
   signed char **buffStrPtr,
   unsigned long *lenBuffULPtr,
   void *samFILE,
   void *outFILE
);

/*-------------------------------------------------------\
| Fun02: phead_samToAln
|   - prints expanded cigar format header for alignment
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct with aligned sequence
|   - refIdStr:
|     o c-string with reference name
|   - scoreSL:
|     o score for the alignment
|   - numAnonSL:
|     o number of anonymous matches
|   - outFILE:
|     o pointer to file to print header to
| Outputs:
|   - Prints:
|     o header for expanded cigar format to outFILE
\-------------------------------------------------------*/
void
phead_samToAln(
   struct samEntry *samSTPtr,
   signed char *refIdStr,
   signed long scoreSL,
   signed long numAnonSL,
   void *outFILE
);

/*-------------------------------------------------------\
| Fun03: startLine_samToAln
|   - setup start of line in buffer
| Input:
|   - pposBl:
|     o 1: add teh positoin to the buffer
|   - refBuffStr:
|     o buffer to set up for reference sequence
|   - qryBuffStr:
|     o buffer to set up for query sequence
|   - cigBuffStr:
|     o buffer to set up for expanded cigar entry
|   - refPosUL:
|     o reference position to add to buffer (i pposBl = 1)
|   - qryPosUL:
|     o query position to add to buffer (i pposBl = 1)
|   - lenLineUS:
|     o pointer to unsigned short to have new buff length
| Output:
|   - Modifies:
|     o refBuffStr to have the reference entry
|     o qryBuffStr to have the query entry
|     o cigBuffStr to have the expand cigar entry
|     o lenbuffUS to length of start of alignment
|     o lenLineUS to length of start of alignment
\-------------------------------------------------------*/
void
startLine_samToAln(
   signed char pposBl,        /*1: print positions*/
   signed char *refBuffStr,
   signed char *qryBuffStr,
   signed char *cigBuffStr,
   unsigned long qryPosUL,    /*qury positoin to add*/
   unsigned long refPosUL,    /*ref position to add*/
   unsigned long *lenLineUL   /*characters in buffer*/
);

/*-------------------------------------------------------\
| Fun04: paln_samToAln
|   - prints an alignment in expanded cigar format
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct with aligned sequence
|   - refSTPtr:
|     o pointer to seqST struct with reference sequence
|   - scoreSL:
|     o score for the alignment
|   - setSTPtr:
|     o pointer to alnSet struct with alignment setSTPtr
|   - outFILE:
|     o pointer to file to print alignment to
| Outputs:
|   - Prints:
|     o header for expanded cigar format to outFILE
\-------------------------------------------------------*/
signed char
paln_samToAln(
   struct samEntry *samSTPtr,
   struct seqST *refSTPtr,
   struct alnSet *setSTPtr,
   void *outFILE
);

#endif

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However,for cases were the public domain is not
:   suitable,such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconvient / not possible,this code is under the
:   MIT license.
: 
: Public domain:
: 
: This is free and unencumbered software released into the
:   public domain.
: 
: Anyone is free to copy,modify,publish,use,compile,
:   sell,or distribute this software,either in source
:   code form or as a compiled binary,for any purpose,
:   commercial or non-commercial,and by any means.
: 
: In jurisdictions that recognize copyright laws,the
:   author or authors of this software dedicate any and
:   all copyright interest in the software to the public
:   domain. We make this dedication for the benefit of the
:   public at large and to the detriment of our heirs and
:   successors. We intend this dedication to be an overt
:   act of relinquishment in perpetuity of all present and
:   future rights to this software under copyright law.
: 
: THE SOFTWARE IS PROVIDED "AS IS\n",WITHOUT WARRANTY OF
:   ANY KIND,EXPRESS OR IMPLIED,INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO
:   EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM,
:   DAMAGES OR OTHER LIABILITY,WHETHER IN AN ACTION OF
:   CONTRACT,TORT OR OTHERWISE,ARISING FROM,OUT OF OR
:   IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
:   DEALINGS IN THE SOFTWARE.
: 
: For more information,please refer to
:   <https://unlicense.org>
: 
: MIT License:
: 
: Copyright (c) 2024 jeremyButtler
: 
: Permission is hereby granted,free of charge,to any
:   person obtaining a copy of this software and
:   associated documentation files (the "Software"),to
:   deal in the Software without restriction,including
:   without limitation the rights to use,copy,modify,
:   merge,publish,distribute,sublicense,and/or sell
:   copies of the Software,and to permit persons to whom
:   the Software is furnished to do so,subject to the
:   following conditions:
: 
: The above copyright notice and this permission notice
:   shall be included in all copies or substantial
:   portions of the Software.
: 
: THE SOFTWARE IS PROVIDED "AS IS\n",WITHOUT WARRANTY OF
:   ANY KIND,EXPRESS OR IMPLIED,INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM,DAMAGES OR OTHER LIABILITY,WHETHER IN
:   AN ACTION OF CONTRACT,TORT OR OTHERWISE,ARISING
:   FROM,OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
