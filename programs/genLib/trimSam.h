/*########################################################
# Name: trimSam
#   - Holds functions for trimming softmasked regoins at
#     the start and end of sequences in an sam file
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - header guards
'   o fun01 trimSeq_trimSam:
'     o Trim soft mask regions off end of sam entry
'   o fun02: trimCoords_trimSam
'     - Trim an sam file entry by coordinates
'   o fun03 trimReads_trimSam:
'     o Trims soft mask regions for all reads with a
'       sequence in a sam file
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - header guards
\-------------------------------------------------------*/

#ifndef TRIMSAM_H
#define TRIMSAM_H

typedef struct samEntry samEntry;

#define def_rangeErr_trimSam 1
#define def_memErr_trimSam 2
#define def_fileErr_trimSam 4
#define def_header_trimSam 8
#define def_noMap_trimSam 16
#define def_noSeq_trimSam 32

/*-------------------------------------------------------\
| Fun01: trimSeq_trimSam
| Use:
|  - Trims off the soft masked regions of a sam entry
| Input:
|  - samVoidST:
|    o samEntry struct with sam entry to trim. Named
|      void so I can cast it as a samEntry to samST
| Output:
|  - Returns:
|    o 0 if suceeded
|    o def_header_trimSam if header (invalid and ignored)
|    o def_noMap_trimSam if unmapped read (no reference)
|    o def_noSeq_trimSam if no sequence line
|  - Modifies:
|    o Trims cigar, sequence, & q-score entries in samST.
\-------------------------------------------------------*/
signed char
trimSeq_trimSam(
   struct samEntry *samST
);

/*-------------------------------------------------------\
| Fun02: trimCoords_trimSam
|   - Trim an sam file entry by coordinates
| Input:
|   - samSTPtr:
|     o Pointer to an sam entry structure with an read to
|       trim
|   - startSI:
|     o Singed integer with the starting coordinate
|   - endSI:
|     o Singed integer with the ending coordinate
|     o < 1 is treated as do not trim
| Output:
|   - Modifies:
|     o samSTPtr to be trimmed
|   - Returns:
|     o 0 for no errors
|     o 1 for coordinates out of range
\-------------------------------------------------------*/
signed char
trimCoords_trimSam(
   struct samEntry *samSTPtr,
   signed int startSI,
   signed int endSI
);

/*-------------------------------------------------------\
| Fun03: trimReads_trimSam
| Use:
|  - Goes though sam file and calls trimSeq_trimSam for each
|    entry
| Input:
|  - samFILE:
|    o Pointer to FILE. Sam file with entries to trim
|  - outFILE:
|    o File to write the trimmed sam file entries to
|  - keepUmapBl:
|    o Also print out entries with unmapped reads
| Output:
|  - Returns:
|    o 0 if suceeded
|    o def_memErr_trimSam for memory errors
|  - Prints:
|    o trimmed sam entries with sequences to outFILE, but
|      ignores sam entries without sequences
\-------------------------------------------------------*/
signed char
trimReads_trimSam(
    void *samFILE,           /*sam file to trim*/
    void *outFILE,           /*file to store output*/
    signed char keepUmapBl   /*1: keep unmapped reads*/
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
