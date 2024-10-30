/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' ampDepth SOF: Start Of File
'   - Had functions for updating and printing an read
'     depth histogram (array of integers)
'   o header:
'     - foward declarations and guards
'   o fun01: addRead_ampDepth
'     - adds a read to an amp depth arrary (histogram)
'   o fun02: phead_ampDepth
'     - prints the header for an ampDepth tsv file
'   o fun03: phist_ampDepth
'     - prints an histogram of read depths as an tsv with
'       individual genes for ampDepth
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - foward declarations and guards
\-------------------------------------------------------*/

#ifndef AMPLICON_DEPTH_H
#define AMPLICON_DEPTH_H

typedef struct samEntry samEntry;
typedef struct geneCoord geneCoord;

/*-------------------------------------------------------\
| Fun01: addRead_ampDepth
|   - adds a read to an amp depth arrary (histogram)
| Input:
|   - samSTPtr:
|     o pointer to an samEntry structure with a read to
|       add to the histogram
|   - geneCoordSTPtr:
|     o pointer to an array of geneCoord structures with
|       gene coordinates for each gene of interest
|   - startSeqUI:
|     o first base in target region of reference 
|   - endSeqUI:
|     o last base in target region of reference
|   - numOffTargSI:
|     o number of reads not mapping to a gene coordiante,
|       but are mapped to the genome
| Output:
|   - Modifies:
|     o ampDepthIAry to be incurmented for each reference
|       base in the read
\-------------------------------------------------------*/
void
addRead_ampDepth(
   struct samEntry *samSTPtr,
   unsigned int startSeqUI,  /*first reference coorinate*/
   unsigned int endSeqUI,    /*last reference coordinate*/
   signed int *depthArySI,
   signed int *numOffTargSI
);

/*-------------------------------------------------------\
| Fun02: phead_ampDepth
|   - prints the header for an ampDepth tsv file
| Input:
|   - outFILE:
|     o file to print the header to
| Output:
|   - Prints:
|     o ampDepth header to outFILE
\-------------------------------------------------------*/
void
phead_ampDepth(
   void *outFILE
);

/*-------------------------------------------------------\
| Fun03: phist_ampDepth
|   - prints an histogram of read depths as an tsv with
|     individual genes for ampDepth
| Input:
|   - histArySI:
|     o integer array with the histogram to print out
|   - minDepthSI:
|     o integer with the min depth to keep an histogram
|       entry
|   - geneCoordSTPtr:
|     o pointer to an geneCoord structure array with the
|       locations of each gene
|   - numGenesSI:
|     o integer with number of genes in genesST
|   - offTargSI:
|     o number off of target reads to print out
|   - noMapSI:
|     o number off of unmapped reads to print out
|   - extraColStr:
|     o c-string with the first column to print out
|   - outVoidFILE:
|     o FILE pointer with output file to print to
| Output:
|   - Prints:
|     o stats for each gene to outVoidFILE
\-------------------------------------------------------*/
void
phist_ampDepth(
   signed int *histArySI,
   signed int minDepthSI,
   struct geneCoord *geneCoordSTPtr,
   signed int numGenesSI,
   signed int offTargSI,
   signed int noMapSI,
   signed char *extraColStr,
   void *outVoidFILE
);

#endif

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconveint / not possible, this code is under the
:   MIT license
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
