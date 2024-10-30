/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' geneCoord SOF: Start Of File
'   - holds the geneCoord struct and its supporting
'     functions I use to read in gene positions form a paf
'     file
'   o header:
'     - Header guards
'   o .h st: geneCoord
'     - Holds the arrays for the gene coordinates and ids 
'       in an gene coordiantes table
'   o fun02: freeStack_geneCoord
'     - Frees the arrays in a geneCoord structure
'   o fun03: freeHeap_geneCoord
'     - Frees the a heap alloacted geneCoord structure
'   o fun01: init_geneCoord
'     - Initializes a geneCoord structure
'   o fun04: mk_geneCoord
'     - Makes a heap allocated geneCoord structure
'   o fun05: getPaf_geneCoord
'     - Get the id and coordinates for a gene from a paf
'       file
'   o .c fun06: swap_geneCoord
'     - Swaps two array items in a geneCoord structure
'       around
'   o fun07: sort_geneCoord
'     - Sorts the arrays in a genesCoord structure by
'       starting positiion with shell short.
'   o fun08: findStart_geneCoord
'     - Does a binary search by starting coordinate for a
'       potentail gene in a geneCoord structure
'   o fun09: sortName_geneCoord
'     - Sorts the arrays in a genesCoord structure by
'       gene name
'   o fun10: findName_geneCoord
'     - Does a binary search to find an gene name in an
'       gene geneCoord structer (must be sorted by name)
'   o fun11: getCoords_geneCoord
'     - Gets the gene coordinates from an gene coordinates
'       table
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Header guards
\-------------------------------------------------------*/

#ifndef GENE_COORDINATE_STRUCTER_H
#define GENE_COORDINATE_STRUCTER_H

#define def_fileErr_geneCoord 2
#define def_invalidEntry_geneCoord 4
#define def_memErr_geneCoord 64

#define def_lenId_geneCoord 64

/*-------------------------------------------------------\
| ST01: geneCoord
|   - Holds the arrays for the gene coordinates and ids 
|     in an gene coordinates table
\-------------------------------------------------------*/
typedef struct geneCoord{
   signed char (*idStrAry)[def_lenId_geneCoord];
      /*id of gene/section to change*/

   signed char (*refAryStr)[def_lenId_geneCoord];
      /*id of reference to convert coordinates to*/

   unsigned int *startAryUI; /*new starting position*/
   unsigned int *endAryUI;   /*new ending position*/
   unsigned char *dirAryUC;  /*direction*/
}geneCoord;

/*#define getIdIndex_geneCoord(indexUI) ((indexUI) * def_lenId_geneCoord)*/

/*-------------------------------------------------------\
| Fun01: init_geneCoord
|    - Initializes a geneCoord structure
| Input:
|   - geneCoordST:
|     o Pointer to a geneCoord structure to initialize
| Output:
|   - Sets:
|     o idStrAry, startAryUI, and endAryUI to 0
\-------------------------------------------------------*/
void
init_geneCoord(
   struct geneCoord *geneCoordST
);

/*-------------------------------------------------------\
| Fun02: freeStack_geneCoord
|    - Frees the arrays in a geneCoord structure
| Input:
|   - geneCoordST:
|     o Pointer to a geneCoord structure with arrays to
|       free
| Output:
|   - Frees:
|     o idStrAry, startAryUI, and endAryUI (and sets to 0)
\-------------------------------------------------------*/
void
freeStack_geneCoord(
   struct geneCoord *geneCoordST
);

/*-------------------------------------------------------\
| Fun03: freeHeap_geneCoord
|    - Frees the a heap alloacted geneCoord structure
| Input:
|   - geneCoordST:
|     o Pointer to a geneCoord structure to free
| Output:
|   - Frees:
|     o geneCoordST
\-------------------------------------------------------*/
void
freeHeap_geneCoord(
   struct geneCoord *geneCoordST
);

/*-------------------------------------------------------\
| Fun04: mk_geneCoord
|    - Makes a heap allocated geneCoord structure
| Input:
|   - numGenesUI:
|     o Length of the arrays in the new geneCoord
|       structure
| Output:
|   - Returns:
|     o geneCoord structure with allocated arrays
| Note:
|   - Do not use init_geneCoord on the returned structure
\-------------------------------------------------------*/
struct geneCoord *
mk_geneCoord(
   unsigned int numGenesUI
);

/*-------------------------------------------------------\
| Fun05: getPaf_geneCoord
|    - Get the id and coordinates for a gene from a paf
|      file
| Input:
|   - geneCoordST:
|     o Pointer to a geneCoordST structure to hold the
|       new gene
|   - posUI:
|     o Position in the arrays in the geneCoord structure
|       to add the new gene at
|   - typeC:
|     o Alignment type for the extracted gene
|   - pafLineStr:
|     o Line from the paf file with the gene to add
| Output:
|   - Modifies:
|     o All arrays in geneCoordST to hold the new gene
\-------------------------------------------------------*/
void
getPaf_geneCoord(
   struct geneCoord *geneCoordST,
   unsigned int posUI,
   signed char *typeC,
   signed char *pafLineStr
); /*getPaf_geneCoord*/

/*-------------------------------------------------------\
| Fun07: sort_geneCoord
|  - Sorts the arrays in a genesCoord structure by
|    starting positiion with shell short.
| Input:
|  - geneCoordST:
|    o Pointer to geneCoord structure with gene
|      coordinates to sort
|  - startUI:
|    o First element to start sorting at
|  - endUI:
|    o Last element to sort
| Output:
|  - Modifies:
|    o Arrays in geneCoordST to be sorted by the gene
|      starting coordinate (lowest first)
\-------------------------------------------------------*/
void
sort_geneCoord(
   struct geneCoord *geneCoordST,
   unsigned int startUI,
   unsigned int endUI
);

/*-------------------------------------------------------\
| Fun08: findStart_geneCoord
|  - Does a binary search by starting coordinate for a
|    potentail gene in a geneCoord structure
| Input:
|  - geneCoordST:
|    o Pointer to geneCoord structure with starting gene
|      coordinates to search
|  - qryUI:
|    o Starting coordinate (query) to search for in
|      geneCoordST
|  - numGenesUI:
|    o Number of genes in geneCoordST (index 1)
| Output:
|  - Returns:
|    o The index of the starting position
|    o -1 if there was no gene
\-------------------------------------------------------*/
int
findStart_geneCoord(
   struct geneCoord *geneST,
   unsigned int qryUI,
   signed int numGenesSI
);

/*-------------------------------------------------------\
| Fun09: sortName_geneCoord
|  - Sorts the arrays in a genesCoord structure by
|    gene name
| Input:
|  - geneCoordST:
|    o Pointer to geneCoord structure with gene
|      coordinates to sort
|  - startUI:
|    o First element to start sorting at
|  - endUI:
|    o Last element to sort (index 0)
| Output:
|  - Modifies:
|    o Arrays in geneCoordST to be sorted by the gene
|      starting coordinate (lowest first)
\-------------------------------------------------------*/
void
sortName_geneCoord(
   struct geneCoord *geneCoordST,
   unsigned int startUI,
   unsigned int endUI
);

/*-------------------------------------------------------\
| Fun10: findName_geneCoord
|  - Does a binary search to find an gene name in an gene
|    geneCoord structer (must be sorted by name)
| Input:
|  - geneCoordST:
|    o Pointer to geneCoord structure with starting gene
|      coordinates to search
|  - nameStr:
|    o c-string with name to search for
|  - numGenesSI:
|    o Number of genes in geneCoordST (index 1)
| Output:
|  - Returns:
|    o The index of gene with the same name
|    o -1 if there was no gene
\-------------------------------------------------------*/
int
findName_geneCoord(
   struct geneCoord *geneST,
   signed char *nameStr,
   signed int numGenesSI
);

/*-------------------------------------------------------\
| Fun11: getCoords_geneCoord
|  - Gets the gene coordinates from a gene table (tsv)
| Input:
|  - geneTblFileStr:
|    o C-string with name of the gene table file to
|      extract the gene coordinates and names from
|  - numGenesSI:
|    o Will hold the Number of genes extracted
|  - errULPtr:
|    o Will hold the error return value
| Output:
|  - Returns:
|    o Pointer to an sorted geneCoord structure with the 
|      gene coordinates
|    o 0 for errors
|  - Modifies:
|    o numGenesI to have the number of genes (index 0)
|      extracted
|    o errULPtr to hold the error
|      - 0 for no errors
|      - def_fileErr_geneCoord for an file opening error
|      - def_memErr_geneCoord for an memor error
|      - (line_number << 8) | def_invalidEntry_geneCoord
|        for an invalid line in the file
|        o Get the line number with (*errULPtr >> 8)
\-------------------------------------------------------*/
struct geneCoord *
getCoords_geneCoord(
   signed char *geneTblFileStr,
   signed int *numGenesSI, /*Number of genes extracted*/
   unsigned long *errULPtr
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
