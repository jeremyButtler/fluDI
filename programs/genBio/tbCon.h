/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' tbCon SOF: Start Of File
'  - holds functions to make a quick consensus, more like
'    ivar. this is a bit memory intensive,
'    (40 + # insertions * (17) + bases in insertions)
'    bytes 
'   o header:
'     - Included libraries
'   o .h st01: insBase
'     - structure to hold an list of insetions for a
'       single position
'   o .h st02: conNt_tbCon
'     - structure to hold an single base in a sequence
'   o .h st03: set_tbCon
'     - holds the settings for tbCon
'   o fun01: initInsBase
'     - initializes an initBase structure
'   o fun02: mkIns_tbCon
'     - makes and initializes an ins_tbCon structure
'   o fun03: freeHeap_ins_tbCon
'     - Frees an heap allocated insBase structure
'   o fun04: blank_conNt_tbCon
'     - sets all variables in conNt_tbCon structure to 0
'   o fun05: init_conNt_tbCon
'     - initializes a conNt_tbCon structure
'   o fun06: mk_conBase_tbCon
'     - Makes an initialized, heap allocated conNt_tbCon
'       structure
'   o fun07: freeStack_conNt_tbCon
'     - frees heap allocated variables in a conNt_tbCon
'       structure (also initializes structure)
'   o fun08: freeHeap_conNt_tbCon
'     - frees a conNt_tbCon structure
'   o fun09: freeHeapAry_conNt_tbCon
'     - frees array of heap allocated conNt_tbCon structs
'   o fun10: init_set_tbCon
'     - initialize a set_tbCon struct to default settings
'   o fun11: freeStack_set_tbCon
'     - frees variables inside a set_tbCon struct
'   o fun12: addRead_tbCon
'     - adds read to a conNt_tbCon array
'   o fun13: collapse_tbCon
'     - Collapses a conNt_tbCon array into an array of
'       samEntry structures
'   o fun14: noFragCollapse_tbCon
'     - collapses a conNt_tbCon struct array into a single
'       samEntry struct (low read depth is masked)
'   o fun15: pvar_tbCon
'     - print entries in a conNt_tbCon structure array
'       that are above the minimum read depth
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Defined variables and header guards
\-------------------------------------------------------*/

#ifndef TB_CONSENSUS_H
#define TB_CONSENSUS_H

typedef struct samEntry samEntry;

/*Settings for consensus buildiing*/

/*-------------------------------------------------------\
| ST01: ins_tbCon
|   - Structure to hold an list of insetions for a single
|     position
\-------------------------------------------------------*/
typedef struct ins_tbCon{
   signed char *insStr;  /*c-string with insertion*/
   signed int lenInsSI;  /*length of insertion*/
   signed int numHitsSI; /*number reads supporting a ins*/

   struct ins_tbCon *nextIns; /*Next insertion*/
}ins_tbCon;

/*-------------------------------------------------------\
| ST02: conNt_tbCon
|   - Structure to hold an single base in a sequence
\-------------------------------------------------------*/
typedef struct conNt_tbCon{
   signed int numASI;     /*number reads supporting a*/
   signed int numTSI;     /*number reads supporting t*/
   signed int numCSI;     /*number reads supporting c*/
   signed int numGSI;     /*number reads supporting g*/
   signed int numDelSI;   /*number reads supporting del*/

   signed int numInsSI;   /*number reads supporting ins*/
   struct ins_tbCon *insList; /*List of insertions*/

   signed int totalNtSI;  /*reads supporting position*/
   signed int ntKeptSI;   /*reads above min q-score*/
}conNt_tbCon;

/*-------------------------------------------------------\
| ST03: set_tbCon
|   - Structure to hold settings for tbCon
\-------------------------------------------------------*/
typedef struct set_tbCon{
   signed int minQSI;       /*min q-score to keep a base*/
   signed int minInsQSI;    /*min q-score to keep ins*/
   unsigned char minMapqUC; /*min mapq to keep a read*/

   /*min percent/depth/length to keep an match/mutation*/
   signed int minLenSI;     /*min fragment length*/
   signed int minDepthSI;   /*min depth to keep base*/
   float minPercSnpF;  /*min % support to keep snp/match*/
   float minPercDelF;  /*min % support to keep insertion*/
   float minPercInsF;  /*min % support to keep deletion*/

   /*to mask low quality/support cosnesus bases with*/
   signed char maskSC;      /*character to mask with*/

   /*min percent/depth to print an snp/match,ins,or del*/
   signed int minPrintDepthSI;
   float printMinSupSnpF;
   float printMinSupInsF;
   float printMinSupDelF;
}set_tbCon;

/*-------------------------------------------------------\
| Fun02: mkIns_tbCon
|   - makes and initializes an ins_tbCon structure
| Input:
| Output:
|   - Returns:
|     o pointer to initialized ins_tbCon structure
|     o 0 for memory errors
\-------------------------------------------------------*/
struct ins_tbCon *
mkIns_tbCon(
);

/*-------------------------------------------------------\
| Fun03: freeHeap_ins_tbCon
|   - frees heap allocated insBase structure
| Input:
|   - insSTPtr:
|     o pointer to insBase structure to free
| Output:
|   - Frees:
|     o insSTPtr
|   - Returns:
|     o next insBase structure (nextIns) in insBase list
\-------------------------------------------------------*/
struct ins_tbCon *
freeHeap_ins_tbCon(
   struct ins_tbCon *insSTPtr
);

/*-------------------------------------------------------\
| Fun04: blank_conNt_tbCon
|   - sets all variables in conNt_tbCon structure to 0
| Input:
|   - conNtSTPtr:
|     o pointer to conNt_tbCon structure to blank
| Output:
|   - Modifies:
|     o all non-pointer values in conNtSTPtr to be 0
\-------------------------------------------------------*/
void
blank_conNt_tbCon(
   struct conNt_tbCon *conNtSTPtr
);

/*-------------------------------------------------------\
| Fun05: init_conNt_tbCon
|   - initializes a conNt_tbCon structure
| Input:
|   - conNtSTPtr:
|     o pointer to an conNt_tbCon structure to initialize
| Output:
|   - Modifies:
|     o all values (including pointers) in conNtSTPtr to
|       be 0
\-------------------------------------------------------*/
void
init_conNt_tbCon(
   struct conNt_tbCon *conNtSTPtr
);

/*-------------------------------------------------------\
| Fun06: mk_conBase_tbCon
|   - makes initialized, heap allocated conNt_tbCon struct
| Input:
| Output:
|   - Returns:
|     o pointer to an initialized conNt_tbCon structure
|     o 0 for memory errors
\-------------------------------------------------------*/
struct conNt_tbCon *
mk_conBase_tbCon(
);

/*-------------------------------------------------------\
| Fun07: freeStack_conNt_tbCon
|   - frees heap allocated variables in a conNt_tbCon
|     structure (also initializes structure)
| Input:
|   - conNtSTPtr:
|     o pointer to conNt_tbCon struct to free variables
| Output:
|   - Frees:
|     o frees the conNt_tbConST->insList 
|   - Sets:
|     o all values in conNt_tbConSTPTr to 0
\-------------------------------------------------------*/
void
freeStack_conNt_tbCon(
   struct conNt_tbCon *conNtSTPtr
);
/*-------------------------------------------------------\
| Fun08: freeHeap_conNt_tbCon
|   - Frees a conNt_tbCon structure
| Input:
|   - conNtSTPtr:
|     o pointer toa conNt_tbCon structure to free
| Output:
|   - Frees:
|     o frees conNt_tbConST
\-------------------------------------------------------*/
void
freeHeap_conNt_tbCon(
   struct conNt_tbCon *conNtSTPtr
);

/*-------------------------------------------------------\
| Fun09: freeHeapAry_conNt_tbCon
|   - frees array of heap allocated conNt_tbCon structs
| Input:
|   - conNtAryST:
|     o pointer conNt_tbCon struct array to free
|   - lenArySI:
|     o number of conNt_tbCon structurs in array
| Output:
|   - Frees:
|     o array in conNtSTPtr
\-------------------------------------------------------*/
void
freeHeapAry_conNt_tbCon(
   struct conNt_tbCon *conNtAryST,
   signed int lenArySI
);

/*-------------------------------------------------------\
| Fun09: freeHeap_conNt_tbConAry
|   - Frees an array of conNt_tbCon structures (assumes heap)
| Input:
|   - conNtSTPtr:
|     o Pointer to an pointer to the start of an conNt_tbCon
|       structure array to free
|   - lenAryI:
|     o Number of conNt_tbCon structurs in the array
| Output:
|   - Frees:
|     o The array in conNtSTPtr
|   - Sets:
|     o conNt_tbConST to 0
\-------------------------------------------------------*/
void
freeHeap_conNt_tbConAry(
   struct conNt_tbCon **conNtSTPtrPtr,
   int lenAryI
);

/*-------------------------------------------------------\
| Fun10: init_set_tbCon
|   - initialize a set_tbCon struct to default settings
| Input:
|   - setSTPtr:
|     o Pointer to the set_tbCon structure with settings
|       for tbCon to set to defaults
| Output:
|   - Modifies:
|     o set_tbConSTPtr to default settings
\-------------------------------------------------------*/
void
init_set_tbCon(
   struct set_tbCon *setSTPtr
);

/*-------------------------------------------------------\
| Fun11: freeStack_set_tbCon
|   - frees variables inside a set_tbCon struct
| Input:
|   - setSTPtr:
|     o pointer to the set_tbCon struct
| Output:
|   - Nothing:
|     o here for future use
\-------------------------------------------------------*/
void
freeStack_set_tbCon(
   struct set_tbCon *setSTPtr
);

/*-------------------------------------------------------\
| Fun12: addRead_tbCon
|   - adds read to a conNt_tbCon array
| Input:
|   - sameEntrySTPtr:
|     o pointer to a samEntry struct with read to add
|   - conNtAryST:
|     o pointer to conNt_tbCon struct array to add read to
|     o the array is resize if to small
|   - lenRefUI:
|     o length of the reference sequence; use 0 for unkown
|     o is updated as conNt_tbConSTAry expands
|   - settings:
|     o pointer to set_tbCon struct with settings
| Output:
|   - Modifies:
|     o conNtAryST to have new read
|     o lenRefUI to have new conNtAryST size if conNtAryST
|       is resized
|   - Returns:
|     o for no errors
|     o def_header_tbConDefs for sam file header entry
|     o def_noMap_tbConDefs for an unmapped read
|     o def_noSeq_tbConDefs for a read missing a sequence
|     o def_memErr_tbConDefs for memory errors
\-------------------------------------------------------*/
signed char
addRead_tbCon(
   samEntry *samSTPtr,       /*read to add to consensus*/
   struct conNt_tbCon *conNtAryST[], /*Consensus array*/
   unsigned int *lenRefUI,       /*Length of reference*/
   struct set_tbCon *settings
);

/*-------------------------------------------------------\
| Fun13: collapse_tbCon
|   - collapses a conNt_tbCon strrucvt array into an array
|     of samEntry structs
| Input:
|   - conNtAryST:
|     o pointer to a conNt_tbCon struct array to collapse
|   - lenConAryUI:
|     o length of the conNt_tbCon array
|   - lenSamArySI:
|     o pointer to integer to have number of samEntry
|       structures made
|   - refIdStr:
|     o c-string with referernce sequence name
|   - settings:
|     o pointer to set_tbCon struct with consensus
|       building settings
|   - errSC:
|     o pointer to character to hold error type
| Output:
|   - Modifies:
|     o lenSamST to hold the returned samEntry arrray
|       length
|     o errSC:
|       - 0 for no errors
|       - def_noSeq_tbConDefs if could not build consensus
|       - def_memErr_tbConDefs if had memory error
|   - Returns:
|     o array of sam entry structures with consensus 
|       fragments.
|       - idealy this would be one structure, but for
|         amplicons it will be an array
|     o 0 for memory erors
\-------------------------------------------------------*/
samEntry *
collapse_tbCon(
   struct conNt_tbCon conNtAryST[], /*to collapse*/
   unsigned int lenConAryUI,   /*length of consensus*/
   signed int *lenSamArySI,    /*set to out array length*/
   signed char *refIdStr,      /*name of reference seq*/
   struct set_tbCon *settings, /*settings for collapsing*/
   signed char *errSC          /*error reports*/
);

/*-------------------------------------------------------\
| Fun14: noFragCollapse_tbCon
|   - collapses a conNt_tbCon struct array into a single
|     samEntry struct (low read depth is masked)
| Input:
|   - conNtAryST:
|     o pointer to a conNt_tbCon struct array to collapse
|   - lenConAryUI:
|     o length of the conNt_tbCon array
|   - refIdStr:
|     o c-string with referernce sequence name
|   - numMaskUIPtr:
|     o pointer to unsigned int to have number of bases
|       masked in consensus
|   - settings:
|     o pointer to set_tbCon struct with consensus
|       building settings
|   - errSC:
|     o pointer to character to hold error type
| Output:
|   - Modifies:
|     o lenSamST to hold the returned samEntry arrray
|       length
|     o numMaskUIPtr to have number of masked bases
|     o errSC:
|       - 0 for no errors
|       - def_noSeq_tbConDefs if could not build consensus
|       - def_memErr_tbConDefs if had memory error
|   - Returns:
|     o array of sam entry structures with consensus 
|       fragments.
|       - idealy this would be one structure, but for
|         amplicons it will be an array
|     o 0 for memory erors
\-------------------------------------------------------*/
samEntry *
noFragCollapse_tbCon(
   struct conNt_tbCon conNtAryST[], /*to collapse*/
   unsigned int lenConAryUI,   /*length of consensus*/
   signed char *refIdStr,      /*name of reference seq*/
   unsigned int *numMaskUIPtr, /*# bases masked*/
   struct set_tbCon *settings, /*settings for collapsing*/
   signed char *errSC          /*error reports*/
);

/*-------------------------------------------------------\
| Fun15: pvar_tbCon
|   - Print entries in a conNt_tbCon structure array that
|     are above the minimum read depth
| Input:
|   - conNtAryST:
|     o pointer to a conNt_tbCon structure array to print
|   - lenConAryUI:
|     o length of the conNt_tbCon array
|   - minDepthSI:
|     o minimum read depth to print out a position
|   - refIdStr:
|     o c-string with referernce sequence name
|   - outFILE:
|     o c-string with name of file to print everything to
|     o will append to file if already exists
| Output:
|   - Prints:
|     o entries in conNtAryST to outFILE
|   - Returns:
|     o 0 for success
|     o def_fileErr_tbConDefs for file errors
\-------------------------------------------------------*/
char
pvar_tbCon(
   struct conNt_tbCon conNtAryST[], /*consensus array*/
   unsigned int lenConAryUI,     /*length of conNtAryST*/
   signed char *refIdStr,        /*name of reference seq*/
   struct set_tbCon *settings,   /*settings for printing*/
   signed char *outStr           /*file to print to*/
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
