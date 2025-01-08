/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' edDist SOF: Start Of File
'   - has functions to find edit distances (does some
'     filtering)
'   o header:
'     - included libraries
'   o fun01: blank_res_edDist
'     - blanks (sets 0) values in a res_edDist struct
'   o fun02: init_res_edDist
'     - initializes and res_edDist struct
'   o fun03: freeStack_res_edDist
'     - frees heap varaibles in a res_edDist struct
'   o fun04: freeHeap_res_edDist
'     - frees heap allocated res_edDist struct
'   o .c fun05: isTransNt_edDist
'     - checks to see if query base is a transversion
'     - not used but want to keep code
'   o fun06: readCmpDist_edDist
'     - gets edit distances between the query & reference
'   o fun07: dist_edDist
'     - gets edit distances for reference
'   o fun08: percDist_edDist
'     - finds probablity of reads matching by chance
'   o fun09: addReadToDepth_edDist
'     - adds a read to a depth profile
'   o fun10: mkDepthProfile_edDist
'     - finds depth of every base in samEntry reference
'   o fun11: phead_edDist
'     - prints header for edDist output tsv
'   o fun12: pdist_edDist
'     - prints edit distance tsv line for edDist
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards and foward declerations
\-------------------------------------------------------*/

#ifndef SAM_FILE_CLUSTERING_H
#define SAM_FILE_CLUSTERING_H

struct samEntry;
struct seqST;

#define def_memErr_edDist -16
#define def_fileErr_edDist -32

#define def_noOverlap_edDist -1
#define def_diffRef_edDist -2
#define def_noSeq_edDist -4
#define def_noMap_edDist -8

/*-------------------------------------------------------\
| ST01: res_edDist
|   - holds results and windows for edDist struct 
\-------------------------------------------------------*/
typedef struct res_edDist
{
   unsigned int *depthAryUI;   /*depth array for profile*/
   unsigned long sizeDepthUL;  /*size of depthAryUI*/

   unsigned int probUI;    /*percent dist/error*/
   unsigned int overlapUI;     /*overlap between reads*/
   unsigned int startUI;       /*first shared base*/
   unsigned int numIndelUI;    /*number indels kept*/
   unsigned int indelEventsUI; /*number indel events*/

   unsigned int maxWinDistUI;  /*maximum dist in window*/
   unsigned int probMaxWinUI;  /*max window % dist/error*/
   unsigned int minWinDistUI;  /*minimum dist in window*/
   float avgWinDistF;          /*average dist in window*/

   unsigned int numFullWinUI;  /*number of full windows*/
   unsigned int extraWinNtUI; 
      /*number bases that did not fit in a window*/
   unsigned int lastWinDiffUI;
      /*number differences in last window*/

   signed long edDistSL;       /*edit distance*/
}res_edDist;

/*-------------------------------------------------------\
| Fun01: blank_res_edDist
|   - blanks (sets 0) values in a res_edDist struct
| Input:
|   - resSTPtr
|     o pointer to res_edDist struct to blank
| Output:
|   - Modifies:
|     o all values in resSTPtr, except depthAryUI and
|       sizeDepthUI to be 0
\-------------------------------------------------------*/
void
blank_res_edDist(
   struct res_edDist *resSTPtr
);

/*-------------------------------------------------------\
| Fun02: init_res_edDist
|   - initializes and res_edDist struct
| Input:
|   - resSTPtr
|     o pointer to res_edDist struct to initialize
| Output:
|   - Modifies:
|     o all values in resSTPtr to be 0
\-------------------------------------------------------*/
void
init_res_edDist(
   struct res_edDist *resSTPtr
);

/*-------------------------------------------------------\
| Fun03: freeStack_res_edDist
|   - frees heap varaibles in a res_edDist struct
| Input:
|   - resSTPtr
|     o pointer to res_edDist struct with vars to free
| Output:
|   - Frees:
|     o depthAryUI; the calls init_res_edDist (fun02)
\-------------------------------------------------------*/
void
freeStack_res_edDist(
   struct res_edDist *resSTPtr
);

/*-------------------------------------------------------\
| Fun04: freeHeap_res_edDist
|   - frees heap allocated res_edDist struct
| Input:
|   - resSTPtr
|     o pointer to res_edDist struct to free
| Output:
|   - Frees:
|     o resSTPtr (you must set to 0)
\-------------------------------------------------------*/
void
freeHeap_res_edDist(
   struct res_edDist *resSTPtr
);

/*-------------------------------------------------------\
| Fun06: readCmpDist_edDist
|   - gets edit distances between two reads mapped to the
|     same reference
|   - deletions and insertions are only counted if they
|     execed a minimum length.
| Input:
|   - qrySTPtr:
|     o pointer to samEntry structure with read (query) to
|       find the edit distance for
|   - refSTPtr:
|     o pointer to samEntry structure with reference to
|       compare query (qrySTPtr) to
|   - indelLenUI:
|     o minimum length for a indel to count as an event
|   - minQUC:
|     o minimum q-score to not discard an snp
|   - minOverlapF:
|     o minimum percent overlap to score
|     o use 0 or negative value to not check overlap
|       - still checks to see reads share bases
|   - minDepthUI:
|     o minimum depth for reference base to keep
|       difference if using a depth profile (depthAryUI)
|   - depthAryUI:
|     o unsigned int array with read depths for each
|       reference base (use 0 to not use)
|   - winSizeUI:
|     o size of window for window error rates
|     o windows shifts by size every time
|   - resSTPtr:
|     o pointer to res_edDist struct to have results
| Output:
|   - Modifies:
|     o numIndelUI in resSTPtr to have number kept indels
|     o indelEventsUI in resSTPtr to have number of indel
|       events
|     o overlapUI in resSTPtr to hold aligned length
|     o maxWinDistUI in resSTPtr to hold maximum number
|       of errors in all windows
|     o minWinDistUI in resSTPtr to hold minimum number
|       of errors in all windows
|     o avgWinDistF in resSTPtr to hold mean number of
|       errors in all widows
|     o edDistSL in resSTPtr to have edit distance
|   - Returns:
|     o edit distance between query and ref
|     o negative value if could not find distance
|       - def_noOveralp_edDist if reads to not overlap
|       - def_diffRef_edDist if reads have different ref
|       - def_noSeq_edDist if one read is missing sequence
|       - def_noMap_edDist if one read is unmapped
\-------------------------------------------------------*/
signed long
readCmpDist_edDist(
   struct samEntry *qrySTPtr, /*read for edit distance*/
   struct samEntry *refSTPtr, /*ref to compare*/
   unsigned int indelLenUI,   /*min indel length*/
   unsigned char minQUC,      /*min Q-score (snp/ins)*/
   float minOverlapF,         /*min % overlap*/
   unsigned int minDepthUI,   /*min depth if profiling*/
   unsigned int winSizeUI,    /*size of one window*/
   struct res_edDist *resSTPtr/*results/depth profile*/
);

/*-------------------------------------------------------\
| Fun07: dist_edDist
|   - gets edit distances for mapped reference
|   - deletions and insertions are only counted if they
|     execed a minimum length.
| Input:
|   - qrySTPtr:
|     o pointer to samEntry structure with read (query) to
|       find the edit distance for
|   - refSTPtr:
|     o reference sequence, if provided checks 'M' cases
|     o use 0 to not check
|   - indelLenUI:
|     o minimum length for a indel to count as an event
|   - minQUC:
|     o minimum q-score to not discard an snp
|   - winSizeUI:
|     o size of window for window error rates
|     o windows shifts by size every time
|   - resSTPtr:
|     o pointer to res_edDist struct to have results
| Output:
|   - Modifies:
|     o numIndelUI in resSTPtr to have number kept indels
|     o indelEventsUI in resSTPtr to have number of indel
|       events
|     o overlapUI in resSTPtr to hold aligned length
|     o maxWinDistUI in resSTPtr to hold maximum number
|       of errors in all windows
|     o minWinDistUI in resSTPtr to hold minimum number
|       of errors in all windows
|     o avgWinDistF in resSTPtr to hold mean number of
|       errors in all widows
|     o edDistSL in resSTPtr to have edit distance
|   - Returns:
|     o edit distance for query and mapped reference
|     o negative value if could not find distance
|       - def_noSeq_edDist if one read is missing sequence
|       - def_noMap_edDist if one read is unmapped
\-------------------------------------------------------*/
signed long
dist_edDist(
   struct samEntry *qrySTPtr, /*read for edit distance*/
   struct seqST *refSTPtr,    /*has reference sequence*/
   unsigned int indelLenUI,   /*min indel length*/
   unsigned char minQUC,      /*min Q-score (snp/ins)*/
   unsigned int winSizeUI,    /*size of one window*/
   struct res_edDist *resSTPtr/*holds results*/
);

/*-------------------------------------------------------\
| Fun08: percDist_edDist
|   - gives a rough idea on precentage of difference from
|     error
|   - not great, but allows lumping reads together
| Input:
|   - resSTPtr:
|     o pointer to res_edDist struct with edit distance
|       and overlap size
|   - winSizeUI:
|     o size of one window
|   - percErrF:
|     o expected percent of errors (0 to 1) in reads
|     o if doing read to read comparsions; double this
| Output:
|   - Returns:
|     o probablity of edit distance due to chance
\-------------------------------------------------------*/
signed int
percDist_edDist(
   struct res_edDist *resSTPtr,
   unsigned int winSizeUI,
   float percErrF
);

/*-------------------------------------------------------\
| Fun09: addReadToDepth_edDist
|   - adds a read to a depth profile
| Input:
|   - refSTPtr:
|     o pointer to samEntry struct with reference
|   - qrySTPtr:
|     o query to add to profile
|   - minQUC:
|     o minimum q-score to keep snp
|   - minOverlapF:
|     o minimum percent overlap to score
|     o use 0 or negative value to not check overlap
|       - still checks to see reads share bases
|   - resSTPtr:
|     o pointer to res_edDist struct with depth array to
|       store read depths in (must be length of reference)
| Output:
|   - Modifies:
|     o depthAryUI in resSTPtr to have query added
|   - Returns:
|     o 0 for no errors
|     o 1 if read was not added (failed filters)
\-------------------------------------------------------*/
signed char
addReadToDepth_edDist(
   struct samEntry *refSTPtr, /*reference samEntry*/
   struct samEntry *qrySTPtr, /*query samEntry*/
   unsigned char minQUC,      /*min Q-score (snp/ins)*/
   float minOverlapF,         /*min % overlap*/
   struct res_edDist *resSTPtr/*has depth profile array*/
);

/*-------------------------------------------------------\
| Fun10: mkDepthProfile_edDist
|   - finds the depth of every base in samEntry reference
| Input:
|   - refSTPtr:
|     o pointer to samEntry struct with reference
|   - minQUC:
|     o minimum q-score to keep snp
|   - minOverlapF:
|     o minimum percent overlap to score
|   - resSTPtr:
|     o pointer to res_edDist struct to have depth profile
|   - samSTPtr:
|     o for reading each line in the sam file
|   - buffStrPtr:
|     o pointer to c-string with buffer for reading file
|   - lenBuffULPtr:
|     o pointer to unsigned long to hold buffStrPtr size
|   - samFILE:
|     o sam file to scan
| Output:
|   - Modifies:
|     o samSTPtr to have last entry in file
|     o samFILE to point to start
|     o buffStrPtr to have last line in sam file and is
|       resized if needed
|     o lenBuffULPtr to new buffStrPtr size (if resized)
|     o depthAryUI in resSTPtr to have depth profile
|   - Returns:
|     o 0 for no errors
|     o def_memErr_edDist for memory errors
|     o def_fileErr_edDist for memory errors
\-------------------------------------------------------*/
signed char
mkDepthProfile_edDist(
   struct samEntry *refSTPtr, /*reference samEntry*/
   unsigned char minQUC,      /*min Q-score (snp/ins)*/
   float minOverlapF,         /*min % overlap*/
   struct res_edDist *resSTPtr,/*has depth array*/
   struct samEntry *samSTPtr, /*for reading sam file*/
   signed char **buffStrPtr,  /*for reading sam file*/
   unsigned long *lenBuffULPtr, /*size of buffStrPtr*/
   void *samFILE
);

/*-------------------------------------------------------\
| Fun11: phead_edDist
|   - prints header for edDist output tsv
| Input:
|   - outFILE:
|     o FILE pointer to print header to
| Output:
|   - Prints:
|     o edDist header to outFILE
\-------------------------------------------------------*/
void
phead_edDist(
   void *outFILE
);

/*-------------------------------------------------------\
| Fun12: pdist_edDist
|   - prints edit distance tsv line for edDist
| Input:
|   - resSTPtr:
|     o pointer to res_edDist struct with distance to
|       print
|   - qryIdStr:
|     o c-string with query id (name) to print
|   - refIdStr:
|     o c-string with reference id (name) to print
|   - outFILE:
|     o FILE pointer to print header to
| Output:
|   - Prints:
|     o edDist header to outFILE
\-------------------------------------------------------*/
void
pdist_edDist(
   struct res_edDist *resSTPtr,
   signed char *qryIdStr,
   signed char *refIdStr,
   void *outFILE
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
