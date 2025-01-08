/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' edClust SOF: Start Of File
'   - clusters reads by similarity using edit distance
'   o header:
'     - included libraries
'   o fun01: getBestRead_edClust
'     - finds read with best score from a score array
'   o fun02: depthProf_edClust
'     - find the depth profile for a read
'   o fun03: findNumMap_edClust
'     - finds number of reads mapping to best read
'   o fun04: getDepth_edClust
'     - gets number of reads that overlap a cluster
'   o fun05: cluster_edClust
'     - clusters reads
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - forward declarations, defined variables, and guards
\-------------------------------------------------------*/

#ifndef EDIT_DISTANCE_CLUSTER_H
#define EDIT_DISTANCE_CLUSTER_H

struct samEntry;
struct set_tbCon;
struct res_edDist;

struct set_clustST ;
struct index_clustST;
struct con_clustST;

#define def_memErr_edClust -1
#define def_fileErr_edClust -2
#define def_noReads_edClust -4

/*-------------------------------------------------------\
| Fun01: getBestRead_edClust
|   - finds read with best score from a score array; fun01
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to have best read
|   - indexSTPtr:
|     o pointer to index_clustST with score array to
|       search
|   - buffStrPtr:
|     o c-string to use reading sam file
|     o this must be big enough to get the entire line
|       use buffer used with mk_index_clustST to get this
|   - samFILE:
|     o sam file to find best read in
| Output:
|   - Modifies:
|     o samFILE to be at start
|     o samSTPtr to have selected read
|   - Returns:
|     o index of selected read in scoreArySC
|     o def_memErr_edClust for memory errors
|     o def_fileErr_edClust for file errors
|     o def_noReads_edClust if no reads in file
\-------------------------------------------------------*/
signed long
getBestRead_edClust(
   struct samEntry *samSTPtr,
   struct index_clustST *indexSTPtr,
   signed char *buffStr,
   void *samFILE
);

/*-------------------------------------------------------\
| Fun02: depthProf_edClust
|   - find the depth profile for a read
| Input:
|   - refSTPtr:
|     o pointer to samEntry struct with reference (read)
|       to get profile for
|   - refNumUI:
|     o reference the reference read was mapped to
|   - clustSI:
|     o cluster working on
|   - indexSTPtr:
|     o pointer to index_clustST struct with file index's
|   - clustSetSTPtr:
|     o pointer to set_clustST struct with settings
|   - resEdDistSTPtr:
|     o pointer to res_edDist struct for holding edit
|       distance results
|   - samSTPtr:
|     o pointer to samEntry struct to use in file reading
|   - buffStr:
|     o pointer to c-string use as buffer in file reading
|     o needs to be able to hold largest line in file
|   - samFILE:
|     o FILE pointer to sam file with reads
| Output:
|   - Modifies:
|     o samSTPtr to have last read in file
|     o buffStr to have last line in file
|     o samFILE to point to the start
|     o errSCPtr to hold error messages:
|  - Returns:
|    o 0 for no errors
|    o def_memErr_edClust for memory error
|    o def_fileErr_edClust for file error
\-------------------------------------------------------*/
signed char
depthProf_edClust(
   struct samEntry *refSTPtr, /*read to get profile for*/
   unsigned int refNumUI,   /*reference read mapped to*/
   signed int clustSI,               /*cluster on*/
   struct index_clustST *indexSTPtr, /*has file index's*/
   struct set_clustST *clustSetSTPtr,/*settings*/
   struct res_edDist *resEdDistSTPtr,/*for edDist result*/
   struct samEntry *samSTPtr,        /*for reading file*/
   signed char *buffStr,             /*for reading file*/
   void *samFILE                     /*sam file*/
);

/*-------------------------------------------------------\
| Fun03: findNumMap_edClust
|   - finds number of reads mapping to best read
| Input:
|   - numReadsSL:
|     o pointer to unsigned long to have number of mapped
|       reads
|   - lenConUI:
|     o pointer to unsigned int to have consensus length
|   - bestSTPtr:
|     o pointer to samEntry struct with best read
|   - bestIndexSL:
|     o index of best read (return from fun01)
|   - clustSetSTPtr:
|     o has settings for filtering and stores the first
|       and last possible base of the consensesus
|   - indexSTPtr:
|     o pointer to index_clustST struct with reads index's
|       and open reads
|   - resEdDistSTPtr:
|     o pointer to res_edDist struct to hold results from
|       edDist comparision
|   - clustSI:
|     o cluster working on (what to assign to targArySI)
|   - conBl:
|     o 1: use consensus to read error rate
|     o 0: use read to read error rate
|   - buffStrPtr:
|     o c-string to use reading sam file
|     o should be able to hold longest line in file
|       (fun19 sets this up)
|   - samSTPtr:
|     o pointer to samEntry struct for reading sam file
|   - samFILE:
|     o sam file with reads to compare to best read
| Output:
|   - Modifies:
|     o numReadsUL to have number of reads that mapped or
|       the error
|       - def_memErr_edClust for memory errors
|       - def_fileErr_edClust for file errors
|     o lenConUI to have length of returend consensus
|     o samSTPtr to be last read in file
|     o samFILE to be at start
|     o buffStrPtr to have last sam file entry
|     o startUI and endUI in clustSetSTPtr to have start
|       and end of the consensus
|     o clustArySI in indexSTPtr to have hits assigned to
|       clustSI
|   - Returns:
|     o conNt_tbCon array with uncollapsed consensus
\-------------------------------------------------------*/
struct conNt_tbCon *
findNumMap_edClust(
   signed long *numReadsSL,         /*num mapped reads*/
   unsigned int *lenConUI,          /*consensus length*/
   struct samEntry *bestSTPtr,      /*best read for map*/
   signed long bestIndexSL,         /*best read index*/
   struct set_clustST *clustSetSTPtr,/*cluster settings*/
   struct set_tbCon *conSetSTPtr,   /*consensus settings*/
   struct index_clustST *indexSTPtr,/*has read index's*/
   struct res_edDist *resEdDistSTPtr, /*for edDist*/
   signed int clustSI,              /*cluster working on*/
   signed char conBl,             /*1: best is consensus*/
   struct samEntry *samSTPtr,     /*for reading sam file*/
   signed char *buffStr,          /*for reading sam file*/
   void *samFILE                  /*has reads*/
);

/*-------------------------------------------------------\
| Fun04: getDepth_edClust
|   - gets number of reads that overlap a cluster
| Input:
|   - conListSTPtr:
|     o pointer to head of con_clustST struct list with
|       cluster consensus to check depths for
|   - indexSTPtr:
|     o pointer to index_clustST struct with each reads
|       assigned cluster, start position, and end position
|   - minOverlapF:
|     o mininum overlap to count read and cluster as
|       overlapping (for total count)
| Output:
|   - Modifies:
|     o maxReadsUL in each con_clustST struct in
|       conListSTPtr to have number of possible overlaping
|       reads
\-------------------------------------------------------*/
void
getMaxDepth_edClust(
   struct con_clustST *conListSTPtr,
   struct index_clustST *indexSTPtr,
   float minOverlapF
);

/*-------------------------------------------------------\
| Fun05: cluster_edClust
|   - clusters reads
| Input:
|   - indexSTPtr:
|     o pointer to index_clustST pointer to set to
|       clusters (should be null if not; is freeded)
|   - clustSetSTPtr:
|     o pointer to set_clustST struct with cluster setting
|   - conSetSTPtr:
|     o pointer to set_tbCon struct with consensus setting
|   - samSTPtr:
|     o pointer to samEntry struct to use in reading file
|   - buffStrPtr:
|     o pointer to c-string to use reading sam files
|   - lenBuffULPtr:
|     o pointer to unsigned long with buffStrPtr length
|   - samFILE:
|     o sam file with reads to cluster
|   - logFILE:
|     o file to print log to (use 0 for no logging)
|   - errSCPtr:
|     o pointer to signed char to hold error type
| Output:
|   - Modifies:
|     o all arries clustSTPtr to have clusters and index's
|     o buffStrPtr to be resized to longest line in file
|     o lenBuffULPtr to hav buffStrPtr's length
|     o samFILE to point to start
|     o clustSI in clustSetSTPtr to be on last cluster
|     o errSC:
|       o 0 for no errors
|       o def_memErr_edClust for memory errors
|       o def_fileErr_edClust for file errors
|       o def_noReads_edClust if to few reads
|   - Prints:
|     o log (stats; results, and errors) to logFILE
|   - Returns:
|     o list of con_clustST structs with consensus
|       clusters
\-------------------------------------------------------*/
struct con_clustST *
cluster_edClust(
   struct index_clustST **indexSTPtr,
   struct set_clustST *clustSetSTPtr,
   struct set_tbCon *conSetSTPtr,
   struct samEntry *samSTPtr,
   signed char **buffStrPtr,
   unsigned long *lenBuffULPtr,
   void *samFILE,
   void *logFILE,
   signed char *errSCPtr
);

/*-------------------------------------------------------\
| Fun09: plist_con_clustST
|   - print the clusters consensuses in a con_clustST list
| Input:
|   - conSTPtr:
|     o pointer to con_clustST struct list with
|       consensuses to print
|   - headerStr:
|     o c-string with header to print 
|   - pgHeadStr:
|     o c-string with program header to print 
|   - buffStrPtr:
|     o to c-string to print consensuses with
|   - lenBuffULPtr:
|     o pointer to unsigned long with current length of
|       buffer
|   - outFILE:
|     o file to print consensuses to
| Output:
|   - Modifies:
|     o buffStrPtr to be resized if needed
|     o lenBuffULPtr to have current buffer size
|   - Prints:
|     o headers and consensus to outFILE as a sam file
|   - Returns:
|     o 0 for no errors
|     o def_memErr_edClust for memeory errors
\-------------------------------------------------------*/
signed char
plist_con_clustST(
   struct con_clustST *conSTPtr, /*consensuses to print*/
   signed char *headerStr,       /*sam file header*/
   signed char *pgHeadStr,       /*program header*/
   signed char **buffStrPtr,     /*for printing*/
   unsigned long *lenBuffULPtr,  /*size of buffStrPtr*/
   void *outFILE                 /*file to print to*/
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
