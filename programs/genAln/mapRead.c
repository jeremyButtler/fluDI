/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' mapRead SOF: Start Of File
'   - maps a read to a reference
'   o .h st01: chains_mapRead
'     - holds kmer chains build from kmerChain_mapRead
'   o fun01: blank_chains_mapRead
'     - blanks a chains_mapRead struct
'   o fun02: init_chains_mapRead
'     - initialize a chains_mapRead struct
'   o fun03: freeStack_chains_mapRead
'     - frees all variables in a chains_mapRead struct
'   o fun04: freeHeap_chains_mapRead
'     - frees a chains_mapRead struct
'   o fun05: setup_chains_mapRead
'     - sets up memory for a chains_mapRead struct
'   o fun06: mk_chains_mapRead
'     - makes a heap allocated chains_mapRead struct
'   o fun07: realloc_chains_mapRead
'     - reallocates memroy in a chains_mapRead struct
'   o .c fun08: uiFind_mapRead
'     - finds first signed int in a range
'   o .c fun09: findChain_chains_mapRead
'     - finds kmer in range from kmerChain_mapRead return
'   o fun10: kmerChain_mapRead
'     - merge kmers into chains (longer kmers)
'   o .c fun11: numQryChains_mapRead
'     - gets the number of chains assigned to a single
'       query kmer
'   o fun12: mergeChains_mapRead
'     - merges kmer chains from fun10 (kmerChain) into
'       the chain the maximizes base count
'   o fun13: mergeToSam_mapRead
'     - convert output from mergeChains_map read to sam
'       file entry (mapq is set to 0)
'   o fun14: flipRef_mapRead
'     - flips samEntry reference to query
'   o fun15: scoreSubAln_mapRead
'     - scores the aligment in samEntry and finds best
'       sub-aligment
'   o fun16: blank_set_mapRead
'     - blank a set_mapRead struct
'   o fun17: init_set_mapRead
'     - initializes a set_mapRead struct
'   o fun18: freeStack_set_mapRead
'     - frees variables in a set_mapRead struct
'   o fun19: freeHeap_set_mapRead
'     - frees a set_mapRead struct
'   o fun20: setup_set_mapRead
'     - allocates memory for structs in a set_mapRead
'   o fun21: addQryKmers_mapRead
'     - get forward and reverse kmer sequences
'   o fun22: aln_mapRead
'     - maps a read to a reference
'   o fun23: qckAln_mapRead
'     - maps a read to a reference; returns best alignment
'     - qckAln_mapRead differs from aln_mapRead in that
'       it uses kmerScan to find best window in entire
'       reference, and then after scan, makes alignment
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "mapRead.h"

#include "../genLib/numToStr.h"
#include "../genLib/ulCp.h"
#include "../genLib/shellSort.h"

#include "../genBio/seqST.h"
#include "../genBio/samEntry.h"

#include "../genAln/alnSet.h"
#include "../genAln/dirMatrix.h"
#include "../genAln/needle.h"
#include "../genAln/water.h"
#include "../genAln/kmerFind.h"

/*.h files only*/
#include "defsMapRead.h" /*default settings*/
#include "../genLib/genMath.h" /*.h max/min macros only*/
#include "../genBio/kmerBit.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries
!   - std #include <stdio.h>
!   - .c  #include "../genLib/base10Str.h"
!   - .c  #include "../genLib/charCp.h"
!   - .c  #include "../genLib/strAry.h"
!   - .c  #include "../genAln/indexToCoord.h"
!   - .c  #include "../genAln/memwater.h"
!   - .h  #include "../genBio/ntTo5Bit.h"
!   - .h  #include "../genAln/alnDefs.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: blank_chains_mapRead
|   - blanks a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct to blank
| Output:
|   - Modifies:
|     o refKmerAryUI and qryKmerAryUI in chainsSTPtr to
|       have all values set to -1
|     o chainLenArySI and scoreArySI in chainsSTPtr to
|       have all values set to 0
|     o indexArySI in chainsSTPtr to be filled with -1's
|     o lenChainSI and lenIndexSL in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
blank_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to blank*/
){
   signed long slElm = 0;

   if(! chainsSTPtr)
      return; /*nothing to blank*/

   if(chainsSTPtr->refArySI)
   { /*If: need to blank chain arrays*/

      slElm = 0;

      while(slElm < chainsSTPtr->lenChainSI)
      { /*Loop: blank chain arrays*/
         chainsSTPtr->refArySI[slElm] = -1;
         chainsSTPtr->qryArySI[slElm] = -1;
         chainsSTPtr->chainArySI[slElm] = 0;
         chainsSTPtr->scoreArySI[slElm] = 0;
         chainsSTPtr->nextArySI[slElm++] = 0;
      } /*Loop: blank chain arrays*/

   } /*If: need to blank chain arrays*/

   if(chainsSTPtr->indexArySI)
   { /*If: need to blank chain arrays*/

      slElm = 0;

      while(slElm < chainsSTPtr->lenIndexSI)
         chainsSTPtr->indexArySI[slElm++] = -1;

   } /*If: need to blank chain arrays*/

   chainsSTPtr->lenIndexSI = 0;
   chainsSTPtr->lenChainSI = 0;
   chainsSTPtr->lastIndexSI = 0;
} /*blank_chains_mapRead*/

/*-------------------------------------------------------\
| Fun02: init_chains_mapRead
|   - initialize a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct to initialize
| Output:
|   - Modifies:
|     o refKmerAryUI, qryKmerAryUI, chainLenArySI,
|       scoreArySI, and indexArySI to be 0 (null)
|     o lenChainSI, sizeChainSI, lenIndexSL, and
|       sizeIndexSI in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
init_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to initialize*/
){
   if(! chainsSTPtr)
      return; /*nothing to initialize*/

   chainsSTPtr->refArySI = 0;
   chainsSTPtr->qryArySI = 0;
   chainsSTPtr->chainArySI = 0;
   chainsSTPtr->scoreArySI = 0;
   chainsSTPtr->nextArySI = 0;
   chainsSTPtr->indexArySI = 0;

   chainsSTPtr->sizeChainSI = 0;
   chainsSTPtr->sizeIndexSI = 0;

   blank_chains_mapRead(chainsSTPtr);
} /*init_chains_mapReads*/

/*-------------------------------------------------------\
| Fun03: freeStack_chains_mapRead
|   - frees all variables in a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead with variables to free
| Output:
|   - Modifies:
|     o frees refKmerAryUI, qryKmerAryUI, chainLenArySI,
|       scoreArySI, and indexArySI and sets to 0 (null)
|     o lenChainSI, sizeChainSI, lenIndexSL, and
|       sizeIndexSI in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
freeStack_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to free vars*/
){
   if(! chainsSTPtr)
      return;

   if(chainsSTPtr->refArySI)
      free(chainsSTPtr->refArySI);

   if(chainsSTPtr->qryArySI)
      free(chainsSTPtr->qryArySI);

   if(chainsSTPtr->chainArySI)
      free(chainsSTPtr->chainArySI);

   if(chainsSTPtr->scoreArySI)
      free(chainsSTPtr->scoreArySI);

   if(chainsSTPtr->nextArySI)
      free(chainsSTPtr->nextArySI);

   if(chainsSTPtr->indexArySI)
      free(chainsSTPtr->indexArySI);

   init_chains_mapRead(chainsSTPtr);
} /*freeStack_chains_mapRead*/

/*-------------------------------------------------------\
| Fun04: freeHeap_chains_mapRead
|   - frees a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to free
| Output:
|   - Modifies:
|     o frees chainsSTPtr (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to free*/
){
   if(! chainsSTPtr)
      return;

   freeStack_chains_mapRead(chainsSTPtr);
   free(chainsSTPtr);
} /*freeHeap_chains_mapRead*/

/*-------------------------------------------------------\
| Fun05: setup_chains_mapRead
|   - sets up memory for a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to setup memory
|   - lenChainsSL:
|     o length of chains arrays
|   - lenRefUI:
|     o length of reference
| Output:
|   - Modifies:
|     o allocates lenChainsSL items of memory to
|       refCkmerArySI, qryKmerArySI, chainLenArySI,
|       and scoreArySI in chainsSTPtr
|     o allocates lenRefUI memory to indexArySI in
|       chainsSTPtr
|     o sizeChainSI in chainsSTPtr to be lenChainSI 
|     o sizeIndexSI in chainsSTPtr to be lenIndexSL 
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*to setup*/
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
){
   signed long slElm = 0;

   if(! chainsSTPtr)
      return 0;

   /*making sure no old memory hanging around*/
   freeStack_chains_mapRead(chainsSTPtr);

   chainsSTPtr->refArySI =
      malloc(lenChainsSL * sizeof(signed int));
   if(! chainsSTPtr->refArySI)
      goto memErr_fun05;

   chainsSTPtr->qryArySI =
      malloc(lenChainsSL * sizeof(signed int));
   if(! chainsSTPtr->qryArySI)
      goto memErr_fun05;

   chainsSTPtr->chainArySI =
      calloc(
         lenChainsSL,
         sizeof(signed int)
      );
   if(! chainsSTPtr->chainArySI)
      goto memErr_fun05;

   chainsSTPtr->scoreArySI =
      calloc(
         lenChainsSL,
         sizeof(signed int)
      );
   if(! chainsSTPtr->scoreArySI)
      goto memErr_fun05;

   chainsSTPtr->nextArySI =
      calloc(
         lenChainsSL,
         sizeof(signed int)
      );
   if(! chainsSTPtr->nextArySI)
      goto memErr_fun05;


   chainsSTPtr->indexArySI =
      calloc(
         lenRefUI,
         sizeof(signed int)
      );
   if(! chainsSTPtr->indexArySI)
      goto memErr_fun05;


   chainsSTPtr->sizeChainSI = lenChainsSL;
   chainsSTPtr->sizeIndexSI = (signed int) lenRefUI;


   slElm = 0;
   while(slElm < lenChainsSL)
   { /*Loop: blank kmer starting coords arrays*/
      chainsSTPtr->refArySI[slElm] = -1;
      chainsSTPtr->qryArySI[slElm++] = -1;
   } /*Loop: blank kmer starting coords arrays*/

   slElm = 0;
   while(slElm < lenRefUI)
      chainsSTPtr->indexArySI[slElm++] = -1;


   return 0;

   memErr_fun05:;
      return def_memErr_mapRead;
} /*setup_chains_mapRead*/

/*-------------------------------------------------------\
| Fun06: mk_chains_mapRead
|   - makes a heap allocated chains_mapRead struct
| Input:
|   - lenChainsSL:
|     o length of chains arrays
|   - lenRefUI:
|     o length of reference
| Output:
|   - Returns:
|     o pointer to initialized and setup chains_mapRead
|       struct
|     o 0 for memory errors
\-------------------------------------------------------*/
struct chains_mapRead *
mk_chains_mapRead(
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
){
   struct chains_mapRead *retSTPtr = 0;

   retSTPtr = malloc(sizeof(struct chains_mapRead));
   if(! retSTPtr)
      goto memErr_fun06;

   init_chains_mapRead(retSTPtr);

   if(
      setup_chains_mapRead(
         retSTPtr,
         lenChainsSL,
         (signed long) lenRefUI
      )
   ) goto memErr_fun06;


   goto ret_fun06;

   memErr_fun06:;
      if(retSTPtr)
         freeHeap_chains_mapRead(retSTPtr);
      retSTPtr = 0;
      goto ret_fun06;
   ret_fun06:;
      return retSTPtr;
} /*mk_chains_mapRead*/

/*-------------------------------------------------------\
| Fun07: realloc_chains_mapRead
|   - reallocates memroy in a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to rellocate memory
|   - lenChainsSL:
|     o new length of chains arrays (0 for do not resize)
|   - lenRefUI:
|     o new length of reference (0 for do not resize)
| Output:
|   - Modifies:
|     o refKmerArySI, qryKmerArySI, chainLenArySI,
|       scoreArySI in chainsSTPtr to be reallocated to
|       lenChainSI if lenChainSI >chainsSTPtr->sizeChainSI
|     o indexArySI chainsSTPtr to be reallocated to
|       lenRefUI if lenRefUI  > chainsSTPtr->sizeIndexSI
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
realloc_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*reallocate*/
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC:
   '   - reallocates memroy in a chains_mapRead struct
   '   o fun07 sec01:
   '     - variable declaration + input check
   '   o fun07 sec02:
   '     - chain arrays reallocs + initialize
   '   o fun07 sec03:
   '     - reference index array rellocs
   '   o fun07 sec04:
   '     - return error or no error
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01:
   ^   - variable declaration + input check
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *scPtr = 0;
   signed long slElm = 0;

   if(! chainsSTPtr)
      return 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^   - chain arrays reallocs + initialize
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(lenChainsSL > chainsSTPtr->sizeChainSI)
   { /*If: reallocating chains array memory*/
      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->refArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->refArySI = (signed int *) scPtr;

      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->qryArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->qryArySI = (signed int *) scPtr;

      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->chainArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->chainArySI = (signed int *) scPtr;

      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->scoreArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->scoreArySI = (signed int *) scPtr;


      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->nextArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->nextArySI = (signed int *) scPtr;


      slElm = chainsSTPtr->sizeChainSI;
      while(slElm < lenChainsSL)
      { /*Loop: blank kmer starting coords arrays*/
         chainsSTPtr->refArySI[slElm] = -1;
         chainsSTPtr->qryArySI[slElm] = -1;
         chainsSTPtr->chainArySI[slElm] = 0;
         chainsSTPtr->scoreArySI[slElm] = 0;
         chainsSTPtr->nextArySI[slElm++] = 0;
      } /*Loop: blank kmer starting coords arrays*/

      chainsSTPtr->sizeChainSI = lenChainsSL;
   } /*If: reallocating chains array memory*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^   - reference index array rellocs
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(lenRefUI > (unsigned int) chainsSTPtr->sizeIndexSI)
   { /*If: reallocating chains array memory*/
      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->indexArySI,
            lenRefUI * sizeof(signed long)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->indexArySI = (signed int *) scPtr;

      slElm = chainsSTPtr->sizeIndexSI;
      while(slElm < lenRefUI)
         chainsSTPtr->indexArySI[slElm++] = 0;

      chainsSTPtr->sizeIndexSI = lenRefUI;
   } /*If: reallocating chains array memory*/


   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec04:
   ^   - return error or no error
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun07_sec04:;
      return def_memErr_mapRead;
} /*realloc_chains_mapRead*/

/*-------------------------------------------------------\
| Fun08: uiFind_mapRead
|   - finds first unsigned int in a range
| Input:
|   - siAry:
|     o array to search for query in
|   - qrySI:
|     o query to find
|   - lenSI:
|     o length of array (index 1)
| Output:
|   - Returns:
|     o index of first item in array
|     o -1 if item not in array
\-------------------------------------------------------*/
signed long
uiFind_mapRead(
   unsigned int *uiAry,   /*array to search*/
   unsigned int qryUI,    /*query to search for*/
   signed int lenSI       /*array length (index 1)*/
){
   signed long midSL = 0;
   signed long rightSL = lenSI - 1;
   signed long leftSL = 0;

   if(lenSI <= 0)
      return -1;

   while(leftSL <= rightSL)
   { /*Loop: find query*/
      midSL = (leftSL + rightSL) >> 1;

      if(qryUI > uiAry[midSL])
         leftSL = midSL + 1;
      else if(qryUI < uiAry[midSL])
         rightSL = midSL - 1;
      else
      { /*Else: is match*/
         if(! midSL)
            return midSL; /*at start of array*/

         if(uiAry[midSL - 1] == qryUI)
            rightSL = midSL - 1; /*still more matches*/
         else
            return midSL;
      } /*Else: is match*/
   } /*Loop: find query*/

   return -1; /*no matches found*/
} /*uiFind_mapRead*/

/*--------------------------------------------------------\
| Fun09: findChain_chains_mapRead
|   - finds chain with input kmer coordiantes in a
|     chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to chains_mapRead struct to search
|   - refSI:
|     o reference starting position of kmer
|   - qrySI:
|     o query starting position of kmer
| Output:
|   - Returns:
|     o index of at least one kmer chain with both
|       refSI and qrySI kmer coordaintes
|     o -1 if did not find kmer
| Note:
|   - data must be sorted by reference start kmer position,
|     then by query start kmer positions
\--------------------------------------------------------*/
signed long
findChain_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*has chains*/
   signed int refSI,        /*reference kmer start*/
   signed int qrySI         /*query kmer start*/
){
   signed long midSL = 0;
   signed long rightSL = chainsSTPtr->lenChainSI - 1;
   signed long leftSL = 0;
   signed long refEndSL = 0;
   signed long qryEndSL = 0;

   if(chainsSTPtr->lenChainSI <= 0)
      return -1;

   while(leftSL <= rightSL)
   { /*Loop: find kmers*/
      midSL = (leftSL + rightSL) >> 1;

      refEndSL = chainsSTPtr->refArySI[midSL];
      refEndSL += chainsSTPtr->chainArySI[midSL];

      qryEndSL = chainsSTPtr->qryArySI[midSL];
      qryEndSL += chainsSTPtr->chainArySI[midSL];

      if(refSI > refEndSL)
         leftSL = midSL + 1;
      else if(qrySI > qryEndSL)
         leftSL = midSL + 1;
      else if(refSI < chainsSTPtr->refArySI[midSL])
         rightSL = midSL - 1;
      else if(qrySI < chainsSTPtr->qryArySI[midSL])
         rightSL = midSL - 1;
      else
         return midSL;
         /*chain includes reference and query kmers*/
   } /*Loop: find kmers*/

   return -1; /*nothing found*/
} /*findChain_chains_mapRead*/

/*-------------------------------------------------------\
| Fun10: kmerChain_mapRead
|   - merge kmers into chains (longer kmers)
| Input:
|   - chainsSTPtr:
|     o pointer to chainsSTPtr struct to add kmer chains
|       to
|   - minChainLenSI:
|     o minimum number of kmers to keep a chain
|   - qryKmerArySI:
|     o signed int array with query kmers (by position)
|   - lenQrySI:
|     o number kmers in qryKmerArySI
|   - refKmerArySI:
|     o signed int array with reference kmers to chain
|   - refSortArySI:
|     o signed int array that is same length as
|       refKmerArySI
|     o modified to have sorted kmers in refKmerArySI
|   - lenRefSI:
|     o number kmers in refKmerArySI
|   - refStartSI:
|     o first refrence base in window
| Output:
|   - Modifies:
|     o refKmerArySI, qryKmerArySI, and chainLenArySI in
|       chainsSTPtr to have chain starting coordinates and
|       length
|     o scoreArySI and nextArySI in chainsSTPtr to have
|       all detected kmer positions filled with 0's
|     o indexArySI in chainsSTPtr to have index or -1 for
|       each query kmer start position in qryKmerArySI
|     o lenChainSI and lenIndexSI to have length of
|       indexArySI (num ref kmers) and number of chains
|     o sizeChainSI and sizeIndexSI in chainsSTPtr if
|       needed to rellocate memory for chain/index arrays
|     o refSortArySI to have kmers in refKmerArySI
|       - kmers sorted, with negatives at end (unsigned)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|   - Note:
|     o sorting order of chains is query coordinate and
|       then reference coordinate
\-------------------------------------------------------*/
signed char
kmerChain_mapRead(
   struct chains_mapRead *chainsSTPtr,/*gets kmer chains*/
   signed int minChainLenSI,          /*min chain length*/
   signed int *qryKmerArySI,          /*query kmers*/
   signed int lenQrySI,               /*query length*/
   signed int *refKmerArySI,          /*reference kmers*/
   signed int *refSortArySI,          /*to sort kmers*/
   signed int lenRefSI,               /*number ref kmers*/
   signed int refStartSI              /*first ref base*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun10 TOC:
   '   - merge kmers to get longer kmers
   '   o fun10 sec01:
   '     - variable declarations
   '   o fun10 sec02:
   '     - initial memory allocation
   '   o fun10 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;          /*for returing errors*/

   signed long indexSL = 0;      /*index of kmer matches*/
   signed int qryIndexSI = 0;      /*qry kmer on*/
   signed int refIndexSI = 0;      /*table kmer on*/
   signed int *refPosHeapArySI = 0;

   signed int lenChainSI = 0;      /*length of 1 chain*/
   signed int siKmer = 0;          /*reference kmer iter*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec02:
   ^   - initial memory allocation
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_chains_mapRead(chainsSTPtr);

   if(chainsSTPtr->sizeIndexSI < lenQrySI)
   { /*If: need to resize index*/
      errSC =
         realloc_chains_mapRead(
            chainsSTPtr,
            0, /*not resizing chains array*/
            lenQrySI
         );
      if(errSC)
         goto memErr_fun10_sec04;
   } /*If: need to resize index*/

   chainsSTPtr->lenIndexSI = lenQrySI;


   /*allows me to quickly get index's for chain extend
   `   step
   */
   refPosHeapArySI= malloc(lenRefSI * sizeof(signed int));
   if(! refPosHeapArySI)
      goto memErr_fun10_sec04;
   for(
      indexSL = 0;
      indexSL < lenRefSI;
      ++indexSL
   ) refPosHeapArySI[indexSL] = (signed int) indexSL;


   cpLen_ulCp(
      (signed char *) refSortArySI,
      (signed char *) refKmerArySI,
      lenRefSI * sizeof(signed int)
   ); /*copy kmers*/

   uiTwinSort_kmerFind(
      (unsigned int *) refSortArySI, /*kmers to sort*/
      (unsigned int *) refPosHeapArySI,
      0,                       /*start at first kmer*/
      lenRefSI - 1             /*last kmer (index 0)*/
   ); /*O[r * log(r)] for quick look up*/
      /* by doing an unsigned int sort, I put the negative
      `  values to the end of the array
      */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec03:
   ^   - find kmer chains
   ^   o fun10 sec03 sub01:
   ^     - find if have kmer chain(s) + start detect loop
   ^   o fun10 sec03 sub02:
   ^     - build all kmer chains for one ref kmer
   ^   o fun10 sec03 sub03:
   ^     - find next possible kmer match to ref kmer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun10 Sec03 Sub01:
   *   - find if have kmer chain(s) + start detect loop
   \*****************************************************/

   /* very hard to estimate O/theta (finding worst case
   `   is hard); this is an off the`cuff guess
   `  O[  q * log(r) qry-find * O[r * log(hits) * O[r|q
   `    + log(r) next-ref
   `   ]
   `   = O[q * r * log(2r + hits) * (r|q)]
   `      but r|q worst case is many duplicate kmers
   `      that can never meet the minimum chain
   `      length
   `   o no chain case, never extend chain beyond
   `     minimum chain length (mChain)
   `   = O[q * r * log(2r + hits) * mChain]
   `   o however, this assumes that for each q, we
   `     are iterating over the entireity of r, but
   `     this would require duplicate kmers, which
   `     would mean the minimum chain length,
   `     so there must be enough kmers that non-mathcing
   `     patterns form.
   `     * for 1 kmer only, it is dulpicate, everything
   `       is merged, so [r * log(hits) * mChain] becomes
   `       [1 * log(1) + mChain]
   `     * for 2 kmers, query is divided by 2, inorder to
   `       be worst case the reference must have a
   `       different kmer at every mChain interval. This
   `       is highly unlikey in biological data. The max
   `       length before duplication is (2 ^ mChain)
   `     * howver, in thuis case we have no hits, so
   `       log(2r + hits) becomes log(2q)
   `   = O[q * (r / 2) * log(2r + hits) * mChain]
   `   o the worst case is highly unlikely, also, it is
   `     very likely kmer diversity will be high, so this
   `     more likely scales closer, but not quite to
   `     kmer-space [r / min(length(r), 4 ^ kmer-length]
   `   = theta[q * kmer-space * log(2r + hits) * mChain]
   `   o for most cases the kmer space will be larger
   `     (7mer to 11mer; 16000 kmers to 4000000 kmers)
   `     than the read length, so only duplications
   `     will have an impact on the theta case, so
   `     log(2r + hits) might scale to log(r + hits)
   `   o the kmer scan step should reduce the chance
   `     of duplications, by limiting windows size
   `   o based on this it may often, very roughly scale to
   `   = theta[q * log(r + hits] * mChain]
   `     * probably my overally optimisitic guess
   */
   for(
      siKmer = 0;
      siKmer < lenQrySI;
      ++siKmer
   ){ /*Loop: find matching kmers*/

      chainsSTPtr->indexArySI[siKmer] = -1;

      if(qryKmerArySI[siKmer] < 0)
         continue; /*anonymous base*/

      indexSL =
         uiFind_mapRead(
            (unsigned int *) refSortArySI,
            (unsigned int) qryKmerArySI[siKmer],
            lenRefSI
         ); /*find first query kmer for ref*/
            /*O[log(r)] ref-find*/
            /*unsigned int forces so negative values
            `  to be at end
            */

      if(indexSL < 0)
         continue; /*no match*/

      /**************************************************\
      * Fun10 Sec03 Sub02:
      *   - build all kmer chains for one ref kmer
      \**************************************************/

      /*  O[r] ref-loop
      ` * O[log(hits)] ref-find
      ` * O[r | q] extend
      */
      while(qryKmerArySI[siKmer] == refSortArySI[indexSL])
      { /*Loop: make kmer chains*/

         qryIndexSI = siKmer;
         refIndexSI = refPosHeapArySI[indexSL];
         lenChainSI = 0;

         if(
            findChain_chains_mapRead(
               chainsSTPtr,
               refIndexSI + refStartSI,
               qryIndexSI
            ) >= 0
         ) goto nextMatch_fun10_sec03_sub03;
           /*already have this kmer chain*/
           /*O[log(hits)] hit-find*/


         while(
                qryKmerArySI[qryIndexSI++]
             == refKmerArySI[refIndexSI++]
         ) ++lenChainSI; /*find kmer chain length*/
           /*O[r | q] extend*/


         if(lenChainSI >= minChainLenSI)
         { /*If: keeping kmer chain*/

             if(
                   chainsSTPtr->lenChainSI
                >= chainsSTPtr->sizeChainSI
             ){ /*If: need more memory*/
                errSC =
                   realloc_chains_mapRead(
                      chainsSTPtr,
                      chainsSTPtr->lenChainSI << 1,
                      0 /*not rellocating index arrray*/
                   );
                if(errSC)
                   goto memErr_fun10_sec04;
             } /*If: need more memory*/

            if(chainsSTPtr->indexArySI[siKmer] < 0)
            { /*If: first chain for this query base*/
               chainsSTPtr->indexArySI[siKmer] =
                  chainsSTPtr->lenChainSI;

               chainsSTPtr->lastIndexSI = siKmer;
            } /*If: first chain for this query base*/

             /*reference kmer starting point*/
             chainsSTPtr->refArySI[
                chainsSTPtr->lenChainSI
             ] =  refPosHeapArySI[indexSL] + refStartSI;

             /*query starting point*/
             chainsSTPtr->qryArySI[
                chainsSTPtr->lenChainSI
             ] = siKmer;

             /*kmer chain length*/
             chainsSTPtr->chainArySI[
                chainsSTPtr->lenChainSI++
             ] = lenChainSI;
               
         } /*If: keeping kmer chain*/

         /***********************************************\
         * Fun10 Sec03 Sub03:
         *   - find next possible kmer match to ref kmer
         \***********************************************/

         nextMatch_fun10_sec03_sub03:;
            ++indexSL; /*get off last match*/
      } /*Loop: make kmer chains*/
   }  /*Loop: find matching kmers*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec04:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun10_sec04;

   memErr_fun10_sec04:;
      errSC = def_memErr_mapRead;
      goto ret_fun10_sec04;

   ret_fun10_sec04:
      if(refPosHeapArySI)
         free(refPosHeapArySI);
      refPosHeapArySI = 0;
      return errSC;
} /*kmerChain_mapRead*/

/*-------------------------------------------------------\
| Fun11: numQryChains_mapRead
|   - gets the number of chains assigned to a single
|     query kmer
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct with reference
|       chaines
|   - index:
|     o index of reference kmer (kmer number) to find
|       number of chains for
| Output:
|   - Returns:
|     o number of chains for reference kmer
\-------------------------------------------------------*/
#define numQryChains_mapRead(chainsSTPtr, index) ( ((chainsSTPtr)->indexArySI[(indexSL)] < 0) ? -1 : ( (index) < (chainsSTPtr)->lenIndexSI ) ? ( (chainsSTPtr)->indexArySI[(index) + 1] - (chainsSTPtr)->indexArySI[(index)] ) : ( (chainsSTPtr)->lenChainSI - (chainsSTPtr)->indexArySI[(index)] ) )
/*Logic:
`   - checks to see if reference kmer has any assigned
`     assigned chains
`   - turnary checks if have more index's (ref chains)
`     after the current index
`     (index < chainsSTPtr->lenIndexSI)
`   - if have more index's; length is found by
`     subtracing next reference starting index
`     (indexArySI[index + 1]) from the current index
`     (indexArySI[index])
`   - if at end; length is found by
`     subtracing the array length (chainsSTPtr->lenChainSI)
`     minus the references starting index
`     (chainsSTPtr->indexArySI[index])
*/

/*-------------------------------------------------------\
| Fun12: mergeChains_mapRead
|   - merges kmer chains from fun10 (kmerChain) into
|     the chain the maximizes base count
|   - not optimal solution, but also reduce complexity
|     by not counting for gaps
| Input:
|   - chainAryUI:
|     o unsigned int array with kmer chains to merge into
|       a "best" chain
| Output:
|   - Modifies:
|     o nextArySI in chains_mapRead to have best scoring
|       trail for each chain
|     o scoreArySI to have most bases possible in a chain
|   - Returns:
|     o index of best merged chain
|     o def_noChains_mapRead if no chains in array
\-------------------------------------------------------*/ 
signed int
mergeChains_mapRead(
   struct chains_mapRead *chainsSTPtr  /*chains to merge*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun12 TOC:
   '   - merges kmer chains from fun10 (kmerChain)
   '   o fun12 sec01:
   '     - variable declarations
   '   o fun12 Sec02:
   '     - memory allocation
   '   o fun12 sec03:
   '     - check few chain cases
   '   o fun12 sec04:
   '     - merge chains
   '   o fun12 sec05:
   '     - return results
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int siChain = 0;    /*main chain on*/
   signed int siSub = 0;      /*subchain on*/
   signed int lastKmerSI = 0; /*index of last ref kmer*/
   signed int lastIndexSI = 0;/*index of last refkmer*/
   signed int subIndexSI = 0; /*index of subchain ref*/

   signed int scoreSI = 0;    /*score of merge*/
   signed int refUniqSI = 0;  /*number uniq bases in ref*/
   signed int refEndSI = 0;   /*ref ending coordinate*/

   signed int qryUniqSI = 0; /*number unique query bases*/
   signed int qryEndSI = 0;  /*query ending coordinate*/

   signed int maxNtSI = 0;    /*most bases in chain*/
   signed int curMaxNtSI = 0; /*current best base count*/
   signed int maxIndexSI = 0; /*index of chain maximized*/
   signed int curMaxIndexSI = 0;

   signed int secMaxNtSI = 0; /*2nd most bases in chain*/
   signed int curSecMaxNtSI = 0; /*current rounds score*/
      /*used as a quick check to reduces chances of a full
      `  scan. The logic is if
      `     (nt in chain) + secMaxNtSI
      `   < maxNtSI + (nt in chain) - overlap
      `  then I can safely say that I can not maximize
      `  the score any more for overlap cases
      */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec02:
   ^   - memory allocation
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec03:
   ^   - check few chain cases
   ^   o fun12 sec03 sub01:
   ^     - no chains and only one chain case
   ^   o fun12 sec03 sub02:
   ^     - only two chain case
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun12 Sec03 Sub01:
   *   - no chains and only one chain case
   \*****************************************************/

   if(chainsSTPtr->lenChainSI <= 0)
      goto noChains_fun12_sec05;

   if(chainsSTPtr->lenChainSI == 1)
   { /*If: only one chain*/
      maxIndexSI = 0;

      chainsSTPtr->scoreArySI[0] =
         chainsSTPtr->chainArySI[0];

      chainsSTPtr->nextArySI[0] = 0; /*nothing after*/

      goto ret_fun12_sec05;
   } /*If: only one chain*/

   /*****************************************************\
   * Fun12 Sec03 Sub02:
   *   - two chain case
   *   o fun12 sec03 sub02 cat01:
   *     - two chains, setup + find ending coordinats
   *   o fun12 sec03 sub02 Cat02:
   *     - two chains, handle non-merge cases
   *   o fun12 sec03 sub02 Cat03:
   *     - two chains, handle merge cases
   *   o fun12 sec03 sub02 Cat04:
   *     - goto return section
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun12 Sec03 Sub02 Cat01:
   +   - two chains, setup + find ending coordinates
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   else if(chainsSTPtr->lenChainSI == 2)
   { /*Else If: might be quick merge*/

      chainsSTPtr->scoreArySI[0] =
         chainsSTPtr->chainArySI[0];
      chainsSTPtr->scoreArySI[1] =
         chainsSTPtr->chainArySI[1];
      maxNtSI = 0;

      /*find end of chains*/
      refUniqSI =
           chainsSTPtr->refArySI[0]
         + chainsSTPtr->chainArySI[0];
         /*ending coordinate of first reference chain*/
      refUniqSI =
           chainsSTPtr->refArySI[1]
         + chainsSTPtr->chainArySI[1];
         /*ending coordinate of second reference chain*/

      qryUniqSI =
           chainsSTPtr->qryArySI[0]
         + chainsSTPtr->chainArySI[0];
         /*ending coordiante of frist query chain*/
      qryEndSI =
           chainsSTPtr->qryArySI[1]
         + chainsSTPtr->chainArySI[1];
         /*ending coordiante of second query chain*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun12 Sec03 Sub02 Cat02:
      +   - two chains, handle non-merge cases
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*the chians are always sorted by reference and then
      `  query kmer starting position
      ` first check is checking if one chain enguffled the
      `   other or if second chain is before first
      */
      if(
            qryUniqSI >= qryEndSI
         || refUniqSI >= refEndSI
         ||
               chainsSTPtr->refArySI[0]
            >= chainsSTPtr->refArySI[1]
         ||   
               chainsSTPtr->qryArySI[0]
            >= chainsSTPtr->qryArySI[1]
      ){ /*If: chains can not be merged*/
         if(
              chainsSTPtr->chainArySI[0]
            > chainsSTPtr->chainArySI[1]
         ) maxIndexSI = 0; /*first chain is best*/

         else
            maxIndexSI = 1; /*second chain is best*/
      }  /*If: chains can not be merged*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun12 Sec03 Sub02 Cat03:
      +   - two chains, handle merge cases
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else
      { /*Else: can merge chains*/
         maxIndexSI = 0;

         /*find if have gap between chains*/
         refUniqSI = chainsSTPtr->refArySI[1] - refUniqSI;
         qryUniqSI = chainsSTPtr->qryArySI[1] - qryUniqSI;

         /*check if have overlap (gap)*/
         if(
               refUniqSI < 0
            || qryUniqSI < 0
         ){ /*If: have overlap*/
            maxNtSI =
               min_genMath(
                  refUniqSI,
                  qryUniqSI
               ); /*find most overlap (most negative)*/

            maxNtSI += chainsSTPtr->chainArySI[0];
            /*find length minus overlap*/
         }  /*If: have overlap*/

         else
            maxNtSI = chainsSTPtr->chainArySI[0];

         /*get chain length*/
         maxNtSI += chainsSTPtr->chainArySI[1];
         chainsSTPtr->scoreArySI[0] = maxNtSI;
         chainsSTPtr->nextArySI[0] = 1;
      } /*Else: can merge chains*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun12 Sec03 Sub02 Cat04:
      +   - goto return section
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      goto ret_fun12_sec05;
   } /*Else If: might be quick merge*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec04:
   ^   - merge chains
   ^   o fun12 sec04 sub01:
   ^     - find last mapped reference kmer
   ^   o fun12 sec04 sub02:
   ^     - merge kmers
   ^   o fun12 sec04 sub03:
   ^     - score last reference kmer chains
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun12 Sec04 Sub01:
   *   - find last mapped reference kmer
   \*****************************************************/

   lastKmerSI = chainsSTPtr->lastIndexSI;

   if(lastKmerSI < 0)
      goto noChains_fun12_sec05;
   else
   { /*Else: have chains to merge*/
      lastIndexSI = lastKmerSI;
      lastKmerSI = chainsSTPtr->indexArySI[lastKmerSI];
   } /*Else: have chains to merge*/

   /*****************************************************\
   * Fun12 Sec04 Sub02:
   *   - score last reference kmer chains
   \*****************************************************/

   siChain = chainsSTPtr->lenChainSI - 1;

   while(
         chainsSTPtr->qryArySI[siChain]
      >= chainsSTPtr->qryArySI[lastKmerSI]
   ){ /*Loop: merge kmer chains*/
      chainsSTPtr->scoreArySI[siChain] =
         chainsSTPtr->chainArySI[siChain];
      chainsSTPtr->nextArySI[siChain] = 0;/*end of merge*/

      /*want best score on first reference kmer, that
      `  way I avoid long gaps (hopefully)
      */
      if(chainsSTPtr->chainArySI[siChain] >= maxNtSI)
      { /*If: new best score*/
         secMaxNtSI = maxNtSI;
         maxNtSI = chainsSTPtr->chainArySI[siChain];
         maxIndexSI = siChain;   
      } /*If: new best score*/

      else if(
         chainsSTPtr->chainArySI[siChain] > secMaxNtSI
      ) secMaxNtSI = chainsSTPtr->chainArySI[siChain];
         /*so always have second best case*/

      --siChain;

      if(siChain < 0)
         goto ret_fun12_sec05;
         /*no chains can be merged*/
   } /*Loop: merge kmer chains*/

   /*****************************************************\
   * Fun12 Sec04 Sub03:
   *   - score chains
   *   o fun12 sec04 sub03 cat01:
   *     - loop though chains & find end coords of chain
   *   o fun12 sec04 sub03 cat02:
   *     - handel no overlap cases (use highest score)
   *   o fun12 sec04 sub03 cat03:
   *     - handel overlap, but still best score cases
   *   o fun12 sec04 sub03 cat04:
   *     - handle cases were need to find best score
   *   o fun12 sec04 sub03 cat05:
   *     - merge valid subchains
   *   o fun12 sec04 sub03 cat06:
   *     - move to next sub chain
   *   o fun12 sec04 sub03 cat07:
   *     - move back one chain
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun12 Sec04 Sub03 Cat01:
   +   - loop though chains & find end coords of chain
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   curMaxNtSI = maxNtSI;
   curSecMaxNtSI = secMaxNtSI;

   /*this is kinda like the pheramon trail for the ant
   `   traveling salesman problem.
   */
   while(siChain >= 0)
   { /*Loop: merge chains*/
      qryEndSI =
           chainsSTPtr->qryArySI[siChain]
         + chainsSTPtr->chainArySI[siChain];

      refEndSI =
           chainsSTPtr->qryArySI[siChain]
         + chainsSTPtr->chainArySI[siChain];

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun12 Sec04 Sub03 Cat02:
      +   - handel no overlap cases (use highest score)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      qryUniqSI =
         qryEndSI - chainsSTPtr->qryArySI[maxIndexSI];
      refUniqSI =
         refEndSI - chainsSTPtr->refArySI[maxIndexSI];

      if(
            qryUniqSI < 0
         && refUniqSI < 0
      ){ /*If: no overlap*/
         chainsSTPtr->scoreArySI[siChain] =
            maxNtSI + chainsSTPtr->chainArySI[siChain];
         chainsSTPtr->nextArySI[siChain] = maxIndexSI;

         if(
               chainsSTPtr->scoreArySI[siChain]
            >= curMaxNtSI
         ){ /*If: new best score*/
            curSecMaxNtSI = curMaxNtSI;
            curMaxNtSI = chainsSTPtr->scoreArySI[siChain];
            curMaxIndexSI = siChain;
         } /*If: new best score*/

         else if(
              chainsSTPtr->scoreArySI[siChain]
            > curSecMaxNtSI
         ) curSecMaxNtSI=chainsSTPtr->scoreArySI[siChain];

         goto mvBack_fun12_sec04_sub03_cat07;
      }  /*If: no overlap*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun12 Sec04 Sub03 Cat03:
      +   - handel overlap, but still best score cases
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
              chainsSTPtr->chainArySI[siChain]
           +  min_genMath(refUniqSI, qryUniqSI)
           +  maxNtSI
         > secMaxNtSI
      ){ /*Else If: with overlap still is best score*/
         chainsSTPtr->scoreArySI[siChain] =
               chainsSTPtr->chainArySI[siChain]
            +  min_genMath(refUniqSI, qryUniqSI)
            + maxNtSI; /*find score*/

         chainsSTPtr->nextArySI[siChain] = maxIndexSI;

         if(
              chainsSTPtr->scoreArySI[siChain]
            >= curMaxNtSI
         ){ /*If: new best score*/
            curSecMaxNtSI = curMaxNtSI;
            curMaxNtSI = chainsSTPtr->scoreArySI[siChain];
            curMaxIndexSI = siChain;
         } /*If: new best score*/

         else if(
              chainsSTPtr->scoreArySI[siChain]
            > curSecMaxNtSI
         ) curSecMaxNtSI=chainsSTPtr->scoreArySI[siChain];

         goto mvBack_fun12_sec04_sub03_cat07;
      }  /*Else If: with overlap still is best score*/
         
      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun12 Sec04 Sub03 Cat04:
      +   - handle cases were need to find best score
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*this is the n^2 part, however, I am reallying
      `  on kmer space to be scarce enough that this step
      `  will be rare
      */

      siSub = lastKmerSI; /*last chain was on*/
      subIndexSI = lastIndexSI;/*qry index of last chain*/

      while(siSub < chainsSTPtr->lenChainSI)
      { /*Loop: scan next query kmer*/

         /*skipping chains I can not merge*/
         if(
                refEndSI
            >=  
                chainsSTPtr->refArySI[siSub]
              + chainsSTPtr->chainArySI[siSub]
         ) ;
           /*reference is before/engulfs chain*/

         else if(
              qryEndSI
            >=
                chainsSTPtr->qryArySI[siSub]
              + chainsSTPtr->chainArySI[siSub]
         ){ /*Else If: complete query overlap*/
             ++subIndexSI;
             siSub = chainsSTPtr->indexArySI[subIndexSI];
         }  /*Else If: complete query overlap*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun12 Sec04 Sub03 Cat05:
         +   - merge valid subchains
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         else
         { /*Else: have gap or overlap (ref and query)*/
            qryUniqSI =
               chainsSTPtr->qryArySI[siSub] - qryEndSI;
            refUniqSI =
               chainsSTPtr->refArySI[siSub] - refEndSI;

            if(
                  refUniqSI < 0
               || qryUniqSI < 0
            ){ /*If: have overlap*/
               scoreSI = 
                  min_genMath( 
                     refUniqSI,
                     qryUniqSI
                  ); /*compinsate for overlap*/
            }  /*If: have overlap*/

            scoreSI = chainsSTPtr->chainArySI[siSub];
            scoreSI += chainsSTPtr->chainArySI[siChain];

            if(scoreSI > chainsSTPtr->scoreArySI[siChain])
            { /*If: new local best score*/
               chainsSTPtr->scoreArySI[siChain] = scoreSI;
               chainsSTPtr->nextArySI[siChain] = siSub;

               if(scoreSI > curMaxNtSI)
               { /*If: new best score*/
                  curSecMaxNtSI = maxNtSI;
                  curMaxNtSI = scoreSI;
                  curMaxIndexSI = siChain;
               } /*If: new best score*/

               else if(scoreSI > curSecMaxNtSI)
                  curSecMaxNtSI = scoreSI;
            } /*If: new local best score*/

            else if(scoreSI > curSecMaxNtSI)
               curSecMaxNtSI = scoreSI;
         } /*Else: have gap or overlap (ref and query)*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun12 Sec04 Sub03 Cat06:
         +   - move to next sub chain
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         ++siSub; /*move to next hit*/

      } /*Loop: scan next query kmer*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun12 Sec04 Sub03 Cat07:
      +   - move back one chain
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      mvBack_fun12_sec04_sub03_cat07:;
         --siChain;

         if(! lastKmerSI)
            ; /*If: on last chain*/
         else if(
             chainsSTPtr->qryArySI[siChain]
           < chainsSTPtr->qryArySI[lastKmerSI - 1]
         ){ /*If: on next reference kmer*/
             /*update to add in new maximum values for the
             `  next reference kmer (round)
             */
             maxNtSI = curMaxNtSI;
             maxIndexSI = curMaxIndexSI;
             secMaxNtSI = curSecMaxNtSI;

            --lastIndexSI;

            lastKmerSI =
               chainsSTPtr->indexArySI[lastIndexSI];

            while(lastKmerSI < 0)
            { /*Loop: find previous kmer*/
               if(lastIndexSI < 0)
                  goto ret_fun12_sec05; /*finished*/

               --lastIndexSI;

               lastKmerSI =
                  chainsSTPtr->indexArySI[lastIndexSI];
            } /*Loop: find previous kmer*/
         } /*If: on next reference kmer*/
   } /*Loop: merge chains*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec05:
   ^   - return results
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun12_sec05;

   noChains_fun12_sec05:;
      maxIndexSI = def_noChains_mapRead;
      goto ret_fun12_sec05;

   ret_fun12_sec05:;
      return maxIndexSI;
} /*mergeChains_mapRead*/

/*-------------------------------------------------------\
| Fun13: mergeToSam_mapRead
|   - convert output from mergeChains_map read to sam
|     file entry (mapq is set to 0)
| Input:
|   - chainsSTPtr:
|     o chains_mapRead struct pointer with chains to use
|       in making the alignment
|   - indexSI:
|     o starting index of merged chain to align
|   - refSTPtr:
|     o tbl_kmerFind struct pointer with reference
|       sequence
|   - qrySTPtr:
|     o ref_kmerFind struct used with kmerChain_mapRead
|       (is query)
|   - revBl:
|     o query direction (same as kmerChain_mapRead)
|       * 1 is reverse complement
|       * 9 is foward
|   - chainLenSI:
|     o minimum length of a chain
|     o used to align ends before chains, alignment size
|       is 1.5 chains (for chainLenSI = 20, is 30 bases)
|       * 1.5 is off the top of my head
|   - settings:
|     o settings for alignment
|   - samSTPtr:
|     o samEntry struct pointer to add alignment to
| Output:
|   - Modifies:
|     o samSTPtr to have alignment; everything in samSTPtr
|       is modified
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
mergeToSam_mapRead(
   struct chains_mapRead *chainsSTPtr,/*chains to use*/
   signed int indexSI,            /*chain index to start*/
   struct tblST_kmerFind *refSTPtr,/*reference sequence*/
   struct refST_kmerFind *qrySTPtr,/*has query sequence*/
   signed char revBl,            /*1 is reverse*/
   signed int chainLenSI,      /*min length of one chain*/
   struct alnSet *settings,      /*aligment settings*/
   struct samEntry *samSTPtr     /*gets final alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun13 TOC:
   '   - convert output form mergeChains_map read to sam
   '     file entry (mapq is set to 0)
   '   o fun13 sec01:
   '     - variable declarations
   '   o fun13 sec02:
   '     - initialize and setup
   '   o fun13 sec03:
   '     - align ending soft masked bases
   '   o fun13 sec04:
   '     - align ending soft masked bases
   '   o fun13 sec05:
   '     - finish up samEntry values (no mapq)
   '   o fun13 sec06:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;
   struct seqST *qrySTDoNotFree = 0;
   signed int lastChainBl = 0;
      /*tells if have done last chain yet*/

   /*coordiantes for aligning chunks*/
   signed int refStartSI = 0;
   signed int refEndSI = 0;
   signed int refLastEndSI = 0;

   signed int qryStartSI = 0;
   signed int qryEndSI = 0;
   signed int qryLastEndSI = 0;

   /*for alignments*/
   unsigned int maskUI = 0; /*throw away*/
   unsigned int anonUI = 0; /*throw away*/
   unsigned int tmpUI = 0;  /*cigar lengths (throw away)*/
   struct dirMatrix matrixStackST;

   /*for sam entry editing*/
   signed char *swapSCPtr = 0; /*for cigar reallocs*/
   signed int siCig = 0;
   signed int newCigPosSI = 0;
   signed long scoreSL = 0; /*for end checking*/

   /*for dealing with soft masked bases at start*/
   signed int *cigHeapArySI = 0;
   signed char *cigHeapStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec02:
   ^   - initialize and setup
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   chainLenSI += (chainLenSI >> 1); /*softmasks at end*/

   cigHeapArySI =
      malloc((chainLenSI +1) * sizeof(signed int));
   if(! cigHeapArySI)
      goto memErr_fun13_sec06;

   cigHeapStr =
      malloc((chainLenSI + 1) * sizeof(signed char));
   if(! cigHeapStr)
      goto memErr_fun13_sec06;


   /*for merged chain gap filing*/
   init_dirMatrix(&matrixStackST);
   blank_samEntry(samSTPtr);


   /*get direction of alignment for query*/
   if(revBl)
   { /*If: reverse complement mapping*/
      qrySTDoNotFree = qrySTPtr->revSeqST;
      samSTPtr->flagUS = 16;
   } /*If: reverse complement mapping*/

   else
   { /*If: foward mapping*/
      qrySTDoNotFree = qrySTPtr->forSeqST;
      samSTPtr->flagUS = 0;
   } /*If: foward mapping*/


   /*setup reference start coords*/
   samSTPtr->refStartUI =
      (unsigned int) chainsSTPtr->refArySI[indexSI];

   if(chainsSTPtr->qryArySI[indexSI] > 0)
   { /*If: have starting softmasking*/
      samSTPtr->cigTypeStr[0] = 'S';

      samSTPtr->cigArySI[0] =
         chainsSTPtr->qryArySI[indexSI];
         /*coordinates in index 0, but cigar index 1*/

      siCig = 1;
      samSTPtr->lenCigUI = 1;
   } /*If: have starting softmasking*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec03:
   ^   - align gaps in merged chains
   ^   o fun13 sec03 sub01:
   ^     - expand cigar array memory + start loop
   ^   o fun13 sec03 sub02:
   ^     - find first and last base in current chain
   ^   o fun13 sec03 sub03:
   ^     - handel easy insertion/deletion cases
   ^   o fun13 sec03 sub04:
   ^     - handel easy (single) snp cases
   ^   o fun13 sec03 sub05:
   ^     - handel gaps that need to be aligned
   ^   o fun13 sec03 sub06:
   ^     - add chain length to matches
   ^   o fun13 sec03 sub07:
   ^     - move to next chain in merge
   ^   o fun13 sec03 sub08:
   ^     - check if have last chain to merge
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun13 Sec03 Sub01:
   *   - expand cigar array memory (if need) + start loop
   \*****************************************************/

   samSTPtr->cigTypeStr[siCig] = '=';
   samSTPtr->cigArySI[siCig] =
      chainsSTPtr->chainArySI[indexSI];
   samSTPtr->cigArySI[siCig] += (refSTPtr->lenKmerUC - 2);
        /*last kmer is always two bases off*/
   ++siCig;

   qryLastEndSI = chainsSTPtr->qryArySI[indexSI];
   qryLastEndSI += chainsSTPtr->chainArySI[indexSI];
   qryLastEndSI += (refSTPtr->lenKmerUC - 2);

   refLastEndSI = chainsSTPtr->refArySI[indexSI];
   refLastEndSI += chainsSTPtr->chainArySI[indexSI];
   refLastEndSI += (refSTPtr->lenKmerUC - 2);

   goto nextChain_fun13_sec03_sub07;
   /*first index will often be 0*/

   while(chainsSTPtr->nextArySI[indexSI] > 0)
   { /*Loop: while have chains to merge*/
      mergeChains_fun13_sec03_sub07:;

      /*doing + 3 because max for this loop is two
      `   entries, the thrid is for a softmask ending
      */
      if(
        (siCig +3) >= (signed int) samSTPtr->lenCigBuffUI)
      { /*If: need more cigar memory*/
         newCigPosSI= (signed int) samSTPtr->lenCigBuffUI;
         newCigPosSI += (newCigPosSI >> 1);

         swapSCPtr =
            realloc(
               samSTPtr->cigTypeStr,
               (newCigPosSI + 1) * sizeof(signed char)
            );
         if(! swapSCPtr)
            goto memErr_fun13_sec06;
         samSTPtr->cigTypeStr = swapSCPtr;

         swapSCPtr =
            (signed char *)
            realloc(
               samSTPtr->cigArySI,
               (newCigPosSI + 1) * sizeof(signed int)
            );
         if(! swapSCPtr)
            goto memErr_fun13_sec06;
         samSTPtr->cigArySI = (signed int *) swapSCPtr;

         samSTPtr->lenCigBuffUI =
            (unsigned int) newCigPosSI;
      } /*If: need more cigar memory*/

      /**************************************************\
      * Fun13 Sec03 Sub02:
      *   - find first and last base in current chain
      \**************************************************/

      /*find start and end of kmer chain*/
      qryStartSI = chainsSTPtr->qryArySI[indexSI];
      qryEndSI = qryStartSI;
      qryEndSI += chainsSTPtr->chainArySI[indexSI];
      qryEndSI += (refSTPtr->lenKmerUC - 2);
        /*last kmer is always one base off*/

      refStartSI = chainsSTPtr->refArySI[indexSI];
      refEndSI = refStartSI;
      refEndSI += chainsSTPtr->chainArySI[indexSI];
      refEndSI += (refSTPtr->lenKmerUC - 2);
        /*last kmer is always one base off*/

      /*account for overlaps in kmer chains*/
      if(qryStartSI < qryLastEndSI)
         qryStartSI = qryLastEndSI;

      if(refStartSI < refLastEndSI)
         refStartSI = refLastEndSI;

      /**************************************************\
      * Fun13 Sec03 Sub03:
      *   - handel easy insertion/deletion cases
      \**************************************************/

      /*check if have easy call*/
      if((qryStartSI - qryLastEndSI) < 2)
      { /*If: no query gap (means reference gap)*/
         samSTPtr->cigTypeStr[siCig] = 'D';
         samSTPtr->cigArySI[siCig] = 1;
         ++siCig;
      } /*If: no query gap (means reference gap)*/

      else if((refStartSI - refLastEndSI) < 2)
      { /*If: no reference gap (means query gap)*/
         samSTPtr->cigTypeStr[siCig] = 'I';
         samSTPtr->cigArySI[siCig] = 1;
         ++siCig;
      } /*If: no reference gap (means query gap)*/

      /**************************************************\
      * Fun13 Sec03 Sub04:
      *   - handel easy (single) snp cases
      \**************************************************/

      else if(
            qryStartSI - qryLastEndSI == 2
         && refStartSI - refLastEndSI == 2
      ){ /*Else If: mismatch*/
         samSTPtr->cigTypeStr[siCig] = 'X';
         samSTPtr->cigArySI[siCig] = 1;
         ++siCig;
      }  /*Else If: mismatch*/

      /**************************************************\
      * Fun13 Sec03 Sub05:
      *   - handel gaps that need to be aligned
      \**************************************************/

      else
      { /*Else: need to align*/
         ++qryLastEndSI;  /*move into gap between chains*/
         qrySTDoNotFree->offsetUL = qryLastEndSI;
         qrySTDoNotFree->endAlnUL =
            (unsigned long) qryStartSI - 1;

         ++refLastEndSI;  /*move into gap between chains*/
         refSTPtr->seqSTPtr->offsetUL = refLastEndSI;
         refSTPtr->seqSTPtr->endAlnUL= 
            (unsigned long) refStartSI -1;

         if(
            needle(
               qrySTDoNotFree,
               refSTPtr->seqSTPtr,
               &matrixStackST,
               settings
            ) == def_memErr_needle
         ) goto memErr_fun13_sec06;

         newCigPosSI =
            getCig_dirMatrix(
              &matrixStackST,        /*has alignment*/
              0,                     /*use matrix index*/
              0,                     /*always forward*/
              qrySTDoNotFree,        /*query sequence*/
              refSTPtr->seqSTPtr,    /*reference sequence*/
              &samSTPtr->cigTypeStr, /*cigar type array*/
              &samSTPtr->cigArySI,  /*has cigar lengths*/
              samSTPtr->lenCigUI,   /*position in cigar*/
              &samSTPtr->lenCigUI,  /*length of cigar*/
              (unsigned int *) &refLastEndSI,/*aln start*/
              &anonUI,               /*ignoring*/
              &maskUI,               /*ignoring*/
              1,                     /*keep ending dels*/
              settings               /*settings for aln*/
            ); /*add aligment to cigar*/

         siCig = newCigPosSI + 1;

         if( newCigPosSI < 0)
            goto memErr_fun13_sec06;

      } /*Else: need to align*/

      /**************************************************\
      * Fun13 Sec03 Sub06:
      *   - add chain length to matches
      \**************************************************/

      samSTPtr->refEndUI = (unsigned int) refEndSI;
         /*mark end of alignment*/
      samSTPtr->readLenUI = (unsigned int) qryEndSI;
         /*so know last query base*/

      /*setup for next round*/
      qryLastEndSI = qryEndSI;
      refLastEndSI = refEndSI;

      /*add chain lengths in*/
      qryEndSI -= qryStartSI; /*find query length*/
      refEndSI -= refStartSI; /*find reference length*/
         /*already accounted for extra bases in last
         `  kmer when found qryEndSI and refEndSI
         */

      samSTPtr->cigArySI[siCig] =
         min_genMath(
            qryEndSI,
            refEndSI
         );
      ++samSTPtr->cigArySI[siCig];/*convert to index 1*/

      samSTPtr->cigTypeStr[siCig++] = '=';
      samSTPtr->lenCigUI = (unsigned int) siCig;

      /**************************************************\
      * Fun13 Sec03 Sub07:
      *   - move to next chain in merge
      \**************************************************/

      nextChain_fun13_sec03_sub07:;
         if(chainsSTPtr->nextArySI[indexSI])
            indexSI = chainsSTPtr->nextArySI[indexSI];
            /*move to next chain in merge*/
   } /*Loop: while have chains to merge*/

   /*****************************************************\
   * Fun13 Sec03 Sub08:
   *   - check if have last chain to merge
   \*****************************************************/

   if(! lastChainBl)
   { /*If: still need to finish the last chain*/
      lastChainBl = 1;
      goto mergeChains_fun13_sec03_sub07;
   } /*If: still need to finish the last chain*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec04:
   ^   - align ending soft masked bases
   ^   o fun13 sec04 sub01:
   ^     - deal with starting softmasking
   ^   o fun13 sec04 sub02:
   ^     - deal with ending softmasking
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun13 Sec04 Sub01:
   *   - deal with starting softmasking
   *   o fun13 sec04 sub01 cat01:
   *     - find softamsking coordinates to align
   *   o fun13 sec04 sub01 cat02:
   *     - do waterman alignment8
   *   o fun13 sec04 sub01 cat03:
   *     - get cigar and check if keeping
   *   o fun13 sec04 sub01 cat04:
   *     - check if keeping cigar entries
   *   o fun13 sec04 sub01 cat05:
   *     - realloc main memory (if more needed)
   *   o fun13 sec04 sub01 cat06:
   *     - add new cigar to samEntry
   *   o fun13 sec04 sub01 cat07:
   *     - adjust softmasking at start
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun13 Sec04 Sub01 Cat01:
   +   - find softamsking coordinates to align
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(samSTPtr->cigTypeStr[0] == 'S')
   { /*If: have softmasking at start*/
      /*adding to end so can do swap*/

      if(samSTPtr->refStartUI < (unsigned int) chainLenSI)
         refStartSI = 0;
      else
         refStartSI =
            (signed int)
            samSTPtr->refStartUI - chainLenSI;

      if(samSTPtr->cigArySI[0] < chainLenSI)
         qryStartSI = 0;
      else
         qryStartSI =
            (signed int)
            samSTPtr->cigArySI[0] - chainLenSI;

      if(samSTPtr->refStartUI > 0)
         refEndSI = samSTPtr->refStartUI - 1;
      else
         refEndSI = 0;

      qryEndSI = samSTPtr->cigArySI[0];


      /*check if both sequences have softmasking*/
      if(refEndSI - refStartSI <= 0)
         goto skipStartMask_fun13_sec04_sub02;
      else if (qryEndSI - qryStartSI <= 0)
         goto skipStartMask_fun13_sec04_sub02;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub01 Cat02:
      +   - do waterman alignment
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      qrySTDoNotFree->offsetUL =
         (unsigned long) qryStartSI;
      qrySTDoNotFree->endAlnUL =
         (unsigned long) qryEndSI;

      refSTPtr->seqSTPtr->offsetUL =
         (unsigned long) refStartSI;
      refSTPtr->seqSTPtr->endAlnUL =
         (unsigned long) refEndSI;

      if(
         needle(
            qrySTDoNotFree,
            refSTPtr->seqSTPtr,
            &matrixStackST,
            settings
         ) == def_memErr_water
      ) goto memErr_fun13_sec06;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub01 Cat03:
      +   - get cigar and check if keeping
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      refLastEndSI = refStartSI;
      tmpUI = chainLenSI;
      cigHeapStr[0] = '\0';

      newCigPosSI =
         getCig_dirMatrix(
           &matrixStackST,        /*has alignment*/
           0,                     /*use matrix index*/
           0,                     /*always forward*/
           qrySTDoNotFree,        /*query sequence*/
           refSTPtr->seqSTPtr,    /*reference sequence*/
           &cigHeapStr,           /*cigar type array*/
           &cigHeapArySI,         /*has cigar lengths*/
           0,                     /*first cigar pos (0)*/
           &tmpUI,                /*cigar buffer size*/
           (unsigned int *) &refLastEndSI,/*start of aln*/
           &anonUI,               /*ignoring*/
           &maskUI,               /*ignoring*/
           1,                     /*keep ending dels*/
           settings               /*settings for aln*/
         ); /*add aligment to cigar*/
      
      if(newCigPosSI < 0)
         goto memErr_fun13_sec06;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub01 Cat04:
      +   - check if keeping cigar entries
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      ++newCigPosSI; /*convert to index 1*/

      maskUI = 0; /*storing cigar size in mask*/
      scoreSL = 0;
      siCig = 0;
      refEndSI = refStartSI;
      qryEndSI = qryStartSI;

      while(siCig < newCigPosSI)
      { /*Loop: score alignment*/
         ++maskUI;

         if(
               cigHeapStr[siCig] == 'D'
            || cigHeapStr[siCig] == 'I'
         ){ /*If: deletion*/

            if(cigTypeStr[siCig] == 'D')
               refEndSI += cigHeapArySI[siCig];
            else
               qryEndSI += cigHeapArySI[siCig];


            scoreSL -= settingsST->gapSS;
            scoreSL -=
                 settingsST->extendSS
               * (cigHeapArySI[siCig] - 1);

            if(scoreSL <= 0)
            { /*If: end of alignment*/
               maskUI = 0;
               qryStartSI = qryEndSI;
               refStartSI = refEndSI;
               scoreSL = 0;
            } /*If: end of alignment*/

         } /*If: deletion*/

         else if(
               cigTypeStr[siCig] != '='
            || cigTypeStr[siCig] != 'X'
         ){ /*If: match/snp*/
            tmpUI = (signed int) cigHeapArySI;

            while(tmpUI > 0)
            { /*Loop: find score*/
               scoreSL +=
                  qryDoNotFree->seqStr[qryEndSI++],
                  refSTPtr->seqSTPtr->seqStr[refEndSI++],
                  settingsST
               );

               --tmpUI;
            } /*Loop: find score*/

            if(scoreSL <= 0)
            { /*If: end of alignment*/
               maskUI = 0;
               qryStartSI = qryEndSI;
               refStartSI = refEndSI;
               scoreSL = 0;
            } /*If: end of alignment*/
         } /*If: match/snp*/
      } /*Loop: score alignment*/

      if(! scoreSL)
         goto skipStartMask_fun13_sec04_sub02;
         /*If: no good subaligment by end*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub01 Cat05:
      +   - realloc main memory (if more needed)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      indexSI = (signed int) samSTPtr->lenCigUI + maskUI;

      if(indexSI >= (signed int) samSTPtr->lenCigBuffUI)
      { /*If: need more memory*/

         swapSCPtr =
            realloc(
               samSTPtr->cigTypeStr,
               (indexSI + 1) * sizeof(signed char)
            );
         if(! swapSCPtr)
            goto memErr_fun13_sec06;
         samSTPtr->cigTypeStr = swapSCPtr;

         swapSCPtr =
            (signed char *)
            realloc(
               samSTPtr->cigArySI,
               (indexSI + 1) * sizeof(signed int)
            );
         if(! swapSCPtr)
            goto memErr_fun13_sec06;
         samSTPtr->cigArySI = (signed int *) swapSCPtr;

         samSTPtr->lenCigBuffUI = (unsigned int) indexSI;
         
      } /*If: need more memory*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub01 Cat06:
      +   - add new cigar to samEntry
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      indexSI = (signed int) samSTPtr->lenCigUI - 1;
      indexSI += (newCigPosSI - (cigHeapStr[0] == 'S'));
      /*removes softmasking at start if needed*/

      siCig = samSTPtr->lenCigUI;
      samSTPtr->lenCigUI = indexSI + 1;

      while(siCig > 0)
      { /*Loop: move cigar up*/
         samSTPtr->cigTypeStr[indexSI] ^=
            samSTPtr->cigTypeStr[siCig];
         samSTPtr->cigTypeStr[siCig] ^=
            samSTPtr->cigTypeStr[indexSI];
         samSTPtr->cigTypeStr[indexSI] ^=
            samSTPtr->cigTypeStr[siCig];

         samSTPtr->cigArySI[indexSI] ^=
            samSTPtr->cigArySI[siCig];
         samSTPtr->cigArySI[siCig] ^=
            samSTPtr->cigArySI[indexSI];
         samSTPtr->cigArySI[indexSI] ^=
            samSTPtr->cigArySI[siCig];

         --siCig;
         --indexSI;
      } /*Loop: move cigar up*/

      siCig = 1;
      indexSI = cigHeapStr[0] != 'S';
      /*1 if no softmasking; 0 if softmasking*/
 
      while(indexSI < newCigPosSI)
      { /*Loop: add new cigar entries in*/
         samSTPtr->cigTypeStr[siCig]= cigHeapStr[indexSI];
         samSTPtr->cigArySI[siCig]=cigHeapArySI[indexSI];

         if(cigHeapStr[indexSI] != 'I')
            samSTPtr->refStartUI -= cigHeapArySI[indexSI];
            /*adjust reference starting position*/

         ++siCig;
         ++indexSI;
      } /*Loop: add new cigar entries in*/


      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub01 Cat07:
      +   - adjust softmasking at start
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*removed aligned bases from softmasked entry*/
      if(samSTPtr->cigArySI[0] >= chainLenSI)
         samSTPtr->cigArySI[0] -= chainLenSI;
      else
         samSTPtr->cigArySI[0] = 0;

      if(cigHeapStr[0] == 'S')
         samSTPtr->cigArySI[0] += cigHeapStr[0];
         /*added aligned softmasking to cigar*/

      if(samSTPtr->cigArySI[0] <= 0)
      { /*If: no softmasking at start*/

         /*bit inefficent, but works*/

         siCig = 0;
         while(siCig < (signed int) samSTPtr->lenCigUI)
         { /*Loop: move cigar entries back one*/
            samSTPtr->cigTypeStr[siCig] ^=
               samSTPtr->cigTypeStr[siCig + 1];
            samSTPtr->cigTypeStr[siCig + 1] ^=
               samSTPtr->cigTypeStr[siCig];
            samSTPtr->cigTypeStr[siCig] ^=
               samSTPtr->cigTypeStr[siCig + 1];

            samSTPtr->cigArySI[siCig] ^=
               samSTPtr->cigArySI[siCig + 1];
            samSTPtr->cigArySI[siCig + 1] ^=
               samSTPtr->cigArySI[siCig];
            samSTPtr->cigArySI[siCig] ^=
               samSTPtr->cigArySI[siCig + 1];

            ++siCig;
         } /*Loop: move cigar entries back one*/

         --samSTPtr->lenCigUI;
      } /*If: no softmasking at start*/
   } /*If: have softmasking at start*/

   /*****************************************************\
   * Fun13 Sec04 Sub02:
   *   - deal with ending softmasking
   *   o fun13 sec04 sub02 cat01:
   *     - find softamsking coordinates to align
   *   o fun13 sec04 sub02 cat02:
   *     - do waterman alignment8
   *   o fun13 sec04 sub02 cat03:
   *     - get cigar and check if keeping
   *   o fun13 sec04 sub02 cat04:
   *     - realloc main memory (if more needed)
   *   o fun13 sec04 sub02 cat05:
   *     - add new cigar to samEntry
   *   o fun13 sec04 sub02 cat06:
   *     - adjust softmasking at start
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun13 Sec04 Sub02 Cat01:
   +   - find softamsking coordinates to align
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   skipStartMask_fun13_sec04_sub02:;

   siCig = samSTPtr->lenCigUI - 1;

   if(samSTPtr->cigTypeStr[siCig] == 'S')
   { /*If: have ending softmasking*/
      /*adding to end so can do swap*/

      refStartSI = (signed int) samSTPtr->refEndUI;

      refEndSI =
         (signed int)
         refSTPtr->seqSTPtr->lenSeqUL - refStartSI;

      qryStartSI = (signed int) qrySTDoNotFree->lenSeqUL;
      qryEndSI =
         qryStartSI - (signed int) samSTPtr->readLenUI;


      if(
            refEndSI < 1
         || qryEndSI < 1
      ) goto skipEndMask_fun13_sec05;
        /*at least on position has no softmasking*/


      if(refEndSI < chainLenSI)
         refEndSI =
            (signed int) refSTPtr->seqSTPtr->lenSeqUL - 1; 
      else
         refEndSI = refStartSI + chainLenSI - 1;
            

      if(qryEndSI < chainLenSI)
         qryEndSI =
           (signed int) qrySTDoNotFree->lenSeqUL - 1;
      else
         qryEndSI = qryStartSI + chainLenSI - 1;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub02 Cat02:
      +   - do waterman alignment
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      qrySTDoNotFree->offsetUL=(unsigned long) qryStartSI;
      qrySTDoNotFree->endAlnUL = (unsigned long) qryEndSI;

      refSTPtr->seqSTPtr->offsetUL =
         (unsigned long) refStartSI;
      refSTPtr->seqSTPtr->endAlnUL= (unsigned long) refEndSI;

      if(
         water(
            qrySTDoNotFree,
            refSTPtr->seqSTPtr,
            &matrixStackST,
            settings
         ) == def_memErr_water
      ) goto memErr_fun13_sec06;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub02 Cat03:
      +   - get cigar and check if keeping
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      refLastEndSI = refStartSI;
      tmpUI = chainLenSI;
      cigHeapStr[0] = '\0';

      newCigPosSI =
         getCig_dirMatrix(
            &matrixStackST,      /*has alignment*/
            0,                     /*use matrix index*/
            0,                   /*always forward*/
            qrySTDoNotFree,      /*query sequence*/
            refSTPtr->seqSTPtr,     /*reference sequence*/
            &cigHeapStr,         /*cigar type array*/
            &cigHeapArySI,       /*has cigar lengths*/
            0,                     /*first cigar pos (0)*/
            &tmpUI,                /*cigar buffer size*/
            (unsigned int *) &refStartSI,/*start of aln*/
            &anonUI,             /*ignoring*/
            &maskUI,             /*ignoring*/
            1,                   /*keep ending dels*/
            settings             /*settings for aln*/
         ); /*add aligment to cigar*/
      
      if(newCigPosSI < 0)
         goto memErr_fun13_sec06;


      if(cigHeapStr[0] == 'S')
         goto skipEndMask_fun13_sec05;
         /*softmasking between alignment end and aligned*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub02 Cat04:
      +   - realloc main memory (if more needed)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      indexSI =
         (signed int) samSTPtr->lenCigUI + newCigPosSI;

      if(indexSI >= (signed int) samSTPtr->lenCigBuffUI)
      { /*If: need more memory*/

         swapSCPtr =
            realloc(
               samSTPtr->cigTypeStr,
               (indexSI + 2) * sizeof(signed char)
            );
         if(! swapSCPtr)
            goto memErr_fun13_sec06;
         samSTPtr->cigTypeStr = swapSCPtr;

         swapSCPtr =
            (signed char *)
            realloc(
               samSTPtr->cigArySI,
               (indexSI + 2) * sizeof(signed int)
            );
         if(! swapSCPtr)
            goto memErr_fun13_sec06;
         samSTPtr->cigArySI = (signed int *) swapSCPtr;

         samSTPtr->lenCigBuffUI = (unsigned int) indexSI;
         
      } /*If: need more memory*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub02 Cat05:
      +   - add new cigar to samEntry
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      siCig = samSTPtr->lenCigUI;
      samSTPtr->lenCigUI = indexSI;
      indexSI = 0;

      while(siCig < (signed int) samSTPtr->lenCigUI)
      { /*Loop: move cigar up*/
         if(cigHeapStr[indexSI] == 'S')
            break;

         if(cigHeapStr[siCig] != 'I')
            samSTPtr->refEndUI += cigHeapArySI[siCig];
            /*reference length*/

         if(cigHeapStr[siCig] != 'D')
            samSTPtr->readLenUI += cigHeapArySI[siCig];
            /*query length*/

         samSTPtr->cigTypeStr[siCig] =
            cigHeapStr[indexSI];

         samSTPtr->cigArySI[siCig] =
            cigHeapArySI[indexSI];

         ++siCig;
         ++indexSI;
      } /*Loop: move cigar up*/

      samSTPtr->cigTypeStr[siCig] = '\0';
      samSTPtr->lenCigUI = siCig;
   } /*If: have ending softmasking*/
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec05:
   ^   - finish up samEntry values (no mapq)
   ^   o fun13 sec05 sub01:
   ^     - add ending softmasking in
   ^   o fun13 sec05 sub02:
   ^     - add read lengths
   ^   o fun13 sec05 sub03:
   ^     - copy read ids
   ^   o fun13 sec05 sub04:
   ^     - copy sequence
   ^   o fun13 sec05 sub05:
   ^     - copy q-score
   ^   o fun13 sec05 sub06:
   ^     - find mutation counts
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun13 Sec05 Sub01:
   *   - add ending softmasking in
   \*****************************************************/

   skipEndMask_fun13_sec05:;

   /*add ending softmasking*/
   if(samSTPtr->readLenUI < qrySTDoNotFree->lenSeqUL)
   { /*If: have softmasking at end*/

      samSTPtr->cigTypeStr[samSTPtr->lenCigUI] = 'S';

      samSTPtr->cigArySI[samSTPtr->lenCigUI] =
         qrySTDoNotFree->lenSeqUL - samSTPtr->readLenUI;

      ++samSTPtr->lenCigUI;
   } /*If: have softmasking at end*/

   /*****************************************************\
   * Fun13 Sec05 Sub02:
   *   - add read lengths
   \*****************************************************/

   /*add in read lengths*/
   samSTPtr->readLenUI =
      (unsigned int) qrySTDoNotFree->lenSeqUL;

   samSTPtr->alnReadLenUI =
      samSTPtr->refEndUI - samSTPtr->refStartUI;

   /*****************************************************\
   * Fun13 Sec05 Sub03:
   *   - copy read ids
   \*****************************************************/

   /*copy read ids*/
   samSTPtr->lenRefIdUC =
      (unsigned int)
      cpWhite_ulCp(
         samSTPtr->refIdStr,
         refSTPtr->seqSTPtr->idStr
      );

   samSTPtr->lenQryIdUC =
      (unsigned int)
      cpWhite_ulCp(
         samSTPtr->qryIdStr,
         qrySTDoNotFree->idStr
      );

   /*****************************************************\
   * Fun13 Sec05 Sub04:
   *   - copy sequence
   \*****************************************************/

   if(samSTPtr->lenSeqBuffUI < qrySTDoNotFree->lenSeqUL)
   { /*If: need more memory for sequence*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         malloc(
              (qrySTDoNotFree->lenSeqUL + 9)
            * sizeof(unsigned char)
         );

      if(! samSTPtr)
         goto memErr_fun13_sec06;

      samSTPtr->lenSeqBuffUI =
        (unsigned int) qrySTDoNotFree->lenSeqUL;
   } /*If: need more memory for sequence*/

   samSTPtr->readLenUI =
      cpStr_ulCp(
         samSTPtr->seqStr,
         qrySTDoNotFree->seqStr
      );

   indexToSeq_alnSet(samSTPtr->seqStr);
     /*sequence is in index format, but I need a human
     `   readable format
     */

   /*****************************************************\
   * Fun13 Sec05 Sub05:
   *   - copy q-score
   \*****************************************************/

   if(qrySTDoNotFree->lenQUL > 0)
   { /*If: have q-score entry*/

      if(samSTPtr->lenQBuffUI < qrySTDoNotFree->lenQUL)
      { /*If: need more memory*/
         free(samSTPtr->qStr);
         samSTPtr->qStr = 0;

         samSTPtr->qStr =
            malloc(
                 (qrySTDoNotFree->lenSeqUL + 9)
               * sizeof(unsigned char)
            );

         if(! samSTPtr)
            goto memErr_fun13_sec06;

         samSTPtr->lenSeqBuffUI =
           (unsigned int) qrySTDoNotFree->lenSeqUL;
      } /*If: need more memory*/

      cpQEntry_samEntry(
         samSTPtr,
         qrySTDoNotFree->qStr,
         0 /*histgram and q-score values already blanked*/
      ); /*copy q-score entries*/

   } /*If: have q-score entry*/

   else
   { /*Else: mark no q-score entry*/
      samSTPtr->qStr[0] = '*';
      samSTPtr->qStr[1] = '\0';
   } /*Else: mark no q-score entry*/

   /*****************************************************\
   * Fun13 Sec05 Sub06:
   *   - find mutation counts
   \*****************************************************/

   /*this would be faster if I did it as I aligned,
   `   however, this is eaiser. It should also add little
   `   overhead
   */

   for(
      siCig = 0;
      siCig < (signed int) samSTPtr->lenCigUI;
      ++siCig
   ){ /*Loop: get mutation counts*/

      switch(samSTPtr->cigTypeStr[siCig])
      { /*Switch: find mutation type*/
         case '=':
         /*Case: match*/
            samSTPtr->numMatchUI +=
               samSTPtr->cigTypeStr[siCig];
            break;
         /*Case: match*/

         case 'X':
         /*Case: mismatch*/
            samSTPtr->numSnpUI +=
               samSTPtr->cigTypeStr[siCig];
            break;
         /*Case: mismatch*/

         case 'I':
         /*Case: insertion*/
            samSTPtr->numInsUI +=
               samSTPtr->cigTypeStr[siCig];
            break;
         /*Case: insertion*/

         case 'D':
         /*Case: deletion*/
            samSTPtr->numDelUI +=
               samSTPtr->cigTypeStr[siCig];
            break;
         /*Case: deletion*/

         case 'S':
         /*Case: soft mask*/
            samSTPtr->numMaskUI +=
               samSTPtr->cigTypeStr[siCig];
            break;
         /*Case: soft mask*/
      } /*Switch: find mutation type*/

   }  /*Loop: get mutation counts*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec06:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun13_sec06;

   memErr_fun13_sec06:;
      errSC = def_memErr_mapRead;
      goto ret_fun13_sec06;
   ret_fun13_sec06:;
      qrySTDoNotFree->offsetUL = 0;
      qrySTDoNotFree->endAlnUL =
         qrySTDoNotFree->lenSeqUL - 1;

      refSTPtr->seqSTPtr->offsetUL = 0;
      refSTPtr->seqSTPtr->endAlnUL =
         refSTPtr->seqSTPtr->lenSeqUL - 1;

      freeStack_dirMatrix(&matrixStackST);

      if(cigHeapArySI)
         free(cigHeapArySI);
      cigHeapArySI = 0;

      if(cigHeapStr)
         free(cigHeapStr);
      cigHeapStr = 0;

      return errSC;
} /*mergeToSam_mapRead*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
| Fun14: flipRef_mapRead
|   - flips samEntry reference to query
| Input:
|   - qrySTPtr:
|     o seqST struct pointer with query sequence to add
|   - samSTPtr:
|     o samEntry struct pointer with sam entry to flip
|       * qrySTPtr must be reference mapped to in samSTPtr
| Output:
|   - Modifies:
|     o samSTPtr to have querySTPtr as sequence
|       * most variables are modified
|       * reverse complements if samSTPtr->flagUS & 16
|         - put in foward query sequence, this function
|           reverse complements the sequence copied to
|           samSTPtr
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
signed char
flipRef_mapRead(
   struct seqST *qrySTPtr,   /*query to add in*/
   struct samEntry *samSTPtr /*to flip*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun14 TOC:
   '   - flips samEntry reference to query
   '   o fun14 sec01:
   '     - variable declarations
   '   o fun14 sec02:
   '     - make sure have enough memory for cigar
   '   o fun14 sec03:
   '     - find reference start from sam alignment
   '   o fun14 sec04:
   '     - add starting softmasking (if needed)
   '   o fun14 sec05:
   '     - swap cigar insertions and deletions (for query)
   '   o fun14 sec06:
   '     - add ending softmasking
   '   o fun14 sec07:
   '     - swap query and reference ids (names)
   '   o fun14 sec08:
   '     - add in query sequence (and length)
   '   o fun14 sec09:
   '     - add in query q-score entry
   '   o fun14 sec10:
   '     - reverse complement
   '   o fun14 sec11:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *swapSCPtr = 0;
   unsigned int qryStartUI = samSTPtr->refStartUI;
   signed int posSI = 0;
   signed int siCig = 0;
   unsigned long lenQryUL = 0;

   unsigned int lenRefUI = samSTPtr->readLenUI;
      /*for reverse complent step*/

   struct seqST swapSeqST;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec02:
   ^   - make sure have enough memory for cigar
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    init_seqST(&swapSeqST);

    /*might have to add softmasking to ends (2 extra)*/

    if(samSTPtr->lenCigUI >= samSTPtr->lenCigBuffUI - 2)
    { /*If: need to get more memory*/
       swapSCPtr =
          realloc(
             samSTPtr->cigTypeStr,
               (samSTPtr->lenCigUI + 11)
             * sizeof(signed char)
          );
       if(! swapSCPtr)
          goto memErr_fun14_sec11;
       samSTPtr->cigTypeStr = swapSCPtr;

       swapSCPtr =
          (signed char *)
          realloc(
             samSTPtr->cigArySI,
               (samSTPtr->lenCigUI + 11)
             * sizeof(signed int)
          );
       if(! swapSCPtr)
          goto memErr_fun14_sec11;
       samSTPtr->cigArySI = (signed int *) swapSCPtr;

       samSTPtr->lenCigBuffUI += 2;
    } /*If: need to get more memory*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec03:
   ^   - find reference start from sam alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*softmasking marks end of alignment*/

   if(samSTPtr->cigTypeStr[0] == 'S')
      samSTPtr->refStartUI = samSTPtr->cigArySI[0];
   else
      samSTPtr->refStartUI = 0;

   samSTPtr->refEndUI = samSTPtr->refStartUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec04:
   ^   - add starting softmasking (if needed)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(qryStartUI != 0)
   { /*If: have query soft masking*/

      if(samSTPtr->cigTypeStr[0] != 'S')
      { /*If: no soft mask entry at start*/
         siCig = samSTPtr->lenCigUI;

         while(siCig > 0)
         { /*Loop: shift cigar one up*/
            samSTPtr->cigTypeStr[siCig] =
               samSTPtr->cigTypeStr[siCig - 1];

            samSTPtr->cigArySI[siCig] =
               samSTPtr->cigArySI[siCig - 1];

            --siCig;
         } /*Loop: shift cigar one up*/

         samSTPtr->cigTypeStr[0] = 'S';
         ++samSTPtr->lenCigUI;
      } /*If: no soft mask entry at start*/

      samSTPtr->cigArySI[0] = qryStartUI - 1;
      lenQryUL += qryStartUI - 1;
   } /*If: have query soft masking*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec05:
   ^   - swap cigar insertions and deletions (for query)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   siCig = 0;

   while(siCig < (signed int) samSTPtr->lenCigUI)
   { /*Loop: flip cigars*/

      if(samSTPtr->cigTypeStr[siCig] == 'I')
      { /*If: true reference is insertion*/
         samSTPtr->cigTypeStr[siCig] = 'D';
         samSTPtr->refEndUI +=
            samSTPtr->cigTypeStr[siCig];
      } /*If: true reference is insertion*/

      else if(samSTPtr->cigTypeStr[siCig] == 'D')
      { /*Else: reference had deletion*/
         samSTPtr->cigTypeStr[siCig] = 'I';
         lenQryUL += samSTPtr->cigTypeStr[siCig];
      } /*Else: reference had deletion*/

      else if(samSTPtr->cigTypeStr[siCig] != 'S')
         lenQryUL += samSTPtr->cigTypeStr[siCig];
         /*keep track of query length*/

      ++siCig;

   } /*Loop: flip cigars*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec06:
   ^   - add ending softmasking
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   --siCig; /*get back on last entry*/

   if(lenQryUL < qrySTPtr->lenSeqUL)
   { /*If: ending soft masking*/

      if(samSTPtr->cigTypeStr[siCig] != 'S')
      { /*If: reference had no end softmasking*/
         ++siCig;
         samSTPtr->cigTypeStr[siCig] = 'S';
      } /*If: reference had no end softmasking*/

      samSTPtr->cigArySI[siCig] =
         (signed int)
         (qrySTPtr->lenSeqUL - lenQryUL);

      ++samSTPtr->lenCigUI;
   } /*If: ending soft masking*/

   else if(samSTPtr->cigTypeStr[siCig] == 'S')
      --samSTPtr->lenCigUI; /*remove entry at end*/
   
   samSTPtr->cigTypeStr[samSTPtr->lenCigUI] = '\0';
   samSTPtr->cigArySI[samSTPtr->lenCigUI] = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec07:
   ^   - swap query and reference ids (names)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   siCig = 0;

   /*swap read ids*/
   swapNull_ulCp(
      samSTPtr->qryIdStr, 
      samSTPtr->refIdStr
   );

   samSTPtr->lenQryIdUC ^= samSTPtr->lenRefIdUC;
   samSTPtr->lenRefIdUC ^= samSTPtr->lenQryIdUC;
   samSTPtr->lenQryIdUC ^= samSTPtr->lenRefIdUC;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec08:
   ^   - add in query sequence (and length)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(qrySTPtr->lenSeqUL >= samSTPtr->lenSeqBuffUI)
   { /*If: need to get more sequence buffer memory*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         malloc(
           (qrySTPtr->lenSeqUL + 9) * sizeof(signed char)
         );
      if(! samSTPtr->seqStr)
         goto memErr_fun14_sec11;
   } /*If: need to get more sequence buffer memory*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      qrySTPtr->seqStr,
      (unsigned int) qrySTPtr->lenSeqUL
   );

   samSTPtr->readLenUI= (unsigned int) qrySTPtr->lenSeqUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec09:
   ^   - add in query q-score entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(qrySTPtr->lenQUL > 0)
   { /*If: have q-score entry*/

      if(samSTPtr->lenQBuffUI < qrySTPtr->lenQUL)
      { /*If: need more memory*/

         free(samSTPtr->qStr);
         samSTPtr->qStr = 0;

         samSTPtr->qStr =
            malloc(
                 (qrySTPtr->lenSeqUL + 9)
               * sizeof(unsigned char)
            );

         if(! samSTPtr)
            goto memErr_fun14_sec11;

         samSTPtr->lenSeqBuffUI =
           (unsigned int) qrySTPtr->lenSeqUL;

      } /*If: need more memory*/

      cpQEntry_samEntry(
         samSTPtr,
         qrySTPtr->qStr,
         1 /*need to blank histgram and q-score values*/
      ); /*copy q-score entries*/

   } /*If: have q-score entry*/

   else
   { /*Else: no q-socre; need to blank*/

      samSTPtr->medianQF = 0;
      samSTPtr->meanQF = 0;
      samSTPtr->sumQUL = 0;
      samSTPtr->qStr[0] = '*';
      samSTPtr->qStr[1] = '\0';

      siCig = 0;
      while(siCig < def_maxQ_samEntry)
         samSTPtr->qHistUI[siCig++] = 0;

   } /*Else: no q-score; need to blank*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec10:
   ^   - reverse complement
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! (samSTPtr->flagUS & 16) )
      goto ret_fun14_sec11;
      /*if was a foward alignment*/
  
   /*adjust reference start and end coordinates*/
   lenRefUI = samSTPtr->refEndUI - samSTPtr->refStartUI;
   samSTPtr->refStartUI = lenRefUI - samSTPtr->refEndUI;
   samSTPtr->refEndUI = samSTPtr->refStartUI + lenQryUL;

   siCig = (signed int) samSTPtr->lenCigUI - 1;
   posSI = 0;

   while(posSI < siCig)
   { /*Loop: invert cigar*/

      /*swap cigar entries*/
      samSTPtr->cigTypeStr[siCig] ^=
         samSTPtr->cigTypeStr[posSI];
      samSTPtr->cigTypeStr[posSI] ^=
         samSTPtr->cigTypeStr[siCig];
      samSTPtr->cigTypeStr[siCig] ^=
         samSTPtr->cigTypeStr[posSI];

      samSTPtr->cigArySI[siCig] ^=
         samSTPtr->cigArySI[posSI];
      samSTPtr->cigArySI[posSI] ^=
         samSTPtr->cigArySI[siCig];
      samSTPtr->cigArySI[siCig] ^=
         samSTPtr->cigArySI[posSI];

      ++siCig;
      ++posSI;

   } /*Loop: invert cigar*/

   /*this is a bit lazy, but it works*/
   swapSeqST.seqStr = samSTPtr->seqStr;
   swapSeqST.lenSeqUL = samSTPtr->readLenUI;
   swapSeqST.qStr = samSTPtr->qStr;

   revComp_seqST(&swapSeqST);

   swapSeqST.seqStr = 0;
   swapSeqST.lenSeqUL = 0;
   swapSeqST.qStr = 0;

   samSTPtr->refEndUI =
      samSTPtr->refEndUI - samSTPtr->refStartUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec11:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ret_fun14_sec11:;
      return 0;

   memErr_fun14_sec11:;
      return def_memErr_mapRead;
} /*flipRef_mapRead*/

/*-------------------------------------------------------\
| Fun15: scoreSubAln_mapRead
|   - scores the aligment in samEntry and finds best
|     sub-aligment
| Input:
|   - samSTPtr:
|     o samEntry struct pointer with alignment to score
|     o do not trim
|   - subBl:
|     o modifiy alignment to best sub-alignment
|       * 1: look for best sub-alignment
|       * 0: do not look for sub-alignments
|   - refSTPtr:
|     o seqST struct pointer with reference sequence used
|       in alignment
|     o seqStr must be as lookup index (seqToIndex_alnSet)
|   - qrySTPtr:
|     o seqST struct pointer with query sequence used in
|       alignment
|     o seqStr must be as lookup index (seqToIndex_alnSet)
|   - alnSetSTPTr:
|     o alnSet struct pointer with scoring settings
|   - maxScoreSLPtr:
|     o signed long pointer to hold maximum score possible
|       for the aligned part of the query
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o extraStr in samSTPtr to have score "AS:i:score"
|       entry
|     o cigTypeStr and cigArySI in samSTPtr to be setup
|       for best sub-alignment (if subBl is 1)
|     o lenCigUI in samSTPtr to be new cigar length if
|       looking for best sub-alignment
|     o maxScoreSLPtr to have maximum possible score for
|       returned alignment (exlcudes unaligned bases)
|     o errSCPtr:
|       * 0 for no errors
|       * def_noScore_mapRead if no sub-score greater
|         than 0 (unlikley)
|       * def_memErr_mapRead for memory errors
|   - Returns:
|     o score for best sub-alignment or if subSl is 0, for
|       the entire samEntry
\-------------------------------------------------------*/
signed long
scoreSub_mapRead(
   struct samEntry *samSTPtr, /*alignment to score*/
   signed char subBl,         /*1:look for sub-alingment*/
   struct seqST *refSTPtr,    /*reference sequence*/
   struct seqST *qrySTPtr,    /*query sequence*/
   struct alnSet *alnSetSTPtr,/*scoring settings*/
   signed long *maxScoreSLPtr,/*gets max possible score*/
   signed char *errSCPtr      /*errors*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun15 TOC:
   '   - scores the aligment in samEntry and finds best
   '     sub-aligment
   '   o fun15 sec01:
   '     - variable declartions
   '   o fun15 sec02:
   '     - allocate memory for assembler score (if need)
   '   o fun15 sec03:
   '     - find score and best sub-score
   '   o fun15 sec04:
   '     - convert cigar to sub-alingment cigar; if subaln
   '   o fun15 sec05:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *refStr = refSTPtr->seqStr;
   signed char *qryStr = qrySTPtr->seqStr;

   signed char *swapSCPtr = 0;
   unsigned int uiCig = 0;
   unsigned int uiPos = 0;
   signed long scoreSL = 0;

   signed long subScoreSL = 0;
   unsigned int subStartUI = 0;
   signed long maxSubScoreSL = 0;
   signed long maxBestScoreSL = 0;

   signed long bestScoreSL = 0;
   unsigned int bestStartUI = 0;
   unsigned int bestEndUI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec02:
   ^   - allocate memory for assembler score (if need)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
      samSTPtr->lenExtraUI >= samSTPtr->lenExtraBuffUI -32
   ){ /*If: need more memory*/
      swapSCPtr =
         realloc(
            samSTPtr->extraStr,
              (samSTPtr->lenExtraUI + 41)
            * sizeof(signed char)
         );
      if(! swapSCPtr)
         goto memErr_fun15_sec05;
      samSTPtr->extraStr = swapSCPtr;

      samSTPtr->lenExtraBuffUI += 32;
   }  /*If: need more memory*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec03:
   ^   - find score and best sub-score
   ^   o fun15 sec03 sub01:
   ^     - get to first aligned reference base
   ^   o fun15 sec03 sub02:
   ^     - score matches/snps + start switch + start loop
   ^   o fun15 sec03 sub03:
   ^     - score insertions
   ^   o fun15 sec03 sub04:
   ^     - score deletions
   ^   o fun15 sec03 sub05:
   ^     - masking, move to next query entry
   ^   o fun15 sec03 sub06:
   ^     - add final sub-score to total score
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun15 Sec03 Sub01:
   *   - get to first aligned reference base
   \*****************************************************/

   refStr += samSTPtr->refStartUI;
   *maxScoreSLPtr = 0;

   /*****************************************************\
   * Fun15 Sec03 Sub02:
   *   - score matches/snps + start switch + start loop
   \*****************************************************/

   for(
      uiCig = 0;
      uiCig < samSTPtr->lenCigUI;
      ++uiCig
   ){ /*Loop: score cigar*/

      switch(samSTPtr->cigTypeStr[uiCig])
      { /*Switch: find cigar type*/

         case '=':
         case 'X':
         case 'M':
         /*Case: match or snp*/
            for(
               uiPos = 0;
               uiPos <
                 (unsigned int) samSTPtr->cigArySI[uiCig];
               ++uiPos
            ){ /*Loop: find score*/
               subScoreSL +=
                  getScore_alnSet(
                     *qryStr,
                     *refStr,
                     alnSetSTPtr
                  );

               maxSubScoreSL +=
                  getScore_alnSet(
                     *qryStr,
                     *qryStr,
                     alnSetSTPtr
                  );

               ++refStr;
               ++qryStr;
            }  /*Loop: find score*/

            if(subScoreSL <= 0)
            { /*If: end of sub-alignment*/
               scoreSL += subScoreSL;
               *maxScoreSLPtr += maxSubScoreSL;

               subScoreSL = 0;
               maxSubScoreSL = 0;

               subStartUI = uiCig;
            } /*If: end of sub-alignment*/

            else if(subScoreSL > bestScoreSL)
            { /*Else If: new best score*/
               bestScoreSL = subScoreSL;
               maxBestScoreSL = maxSubScoreSL;
               bestStartUI = subStartUI;
               bestEndUI = uiCig;
            } /*Else If: new best score*/

            break;
         /*Case: match or snp*/

         /***********************************************\
         * Fun15 Sec03 Sub03:
         *   - score insertions
         \************************************************/

         case 'I':
         /*Case: insertion*/

            /*get score and subscore*/
            subScoreSL += alnSetSTPtr->gapSS;

            subScoreSL +=
               (
                    (samSTPtr->cigArySI[uiCig] - 1)
                  * alnSetSTPtr->extendSS
               );

            /*get the maximum possible score*/
            for(
               uiPos = 0;
               uiPos <
                 (unsigned int) samSTPtr->cigArySI[uiCig];
               ++uiPos
            ){ /*Loop: find maximum possible score*/
               maxSubScoreSL +=
                  getScore_alnSet(
                     *qryStr,
                     *qryStr,
                     alnSetSTPtr
                  );

               ++qryStr;
            }  /*Loop: find maximum possible score*/


            /*check if end of sub-alignment*/
            if(subScoreSL <= 0)
            { /*If: end of sub-alignment*/
               scoreSL += subScoreSL;
               *maxScoreSLPtr += maxSubScoreSL;

               subScoreSL = 0;
               maxSubScoreSL = 0;

               subStartUI = uiCig;
            } /*If: end of sub-alignment*/


            break;
         /*Case: insertion*/

         /***********************************************\
         * Fun15 Sec03 Sub04:
         *   - score deletions
         \************************************************/

         case 'D':
         /*Case: deletion*/
            refStr += samSTPtr->cigArySI[uiCig];
            subScoreSL += alnSetSTPtr->gapSS;

            subScoreSL +=
               (
                    (samSTPtr->cigArySI[uiCig] - 1)
                  * alnSetSTPtr->extendSS
               );

            if(subScoreSL <= 0)
            { /*If: end of sub-alignment*/
               scoreSL += subScoreSL;
               *maxScoreSLPtr += maxSubScoreSL;

               subScoreSL = 0;
               maxSubScoreSL = 0;

               subStartUI = uiCig;
            } /*If: end of sub-alignment*/

            break;
         /*Case: deletion*/

         /***********************************************\
         * Fun15 Sec03 Sub05:
         *   - masking, move to next query entry
         \************************************************/

         default:
            qryStr += samSTPtr->cigArySI[uiCig];
      } /*Switch: find cigar type*/

   } /*Loop: score cigar*/

   /*****************************************************\
   * Fun15 Sec03 Sub06:
   *   - add final sub-score to total score
   \*****************************************************/

   scoreSL += subScoreSL;
   *maxScoreSLPtr += maxSubScoreSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec04:
   ^   - convert cigar to sub-alingment cigar (if subaln)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! subBl)
      goto ret_fun15_sec05;
   else if(bestStartUI == bestEndUI)
      goto noScore_fun15_sec05;
   else if(bestStartUI >= samSTPtr->lenCigUI)
      goto noScore_fun15_sec05;

   if(bestStartUI)
   { /*If: need to shift cigar*/
      cpLen_ulCp(
         samSTPtr->cigTypeStr,
         &samSTPtr->cigTypeStr[bestStartUI],
         bestStartUI
      );

      cpLen_ulCp(
         (signed char *) samSTPtr->cigArySI,
         (signed char *) &samSTPtr->cigArySI[bestStartUI],
         bestStartUI * sizeof(signed int)
      );

      samSTPtr->lenCigUI -= bestStartUI;
   } /*If: need to shift cigar*/

   /*add in ending coordinates; should be safe if entire
   `  alignment was best sub-alignment
   */
   samSTPtr->cigTypeStr[bestEndUI - bestStartUI +1]= '\0';
   samSTPtr->cigArySI[bestEndUI - bestStartUI + 1] = 0;
   scoreSL = bestScoreSL;
   *maxScoreSLPtr = maxBestScoreSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec05:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ret_fun15_sec05:;
      *errSCPtr = 0;

      samSTPtr->extraStr[samSTPtr->lenExtraUI++] = 'A';
      samSTPtr->extraStr[samSTPtr->lenExtraUI++] = 'S';
      samSTPtr->extraStr[samSTPtr->lenExtraUI++] = ':';
      samSTPtr->extraStr[samSTPtr->lenExtraUI++] = 'i';
      samSTPtr->extraStr[samSTPtr->lenExtraUI++] = ':';

      bestScoreSL = scoreSL; /*so can  update AS:i score*/

      if(scoreSL < 0)
      { /*If: had negative score*/
         samSTPtr->extraStr[samSTPtr->lenExtraUI++] = '-';
         bestScoreSL *= -1;
      } /*If: had negative score*/
      
      samSTPtr->lenExtraUI +=
         numToStr(
            &samSTPtr->extraStr[samSTPtr->lenExtraUI],
            (unsigned long) bestScoreSL
         );

      return scoreSL;

   noScore_fun15_sec05:;
      *errSCPtr = def_noScore_mapRead;
      return scoreSL;

   memErr_fun15_sec05:;
      *errSCPtr = def_memErr_mapRead;
      return scoreSL;
} /*scoreSubAln_mapRead*/

/*-------------------------------------------------------\
| Fun16: blank_set_mapRead
|   - blank a set_mapRead struct
| Input:
|   - setSTPtr:
|     o pointer to set_mapRead struct to blank
|   - quickBl:
|     o tells if doing quick blank (used kmers only)
|       * 1: quick blank; only blank kmers in tblSTPtr
|       * 0: slow blank; blank all of tblSTPtr
|   - tblSeqBlankBl:
|     o tells if blanking sequence in table
|       * 1: blank sequence in tblSTPtr
|       * 0: do not blank sequece in tblSTPtr
| Output:
|   - Modifies:
|     o blanks structs in setSTPtr
\-------------------------------------------------------*/
void
blank_set_mapRead(
   struct set_mapRead *setSTPtr,
   signed char quickBl,/*1: do quick blank of tblSTPtr*/
   signed char tblSeqBlankBl /*blank sequences*/
){
   
   if(! setSTPtr)
      return;

   if(setSTPtr->tblSTPtr)
   { /*If: have table to blank*/
      if(quickBl)
         qckBlank_tblST_kmerFind(
            setSTPtr->tblSTPtr,
            tblSeqBlankBl
         );
      else
         blank_tblST_kmerFind(
            setSTPtr->tblSTPtr,
            tblSeqBlankBl
         );
   } /*If: have table to blank*/

   if(setSTPtr->qrySTPtr)
      blank_refST_kmerFind(setSTPtr->qrySTPtr);

   if(setSTPtr->chainsSTPtr)
      blank_chains_mapRead(setSTPtr->chainsSTPtr);

} /*blank_set_mapRead*/

/*-------------------------------------------------------\
| Fun17: init_set_mapRead
|   - initializes a set_mapRead struct
| Input:
|   - setSTPtr:
|     o pointer to set_mapRead struct to initializes
| Output:
|   - Modifies:
|     o sets variables to defaults and pointers to null
\-------------------------------------------------------*/
void
init_set_mapRead(
   struct set_mapRead *setSTPtr
){
   if(! setSTPtr)
      return;

   setSTPtr->percExtraNtF = def_percExtraNt_defsMapRead;
   setSTPtr->percShiftF = def_percShift_defsMapRead;

   setSTPtr->lenArySI[0] = def_minLen_defsMapRead;
   setSTPtr->lenArySI[1] = def_secLen_defsMapRead;
   setSTPtr->lenArySI[2] = def_thirdLen_defsMapRead;
   setSTPtr->lenArySI[3] = def_fourthLen_defsMapRead;
   setSTPtr->lenArySI[4] = def_fifthLen_defsMapRead;
   setSTPtr->lenArySI[5] = def_sixthLen_defsMapRead;
   setSTPtr->lenArySI[6] = def_finalLen_defsMapRead;

   setSTPtr->kmerAryUC[0] = def_minKmer_defsMapRead;
   setSTPtr->kmerAryUC[1] = def_secKmer_defsMapRead;
   setSTPtr->kmerAryUC[2] = def_thirdKmer_defsMapRead;
   setSTPtr->kmerAryUC[3] = def_fourthKmer_defsMapRead;
   setSTPtr->kmerAryUC[4] = def_fifthKmer_defsMapRead;
   setSTPtr->kmerAryUC[5] = def_sixthKmer_defsMapRead;
   setSTPtr->kmerAryUC[6] = def_finalKmer_defsMapRead;

   setSTPtr->minKmersAryF[0] = def_minPerc_defsMapRead;
   setSTPtr->minKmersAryF[1] = def_secPerc_defsMapRead;
   setSTPtr->minKmersAryF[2] = def_thirdPerc_defsMapRead;
   setSTPtr->minKmersAryF[3] = def_fourthPerc_defsMapRead;
   setSTPtr->minKmersAryF[4] = def_fifthPerc_defsMapRead;
   setSTPtr->minKmersAryF[5] = def_sixthPerc_defsMapRead;
   setSTPtr->minKmersAryF[6] = def_finalPerc_defsMapRead;

   setSTPtr->chainLenArySI[0] = def_minChain_defsMapRead;
   setSTPtr->chainLenArySI[1] = def_secChain_defsMapRead;
   setSTPtr->chainLenArySI[2] =
        def_thirdChain_defsMapRead;
   setSTPtr->chainLenArySI[3] =
       def_fourthChain_defsMapRead;
   setSTPtr->chainLenArySI[4] =
      def_fifthChain_defsMapRead;
   setSTPtr->chainLenArySI[5] =
      def_sixthChain_defsMapRead;
   setSTPtr->chainLenArySI[6] =
      def_finalChain_defsMapRead;

   setSTPtr->lenKmersUC = 7;


   setSTPtr->subBl = def_subAln_defsMapRead;
   setSTPtr->minScoreF = def_minScore_defsMapRead;
   setSTPtr->minPercLenF = def_minPercLen_defsMapRead;

   setSTPtr->alnSetST = 0;
   setSTPtr->tblSTPtr = 0;
   setSTPtr->qrySTPtr = 0;
   setSTPtr->chainsSTPtr = 0;

   setSTPtr->forQryKmersArySI = 0;
   setSTPtr->revQryKmersArySI = 0;
   setSTPtr->lenQrySI = 0;
   setSTPtr->sizeQrySI = 0;

   setSTPtr->refKmersArySI = 0;
   setSTPtr->refSortArySI = 0;
   setSTPtr->lenRefSI = 0;
   setSTPtr->sizeRefSI = 0;

   /*currently adds nothing, but for future*/
   blank_set_mapRead(
      setSTPtr,
      0,        /*do slower blank (nothing to blank)*/
      1         /*blank sequences (none to blank)*/
   );
} /*init_set_mapRead*/

/*-------------------------------------------------------\
| Fun18: freeStack_set_mapRead
|   - frees variables in a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer with variables to free
| Output:
|   - Modifies:
|     o frees all pointers and intitializes; defaults/null
\-------------------------------------------------------*/
void
freeStack_set_mapRead(
   struct set_mapRead *setSTPtr
){
   if(! setSTPtr)
      return;

   if(setSTPtr->alnSetST)
      freeHeap_alnSet(setSTPtr->alnSetST);

   if(setSTPtr->tblSTPtr)
      freeHeap_tblST_kmerFind(setSTPtr->tblSTPtr);

   if(setSTPtr->qrySTPtr)
      freeHeap_refST_kmerFind(setSTPtr->qrySTPtr);

   if(setSTPtr->chainsSTPtr)
      freeHeap_chains_mapRead(setSTPtr->chainsSTPtr);


   if(setSTPtr->forQryKmersArySI)
      free(setSTPtr->forQryKmersArySI);

   if(setSTPtr->revQryKmersArySI)
      free(setSTPtr->revQryKmersArySI);

   if(setSTPtr->refKmersArySI)
      free(setSTPtr->refKmersArySI);

   if(setSTPtr->refSortArySI)
      free(setSTPtr->refSortArySI);

   init_set_mapRead(setSTPtr);
} /*freeStack_set_mapRead*/

/*-------------------------------------------------------\
| Fun19: freeHeap_set_mapRead
|   - frees a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o frees setSTPtr (you must set to null/0)
\-------------------------------------------------------*/
void
freeHeap_set_mapRead(
   struct set_mapRead *setSTPtr
){
   if(! setSTPtr)
      return;

   freeStack_set_mapRead(setSTPtr);
   free(setSTPtr);
} /*freeHeap_set_mapRead*/

/*-------------------------------------------------------\
| Fun20: setup_set_mapRead
|   - allocates memory for structs in a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer to free
|   - kmerUC:
|     o size of one kmer
|   - maxRefLenSL:
|     o size of expected longest reference
| Output:
|   - Modifies:
|     o all arrays in setSTPtr to have memory
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_set_mapRead(
   struct set_mapRead *setSTPtr, /*to setup*/
   unsigned char kmerUC,         /*size of one kmer*/
   signed long maxRefLenSL       /*size of longest ref*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun20 TOC:
   '   - allocates memory for structs in set_mapRead
   '   o fun20 sec01:
   '     - setup alnset (alnsetsTPtr) struct + input check
   '   o fun20 sec02:
   '     - setup tblsT_kmero fInd (tblsTPtr) struct
   '   o fun20 sec03:
   '     - setup refsT_kmero fInd (qrysTPtr) struct
   '   o fun20 sec04:
   '     - setup chains_mapRead (chainssTPtr) struct
   '   o fun20 sec05:
   '     - setup query kmer arrays
   '   o fun20 sec06:
   '     - setup reference kmer arrays
   '   o fun20 sec07:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec01:
   ^   - setup alnSet (alnSetSTPtr) struct + input check
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! setSTPtr)
      return 0;

   if(! setSTPtr->alnSetST)
   { /*If: do not have alignment settings*/

      setSTPtr->alnSetST = malloc(sizeof(struct alnSet));
      if(! setSTPtr->alnSetST)
         goto memErr_fun20_sec07;
      init_alnSet(setSTPtr->alnSetST);

   } /*If: do not have alignment settings*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec02:
   ^   - setup tblST_kmerfInd (tblSTPtr) struct
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(setSTPtr->tblSTPtr)
      freeHeap_tblST_kmerFind(setSTPtr->tblSTPtr);
   setSTPtr->tblSTPtr = 0;

   setSTPtr->tblSTPtr =
      malloc(sizeof(struct tblST_kmerFind));
   if(! setSTPtr->tblSTPtr)
       goto memErr_fun20_sec07;
   init_tblST_kmerFind(setSTPtr->tblSTPtr);

   if(
      setup_tblST_kmerFind(
         setSTPtr->tblSTPtr,
         kmerUC
      )
   ) goto memErr_fun20_sec07;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec03:
   ^   - setup refST_kmerFind (qrySTPtr) struct
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(setSTPtr->qrySTPtr)
      freeHeap_refST_kmerFind(setSTPtr->qrySTPtr);
   setSTPtr->qrySTPtr = 0;

   setSTPtr->qrySTPtr =
      malloc(sizeof(struct refST_kmerFind));
   if(! setSTPtr->qrySTPtr)
       goto memErr_fun20_sec07;
   init_refST_kmerFind(setSTPtr->qrySTPtr);

   if(
      setup_refST_kmerFind(
         setSTPtr->qrySTPtr,
         kmerUC
      )
   ) goto memErr_fun20_sec07;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec04:
   ^   - setup chains_mapRead (chainsSTPtr) struct
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(setSTPtr->chainsSTPtr)
      freeHeap_chains_mapRead(setSTPtr->chainsSTPtr);
   setSTPtr->chainsSTPtr = 0;

   setSTPtr->chainsSTPtr =
      malloc(sizeof(struct chains_mapRead));
   if(! setSTPtr->chainsSTPtr)
      goto memErr_fun20_sec07;
   init_chains_mapRead(setSTPtr->chainsSTPtr);

   if( 
      setup_chains_mapRead(
         setSTPtr->chainsSTPtr,
         32,
         maxRefLenSL
      )
   ) goto memErr_fun20_sec07;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec05:
   ^   - setup query kmer arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(setSTPtr->forQryKmersArySI)
       free(setSTPtr->forQryKmersArySI);
   setSTPtr->forQryKmersArySI = 0;

   setSTPtr->forQryKmersArySI =
      malloc((maxRefLenSL + 8) * sizeof(signed int));
   if(! setSTPtr->forQryKmersArySI)
      goto memErr_fun20_sec07;


   if(setSTPtr->revQryKmersArySI)
       free(setSTPtr->revQryKmersArySI);
   setSTPtr->revQryKmersArySI = 0;

   setSTPtr->revQryKmersArySI =
      malloc((maxRefLenSL + 8) * sizeof(signed int));
   if(! setSTPtr->revQryKmersArySI)
      goto memErr_fun20_sec07;


   setSTPtr->lenQrySI = 0;
   setSTPtr->sizeQrySI = (signed int) maxRefLenSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec06:
   ^   - setup reference kmer arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(setSTPtr->refKmersArySI)
       free(setSTPtr->refKmersArySI);
   setSTPtr->refKmersArySI = 0;

   setSTPtr->refKmersArySI =
      malloc((maxRefLenSL + 8) * sizeof(signed int));
   if(! setSTPtr->refKmersArySI)
      goto memErr_fun20_sec07;


   if(setSTPtr->refSortArySI)
       free(setSTPtr->refSortArySI);
   setSTPtr->refSortArySI = 0;

   setSTPtr->refSortArySI =
      malloc((maxRefLenSL + 8) * sizeof(signed int));
   if(! setSTPtr->refSortArySI)
      goto memErr_fun20_sec07;


   setSTPtr->lenRefSI = 0;
   setSTPtr->sizeRefSI = (signed int) maxRefLenSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec07:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun20_sec07:;
      return def_memErr_mapRead;
} /*setup_set_mapRead*/


/*-------------------------------------------------------\
| Fun21: addQryKmers_mapRead
|   - get forward and reverse kmer sequences
| Input:
|   - forSeqStr:
|     o c-string with forward query sequence
|   - refSeqStr:
|     o c-string with reverse query sequence
|   - lenSeqSI:
|     o sequence length
|   - forKmerArySI:
|     o signed int array pionter to hold forward kmers
|   - revKmerArySI:
|     o signed int array pionter to hold reverse kmers
|   - sizeSIPtr:
|     o signed int pointer with current for/rev array size
|   - kmerLenUC:
|     o length of one kmer (index 1)
| Output:
|   - Modifies:
|     o forKmerArySI to have forward kmers (resized if
|       needed)
|     o revKmerArySI to have reverse kmers (resized if
|       needed)
|     o sizeSIPtr to have resized array length
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
addQryKmers_mapRead(
   signed char *forSeqStr,     /*forward sequence*/
   signed char *revSeqStr,     /*reverse sequence*/
   signed int lenSeqSI,        /*sequence length*/
   signed int **forKmerArySI,  /*gets forward kmers*/
   signed int **revKmerArySI,  /*gets reverse kmers*/
   signed int *sizeSIPtr,      /*current array sizes*/
   unsigned char kmerLenUC     /*size of one kmer*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun21 TOC:
   '   - get forward and reverse kmer sequences
   '   o fun21 sec01:
   '     - variable declarations
   '   o fun21 sec02:
   '     - allocate memory for arrays
   '   o fun21 sec03:
   '     - add kmers to arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned char ntUC = 0;
   signed int siNt = 0;

   unsigned long forKmerUL = 0; /*forward kmer*/
   unsigned long lenForUL = 0;  /*number kmers*/

   unsigned long revKmerUL = 0;
   unsigned long lenRevUL = 0;  /*number kmers*/

   unsigned long maskUL = mkMask_kmerBit(kmerLenUC);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec02:
   ^   - allocate memory for arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
         ! *forKmerArySI
      || ! *revKmerArySI
      || *sizeSIPtr < lenSeqSI
   ){ /*If: need more memory*/
      if(*forKmerArySI)
         free(*forKmerArySI);
      *forKmerArySI = 0;

      *forKmerArySI =
         malloc( (lenSeqSI + 8) * sizeof(signed int) );
      if(! *forKmerArySI)
         goto memErr_fun21_sec04;


      if(*revKmerArySI)
         free(*revKmerArySI);
      *revKmerArySI = 0;

      *revKmerArySI =
         malloc( (lenSeqSI + 8) * sizeof(signed int) );
      if(! *revKmerArySI)
         goto memErr_fun21_sec04;


      *sizeSIPtr = lenSeqSI;
   }  /*If: need more memory*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec03:
   ^   - add kmers to arrays
   ^   o fun21 sec01 sub01:
   ^     - add forward kmers + start loop
   ^   o fun21 sec01 sub02:
   ^     - add reverse kmers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun21 Sec01 Sub01:
   *   - add forward kmers + start loop
   \*****************************************************/

   --kmerLenUC; /*convert to index 0*/

   for(
      siNt = 0;
      siNt < lenSeqSI;
      ++siNt
   ){ /*Loop: add kmers*/

      ntUC =
         alnNtTo_kmerBit[
           (unsigned char) forSeqStr[ siNt ]
         ];

      forKmerUL =
         ntBitToKmer_kmerBit(
            ntUC,
            forKmerUL,
            maskUL
         );

      if(ntUC < def_anonNt_kmerBit)
      { /*If: no errors/anonymous bases*/
         ++lenForUL; /*kmer length*/

         if(siNt < kmerLenUC)
            ; /*on first few bases*/
         else if(lenForUL < kmerLenUC)
            (*forKmerArySI)[siNt - kmerLenUC] =
                def_noKmer_kmerFind;
            /*had anonymous base*/
         else
            (*forKmerArySI)[siNt - kmerLenUC] =
               (signed int) forKmerUL;
      } /*If: no errors/anonymous bases*/

      else
         lenForUL = 0;

      /**************************************************\
      * Fun21 Sec01 Sub02:
      *   - add reverse kmers
      \**************************************************/

      ntUC =
         alnNtTo_kmerBit[
           (unsigned char) revSeqStr[ siNt ]
         ];

      revKmerUL =
         ntBitToKmer_kmerBit(
            ntUC,
            revKmerUL,
            maskUL
         );

      if(ntUC < def_anonNt_kmerBit)
      { /*If: no errors/anonymous bases*/
         ++lenRevUL; /*kmer length*/

         if(siNt < kmerLenUC)
            ; /*on first few bases*/
         else if(lenRevUL < kmerLenUC)
            (*revKmerArySI)[siNt - kmerLenUC] =
                def_noKmer_kmerFind;
            /*had anonymous base*/
         else
            (*revKmerArySI)[siNt - kmerLenUC] =
               (signed int) revKmerUL;
      } /*If: no errors/anonymous bases*/

      else
         lenRevUL = 0;

   }  /*Loop: add kmers*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun21_sec04:;
      return def_memErr_mapRead;
} /*addQryKmers_mapRead*/


/*TODO:
`   not doing secondary prints, so mabye kmer scans to
`   find best positions, then alignment.
`   Need to think this through
*/
/*-------------------------------------------------------\
| Fun22: aln_mapRead
|   - maps a read to a reference
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to hold alignment
|   - qrySTPtr:
|     o seqST struct pointer with query sequence
|   - refSTPtr:
|     o seqST struct pointer with reference sequence
|   - setSTPtr:
|     o pointer to set_mapRead struct with settings
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o samSTPtr to have alignment (everything modified)
|     o errSCPtr to have error message
|       * 0 for no errors
|       * def_memErr_mapRead for no errors
|       * def_noChains_mapRead if chould not get chains
|       * def_noAln_mapRead if could not find alignment
|     o these are input to avoid making each time (reduce
|       memory allcation/deallocation), so no need to
|       check
|       * everything tblSTPtr
|       * everything in refKmerFindSTPtr
|       * everything in chainsSTPtr
|   - Returns:
|     o score for alignment
\-------------------------------------------------------*/
signed char
aln_mapRead(
   struct samEntry *samSTPtr,     /*gets aligment*/
   struct seqST *qrySTPtr,        /*query sequence*/
   struct seqST *refSTPtr,        /*reference sequence*/
   struct set_mapRead *setSTPtr,  /*settings*/
   signed char *errSCPtr          /*gets errors*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun22 TOC:
   '   - maps a read to a reference
   '   o fun10 sec01:
   '     - variable declarations
   '   o fun22 sec02:
   '     - setup table and query
   '   o fun22 sec03:
   '     - kmer scan
   '   o fun22 sec04:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;
   signed char *swapSCPtr = 0;
   signed char flipBl =
      qrySTPtr->lenSeqUL > refSTPtr->lenSeqUL;
   unsigned char kmerUC = 0;
   float minKmerPercF = 0; /*min % kmers for scan*/

   signed long scoreSL = 0;
   signed long maxScoreSL = 0;
   signed long maxPossibleScoreSL = 0;

   /*holding best reference windows kmers*/
   unsigned int refStartUI = 0;

   signed long numSecSL = 0; /*number secondary alns*/

   signed long maxSL = 0; /*maxium kmers for window*/
   signed char revBl = 0; /*reverse direction best?*/
   signed char firstBl = 1; /*first kmer set*/
   signed char endBl = 0;  /*at end; terminate loop*/
   float percScoreF = 0;
   float percLenF = 0;

   signed long forKmersSL = 0;
   signed long bestForKmersSL = 0;

   signed long revKmersSL = 0;
   signed long bestRevKmersSL = 0;

   signed char winSC = 0;     /*window on*/
   signed int chainLenSI = 0; /*min chain length*/

   struct samEntry tmpSamStackST;
   struct seqST *tblSeqST = 0; /*temporary copy*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec02:
   ^   - setup table and query
   ^   o fun22 sec02 sub01:
   ^     - find best way to align sequence
   ^   o fun22 sec02 sub02:
   ^     - find kmer length for alignment
   ^   o fun22 sec02 sub03:
   ^     - table setup
   ^   o fun22 sec02 sub04:
   ^     - refST_kmerFind (query) setup
   ^   o fun22 sec02 sub05:
   ^     - reference target kmers memory allocation
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun22 Sec02 Sub01:
   *   - find best way to align sequence
   \*****************************************************/

   init_samEntry(&tmpSamStackST);

   if( setup_samEntry(&tmpSamStackST) )
      goto memErr_fun22_sec04;

   flipBl = qrySTPtr->lenSeqUL > refSTPtr->lenSeqUL;

   if(flipBl)
   { /*If swaping query and reference sequence*/
      swapSCPtr = (signed char *) qrySTPtr;
      qrySTPtr = refSTPtr;
      refSTPtr = (struct seqST *) qrySTPtr;
   } /*If swaping query and reference sequence*/

   /*****************************************************\
   * Fun22 Sec02 Sub02:
   *   - find kmer length for alignment
   \*****************************************************/

   /*not effiicent, but should be nothing compared to
   `  kmer scan step
   */
   while(
         kmerUC < setSTPtr->lenKmersUC
      &&
           setSTPtr->lenArySI[kmerUC]
         < (signed int) qrySTPtr->lenSeqUL
   ) ++kmerUC;

   /*only off one index if all lengths were less than*/
   kmerUC -= (kmerUC > 0);
   minKmerPercF = setSTPtr->minKmersAryF[kmerUC];
   chainLenSI = setSTPtr->chainLenArySI[kmerUC];
   kmerUC = setSTPtr->kmerAryUC[kmerUC];
  
   setSTPtr->qrySTPtr->lenKmerUC = kmerUC;

   /*****************************************************\
   * Fun22 Sec02 Sub03:
   *   - table setup
   \*****************************************************/

   /*needed to avoid setup call, no desier to do a full
   `  blank for a 4mb array
   */

   qckBlank_tblST_kmerFind(
      setSTPtr->tblSTPtr,
      0 /*not blanking sequence*/
   ); /*blank kmertable*/

   setSTPtr->tblSTPtr->lenKmerUC = kmerUC;
   setSTPtr->tblSTPtr->kmerMaskUL= mkMask_kmerBit(kmerUC);


   errSC =
      prep_tblST_kmerFind(
         setSTPtr->tblSTPtr,     /*has table*/
         setSTPtr->percExtraNtF, /*% extra nt in window*/
         setSTPtr->percShiftF,   /*% to shift window by*/
         qrySTPtr->lenSeqUL
      ); /*prepare table for alignment*/
      /*this will remove all old sequene kmers, but not
      `  modify the kmer table. techincally blanking
      `  kmerArySI twice, but should be minimal compare
      `  to scan
      */


   tblSeqST = setSTPtr->tblSTPtr->seqSTPtr;
   setSTPtr->tblSTPtr->seqSTPtr = refSTPtr;

   /*****************************************************\
   * Fun22 Sec02 Sub04:
   *   - refST_kmerFind (query) setup
   \*****************************************************/

   /*this step also converts sequences to index's*/
   if(
     ! addSeqToRefST_kmerFind(
        setSTPtr->tblSTPtr, /*tblST_kmerFind struct*/
        setSTPtr->qrySTPtr, /*refST_kmerFind struct*/
        qrySTPtr,           /*query sequence*/
        minKmerPercF,       /*minim % kmers to align*/
        qrySTPtr->lenSeqUL, /*longest sequence*/
        setSTPtr->alnSetST  /*scoring matrix*/
     )
   ) goto memErr_fun22_sec04;

   /*****************************************************\
   * Fun22 Sec02 Sub05:
   *   - reference target kmers memory allocation
   \*****************************************************/

   /*this is to hold best window kmers, so do not loose
   `  them or have to back shift
   */

   setSTPtr->lenRefSI =
      (signed int) setSTPtr->tblSTPtr->numKmerUI;


   if(setSTPtr->sizeRefSI < setSTPtr->lenRefSI)
   { /*If: need to free arrays*/

      if(setSTPtr->refKmersArySI)
         free(setSTPtr->refKmersArySI); 
      setSTPtr->refKmersArySI = 0;


      if(setSTPtr->refSortArySI)
         free(setSTPtr->refSortArySI); 
      setSTPtr->refSortArySI = 0;


      setSTPtr->sizeRefSI = setSTPtr->lenRefSI;

   } /*If: need to free arrays*/
   

   if(! setSTPtr->refKmersArySI)
   { /*If: need to allocate memory*/
      setSTPtr->refKmersArySI =
         malloc(
             (setSTPtr->lenRefSI + 8)
           * sizeof(signed int)
         );

      if(! setSTPtr->refKmersArySI)
         goto memErr_fun22_sec04;
   } /*If: need to allocate memory*/


   if(! setSTPtr->refSortArySI)
   { /*If: need to allocate memory*/
      setSTPtr->refSortArySI =
         malloc(
             (setSTPtr->lenRefSI + 8)
           * sizeof(signed int)
         );

      if(! setSTPtr->refSortArySI)
         goto memErr_fun22_sec04;
   } /*If: need to allocate memory*/

   /*****************************************************\
   * Fun22 Sec02 Sub06:
   *   - get query kmer sequence
   \*****************************************************/

   setSTPtr->lenQrySI =
      (signed int) setSTPtr->qrySTPtr->forSeqST->lenSeqUL;

   errSC =
      addQryKmers_mapRead(
         setSTPtr->qrySTPtr->forSeqST->seqStr,
         setSTPtr->qrySTPtr->revSeqST->seqStr,
         setSTPtr->lenQrySI,
         &setSTPtr->forQryKmersArySI,
         &setSTPtr->revQryKmersArySI,
         &setSTPtr->sizeQrySI,
         kmerUC
      ); /*add query kmer sequence to arrays*/

   if(errSC)
      goto memErr_fun22_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec03:
   ^   - kmer scan
   ^   o fun22 sec03 sub01:
   ^     - get first window (chunk)
   ^   o fun22 sec03 sub02:
   ^     - get and check matchig kmer counts + start loop
   ^   o fun22 sec03 sub03:
   ^     - may have match, find best window to align
   ^   o fun22 sec03 sub04:
   ^     - align kmers
   ^   o fun22 sec03 sub05:
   ^     - move to next window
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun22 Sec03 Sub01:
   *   - get first window (chunk)
   \*****************************************************/

   endBl = 0;

   errSC =
      nextSeqChunk_tblST_kmerFind(
         setSTPtr->tblSTPtr,
         &firstBl
      );

   if(errSC)
   { /*If: have no more windows*/
      endBl = 1;
      goto kmerCount_fun22_sec03_sub0x;
   } /*If: have no more windows*/

   /*****************************************************\
   * Fun22 Sec03 Sub02:
   *   - get and check matchig kmer counts + start loop
   \*****************************************************/
   
   while(! errSC)
   { /*Loop: find alignments*/

      kmerCount_fun22_sec03_sub0x:;

      forKmersSL =
         forCntMatchs_kmerFind(
            setSTPtr->tblSTPtr,
            setSTPtr->qrySTPtr
         ); /*forward kmer count*/

      revKmersSL =
         revCntMatchs_kmerFind(
            setSTPtr->tblSTPtr,
            setSTPtr->qrySTPtr
         ); /*reverse kmer count*/


      /*realing on max here being unpredictiable, so
      `   that branch prediction is bad here
      */
      maxSL =
         max_genMath(      
            forKmersSL,
            revKmersSL
         ); /*direction with most kmers*/
      revBl = (forKmersSL < revKmersSL);

      if(
           maxSL
         < (signed long) setSTPtr->qrySTPtr->minKmersUI
      ) goto nextWin_fun22_sec03_sub05;
         /*nothing here*/

      /**************************************************\
      * Fun22 Sec03 Sub03:
      *   - may have match, find best window to align
      *   o fun22 sec03 sub03 cat01:
      *     - get first windows kmers + scores
      *   o fun22 sec03 sub03 cat02:
      *     - check if next windows are better
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun22 Sec03 Sub03 Cat01:
      +   - get first windows kmers + scores
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      bestForKmersSL = forKmersSL;
      bestRevKmersSL = revKmersSL;

      cpLen_ulCp(
         (signed char *) setSTPtr->refKmersArySI,
         (signed char *) setSTPtr->tblSTPtr->kmerArySI,
         (setSTPtr->tblSTPtr->ntInWinUI + 1)
           * sizeof(signed int)
          /*+1 to get -2 kmer*/
      );

      /*get reference starting coordiante*/
      refStartUI = setSTPtr->tblSTPtr->seqPosUL;
      refStartUI -= setSTPtr->tblSTPtr->ntInWinUI;
         /*account for seqPosUL being end of window*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun22 Sec03 Sub03 Cat02:
      +   - check next windows are better
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(! errSC)
      { /*If: not at end*/

         for(
            winSC = 0;
            winSC < setSTPtr->numWinSC;
            ++winSC
         ){ /*Loop: find best overlapping window*/

            errSC =
               nextSeqChunk_tblST_kmerFind(
                  setSTPtr->tblSTPtr,
                  &firstBl
               ); /*move to next window*/

            bestForKmersSL =
               max_genMath(      
                  forKmersSL,
                  bestForKmersSL
               );

            bestRevKmersSL =
               max_genMath(      
                  revKmersSL,
                  bestRevKmersSL
               );
  
            forKmersSL =
               max_genMath(
                  bestForKmersSL,
                  bestRevKmersSL
               );

            if(forKmersSL > maxSL)
            { /*If: new best*/
               maxSL = forKmersSL;
               revBl = (maxSL == bestRevKmersSL);

               cpLen_ulCp(
                  (signed char *) setSTPtr->refKmersArySI,
                  (signed char *)
                     setSTPtr->tblSTPtr->kmerArySI,
                  (setSTPtr->tblSTPtr->ntInWinUI + 1)
                    * sizeof(signed int)
                   /*+1 to get -2 kmer*/
               );
                  
               /*get reference starting coordiante*/
               refStartUI = setSTPtr->tblSTPtr->seqPosUL;
               refStartUI -=
                   setSTPtr->tblSTPtr->ntInWinUI;
                  /*seqPosUL is at end of window*/
            } /*If: new best*/

            if(errSC)
                break; /*no more windows to scan*/

         } /*Loop: find best overlapping window*/
      } /*If: not at end*/

      /**************************************************\
      * Fun22 Sec03 Sub04:
      *   - align kmers
      *   o fun22 sec03 sub04 cat01:
      *     - get chains and merge chains
      *   o fun22 sec03 sub04 cat02:
      *     - fill in merged chain gaps
      *   o fun22 sec03 sub04 cat03:
      *     - adjust coordinates, flip if needed, & score
      *   o fun22 sec03 sub04 cat04:
      *     - check if keep chain (min score + new best)
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun22 Sec03 Sub04 Cat01:
      +   - get chains and merge chains
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(revBl)
      { /*If: reverse complement best alignment*/
         if(
            kmerChain_mapRead(
               setSTPtr->chainsSTPtr,
               chainLenSI,
               setSTPtr->revQryKmersArySI,
               setSTPtr->lenQrySI,
               setSTPtr->refKmersArySI,
               setSTPtr->refSortArySI,
               setSTPtr->lenRefSI,
               (signed int) refStartUI
            )
         ) goto memErr_fun22_sec04;
      } /*If: reverse complement best alignment*/

      else
      { /*Else: forward best alignment*/
         if(
            kmerChain_mapRead(
               setSTPtr->chainsSTPtr,
               chainLenSI,
               setSTPtr->forQryKmersArySI,
               setSTPtr->lenQrySI,
               setSTPtr->refKmersArySI,
               setSTPtr->refSortArySI,
               setSTPtr->lenRefSI,
               (signed int) refStartUI
            )
         ) goto memErr_fun22_sec04;
      } /*Else: foward best alignment*/



      scoreSL= mergeChains_mapRead(setSTPtr->chainsSTPtr);

      if(scoreSL == def_noChains_mapRead)
         goto nextWin_fun22_sec03_sub05;
         /*no chains to merge*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun22 Sec03 Sub04 Cat02:
      +   - fill in merged chain gaps
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/
      
      if(
         mergeToSam_mapRead(
            setSTPtr->chainsSTPtr,
            (signed int) scoreSL,
            setSTPtr->tblSTPtr,
            setSTPtr->qrySTPtr,
            revBl,
            chainLenSI,
            setSTPtr->alnSetST,
            &tmpSamStackST
         ) /*make merged chains into an alignment*/
      ) goto memErr_fun22_sec04;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun22 Sec03 Sub04 Cat03:
      +   - adjust coordinates, flip if needed, and score
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/
      
      /*adjust for window not being reference start*/
      samSTPtr->refStartUI += refStartUI;
      samSTPtr->refEndUI += refStartUI;


      if(flipBl)
      { /*If: need to flip reference and query*/
         if(
            flipRef_mapRead(
               refSTPtr, /*is really query in this case*/
               &tmpSamStackST
            )
         ) goto memErr_fun22_sec04;

         swapSCPtr = (signed char *) refSTPtr;
         refSTPtr = qrySTPtr;
         qrySTPtr = (struct seqST *) swapSCPtr;

         if(samSTPtr->flagUS & 16)
         { /*If: is reverse complement alignment*/
            revCmpIndex_alnSet(
               qrySTPtr->seqStr,
               qrySTPtr->qStr,
               qrySTPtr->lenSeqUL
            );
         } /*If: is reverse complement alignment*/
      } /*If: need to flip reference and query*/


      scoreSL =
         scoreSub_mapRead(
            &tmpSamStackST,
            setSTPtr->subBl,  /*if doing sub-alignment*/
            refSTPtr,
            qrySTPtr,
            setSTPtr->alnSetST,
            &maxPossibleScoreSL,
            &errSC
         );

      if(flipBl)
      { /*If: need to swap qeury and reference*/
         percLenF = (float) samSTPtr->numMatchUI;
         percLenF += (float) samSTPtr->numSnpUI;
         percLenF += (float) samSTPtr->numInsUI;
         percLenF /= (float) qrySTPtr->lenSeqUL;
        

         swapSCPtr = (signed char *) refSTPtr;
         refSTPtr = qrySTPtr;
         qrySTPtr = (struct seqST *) swapSCPtr;

         if(samSTPtr->flagUS & 16)
         { /*If: is reverse complement alignment*/
            revCmpIndex_alnSet(
               qrySTPtr->seqStr,
               qrySTPtr->qStr,
               qrySTPtr->lenSeqUL
            );
         } /*If: is reverse complement alignment*/
      } /*If: need to swap qeury and reference*/

      if(errSC)
      { /*If: had error*/
         if(errSC == def_noScore_mapRead)
            goto nextWin_fun22_sec03_sub05;
         else
            goto memErr_fun22_sec04;
      } /*If: had error*/

      if(! maxPossibleScoreSL)
         goto nextWin_fun22_sec03_sub05;
         /*should never happen, but in case*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun22 Sec03 Sub04 Cat04:
      +   - check if keep chain (min score + new best)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      percScoreF = (float) scoreSL;
      percScoreF /= (float) maxPossibleScoreSL;

      /*direction does not matter here*/
      if(percScoreF < setSTPtr->minScoreF)
         ; /*do nothing; score beneath min value*/

      else if(percLenF < setSTPtr->minPercLenF)
         ; /*do nothing; score beneath min value*/
    
      else if(maxScoreSL < scoreSL)
      { /*Else If: better score*/
         ++numSecSL;

         swap_samEntry(
            &tmpSamStackST,
            samSTPtr
         );
      } /*Else If: better score*/

      else
         ++numSecSL; /*secondary alignment, not better*/

      blank_samEntry(&tmpSamStackST);

      /**************************************************\
      * Fun22 Sec03 Sub05:
      *   - move to next window
      \**************************************************/

      nextWin_fun22_sec03_sub05:;
         if(endBl)
            break;

         errSC =
            nextSeqChunk_tblST_kmerFind(
               setSTPtr->tblSTPtr,
               &firstBl
            ); /*move to next window*/
   }  /*Loop: find alignments*/

   if(! endBl)
   { /*If: on last window*/
      endBl = 1;
      goto kmerCount_fun22_sec03_sub0x;
   } /*If: on last window*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! numSecSL)
      *errSCPtr = def_noAln_mapRead;
   else
      *errSCPtr = 0;

   goto ret_fun22_sec04;

   memErr_fun22_sec04:;
      *errSCPtr = def_memErr_mapRead;
      goto ret_fun22_sec04;

   ret_fun22_sec04:;
      freeStack_samEntry(&tmpSamStackST);
      setSTPtr->tblSTPtr->seqSTPtr = tblSeqST;
      return maxScoreSL;
} /*aln_mapRead*/

/*-------------------------------------------------------\
| Fun23: qckAln_mapRead
|   - maps a read to a reference; returns best alignment
|   - qckAln_mapRead differs from aln_mapRead in that
|     it uses kmerScan to find best window in entire
|     reference, and then after scan, makes alignment
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to hold alignment
|   - qrySTPtr:
|     o seqST struct pointer with query sequence
|     o needs to be converted with seqToIndex_alnSet()
|   - refSTPtr:
|     o seqST struct pointer with reference sequence
|     o needs to be converted with seqToIndex_alnSet()
|   - setSTPtr:
|     o pointer to set_mapRead struct with settings
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o samSTPtr to have alignment (everything modified)
|     o errSCPtr to have error message
|       * 0 for no errors
|       * def_memErr_mapRead for no errors
|       * def_noChains_mapRead if chould not get chains
|       * def_noAln_mapRead if could not find alignment
|     o these are input to avoid making each time (reduce
|       memory allcation/deallocation), so no need to
|       check
|       * everything tblSTPtr
|       * everything in refKmerFindSTPtr
|       * everything in chainsSTPtr
|   - Returns:
|     o score for alignment
\-------------------------------------------------------*/
signed long
qckAln_mapRead(
   struct samEntry *samSTPtr,     /*gets aligment*/
   struct seqST *qrySTPtr,        /*query sequence*/
   struct seqST *refSTPtr,        /*reference sequence*/
   struct set_mapRead *setSTPtr,  /*settings*/
   signed char *errSCPtr          /*gets errors*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun23 TOC:
   '   - maps a read to a reference
   '   o fun10 sec01:
   '     - variable declarations
   '   o fun23 sec02:
   '     - setup table and query
   '   o fun23 sec03:
   '     - kmer scan
   '   o fun23 sec04:
   '     - align kmers
   '   o fun23 sec05:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;
   signed char *swapSCPtr = 0;
   signed char flipBl =
      qrySTPtr->lenSeqUL > refSTPtr->lenSeqUL;
   unsigned char kmerUC = 0;
   float minKmerPercF = 0; /*min % kmers for scan*/

   signed long scoreSL = 0;
   signed long maxScoreSL = 0; /*maximum possible score*/

   /*holding best reference windows kmers*/
   unsigned int refStartUI = 0;
   signed char firstBl = 1; /*first chunk*/

   signed long maxKmerSL = 0;
   signed long forKmersSL = 0; /*forward win kmer count*/
   signed long revKmersSL = 0; /*reverse win kmer count*/

   signed long maxSL = 0; /*maxium kmers for window*/
   signed char revBl = 0; /*reverse direction best?*/
   signed char endBl = 0;  /*at end; terminate loop*/
   signed long chainLenSI = 0; /*minimum chain length*/
   float percScoreF = 0;

   struct seqST *tblSeqST = 0; /*temporary copy*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec02:
   ^   - setup table and query
   ^   o fun23 sec02 sub01:
   ^     - find best way to align sequence
   ^   o fun23 sec02 sub02:
   ^     - find kmer length for alignment
   ^   o fun23 sec02 sub03:
   ^     - table setup
   ^   o fun23 sec02 sub04:
   ^     - refST_kmerFind (query) setup
   ^   o fun23 sec02 sub05:
   ^     - reference target kmers memory allocation
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun23 Sec02 Sub01:
   *   - find best way to align sequence
   \*****************************************************/

   flipBl = qrySTPtr->lenSeqUL > refSTPtr->lenSeqUL;

   if(flipBl)
   { /*If swaping query and reference sequence*/
      swapSCPtr = (signed char *) qrySTPtr;
      qrySTPtr = refSTPtr;
      refSTPtr = (struct seqST *) qrySTPtr;
   } /*If swaping query and reference sequence*/

   /*****************************************************\
   * Fun23 Sec02 Sub02:
   *   - find kmer length for alignment
   \*****************************************************/

   /*not effiicent, but should be nothing compared to
   `  kmer scan step
   */
   while(
         kmerUC < setSTPtr->lenKmersUC
      &&
           setSTPtr->lenArySI[kmerUC]
         < (signed int) qrySTPtr->lenSeqUL
   ) ++kmerUC;

   /*always one kmer off*/
   kmerUC -= (kmerUC > 0);
   minKmerPercF = setSTPtr->minKmersAryF[kmerUC];
   chainLenSI = setSTPtr->chainLenArySI[kmerUC];
   kmerUC = setSTPtr->kmerAryUC[kmerUC];
  
   setSTPtr->qrySTPtr->lenKmerUC = kmerUC;

   /*****************************************************\
   * Fun23 Sec02 Sub03:
   *   - table setup
   \*****************************************************/

   /*needed to avoid setup call, no desier to do a full
   `  blank for a 4mb array
   */

   qckBlank_tblST_kmerFind(
      setSTPtr->tblSTPtr,
      0 /*not blanking sequence*/
   ); /*blank kmertable*/

   setSTPtr->tblSTPtr->lenKmerUC = kmerUC;
   setSTPtr->tblSTPtr->kmerMaskUL= mkMask_kmerBit(kmerUC);


   errSC =
      prep_tblST_kmerFind(
         setSTPtr->tblSTPtr,     /*has table*/
         setSTPtr->percExtraNtF, /*% extra nt in window*/
         setSTPtr->percShiftF,   /*% to shift window by*/
         qrySTPtr->lenSeqUL
      ); /*prepare table for alignment*/
      /*this will remove all old sequene kmers, but not
      `  modify the kmer table. techincally blanking
      `  kmerArySI twice, but should be minimal compare
      `  to scan
      `  It also allocates memory if table is to small
      */

   tblSeqST = setSTPtr->tblSTPtr->seqSTPtr;
   setSTPtr->tblSTPtr->seqSTPtr = refSTPtr;

   /*****************************************************\
   * Fun23 Sec02 Sub04:
   *   - refST_kmerFind (query) setup
   \*****************************************************/

   /*this step also converts sequences to index's*/
   if(
     ! addSeqToRefST_kmerFind(
        setSTPtr->tblSTPtr, /*tblST_kmerFind struct*/
        setSTPtr->qrySTPtr, /*refST_kmerFind struct*/
        qrySTPtr,           /*query sequence*/
        minKmerPercF,       /*minim % kmers to align*/
        qrySTPtr->lenSeqUL, /*longest sequence*/
        setSTPtr->alnSetST  /*scoring matrix*/
     )
   ) goto memErr_fun23_sec05;

   /*****************************************************\
   * Fun23 Sec02 Sub05:
   *   - reference target kmers memory allocation
   \*****************************************************/

   /*this is to hold best window kmers, so do not loose
   `  them or have to back shift
   */

   setSTPtr->lenRefSI =
      (signed int) setSTPtr->tblSTPtr->numKmerUI;


   if(setSTPtr->sizeRefSI < setSTPtr->lenRefSI)
   { /*If: need to free arrays*/

      if(setSTPtr->refKmersArySI)
         free(setSTPtr->refKmersArySI); 
      setSTPtr->refKmersArySI = 0;


      if(setSTPtr->refSortArySI)
         free(setSTPtr->refSortArySI); 
      setSTPtr->refSortArySI = 0;


      setSTPtr->sizeRefSI = setSTPtr->lenRefSI;

   } /*If: need to free arrays*/
   

   if(! setSTPtr->refKmersArySI)
   { /*If: need to allocate memory*/
      setSTPtr->refKmersArySI =
         malloc(
             (setSTPtr->lenRefSI + 8)
           * sizeof(signed int)
         );

      if(! setSTPtr->refKmersArySI)
         goto memErr_fun23_sec05;
   } /*If: need to allocate memory*/


   if(! setSTPtr->refSortArySI)
   { /*If: need to allocate memory*/
      setSTPtr->refSortArySI =
         malloc(
             (setSTPtr->lenRefSI + 8)
           * sizeof(signed int)
         );

      if(! setSTPtr->refSortArySI)
         goto memErr_fun23_sec05;
   } /*If: need to allocate memory*/

   /*****************************************************\
   * Fun23 Sec02 Sub06:
   *   - get query kmer sequence
   \*****************************************************/

   setSTPtr->lenQrySI =
      (signed int) setSTPtr->qrySTPtr->forSeqST->lenSeqUL;

   errSC =
      addQryKmers_mapRead(
         setSTPtr->qrySTPtr->forSeqST->seqStr,
         setSTPtr->qrySTPtr->revSeqST->seqStr,
         setSTPtr->lenQrySI,
         &setSTPtr->forQryKmersArySI,
         &setSTPtr->revQryKmersArySI,
         &setSTPtr->sizeQrySI,
         kmerUC
      ); /*add query kmer sequence to arrays*/

   if(errSC)
      goto memErr_fun23_sec05;


   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec03:
   ^   - kmer scan
   ^   o fun23 sec03 sub01:
   ^     - get first window (chunk)
   ^   o fun23 sec03 sub02:
   ^     - get and check matchig kmer counts + start loop
   ^   o fun23 sec03 sub03:
   ^     - have new best, record
   ^   o fun23 sec03 sub04:
   ^     - align kmers
   ^   o fun23 sec03 sub05:
   ^     - move to next window
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun23 Sec03 Sub01:
   *   - get first window (chunk)
   \*****************************************************/

   endBl = 0;

   errSC =
      nextSeqChunk_tblST_kmerFind(
         setSTPtr->tblSTPtr,
         &firstBl
      );

   if(errSC)
   { /*If: have no more windows*/
      endBl = 1;
      goto kmerCount_fun23_sec03_sub0x;
   } /*If: have no more windows*/

   /*****************************************************\
   * Fun23 Sec03 Sub02:
   *   - get and check matchig kmer counts + start loop
   \*****************************************************/
   
   while(! errSC)
   { /*Loop: find alignments*/

      kmerCount_fun23_sec03_sub0x:;

      forKmersSL =
         forCntMatchs_kmerFind(
            setSTPtr->tblSTPtr,
            setSTPtr->qrySTPtr
         ); /*forward kmer count*/

      revKmersSL =
         revCntMatchs_kmerFind(
            setSTPtr->tblSTPtr,
            setSTPtr->qrySTPtr
         ); /*reverse kmer count*/


      /*realing on max here being unpredictiable, so
      `   that branch prediction is bad here
      */
      maxSL =
         max_genMath(      
            forKmersSL,
            revKmersSL
         ); /*direction with most kmers*/
      revBl = (forKmersSL < revKmersSL);

      if(
           maxSL
         < (signed long) setSTPtr->qrySTPtr->minKmersUI
      ) goto nextWin_fun23_sec03_sub05;
         /*nothing here*/

      /**************************************************\
      * Fun23 Sec03 Sub03:
      *   - have new best, record
      \**************************************************/

      if(maxSL <= maxKmerSL)
         goto nextWin_fun23_sec03_sub05;

      maxKmerSL = maxSL;

      cpLen_ulCp(
         (signed char *) setSTPtr->refKmersArySI,
         (signed char *) setSTPtr->tblSTPtr->kmerArySI,
         (setSTPtr->tblSTPtr->ntInWinUI + 1)
           * sizeof(signed int)
          /*+1 to get -2 kmer*/
      );

      /*get reference starting coordiante*/
      refStartUI = setSTPtr->tblSTPtr->seqPosUL;
      refStartUI -= setSTPtr->tblSTPtr->ntInWinUI;
         /*account for seqPosUL being end of window*/

      /**************************************************\
      * Fun23 Sec03 Sub05:
      *   - move to next window
      \**************************************************/

      nextWin_fun23_sec03_sub05:;
         if(endBl)
            break;

         errSC =
            nextSeqChunk_tblST_kmerFind(
               setSTPtr->tblSTPtr,
               &firstBl
            ); /*move to next window*/
   }  /*Loop: find alignments*/

   if(! endBl)
   { /*If: on last window*/
      endBl = 1;
      goto kmerCount_fun23_sec03_sub0x;
   } /*If: on last window*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec04:
   ^   - align kmers
   ^   o fun23 sec03 sub01:
   ^     - get chains and merge chains
   ^   o fun23 sec03 sub02:
   ^     - fill in merged chain gaps
   ^   o fun23 sec03 sub03:
   ^     - adjust coordinates, flip if needed, & score
   ^   o fun23 sec03 sub04:
   ^     - check if keep chain (min score + new best)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*+++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun23 Sec04 Sub01:
   +   - get chains and merge chains
   \+++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(! maxKmerSL)
      goto noAln_fun23_sec05;
      /*no match found*/

   if(revBl)
   { /*If: reverse complement best alignment*/
      if(
         kmerChain_mapRead(
            setSTPtr->chainsSTPtr,
            chainLenSI,
            setSTPtr->revQryKmersArySI,
            setSTPtr->lenQrySI,
            setSTPtr->refKmersArySI,
            setSTPtr->refSortArySI,
            setSTPtr->lenRefSI,
            (signed int) refStartUI
         )
      ) goto memErr_fun23_sec05;
   } /*If: reverse complement best alignment*/

   else
   { /*Else: forward best alignment*/
      if(
         kmerChain_mapRead(
            setSTPtr->chainsSTPtr,
            chainLenSI,
            setSTPtr->forQryKmersArySI,
            setSTPtr->lenQrySI,
            setSTPtr->refKmersArySI,
            setSTPtr->refSortArySI,
            setSTPtr->lenRefSI,
            (signed int) refStartUI
         )
      ) goto memErr_fun23_sec05;
   } /*Else: foward best alignment*/

   scoreSL= mergeChains_mapRead(setSTPtr->chainsSTPtr);

   if(scoreSL == def_noChains_mapRead)
      goto noAln_fun23_sec05;
      /*no chains to merge*/

   /**************************************************\
   * Fun23 Sec03 Sub02:
   *   - fill in merged chain gaps
   \**************************************************/
   
   if(
      mergeToSam_mapRead(
         setSTPtr->chainsSTPtr,
         (signed int) scoreSL,
         setSTPtr->tblSTPtr,
         setSTPtr->qrySTPtr,
         revBl,
         chainLenSI,
         setSTPtr->alnSetST,
         samSTPtr
      ) /*make merged chains into an alignment*/
   ) goto memErr_fun23_sec05;

   /**************************************************\
   * Fun23 Sec03 Sub03:
   *   - adjust coordinates, flip if needed, and score
   \**************************************************/
   
   /*adjust for window not being reference start*/
   samSTPtr->refStartUI += refStartUI;
   samSTPtr->refEndUI += refStartUI;


   if(flipBl)
   { /*If: need to flip reference and query*/
      if(
         flipRef_mapRead(
            refSTPtr, /*is really query in this case*/
            samSTPtr
         )
      ) goto memErr_fun23_sec05;

      swapSCPtr = (signed char *) refSTPtr;
      refSTPtr = qrySTPtr;
      qrySTPtr = (struct seqST *) swapSCPtr;

      if(samSTPtr->flagUS & 16)
      { /*If: is reverse complement alignment*/
         revCmpIndex_alnSet(
            qrySTPtr->seqStr,
            qrySTPtr->qStr,
            qrySTPtr->lenSeqUL
         );
      } /*If: is reverse complement alignment*/
   } /*If: need to flip reference and query*/


   scoreSL =
      scoreSub_mapRead(
         samSTPtr,
         setSTPtr->subBl,  /*if doing sub-alignment*/
         refSTPtr,
         qrySTPtr,
         setSTPtr->alnSetST,
         &maxScoreSL,
         &errSC
      );


   if(samSTPtr->flagUS & 16)
   { /*If: is reverse complement alignment*/
      revCmpIndex_alnSet(
         qrySTPtr->seqStr,
         qrySTPtr->qStr,
         qrySTPtr->lenSeqUL
      );
   } /*If: is reverse complement alignment*/

   if(errSC)
   { /*If: had error*/
      if(errSC == def_noScore_mapRead)
         goto noAln_fun23_sec05;
      else
         goto memErr_fun23_sec05;
   } /*If: had error*/

   if(! maxScoreSL)
      goto noAln_fun23_sec05;
      /*should never happen, but in case*/

   /**************************************************\
   * Fun23 Sec03 Sub04:
   *   - check if keep chain (min score * new best)
   \**************************************************/

   percScoreF = (float) scoreSL;
   percScoreF /= (float) maxScoreSL;

   /*direction does not matter here*/
   if(percScoreF < setSTPtr->minScoreF)
      goto noAln_fun23_sec05;
      /*do nothing; score beneath min value*/


   percScoreF = (float) samSTPtr->numMatchUI;
   percScoreF += (float) samSTPtr->numSnpUI;
   percScoreF += (float) samSTPtr->numInsUI;
   percScoreF /= (float) qrySTPtr->lenSeqUL;

   if(percScoreF < setSTPtr->minPercLenF)
      goto noAln_fun23_sec05;
      /*do nothing; length beneath min value*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec05:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *errSCPtr = 0;
   goto ret_fun23_sec05;

   noAln_fun23_sec05:;
      *errSCPtr = def_noAln_mapRead;
      goto ret_fun23_sec05;

   memErr_fun23_sec05:;
      *errSCPtr = def_memErr_mapRead;
      goto ret_fun23_sec05;

   ret_fun23_sec05:;
      setSTPtr->tblSTPtr->seqSTPtr = tblSeqST;
      return scoreSL;
} /*aln_mapRead*/

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
