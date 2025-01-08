/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' clustST SOF: Start Of File
'   - holds structures used in clustST and clustST
'     (for clustering)
'   o TOF: Table Of Functions
'     - tof00: header entry
'     - tof01: structures in this file
'     - tof02: set_clustST support fuctions
'     - tof03: index_clustST general functions
'     - tof04: con_clust general functions
'     - tof05: build/use index_clustST; consensus to
'     - tof06: consensus and cluster functions
'     - tof07: var_clust general
'     - tof08: hist_clust support + probability functions
'   - TOF00: header entry
'     o header:
'       - included libraries
'   - TOF01: structures
'     o .h st01: set_clustST
'       - has settings for clustering
'     o .h st02: index_clustST
'       - holds scores, clusters, and index's for sam file
'       - ~ 29 bytes per line
'     o .h st03: con_clustST
'       - holds a single consensus/profile for a cluster
'     o .h st04: var_clustST
'       - holds depth of variants in reference (clustST)
'     o .h st05: hist_clustST
'       - holds profile history for reads (clustST)
'   - TOF02: set_clustST support fuctions
'     o fun01: blank_set_clustST
'       - sets start and end to 0
'     o fun02: init_set_clustST
'       - sets settings in set_clustST struct to defaults
'     o fun03: freeStack_set_clustST
'       - no internal variables to free (here for support)
'     o fun04: freeHeap_set_clustST
'       - frees a heap allocated set_clustST struct
'   - TOF03: index_clustST general functions
'     o fun05: blank_index_clustST
'       - blanks an index_clustST struct
'     o fun06: init_index_clustST
'       - initializes an index_clustST struct
'     o fun07: freeStack_index_clustST
'       - frees and initializes variables in index_clustST
'     o fun09: freeHeap_index_clustST
'       - frees index_clustST struct
'     o fun10: setup_index_clustST
'       - setsup memory for index_clustST struct
'     o fun11: realloc_index_clustST
'       - reallocates memory for an index_clustST struct
'   - TOF04: con_clust general functions
'     o fun12: blank_con_clustST
'       - blank an clustST con_clustST struct
'     o fun13: init_con_clustST
'       - initialzie an clustST con_clustST struct
'     o fun14: freeStack_con_clustST
'       - frees variables inside an con_clustST struct
'     o fun15: freeHeap_con_clustST
'       - frees a con_clustST struct
'     o fun16: freeHeapList_con_clustST
'       - frees list of heap allocated con_clustST structs
'     o fun17: mk_con_clustST
'       - makes a con_clustST struct on the heap
'     o fun18: swap_con_clustST
'       - swaps to con_clustST structs (except nextST ptr)
'   - TOF05: build/use index_clustST struct
'     o fun19: getNumLines_clustST
'       - finds number of lines in a file
'     o fun20: mk_index_clustST
'       - get read scores for a sam file
'     o fun21: getRead_clustST
'       - gets a read by index from a sam file
'   - TOF06: consensus and cluster functions
'     o fun22: getCon_clustST
'       - builds a consensus for a cluster
'     o fun23: cmpCons_clustST
'       - compares two consensus clusters to see if same
'     o fun24: plist_con_clustST
'      - print clusters consensuses in a con_clustST list
'     o fun25: getClust_clustST
'       - extracts reads for a signle cluster
'     o fun26: pbins_clustST
'       - extracts reads for each cluster to a sam file
'   - TOF07: var_clust general
'     o fun27: blank_var_clustST
'       - blanks var_clustST for use with a new reference
'     o fun28: init_var_clustST
'       - initializes a var_clustST array
'     o fun29: freeStack_var_clustST
'       - frees arrays in a var_clustST struct
'     o fun30: freeHeap_var_clustST
'       - frees a heap allocated var_clustST struct
'     o fun31: freeHeapAry_var_clustST
'       - frees an array of heap allocated var_clustST
'     o fun32: setupRef_var_clustST
'       - setsup a var_clustST for a reference
'     o fun33: rmLowDepth_var_clustST
'       - remove low read depth variants in a var_clustST
'   - TOF08: hist_clust support + probability functions
'     o fun34: blank_hist_clustST
'       - blanks hist_clustST for use
'     o fun35: init_hist_clustST
'       - initializes hist_clustST for use
'     o fun36: freeStack_hist_clustST
'       - frees variables inside a hist_clustST
'     o fun37: freeHeap_hist_clustST
'       - frees a hist_clustST struct
'     o fun38: prob_clustST
'       - find probability that a profile is supported
'     o fun39: roughProb_clustST
'       - get expected power of 10 of profile being unique
'     o fun40: findBestProf_hist_clustST
'       - finds best profile in a hist_clustST history
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "clustST.h"

#include <stdio.h>

#include "../genLib/numToStr.h"
#include "../genLib/ulCp.h"
#include "../genLib/strAry.h"

#include "../genBio/samEntry.h"
#include "../genBio/edDist.h"
#include "../genBio/tbCon.h"

/*.h files only (no .c files)*/
#include "../genLib/dataTypeShortHand.h"
#include "../genLib/genMath.h" /*only using .h min macro*/
#include "../genBio/tbConDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../genLib/base10str.h"
!   - .c  #include "../genLib/strAry.h"
!   - .h  #include "../genBio/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: blank_set_clustST
|   - sets start and end to 0
| Input:
|   - clustSetSTPtr;
|     o pointer to set_clustST struct to blank
| Output:
|   - Modifies:
|     o startUI and endUI in clustSetSTPtr to be 0
\-------------------------------------------------------*/
void
blank_set_clustST(
   struct set_clustST *clustSetSTPtr
){
   clustSetSTPtr->startUI = 0;
   clustSetSTPtr->endUI = 0;
   clustSetSTPtr->clustSI = 0;
} /*blank_set_clustST*/

/*-------------------------------------------------------\
| Fun02: init_set_clustST
|   - sets all settings in a set_clustST struct to
|     defaults
| Input:
|   - clustSetSTPtr;
|     o pointer to set_clustST struct to initialize
| Output:
|   - Modifies:
|     o all variables in clustSetSTPtr to be defaults
\-------------------------------------------------------*/
void
init_set_clustST(
   struct set_clustST *clustSetSTPtr
){
   clustSetSTPtr->repIntervalSL = def_repInterval_clustST;

   /*unique to clustST*/
   clustSetSTPtr->minProbF = def_minProb_clustST;
   clustSetSTPtr->minVarUI = def_minVar_clustST;

   /*unique to clustST*/
   clustSetSTPtr->conRebuildUC = def_conRebuild_clustST;
   clustSetSTPtr->depthProfBl = def_depthProfile_clustST;
   clustSetSTPtr->readErrRateF = def_readErrRate_clustST;
   clustSetSTPtr->conErrRateF = def_conErrRate_clustST;
   clustSetSTPtr->varToErrUI = def_maxVarErrRatio_clustST;
   clustSetSTPtr->lenWeightF = def_lenWeight_clustST;

   /*initial read filtering settings*/
   clustSetSTPtr->minLenUI = def_minLen_clustST;
   clustSetSTPtr->minMapqUC = def_minMapq_clustST;
   clustSetSTPtr->minMedQUI = def_minMedQ_clustST;
   clustSetSTPtr->minAvgQUI = def_minAvgQ_clustST;

   clustSetSTPtr->minPercDepthF =def_minPercDepth_clustST;
   clustSetSTPtr->percOverlapF = def_percOverlap_clustST;

   clustSetSTPtr->winSizeUI = def_window_clustST;
   clustSetSTPtr->winErrUI = def_windowError_clustST;

   /*consensus filtering*/
   clustSetSTPtr->maxConSimF = def_maxConSim_clustST;
   clustSetSTPtr->maxNPercF = def_maxNPerc_clustST;

   /*variant filtering*/
   clustSetSTPtr->minDepthUI = def_minDepth_clustST;
   clustSetSTPtr->indelLenUI = def_indelLen_clustST;
   clustSetSTPtr->minSnpQUC = def_minSnpQ_clustST;

   blank_set_clustST(clustSetSTPtr);
} /*init_set_clustST*/

/*-------------------------------------------------------\
| Fun03: freeStack_set_clustST
|   - no internal variables to free (here for support)
| Input:
|   - clustSetSTPtr;
|     o pointer to set_clustST struct with variables to
|       free
| Output:
|   - Modifies:
|     o initializes all variables (nothing to free)
\-------------------------------------------------------*/
void
freeStack_set_clustST(
   struct set_clustST *clustSetSTPtr
){
   if(clustSetSTPtr)
      init_set_clustST(clustSetSTPtr);
} /*freeStack_set_clustST*/

/*-------------------------------------------------------\
| Fun04: freeHeap_set_clustST
|   - frees a heap allocated set_clustST struct
| Input:
|   - clustSetSTPtr;
|     o pointer to set_clustST struct to free
| Output:
|   - Frees:
|     o clustSetSTPtr (you must set to 0)
\-------------------------------------------------------*/
void
freeHeap_set_clustST(
   struct set_clustST *clustSetSTPtr
){
   if(clustSetSTPtr)
   { /*If: have someting to free*/
      freeStack_set_clustST(clustSetSTPtr);
      free(clustSetSTPtr);
   } /*If: have someting to free*/
} /*freeHeap_set_clustST*/

/*-------------------------------------------------------\
| Fun05: blank_index_clustST
|   - blanks an index_clustST struct
| Input:
|   - indexSTPtr:
|     o pointer to index_clustST struct to blank
| Output:
|   - Modifies:
|     o lenSamUL, keptUL, and numRefUI in struct be 0
\-------------------------------------------------------*/
void
blank_index_clustST(
   struct index_clustST *indexSTPtr
){
   indexSTPtr->lenSamUL = 0;
   indexSTPtr->keptSL = 0;
   indexSTPtr->numRefUI = 0;
} /*blank_index_clustST*/

/*-------------------------------------------------------\
| Fun06: init_index_clustST
|   - initializes an index_clustST struct
| Input:
|   - indexSTPtr:
|     o pointer to index_clustST struct to initializes
| Output:
|   - Modifies:
|     o all arrays to be null and lengths to be 0
\-------------------------------------------------------*/
void
init_index_clustST(
   struct index_clustST *indexSTPtr
){
   indexSTPtr->scoreArySC = 0;
   indexSTPtr->startAryUI = 0;
   indexSTPtr->endAryUI = 0;
   indexSTPtr->refAryUI = 0;
   indexSTPtr->clustArySI = 0;

   indexSTPtr->indexAryUL = 0;
   indexSTPtr->lenLineAryUI = 0;

   indexSTPtr->lenUL = 0;

   indexSTPtr->refIdAryStr = 0;
   indexSTPtr->refNumAryUI = 0;
   indexSTPtr->lenRefUI = 0;

   blank_index_clustST(indexSTPtr);
} /*init_index_clustST*/

/*-------------------------------------------------------\
| Fun07: freeStack_index_clustST
|   - frees and initializes variables in index_clustST
| Input:
|   - indexSTPtr:
|     o pointer to index_clustST struct with variables to
|       free
| Output:
|   - Frees:
|     o all arrays, then initializes
\-------------------------------------------------------*/
void
freeStack_index_clustST(
   struct index_clustST *indexSTPtr
){
   if(indexSTPtr)
   { /*If: have structure*/
      if(indexSTPtr->scoreArySC)
         free(indexSTPtr->scoreArySC);

      if(indexSTPtr->startAryUI)
         free(indexSTPtr->startAryUI);

      if(indexSTPtr->endAryUI)
         free(indexSTPtr->endAryUI);

      if(indexSTPtr->refAryUI)
         free(indexSTPtr->refAryUI);

      if(indexSTPtr->clustArySI)
         free(indexSTPtr->clustArySI);


      if(indexSTPtr->indexAryUL)
         free(indexSTPtr->indexAryUL);
      indexSTPtr->indexAryUL = 0;

      if(indexSTPtr->lenLineAryUI)
         free(indexSTPtr->lenLineAryUI);


      if(indexSTPtr->refIdAryStr)
         free(indexSTPtr->refIdAryStr);

      if(indexSTPtr->refNumAryUI)
         free(indexSTPtr->refNumAryUI);

      init_index_clustST(indexSTPtr);
   } /*If: have structure*/
} /*freeStack_index_clustST*/

/*-------------------------------------------------------\
| Fun09: freeHeap_index_clustST
|   - frees index_clustST struct
| Input:
|   - indexSTPtr:
|     o pointer to index_clustST struct to free
| Output:
|   - Frees:
|     o indexSTPtr (you must set to 0)
\-------------------------------------------------------*/
void
freeHeap_index_clustST(
   struct index_clustST *indexSTPtr
){
   if(indexSTPtr)
   { /*If: have struct to free*/
      freeStack_index_clustST(indexSTPtr);
      free(indexSTPtr);
   } /*If: have struct to free*/
} /*freeHeap_index_clustST*/

/*-------------------------------------------------------\
| Fun10: setup_index_clustST
|   - setsup memory for index_clustST struct
| Input:
|   - indexSTPtr:
|     o pointer to index_clustST struct to free
|   - sizeUL:
|     o new size of arrays
| Output:
|   - Modifies:
|     o arrays in index_clustST to have memory allocated
|   - Returns:
|     o 0 for no errors
|     o def_memErr_clustST for memory errors
\-------------------------------------------------------*/
signed char
setup_index_clustST(
   struct index_clustST *indexSTPtr,
   unsigned long sizeUL
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun10 TOC:
   '   o fun10 sec01:
   '     - check if need to set default array size
   '   o fun10 sec02:
   '     - scoring and read mapping arrays
   '   o fun10 sec03:
   '     - clustering assignment array
   '   o fun10 sec04:
   '     - file coordinates arrays
   '   o fun10 sec05:
   '     - reference id array
   '   o fun10 sec06:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec01:
   ^   - check if need to set default array size
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(sizeUL == 0)
      sizeUL = 4096;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec02:
   ^   - scoring and read mapping arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   indexSTPtr->lenUL = sizeUL;

   if(indexSTPtr->scoreArySC)
      free(indexSTPtr->scoreArySC);
   indexSTPtr->scoreArySC = 0;

   indexSTPtr->scoreArySC =
      malloc(
         sizeUL *
         sizeof(schar)
      );

   if(! indexSTPtr->scoreArySC)
      goto memErr_fun10_sec06;


   if(indexSTPtr->startAryUI)
      free(indexSTPtr->startAryUI);
   indexSTPtr->startAryUI = 0;

   indexSTPtr->startAryUI =
      malloc(
         sizeUL *
         sizeof(uint)
      );

   if(! indexSTPtr->startAryUI)
      goto memErr_fun10_sec06;


   if(indexSTPtr->refAryUI)
      free(indexSTPtr->refAryUI);
   indexSTPtr->refAryUI = 0;

   indexSTPtr->refAryUI =
      malloc(
         sizeUL *
         sizeof(uint)
      );

   if(! indexSTPtr->refAryUI)
      goto memErr_fun10_sec06;


   if(indexSTPtr->endAryUI)
      free(indexSTPtr->endAryUI);
   indexSTPtr->endAryUI = 0;

   indexSTPtr->endAryUI =
      malloc(
         sizeUL *
         sizeof(uint)
      );

   if(! indexSTPtr->endAryUI)
      goto memErr_fun10_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec03:
   ^   - clustering assignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(indexSTPtr->clustArySI)
      free(indexSTPtr->clustArySI);
   indexSTPtr->clustArySI = 0;

   indexSTPtr->clustArySI =
      malloc(
         sizeUL *
         sizeof(uint)
      );

   if(! indexSTPtr->clustArySI)
      goto memErr_fun10_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec04:
   ^   - file coordinates arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(indexSTPtr->indexAryUL)
      free(indexSTPtr->indexAryUL);
   indexSTPtr->indexAryUL = 0;

   indexSTPtr->indexAryUL =
      malloc(
         sizeUL *
         sizeof(ulong)
      );

   if(! indexSTPtr->indexAryUL)
      goto memErr_fun10_sec06;


   if(indexSTPtr->lenLineAryUI)
      free(indexSTPtr->lenLineAryUI);
   indexSTPtr->lenLineAryUI = 0;

   indexSTPtr->lenLineAryUI =
      malloc(
         sizeUL *
         sizeof(uint)
      );

   if(! indexSTPtr->lenLineAryUI)
      goto memErr_fun10_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec05:
   ^   - reference id array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(indexSTPtr->refIdAryStr)
      free(indexSTPtr->refIdAryStr);
   indexSTPtr->refIdAryStr = 0;

   indexSTPtr->refIdAryStr = mk_strAry(16);

   if(! indexSTPtr->refIdAryStr)
      goto memErr_fun10_sec06;


   if(indexSTPtr->refNumAryUI)
      free(indexSTPtr->refNumAryUI);
   indexSTPtr->refNumAryUI = 0;

   indexSTPtr->refNumAryUI = malloc(16 * sizeof(uint));

   if(! indexSTPtr->refNumAryUI)
      goto memErr_fun10_sec06;


   indexSTPtr->lenRefUI = 16;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec06:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun10_sec06:;
   return def_memErr_clustST;
} /*setup_index_clustST*/

/*-------------------------------------------------------\
| Fun11: realloc_index_clustST
|   - reallocates memory for an index_clustST struct
| Input:
|   - indexSTPtr:
|     o pointer to index_clustST struct to assign more
|       memory to
|   - sizeUL:
|     o new size of arrays; use indexSTPtr->lenUL to not
|       resize arrays (except refIdAryStr)
|   - numRefUI:
|     o new number of references; use indexSTPtr->lenRefUI
|       to not resize reference array
| Output:
|   - Modifies:
|     o arrays in index_clustST to be size of sizeUL
|   - Returns:
|     o 0 for no errors
|     o def_memErr_clustST for memory errors
\-------------------------------------------------------*/
signed char
realloc_index_clustST(
   struct index_clustST *indexSTPtr,
   unsigned long sizeUL,
   unsigned int numRefUI
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC:
   '   o fun11 sec01:
   '     - variable declarations
   '   o fun11 sec02:
   '     - check if need to resize + resize read stats
   '   o fun11 sec03:
   '     - resize cluster assignment array
   '   o fun11 sec04:
   '     - resize file coordinate array
   '   o fun11 sec05:
   '     - check if need to resize reference id array
   '   o fun11 sec06:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned int *tmpUIPtr = 0;
   signed char *tmpStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec02:
   ^   - check if need to resize + resize read stats
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(sizeUL != indexSTPtr->lenUL)
   { /*If: need to resize arrays*/
      tmpUIPtr =
         (uint *)
         realloc(
            indexSTPtr->scoreArySC,
            sizeUL * sizeof(schar)
         );

      if(! tmpUIPtr)
         goto memErr_fun11_sec06;

      indexSTPtr->scoreArySC = (schar *) tmpUIPtr;



      tmpUIPtr =
         realloc(
            indexSTPtr->startAryUI,
            sizeUL * sizeof(uint)
         );

      if(! tmpUIPtr)
         goto memErr_fun11_sec06;

      indexSTPtr->startAryUI = tmpUIPtr;


      tmpUIPtr =
         realloc(
            indexSTPtr->endAryUI,
            sizeUL * sizeof(uint)
         );

      if(! tmpUIPtr)
         goto memErr_fun11_sec06;

      indexSTPtr->endAryUI = tmpUIPtr;


      tmpUIPtr =
         realloc(
            indexSTPtr->refAryUI,
            sizeUL * sizeof(uint)
         );

      if(! tmpUIPtr)
         goto memErr_fun11_sec06;

      indexSTPtr->refAryUI = tmpUIPtr;

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun11 Sec03:
      ^   - resize cluster assignment array
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      tmpUIPtr =
         (uint *)
         realloc(
            indexSTPtr->clustArySI,
            sizeUL * sizeof(sint)
         );

      if(! tmpUIPtr)
         goto memErr_fun11_sec06;

      indexSTPtr->clustArySI = (sint *) tmpUIPtr;


      indexSTPtr->lenUL = sizeUL;

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun11 Sec04:
      ^   - resize file coordinate array
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      tmpUIPtr =
         (uint *)
         realloc(
            indexSTPtr->indexAryUL,
            sizeUL * sizeof(ulong)
         );

      if(! tmpUIPtr)
         goto memErr_fun11_sec06;

      indexSTPtr->indexAryUL = (ulong *) tmpUIPtr;


      tmpUIPtr =
         realloc(
            indexSTPtr->lenLineAryUI,
            sizeUL * sizeof(uint)
         );

      if(! tmpUIPtr)
         goto memErr_fun11_sec06;

      indexSTPtr->lenLineAryUI = tmpUIPtr;

   } /*If: need to resize arrays*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec05:
   ^   - check if need to resize reference id array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numRefUI != indexSTPtr->lenRefUI)
   { /*If: resizing reference array*/
      tmpStr =
         realloc_strAry(
            indexSTPtr->refIdAryStr,
            numRefUI
         );

      if(! tmpStr)
         goto memErr_fun11_sec06;

      indexSTPtr->refIdAryStr = tmpStr;


      tmpUIPtr =
         realloc(
            indexSTPtr->refNumAryUI,
            numRefUI
         );

      if(! tmpUIPtr)
         goto memErr_fun11_sec06;

      indexSTPtr->refNumAryUI = tmpUIPtr;

      indexSTPtr->lenRefUI = numRefUI;
   } /*If: resizing reference array*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec06:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun11_sec06:;
   return def_memErr_clustST;
} /*realloc_index_clustST*/

/*-------------------------------------------------------\
| Fun12: blank_con_clustST
|   - blank an clustST con_clustST struct
| Input:
|   - conSTPtr:
|     o pointer to con_clustST struct to blank
| Output:
|   - Modifies:
|     o conSTPtr to be blank 
\-------------------------------------------------------*/
void
blank_con_clustST(
   struct con_clustST *conSTPtr
){
   if(conSTPtr->samSTPtr)
      blank_samEntry(conSTPtr->samSTPtr);

   conSTPtr->clustSI = 0;
   conSTPtr->numReadsUL = 0;
   conSTPtr->maxReadsUL = 0;
   conSTPtr->numVarUI = 0;
   conSTPtr->startUI = 0;
   conSTPtr->endUI = 0;
} /*blank_con_clustST*/

/*-------------------------------------------------------\
| Fun13: init_con_clustST
|   - initialzie an clustST con_clustST struct
| Input:
|   - conSTPtr:
|     o pointer to con_clustST struct to initialize
| Output:
|   - Modifies:
|     o conSTPtr to intialize
\-------------------------------------------------------*/
void
init_con_clustST(
   struct con_clustST *conSTPtr
){
   conSTPtr->samSTPtr = 0;
   conSTPtr->nextST = 0;

   conSTPtr->typeArySC = 0;
   conSTPtr->posAryUI = 0;
   conSTPtr->lenProfUI = 0;

   blank_con_clustST(conSTPtr);
} /*init_con_clustST*/

/*-------------------------------------------------------\
| Fun14: freeStack_con_clustST
|   - frees variables inside an con_clustST struct
| Input:
|   - conSTPtr:
|     o pointer to con_clustST struct to free varaibles
| Output:
|   - Frees:
|     o samSTPtr in conSTPtr (nextST not freeded)
|   - Modifies:
|     o blanks conSTPtr
\-------------------------------------------------------*/
void
freeStack_con_clustST(
   struct con_clustST *conSTPtr
){
   if(conSTPtr)
   { /*If: have a struct*/
      if(conSTPtr->samSTPtr)
         freeHeap_samEntry(conSTPtr->samSTPtr);

      if(conSTPtr->typeArySC)
         free(conSTPtr->typeArySC);

      if(conSTPtr->posAryUI)
         free(conSTPtr->posAryUI);

      init_con_clustST(conSTPtr);
   } /*If: have a struct*/
} /*freeStack_con_clustST*/

/*-------------------------------------------------------\
| Fun15: freeHeap_con_clustST
|   - frees a con_clustST struct
| Input:
|   - conSTPtr:
|     o pointer to con_clustST struct to free
| Output:
|   - Frees:
|     o conSTPtr
|   - Returns:
|     o next con_clustST struct in list
\-------------------------------------------------------*/
con_clustST *
freeHeap_con_clustST(
   struct con_clustST *conSTPtr
){
   struct con_clustST *nextSTPtr = 0;

   if(conSTPtr)
   { /*If: have a struct*/
      nextSTPtr = conSTPtr->nextST;

      freeStack_con_clustST(conSTPtr);
      free(conSTPtr);
   } /*If: have a struct*/

   return nextSTPtr;
} /*freeHeap_con_clustST*/

/*-------------------------------------------------------\
| Fun16: freeHeapList_con_clustST
|   - frees a list of heap allocated con_clustST structs
| Input:
|   - conSTPtr:
|     o pointer to con_clustST struct list to free
| Output:
|   - Frees:
|     o conSTPtr and all other structs in the list
\-------------------------------------------------------*/
void
freeHeapList_con_clustST(
   struct con_clustST *conSTPtr
){
   while(conSTPtr)
      conSTPtr = freeHeap_con_clustST(conSTPtr);
} /*freeHeapList_con_clustST*/

/*-------------------------------------------------------\
| Fun17: mk_con_clustST
|   - makes a con_clustST struct on the heap
| Input:
|   - samSTPtr:
|     o pointer to sam entry struct to store (0 for none)
|     o only a shallow copy is done (pointer saved only)
|   - clustSI:
|     o cluster number to assign
|   - numReadsUL:
|     o number of reads in cluster
| Output:
|   - Returns:
|     o pointer con_clustST struct with input
|     o 0 for memory error
\-------------------------------------------------------*/
con_clustST *
mk_con_clustST(
   struct samEntry *samSTPtr, /*consensus to add*/
   signed int clustSI,        /*cluster number*/
   unsigned long numReadsUL   /*number reads in cluster*/
){
   struct con_clustST *retSTPtr = 0;

   retSTPtr = malloc(sizeof(struct con_clustST));

   if(! retSTPtr)
      goto memErr_fun17_sec0x;

   init_con_clustST(retSTPtr);

   retSTPtr->samSTPtr = samSTPtr;
   retSTPtr->clustSI = clustSI;
   retSTPtr->numReadsUL = numReadsUL;

   return retSTPtr;

   memErr_fun17_sec0x:;
   if(retSTPtr)
      freeHeap_con_clustST(retSTPtr);
   retSTPtr = 0;
   return 0;
} /*mk_con_clustST*/

/*-------------------------------------------------------\
| Fun18: swap_con_clustST
|   - swaps to con_clustST structs (except nextST ptr)
| Input:
|   - firstSTPtr:
|     o first struct to swap
|   - secSTPtr:
|     o second struct to swap
| Output:
|   - Modifies:
|     o firstSTPtr to have secSTPtr variables
|     o secSTPtr to have frstSTPtr variables
|   - Note: does not swap nextST pointer
\-------------------------------------------------------*/
void
swap_con_clustST(
   struct con_clustST *firstSTPtr,
   struct con_clustST *secSTPtr
){
   struct samEntry *samSwapSTPtr = 0;
   signed char *swapSCPtr = 0;
   
   samSwapSTPtr = firstSTPtr->samSTPtr;
   firstSTPtr->samSTPtr = secSTPtr->samSTPtr;
   secSTPtr->samSTPtr = samSwapSTPtr;

   firstSTPtr->clustSI ^= secSTPtr->clustSI;
   secSTPtr->clustSI ^= firstSTPtr->clustSI;
   firstSTPtr->clustSI ^= secSTPtr->clustSI;

   firstSTPtr->numReadsUL ^= secSTPtr->numReadsUL;
   secSTPtr->numReadsUL ^= firstSTPtr->numReadsUL;
   firstSTPtr->numReadsUL ^= secSTPtr->numReadsUL;

   firstSTPtr->startUI ^= secSTPtr->startUI;
   secSTPtr->startUI ^= firstSTPtr->startUI;
   firstSTPtr->startUI ^= secSTPtr->startUI;

   firstSTPtr->endUI ^= secSTPtr->endUI;
   secSTPtr->endUI ^= firstSTPtr->endUI;
   firstSTPtr->endUI ^= secSTPtr->endUI;

   /*swap the profile*/

   swapSCPtr = firstSTPtr->typeArySC;
   firstSTPtr->typeArySC = secSTPtr->typeArySC;
   secSTPtr->typeArySC = swapSCPtr;

   swapSCPtr = (schar *) firstSTPtr->posAryUI;
   firstSTPtr->posAryUI = secSTPtr->posAryUI;
   secSTPtr->posAryUI = (uint *) swapSCPtr;

   firstSTPtr->numVarUI ^= secSTPtr->numVarUI;
   secSTPtr->numVarUI ^= firstSTPtr->numVarUI;
   firstSTPtr->numVarUI ^= secSTPtr->numVarUI;

   firstSTPtr->lenProfUI ^= secSTPtr->lenProfUI;
   secSTPtr->lenProfUI ^= firstSTPtr->lenProfUI;
   firstSTPtr->lenProfUI ^= secSTPtr->lenProfUI;
} /*swap_con_clustST*/

/*-------------------------------------------------------\
| Fun19: getNumLines_clustST
|   - finds number of lines in a file
| Input:
|   - inFILE:
|     o file to find number of lines in
| Output:
|   - Modifies:
|     o inFILE to point to start of file
|   - Returns:
|     o number of lines in file
\-------------------------------------------------------*/
unsigned long
getNumLines_clustST(
   void *inFILE
){
   #define len_fun18 4096
   schar buffStr[len_fun18 + 1];
   schar *tmpStr = 0;
   ulong bytesUL = 0;
   ulong lineUL = 0;

   bytesUL =
      fread(
         buffStr,
         sizeof(schar),
         len_fun18,
         inFILE
      ); /*read first line*/

   buffStr[bytesUL] = '\0';
   tmpStr = buffStr;

   while(
         bytesUL
      || *tmpStr != '\0'
   ){ /*Loop: find number of new lines*/
      if(*tmpStr == '\n')
      { /*If: have new line*/
         ++lineUL;
         ++tmpStr;
      } /*If: have new line*/

      else if(*tmpStr == '\0')
      { /*Else: read in more file*/
         bytesUL =
            fread(
               buffStr,
               sizeof(schar),
               len_fun18,
               inFILE
            ); /*read first line*/

         if(bytesUL > 0)
            buffStr[bytesUL] = '\0';

         tmpStr = buffStr;
      } /*Else: read in more file*/

      tmpStr += endLineUnix_ulCp(tmpStr);
         /*in this case do not care about '\r'*/
   } /*Loop: find number of new lines*/

   if(*(tmpStr - 1) != '\n')
      ++lineUL; /*account for last line being '\0'*/

   fseek(
      inFILE,
      0,
      SEEK_SET
   ); /*move to start of file*/

   return lineUL;
} /*getNumLines_clustST*/

/*-------------------------------------------------------\
| Fun20: mk_index_clustST
|   - makes an index_clustST struct for a sam file
| Input:
|   - sizeUL:
|     o size to make arrays at start
|     o any value > 0 will do since arrays are resized
|       when needed, however, for best memory usage, input
|       number of lines in sam file
|   - clustSetSTPtr:
|     o settings for filtering reads
|   - samSTPtr:
|     o pointer to sam struct to use in read hunt
|   - buffStrPtr:
|     o pointer to c-string to use reading sam file
|   - lenBuffULPtr:
|     o pointer to unsigned long with size of buffStrPtr
|   - samFILE:
|     o sam file to convert to stats
| Output:
|   - Modifies:
|     o samSTPtr to hold last entry
|     o buffStrPtr to be larger if needed
|     o lenBuffULPtr to be new buffStrPtr size (if resized)
|     o samFILE to point to start
|   - Returns:
|     o 0 for memory error
|     o pointer to index_clustST struct with scores and
|       index's (index_clustST->keptSL is 0 if no reads)
|       - in clustArySI:
|         i. open reads set to 0
|         ii. headers set to def_header_clustST
|         iii. discarded reads set to def_discard_clustST
\-------------------------------------------------------*/
struct index_clustST *
mk_index_clustST(
   unsigned long sizeUL,
   struct set_clustST *clustSetSTPtr,
   struct samEntry *samSTPtr,
   signed char **buffStrPtr,
   unsigned long *lenBuffUL,
   void *samFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun20 TOC:
   '   - get read scores for a sam file
   '   o fun20 sec01:
   '     - variable declarations
   '   o fun20 sec02:
   '     - setup; move to file start and memory allocation
   '   o fun20 sec03:
   '     - get scores
   '   o fun20 sec04:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ulong lineUL = 0;
   schar errSC = 0;
   uchar scoreUC = 0;  /*used in scoring step*/
   uint lenUI = 0;     /*used in scoring (log10) step*/

   slong filePosSL = 0; /*position at in file*/

   slong indexSL = 0;        /*find reads reference*/
   slong slCnt = 0;          /*shifting indexs*/

   struct index_clustST *retHeapST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec02:
   ^   - setup (move to file start and memory allocation)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fseek(
      samFILE,
      0,
      SEEK_SET
   ); /*make sure at start of file*/

   /*allocate memory for scores, index's, and clusters*/
   retHeapST = malloc(sizeof(struct index_clustST));

   if(! retHeapST)
      goto memErr_fun20_sec04;
   
   init_index_clustST(retHeapST);

   if(sizeUL == 0)
      sizeUL = 4096;

   errSC =
      setup_index_clustST(
         retHeapST,
         sizeUL
      );

   if(errSC)
      goto memErr_fun20_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec03:
   ^   - get scores
   ^   o fun20 sec03 sub01:
   ^     - get first line of sam file
   ^   o fun20 sec03 sub02:
   ^     - start loop & handel memory resizing (if needed)
   ^   o fun20 sec03 sub03:
   ^     - get score or mark as discarded
   ^   o fun20 sec03 sub04:
   ^     - get next entry
   ^   o fun20 sec03 sub05:
   ^     - final error check
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun20 Sec03 Sub01:
   *   - get first line of sam file
   \*****************************************************/

   errSC =
      get_samEntry(
         samSTPtr,
         buffStrPtr,
         lenBuffUL,
         samFILE
      );

   if(errSC == def_memErr_samEntry)
      goto memErr_fun20_sec04;

   else if(errSC)
      goto noReads_fun20_sec04;

   /*****************************************************\
   * Fun20 Sec03 Sub02:
   *   - start loop and handel memory resizing (if needed)
   \*****************************************************/

   while(! errSC)
   { /*Loop: print out stats*/
      if(retHeapST->lenSamUL >= retHeapST->lenUL)
      { /*If: need to get more memory*/
         errSC =
            realloc_index_clustST(
               retHeapST,
               retHeapST->lenUL << 1,
               retHeapST->lenRefUI    /*not resized*/
            );

         if(errSC)
            goto memErr_fun20_sec04;
      } /*If: need to get more memory*/

      /**************************************************\
      * Fun20 Sec03 Sub03:
      *   - detect entires to discard
      \**************************************************/

      if(samSTPtr->extraStr[0] == '@')
      { /*If: comment*/
         retHeapST->scoreArySC[lineUL] =
            def_header_clustST;

         retHeapST->clustArySI[lineUL] =
            def_header_clustST;

         retHeapST->indexAryUL[lineUL] =(ulong) filePosSL;

         retHeapST->startAryUI[lineUL] = 0;
         retHeapST->endAryUI[lineUL] = 0;
         retHeapST->refAryUI[lineUL] = 0;
      } /*If: comment*/

      else if(samSTPtr->flagUS & (4 | 256 | 2048))
      { /*Else If: unmapped, secondary, or supplemental*/
         retHeapST->scoreArySC[lineUL] =
            def_discard_clustST;

         retHeapST->clustArySI[lineUL] =
            def_discard_clustST;

         retHeapST->indexAryUL[lineUL] =(ulong) filePosSL;

         retHeapST->startAryUI[lineUL] = 0;
         retHeapST->endAryUI[lineUL] = 0;
         retHeapST->refAryUI[lineUL] = 0;
      } /*Else If: unmapped, secondary, or supplemental*/

      else if(samSTPtr->mapqUC < clustSetSTPtr->minMapqUC)
      { /*Else If: mapping quality to low*/
         retHeapST->scoreArySC[lineUL] =
            def_discard_clustST;

         retHeapST->clustArySI[lineUL] =
            def_discard_clustST;

         retHeapST->indexAryUL[lineUL] =(ulong) filePosSL;

         retHeapST->startAryUI[lineUL] = 0;
         retHeapST->endAryUI[lineUL] = 0;
         retHeapST->refAryUI[lineUL] = 0;
      } /*Else If: mapping quality to low*/

      else if(
         samSTPtr->alnReadLenUI < clustSetSTPtr->minLenUI
      ){ /*Else If: aligned length to short*/
         retHeapST->scoreArySC[lineUL] =
            def_discard_clustST;

         retHeapST->clustArySI[lineUL] =
            def_discard_clustST;

         retHeapST->indexAryUL[lineUL] =(ulong) filePosSL;

         retHeapST->startAryUI[lineUL] = 0;
         retHeapST->endAryUI[lineUL] = 0;
         retHeapST->refAryUI[lineUL] = 0;
      } /*Else If: aligned length to short*/

      else if(
         samSTPtr->medianQF < clustSetSTPtr->minMedQUI
      ){ /*Else If: median q-score to low*/
         retHeapST->scoreArySC[lineUL] =
            def_discard_clustST;

         retHeapST->clustArySI[lineUL] =
            def_discard_clustST;

         retHeapST->refAryUI[lineUL] =
            def_discard_clustST;

         retHeapST->indexAryUL[lineUL] =(ulong) filePosSL;

         retHeapST->startAryUI[lineUL] = 0;
         retHeapST->endAryUI[lineUL] = 0;
         retHeapST->refAryUI[lineUL] = 0;
      } /*Else If: median q-score to low*/

      else if(samSTPtr->meanQF < clustSetSTPtr->minAvgQUI)
      { /*Else If: mean q-score to low*/
         retHeapST->scoreArySC[lineUL] =
            def_discard_clustST;

         retHeapST->clustArySI[lineUL] =
            def_discard_clustST;

         retHeapST->indexAryUL[lineUL] =(ulong) filePosSL;

         retHeapST->startAryUI[lineUL] = 0;
         retHeapST->endAryUI[lineUL] = 0;
         retHeapST->refAryUI[lineUL] = 0;
      } /*Else If: mean q-score to low*/

      /**************************************************\
      * Fun20 Sec03 Sub04:
      *   - score and index reads
      *   o fun20 sec03 sub04 cat01:
      *     - assign ref number; check if new ref
      *   o fun20 sec03 sub04 cat02:
      *     - assign ref number; add new ref
      *   o fun20 sec03 sub04 cat03:
      *     - get/assign score and index
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun20 Sec03 Sub04 Cat01:
      +   - assign ref number; check if new ref
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else
      { /*Else: keep read*/
         ++retHeapST->keptSL;

         if(retHeapST->numRefUI == 0)
         { /*If: first reference*/
            add_strAry(
               samSTPtr->refIdStr,
               retHeapST->refIdAryStr,
               0
            );

            retHeapST->refNumAryUI[0] = 1;
            retHeapST->numRefUI = 1;

            retHeapST->refAryUI[lineUL] =
               retHeapST->numRefUI;
         } /*If: first reference*/

         else
         { /*Else: get reference*/
            indexSL =
               find_strAry(
                  retHeapST->refIdAryStr,
                  samSTPtr->refIdStr,
                  retHeapST->numRefUI
               );

            if(indexSL >= 0)
               retHeapST->refAryUI[lineUL] = 
                  retHeapST->refNumAryUI[indexSL];
                  /*assigned reference number*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun20 Sec03 Sub04 Cat02:
            +   - assign ref number; add new ref
            \+++++++++++++++++++++++++++++++++++++++++++*/

            else
            { /*Else: is a new reference*/
               if(
                     retHeapST->numRefUI + 1
                  >= retHeapST->lenRefUI
               ){ /*If: need more memory*/

                  errSC =
                     realloc_index_clustST(
                        retHeapST,
                        retHeapST->lenUL, /*not resizing*/
                        (retHeapST->numRefUI + 1) << 1
                     );
                     /*we have two arrays, the reference
                     `   id array and the read array, we
                     `   only need to resize reference id
                     */

                  if(errSC)
                     goto memErr_fun20_sec04;

               } /*If: need more memory*/

               indexSL =
                  (signed long)
                  addSort_strAry(
                     samSTPtr->refIdStr,
                     retHeapST->refIdAryStr,
                     retHeapST->numRefUI
                  );

               /*need to account for inserted index*/
               for(
                  slCnt =(slong) retHeapST->numRefUI;
                  slCnt > indexSL;
                  --slCnt
               ) retHeapST->refNumAryUI[slCnt] =
                    retHeapST->refNumAryUI[slCnt - 1];

               ++retHeapST->numRefUI;

               retHeapST->refNumAryUI[indexSL] =
                  retHeapST->numRefUI;

               retHeapST->refAryUI[lineUL] =
                  retHeapST->numRefUI;
            } /*Else: is a new reference*/
         } /*Else: get reference*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun20 Sec03 Sub04 Cat03:
         +   - get/assign score and index
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         lenUI = samSTPtr->alnReadLenUI;

         scoreUC = 0;

         while(lenUI > 0)
         { /*Loop: log10(length)*/
            ++scoreUC;
            lenUI /= 10;
         } /*Loop: log10(length)*/

         scoreUC *= clustSetSTPtr->lenWeightF;

         retHeapST->scoreArySC[lineUL] = 
            scoreUC + (uchar) samSTPtr->medianQF;

         retHeapST->startAryUI[lineUL] =
            samSTPtr->refStartUI;

         retHeapST->endAryUI[lineUL] =
            samSTPtr->refEndUI;

         retHeapST->clustArySI[lineUL] = 0;
         retHeapST->indexAryUL[lineUL] =(ulong) filePosSL;
      } /*Else: keep read*/
 
      /**************************************************\
      * Fun20 Sec03 Sub04:
      *   - get next entry
      \**************************************************/

      /*moving to end of line, so want to skip '\r'*/
      retHeapST->lenLineAryUI[lineUL] =
         endLineUnix_ulCp(*buffStrPtr);

      ++retHeapST->lenLineAryUI[lineUL]; /*for '\n'*/

      filePosSL += retHeapST->lenLineAryUI[lineUL];

      errSC =
         get_samEntry(
            samSTPtr,
            buffStrPtr,
            lenBuffUL,
            samFILE
         );

      ++lineUL;
   } /*Loop: print out stats*/

   /*****************************************************\
   * Fun20 Sec03 Sub05:
   *   - final error check
   \*****************************************************/

   retHeapST->lenSamUL = lineUL;

   if(errSC == def_memErr_samEntry)
      goto memErr_fun20_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun20_sec04;

   memErr_fun20_sec04:;

   if(retHeapST)
      freeHeap_index_clustST(retHeapST);
   retHeapST = 0;

   goto ret_fun20_sec04;

   noReads_fun20_sec04:;
   retHeapST->keptSL = 0;
   goto ret_fun20_sec04;

   ret_fun20_sec04:;

   fseek(
      samFILE,
      0,
      SEEK_SET
   );

   return retHeapST;
} /*mk_index_clustST*/

/*-------------------------------------------------------\
| Fun21: getRead_clustST
|   - gets a read by index from a sam file
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to hold read
|   - lineUL:
|     o line number of read in sam file
|   - offsetUL:
|     o number of bytes from current position to read
|   - indexSTPtr:
|     o pointer to index_clustST struct with line lengths
|   - buffStr:
|     o c-string to read in line with
|     o should be length of target reads line + 2
|   - samFILE:
|     o FILE pointer to sam file to get read from
| Output:
|   - Modifies:
|     o samSTPtr to have new read
|     o buffStr to have reads line
|     o samFILE to be on read after current read
|   - Returns:
|     o 0 for no errors
|     o def_memErr_clustST for memory errors
|     o def_fileErr_clustST for file errors
\-------------------------------------------------------*/
signed char
getRead_clustST(
   struct samEntry *samSTPtr,       /*will have read*/
   unsigned long lineUL,            /*line of read*/
   unsigned long offsetUL,          /*offset to read*/
   struct index_clustST *indexSTPtr,/*index's of reads*/
   signed char *buffStr,            /*for getting line*/
   void *samFILE                    /*file with read*/
){
   signed char errSC = 0;

   if(offsetUL)
   { /*If: skipping entries*/
      fseek(
         samFILE,
         offsetUL,
         SEEK_CUR
      );

      offsetUL = 0;
   } /*If: skipping entries*/

   offsetUL =
      fread(
         buffStr,
         sizeof(schar),
         indexSTPtr->lenLineAryUI[lineUL],
         samFILE
      );

   if(! offsetUL)
      return def_fileErr_clustST;

   buffStr[offsetUL] = '\0';

   errSC =
      lineTo_samEntry(
         samSTPtr,
         buffStr
      );

   if(errSC)
      return def_memErr_clustST;

   return 0;
} /*getRead_clustST*/

/*-------------------------------------------------------\
| Fun22: getCon_clustST
|   - builds a consensus for a cluster
| Input:
|   - indexSTPtr:
|     o pointer to index_clustST struct with cluster array
|       to use
|     o clusters assigned with fun03 (0 is no cluster)
|   - clustSI:
|     o cluster number to build consensus for
|   - conSetSTPtr:
|     o pointer to set_tbCon struct with setting for
|       consensus building
|   - clustSetSTPtr:
|     o has settings for filtering and stores the first
|       and last possible base of the consensesus
|   - samSTPtr:
|     o pointer to samEntry struct for reading sam file
|   - buffStr:
|     o c-string to use reading sam file
|   - samFILE:
|     o sam file with reads to add to consensus
|   - errSCPtr:
|     o pointer to signed char to hold error type
| Output:
|   - Modifies:
|     o samSTPtr to be last read in file
|     o samFILE to be at start
|     o buffStrPtr to have last sam file entry
|     o errSCPtr:
|       - 0 for no errors
|       - def_memErr_clustST for any other error
|       - def_noCon_clustST if consensus had to many N's
|         or could not be built
|   - Returns:
|     o 0 if could not build consenensus or had error
|     o samEntry strucure array with consensus fragments
\-------------------------------------------------------*/
struct samEntry *
getCon_clustST(
   struct index_clustST *indexSTPtr,
   signed int clustSI,
   struct set_clustST *clustSetSTPtr,
   struct set_tbCon *conSetSTPtr,
   struct samEntry *samSTPtr,
   signed char *buffStr,
   void *samFILE,
   signed char *errSCPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun22 TOC:
   '   - finds number of reads mapping to best read
   '   o fun22 sec01:
   '     - variable declarations
   '   o fun22 sec02:
   '     - get first entry in sam file
   '   o fun22 sec03:
   '     - build consensus
   '   o fun22 sec04:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   slong lineSL = 0;

   struct samEntry *samHeapST = 0;
   struct conNt_tbCon *conHeapAryST = 0;
   uint lenConUI = 0;
   uint numMaskUI = 0; /*number of masked bases in con*/
   float percF = 0;

   ulong offsetUL = 0; /*offset to get target read*/

   schar refIdStr[64];

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec02:
   ^   - get first entry in sam file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   refIdStr[0] = '\0';

   lenConUI =
      clustSetSTPtr->endUI - clustSetSTPtr->startUI;

   ++lenConUI;
   lenConUI += 128; /*allow for some offset*/

   fseek(
      samFILE,
      0,
      SEEK_SET
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec03:
   ^   - build consensus
   ^   o fun22 sec03 sub01:
   ^     - filter reads + start loop
   ^   o fun22 sec03 sub02:
   ^     - get read
   ^   o fun22 sec03 sub03:
   ^     - add reads to consensus
   ^   o fun22 sec03 sub04:
   ^     - collapse consensus and check for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun22 Sec03 Sub01:
   *   - filter reads + start loop
   \*****************************************************/

   for(
      lineSL = 0;
      lineSL < (slong) indexSTPtr->lenSamUL;
      ++lineSL
   ){ /*Loop: add reads to consensus*/

      if(indexSTPtr->clustArySI[lineSL] != clustSI)
      { /*If: read not in cluster*/
         offsetUL += indexSTPtr->lenLineAryUI[lineSL];
         continue;
      } /*If: read not in cluster*/

      /**************************************************\
      * Fun22 Sec03 Sub02:
      *   - get read
      \**************************************************/

      errSC =
         getRead_clustST(
            samSTPtr,
            (ulong) lineSL,
            offsetUL,
            indexSTPtr,
            buffStr,
            samFILE
         ); /*get read from sam file*/

      offsetUL = 0;

      if(errSC)
         goto err_fun22_sec04;

      /**************************************************\
      * Fun22 Sec03 Sub03:
      *   - add reads to conensus
      \**************************************************/

      if(refIdStr[0] == '\0')
      { /*If: have not recorded reference name*/
         cpLen_ulCp(
            refIdStr,
            samSTPtr->refIdStr,
            samSTPtr->lenRefIdUC
         );
      } /*If: have not recorded reference name*/

      /*adjust to local coordiantes*/

      numMaskUI =
         min_genMath(
            samSTPtr->refStartUI,
            clustSetSTPtr->startUI
         );

      samSTPtr->refStartUI -= numMaskUI;
      samSTPtr->refEndUI -= numMaskUI;

      /*add read to consensus*/

      errSC =
         addRead_tbCon(
            samSTPtr,
            &conHeapAryST,
            &lenConUI,
            conSetSTPtr
         );

      if(errSC == def_memErr_tbConDefs)
         goto err_fun22_sec04;
   } /*Loop: add reads to consensus*/

   /*****************************************************\
   * Fun22 Sec03 Sub04:
   *   - collapse consensus and check for errors
   \*****************************************************/

   samHeapST =
      noFragCollapse_tbCon(
         conHeapAryST,
         lenConUI,
         refIdStr,
         &numMaskUI,
         conSetSTPtr,
         &errSC
      );

   if(errSC == def_memErr_tbConDefs)
      goto err_fun22_sec04;

   if(! samHeapST)
      goto noCon_fun22_sec04;

   percF = (float) numMaskUI;
   percF /= (float) samHeapST->readLenUI;

   if(percF > clustSetSTPtr->maxNPercF)
      goto noCon_fun22_sec04;

   /*set up global coordinates*/
   samHeapST->refStartUI += clustSetSTPtr->startUI;
   samHeapST->refEndUI += clustSetSTPtr->startUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec04:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *errSCPtr = 0;
   goto ret_fun22_sec04;


   noCon_fun22_sec04:;
   *errSCPtr = def_noCon_clustST;

   if(samHeapST)
      freeHeap_samEntry(samHeapST);
   samHeapST = 0;
   goto ret_fun22_sec04;


   err_fun22_sec04:;
   *errSCPtr = def_memErr_clustST;

   if(samHeapST)
      freeHeap_samEntry(samHeapST);
   samHeapST = 0;

   goto ret_fun22_sec04;



   ret_fun22_sec04:;

   if(conHeapAryST)
      freeHeapAry_conNt_tbCon(
         conHeapAryST,
         lenConUI
      );

   lenConUI = 0;
   conHeapAryST = 0;

   fseek(
      samFILE,
      0,
      SEEK_SET
   );

   return samHeapST;
} /*getCon_clustST*/

/*-------------------------------------------------------\
| Fun23: cmpCons_clustST
|   - compares two consensus clusters to see if same
| Input:
|   - firstConSTPtr:
|     o pointer to con_clustST struct with first consensus
|       to compare
|   - secConSTPtr:
|     o pointer to con_clustST struct with second
|       consensus to compare
|   - noBoundCheckBl:
|     o 1: skip bounds check (one read overlaps other)
|     o 0: do a bounds check
|   - indexSTPtr:
|     o pointer to index_clustST with clusters
|   - clustSetSTPtr:
|     o pointer to set_clustST struct with settings
|   - edDistResSTPtr:
|     o pointer to edDistResSTPtr struct to hold results
|       from edDist
| Output:
|   - Modifies:
|     o clustArySI in indexSTPtr to have clusters merged
|       if consensus are to similar (lowest cluster
|       assigned)
|     o merges second cluster into first cluster if
|       clusters were to similar
|       - you must free the second cluster (always has
|         discarded cluster)
|   - Returns:
|     o 0 if clusters are different
|     o 1 if merged clusters
\-------------------------------------------------------*/
signed char
cmpCons_clustST(
   struct con_clustST *firstConSTPtr,
   struct con_clustST *secConSTPtr,
   signed char noBoundCheckBl, /*1: no bounds check*/
   struct index_clustST *indexSTPtr,
   struct set_clustST *clustSetSTPtr,
   struct res_edDist *edDistResSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun23: cmpCons_clustST
   '   - compares two consensus clusters to see if same
   '   o fun23 sec01:
   '     - variable declarations   
   '   o fun23 sec02:
   '     - get edit distance and filter distance cons out
   '   o fun23 sec03:
   '     - merge similar consensus
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   slong distSL = 0;

   schar retSC = 0;

   /*changing cluster numbers*/
   sint keepClustSI = 0;
   sint changeClustSI = 0;

   float distPercF = 0;
   uint *tmpUIPtr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec02:
   ^   - get edit distance and filter distance cons out
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*disable depth profiling for consensuses*/
   tmpUIPtr = edDistResSTPtr->depthAryUI;
   edDistResSTPtr->depthAryUI = 0;

   if(! noBoundCheckBl)
   { /*If: not doing bounds check*/
      if(
              firstConSTPtr->samSTPtr->refStartUI 
            < secConSTPtr->samSTPtr->refStartUI
         &&
              firstConSTPtr->samSTPtr->refEndUI 
            < secConSTPtr->samSTPtr->refEndUI
      ) return 0; /*one consensus does not cover other*/

      if(
              firstConSTPtr->samSTPtr->refStartUI 
            > secConSTPtr->samSTPtr->refStartUI
         &&
              firstConSTPtr->samSTPtr->refEndUI 
            > secConSTPtr->samSTPtr->refEndUI
      ) return 0; /*one consensus does not cover other*/
   } /*If: not doing bounds check*/
   
   distSL =
      readCmpDist_edDist(
         firstConSTPtr->samSTPtr,
         secConSTPtr->samSTPtr,
         clustSetSTPtr->indelLenUI,
         clustSetSTPtr->minSnpQUC,
         -1,              /*will check overlaps later*/
         0,               /*not doing depth profile*/
         clustSetSTPtr->winSizeUI,
         edDistResSTPtr
      ); /*get edit distance*/

   edDistResSTPtr->depthAryUI = tmpUIPtr;

   distPercF = (float) distSL;
   distPercF /= (float) edDistResSTPtr->overlapUI;

   if(1 - distPercF <= clustSetSTPtr->maxConSimF)
      return 0; /*to much difference*/

   distSL =
      percDist_edDist(
        edDistResSTPtr,
        clustSetSTPtr->winSizeUI,
        clustSetSTPtr->conErrRateF
      );

   if(
           edDistResSTPtr->overlapUI
        >= clustSetSTPtr->winSizeUI
     &&
           edDistResSTPtr->probMaxWinUI
        > clustSetSTPtr->winErrUI
   ) return 0; /*window scan found difference region*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec03:
   ^   - merge similar consensus
   ^   o fun23 sec03 sub01:
   ^     - find lowest cluster number
   ^   o fun23 sec03 sub02:
   ^     - find consensus to return as keep
   ^   o fun23 sec03 sub03:
   ^     - merge clusters and return keep consensus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun23 Sec03 Sub01:
   *   - find lowest cluster number
   \*****************************************************/

   if(
        firstConSTPtr->clustSI
      < secConSTPtr->clustSI
   ){ /*If: first cluster has lowest number*/
      keepClustSI = firstConSTPtr->clustSI;
      changeClustSI = secConSTPtr->clustSI;

      secConSTPtr->clustSI = firstConSTPtr->clustSI;
   } /*If: first cluster has lowest number*/

   else
   { /*Else: second cluster has lowest number*/
      keepClustSI = secConSTPtr->clustSI;
      changeClustSI = firstConSTPtr->clustSI;

      firstConSTPtr->clustSI = secConSTPtr->clustSI;
   } /*Else: second cluster has lowest number*/

   /*****************************************************\
   * Fun23 Sec03 Sub02:
   *   - find consensus to return as keep
   \*****************************************************/

   firstConSTPtr->numReadsUL += secConSTPtr->numReadsUL;
   secConSTPtr->numReadsUL = firstConSTPtr->numReadsUL;

   if(
        firstConSTPtr->samSTPtr->alnReadLenUI
      > secConSTPtr->samSTPtr->alnReadLenUI
   ) retSC = 1; /*fist cluster is longer*/

   else if(
        firstConSTPtr->samSTPtr->alnReadLenUI
      < secConSTPtr->samSTPtr->alnReadLenUI
   ){ /*Else If: second cluster is longer*/
      swap_con_clustST(
         firstConSTPtr,
         secConSTPtr
      );

      retSC = -1;
   } /*Else If: second cluster is longer*/

   else if(
        firstConSTPtr->numReadsUL
      < secConSTPtr->numReadsUL
   ){ /*Else If: second cluster has more reads*/
      swap_con_clustST(
         firstConSTPtr,
         secConSTPtr
      );

      retSC = -1;
   } /*Else If: second cluster has more reads*/

   else
      retSC = 1; /*first cluster is equal or better*/

   firstConSTPtr->startUI =
      min_genMath(
         firstConSTPtr->startUI,
         secConSTPtr->startUI
      );

   firstConSTPtr->endUI =
      max_genMath(
         firstConSTPtr->endUI,
         secConSTPtr->endUI
      );


   /*****************************************************\
   * Fun23 Sec03 Sub03:
   *   - merge clusters and return keep consensus
   \*****************************************************/

   for(
      distSL = 0;
      distSL < (slong) indexSTPtr->lenSamUL;
      ++distSL
   ){ /*Loop: merge clusters*/
      if(indexSTPtr->clustArySI[distSL] == changeClustSI)
         indexSTPtr->clustArySI[distSL] = keepClustSI;
   } /*Loop: merge clusters*/

   return retSC;
} /*cmpCons_clustST*/

/*-------------------------------------------------------\
| Fun24: plist_con_clustST
|   - print the clusters consensuses in a con_clustST list
| Input:
|   - conSTPtr:
|     o pointer to con_clustST struct list with
|       consensuses to print
|   - headerStr:
|     o c-string with header to print (null = no header)
|   - pgHeadStr:
|     o c-string with program header to print (null = no
|       header)
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
|     o def_memErr_clustST for memeory errors
\-------------------------------------------------------*/
signed char
plist_con_clustST(
   struct con_clustST *conSTPtr, /*consensuses to print*/
   signed char *headerStr,       /*sam file header*/
   signed char *pgHeadStr,       /*program header*/
   signed char **buffStrPtr,     /*for printing*/
   unsigned long *lenBuffULPtr,  /*size of buffStrPtr*/
   void *outFILE                 /*file to print to*/
){
   schar errSC = 0;
   schar *tmpStr = 0;
   struct samEntry *samSTPtr = 0;

   if(headerStr)
      fprintf(
         outFILE,
         "%s",
         headerStr
      ); /*print headers for sam file*/

   if(pgHeadStr)
      fprintf(
         outFILE,
         "%s",
         pgHeadStr
      ); /*prints the program header for sam file*/

   while(conSTPtr)
   { /*Loop: print consensuses*/
      samSTPtr = conSTPtr->samSTPtr;

      if(! samSTPtr)
      { /*If: cluster is empty*/
         conSTPtr = conSTPtr->nextST;
         continue;
      } /*If: cluster is empty*/

      /*most of time query id not informative*/
      tmpStr = samSTPtr->qryIdStr;
      samSTPtr->lenQryIdUC = 0;

      *tmpStr++ = 'c';
      ++samSTPtr->lenQryIdUC;
      *tmpStr++ = 'l';
      ++samSTPtr->lenQryIdUC;
      *tmpStr++ = 'u';
      ++samSTPtr->lenQryIdUC;
      *tmpStr++ = 's';
      ++samSTPtr->lenQryIdUC;
      *tmpStr++ = 't';
      ++samSTPtr->lenQryIdUC;
      *tmpStr++ = ':';
      ++samSTPtr->lenQryIdUC;

      samSTPtr->lenQryIdUC +=
         numToStr(
            tmpStr,
            conSTPtr->clustSI
         );

      errSC =
         p_samEntry(
            samSTPtr,
            buffStrPtr,
            lenBuffULPtr,
            1,
            outFILE
         );

      if(errSC)
         goto memErr_fun24;

      fprintf(
         outFILE,
         "\tclust:i:%i\n",
         conSTPtr->clustSI
      );

      conSTPtr = conSTPtr->nextST;
      samSTPtr = 0;
   } /*Loop: print consensuses*/

   fflush(outFILE);
   return 0;

   memErr_fun24:;
   return def_memErr_clustST;
} /*plist_con_edCLust*/

/*-------------------------------------------------------\
| Fun25: getClust_clustST
|   - extracts reads for a signle cluster
| Input:
|   - indexSTPtr:
|     o pointer to index_clustST struct with clusters
|   - clustSI:
|     o cluster number to extract reads for
|   - headStr:
|     o program header to print
|   - samSTPtr:
|     o pointer to samEntry struct for reading sam file
|   - buffStrPtr:
|     o pointer to c-string to use reading sam file
|   - lenBuffULPtr:
|     o pointer to unsigned long with size of buffStrPtr
|   - samFILE:
|     o sam file to extract reads from
|   - outFILE:
|     o file to save reads to
| Output:
|   - Modifies:
|     o samSTPtr to be last read in file
|     o samFILE to be at start
|     o buffStrPtr to have last sam file entry
|     o lenBuffULPtr to have buffStrPtr size (if resized)
|   - Returns:
|     o 0 for no errors
|     o def_memErr for memory errors
\-------------------------------------------------------*/
signed char
getClust_clustST(
   struct index_clustST *indexSTPtr,
   signed int clustSI,
   signed char *headStr,
   struct samEntry *samSTPtr,
   signed char **buffStrPtr,
   unsigned long *lenBuffULPtr,
   void *samFILE,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun25 TOC:
   '   - finds number of reads mapping to best read
   '   o fun25 sec01:
   '     - variable declarations
   '   o fun25 sec02:
   '     - move to start of sam file
   '   o fun25 sec03:
   '     - extract reads
   '   o fun25 sec04:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun25 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   slong lineSL = 0;
   ulong offsetUL = 0;

   schar headBl = 1; /*need to print out header*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun25 Sec02:
   ^   - move to start of sam file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fseek(
      samFILE,
      0,
      SEEK_SET
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun25 Sec03:
   ^   - extract reads
   ^   o fun25 sec03 sub01:
   ^     - check if read in cluster (or is header)
   ^   o fun25 sec03 sub02:
   ^     - get read
   ^   o fun25 sec03 sub03:
   ^     - print read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun25 Sec03 Sub01:
   *   - check if read in cluster (or is header)
   \*****************************************************/

   for(
      lineSL = 0;
      lineSL < (slong) indexSTPtr->lenSamUL;
      ++lineSL
   ){ /*Loop: find reads assigned to cluster*/

      if(
             indexSTPtr->clustArySI[lineSL]
          == def_header_clustST
      ) ; /*is a header; print out*/

      else if(indexSTPtr->clustArySI[lineSL] != clustSI)
      { /*If: read not in cluster*/
         offsetUL += indexSTPtr->lenLineAryUI[lineSL];
         continue;
      } /*If: read not in cluster*/

      /**************************************************\
      * Fun25 Sec03 Sub02:
      *   - get read
      \**************************************************/

      errSC =
         getRead_clustST(
            samSTPtr,
            (ulong) lineSL,
            offsetUL,
            indexSTPtr,
            *buffStrPtr,
            samFILE
         ); /*get read from sam file*/

      offsetUL = 0;

      if(errSC)
         goto err_fun25_sec04;

      /**************************************************\
      * Fun25 Sec03 Sub03:
      *   - print read
      \**************************************************/
 
      if(samSTPtr->extraStr[0] == '@')
      { /*If: comment*/
         fprintf(
            (FILE *) outFILE,
            "%s\n",
            samSTPtr->extraStr
         );

         continue;
      } /*If: comment*/

      if(headBl)
      { /*If: need to print program header*/
         headBl = 0;

         fprintf(
            (FILE *) outFILE,
            "%s\n",
            headStr
         );
      } /*If: need to print program header*/

      if(indexSTPtr->clustArySI[lineSL] == clustSI)
      { /*Else If: printing read*/
         errSC =
            p_samEntry(
               samSTPtr,
               buffStrPtr,
               lenBuffULPtr,
               1,            /*want to print cluster #*/
               outFILE
            ); /*print out the entry*/

         if(errSC == def_memErr_samEntry)
            goto err_fun25_sec04;

         fprintf(
            (FILE *) outFILE,
            "\tclust:i:%i\n",
            clustSI
         );
      } /*Else If: printing read*/
   } /*Loop: find number of matches*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun25 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun25_sec04;


   err_fun25_sec04:;
   errSC = def_memErr_clustST;
   goto ret_fun25_sec04;


   ret_fun25_sec04:;

   fseek(
      samFILE,
      0,
      SEEK_SET
   );

   return errSC;
} /*getClust_clustST*/

/*-------------------------------------------------------\
| Fun26: pbins_clustST
|   - extracts reads for each cluster to a sam file
| Input:
|   - prefixStr:
|     o c-string with prefix for output file name
|   - clustSI:
|     o last cluster to extract reads for, should be last
|       cluster in targArySI, but does not need to be
|   - indexSTPtr:
|     o pointer to index_clustST struct with clusters
|   - pgHeadStr:
|     o program header to print (null = no print)
|   - samSTPtr:
|     o pointer to samEntry struct for reading sam file
|   - buffStrPtr:
|     o pointer to c-string to use reading sam file
|   - lenBuffULPtr:
|     o pointer to unsigned long with size of buffStrPtr
|   - samFILE:
|     o sam file to extract reads from
| Output:
|   - Modifies:
|     o samSTPtr to be last read in file
|     o samFILE to be at start
|     o buffStrPtr to have last sam file entry
|     o lenBuffULPtr to have buffStrPtr size (if resized)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_clustST for memory errors
|     o def_fileErr_clustST for file errors
\-------------------------------------------------------*/
signed char
pbins_clustST(
   signed char *prefixStr,
   signed int clustSI,
   struct index_clustST *indexSTPtr,
   signed char *pgHeadStr,
   struct samEntry *samSTPtr,
   signed char **buffStrPtr,
   unsigned long *lenBuffULPtr,
   void *samFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun26 TOC:
   '   - extracts reads for each cluster to a sam file
   '   o fun26 sec01:
   '     - variable declarations
   '   o fun26 sec02:
   '     - allocate memory for header
   '   o fun26 sec03:
   '     - get header
   '   o fun26 sec04:
   '     - extract reads
   '   o fun26 sec05:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   slong lineSL = 0;

   schar *tmpStr = 0;
   schar *headHeapStr = 0;
   uint lenHeadUI = 4096;
   uint bytesUI = 0;

   schar *clustHeapAryBl = 0;  /*tells if printed header*/
   schar outFileStr[2048]; /*holds output file name*/
   uint lenPrefixUI = endStr_ulCp(prefixStr);
   FILE *outFILE = 0;

   ulong offsetUL = 0; /*offset to next read*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec02:
   ^   - allocate memory for header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fseek(
      samFILE,
      0,
      SEEK_SET
   );


   for(
      lineSL = 0;
      lineSL < (slong) indexSTPtr->lenSamUL;
      ++lineSL
   ){ /*Loop: get header length*/
      if(
             indexSTPtr->clustArySI[lineSL]
          != def_header_clustST
      ) break;; /*on reads*/

      lenHeadUI += indexSTPtr->lenLineAryUI[lineSL];
   } /*Loop: get header length*/


   headHeapStr = malloc((lenHeadUI + 1) * sizeof(schar));

   if(! headHeapStr)
      goto memErr_fun26_sec05;

   if(clustSI > 0)
      clustHeapAryBl =
         calloc(
            clustSI + 1,
            sizeof(schar)
         );
   else
      clustHeapAryBl =
         calloc(
            2,
            sizeof(schar)
         );

   if(! clustHeapAryBl)
      goto memErr_fun26_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec03:
   ^   - get header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      lineSL = 0;
      lineSL < (slong) indexSTPtr->lenSamUL;
      ++lineSL
   ){ /*Loop: get header*/

      if(
             indexSTPtr->clustArySI[lineSL]
          != def_header_clustST
      ) break;; /*on reads*/

      offsetUL =
         fread(
            *buffStrPtr,
            sizeof(schar),
            indexSTPtr->lenLineAryUI[lineSL],
            samFILE
         );

      (*buffStrPtr)[offsetUL] = '\0';

      if(! offsetUL)
         goto fileErr_fun26_sec05;

      cpLen_ulCp(
         &headHeapStr[bytesUI],
         *buffStrPtr,
         indexSTPtr->lenLineAryUI[lineSL]
      );

      bytesUI += indexSTPtr->lenLineAryUI[lineSL];
   } /*Loop: get header*/

   offsetUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec04:
   ^   - extract reads
   ^   o fun26 sec04 sub01:
   ^     - filter out non-cluster reads + start loop
   ^   o fun26 sec04 sub02:
   ^     - get read
   ^   o fun26 sec04 sub03:
   ^     - build file name for the reads cluster
   ^   o fun26 sec04 sub04:
   ^     - open the clusters sam file (bin)
   ^   o fun26 sec04 sub05:
   ^     - print read to clusters sam file
   ^   o fun26 sec04 sub06:
   ^     - move to next read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun26 Sec04 Sub01:
   *   - filter out non-cluster reads + start loop
   \*****************************************************/

   while(lineSL < (slong) indexSTPtr->lenSamUL)
   { /*Loop: print clusters*/
      if(indexSTPtr->clustArySI[lineSL] <= 0)
      { /*If: not a cluster*/
         offsetUL += indexSTPtr->lenLineAryUI[lineSL];
         ++lineSL;
         continue;
      } /*If: not a cluster*/

      if(indexSTPtr->clustArySI[lineSL] > clustSI)
      { /*If: beyond cluster range*/
         offsetUL += indexSTPtr->lenLineAryUI[lineSL];
         ++lineSL;
         continue;
      } /*If: beyond cluster range*/

      /**************************************************\
      * Fun25 Sec03 Sub02:
      *   - get read
      \**************************************************/

      errSC =
         getRead_clustST(
            samSTPtr,
            (ulong) lineSL,
            offsetUL,
            indexSTPtr,
            *buffStrPtr,
            samFILE
         ); /*get read from sam file*/

      offsetUL = 0;

      if(errSC)
      { /*If: had error*/
         if(errSC == def_memErr_samEntry)
            goto memErr_fun26_sec05;
         if(errSC == def_fileErr_samEntry)
            goto fileErr_fun26_sec05;
      } /*If: had error*/

      /**************************************************\
      * Fun26 Sec04 Sub03:
      *   - build file name for the reads cluster
      \**************************************************/

       tmpStr = outFileStr;

       cpLen_ulCp(
          tmpStr,
          prefixStr,
          lenPrefixUI
       ); /*ad prefix to file name*/

       tmpStr += lenPrefixUI;
       *tmpStr++ = '-';

       cpLen_ulCp(
          tmpStr,
          samSTPtr->refIdStr,
          samSTPtr->lenRefIdUC
       ); /*add reference id to file name*/

       tmpStr += samSTPtr->lenRefIdUC;
       *tmpStr++ = '-';

       tmpStr +=
          numToStr(
             tmpStr,
             indexSTPtr->clustArySI[lineSL]
          ); /*add cluster number to file name*/

       *tmpStr++ = '.';
       *tmpStr++ = 's';
       *tmpStr++ = 'a';
       *tmpStr++ = 'm';
       *tmpStr = '\0';

      /**************************************************\
      * Fun26 Sec04 Sub04:
      *   - open the clusters sam file (bin)
      \**************************************************/

       if(
            clustHeapAryBl[indexSTPtr->clustArySI[lineSL]]
          == 0
       ){ /*If: first read in cluster*/

          outFILE =
             fopen(
                (char *) outFileStr,
                "w"
             );

          if(! outFILE)
             goto fileErr_fun26_sec05;

          fprintf(
             outFILE,
             "%s",
             headHeapStr /*has new line at end*/
          ); /*print header for sam file*/

          if(pgHeadStr)
             fprintf(
                outFILE,
                "%s",
                pgHeadStr
             ); /*print header for sam file*/

          clustHeapAryBl[
             indexSTPtr->clustArySI[lineSL]
          ] = 1;
       } /*If: first read in cluster*/

       else
       { /*Else: not a new file*/
          outFILE =
             fopen(
                (char *) outFileStr,
                "a"
             );

          if(! outFILE)
             goto fileErr_fun26_sec05;
       } /*Else: not a new file*/

      /**************************************************\
      * Fun26 Sec04 Sub05:
      *   - print read to clusters sam file
      \**************************************************/

       errSC =
          p_samEntry(
             samSTPtr,
             buffStrPtr,
             lenBuffULPtr,
             1,            /*want to add cluster number*/
             outFILE
          ); /*print out the entry*/

       if(errSC == def_memErr_samEntry)
          goto memErr_fun26_sec05;

       fprintf(
          outFILE,
          "\tclust:i:%i\n",
          indexSTPtr->clustArySI[lineSL]
       );

      fclose(outFILE);
      outFILE = 0;

      /**************************************************\
      * Fun26 Sec04 Sub06:
      *   - move to next read
      \**************************************************/

      ++lineSL;
   } /*Loop: print clusters*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec05:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
      goto ret_fun26_sec05;


   memErr_fun26_sec05:;
      errSC = def_memErr_clustST;
      goto ret_fun26_sec05;


   fileErr_fun26_sec05:;
      errSC = def_fileErr_clustST;
      goto ret_fun26_sec05;


   ret_fun26_sec05:;
      if(outFILE)
         fclose(outFILE);
      outFILE = 0;

      if(headHeapStr)
         free(headHeapStr);
      headHeapStr = 0;

      if(clustHeapAryBl)
         free(clustHeapAryBl);
      clustHeapAryBl = 0;

      fseek(
         samFILE,
         0,
         SEEK_SET
      );

      return errSC;
} /*pbins_clustST*/

/*-------------------------------------------------------\
| Fun27: blank_var_clustST
|   - blanks var_clustST for use with a new reference
| Input:
|   - varSTPtr:
|     o pointer to var_clustST struct to blank
| Output:
|   - Modifies:
|     o sets all values in aAryUI, cAryUI, gAryUI,
|       tAryUI, delAryUI, and insAryUI to 0
|     o sets lenRefUI to 0
\-------------------------------------------------------*/
void
blank_var_clustST(
   struct var_clustST *varSTPtr
){
   unsigned int uiBase = 0;

   varSTPtr->refNumUI = 0;
   varSTPtr->lenRefUI = 0;

   for(
      uiBase = 0;
      uiBase < varSTPtr->lenUI;
      ++uiBase
   ){ /*Loop: blank arrays*/
      if(varSTPtr->aAryUI)
         varSTPtr->aAryUI[uiBase] = 0;

      if(varSTPtr->cAryUI)
         varSTPtr->cAryUI[uiBase] = 0;

      if(varSTPtr->gAryUI)
         varSTPtr->gAryUI[uiBase] = 0;

      if(varSTPtr->tAryUI)
         varSTPtr->tAryUI[uiBase] = 0;

      if(varSTPtr->delAryUI)
         varSTPtr->delAryUI[uiBase] = 0;

      if(varSTPtr->insAryUI)
         varSTPtr->insAryUI[uiBase] = 0;
   } /*Loop: blank arrays*/
} /*blank_var_clustST*/

/*-------------------------------------------------------\
| Fun28: init_var_clustST
|   - initializes a var_clustST array
| Input:
|   - varSTPtr:
|     o pointer to var_clustST struct to initialize
| Output:
|   - Modifies:
|     o sets all values to 0
\-------------------------------------------------------*/
void
init_var_clustST(
   struct var_clustST *varSTPtr
){
   varSTPtr->aAryUI = 0;
   varSTPtr->cAryUI = 0;
   varSTPtr->gAryUI = 0;
   varSTPtr->tAryUI = 0;
   varSTPtr->delAryUI = 0;
   varSTPtr->insAryUI = 0;

   varSTPtr->lenUI = 0;
   blank_var_clustST(varSTPtr);
} /*init_var_clustST*/

/*-------------------------------------------------------\
| Fun29: freeStack_var_clustST
|   - frees arrays in a var_clustST struct
| Input:
|   - varSTPtr:
|     o pointer to var_clustST with arrays to free
| Output:
|   - Frees:
|     o all arrays (and initializes)
\-------------------------------------------------------*/
void
freeStack_var_clustST(
   struct var_clustST *varSTPtr
){
   if(! varSTPtr)
      return; /*no struct input*/

   if(varSTPtr->aAryUI)
      free(varSTPtr->aAryUI);

   if(varSTPtr->cAryUI)
      free(varSTPtr->cAryUI);

   if(varSTPtr->gAryUI)
      free(varSTPtr->gAryUI);

   if(varSTPtr->tAryUI)
      free(varSTPtr->tAryUI);

   if(varSTPtr->delAryUI)
      free(varSTPtr->delAryUI);

   if(varSTPtr->insAryUI)
      free(varSTPtr->insAryUI);

   init_var_clustST(varSTPtr);
} /*freeStack_var_clustST*/

/*-------------------------------------------------------\
| Fun30: freeHeap_var_clustST
|   - frees a heap allocated var_clustST struct
| Input:
|   - varSTPtr:
|     o pointer to var_clustST to free
| Output:
|   - Frees:
|     o varSTPtr (you must set to 0)
\-------------------------------------------------------*/
void
freeHeap_var_clustST(
   struct var_clustST *varSTPtr
){
   if(varSTPtr)
   { /*If: have struct to free*/
      freeStack_var_clustST(varSTPtr);
      free(varSTPtr);
   } /*If: have struct to free*/
} /*freeHeap_var_clustST*/

/*-------------------------------------------------------\
| Fun31: freeHeapAry_var_clustST
|   - frees an array of heap allocated var_clustST
| Input:
|   - varSTPtr:
|     o pointer to var_clustST array to free
|   - lenUI:
|     o length of var_clustST array (index 1)
| Output:
|   - Frees:
|     o all var_clustST structs in varSTPtr
|      (you must set varSTPtr to 0)
\-------------------------------------------------------*/
void
freeHeapAry_var_clustST(
   struct var_clustST *varSTPtr,
   signed int lenSI
){
   if(! varSTPtr)
      return; /*not structs input*/

   --lenSI;

   while(lenSI >= 0)
      freeHeap_var_clustST(&varSTPtr[lenSI]);

   free(varSTPtr);
} /*freeHeapAry_var_clustST*/

/*-------------------------------------------------------\
| Fun32: setupRef_var_clustST
|   - setsup a var_clustST for a reference
| Input:
|   - varSTPtr:
|     o pointer to var_clustST to set up
|   - lenRefUI:
|     o length of reference to add (index 1)
| Output:
|   - Modifies:
|     o arrays in varSTPtr to have allocated memory
|     o lenRefUI in varSTPtr to be lenRefUI
|     o lenUI in varSTPtr to be lenRefUI (if shorter)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_clustST for memory errors
\-------------------------------------------------------*/
signed char
setupRef_var_clustST(
   struct var_clustST *varSTPtr,
   unsigned int lenRefUI
){
   varSTPtr->lenRefUI = lenRefUI;

   if(varSTPtr->lenUI < lenRefUI)
   { /*If: need to resiaze arrays*/
      freeStack_var_clustST(varSTPtr);
      varSTPtr->lenUI = lenRefUI;
   } /*If: need to resiaze arrays*/

   if(! varSTPtr->aAryUI)
   { /*If: a array needs memory*/
      varSTPtr->aAryUI = malloc(lenRefUI * sizeof(uint));

      if(! varSTPtr->aAryUI)
         goto memErr_fun06;
   } /*If: a array needs memory*/

   if(! varSTPtr->cAryUI)
   { /*If: c array needs memory*/
      varSTPtr->cAryUI = malloc(lenRefUI * sizeof(uint));

      if(! varSTPtr->cAryUI)
         goto memErr_fun06;
   } /*If: c array needs memory*/

   if(! varSTPtr->gAryUI)
   { /*If: g array needs memory*/
      varSTPtr->gAryUI = malloc(lenRefUI * sizeof(uint));

      if(! varSTPtr->gAryUI)
         goto memErr_fun06;
   } /*If: g array needs memory*/

   if(! varSTPtr->tAryUI)
   { /*If: t array needs memory*/
      varSTPtr->tAryUI = malloc(lenRefUI * sizeof(uint));

      if(! varSTPtr->tAryUI)
         goto memErr_fun06;
   } /*If: t array needs memory*/

   if(! varSTPtr->delAryUI)
   { /*If: del array needs memory*/
      varSTPtr->delAryUI =malloc(lenRefUI * sizeof(uint));

      if(! varSTPtr->delAryUI)
         goto memErr_fun06;
   } /*If: del array needs memory*/

   if(! varSTPtr->insAryUI)
   { /*If: ins array needs memory*/
      varSTPtr->insAryUI =malloc(lenRefUI * sizeof(uint));

      if(! varSTPtr->insAryUI)
         goto memErr_fun06;
   } /*If: ins array needs memory*/

   blank_var_clustST(varSTPtr);

   return 0;

   memErr_fun06:;
   return def_memErr_clustST;
} /*setupRef_var_clustST*/

/*-------------------------------------------------------\
| Fun33: rmLowDepth_var_clustST
|   - remove low read depth variants in a var_clustST
| Input:
|   - varSTPtr:
|     o pointer to var_clustST to set up
|   - minDepthUI:
|     o min depth to keep a variant
| Output:
|   - Modifies:
|     o arrays in varSTPtr to have low read depth
|       positions set to 0
\-------------------------------------------------------*/
void
rmLowDepth_var_clustST(
   struct var_clustST *varSTPtr,
   unsigned int minDepthUI
){
   unsigned int uiBase = 0;

   for(
      uiBase = 0;
      uiBase < varSTPtr->lenRefUI;
      ++uiBase
   ){ /*Loop: remove low read depth variants*/
      if(varSTPtr->aAryUI[uiBase] < minDepthUI)
         varSTPtr->aAryUI[uiBase] = 0;

      if(varSTPtr->cAryUI[uiBase] < minDepthUI)
         varSTPtr->cAryUI[uiBase] = 0;

      if(varSTPtr->gAryUI[uiBase] < minDepthUI)
         varSTPtr->gAryUI[uiBase] = 0;

      if(varSTPtr->tAryUI[uiBase] < minDepthUI)
         varSTPtr->tAryUI[uiBase] = 0;

      if(varSTPtr->delAryUI[uiBase] < minDepthUI)
         varSTPtr->delAryUI[uiBase] = 0;

      if(varSTPtr->insAryUI[uiBase] < minDepthUI)
         varSTPtr->insAryUI[uiBase] = 0;
   } /*Loop: remove low read depth variants*/
} /*rmLowDepth_var_clustST*/

/*-------------------------------------------------------\
| Fun34: blank_hist_clustST
|   - blanks hist_clustST for use
| Input:
|   - histSTPtr:
|     o pointer to hist_clustST struct to blank
| Output:
|   - Modifies:
|     o sets all values (except lenUI) to 0
|     o sets all arrays to be filled with 0's
\-------------------------------------------------------*/
void
blank_hist_clustST(
   struct hist_clustST *histSTPtr
){
   uint uiHist = 0;

   histSTPtr->startUI = 0;
   histSTPtr->endUI = 0;
   histSTPtr->numVarUI = 0;
   histSTPtr->depthUL = 0;

   for(
      uiHist = 0;
      uiHist < histSTPtr->lenUI;
      ++uiHist
   ){ /*Loop: blank arrays*/
      if(histSTPtr->indexAryUI)
         histSTPtr->indexAryUI[uiHist] = 0;

      if(histSTPtr->typeAryUC)
         histSTPtr->typeAryUC[uiHist] = 0;

      if(histSTPtr->profAryUI)
         histSTPtr->profAryUI[uiHist] = 0;

      if(histSTPtr->histAryUI)
         histSTPtr->histAryUI[uiHist] = 0;

      if(histSTPtr->depthAryUL)
         histSTPtr->depthAryUL[uiHist] = 0;
   } /*Loop: blank arrays*/
} /*blank_hist_clustST*/

/*-------------------------------------------------------\
| Fun35: init_hist_clustST
|   - initializes hist_clustST for use
| Input:
|   - histSTPtr:
|     o pointer to hist_clustST struct to initialize
| Output:
|   - Modifies:
|     o sets all values to 0
\-------------------------------------------------------*/
void
init_hist_clustST(
   struct hist_clustST *histSTPtr
){
   histSTPtr->lenUI = 0;

   histSTPtr->indexAryUI = 0;
   histSTPtr->typeAryUC = 0;
   histSTPtr->profAryUI = 0;
   histSTPtr->histAryUI = 0;
   histSTPtr->depthAryUL = 0;

   blank_hist_clustST(histSTPtr);
} /*init_hist_clustST*/

/*-------------------------------------------------------\
| Fun36: freeStack_hist_clustST
|   - frees variables inside a hist_clustST
| Input:
|   - histSTPtr:
|     o pointer to hist_clustST with variables to free
| Output:
|   - Modifies:
|     o frees all arrays, then initializes
\-------------------------------------------------------*/
void
freeStack_hist_clustST(
   struct hist_clustST *histSTPtr
){
   if(! histSTPtr)
      return; /*nothing to free*/

   if(histSTPtr->indexAryUI)
      free(histSTPtr->indexAryUI);

   if(histSTPtr->typeAryUC)
      free(histSTPtr->typeAryUC);

   if(histSTPtr->profAryUI)
      free(histSTPtr->profAryUI);

   if(histSTPtr->histAryUI)
      free(histSTPtr->histAryUI);

   if(histSTPtr->depthAryUL)
      free(histSTPtr->depthAryUL);

   init_hist_clustST(histSTPtr);
} /*freeStack_hist_clustST*/

/*-------------------------------------------------------\
| Fun37: freeHeap_hist_clustST
|   - frees a hist_clustST struct
| Input:
|   - histSTPtr:
|     o pointer to hist_clustST struct to free
| Output:
|   - Frees:
|     o histSTPtr (you must set to 0)
\-------------------------------------------------------*/
void
freeHeap_hist_clustST(
   struct hist_clustST *histSTPtr
){
   if(histSTPtr)
   { /*If: have struct to free*/
      freeStack_hist_clustST(histSTPtr);
      free(histSTPtr);
   } /*If: have struct to free*/
} /*freeHeap_hist_clustST*/

/*-------------------------------------------------------\
| Fun38: prob_clustST
|   - find probability that a profile is supported
|   - probablity calcualtion is from hairSplitter
|     Roland Faure, Dominique Lavenier,
|     and Jean-Francios Flot; HairSplittler: haplotype
|     assembly from long, nosiy reads; preprint; 2024
|     https://doi.org/10.1101/2024.02/13.580067
| Input:
|   - numVarUI:
|     o number of variants in profile
|   - lenTargUI:
|     o number of bases in target region
|   - depthUL:
|     o depth of profile
|   - maxDepthUL:
|     o maximum number of reads
|   - errF:
|     o expected error rate (as percentage)
| Output:
|   - Returns:
|     o probability of profile
\-------------------------------------------------------*/
double
prob_clustST(
   unsigned int numVarUI,
   unsigned int lenTargUI,
   unsigned long depthUL,
   unsigned long maxDepthUL,
   float errF
){
   uint uiProb = 0;
   double probDbl = 0;

   /*b choose m (variants)*/
   for(
      uiProb = 0;
      uiProb < numVarUI;
      ++uiProb
   ){ /*Loop: find number of variations*/
      probDbl *=
        (
             (double) (lenTargUI - uiProb)
           / (double) (uiProb + 1)
        ) * errF;
   } /*Loop: find number of variations*/

   /*(b choose m)*/
   for(
      uiProb = 0;
      uiProb < depthUL;
      ++uiProb
   ){ /*Loop: find combinations for read depth*/
      probDbl *=
        (
             (double) (maxDepthUL - uiProb)
           / (double) (uiProb + 1)
        ) * errF;
   } /*Loop: find combinations for read depth*/

   /*finish off e^[a * b]*/
   for(
      uiProb = numVarUI + depthUL;
      uiProb < numVarUI * depthUL;
      ++uiProb
   ) probDbl *= errF;

   return probDbl;
} /*prob_clustST*/

/*-------------------------------------------------------\
| Fun39: roughProb_clustST
|   - finds expected power of 10 of profile being unique
|   - always use prob_clustST to verify
| Input:
|   - numVarUI:
|     o number of variants in profile
|   - lenTargUI:
|     o number of bases in target region
|   - depthUL:
|     o depth of profile
|   - maxDepthUL:
|     o maximum number of reads
|   - errF:
|     o expected error rate (as percentage)
| Output:
|   - Returns:
|     o a rough estimate of the power 10 of probability
\-------------------------------------------------------*/
signed long
roughProb_clustST(
   unsigned int numVarUI,
   unsigned int lenTargUI,
   unsigned long depthUL,
   unsigned long maxDepthUL,
   float errF
){
   uint errDigUI = 0;
   uint lenDigUI = 0;
   uint depthDigUI = 0;
   uint tmpSI = 0;

   do{ /*Loop: find power of error*/
      tmpSI = errF * 10;
      ++errDigUI;
      errF *= 10;
   } while(tmpSI == 0);
     /*Loop: find power of error*/

   do{ /*Loop: find power of length*/
      tmpSI = lenTargUI / 10;
      ++lenDigUI;
      lenTargUI /= 10;
   } while(tmpSI != 0);
     /*Loop: find power of length*/

   do{ /*Loop: find power of depth*/
      tmpSI = maxDepthUL / 10;
      ++depthDigUI;
      maxDepthUL /= 10;
   } while(tmpSI != 0);
     /*Loop: find power of depth*/

   return
      (slong)
        ( (depthDigUI * depthUL) + (lenDigUI * numVarUI) )
      - (errDigUI * numVarUI * depthUL);
      /*this is very rough, but it should give me an idea
      `   if two profiles are in the same power of 10
      */
} /*roughProb_clustST*/

/*-------------------------------------------------------\
| Fun40: findBestProf_hist_clustST
|   - finds best profile in a hist_clustST history
| Input:
|   - histSTPtr:
|     o pointer to hist_clustST with history to search
|   - probDblPtr:
|     o pointer to double to hold probability
|   - errRateF:
|     o expected error rate (%)
| Output:
|   - Modifies:
|     o probDblPtr to have probability of best history
|   - Returns:
|     o index of history with best probability
\-------------------------------------------------------*/
unsigned int
findBestProf_hist_clustST(
   struct hist_clustST *histSTPtr, /*has history*/
   double *probDblPtr,             /*will have odds*/
   float errRateF                  /*expected error*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun40 TOC:
   '   - finds best profile in a hist_clustST history
   '   o fun40 sec01:
   '     - variable declarations
   '   o fun40 sec02:
   '     - find best profile in history
   '   o fun40 sec03:
   '     - find probability of best profile in history
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun40 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint uiHist = 0;
   uint lenUI = histSTPtr->endUI - histSTPtr->startUI + 1;

   signed long probSL = 0;
   double probDbl = 0;

   schar haveBestProbBl = 0;

   unsigned int bestIndexUI = 0;
   signed long bestProbSL = 0;
   double bestProbDbl = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun40 Sec02:
   ^   - find best profile in history
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      uiHist = 0;
      uiHist < histSTPtr->lenUI;
      ++uiHist
   ){ /*Loop: find probabilities*/
      probSL =
        roughProb_clustST(
           histSTPtr->histAryUI[uiHist], /*# variants*/
           lenUI,                        /*overlap len*/
           histSTPtr->depthAryUL[uiHist],/*depth of prof*/
           histSTPtr->depthUL,           /*read depth*/
           errRateF                      /*error rate*/
        );

      if(bestProbSL > probSL)
      { /*If: new best*/
         bestProbSL = probSL;
         bestIndexUI = uiHist;
         *probDblPtr = 0;
         haveBestProbBl = 0;
         continue;
      } /*If: new best*/

      else if(bestProbSL < probSL)
          continue; /*old probability is better*/

      /*this is the acutal probability. I use rough
      `  becuase it is quicker
      */
      if(! haveBestProbBl)
      { /*If: have not found probability yet*/
         *probDblPtr =
            prob_clustST(
              histSTPtr->histAryUI[bestIndexUI],
              lenUI,                       /*overlap len*/
              histSTPtr->depthAryUL[bestIndexUI],
              histSTPtr->depthUL,           /*read depth*/
              errRateF                      /*error rate*/
            );

         haveBestProbBl = 1;
      } /*If: have not found probability yet*/
      
      probDbl =
        prob_clustST(
           histSTPtr->histAryUI[uiHist], /*# variants*/
           lenUI,                        /*overlap len*/
           histSTPtr->depthAryUL[uiHist],/*depth of prof*/
           histSTPtr->depthUL,           /*read depth*/
           errRateF                      /*error rate*/
        );

      if(bestProbDbl > probDbl)
      { /*If: new best probability*/
         bestProbSL = probSL;
         bestIndexUI = uiHist;
         *probDblPtr = probDbl;
      } /*If: new best probability*/
   } /*Loop: find probabilities*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun40 Sec03:
   ^   - find probability of best profile in history
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! haveBestProbBl)
   { /*If: have not found probability yet*/
      *probDblPtr =
         prob_clustST(
           histSTPtr->histAryUI[bestIndexUI],
           lenUI,                       /*overlap len*/
           histSTPtr->depthAryUL[bestIndexUI],
           histSTPtr->depthUL,           /*read depth*/
           errRateF                      /*error rate*/
         );
   } /*If: have not found probability yet*/

   return bestIndexUI;
} /*findBestProf_hist_clustST*/
