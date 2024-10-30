/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' illNano SOF: Start Of File
'   - extracs ONT reads with Illumina variants
'   o header:
'     - included libraries
'   o .h st01: prof_illNano
'     - holds a unique profile from illNano
'   o .h st02: profList_illNano
'     - head of a prof_illNano struct list
'   o fun01: blank_prof_illNano
'     - blanks a prof_illNano struct
'   o fun02: init_prof_illNano
'     - initialize a prof_illNano struct
'   o fun03: freeStack_prof_illNano
'     - frees arrays in prof_illNano (heap allocated)
'   o fun04: freeHeap_prof_illNano
'     - a heap allocated prof_illNano struct
'   o fun05: freeHeapList_prof_illNano
'     - a list of heap allocated prof_illNano structs
'   o fun06: mk_prof_illNano
'     - makes a prof_illNano struct
'   o fun07: cmp_prof_illNano
'     - compares two prof_illNano structs
'   o fun08: swap_prof_illNano
'     - swaps variables in two prof_illNano structs
'   o fun09: merge_prof_illNano
'     - merges two prof_illNano structs into one profile
'   o fun10: blank_profList_illNano
'     - blanks a profList_illNano struct (does nothing)
'   o fun11: init_profList_illNano
'     - initializes a profList_illNano struct
'   o fun12: freeStack_profList_illNano
'     - frees variables (list) inside a profList_illNano
'   o fun13: freeHeap_profList_illNano
'     - frees a profList_illNano struct
'   o fun14: getIllSnp_illNano
'     - gets snps from tbCon tsv file for Illumina data
'   o fun15: getVarNano_illNano
'     - identifies variants in a nanopore read using an
'       Illumina profile
'   o fun16: getNanoReads_illNano
'     - gets nanopore reads using Illumina profile
'   o fun17: phead_profList_illNano
'     - prints header for unique profile tsv
'   o fun18: p_profList_illNano
'     - prints unique profiles in a profList struct
'   o fun19: run_illNano
'     - runs the illNano algorithim
'   o license:
'     - licensing for this code (public domain / mit)
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

#include "illNano.h"

#include <stdio.h>

#include "../genLib/base10str.h"
#include "../genLib/numToStr.h"
#include "../genLib/ulCp.h"
#include "../genBio/samEntry.h"

/*.h files only*/
#include "../genLib/dataTypeShortHand.h"
#include "../genLib/genMath.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries
!   o .c   #include "../genLib/numToStr.h"
!   o .c   #include "../genLib/strAry.h"
!   o .h   #include "ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: blank_prof_illNano
|   - blanks a prof_illNano struct
| Input:
|   - profSTPtr:
|     o pointer to prof_illNano struct to blank
| Output:
|   - Modifies:
|     o varInProfUI, maxDiffUI, minDiffUI, and avgDiffF
|       to be 0
|     o all arrays to be filled with 0 (if have)
\-------------------------------------------------------*/
void
blank_prof_illNano(
   struct prof_illNano *profSTPtr
){
   uint uiPos = 0;

   profSTPtr->varInProfUI = 0;

   profSTPtr->maxDiffUI = 0;
   profSTPtr->minDiffUI = 0;

   profSTPtr->sumDiffUL = 0;
   profSTPtr->overlapUI = 0;
   profSTPtr->avgDiffF = 0;

   if(
         profSTPtr->posAryUI
      || profSTPtr->ntArySC
      || profSTPtr->depthAryUI
      || profSTPtr->xDepthAryUI
   ){ /*If: at least on array to blank*/
      for(
         uiPos = 0;
         uiPos < profSTPtr->sizeProfUI;
         ++uiPos
      ){ /*Loop: blank arrays*/
         if(profSTPtr->posAryUI)
            profSTPtr->posAryUI = 0;

         if(profSTPtr->ntArySC)
            profSTPtr->ntArySC = (schar) '\0';

         if(profSTPtr->depthAryUI)
            profSTPtr->depthAryUI = 0;

         if(profSTPtr->xDepthAryUI)
            profSTPtr->xDepthAryUI = 0;
      } /*Loop: blank arrays*/
   } /*If: at least on array to blank*/
} /*blank_prof_illNano*/

/*-------------------------------------------------------\
| Fun02: init_prof_illNano
|   - initialize a prof_illNano struct
| Input:
|   - profSTPtr:
|     o pointer to prof_illNano struct to initialize
| Output:
|   - Modifies:
|     o all variants in profSTPtr to be 0
\-------------------------------------------------------*/
void
init_prof_illNano(
   struct prof_illNano *profSTPtr
){
   profSTPtr->posAryUI = 0;
   profSTPtr->ntArySC = 0;
   profSTPtr->depthAryUI = 0;
   profSTPtr->xDepthAryUI = 0;

   profSTPtr->sizeProfUI = 0;

   profSTPtr->nextST = 0;
   profSTPtr->lastST = 0;

   blank_prof_illNano(profSTPtr);
} /*init_prof_illNano*/

/*-------------------------------------------------------\
| Fun03: freeStack_prof_illNano
|   - frees arrays in prof_illNano (assumes heap allocate)
| Input:
|   - profSTPtr:
|     o pointer to prof_illNano struct with arrays to free
| Output:
|   - Frees:
|     o posAryUI, ntArySC, depthAryUI, and xDepthAryUI in
|       profSTPtr
\-------------------------------------------------------*/
void
freeStack_prof_illNano(
   struct prof_illNano *profSTPtr
){
   if(! profSTPtr)
      return;

   if(profSTPtr->ntArySC)
      free(profSTPtr->ntArySC);

   if(profSTPtr->posAryUI)
      free(profSTPtr->posAryUI);

   if(profSTPtr->depthAryUI)
      free(profSTPtr->depthAryUI);

   if(profSTPtr->xDepthAryUI)
      free(profSTPtr->xDepthAryUI);

   init_prof_illNano(profSTPtr);
} /*freeStack_prof_illNano*/

/*-------------------------------------------------------\
| Fun04: freeHeap_prof_illNano
|   - a heap allocated prof_illNano struct
| Input:
|   - profSTPtr:
|     o pointer to prof_illNano struct to free
| Output:
|   - Frees:
|     o profSTPtr (you must set to 0)
|   - Returns:
|     o next prof_illNano struct in list
\-------------------------------------------------------*/
struct prof_illNano *
freeHeap_prof_illNano(
   struct prof_illNano *profSTPtr
){
   struct prof_illNano *nextSTPtr = 0;

   if(! profSTPtr)
      return 0;

   nextSTPtr = profSTPtr->nextST;
   freeStack_prof_illNano(profSTPtr);
   free(profSTPtr);
   return nextSTPtr;
} /*freeHeap_prof_illNano*/

/*-------------------------------------------------------\
| Fun05: freeHeapList_prof_illNano
|   - a list of heap allocated prof_illNano structs
| Input:
|   - profSTPtr:
|     o pointer to list of prof_illNano struct to free
| Output:
|   - Frees:
|     o every prof_illNano struct in profSTPtr
|       (you must set to profSTPtr to 0)
\-------------------------------------------------------*/
void
freeHeapList_prof_illNano(
   struct prof_illNano *profSTPtr
){
   while(profSTPtr)
      profSTPtr = freeHeap_prof_illNano(profSTPtr);
} /*freeHeapList_prof_illNano*/

/*-------------------------------------------------------\
| Fun06: mk_prof_illNano
|   - makes a prof_illNano struct
| Input:
|   - posAryUI:
|     o unsigned int array with positions to copy
|     o use 0 (null) to not use
|   - ntArySC:
|     o signed char array with nucleotides to copy
|     o use 0 (null) to not use
|   - numVarUI:
|     o number of variants in ntArySC and posAryUI
| Output:
|   - Returns:
|     o heap allocated prof_illNano struct with copied
|       arrays
|     o 0 for errors
\-------------------------------------------------------*/
struct prof_illNano *
mk_prof_illNano(
   unsigned int *posAryUI, /*has positions to copy*/
   signed char *ntArySC,   /*has nucleotides to copy*/
   unsigned int numVarUI   /*number variants in arrays*/
){
   uint uiVar = 0;
   struct prof_illNano *retHeapST = 0;


   /*create struct*/
   retHeapST = malloc(sizeof(struct prof_illNano));

   if(! retHeapST)
      goto memErr_fun06;

   init_prof_illNano(retHeapST);


   /*allocate array memory*/

   retHeapST->posAryUI =
      malloc((numVarUI + 9) * sizeof(uint));

   if(! retHeapST->posAryUI)
      goto memErr_fun06;


   retHeapST->ntArySC =
      malloc((numVarUI + 9) * sizeof(schar));

   if(! retHeapST->ntArySC)
      goto memErr_fun06;


   retHeapST->depthAryUI =
      malloc((numVarUI + 9) * sizeof(uint));

   if(! retHeapST->depthAryUI)
      goto memErr_fun06;


   retHeapST->xDepthAryUI =
      malloc((numVarUI + 9) * sizeof(uint));

   if(! retHeapST->xDepthAryUI)
      goto memErr_fun06;

   for(
      uiVar = 0;
      uiVar < numVarUI;
      ++uiVar
   ){ /*Loop: initialize depth to 1 read*/
      if(ntArySC[uiVar] == def_unkown_illNano)
      { /*If: unknown base*/
         retHeapST->depthAryUI[uiVar] = 0;
         retHeapST->xDepthAryUI[uiVar] = 1;
      } /*If: unknown base*/

      else
      { /*Else: actual variant base*/
         retHeapST->depthAryUI[uiVar] = 1;
         retHeapST->xDepthAryUI[uiVar] = 0;
      } /*Else: actual variant base*/
   } /*Loop: initialize depth to 1 read*/


   /*copy array values*/
   retHeapST->varInProfUI = numVarUI;
   retHeapST->sizeProfUI = numVarUI;

   cpLen_ulCp(
      (schar *) retHeapST->posAryUI,
      (schar *) posAryUI,
      numVarUI * sizeof(uint)
   );

   cpLen_ulCp(
      retHeapST->ntArySC,
      ntArySC,
      numVarUI
   );

   /*return result/clean up if needed*/

   goto ret_fun06;

   memErr_fun06:;
      if(retHeapST)
         freeHeap_prof_illNano(retHeapST);
      retHeapST = 0;
      goto ret_fun06;

   ret_fun06:;
      return retHeapST;
} /*mk_prof_illNano*/

/*-------------------------------------------------------\
| Fun07: cmp_prof_illNano
|   - compares two prof_illNano structs
| Input:
|   - firstProfSTPtr;
|     o pointer to first prof_illNano struct to compare
|   - secProfSTPtr;
|     o pointer to second prof_illNano struct to compare
|   - cmpSmallLapBl:
|     o 1: compare event when only a few variants overlap
|     o 0: require at least one profile to have at least
|          the variants present in the other profile
|          (can have more)
| Output:
|   - Returns:
|     o > -1, number of differences between profiles
|       - def_unkown_illNano bases are treated as wild
|         cards
|     o 0 no differences
|     o < 0 no overlap
\-------------------------------------------------------*/
signed int
cmp_prof_illNano(
   struct prof_illNano *firstProfSTPtr, /*first struct*/
   struct prof_illNano *secProfSTPtr,   /*second struct*/
   signed char cmpSmallLapBl/*1 do not worry about exact*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC: cmp_prof_illNano
   '   - compares two prof_illNano structs
   '   o fun07 sec01:
   '     - variable declarations
   '   o fun07 sec02:
   '     - check if have variants
   '   o fun07 sec03:
   '     - move to first shared base in profile
   '   o fun07 sec04:
   '     - count number of differences
   '   o fun07 sec05:
   '     - return the number of differences
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint uiFirst = 0;
   uint uiSec = 0;
   sint diffSI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^   - check if have variants
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! firstProfSTPtr->varInProfUI)
      goto ret_fun07_sec05;
      /*no variants in first profile*/

   if(! secProfSTPtr->varInProfUI)
      goto ret_fun07_sec05;
      /*no variants in second profile*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^   - move to first shared base in profile
   ^   o fun07 sec03 sub01:
   ^     - move till 1st profile position is not > 2nd
   ^   o fun07 sec03 sub02:
   ^     - move till 2nd profile position is not > first
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun07 Sec03 Sub01:
   *   - move till 1st profile position is not > 2nd
   \*****************************************************/

   while(
        firstProfSTPtr->posAryUI[uiFirst]
      > secProfSTPtr->posAryUI[uiSec]
   ){ /*Loop: find first matching position*/
      ++uiSec;

      if(uiFirst >= firstProfSTPtr->varInProfUI)
      { /*If: no ovlerapping variants*/
         diffSI = -1 * firstProfSTPtr->varInProfUI;
         goto ret_fun07_sec05;
      } /*If: no ovlerapping variants*/

      if(uiSec >= secProfSTPtr->varInProfUI)
      { /*If: no ovlerapping variants*/
         diffSI = -1 * secProfSTPtr->varInProfUI;
         goto ret_fun07_sec05;
      } /*If: no ovlerapping variants*/
   } /*Loop: find first matching position*/

   /*****************************************************\
   * Fun07 Sec03 Sub02:
   *   - move till 2nd profile position is not > first
   \*****************************************************/

   while(
        firstProfSTPtr->posAryUI[uiFirst]
      < secProfSTPtr->posAryUI[uiSec]
   ){ /*Loop: find first matching position*/
      ++uiFirst;

      if(uiFirst >= firstProfSTPtr->varInProfUI)
      { /*If: no ovlerapping variants*/
         diffSI = -1 * firstProfSTPtr->varInProfUI;
         goto ret_fun07_sec05;
      } /*If: no ovlerapping variants*/

      if(uiSec >= secProfSTPtr->varInProfUI)
      { /*If: no ovlerapping variants*/
         diffSI = -1 * secProfSTPtr->varInProfUI;
         goto ret_fun07_sec05;
      } /*If: no ovlerapping variants*/
   } /*Loop: find first matching position*/

   if(
         uiFirst
      && uiSec
      && ! cmpSmallLapBl
   ){ /*If: only comparing when have full coverage*/
      diffSI = -1 * firstProfSTPtr->varInProfUI;
      goto ret_fun07_sec05;
   } /*If: only comparing when have full coverage*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec04:
   ^   - count number of differences
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*profiles are order by position (lowest first) and
   `  either have variant base or consensus base, so
   `  at this point profiles will be lined up and will
   `  have 0 differences (though one might be missing a
   `  few positions)
   */

   while(uiSec < secProfSTPtr->varInProfUI)
   { /*Loop: find number of differences*/
      if(uiFirst >= firstProfSTPtr->varInProfUI)
      { /*If: had end of second profile*/
         if(! cmpSmallLapBl)
            diffSI *= -1;
            /*second reference has extra bases*/

         goto ret_fun07_sec05;
         /*at end or profile only covers part of first*/
      } /*If: had end of second profile*/

      if(
           firstProfSTPtr->ntArySC[uiFirst]
        == def_unkown_illNano
      ) ;

      else if(
            secProfSTPtr->ntArySC[uiSec]
         == def_unkown_illNano
      ) ;

      else 
      { /*Else: difference between profiles*/
         diffSI +=
            (     firstProfSTPtr->ntArySC[uiFirst]
               != secProfSTPtr->ntArySC[uiSec]
            ); /*check if difference*/
      } /*Else: difference between profiles*/

      ++uiFirst;
      ++uiSec;
   } /*Loop: find number of differences*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec05:
   ^   - return the number of differences
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ret_fun07_sec05:;
      return diffSI;
} /*cmp_prof_illNano*/

/*-------------------------------------------------------\
| Fun08: swap_prof_illNano
|   - swaps variables in two prof_illNano structs
| Input:
|   - firstProfSTPtr;
|     o pointer to first prof_illNano struct to swap
|   - secProfSTPtr;
|     o pointer to second prof_illNano struct to swap
| Output:
|   - Modifies:
|     o firstProfSTPtr to have values of secProfSTPtr
|     o secProfSTPtr to have values of firstProfSTPtr
| WARNING:
|   - does not swap nextST
\-------------------------------------------------------*/
void
swap_prof_illNano(
   struct prof_illNano *firstProfSTPtr, /*first struct*/
   struct prof_illNano *secProfSTPtr    /*second struct*/
){
   schar *swapSCPtr = 0;
   float swapF = 0;

   /*swap the arrays*/
   swapSCPtr = (schar *) firstProfSTPtr->posAryUI;
   firstProfSTPtr->posAryUI = secProfSTPtr->posAryUI;
   secProfSTPtr->posAryUI = (uint *) swapSCPtr;

   swapSCPtr = firstProfSTPtr->ntArySC;
   firstProfSTPtr->ntArySC = secProfSTPtr->ntArySC;
   secProfSTPtr->ntArySC = swapSCPtr;

   swapSCPtr = (schar *) firstProfSTPtr->depthAryUI;
   firstProfSTPtr->depthAryUI = secProfSTPtr->depthAryUI;
   secProfSTPtr->depthAryUI = (uint *) swapSCPtr;

   swapSCPtr = (schar *) firstProfSTPtr->xDepthAryUI;
   firstProfSTPtr->xDepthAryUI=secProfSTPtr->xDepthAryUI;
   secProfSTPtr->xDepthAryUI = (uint *) swapSCPtr;

   /*swap the array sizes, max, min, and mean*/
   firstProfSTPtr->sizeProfUI ^= secProfSTPtr->sizeProfUI;
   secProfSTPtr->sizeProfUI ^= firstProfSTPtr->sizeProfUI;
   firstProfSTPtr->sizeProfUI ^= secProfSTPtr->sizeProfUI;

   firstProfSTPtr->varInProfUI ^=
      secProfSTPtr->varInProfUI;
   secProfSTPtr->varInProfUI ^=
      firstProfSTPtr->varInProfUI;
   firstProfSTPtr->varInProfUI ^=
      secProfSTPtr->varInProfUI;

   firstProfSTPtr->maxDiffUI ^= secProfSTPtr->maxDiffUI;
   secProfSTPtr->maxDiffUI ^= firstProfSTPtr->maxDiffUI;
   firstProfSTPtr->maxDiffUI ^= secProfSTPtr->maxDiffUI;

   firstProfSTPtr->minDiffUI ^= secProfSTPtr->minDiffUI;
   secProfSTPtr->minDiffUI ^= firstProfSTPtr->minDiffUI;
   firstProfSTPtr->minDiffUI ^= secProfSTPtr->minDiffUI;

   firstProfSTPtr->sumDiffUL ^= secProfSTPtr->sumDiffUL;
   secProfSTPtr->sumDiffUL ^= firstProfSTPtr->sumDiffUL;
   firstProfSTPtr->sumDiffUL ^= secProfSTPtr->sumDiffUL;


   swapF = firstProfSTPtr->avgDiffF;
   firstProfSTPtr->avgDiffF = secProfSTPtr->avgDiffF;
   secProfSTPtr->avgDiffF = swapF;
} /*swap_prof_illNano*/

/*-------------------------------------------------------\
| Fun09: merge_prof_illNano
|   - merges two prof_illNano structs into one profile
| Input:
|   - firstProfSTPtr;
|     o pointer to first prof_illNano struct to merge
|     o should have link to list
|   - secProfSTPtr;
|     o pointer to second prof_illNano struct to merge
| Output:
|   - Modifies:
|     o firstProfSTPtr to have merged profile
|     o secProfSTPtr to have umerged values to free
| WARNING:
|   - does not free the second structure
|   - assumes one struct has all varaints in both structs
|   - will swap first and second struct if needed, but
|     will not swap nextST (list order)
\-------------------------------------------------------*/
void
merge_prof_illNano(
   struct prof_illNano *firstProfSTPtr, /*first struct*/
   struct prof_illNano *secProfSTPtr    /*second struct*/
){
   uint uiFirst = 0;
   uint uiSec = 0;

   if(
        firstProfSTPtr->varInProfUI
      < secProfSTPtr->varInProfUI
   ){ /*If: need to swap*/
      swap_prof_illNano(
         firstProfSTPtr,
         secProfSTPtr
      );
   } /*If: need to swap*/


   while(
         firstProfSTPtr->posAryUI[uiFirst]
      < secProfSTPtr->posAryUI[uiSec]
   ) ++uiFirst; /*find first overlapping position*/


   while(uiSec < secProfSTPtr->varInProfUI)
   { /*Loop: merge profiles*/
      if(uiFirst >= firstProfSTPtr->varInProfUI)
         break;

      if(
            firstProfSTPtr->ntArySC[uiFirst]
         == def_unkown_illNano
      ){ /*If: no idea about base*/
         firstProfSTPtr->ntArySC[uiFirst] =
            secProfSTPtr->ntArySC[uiSec];
      } /*If: no idea about base*/

      firstProfSTPtr->depthAryUI[uiFirst] +=
         secProfSTPtr->depthAryUI[uiSec];

      firstProfSTPtr->xDepthAryUI[uiFirst] +=
         secProfSTPtr->xDepthAryUI[uiSec];

      ++uiSec;
      ++uiFirst;
   } /*Loop: merge profiles*/
} /*merge_prof_illNano*/

/*-------------------------------------------------------\
| Fun10: blank_profList_illNano
|   - blanks a profList_illNano struct (does nothing)
| Input:
|   - profListSTPtr;
|     o pointer to profList_illNano struct to blank
| Output:
\-------------------------------------------------------*/
void
blank_profList_illNano(
   struct profList_illNano *profListSTPtr
){
   profListSTPtr=profListSTPtr; /*stop compiler warnings*/
} /*blank_profList_illNano*/

/*-------------------------------------------------------\
| Fun11: init_profList_illNano
|   - initializes a profList_illNano struct
| Input:
|   - profListSTPtr;
|     o pointer to profList_illNano struct to initialize
| Output:
|   - Modifies:
|     o all values in profListSTPtr to be 0 (null)
\-------------------------------------------------------*/
void
init_profList_illNano(
   struct profList_illNano *profListSTPtr
){
   profListSTPtr->numProfUL = 0;
   profListSTPtr->listST = 0;
} /*init_profList_illNano*/

/*-------------------------------------------------------\
| Fun12: freeStack_profList_illNano
|   - frees variables (list) inside a profList_illNano
| Input:
|   - profListSTPtr;
|     o pointer to profList_illNano with variables to free
| Output:
|   - Frees:j
|     o heap allocated (list) variables in profListSTPtr
\-------------------------------------------------------*/
void
freeStack_profList_illNano(
   struct profList_illNano *profListSTPtr
){
   if(! profListSTPtr)
      return;

   freeHeapList_prof_illNano(profListSTPtr->listST);
   init_profList_illNano(profListSTPtr);
} /*freeStack_profList_illNano*/

/*-------------------------------------------------------\
| Fun13: freeHeap_profList_illNano
|   - frees a profList_illNano struct
| Input:
|   - profListSTPtr;
|     o pointer to profList_illNano struct to free
| Output:
|   - Frees:
|     o profListSTPtr (you must set to 0 [null])
\-------------------------------------------------------*/
void
freeHeap_profList_illNano(
   struct profList_illNano *profListSTPtr
){
   if(! profListSTPtr)
      return;

   freeStack_profList_illNano(profListSTPtr);
   free(profListSTPtr);
} /*freeHeap_profList_illNano*/

/*-------------------------------------------------------\
| Fun14: getIllSnp_illNano
|   - gets snps from tbCon tsv file for Illumina data
| Input:
|   - minDepthUI:
|     o min read depth to keep varaint
|   - minPercDepthF:
|     o min percent read depth to keep varaint
|   - minPercDelF:
|     o min percent read depth to keep deletion; you
|       should use at least 50%, so deletion needs to be
|       major choice
|   - numVarUIPtr:
|     o pointer to unsigned int to hold the number of
|       variant positions found
|   - lenProfUI:
|     o pointer to unsigned int to hold length (final
|       reference position) of profile
|   - tsvFILE:
|     o FILE pointer to tsv file from tbCon
|   - errSCPtr:
|     o pointer to signed char to hold returned error
| Output:
|   - Modifies:
|     o numVarUIPtr to have number of variant positions
|     o lenProfUI to have reference position of last
|       variant
|     o tsvFILE to point to end of file (EOF)
|     o errSCPtr:
|       - 0 for no errors
|       - def_memErr_illNano for memory errors
|       - def_fileErr_illNano for file errors
|   - Returns:
|     o signed char array with variants
|       o def_A_illNano if A supported
|       o def_C_illNano if C supported
|       o def_G_illNano if G supported
|       o def_T_illNano if T supported
|       o def_del_illNano if del supported
|       o def_noBase_illNano if no variant for position
|       o '\0' for end of variants (no more variants)
|     o 0 for memory error
\-------------------------------------------------------*/
signed char *
getIllSnp_illNano(
   unsigned int minDepthUI, /*min depth to keep snp*/
   float minPercDepthF,     /*min % depth to keep snp*/
   float minPercDelF,       /*min % depth to keep del*/
   unsigned int *numVarUIPtr,/*number variants found*/
   unsigned int *lenProfUI, /*length of final profile*/
   void *tsvFILE,           /*tsv with snps*/
   signed char *errSCPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun14 TOC:
   '   - gets snps from tbCon tsv file for Illumina data
   '   o fun14 sec01:
   '     - variable declarations
   '   o fun14 sec02:
   '     - allocate profile memory (initial)
   '   o fun14 sec03:
   '     - read in variants
   '   o fun14 sec04:
   '     - remove positons with no varaints (one base)
   '   o fun14 sec05:
   '     - clean up and return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define lenBuff_fun14 4096
   schar buffStr[lenBuff_fun14];
   uint uiPos = 0;

   uint refPosUI = 0;
   schar baseSC = 0;
   uint depthUI = 0;
   float percDepthF = 0;
   schar *profileStr = 0;
   uint sizeProfUI = 0;
   uint tmpUI = 0;

   schar *tmpStr = 0;
   schar delBl = 0;   /*1: deletion entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec02:
   ^   - allocate profile memory (initial)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *numVarUIPtr = 0;
   *lenProfUI = 0;

   profileStr = malloc((4096 + 9) * sizeof(schar));
   sizeProfUI = 4095;

   if(! profileStr)
      goto memErr_fun14_sec05;

   for(
      tmpUI = 0;
      tmpUI <= sizeProfUI;
      ++tmpUI
   ) profileStr[tmpUI] = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec03:
   ^   - read in variants
   ^   o fun14 sec03 sub01:
   ^     - get past reference id/name + start loop
   ^   o fun14 sec03 sub02:
   ^     - get reference position
   ^   o fun14 sec03 sub03:
   ^     - check mutation type
   ^   o fun14 sec03 sub04:
   ^     - get nucleotide (if snp)
   ^   o fun14 sec03 sub05:
   ^     - get read depth (do minimum depth checks)
   ^   o fun14 sec03 sub06:
   ^     - get percent read depth (do min/max check)
   ^   o fun14 sec03 sub07:
   ^     - check if profile needs more memory
   ^   o fun14 sec03 sub08:
   ^     - assign variant to profile
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun14 Sec03 Sub01:
   *   - get past reference id/name + start loop
   \*****************************************************/

   tmpStr =
      (schar *)
      fgets(
         (char *) buffStr,
         lenBuff_fun14,
         (FILE *) tsvFILE
      ); /*get past header*/

   while(
      fgets(
         (char *) buffStr,
         lenBuff_fun14,
         (FILE *) tsvFILE
      )
   ){ /*Loop: get varaints*/
      uiPos = 0;

      while(buffStr[uiPos++] > 32) ;

      if(
            buffStr[uiPos - 1] != '\t'
         && buffStr[uiPos - 1] != ' '
      ) goto fileErr_fun14_sec05;

      while(buffStr[uiPos] < 33)
      { /*Loop: get off white space*/
         if(buffStr[uiPos] == '\0')
            goto fileErr_fun14_sec05; /*end of line*/
         ++uiPos;
      } /*Loop: get off white space*/

      /**************************************************\
      * Fun14 Sec03 Sub02:
      *   - get reference position
      \**************************************************/

      uiPos +=
         strToUI_base10str(
            &buffStr[uiPos],
            &refPosUI
         ); /*get position*/

      --refPosUI; /*convert to index 0*/

      if(
            buffStr[uiPos] != '\t'
         && buffStr[uiPos] != ' '
      ) goto fileErr_fun14_sec05;

      /*move to next entry*/

      while(buffStr[uiPos] < 33)
      { /*Loop: get off white space*/
         if(buffStr[uiPos] == '\0')
            goto fileErr_fun14_sec05; /*end of line*/
         ++uiPos;
      } /*Loop: get off white space*/

      /**************************************************\
      * Fun14 Sec03 Sub03:
      *   - check mutation type
      \**************************************************/

      /*get past type entry*/
      if(
            buffStr[uiPos] == 'i'
         && buffStr[uiPos + 1] == 'n'
         && buffStr[uiPos + 2] == 's'
      ) continue; /*is insertion entry*/

      if(
            buffStr[uiPos] == 'd'
         && buffStr[uiPos + 1] == 'e'
         && buffStr[uiPos + 2] == 'l'
      ) delBl = 1; /*deletion entry*/

      else
         delBl = 0;

      while(buffStr[uiPos++] > 32) ;

      if(
            buffStr[uiPos - 1] != '\t'
         && buffStr[uiPos - 1] != ' '
      ) goto fileErr_fun14_sec05;

      while(buffStr[uiPos] < 33)
      { /*Loop: get off white space*/
         if(buffStr[uiPos] == '\0')
            goto fileErr_fun14_sec05; /*end of line*/
         ++uiPos;
      } /*Loop: get off white space*/

      /**************************************************\
      * Fun14 Sec03 Sub04:
      *   - get nucleotide (if snp)
      \**************************************************/

      if(! delBl)
         baseSC = buffStr[uiPos]; /*nucleotide*/
      else
         baseSC = 'D';

      /*move to next entry*/

      while(buffStr[uiPos++] > 32) ;

      if(
            buffStr[uiPos - 1] != '\t'
         && buffStr[uiPos - 1] != ' '
      ) goto fileErr_fun14_sec05;

      while(buffStr[uiPos] < 33)
      { /*Loop: get off white space*/
         if(buffStr[uiPos] == '\0')
            goto fileErr_fun14_sec05; /*end of line*/
         ++uiPos;
      } /*Loop: get off white space*/

      /**************************************************\
      * Fun14 Sec03 Sub05:
      *   - get read depth (do minimum depth checks)
      \**************************************************/

      uiPos +=
         strToUI_base10str(
            &buffStr[uiPos],
            &depthUI
         ); /*get read depth*/

      if(
            buffStr[uiPos] != '\t'
         && buffStr[uiPos] != ' '
      ) goto fileErr_fun14_sec05;

      if(buffStr[uiPos] > 32)
         goto fileErr_fun14_sec05; /*end of line*/

      if(depthUI < minDepthUI)
         continue; /*variant not worth keeping*/

      /*move to next entry*/

      while(buffStr[uiPos] < 33)
      { /*Loop: get off white space*/
         if(buffStr[uiPos] == '\0')
            goto fileErr_fun14_sec05; /*end of line*/
         ++uiPos;
      } /*Loop: get off white space*/

      /**************************************************\
      * Fun14 Sec03 Sub06:
      *   - get percent read depth (do min/max check)
      \**************************************************/

      tmpUI = uiPos;

      while(buffStr[tmpUI++] > 32) ;
      --tmpUI;
      buffStr[tmpUI] = '\0';
      tmpUI = 0;

      percDepthF = atof((char *) &buffStr[uiPos]);

      if(
            ! delBl
         && percDepthF < minPercDepthF
      ) continue; /*variant not worth keeping*/

      else if(
            delBl
         && percDepthF < minPercDelF
      ) continue; /*deletion has to low support*/

      /**************************************************\
      * Fun14 Sec03 Sub07:
      *   - check if profile needs more memory
      \**************************************************/

      if(sizeProfUI < refPosUI)
      { /*If: need to resize array*/
         tmpStr =
            realloc(
               profileStr,
               ((refPosUI << 1) + 9) * sizeof(schar)
            );

         if(! tmpStr)
            goto memErr_fun14_sec05;

         profileStr = tmpStr;
         tmpStr = 0;

         for(
            tmpUI = sizeProfUI;
            tmpUI <= sizeProfUI << 1;
            ++tmpUI
         ) profileStr[tmpUI] = 0;

         sizeProfUI <<= 1;
      } /*If: need to resize array*/

      /**************************************************\
      * Fun14 Sec03 Sub08:
      *   - assign variant to profile
      \**************************************************/

      *lenProfUI = refPosUI;

      if(delBl)
      { /*If: was deletion*/
         profileStr[refPosUI] |= def_del_illNano;
         continue;
      } /*If: was deletion*/

      switch(baseSC & ~32)
      { /*Switch: find base*/
         case 'A':
            profileStr[refPosUI] |= def_A_illNano;
            break;
         case 'C':
            profileStr[refPosUI] |= def_C_illNano;
            break;
         case 'G':
            profileStr[refPosUI] |= def_G_illNano;
            break;
         case 'T':
            profileStr[refPosUI] |= def_T_illNano;
            break;
      } /*Switch: find base*/
   } /*Loop: get varaints*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec04:
   ^   - remove positons with no varaints (one base)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ++(*lenProfUI);

   for(
      refPosUI = 0;
      refPosUI < *lenProfUI;
      ++refPosUI
   ){ /*Loop: remove positons with one varaint*/
      baseSC = profileStr[refPosUI];

      /*from stanford bit twiddeling guie; Brain Kernihn's
      `  method for counting bits
      */
      for(
         depthUI = 0;
         baseSC;
         ++depthUI
      ) baseSC &= (baseSC - 1);
        /*baseSC &= (baseSC - 1) removes last set bit
           baseSC &= 0x1111 - 0x1
             goes to 0x1110
           baseSC &= 0x1110 - 0x1
             goes to 0x1101 - 0x1
             goes to 0x1100
        */ 

      if(depthUI < 2)
         profileStr[refPosUI] = def_noBase_illNano;
         /*this case is true most of the time, so
         `   branch prediction is good here
         */
      else
         ++(*numVarUIPtr); /*kept variant position*/
   } /*Loop: remove positons with one varaint*/

   profileStr[*lenProfUI] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec05:
   ^   - clean up and return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *errSCPtr = 0;
   goto ret_fun14_sec05;

   memErr_fun14_sec05:;
     *errSCPtr = def_memErr_illNano;
     goto errCleanUp_fun14_sec05;

   fileErr_fun14_sec05:;
      *errSCPtr = def_fileErr_illNano;
      goto errCleanUp_fun14_sec05;

   errCleanUp_fun14_sec05:;
      if(profileStr)
         free(profileStr);
      profileStr=0;
      goto ret_fun14_sec05;

   ret_fun14_sec05:;
      return profileStr;
} /*getIllSnp_illNano*/

/*-------------------------------------------------------\
| Fun15: getVarNano_illNano
|   - identifies variants in a nanopore read using an
|     Illumina profile
| Input:
|   - samSTPtr:
|     o pointer to samEntry stuct with read to scan
|   - profileStr:
|     o c-strint with variants (profile) to extract
|   - mapPosAryUI:
|     o unsigned int array to store variant positions in
|       (size needs to be number variants in profileStr)
|   - mapNtArySC:
|     o signed char array to store variant nucleotides in
|       (size needs to be number variants in profileStr)
|   - tranBl:
|     o 1: for unkown entries; assigned variant to a
|          transversion
|     o 0: all unidentifed entries are def_unkown_illNano
| Output:
|   - Modifies:
|     o mapPosAryUI to have variant positions
|     o mapNtArySC to have variant nucleotdies
|   - Returns:
|     o number of detected variant positions
\-------------------------------------------------------*/
unsigned int
getVarNano_illNano(  
   struct samEntry *samSTPtr,/*read to scan*/
   signed char *profileStr,  /*illumina profile*/
   signed char tranBl,       /*1: is uknown transition*/
   unsigned int *mapPosAryUI,/*will have variant pos*/
   signed char *mapNtAryUI  /*will have variant bases*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun15 TOC: getVarNano_illNano
   '   - identifies variants in a nanopore read using an
   '     Illumina profile
   '   o fun15 sec01:
   '     - allocate memory and get first line of sam file
   '   o fun15 sec02:
   '     - get first line of sam file
   '   o fun15 sec03:
   '     - scan read for variants
   '   o fun15 sec04:
   '     - return number detected variants
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint refPosUI = 0;
   uint readPosUI = 0;
   uint cigPosUI = 0;
   uint cigNtUI = 0;

   uint numVarUI = 0;     /*number of variants*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec03:
   ^   - scan read for variants
   ^   o fun15 sec03 sub01:
   ^     - prepare for scan + start loop + cigar switch
   ^   o fun15 sec03 sub02:
   ^     - scan snp/matches for variants
   ^   o fun15 sec03 sub03:
   ^     - check if have deletion variants
   ^   o fun15 sec03 sub04:
   ^     - ignore indel variants
   ^   o fun15 sec03 sub05:
   ^     - move to next cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun15 Sec03 Sub01:
   *   - prepare for scan + start loop + cigar switch
   \*****************************************************/

   refPosUI = samSTPtr->refStartUI;

   if(samSTPtr->cigTypeStr[0] == 'S')
   { /*If: first bases soft masked*/
      cigPosUI = 1;
      readPosUI = samSTPtr->cigArySI[0];
   } /*If: first bases soft masked*/

   if(samSTPtr->cigTypeStr[0] == 'H')
      cigPosUI = 1;

   while(cigPosUI < samSTPtr->lenCigUI)
   { /*Loop: go through cigar*/

      if(profileStr[refPosUI] == '\0')
         goto ret_fun15_sec04;

      switch(samSTPtr->cigTypeStr[cigPosUI])
      { /*Switch: find mutation type*/

         /***********************************************\
         * Fun15 Sec03 Sub02:
         *   - scan snp/matches for variants
         *   o fun15 sec03 sub02 cat01:
         *     - found match/snp; loop though bases
         *   o fun15 sec03 sub02 cat02:
         *     - A variant cases
         *   o fun15 sec03 sub02 cat03:
         *     - C variant cases
         *   o fun15 sec03 sub02 cat04:
         *     - G variant cases
         *   o fun15 sec03 sub02 cat05:
         *     - T variant cases
         *   o fun15 sec03 sub02 cat06:
         *     - move to next base
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun15 Sec03 Sub02 Cat01:
         +   - found match/snp; loop though bases
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'M':
         case '=':
         case 'X':
         /*Case: match or snp*/
            cigNtUI = samSTPtr->cigArySI[cigPosUI];

            while(cigNtUI > 0)
            { /*Loop: see if have variant*/
               if(
                    profileStr[refPosUI]
                  & def_noBase_illNano
               ) goto nextNt_fun15_sec03_sub02_cat06;
                  /*no variants kept for position*/

               if(profileStr[refPosUI] == '\0')
                  goto ret_fun15_sec04;
               
               mapPosAryUI[numVarUI] = refPosUI;

               switch(
                  samSTPtr->seqStr[readPosUI] & ~32
               ){ /*Switch: find sequence base*/

                  /*+++++++++++++++++++++++++++++++++++++\
                  + Fun15 Sec03 Sub02 Cat02:
                  +  - A variant cases
                  \+++++++++++++++++++++++++++++++++++++*/

                  case 'A':
                  /*Case: read has A*/
                     if(
                          profileStr[refPosUI]
                        & def_A_illNano
                     ) mapNtAryUI[numVarUI] = 'A';
                       /*if variant dected*/

                     else if(tranBl)
                     { /*Else If: unkown; check transit*/
                        if(
                             profileStr[refPosUI]
                           & def_G_illNano
                        ) mapNtAryUI[numVarUI] = 'A';
                          /*is transition of varaint*/

                        else
                           mapNtAryUI[numVarUI] =
                              def_unkown_illNano;
                     } /*Else If: unkown; check transit*/

                     else
                        mapNtAryUI[numVarUI] =
                           def_unkown_illNano;
                           /*no idea what is*/

                     break;
                  /*Case: read has A*/

                  /*+++++++++++++++++++++++++++++++++++++\
                  + Fun15 Sec03 Sub02 Cat03:
                  +  - C variant cases
                  \+++++++++++++++++++++++++++++++++++++*/

                  case 'C':
                  /*Case: read has C*/
                     if(
                          profileStr[refPosUI]
                        & def_C_illNano
                     ) mapNtAryUI[numVarUI] = 'C';
                       /*variant detected*/

                     else if(tranBl)
                     { /*Else If: unkown; check transit*/
                        if(
                             profileStr[refPosUI]
                           & def_T_illNano
                        ) mapNtAryUI[numVarUI] = 'C';
                          /*is transition of varaint*/

                        else
                           mapNtAryUI[numVarUI] =
                              def_unkown_illNano;
                     } /*Else If: unkown; check transit*/

                     else
                        mapNtAryUI[numVarUI] =
                           def_unkown_illNano;
                           /*no idea what is*/

                     break;
                  /*Case: read has C*/

                  /*+++++++++++++++++++++++++++++++++++++\
                  + Fun15 Sec03 Sub02 Cat04:
                  +  - G variant cases
                  \+++++++++++++++++++++++++++++++++++++*/

                  case 'G':
                  /*Case: read has G*/
                     if(
                          profileStr[refPosUI]
                        & def_G_illNano
                     ) mapNtAryUI[numVarUI] = 'G';
                       /*matches varaint*/

                     else if(tranBl)
                     { /*Else If: unkown; check transit*/
                        if(
                             profileStr[refPosUI]
                           & def_A_illNano
                        ) mapNtAryUI[numVarUI] = 'G';
                          /*is transition of varaint*/

                        else
                           mapNtAryUI[numVarUI] =
                              def_unkown_illNano;
                     } /*Else If: unkown; check transit*/

                     else
                        mapNtAryUI[numVarUI] =
                           def_unkown_illNano;
                           /*no idea what is*/
                     break;
                  /*Case: read has G*/

                  /*+++++++++++++++++++++++++++++++++++++\
                  + Fun15 Sec03 Sub02 Cat05:
                  +  - T variant cases
                  \+++++++++++++++++++++++++++++++++++++*/

                  case 'T':
                  /*Case: read has T*/
                     if(
                          profileStr[refPosUI]
                        & def_T_illNano
                     ){ /*If: have variant*/
                        mapNtAryUI[numVarUI] = 'T';
                     } /*If: have variant*/

                     else if(tranBl)
                     { /*Else If: unkown; check transit*/
                        if(
                             profileStr[refPosUI]
                           & def_C_illNano
                        ) mapNtAryUI[numVarUI] = 'T';
                          /*is transition of varaint*/

                        else
                           mapNtAryUI[numVarUI] =
                              def_unkown_illNano;
                     } /*Else If: unkown; check transit*/

                     else
                        mapNtAryUI[numVarUI] =
                           def_unkown_illNano;
                           /*no idea what is*/

                     break;
                  /*Case: read has T*/
               } /*Switch: find sequence base*/

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun15 Sec03 Sub02 Cat06:
               +   - move to next base
               \++++++++++++++++++++++++++++++++++++++++*/

               ++numVarUI; /*only fires if profile pos*/

               nextNt_fun15_sec03_sub02_cat06:;

               --cigNtUI;
               ++readPosUI;
               ++refPosUI;
            } /*Loop: see if have variant*/

            break;
         /*Case: match or snp*/

         /***********************************************\
         * Fun15 Sec03 Sub03:
         *   - check if have deletion variants
         \***********************************************/

         case 'D':
         /*Case: deletion*/
            cigNtUI = samSTPtr->cigArySI[cigPosUI];

            while(cigNtUI > 0)
            { /*Loop: see if have variant*/
               if(profileStr[refPosUI] == '\0')
                  goto ret_fun15_sec04;
                  /*no more variants in profile*/;

               if(
                     profileStr[refPosUI]
                  != def_noBase_illNano
               ){ /*If: have variant at position*/
                  mapPosAryUI[numVarUI] = refPosUI;

                  if(
                     profileStr[refPosUI] &def_del_illNano
                  ) mapNtAryUI[numVarUI] =
                     def_delChar_illNano;
                     /*deletion variant*/
                  else
                     mapNtAryUI[numVarUI] =
                        def_unkown_illNano;
                        /*no idea what is*/

                  ++numVarUI;
               } /*If: have variant at position*/

               ++refPosUI;
               --cigNtUI;
            } /*Loop: see if have variant*/

            break;
         /*Case: deletion*/

         /***********************************************\
         * Fun15 Sec03 Sub04:
         *   - ignore indel variants
         \***********************************************/

         case 'I':
         /*Case: insertion*/
            readPosUI += samSTPtr->cigArySI[cigPosUI];
            break;
         /*Case: insertion*/
      } /*Switch: find mutation type*/

      /**************************************************\
      * Fun15 Sec03 Sub05:
      *   - move to next cigar entry
      \**************************************************/

      ++cigPosUI;
   } /*Loop: go through cigar*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec04:
   ^   - return number detected variants
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   mapPosAryUI[numVarUI] = 0;
   mapNtAryUI[numVarUI] = '\0';
   
   ret_fun15_sec04:;
   return numVarUI;
} /*getVarNano_illNano*/

/*-------------------------------------------------------\
| Fun16: getNanoReads_illNano
|   - gets nanopore reads using Illumina profile
| Input:
|   - profileStr:
|     o c-strint with variants (profile) to extract
|   - lenProfileUI:
|     o position of last reference base in profile
|   - uniqProfSTPtr:
|     o pointer to profList_illNano struct to store unqiue
|       profile list (nothting done if 0)
|   - tranBl:
|     o 1: for unkown entries; assigned variant to a
|          transversion
|     o 0: all unidentifed entries are def_unkown_illNano
|   - mergeOverBl:
|     o 1: merge overlaps were a read does not completly
|          cover, but their is no difference
|     o 0: only merge when at least on read completly
|          overlaps other
|   - buffStrPtr:
|     o pointer to c-string to use in reading sam file
|   - lenBuffULPtr:
|     o length of buffStrPtr
|   - samFILE:
|     o FILE pointer to sam file with reads to get
|   - outFILE:
|     o FILE pointer to file to print reads to
| Output:
|   - Modifies:
|     o samSTPtr to have last sam file entry
|     o buffStrPtr to have last line in sam flie
|     o lenBuffULPtr to have buffStrPtr size (if resized)
|     o samFILE to point to end of file (EOF)
|   - Prints:
|     o read ids + profile to outFILE as tsv
|   - Returns:
|     o 0 for no errors
|     o def_memErr_illNano for memory errors
|     o def_fileErr_illNano for file errors
\-------------------------------------------------------*/
signed char
getNanoReads_illNano(  
   signed char *profileStr,  /*illumina profile*/
   unsigned int lenProfileUI,/*last ref position in prof*/
   struct profList_illNano *profListSTPtr,
   signed char tranBl,       /*1: is uknown transition*/
   signed char mergeOverBl,/*1 merge incompelte overlaps*/
   struct samEntry *samSTPtr,/*sam file reading*/
   signed char **buffStrPtr, /*buffer for file reading*/
   unsigned long *lenBuffULPtr,/*length of buffStrPtr*/
   void *samFILE,            /*sam file with ONT reads*/
   void *outFILE             /*file to save reads to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun16 TOC: getNanoReads_illNano
   '   - gets nanopore reads using Illumina profile
   '   o fun16 sec01:
   '     - allocate memory and get first line of sam file
   '   o fun16 sec02:
   '     - get first line of sam file
   '   o fun16 sec03:
   '     - get reads from sam file
   '   o fun16 sec04:
   '     - clean up and return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;

   uint *mapPosHeapAryUI = 0; /*position of mapped var*/
   schar *mapNtHeapArySC = 0; /*nucleotide of mapped var*/
   uint numVarUI = 0;     /*number of variants*/
   uint uiVar = 0;        /*for printing variants*/
   uint uiPos = 0;        /*position in output buffer*/

   struct prof_illNano *tmpHeapProfST = 0;
   struct prof_illNano *nodeSTPtr = 0; /*node in list*/

   schar mergeBl = 0;
   sint diffSI = 0;  /*# differences between profiles*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec02:
   ^   - allocate memory and get first line of sam file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   mapPosHeapAryUI = malloc(lenProfileUI * sizeof(uint));

   if(! mapPosHeapAryUI)
      goto memErr_fun16_sec04;

   mapNtHeapArySC = malloc(lenProfileUI * sizeof(schar));

   if(! mapNtHeapArySC)
      goto memErr_fun16_sec04;

   errSC =
      get_samEntry(
         samSTPtr,
         buffStrPtr,
         lenBuffULPtr,
         samFILE
       );

   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_samEntry)
        goto memErr_fun16_sec04;
      else
        goto fileErr_fun16_sec04;
   } /*If: had error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec03:
   ^   - get reads from sam file
   ^   o fun16 sec03 sub01:
   ^     - filter reads + start loop
   ^   o fun16 sec03 sub02:
   ^     - find Illumina variants in read
   ^   o fun16 sec03 sub03:
   ^     - see if profile is unique (if asked to do)
   ^   o fun16 sec03 sub04:
   ^     - print kept read (cluster is cigar like format)
   ^   o fun16 sec03 sub05:
   ^     - move to next read
   ^   o fun16 sec03 sub06:
   ^     - final error check
   ^   o fun16 Sec03 Sub07:
   ^     - find stats for profile list
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun16 Sec03 Sub01:
   *   - filter reads + start loop
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "id\tref\tnum_var\tprofile\n"
   );

   while(! errSC)
   { /*Loop: get reads*/
      if(samSTPtr->extraStr[0] == '@')
         goto nextRead_fun16_sec03_sub05;

      if(samSTPtr->flagUS & (4 | 256 | 2048))
         goto nextRead_fun16_sec03_sub05;

      /**************************************************\
      * Fun16 Sec03 Sub02:
      *   - find Illumina variants in read
      \**************************************************/

      numVarUI =
         getVarNano_illNano(
            samSTPtr,
            profileStr,
            tranBl,
            mapPosHeapAryUI,
            mapNtHeapArySC
         );

      if(numVarUI == 0)
         goto nextRead_fun16_sec03_sub05;
         /*no variant positions in read; ignore*/

      /**************************************************\
      * Fun16 Sec03 Sub03:
      *   - see if profile is unique (if asked to do)
      \**************************************************/

      if(profListSTPtr)
      { /*If: finding unique profiles*/
         nodeSTPtr = profListSTPtr->listST;

         tmpHeapProfST =
            mk_prof_illNano(
               mapPosHeapAryUI,
               mapNtHeapArySC,
               numVarUI
            ); /*make a new node*/
    
         if(! tmpHeapProfST)
            goto memErr_fun16_sec04;

         if(! nodeSTPtr)
         { /*If: first profile*/
            profListSTPtr->numProfUL = 1;
            profListSTPtr->listST = tmpHeapProfST;
            tmpHeapProfST = 0;
            goto pprofile_sec03_sub04;
         } /*If: first profile*/

         mergeBl = 0;

         while(nodeSTPtr)
         { /*Loop: compare profiles*/
            diffSI =
               cmp_prof_illNano(
                  tmpHeapProfST,
                  nodeSTPtr,
                  mergeOverBl /*overlap limits*/
               );

            if(! diffSI)
            { /*If: profiles same, or one covers other*/
               merge_prof_illNano(
                  nodeSTPtr,  /*want to merge into*/
                  tmpHeapProfST
               );

               freeHeap_prof_illNano(tmpHeapProfST);
               tmpHeapProfST = 0;
               mergeBl = 1;
               break;
            } /*If: profiles same, or one covers other*/

            if(! nodeSTPtr->nextST)
               break;

            nodeSTPtr = nodeSTPtr->nextST;
         } /*Loop: compare profiles*/

         if(! mergeBl)
         { /*If: new profile is unique*/
            ++profListSTPtr->numProfUL;
            nodeSTPtr->nextST = tmpHeapProfST;
            tmpHeapProfST->lastST = nodeSTPtr;
         } /*If: new profile is unique*/

         tmpHeapProfST = 0;
      } /*If: finding unique profiles*/

      /**************************************************\
      * Fun16 Sec03 Sub04:
      *   - print kept read (cluster is cigar like format)
      \**************************************************/

      pprofile_sec03_sub04:;

      fprintf(
         (FILE *) outFILE,
         "%s\t%s\t%u\t",
         samSTPtr->qryIdStr,
         samSTPtr->refIdStr,
         numVarUI
      );

      uiPos = 0;

      for(
         uiVar = 0;
         uiVar < numVarUI;
         ++uiVar
      ){ /*Loop: print out each variant position*/
         if(uiVar > 0)
            (*buffStrPtr)[uiPos++] = '_'; /*for humans*/
         
         uiPos +=
            numToStr(
               &(*buffStrPtr)[uiPos],
               mapPosHeapAryUI[uiVar] + 1  /*to index 1*/
            );
               
         (*buffStrPtr)[uiPos++] = mapNtHeapArySC[uiVar];
      } /*Loop: print out each variant position*/

      (*buffStrPtr)[uiPos++] = '\n';
      (*buffStrPtr)[uiPos] = '\0';

      fprintf(
         (FILE *) outFILE,
         "%s",
         *buffStrPtr
      ); /*add new line ot end of line*/

      /**************************************************\
      * Fun16 Sec03 Sub05:
      *   - move to next read
      \**************************************************/

      nextRead_fun16_sec03_sub05:;

      errSC =
         get_samEntry(
            samSTPtr,
            buffStrPtr,
            lenBuffULPtr,
            samFILE
          );
   } /*Loop: get varaints*/

   /*****************************************************\
   * Fun16 Sec03 Sub06:
   *   - move to next read
   \*****************************************************/

   if(errSC == def_memErr_samEntry)
      goto memErr_fun16_sec04;

   /*****************************************************\
   * Fun16 Sec03 Sub07:
   *   - find stats for profile list
   *   o fun16 sec03 sub07 cat01:
   *     - do scan to merge overlapping profiles (n^2)
   *   o fun16 sec03 sub07 cat02:
   *     - blank stats
   *   o fun16 sec03 sub07 cat03:
   *     - get differences (n^2 loop)
   *   o fun16 sec03 sub07 cat04:
   *     - find means
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun16 Sec03 Sub07 Cat01:
   +   - do scan to merge overlapping profiles (n^2)
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(profListSTPtr)
   { /*If: building profile list*/

      nodeSTPtr = profListSTPtr->listST;
      tmpHeapProfST = nodeSTPtr->nextST;

      while(tmpHeapProfST)
      { /*Loop: finish merging*/
         diffSI =
            cmp_prof_illNano(
               nodeSTPtr,
               tmpHeapProfST,
               mergeOverBl /*overlap limits*/
            );

            if(! diffSI)
            { /*If: merging old nodes*/
               --profListSTPtr->numProfUL;

               merge_prof_illNano(
                  nodeSTPtr,  /*want to merge into*/
                  tmpHeapProfST
               );

                tmpHeapProfST->lastST->nextST =
                   tmpHeapProfST->nextST;

               if(tmpHeapProfST->nextST)
                  tmpHeapProfST->nextST->lastST =
                     tmpHeapProfST->lastST;

               freeHeap_prof_illNano(tmpHeapProfST);
               goto nextNode_fun16_sec03_usb07_cat02;
            } /*If: merging old nodes*/

         if(tmpHeapProfST->nextST)
            tmpHeapProfST = tmpHeapProfST->nextST;
         else
         { /*Else: moving to next node*/
            nextNode_fun16_sec03_usb07_cat02:;

            nodeSTPtr = nodeSTPtr->nextST;
            tmpHeapProfST = nodeSTPtr->nextST;

            if(! tmpHeapProfST)
               break; /*no next node*/
         } /*Else: moving to next node*/
      } /*Loop: finish merging*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun16 Sec03 Sub07 Cat02:
      +   - blank stats
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      nodeSTPtr = profListSTPtr->listST;

      while(nodeSTPtr)
      { /*Loop: blank stats*/
         nodeSTPtr->overlapUI = 0;

         nodeSTPtr->minDiffUI =
            ((uint) 1 << (((sizeof(uint)) << 3) - 1)) - 1;
            /*maximum value possible for sint*/

         nodeSTPtr->maxDiffUI = 0;

         nodeSTPtr->sumDiffUL = 0;
         nodeSTPtr->avgDiffF = 0;

         nodeSTPtr = nodeSTPtr->nextST;
      } /*Loop: blank stats*/


      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun16 Sec03 Sub07 Cat03:
      +   - get differences (n^2 loop)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      nodeSTPtr = profListSTPtr->listST;
      tmpHeapProfST = nodeSTPtr->nextST;

      /*n^2 loop*/
      while(tmpHeapProfST)
      { /*Loop: compare profiles*/

         diffSI =
            cmp_prof_illNano(
               nodeSTPtr,
               tmpHeapProfST,
               1          /*want difference*/
            );

         if(diffSI >= 0)
         { /*If: profiles overlap*/
            ++nodeSTPtr->overlapUI;
            ++tmpHeapProfST->overlapUI;

            nodeSTPtr->sumDiffUL += diffSI;
            tmpHeapProfST->sumDiffUL += diffSI;


            nodeSTPtr->minDiffUI =
               min_genMath(
                  nodeSTPtr->minDiffUI,
                  (uint) diffSI
               );

            tmpHeapProfST->minDiffUI =
               min_genMath(
                  tmpHeapProfST->minDiffUI,
                  (uint) diffSI
               );


            nodeSTPtr->maxDiffUI =
               max_genMath(
                  nodeSTPtr->maxDiffUI,
                  (uint) diffSI
               );

            tmpHeapProfST->maxDiffUI =
               max_genMath(
                  tmpHeapProfST->maxDiffUI,
                  (uint) diffSI
               );
         } /*If: profiles overlap*/

         if(tmpHeapProfST->nextST)
            tmpHeapProfST = tmpHeapProfST->nextST;
         else
         { /*Else: moving to next node*/
            nodeSTPtr = nodeSTPtr->nextST;
            tmpHeapProfST = nodeSTPtr->nextST;

            if(! tmpHeapProfST)
               break; /*no next node*/
         } /*Else: moving to next node*/
      } /*Loop: compare profiles*/

      tmpHeapProfST = 0;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun16 Sec03 Sub07 Cat04:
      +   - find means
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      nodeSTPtr = profListSTPtr->listST;

      while(nodeSTPtr->nextST)
      { /*Loop: blank stats*/
         if(! nodeSTPtr->overlapUI)
           nodeSTPtr->minDiffUI = 0;

         else
         { /*Else: have overlaping profiles*/
            nodeSTPtr->avgDiffF = (float)
                 (float) nodeSTPtr->sumDiffUL
               / (float) (nodeSTPtr->overlapUI);
         } /*Else: have overlaping profiles*/

         nodeSTPtr = nodeSTPtr->nextST;
      } /*Loop: blank stats*/
   } /*If: building profile list*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec04:
   ^   - clean up and return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun16_sec04;

   memErr_fun16_sec04:;
      errSC = def_memErr_illNano;
      goto ret_fun16_sec04;

   fileErr_fun16_sec04:;
      errSC = def_fileErr_illNano;
      goto ret_fun16_sec04;

   ret_fun16_sec04:;
      if(tmpHeapProfST)
         freeHeap_prof_illNano(tmpHeapProfST);
      tmpHeapProfST = 0;

      if(mapPosHeapAryUI)
         free(mapPosHeapAryUI);
      mapPosHeapAryUI = 0;

      if(mapNtHeapArySC)
         free(mapNtHeapArySC);
      mapNtHeapArySC = 0;

      return errSC;
} /*getNanoReads_illNano*/

/*-------------------------------------------------------\
| Fun17: phead_profList_illNano
|   - prints header for unique profile tsv
| Input
|   - outFILE:
|     o FILE pointer to print unique profiler header to
| Output:
|   - Prints:
|     o header (unique profile) to outFILE
\-------------------------------------------------------*/
void
phead_profList_illNano(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "num_var\tnum_overlap\tmin_diff\tavg_diff"
   );

   fprintf(
      (FILE *) outFILE,
      "\tmax_diff\tprofile\tvar_depth\tx_var_depth\n"
   );
} /*phead_profList_illNano*/

/*-------------------------------------------------------\
| Fun18: p_profList_illNano
|   - prints unique profiles in a profList struct
| Input
|   - profListSTPtr:
|     o pointer to profListSTPtr with profiles to print
|   - outFILE:
|     o FILE pointer to print unique profiler header to
| Output:
|   - Prints:
|     o unique profiles to outFILE
\-------------------------------------------------------*/
void
p_profList_illNano(
   struct profList_illNano *profListSTPtr,
   void *outFILE
){
   uint uiVar = 0;
   struct prof_illNano *nodeST = 0;

   nodeST = profListSTPtr->listST;

   while(nodeST)
   { /*Loop: print profile list*/
      fprintf(
         (FILE *) outFILE,
         "%u\t%u\t%u\t%0.3f\t%u\t",
         nodeST->varInProfUI,
         nodeST->overlapUI,
         nodeST->minDiffUI,
         nodeST->avgDiffF,
         nodeST->maxDiffUI
      ); /*print out stats*/

      for(
         uiVar = 0;
         uiVar < nodeST->varInProfUI;
         ++uiVar
      ){ /*Loop: print profile*/
         if(uiVar)
            fprintf(
               (FILE *) outFILE,
               "_"
            ); /*separator for new item*/

         fprintf(
            (FILE *) outFILE,
            "%u%c",
            nodeST->posAryUI[uiVar],
            nodeST->ntArySC[uiVar]
         ); /*print out variant*/
      } /*Loop: print profile*/

      fprintf(
         (FILE *) outFILE,
         "\t"
      ); /*tab for depth entry*/

      for(
         uiVar = 0;
         uiVar < nodeST->varInProfUI;
         ++uiVar
      ){ /*Loop: print depth for profile*/
         if(uiVar)
            fprintf(
               (FILE *) outFILE,
               "_"
            ); /*separator for new item*/

         fprintf(
            (FILE *) outFILE,
            "%u:%u",
            nodeST->posAryUI[uiVar],
            nodeST->depthAryUI[uiVar]
         ); /*print out variant*/
      } /*Loop: print depth for profile*/

      fprintf(
         (FILE *) outFILE,
         "\t"
      ); /*tab for depth entry*/

      for(
         uiVar = 0;
         uiVar < nodeST->varInProfUI;
         ++uiVar
      ){ /*Loop: print unkown depth for profile*/
         if(uiVar)
            fprintf(
               (FILE *) outFILE,
               "_"
            ); /*separator for new item*/

         fprintf(
            (FILE *) outFILE,
            "%u:%u",
            nodeST->posAryUI[uiVar],
            nodeST->xDepthAryUI[uiVar]
         ); /*print out variant*/
      } /*Loop: print unkown depth for profile*/

      fprintf(
         (FILE *) outFILE,
         "\n"
      );

      nodeST = nodeST->nextST;
   } /*Loop: print profile list*/
} /*p_profList_illNano*/

/*-------------------------------------------------------\
| Fun19: run_illNano
|   - runs the illNano algorithim
| Input:
|   - minDepthUI:
|     o minimum depth to keep illumina variant
|   - minPercDepthF:
|     o minimum percent read depth to keep Illumina SNP
|   - minPercDelF:
|     o min percent read depth to keep Illumina deletion
|   - tranBl:
|     o 1: for unkown entries; assigned variant to a
|          transversion
|     o 0: all unidentifed entries are def_unkown_illNano
|   - mergeOverBl:
|     o 1: merge overlaps were a read does not completly
|          cover, but their is no difference
|     o 0: only merge when at least on read completly
|          overlaps other
|   - illTsvFILE:
|     o FILE pointer to tbCon tsv file with Illumina vars
|   - samFILE:
|     o FILE pointer to sam file with reads to get
|   - outFILE:
|     o FILE pointer to file to print reads to
|   - uniqFILE:
|     o FILE pointer to file to print uniqe profiles to
| Output:
|   - Modifies:
|     o illTsvFILE to point to end of file (EOF)
|     o samFILE to point to end of file (EOF)
|   - Prints:
|     o read ids + profile to outFILE as tsv
|     o unique profiels to uniqFILE (if provided)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_illNano for memory errors
|     o def_fileErr_illNano for file errors
\-------------------------------------------------------*/
signed char
run_illNano(  
   unsigned int minDepthUI,  /*min depth for illumina*/
   float minPercDepthF,      /*min % depth for ill snp*/
   float minPercDelF,        /*min % depth for ill del*/
   signed char tranBl,       /*1: is uknown transition*/
   signed char mergeOverBl,/*1 merge incompelte overlaps*/
   void *illTsvFILE,         /*tbCon tsv file for ill*/
   void *samFILE,            /*sam file with ONT reads*/
   void *outFILE,            /*file to save reads to*/
   void *uniqFILE            /*file for unique profiles*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun19 TOC:
   '   - runs the illNano algorithim
   '   o fun19 sec01:
   '     - variable declarations
   '   o fun19 sec02:
   '     - set up (initialize and memory allocation)
   '   o fun19 sec03:
   '     - get illumina profile
   '   o fun19 sec04:
   '     - find reads with variant positions in profile
   '   o fun19 sec05:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;

   schar *profHeapStr = 0;
   uint lenProfUI = 0;
   uint varPosUI = 0;  /*number of variant positions*/

   struct samEntry samStackST;
   schar *buffHeapStr = 0;
   ulong lenBuffUL = 0;

   struct profList_illNano profListStackST;
   struct profList_illNano *profListSTPtr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec02:
   ^   - set up (initialize and memory allocation)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_samEntry(&samStackST);
   init_profList_illNano(&profListStackST);

   errSC = setup_samEntry(&samStackST);

   if(errSC)
      goto memErr_fun19_sec05;

   if(uniqFILE)
      profListSTPtr = &profListStackST;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec03:
   ^   - get illumina profile
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   profHeapStr =
      getIllSnp_illNano(
         minDepthUI,
         minPercDepthF,
         minPercDelF,
         &varPosUI,
         &lenProfUI,
         illTsvFILE,
         &errSC
      ); /*get profile for Illumina varaints*/

   if(errSC)
   { /*If: have error*/
      if(errSC == def_memErr_illNano)
         goto memErr_fun19_sec05;
      else
         goto fileErr_fun19_sec05;
   } /*If: have error*/

   if(! varPosUI)
      goto noVars_fun19_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec04:
   ^   - find reads with variant positions in profile
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC =
      getNanoReads_illNano(
         profHeapStr,
         lenProfUI,
         profListSTPtr,
         tranBl,
         mergeOverBl,
         &samStackST,
         &buffHeapStr,
         &lenBuffUL,
         samFILE,
         outFILE
      ); /*extract reads with variant positions*/

   profListSTPtr = 0;

   if(errSC == def_memErr_illNano)
      goto memErr_fun19_sec05;

   if(errSC == def_fileErr_illNano)
      goto fileErr_fun19_sec05;

   if(uniqFILE)
   { /*If: printing unique profile ids*/
      phead_profList_illNano(uniqFILE);

      p_profList_illNano(
         &profListStackST,
         uniqFILE
      );
   } /*If: printing unique profile ids*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec05:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun19_sec05;

   memErr_fun19_sec05:;
      errSC = def_memErr_illNano;
      goto ret_fun19_sec05;

   fileErr_fun19_sec05:;
      errSC = def_fileErr_illNano;
      goto ret_fun19_sec05;

   noVars_fun19_sec05:;
      errSC = def_noVar_illNano;
      goto ret_fun19_sec05;

   ret_fun19_sec05:;
      if(profHeapStr)
         free(profHeapStr);
      profHeapStr = 0;
      lenProfUI = 0;

      if(buffHeapStr)
         free(buffHeapStr);
      buffHeapStr = 0;
      lenBuffUL = 0;

      freeStack_samEntry(&samStackST);
      freeStack_profList_illNano(&profListStackST);

      return errSC;
} /*run_illNano*/

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
