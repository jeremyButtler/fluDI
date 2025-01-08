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
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "edClust.h"

#include <stdio.h>

#include "../genBio/samEntry.h"
#include "../genBio/tbCon.h"
#include "../genBio/edDist.h"

#include "clustST.h"

/*.h files only (no .c files)*/
#include "../genLib/dataTypeShortHand.h"
#include "../genLib/genMath.h" /*only using .h min macro*/

#include "../genBio/tbConDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../genLib/base10str.h"
!   - .c  #include "../genLib/numToStr.h"
!   - .c  #include "../genLib/ulCp.h"
!   - .c  #include "../genLib/charCp.h"
!   - .c  #include "../genLib/strAry.h"
!   - .h  #include "../genBio/ntTo5Bit.h"
!   - .h  #include "../genBio/tbConDefs.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - finds read with best score from score array
   '   o fun01 sec01:
   '     - variable declarations
   '   o fun01 sec02:
   '     - find highest score
   '   o fun01 sec03:
   '     - get highest scoring read
   '   o fun01 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   sint bestScoreSI = 0;
   slong entrySL = 0;
   slong bestIndexSL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - find highest score
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*I pay O(n) cost here, but O(1) cost for setting
   `   reads to -1 (discarded) later. Also, allows me to
   `   use the score index as the line number in the sam
   `   file.
   */

   for(
      entrySL = 0;
      entrySL < (slong) indexSTPtr->lenSamUL;
      ++entrySL
   ){ /*Loop: find best score*/

      if(indexSTPtr->clustArySI[entrySL])
         continue; /*not read or read already assigned*/

      if(bestScoreSI < indexSTPtr->scoreArySC[entrySL])
      { /*If: new best score*/
         bestScoreSI = indexSTPtr->scoreArySC[entrySL];
         bestIndexSL = entrySL;
      } /*If: new best score*/

   } /*Loop: find best score*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - get highest scoring read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(indexSTPtr->clustArySI[bestIndexSL])
      goto noReadErr_fun01_sec04; /*no reads in file*/
      
   fseek(
      samFILE,
      indexSTPtr->indexAryUL[bestIndexSL],
      SEEK_SET
   );

   errSC =
      getRead_clustST(
         samSTPtr,
         (unsigned long) bestIndexSL,
         0,
         indexSTPtr,
         buffStr,
         samFILE
      );
      
   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_clustST)
         goto memErr_fun01_sec04;
      if(errSC == def_fileErr_clustST)
         goto fileErr_fun01_sec04;
   } /*If: had error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun01_sec04;


   memErr_fun01_sec04:;
      bestIndexSL = def_memErr_edClust;
      goto ret_fun01_sec04;


   fileErr_fun01_sec04:;
      bestIndexSL = def_fileErr_edClust;
      goto ret_fun01_sec04;


   noReadErr_fun01_sec04:;
      bestIndexSL = def_noReads_edClust;
      goto ret_fun01_sec04;


   ret_fun01_sec04:;
      fseek(
         samFILE,
         0,
         SEEK_SET
      );

      return bestIndexSL;
} /*getBestRead_edClust*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - find the depth profile for a read
   '   o fun02 sec01:
   '     - varaible declarations
   '   o fun02 sec02:
   '     - allocate memory and get start of file
   '   o fun02 sec03:
   '     - build depth profile
   '   o fun02 sec04:
   '     - return results
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - varaible declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;

   ulong offsetUL = 0;
   ulong lineUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - allocate memory and get start of file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
         ! resEdDistSTPtr->depthAryUI
      || resEdDistSTPtr->sizeDepthUL < refSTPtr->readLenUI
   ){ /*If: need more memory*/

      if(resEdDistSTPtr->depthAryUI)
         free(resEdDistSTPtr->depthAryUI);
      resEdDistSTPtr->depthAryUI = 0;

      resEdDistSTPtr->depthAryUI =
         calloc(
            refSTPtr->readLenUI + 9,
            sizeof(uint)
         ); /*allocate memory for depth array*/

      if(! resEdDistSTPtr->depthAryUI)
         goto memErr_fun02_sec04;
   } /*If: need more memory*/

   fseek(
      samFILE,
      0,
      SEEK_SET
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - build depth profile
   ^   o fun02 sec03 sub01:
   ^     - filter reads
   ^   o fun02 sec03 sub02:
   ^     - get read and add to profile
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec03 Sub01:
   *   - filter reads
   \*****************************************************/

   for(
      lineUL = 0;
      lineUL < indexSTPtr->lenSamUL;
      ++lineUL
   ){ /*Loop: add reads to depth profile*/

      if(indexSTPtr->clustArySI[lineUL] == 0) ;

      else if(indexSTPtr->clustArySI[lineUL] == clustSI) ;

      else
      { /*Else: already assigned or discarded*/
         offsetUL += indexSTPtr->lenLineAryUI[lineUL];
         continue;
      } /*Else: already assigned or discarded*/

      if(refNumUI != indexSTPtr->refAryUI[lineUL])
      { /*If: references are different*/
         offsetUL += indexSTPtr->lenLineAryUI[lineUL];
         continue;
      } /*If: references are different*/

      if(
           indexSTPtr->startAryUI[lineUL]
         < refSTPtr->refStartUI
      ){ /*If: read has bases before reference*/
         offsetUL += indexSTPtr->lenLineAryUI[lineUL];
         continue;
      } /*If: read has bases before reference*/

      if(indexSTPtr->endAryUI[lineUL] >refSTPtr->refEndUI)
      { /*If: read has bases after references*/
         offsetUL += indexSTPtr->lenLineAryUI[lineUL];
         continue;
      } /*If: read has bases after references*/

      /**************************************************\
      * Fun02 Sec03 Sub02:
      *   - get read and add to profile
      \**************************************************/

      errSC =
         getRead_clustST(
            samSTPtr,
            lineUL,
            offsetUL,
            indexSTPtr,
            buffStr,
            samFILE
         ); /*get read from sam file*/

      offsetUL = 0;

      if(errSC)
      { /*If: had error*/
         if(errSC == def_memErr_clustST)
            goto memErr_fun02_sec04;
         if(errSC == def_fileErr_clustST)
            goto fileErr_fun02_sec04;
      } /*If: had error*/


      addReadToDepth_edDist(
         refSTPtr,
         samSTPtr,
         clustSetSTPtr->minSnpQUC,
         -1,                   /*already checked overlap*/
         resEdDistSTPtr
      );
   } /*Loop: add reads to depth profile*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^   - return results
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun02_sec04;

   memErr_fun02_sec04:;
      errSC = def_memErr_edClust;
      goto ret_fun02_sec04;

   fileErr_fun02_sec04:;
      errSC = def_fileErr_edClust;
      goto ret_fun02_sec04;

   ret_fun02_sec04:;
      fseek(
         samFILE,
         0,
         SEEK_SET
      );

      return errSC;
} /*depthProf_edClust*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - finds number of reads mapping to best read
   '   o fun03 sec01:
   '     - variable declarations
   '   o fun03 sec02:
   '     - get first entry in sam file
   '   o fun03 sec03:
   '     - find number of hits
   '   o fun03 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   slong distSL = 0;
   slong lineSL = 0;

   struct conNt_tbCon *conHeapAryST = 0;
   uint firstBaseUI = 0;   /*offset for start*/
   uint lastBaseUI = 0; /*offset for start*/

   float errRateF = 0;

   ulong offsetUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - get first entry in sam file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *numReadsSL = 0;
   firstBaseUI = bestSTPtr->refStartUI;
   lastBaseUI = bestSTPtr->refEndUI;
   *lenConUI = bestSTPtr->alnReadLenUI;

   conHeapAryST = 
      malloc(
        (*lenConUI) * sizeof(struct conNt_tbCon)
      );

   if(! conHeapAryST)
      goto memErr_fun03_sec04;

   for(
      distSL = 0;
      distSL < (slong) *lenConUI;
      ++distSL
   ) init_conNt_tbCon(&conHeapAryST[distSL]);
     /*initialize conensus array*/

   if(conBl)
      errRateF = clustSetSTPtr->conErrRateF;
   else
   { /*Else: read input for best read*/
      errRateF = clustSetSTPtr->readErrRateF;

   } /*Else: read input for best read*/

   clustSetSTPtr->startUI = 0;
   clustSetSTPtr->endUI = 0;

   fseek(
      samFILE,
      0,
      SEEK_SET
   );

   if(clustSetSTPtr->depthProfBl)
   { /*If: doing depth profiling*/
      errSC =
         depthProf_edClust(
            bestSTPtr,
            indexSTPtr->refAryUI[bestIndexSL],
            clustSetSTPtr->clustSI,
            indexSTPtr,
            clustSetSTPtr,
            resEdDistSTPtr,
            samSTPtr,
            buffStr,
            samFILE
         ); /*get read depth profile*/

      if(errSC)
      { /*If: had error*/
         if(errSC == def_memErr_edClust)
            goto memErr_fun03_sec04;
         if(errSC == def_fileErr_edClust)
            goto fileErr_fun03_sec04;
      } /*If: had error*/
   } /*If: doing depth profiling*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - find number of hits
   ^   o fun03 sec03 sub01:
   ^     - filter reads + start loop
   ^   o fun03 sec03 sub02:
   ^     - get read
   ^   o fun03 sec03 sub03:
   ^     - find edit distance and decide if keep
   ^   o fun03 sec03 sub04:
   ^     - move to next read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - filter reads + start loop
   \*****************************************************/

   for(
      lineSL = 0;
      lineSL < (slong) indexSTPtr->lenSamUL;
      ++lineSL
   ){ /*Loop: find number of matches*/
      if(indexSTPtr->clustArySI[lineSL] == 0) ;

      else if(indexSTPtr->clustArySI[lineSL] == clustSI) ;

      else
      { /*Else: already assigned or discarded*/
         offsetUL += indexSTPtr->lenLineAryUI[lineSL];
         continue;
      } /*Else: already assigned or discarded*/

      if(
            indexSTPtr->refAryUI[bestIndexSL]
         != indexSTPtr->refAryUI[lineSL]
      ) goto noMatch_fun03_sec03_sub01;
        /*different references*/

      if(
           indexSTPtr->startAryUI[lineSL]
         < firstBaseUI
      ) goto noMatch_fun03_sec03_sub01;
        /*starts before the best read*/

      if(
           indexSTPtr->endAryUI[lineSL]
         > lastBaseUI
      ){ /*If: goes out of best reads bounds*/
         noMatch_fun03_sec03_sub01:;

         indexSTPtr->clustArySI[lineSL] = 0;
         offsetUL += indexSTPtr->lenLineAryUI[lineSL];
         continue;
      } /*If: goes out of best reads bounds*/

      /**************************************************\
      * Fun03 Sec03 Sub02:
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
      { /*If: had error*/
         if(errSC == def_memErr_clustST)
            goto memErr_fun03_sec04;
         if(errSC == def_fileErr_clustST)
            goto fileErr_fun03_sec04;
      } /*If: had error*/

      /**************************************************\
      * Fun03 Sec03 Sub02:
      *   - find edit distance and decide if keep
      \**************************************************/

      distSL =
         readCmpDist_edDist(
             samSTPtr,
             bestSTPtr,
             clustSetSTPtr->indelLenUI,
             clustSetSTPtr->minSnpQUC,
             -1,
             clustSetSTPtr->minDepthUI,
             clustSetSTPtr->winSizeUI,
             resEdDistSTPtr
         );

      if(distSL < 0)
      { /*If: distance to small*/
         indexSTPtr->clustArySI[lineSL] = 0;
         continue;
      } /*If: distance to small*/

      distSL =
         percDist_edDist(
            resEdDistSTPtr,
            clustSetSTPtr->winSizeUI,
            errRateF
         ); /*find variant:error ratio*/

      if(distSL > clustSetSTPtr->varToErrUI)
      { /*If: reads are to distant*/
         indexSTPtr->clustArySI[lineSL] = 0;
         continue;
      } /*If: reads are to distant*/

      if(
               samSTPtr->alnReadLenUI
            >= clustSetSTPtr->winSizeUI
         &&
              resEdDistSTPtr->probMaxWinUI
            > clustSetSTPtr->winErrUI
      ){ /*If: reads are to distant*/
         indexSTPtr->clustArySI[lineSL] = 0;
         continue;
      } /*If: reads are to distant*/

      indexSTPtr->clustArySI[lineSL] = clustSI;
      ++(*numReadsSL);

      if(clustSetSTPtr->startUI == 0)
         clustSetSTPtr->startUI = samSTPtr->refStartUI;

      clustSetSTPtr->startUI =
         min_genMath(
            samSTPtr->refStartUI,
            clustSetSTPtr->startUI
         );

      clustSetSTPtr->endUI =
         max_genMath(
            samSTPtr->refEndUI,
            clustSetSTPtr->endUI
         );

      /*set offsets for alignment*/
      samSTPtr->refStartUI -= firstBaseUI;
      samSTPtr->refEndUI -= firstBaseUI;

      errSC =
         addRead_tbCon(
            samSTPtr,
            &conHeapAryST,
            lenConUI,
            conSetSTPtr
         ); /*build the conensus*/

      if(errSC)
         goto memErr_fun03_sec04;
   } /*Loop: find number of matches*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun03_sec04;

   memErr_fun03_sec04:;
      *numReadsSL = def_memErr_edClust;
      goto errCleanUp_fun03_sec04;

   fileErr_fun03_sec04:;
     *numReadsSL = def_fileErr_edClust;
     goto errCleanUp_fun03_sec04;

   errCleanUp_fun03_sec04:;
      if(conHeapAryST)
         freeHeapAry_conNt_tbCon(
            conHeapAryST,
            *lenConUI
         );

      conHeapAryST = 0;
      *lenConUI = 0;
      goto ret_fun03_sec04;

   ret_fun03_sec04:;
      fseek(
         samFILE,
         0,
         SEEK_SET
      );

   return conHeapAryST;
} /*findNumMap_edClust*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun04 TOC:
   '   - gets number of reads that overlap a cluster
   '   o fun04 sec01:
   '     - variable declarations
   '   o fun04 sec02:
   '     - find coordinates for cluster + start main loop
   '   o fun04 sec03:
   '     - find number reads within coordinates
   '   o fun04 sec04:
   '     - move to next cluster
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint startUI = (uint) -1; /*maximum value for uint*/
   uint endUI = 0;
   ulong numReadsUL = 0;
   ulong indexUL = 0;

   float percOverlapF = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec02:
   ^   - find coordinates for cluster + start main loop
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(conListSTPtr)
   { /*Loop: while have clusters*/

      startUI = (uint) -1; /*maximum value for uint*/
      endUI = 0;
      numReadsUL = 0;

      for(
         indexUL = 0;
         indexUL < indexSTPtr->lenSamUL;
         ++indexUL
      ){ /*Loop: find start and end coordinates*/
         if(
               conListSTPtr->clustSI
            == indexSTPtr->clustArySI[indexUL]
         ){ /*If: reads in same cluster*/
            startUI =
               min_genMath(
                  startUI,
                  indexSTPtr->startAryUI[indexUL]
               );

            endUI =
               max_genMath(
                  endUI,
                  indexSTPtr->endAryUI[indexUL]
               );

            ++numReadsUL;
         } /*If: reads in same cluster*/
      } /*Loop: find start and end coordinates*/

      if(! numReadsUL)
         goto nextClust_fun04_sec04;

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun04 Sec03:
      ^   - find number reads within coordinates
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      numReadsUL = 0;

      for(
         indexUL = 0;
         indexUL < indexSTPtr->lenSamUL;
         ++indexUL
      ){ /*Loop: find number overlaping reads*/
         if(indexSTPtr->clustArySI[indexUL] <= 0)
            continue; /*read not in cluster*/

         /*find overlap length (read overlapping coords)*/

         percOverlapF =
            (float)
            min_genMath(
               endUI,
               indexSTPtr->endAryUI[indexUL]
            );

         percOverlapF -=
            (float)
            max_genMath(
               startUI,
               indexSTPtr->startAryUI[indexUL]
            );

         /*find % of read in coordinates*/

         percOverlapF =
              (float)
                 (
                     indexSTPtr->endAryUI[indexUL]
                   - indexSTPtr->startAryUI[indexUL]
                 )
            / percOverlapF; /*get perce overlap*/

         /*check if read was within minimum overlap*/

         if(percOverlapF >= minOverlapF)
            ++numReadsUL;
      } /*Loop: find number overlaping reads*/
      
      conListSTPtr->maxReadsUL = numReadsUL;

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun04 Sec04:
      ^   - move to next cluster
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      nextClust_fun04_sec04:;
         conListSTPtr = conListSTPtr->nextST;
   } /*Loop: while have clusters*/
} /*getMaxDepth_edClust*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun05 TOC:
   '   - cluster reads
   '   o fun05 sec01:
   '     - variable declarations
   '   o fun05 sec02:
   '     - initialize and set up variables
   '   o fun05 sec03:
   '     - set up index array (and score)
   '   o fun05 sec04:
   '     - cluster reads
   '   o fun05 sec05:
   '     - final consensus comparison
   '   o fun05 sec06:
   '     - remove low read depth clusters
   '   o fun05 sec07:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;         /*for error checking*/

   slong numReadsSL = 0;    /*number of reads left*/
   slong lastReadsSL = 0;   /*last number of reads left*/
   ulong clustReadsUL = 0;  /*number reads clustered*/
   slong depthSL = 0;       /*read depth for cluster*/
    
   uchar ucCon = 0;    /*iterator for consensus building*/

   struct samEntry bestSamStackST;/*read to cluster with*/
   slong bestIndexSL = 0;         /*index of best read*/

   struct samEntry *conSamHeapST = 0; /*for consensus*/

   struct conNt_tbCon *conHeapAryST = 0;
   uint lenConUI = 0;
   uint numMaskUI = 0; /*throw away*/

   struct con_clustST *retHeapST = 0;
   struct con_clustST *conNodeST = 0;
   struct con_clustST *cmpNodeST = 0;
   struct con_clustST *lastNodeST = 0; /*final compare*/

   struct res_edDist resEdDistStackST;

   float qryPercOverlapF = 0;
   float refPercOverlapF = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec02:
   ^   - initialize and set up variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(*indexSTPtr)
      freeHeap_index_clustST(*indexSTPtr);
   *indexSTPtr = 0;

   init_samEntry(&bestSamStackST);
   init_res_edDist(&resEdDistStackST);

   errSC = setup_samEntry(&bestSamStackST);

   if(errSC)
      goto memErr_fun05_sec07_sub02;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec03:
   ^   - set up index array (and score)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   numReadsSL = getNumLines_clustST(samFILE);

   if(! numReadsSL)
      goto fileErr_fun05_sec07_sub03;

   *indexSTPtr =
      mk_index_clustST(
         numReadsSL + 1, /*number lines in sam file*/
         clustSetSTPtr,  /*has settings for indexing*/
         samSTPtr,       /*for reading sam file*/
         buffStrPtr,     /*for reading sam file*/
         lenBuffULPtr,   /*length of buffStrPtr*/
         samFILE 
      );

   if(! indexSTPtr)
      goto memErr_fun05_sec07_sub02;

   if(
        (*indexSTPtr)->keptSL
      < (slong) clustSetSTPtr->minDepthUI
   ) goto lowReads_fun05_sec07_sub04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec04:
   ^   - cluster reads
   ^   o fun05 sec04 sub01:
   ^     - print progress + start loop
   ^   o fun05 sec04 sub02:
   ^     - get best read and its cluster
   ^   o fun05 sec04 sub03:
   ^     - build consensus and get depth
   ^   o fun05 sec04 sub04:
   ^     - compare consensus (and update read depth)
   ^   o fun05 sec04 sub05:
   ^     - print comparision result to log
   ^   o fun05 sec04 sub06:
   ^     - if merged, try to rebuild consensus
   ^   o fun05 sec04 sub07:
   ^     - zero extra pointers before moving on
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec04 Sub01:
   *   - print progress + start loop
   \*****************************************************/

   numReadsSL = (*indexSTPtr)->keptSL;
   lastReadsSL = numReadsSL;

   if(logFILE)
      fprintf(
         logFILE,
         "kept reads: %li\n\nClusters:\n",
         numReadsSL
      );

   if(clustSetSTPtr->repIntervalSL > 0)
   { /*If: user wanted logging*/
      fprintf(
         stderr,
         "kept reads: %li\n",
         numReadsSL
      ); /*If reporting read counts*/

      fflush(stderr);
   } /*If: user wanted logging*/

   clustSetSTPtr->clustSI = 1; /*on first cluster*/

   while(numReadsSL >= clustSetSTPtr->minDepthUI)
   { /*Loop: cluster reads*/

      if(clustSetSTPtr->repIntervalSL > 0)
      { /*If: reporting status*/
         if(
               (lastReadsSL - numReadsSL)
            >= clustSetSTPtr->repIntervalSL
         ){ /*If: logging interval*/
            fprintf(
               stderr,
               "Reads left: %li\n",
               numReadsSL
            ); /*update on how many reads left*/

            lastReadsSL = numReadsSL;
         } /*If: logging interval*/

         fflush(stderr);
      } /*If: reporting status*/

      /**************************************************\
      * Fun05 Sec04 Sub02:
      *   - get best read and its cluster
      \**************************************************/

      bestIndexSL =
         getBestRead_edClust(
            &bestSamStackST,/*get un-clusetred best read*/
            *indexSTPtr,    /*has cluster assignments*/
            *buffStrPtr,    /*for getting best read*/
            samFILE
         ); /*find best read to cluster with*/

      if(bestIndexSL < 0)
      { /*If: had problem*/
         if(bestIndexSL == def_noReads_edClust)
            break; /*finished*/

         else if(bestIndexSL == def_memErr_edClust)
         { /*Else If: memory error for best read*/
            fprintf(
               (FILE *) logFILE,
               "\nmemory error getting best read\n"
            );

            goto memErr_fun05_sec07_sub02;
         } /*Else If: memory error for best read*/

         else
            goto fileErr_fun05_sec07_sub03;
      } /*If: had problem*/

      if(conHeapAryST)
         freeHeapAry_conNt_tbCon(
            conHeapAryST,
            lenConUI
         );
      conHeapAryST = 0;
      lenConUI = 0;

      conHeapAryST =
         findNumMap_edClust(
            &depthSL,
            &lenConUI,
            &bestSamStackST, /*read to cluster around*/
            bestIndexSL,  /*read index to cluster around*/
            clustSetSTPtr, /*settings for clustering*/
            conSetSTPtr,   /*consensus settings*/
            *indexSTPtr,   /*has cluster array to update*/
            &resEdDistStackST,/*getting edit distance*/
            clustSetSTPtr->clustSI, /*cluster on*/
            0,             /*not a consensus*/
            samSTPtr,      /*for getting reads*/
            *buffStrPtr, /*for getting reads*/
            samFILE
         ); /*find cluster depth for best read*/

      if(depthSL == def_memErr_edClust)
      { /*If: memory error*/
         fprintf(
               (FILE *) logFILE,
               "\nerror mapping reads to best read\n"
            );

         goto memErr_fun05_sec07_sub02;
      } /*If: memory error*/

      if(depthSL == def_fileErr_edClust)
         goto fileErr_fun05_sec07_sub03;

      if(depthSL < (sint) clustSetSTPtr->minDepthUI)
      { /*If: read depth to low*/
         (*indexSTPtr)->clustArySI[bestIndexSL] =
            def_discard_clustST;/*mark as discarded read*/

         freeHeapAry_conNt_tbCon(
            conHeapAryST,
            lenConUI
         );

         conHeapAryST = 0;
         lenConUI = 0;

         --numReadsSL; /*1 less read*/
         continue;     /*find next best read*/
      } /*If: read depth to low*/

      /**************************************************\
      * Fun05 Sec04 Sub03:
      *   - build consensus and get depth
      \**************************************************/
      
      for(
         ucCon = 0;
         ucCon <= clustSetSTPtr->conRebuildUC;
         ++ucCon
      ){ /*Loop: build cosensus*/
         if(conSamHeapST)
            freeHeap_samEntry(conSamHeapST);
         conSamHeapST = 0;

         conSamHeapST =
            noFragCollapse_tbCon(
               conHeapAryST,  /*uncollapsed consensus*/
               lenConUI,      /*length of conHeapAryST*/
               bestSamStackST.refIdStr, /*reference id*/
               &numMaskUI,    /*do not care about*/
               conSetSTPtr,   /*consensus settings*/
               &errSC         /*holds error type*/
            ); /*build the consensus*/

         freeHeapAry_conNt_tbCon(
            conHeapAryST,
            lenConUI
         );

         conHeapAryST = 0;
         lenConUI = 0;

         if(errSC)
         { /*If: memory error*/
            if(errSC == def_noSeq_tbConDefs)
            { /*If: no consensus made*/
               errSC = def_noCon_clustST;
               break;
            } /*If: no consensus made*/

            fprintf(
                  (FILE *) logFILE,
                  "\nerror building consensus\n"
               );

            goto memErr_fun05_sec07_sub02;
         } /*If: memory error*/

         if(
               ! conSamHeapST
            ||
                   (float) numMaskUI
                 / (float) conSamHeapST->alnReadLenUI
               > clustSetSTPtr->maxNPercF
         ){ /*If: no cosensus made*/
            errSC = def_noCon_clustST;
            break;
         } /*If: no cosensus made*/

         conSamHeapST->refStartUI +=
            bestSamStackST.refStartUI;

         conSamHeapST->refEndUI =
            bestSamStackST.refStartUI;

         conSamHeapST->refEndUI +=
            conSamHeapST->alnReadLenUI;

         --conSamHeapST->refEndUI; /*convert to index 0*/

         conHeapAryST = 
            findNumMap_edClust(
               &depthSL,
               &lenConUI,
               conSamHeapST,  /*consensus to cluster*/
               bestIndexSL,   /*original read index*/
               clustSetSTPtr, /*settings for clustering*/
               conSetSTPtr,   /*consensus settings*/
               *indexSTPtr,   /*cluster array to update*/
               &resEdDistStackST,/*getting edit distance*/
               clustSetSTPtr->clustSI,  /*cluster on*/
               1,             /*is a consensus*/
               samSTPtr,      /*for getting reads*/
               *buffStrPtr,    /*for getting reads*/
               samFILE
            ); /*find cluster depth for best read*/

         if(depthSL == def_memErr_edClust)
         { /*If: memory error*/
            fprintf(
                  (FILE *) logFILE,
                  "\nerror mapping reads to consensus\n"
               );

            goto memErr_fun05_sec07_sub02;
         } /*If: memory error*/

         if(depthSL < (sint) clustSetSTPtr->minDepthUI)
         { /*If: read depth to low*/
            errSC = def_noCon_clustST;
            break;
         } /*If: read depth to low*/
      } /*Loop: build cosensus*/

      if(conHeapAryST)
         freeHeapAry_conNt_tbCon(
            conHeapAryST,
            lenConUI
         );

      conHeapAryST = 0;
      lenConUI = 0;

      if(
            errSC == def_noCon_clustST
         || depthSL < (sint) clustSetSTPtr->minDepthUI
      ){ /*If: could not build consensus or low depth*/
         (*indexSTPtr)->clustArySI[bestIndexSL] =
            def_discard_clustST;/*mark as discarded read*/

         if(conSamHeapST)
            freeHeap_samEntry(conSamHeapST);
         conSamHeapST = 0;

         --numReadsSL; /*1 less read*/
         continue;     /*find next best read*/
      } /*If: could not build consensus or low depth*/

      /**************************************************\
      * Fun05 Sec04 Sub04:
      *   - compare consensus (and update read depth)
      \**************************************************/

      conNodeST =
         mk_con_clustST(
            conSamHeapST,
            clustSetSTPtr->clustSI,
            (ulong) depthSL
         ); /*build node con_clustST node for consensus*/

      conNodeST->startUI = clustSetSTPtr->startUI;
      conNodeST->endUI = clustSetSTPtr->endUI;

      if(! conNodeST)
      { /*If: memory error*/
         fprintf(
               (FILE *) logFILE,
               "\nerror making consensus struct\n"
            );

         goto memErr_fun05_sec07_sub02;
      } /*If: memory error*/

      clustReadsUL += depthSL;
      numReadsSL -= depthSL;
         /*at this points reads are in a cluster*/

      if(logFILE)
         fprintf(
            logFILE,
            "   clust %i: %lu reads\t%s\n",
            clustSetSTPtr->clustSI,
            depthSL,
            bestSamStackST.qryIdStr
         );

      if(! retHeapST)
         retHeapST = conNodeST;

      else
      { /*Else: not first cosensus*/
         cmpNodeST = retHeapST; /*first cosensus*/

         depthSL =
            cmpCons_clustST(
               cmpNodeST,
               conNodeST,
               0,          /*want bounds check*/
               *indexSTPtr,
               clustSetSTPtr,
               &resEdDistStackST /*getting edit distance*/
            ); /*compare cons and merge similar clusters*/

         while(cmpNodeST->nextST)
         { /*Loop: compare all consensus*/
            if(depthSL)
            { /*If: merged clusters*/
               freeHeap_con_clustST(conNodeST);
               conNodeST = 0;
               break;
            } /*If: merged clusters*/

            cmpNodeST = cmpNodeST->nextST;

            depthSL =
               cmpCons_clustST(
                  cmpNodeST,
                  conNodeST,
                  0,          /*want bounds check*/
                  *indexSTPtr,
                  clustSetSTPtr,
                  &resEdDistStackST
               ); /*compare and merge similar clusters*/
         } /*Loop: compare all consensus*/

         if(depthSL)
         { /*If: merged clusters*/
            freeHeap_con_clustST(conNodeST);
            conNodeST = 0;
         } /*If: merged clusters*/

         else
            cmpNodeST->nextST = conNodeST; /*add to list*/
      } /*Else: not first cosensus*/

      /**************************************************\
      * Fun05 Sec04 Sub05:
      *   - print comparison result to log
      \**************************************************/

      if(conNodeST)
      { /*If: kept cluster*/
         conNodeST = 0;
         ++clustSetSTPtr->clustSI;/*move to next cluster*/
      } /*If: kept cluster*/

      else
      { /*Else: merged clusters*/
         if(logFILE)
            fprintf(
              logFILE,
              "   clust %i->%i: clust %i has %lu reads\n",
              clustSetSTPtr->clustSI,
              cmpNodeST->clustSI,
              cmpNodeST->clustSI,
              cmpNodeST->numReadsUL
            );

         /***********************************************\
         * Fun05 Sec04 Sub06:
         *   - if merged, try to rebuild consensus
         \***********************************************/

         clustSetSTPtr->startUI = cmpNodeST->startUI;
         clustSetSTPtr->endUI = cmpNodeST->endUI;

         conSamHeapST =
            getCon_clustST(
               *indexSTPtr,
               cmpNodeST->clustSI,
               clustSetSTPtr,
               conSetSTPtr,
               samSTPtr,
               *buffStrPtr,
               samFILE,
               &errSC
            ); /*rebuild the consensus*/

         if(errSC == def_memErr_clustST)
         { /*If: memory error*/
            fprintf(
                  (FILE *) logFILE,
                  "\nerror rebuilding con after clust\n"
               );

            goto memErr_fun05_sec07_sub02;
         } /*If: memory error*/

         if(errSC != def_noCon_clustST)
         { /*If: rebuilt consensus*/
            freeHeap_samEntry(cmpNodeST->samSTPtr);
            cmpNodeST->samSTPtr = conSamHeapST;
         } /*If: rebuilt consensus*/
      } /*Else: merged clusters*/

      /**************************************************\
      * Fun05 Sec04 Sub07:
      *   - zero extra pointers before moving on
      \**************************************************/

      if(logFILE)
         fflush(logFILE);

      cmpNodeST = 0;
      conSamHeapST = 0;
   } /*Loop: cluster reads*/

   --clustSetSTPtr->clustSI; /*no longer on last cluster*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec05:
   ^   - final consensus comparison
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) logFILE,
      "\nOverlap Merging:\n"
   );

   if(! retHeapST)
   { /*If: no consensus detected*/
      fprintf(
         (FILE *) logFILE,
         "\nNo clusters found\n"
      );

      goto lowReads_fun05_sec07_sub04;
   } /*If: no consensus detected*/

   cmpNodeST = retHeapST; /*first cosensus*/
   lastNodeST = cmpNodeST;
   conNodeST = cmpNodeST->nextST;

   while(cmpNodeST)
   { /*Loop: compare all consensus*/
      if(! conNodeST)
         break;

      depthSL = 0; 

      qryPercOverlapF =
         (float)
         min_genMath(
            cmpNodeST->endUI,
            conNodeST->endUI
         );

      qryPercOverlapF -=
         (float)
         max_genMath(
            cmpNodeST->startUI,
            conNodeST->startUI
         ); /*find overlap length*/

      if(qryPercOverlapF <= 0)
         goto nextCon_fun05_sec05;
      
      refPercOverlapF =
           qryPercOverlapF
         / (float) cmpNodeST->samSTPtr->alnReadLenUI;
         /*find reference % coverate*/

      qryPercOverlapF =
           qryPercOverlapF
         / (float) conNodeST->samSTPtr->alnReadLenUI;
        /*find query % coverage*/

      if(  
            qryPercOverlapF < clustSetSTPtr->percOverlapF
         && refPercOverlapF < clustSetSTPtr->percOverlapF
      ) goto nextCon_fun05_sec05;
        /*not enough overlap*/

      lenConUI = (uint) conNodeST->clustSI;
      numMaskUI = (uint) cmpNodeST->clustSI;

      depthSL =
         cmpCons_clustST(
            cmpNodeST,
            conNodeST,
            1,          /*checking bounds before call*/
            *indexSTPtr,
            clustSetSTPtr,
            &resEdDistStackST
         ); /*compare and merge similar clusters*/

         if(depthSL)
         { /*If: merged clusters*/

            if(lenConUI == (uint) cmpNodeST->clustSI)
               lenConUI = numMaskUI;

            if(logFILE)
             fprintf(
              logFILE,
              "   clust %u->%i: clust %i has %lu reads\n",
              lenConUI,
              cmpNodeST->clustSI,
              cmpNodeST->clustSI,
              cmpNodeST->numReadsUL
             );

            lastNodeST->nextST =
               freeHeap_con_clustST(conNodeST);

            conNodeST = lastNodeST;
         } /*If: merged clusters*/

      nextCon_fun05_sec05:;
         if(
               ! conNodeST
            || ! conNodeST->nextST
         ){ /*If: finished one set of comparisons*/
            if(! cmpNodeST->nextST)
               break; /*finished*/

            cmpNodeST = cmpNodeST->nextST;
            conNodeST = cmpNodeST->nextST;
            lastNodeST = cmpNodeST;
         } /*If: finished one set of comparisons*/

         else
         { /*Else: consensuses were different*/
            if(conNodeST == lastNodeST)
               conNodeST = lastNodeST->nextST;
               /*merged clusters*/
 
            else
            { /*Else: did not merge clusters*/
               conNodeST = conNodeST->nextST;
               lastNodeST = lastNodeST->nextST;
            } /*Else: did not merge clusters*/
         } /*Else: consensuses were different*/
   } /*Loop: compare all consensus*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec06:
   ^   - remove low read depth clusters
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) logFILE,
      "\nDepths:\n"
   );

   if(clustSetSTPtr->minPercDepthF > 0)
   { /*If: have minimum read depth to filter with*/

      getMaxDepth_edClust(
         retHeapST,
         *indexSTPtr,
         clustSetSTPtr->percOverlapF
      ); /*get maximum possible depths for each cluster*/

      cmpNodeST = retHeapST; /*first cosensus*/
      lastNodeST = 0;

      while(cmpNodeST)
      { /*Loop: remove low read depth clusters*/
         qryPercOverlapF =
            (float)
              (
                  (double) cmpNodeST->numReadsUL
                / (double) cmpNodeST->maxReadsUL
              );

         if(qryPercOverlapF <clustSetSTPtr->minPercDepthF)
         { /*If: discarding cluster*/

            fprintf(
               (FILE *) logFILE,
               "   clust %i; low read depth (%%%0.4f)\n",
               cmpNodeST->clustSI,
               qryPercOverlapF
            );
               
            if(cmpNodeST == retHeapST)
            { /*If: removing head of list*/
               retHeapST =
                  freeHeap_con_clustST(cmpNodeST);

               cmpNodeST = retHeapST;
            } /*If: removing head of list*/

            else
            { /*Else: freeing non-head node*/
               lastNodeST->nextST =
                  freeHeap_con_clustST(cmpNodeST);

               cmpNodeST = lastNodeST->nextST;
            } /*Else: freeing non-head node*/
         } /*If: discarding cluster*/

         else
         { /*Else: keeping cluster*/
          fprintf(
           (FILE *) logFILE,
           "   clust %i has %lu of %lu (%0.4f%%) reads\n",
           cmpNodeST->clustSI,
           cmpNodeST->numReadsUL,
           cmpNodeST->maxReadsUL,
           qryPercOverlapF * 100
          );

            lastNodeST = cmpNodeST;
            cmpNodeST = lastNodeST->nextST;
         } /*Else: keeping cluster*/
      } /*Loop: remove low read depth clusters*/
   } /*If: have minimum read depth to filter with*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec07:
   ^   - clean up
   ^   o fun05 sec07 sub01:
   ^     - no error clean up
   ^   o fun05 sec07 sub02:
   ^     - memory error clean up
   ^   o fun05 sec07 sub03:
   ^     - file error clean up
   ^   o fun05 sec07 sub04:
   ^     - to few reads error clean up
   ^   o fun05 sec07 sub05:
   ^     - error clean up (all errors call)
   ^   o fun05 sec07 sub06:
   ^     - general clean up (everything calls)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec07 Sub01:
   *   - no error clean up
   \*****************************************************/

   if(logFILE)
   { /*If: logging*/
      fprintf(
         logFILE,
         "\nkept reads: %li\n",
         (*indexSTPtr)->keptSL
         
      );

      fprintf(
         logFILE,
         "clustered reads: %lu\n",
         clustReadsUL
      );

      fprintf(
         logFILE,
         "percent clustered:%0.2f\n",
         100
           * (
                (float) clustReadsUL
              / (float) (*indexSTPtr)->keptSL
            )
      );

      fflush(logFILE);
   } /*If: logging*/

   fprintf(
      stderr,
      "\nkept reads: %li\n",
      (*indexSTPtr)->keptSL
      
   );

   fprintf(
      stderr,
      "clustered reads: %lu\n",
      clustReadsUL
   );

   fprintf(
      stderr,
      "percent clustered:%0.2f\n",
      100
        * (
             (float) clustReadsUL
           / (float) (*indexSTPtr)->keptSL
         )
   );

   fflush(stderr);

   *errSCPtr = 0;
   goto ret_fun05_sec07_sub06;

   /*****************************************************\
   * Fun05 Sec07 Sub02:
   *   - memory error clean up
   \*****************************************************/

   memErr_fun05_sec07_sub02:;
      if(logFILE)
         fprintf(
            logFILE,
            "memory error when (%li) reads left\n",
            numReadsSL
         );

      *errSCPtr = def_memErr_edClust;
      goto errCleanup_fun05_sec07_sub05;

   /*****************************************************\
   * Fun05 Sec07 Sub03:
   *   - file error clean up
   \*****************************************************/

   fileErr_fun05_sec07_sub03:;
      *errSCPtr = def_fileErr_edClust;

      if(logFILE)
         fprintf(
            logFILE,
            "file error when (%li) reads left\n",
            numReadsSL
         );

      goto errCleanup_fun05_sec07_sub05;

   /*****************************************************\
   * Fun05 Sec07 Sub04:
   *   - to few reads error clean up
   \*****************************************************/

   lowReads_fun05_sec07_sub04:;
      *errSCPtr = def_noReads_edClust;

      if(logFILE)
         fprintf(
           logFILE,
           "to few reads (%lu) after filter to cluster\n",
           (*indexSTPtr)->keptSL
         );

      goto errCleanup_fun05_sec07_sub05;

   /*****************************************************\
   * Fun05 Sec07 Sub05:
   *   - error clean up (all errors call)
   \*****************************************************/

   errCleanup_fun05_sec07_sub05:;
      if(retHeapST)
         freeHeapList_con_clustST(retHeapST);
      retHeapST = 0;
      goto ret_fun05_sec07_sub06;

   /*****************************************************\
   * Fun05 Sec07 Sub06:
   *   - general clean up (everything calls)
   \*****************************************************/

   ret_fun05_sec07_sub06:;
      fseek(
         samFILE,
         0,
         SEEK_SET
      );

      if(conHeapAryST)
         freeHeapAry_conNt_tbCon(
            conHeapAryST,
            lenConUI
         );
      conHeapAryST = 0;
      lenConUI = 0;

      freeStack_samEntry(&bestSamStackST);
      freeStack_res_edDist(&resEdDistStackST);

      return retHeapST;
} /*cluster_edClust*/

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
