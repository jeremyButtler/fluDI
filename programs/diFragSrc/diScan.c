/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   - diScan:
'     o had functions to scan for DI fragments
'   o header:
'     - included libraries
'   o fun01: findSeg_diScan
'     - uses kmer profiling to find the best matching
'       segment
'   o fun02: waterScan_diScan
'     - scan for DI fragments using kmer profiling and
'       Watermen alignment
'   o fun03: phead_diScan
'     - print out header for fragment tsv
'   o fun04: pfrag_diScan
'     - print out the fragment
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

#include "diScan.h"

#include <stdio.h>

#include "../genLib/genMath.h"
#include "../genLib/kmerCnt.h"
#include "../genLib/seqST.h"
#include "../genLib/samEntry.h"

#include "../genAln/alnSet.h"
#include "../genAln/dirMatrix.h"
#include "../genAln/water.h"

#include "../diCoordsSrc/diCoords.h"

/*.h only files*/
#include "../genLib/dataTypeShortHand.h"
#include "../genAln/alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #incldue "../genLib/ulCp.h"
!   - .c  #include "../genAln/indexToCoord.h"
!   - .h  #include "../genLib/ntTo2Bit.h"
!   - .h  #include "../genLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: findSeg_diScan
|   - uses kmer profiling to find the best matching
|     segment
| Input:
|   - seqSTPtr:
|     o pointer to seqST with sequence to find segment
|       for
|   - refAryST:
|     o pointer to kmerCnt structure array with segments
|       to compare againts
|   - numSegUI:
|     o number segments (references) in refAryST
|   - kmerArySI:
|     o signed int array to hold the kmers in seqSTPtr
|   - cntArySI:
|     o signed int array to hold kmer counts for seqSTPtr
|   - lenKmerUC:
|     o length of one kmer
|   - minKmerPercF:
|     o min percent of shared kmers to keep a read
|   - maxKmersSI:
|     o pointer to signed long to hold number of kmers
|       that mapped (- for reverse, + for forward)
| Output:
|   - Modifies:
|     o maxKmersSI to have the number of matches
|       - is positive if foward reference was best
|       - is negative if reverse reference was best
|   - Returns:
|     o -1 if no segment mapped
|     o > 0 (segment index) if a segment mapped
\-------------------------------------------------------*/
signed int
findSeg_diScan(
   struct seqST *seqSTPtr,    /*sequence to check*/
   struct kmerCnt *refAryST,  /*array of references*/
   unsigned int numSegUI,      /*length of ref array*/
   signed int *kmerArySI,      /*holds sequence kmers*/
   signed int *cntArySI,       /*holds kmer counts*/
   unsigned char lenKmerUC,    /*length of one kmer*/
   float minKmerPercF,         /*min perc kmers to keep*/
   signed int *maxKmersSI      /*will hold kmer count*/
){
   sint siSeg = 0;
   sint totalKmersSI = 0;
   sint numKmersSI = 0;
   sint maxSegSI = 0;

   float scoreF = 0;

   *maxKmersSI = 0;

   totalKmersSI =
      ntToKmerAry_kmerCnt(
         seqSTPtr,
         lenKmerUC,
         kmerArySI,
         cntArySI
      ); /*set up kmer count array*/

   for(
      siSeg = 0;
      siSeg < (sint) numSegUI;
      ++siSeg
   ){ /*Loop: find the best segment*/
      numKmersSI =
         get_kmerCnt(
            &refAryST[siSeg],
            kmerArySI,
            cntArySI
         ); /*get kmer counts for segments*/

      if(
           ab_genMath(numKmersSI)
         > ab_genMath(*maxKmersSI)
      ){ /*If: found a better segment*/
         *maxKmersSI = numKmersSI;
         maxSegSI = siSeg;
      } /*If: found a better segment*/
   } /*Loop: find the best segment*/

   /*check if enough kmers were present to align*/
   scoreF = (float) ab_genMath(*maxKmersSI);
   scoreF /= (float) totalKmersSI;

   return maxSegSI | (sint) ( -(scoreF < minKmerPercF) );
      /*-1 if score is to low, else number*/
} /*findSeg_diScan*/

/*-------------------------------------------------------\
| Fun02: waterScan_diScan
|   - scan for DI fragments using kmer profiling and
|     Watermen alignment
| Input:
|    - seqSTPtr:
|      o pointer to seqST with sequence to scan
|   - refAryST:
|     o pointer to kmerCnt structure array with segments
|       to compare againts
|   - lenRefUI:
|     o number of kmerCnt structures to scan in refAryST
|   - kmerArySI:
|     o signed int array to hold the kmers in seqSTPtr
|   - cntArySI:
|     o signed int array to hold kmer counts for seqSTPtr
|   - minPercScoreF:
|     o minimum percent score from waterman alingment to
|       check for DI (count as mapped)
|   - minKmerPercF:
|     o min percent of shared kmers to keep a read
|   - lenKmerUC:
|     o length of one kmer
|   - samSTPtr:
|     o pointer to samEntry structure to hold aligned
|       sequence
|   - segSIPtr:
|     o pointer to signed int ot hold the segment mapped
|       to or a negative number
|   - minDIDelUI:
|     o minimum deletion size to flag as a DI sequence
|   - minEndNtUI:
|     o how many bases in a DI event must be to be a DI
|       event
|   - numKmersSIPtr:
|     o pointer to signed int to hold the number of kmers
|       shared with the reference
|   - alnSetSTPtr:
|     o pointer to alnSet struct with alignment settings
|   - matrixSTPtr:
|     o pointer to dirMatrix struct to use in alignment
| Output:
|   - Modifies:
|     o samSTPtr to hold the new alignment
|     o matrixSTPtr to have the direction matrix, score,
|       and index from alignment
|     o kmerArySI to have the sequences kmers list
|     o cntArySI to have the sequences kmer counts
|     o segSIPtr to have the segment mapped to or -1
|   - Returns:
|     o 0 if there were no DI events
|     o > 0 if found DI events (number of events returned)
|     o def_noMatch_diScan (-1) if not a DI sequence
|     o def_memErr_diScan (-2) if memory error
\-------------------------------------------------------*/
signed int
waterScan_diScan(
   struct seqST *seqSTPtr,/*sequence to check*/
   struct kmerCnt *refAryST,  /*array of references*/
   unsigned int lenRefUI,     /*# of references*/
   signed int *kmerArySI,     /*holds sequence kmers*/
   signed int *cntArySI,      /*holds kmer counts*/
   float minPercScoreF,       /*min % score to check DIs*/
   float minKmerPercF,         /*min perc kmers to keep*/
   unsigned char lenKmerUC,    /*length of one kmer*/
   unsigned int minDIDelUI,   /*min del size in DI*/
   unsigned int minPadNtUI,   /*min start/end length*/
   struct samEntry *samSTPtr, /*holds alignment*/
   signed int *segSIPtr,      /*segment mapped to*/
   signed int *numKmersSIPtr, /*number kmers shared*/
   struct alnSet *alnSetSTPtr,/*alignment settings*/
   struct dirMatrix *matrixSTPtr /*matrix for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - scan for DI fragments using kmer profiling and
   '     Watermen alignment
   '   o fun02 sec01:
   '     - variable declerations
   '   o fun02 sec02:
   '     - find the best segment
   '   o fun02 sec03:
   '     - get alignment
   '   o fun02 sec04:
   '     - filter out low scores and find DI count
   '   o fun02 sec05:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;       /*for detecting errors*/

   sint numDISI = 0;     /*number of DI events*/

   float scoreF = 0;     /*score from alignment*/
   float maxScoreF = 0;  /*maximum score possible*/

   uint numAnonUI = 0;   /*# anonymous bases (ignore)*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - find the best segment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *numKmersSIPtr = 0;

   *segSIPtr = 
      findSeg_diScan(
         seqSTPtr,
         refAryST,          /*array of references*/
         lenRefUI,          /*number of references*/
         kmerArySI,         /*holds sequence kmers*/
         cntArySI,          /*holds kmer counts*/
         lenKmerUC,         /*length of one kmer*/
         minKmerPercF,      /*min % shared kmers to keep*/
         numKmersSIPtr      /*-1 or number matched kmers*/
      );

   if(*segSIPtr < 0)
      goto noMatch_fun02_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - get alignment
   ^   o fun02 sec03 sub01:
   ^     - align to foward reference
   ^   o fun02 sec03 sub02:
   ^     - align to reverse complement reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec03 Sub01:
   *   - align to foward reference
   \*****************************************************/

   seqSTPtr->offsetUL = 0;
   seqSTPtr->endAlnUL = seqSTPtr->lenSeqUL - 1;

   seqToIndex_alnSet(seqSTPtr->seqStr);

   if(*numKmersSIPtr > 0)
   { /*If: best match forward segment*/
      seqToIndex_alnSet(
         refAryST[*segSIPtr].forSeqST->seqStr
      );

      refAryST[*segSIPtr].forSeqST->offsetUL = 0;

      refAryST[*segSIPtr].forSeqST->endAlnUL =
         refAryST[*segSIPtr].forSeqST->lenSeqUL - 1;

      scoreF =
         (float)
         water(
            seqSTPtr,
            refAryST[*segSIPtr].forSeqST,
            matrixSTPtr,
            alnSetSTPtr
         ); /*align the sequence*/

      if(scoreF == 0)
      { /*If: low score*/
         indexToSeq_alnSet(
            refAryST[*segSIPtr].forSeqST->seqStr
         );

         goto memErr_fun02_sec06; /*memory error*/
      } /*If: low score*/

      errSC =
         getAln_dirMatrix(
            matrixSTPtr,
            0,                 /*use index from matrix*/
            0,                 /*forward alignment*/
            seqSTPtr,
            refAryST[*segSIPtr].forSeqST,
            samSTPtr,
            &numAnonUI,        /*discarding this*/
            alnSetSTPtr
         ); /*get the alignment*/


      indexToSeq_alnSet(
         refAryST[*segSIPtr].forSeqST->seqStr
      );

      if(errSC)
         goto memErr_fun02_sec06; /*memory error*/
   } /*If: best match forward segment*/

   /*****************************************************\
   * Fun02 Sec03 Sub02:
   *   - align to reverse complement reference
   \*****************************************************/

   else
   { /*Else: best match reverse segment*/
      seqToIndex_alnSet(
         refAryST[*segSIPtr].revSeqST->seqStr
      );

      refAryST[*segSIPtr].revSeqST->offsetUL = 0;

      refAryST[*segSIPtr].revSeqST->endAlnUL = 
         refAryST[*segSIPtr].revSeqST->lenSeqUL - 1;

      scoreF =
         (float)
         water(
            seqSTPtr,
            refAryST[*segSIPtr].revSeqST,
            matrixSTPtr,
            alnSetSTPtr
         ); /*align the sequences*/

      if(scoreF == 0)
      { /*If: low score*/
         indexToSeq_alnSet(
            refAryST[*segSIPtr].revSeqST->seqStr
         );

         goto memErr_fun02_sec06; /*memory error*/
      } /*If: low score*/

      errSC =
         getAln_dirMatrix(
            matrixSTPtr,
            0,                 /*use index from matrix*/
            1,                 /*revese alignment*/
            seqSTPtr,
            refAryST[*segSIPtr].revSeqST,
            samSTPtr,
            &numAnonUI,        /*discarding this*/
            alnSetSTPtr
         );

      indexToSeq_alnSet(
         refAryST[*segSIPtr].revSeqST->seqStr
      );

      if(errSC)
         goto memErr_fun02_sec06; /*memory error*/
   } /*Else: best match reverse segment*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^   - filter out low scores and find DI count
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   indexToSeq_alnSet(seqSTPtr->seqStr);

   /*account for using integer values for scores*/
   scoreF = (float) matrixSTPtr->scoreSL;
   scoreF /= def_scoreAdj_alnDefs;

   maxScoreF = seqSTPtr->lenSeqUL;
   maxScoreF *= def_matchScore_diScan;

   if(scoreF / maxScoreF < minPercScoreF)
      goto noMatch_fun02_sec06;

   numDISI =
      scan_diCoords(
         samSTPtr,
         minDIDelUI,
         minPadNtUI
      ); /*find the number of DI events*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec05:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return numDISI;

   noMatch_fun02_sec06:;
   return (sint) def_noMatch_diScan;

   memErr_fun02_sec06:;
   return (sint) def_memErr_diScan;
} /*waterScan_diScan*/

/*-------------------------------------------------------\
| Fun03: phead_diScan
|   - print out header for fragment tsv
| Input:
|   - outFILE:
|     o file to print header to
|   - lenKmerUC:
|     o length of one kmer
| Output:
|   - Prints:
|     o header to outFILE
\-------------------------------------------------------*/
void
phead_diScan(
   void * outFILE,
   unsigned char lenKmerUC
){
   fprintf(
      (FILE *) outFILE,
      "id\tref\tclass\tdi_events\tdir\tread_len\taln_len"
   );

   fprintf(
      (FILE *) outFILE,
      "\tref_start\tref_end\tref_len\tscore\tmax_score"
   );

   fprintf(
      (FILE *) outFILE,
      "\tnum_match\tnum_snp\tnum_ins\tnum_del\tnum_mask"
   );

   fprintf(
      (FILE *) outFILE,
      "\tshared_%umers\tmed_q\tmean_q\n",
      lenKmerUC
   );
} /*phead_diScan*/

/*-------------------------------------------------------\
| Fun04: pfrag_diScan
|   - print out the fragment
| Input:
|   - samSTPtr:
|     o pointer to samEntry structure with sequence to
|       print
|   - numDISI:
|     o number of DI events detected
|   - segLenSI:
|     o the length of the mapped segment
|   - scoreSL:
|     o score for waterman alignment
|   - numKmersSI:
|     o number of kmers shared between ref and read
|   - outFILE:
|     o file to print read to
| Output:
|   - Prints:
|     o read id, mapped segment, classifaction, and other
|       information to a tsv file
\-------------------------------------------------------*/
void
pfrag_diScan(
   struct samEntry *samSTPtr,
   signed int numDISI,
   signed int segLenSI,
   signed long scoreSL,
   signed int numKmersSI,
   void * outFILE
){
   schar *classStr = 0;
   schar *dirStr = 0;

   numKmersSI = ab_genMath(numKmersSI); /*make sure +*/

   if(numDISI > 0)
      classStr = (schar *) "diRNA";
   else
      classStr = (schar *) "vRNA";

   if(samSTPtr->flagUS & 16)
      dirStr = (schar *) "R";
   else
      dirStr = (schar *) "F";

   fprintf(
      (FILE *) outFILE,
      "%s\t%s\t%s\t%i\t%s\t%u\t%u\t%u\t%u\t%i\t%0.2f",
      samSTPtr->qryIdStr, /*read id*/
      samSTPtr->refIdStr, /*segment/ref*/
      classStr,           /*classification*/
      numDISI,            /*number DI events*/
      dirStr,             /*direction*/
      samSTPtr->readLenUI,
      samSTPtr->alnReadLenUI,
      samSTPtr->refStartUI + 1,
      samSTPtr->refEndUI + 1,
      segLenSI,
      (float) scoreSL / (float) def_scoreAdj_alnDefs
   ); /*print out general read stats*/

   fprintf(
      (FILE *) outFILE,
      "\t%0.2f\t%u\t%u\t%u\t%u\t%u\t%i\t%0.2f\t%0.2f\n",
      (float)
           maxScore_alnDefs(samSTPtr->readLenUI)
         / def_scoreAdj_alnDefs,
      samSTPtr->numMatchUI,
      samSTPtr->numSnpUI,
      samSTPtr->numInsUI,
      samSTPtr->numDelUI,
      samSTPtr->numMaskUI,
      numKmersSI,
      samSTPtr->medianQF,
      samSTPtr->meanQF
   ); /*print out read stats*/
} /*pfrag_diScan*/
