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
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "edDist.h"

#include <stdio.h>

#include "../genLib/ulCp.h"
#include "seqST.h"
#include "samEntry.h"

/*.h files only*/
#include "../genLib/dataTypeShortHand.h"
#include "../genLib/genMath.h" /*min .h macro only*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .c   #include "../genLib/base10str.h"
!   o .c   #include "../genLib/numToStr.h"
!   o .c   #include "../genLib/strAry.h"
!   o .h   #include "ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
){
   resSTPtr->probUI = 0;

   resSTPtr->overlapUI = 0;
   resSTPtr->startUI = 0;

   resSTPtr->numIndelUI = 0;
   resSTPtr->indelEventsUI = 0;

   resSTPtr->maxWinDistUI = 0;
   resSTPtr->probMaxWinUI = 0;
   resSTPtr->minWinDistUI = 0;
   resSTPtr->avgWinDistF = 0;

   resSTPtr->numFullWinUI = 0;
   resSTPtr->extraWinNtUI = 0;
   resSTPtr->lastWinDiffUI = 0;

   resSTPtr->edDistSL = 0;
} /*blank_res_edDist*/

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
){
   resSTPtr->depthAryUI = 0;
   resSTPtr->sizeDepthUL = 0;

   blank_res_edDist(resSTPtr);
} /*init_res_edDist*/

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
){
   if(! resSTPtr)
      return;

   if(resSTPtr->depthAryUI)
      free(resSTPtr->depthAryUI);

   init_res_edDist(resSTPtr);
} /*freeStack_res_edDist*/

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
){
   if(! resSTPtr)
      return;

   freeStack_res_edDist(resSTPtr);
   free(resSTPtr);
   resSTPtr = 0;
} /*freeHeap_res_edDist*/ 

/*-------------------------------------------------------\
| Fun05: isTransNt_edDist
|   - checks to see if query base is a transversion
| Input:
|   - qryUC:
|     o query nucleotide to check
|   - refUC:
|     o reference nucleotide to check
| Output:
|   - Returns:
|     o -1 if was a transversion
|     o 0 if was a transition
\-------------------------------------------------------*/
signed char
isTransNt_edDist(
   unsigned char qryUC,
   unsigned char refUC
){
    switch(qryUC & (~32))
    { /*switch, reverse complement*/
        case 'A':
        case 'G':
        case 'R': /*A/G*/
        /*Case: query is A or G*/
           switch(refUC & (~32))
           { /*Switch: query is A, G or R (A/G)*/
              case 'C':
                 return -1;
              case 'T':
                 return -1;
              case 'Y': /*C/T*/
                 return -1;
              default:
                 return 0;
           } /*Switch: query is A, G or Y (C/T)*/

           return 0;
        /*Case: query is G*/

        case 'T':
        case 'U':
        case 'C':
        case 'Y': /*T/C*/
        /*Case: query is C, T/U, or Y (T/C)*/ 
           switch(refUC & (~32))
           { /*Switch: query is C*/
              case 'A':
                 return -1;
              case 'G':
                 return -1;
              case 'R': /*A/G*/
                 return -1;
              default:
                 return 0;
           } /*Switch: query is C*/

           return 0;
        /*Case: query is C, T/U, or R (A/G)*/

        default:
           return 0;
    } /*switch, reverse complement*/

    return 0;
} /*isTransNt_edDist*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun06 TOC:
   '   - gets edit distances between the query & reference
   '   o fun06 sec01:
   '     - variable declerations
   '   o fun06 sec02:
   '     - check if can get edit distance
   '   o fun06 sec03:
   '     - initialize and check if have q-score entreis
   '   o fun06 sec04:
   '     - find start of overlap between query & reference
   '   o fun06 sec05:
   '     - find edit distance
   '   o fun06 sec08:
   '     - add last window to data
   '   o fun06 sec07:
   '     - return the edit distance
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint siQryCig = 0;  /*query cigar position*/
   sint siRefCig = 0;  /*reference cigar position*/
   sint refValSI = 0;
   sint qryValSI = 0;

   sint tmpSI = 0;

   uint uiQry = 0;     /*query nucleotide on*/
   uint uiRef = 0;     /*reference nucleotide on*/

   float overlapF = 0;
   slong retSL = 0;    /*edit distance*/

   schar qryQBl = 0;
   schar refQBl = 0;
   uchar qUC = 0;

   schar indelBl = 0;     /*for adding indesl to windows*/
   uint baseInWinUI = 0;  /*number of bases in window*/
   uint winDistUI = 0;    /*edit distance of window*/
   uint indelAddUI = 0;   /*find offset for indels*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec02:
   ^   - check if can get edit distance
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_res_edDist(resSTPtr);

   if(qrySTPtr->flagUS & 4)
      goto noMap_fun06_sec07;

   if(refSTPtr->flagUS & 4)
      goto noMap_fun06_sec07;

   if(qrySTPtr->refEndUI < refSTPtr->refStartUI)
      goto noOverlap_fun06_sec07;

   if(qrySTPtr->refStartUI > refSTPtr->refEndUI)
      goto noOverlap_fun06_sec07;

   if(minOverlapF > 0)
   { /*If: user wanted overlap check*/
      overlapF =
        (float)
        min_genMath(
           qrySTPtr->refEndUI,
           refSTPtr->refEndUI
        );

      overlapF -=
        (float)
        max_genMath(
           qrySTPtr->refStartUI,
           refSTPtr->refStartUI
        );

      if(
                overlapF / (float) refSTPtr->alnReadLenUI
              < minOverlapF
         &&
                overlapF / (float) qrySTPtr->alnReadLenUI
              < minOverlapF
      ) goto noOverlap_fun06_sec07;
   } /*If: user wanted overlap check*/

   if(
         ! qrySTPtr->seqStr
      || *qrySTPtr->seqStr == '*' 
      || *qrySTPtr->seqStr == '\0'
   ) goto noSeq_fun06_sec07;

   if(
         ! refSTPtr->seqStr
      || *refSTPtr->seqStr == '*' 
      || *refSTPtr->seqStr == '\0'
   ) goto noSeq_fun06_sec07;

   if(
      eql_ulCp(
         qrySTPtr->refIdStr,
         refSTPtr->refIdStr,
         0,
         '\0'
      )
   ) goto diffRef_fun06_sec07;
     /*reads are mapped to different references*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec03:
   ^   - initialize and check if have q-score entreis
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/


   if(
         qrySTPtr->qStr[0] == '\0'
      ||
         (
             qrySTPtr->qStr[0] == '*'
          && qrySTPtr->qStr[1] == '\0'
         )
   ) qryQBl = 0;

   else
      qryQBl = 1;

   if(
         refSTPtr->qStr[0] == '\0'
      ||
         (
             refSTPtr->qStr[0] == '*'
          && refSTPtr->qStr[1] == '\0'
         )
   ) refQBl = 0;

   else
      refQBl = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec04:
   ^   - find start of overlap between query & reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/


   refValSI = refSTPtr->cigArySI[0];
   qryValSI = qrySTPtr->cigArySI[0];

   if(qrySTPtr->refStartUI > refSTPtr->refStartUI)
   { /*If: move the reference foward*/
      uiQry = refSTPtr->refStartUI;

      findRefPos_samEntry(
         refSTPtr,
         &siRefCig,
         &refValSI,
         qrySTPtr->refStartUI,
         (sint *) &uiQry, /*reference pos (discard)*/
         (sint *) &uiRef
      ); /*set reference to first query base*/

      uiQry = 0;
   } /*If: move reference foward*/

   else if(qrySTPtr->refStartUI < refSTPtr->refStartUI)
   { /*Else If: move query foward*/
      uiRef = qrySTPtr->refStartUI;

      findRefPos_samEntry(
         qrySTPtr,
         &siQryCig,
         &qryValSI,
         refSTPtr->refStartUI, /*end query position*/
         (sint *) &uiRef, /*query pos (discard)*/
         (sint *) &uiQry  /*end query position*/
      ); /*set reference to first query base*/

      uiRef = 0;
   } /*Else If: move query foward*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec05:
   ^   - find edit distance
   ^   o fun06 sec05 sub01:
   ^     - get past masking at start
   ^   o fun06 sec05 sub02:
   ^     - start loop and check if at end (hard/soft mask)
   ^   o fun06 sec05 note01:
   ^     - table of cigar entries added together
   ^   o fun06 sec05 sub03:
   ^     - start switch/check match/possible snp combos
   ^   o fun06 sec05 sub04:
   ^     - check matches or same error type
   ^   o fun06 sec05 sub05:
   ^     - handel deletion cases
   ^   o fun06 sec05 sub06:
   ^     - handel insertion cases
   ^   o fun06 sec05 sub07:
   ^     - check if move to next cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun06 Sec05 Sub01:
   *   - get past masking at start
   \*****************************************************/

   if(refSTPtr->cigTypeStr[siRefCig] == 'S')
   { /*If: have softmasked region*/
      uiRef += refValSI;
      ++siRefCig;
      refValSI = refSTPtr->cigArySI[siRefCig];
   } /*If: have softmasked region*/

   if(qrySTPtr->cigTypeStr[siQryCig] == 'S')
   { /*If: have softmasked region*/
      uiQry += qryValSI;
      ++siQryCig;
      qryValSI = qrySTPtr->cigArySI[siQryCig];
   } /*If: have softmasked region*/

   if(refSTPtr->cigTypeStr[siRefCig] == 'H')
   { /*If: have hardmasked region*/
      ++siRefCig;
      refValSI = refSTPtr->cigArySI[siRefCig];
   } /*If: have hardmasked region*/

   if(qrySTPtr->cigTypeStr[siQryCig] == 'H')
   { /*If: have hardmaksed region*/
      ++siQryCig;
      qryValSI = qrySTPtr->cigArySI[siQryCig];
   } /*If: have hardmasked region*/

   /*****************************************************\
   * Fun06 Sec05 Sub02:
   *   - start loop and check if at end (hard/soft mask)
   \*****************************************************/

   while(
         siQryCig < (sint) qrySTPtr->lenCigUI
      && siRefCig < (sint) refSTPtr->lenCigUI
   ){ /*Loop: get edit distance*/
      if(refSTPtr->cigTypeStr[siRefCig] == 'S')
         break; /*soft masking only at ends (finished)*/

      if(qrySTPtr->cigTypeStr[siQryCig] == 'S')
         break; /*soft masking only at ends (finished)*/

      if(refSTPtr->cigTypeStr[siRefCig] == 'H')
         break; /*hard masking only at ends (finished)*/

      if(qrySTPtr->cigTypeStr[siQryCig] == 'H')
         break; /*hard masking only at ends (finished)*/

      /*`````````````````````````````````````````````````\
      ` Fun06 Sec05 Note01:
      `   - table of cigar entries added together
      \`````````````````````````````````````````````````*/

      /*Values:
      ` =: 61 (match)
      ` X: 88 (snp)
      ` M: 77 (match/snp)
      ` D: 68 (deletion)
      ` I: 73 (insertion)
      ` M + M: 154
      ` = + =: 122
      ` X + X: 176 (need to check if snps agree/disagree)
      ` I + I: 146 (check size; rm low q-scores)
      ` D + D: 136 (check size)
      `
      ` X + M: 165 (need to check if agree/disagree)
      ` = + M: 138 (need to check if agree/disagree)
      ` = + X: 149 (disagree, check q-scores)
      `
      ` D + M: 145 (disagree, check size)
      ` D + =: 129 (disagree, check size)
      ` D + X: 156 (disagree, check siz)
      `
      ` I + M: 150
      ` I + =: 134
      ` I + X: 161
      ` I + D: 141
      */

      /**************************************************\
      * Fun06 Sec05 Sub03:
      *   - start switch/check match/possible snp combos
      *   o fun06 sec05 sub03 cat01:
      *     - get number shared snps/matchs + start loop
      *   o fun06 sec05 sub03 cat02:
      *     - match, anonymous and q-score filter
      *   o fun06 sec05 sub03 cat03:
      *     - depth array filtering (if using)
      *   o fun06 sec05 sub03 cat05:
      *     - add base to window and edit distance
      *   o fun06 sec05 sub03 cat04:
      *     - move to next nucleotide
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun06 Sec05 Sub03 Cat01:
      +   - get number shared snps/matchs + start loop
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      switch(
           refSTPtr->cigTypeStr[siRefCig]
         + qrySTPtr->cigTypeStr[siQryCig]
      ){ /*Switch: find mutation combination*/
         /*cases were I have to compare sequences*/
         case 154: /*M/M*/
         case 176: /*X/X*/
         case 165: /*X/M*/
         case 138: /*=/M*/
         case 149: /*=/X*/
         /*Case: known snp or could be snp*/
            tmpSI =
               min_genMath(
                  refValSI,
                  qryValSI
               );

            resSTPtr->overlapUI += tmpSI;

            while(tmpSI > 0)
            { /*Loop: count number of SNPs*/

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun06 Sec05 Sub03 Cat02:
               +   - match, anonymous and q-score filter
               \++++++++++++++++++++++++++++++++++++++++*/

               if(
                     refSTPtr->seqStr[uiRef]
                  == qrySTPtr->seqStr[uiQry]
               ) goto nextNt_fun06_sec05_sub03_cat05;

               else if(refSTPtr->seqStr[uiRef] == 'N')
                  goto nextNt_fun06_sec05_sub03_cat05;

               else if(qrySTPtr->seqStr[uiQry] == 'N')
                  goto nextNt_fun06_sec05_sub03_cat05;

               if(qryQBl)
               { /*If: have query q-score entry*/
                  qUC =
                       qrySTPtr->qStr[uiQry]
                     - def_adjQ_samEntry;

                  if(qUC < minQUC)
                     goto nextNt_fun06_sec05_sub03_cat05;
               } /*If: have query q-score entry*/

               if(refQBl)
               { /*If: have reference q-score entry*/
                  qUC =
                       refSTPtr->qStr[uiRef]
                     - def_adjQ_samEntry;

                  if(qUC < minQUC)
                     goto nextNt_fun06_sec05_sub03_cat05;
               } /*If: have reference q-score entry*/

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun06 Sec05 Sub03 Cat03:
               +   - depth array filtering (if using)
               \++++++++++++++++++++++++++++++++++++++++*/

               if(resSTPtr->depthAryUI)
               { /*If: using a depth profile*/
                  if(
                       resSTPtr->depthAryUI[uiRef]
                     < minDepthUI
                  ) goto nextNt_fun06_sec05_sub03_cat05;
               } /*If: using a depth profile*/

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun06 Sec05 Sub03 Cat04:
               +   - add base to window and edit distance
               \++++++++++++++++++++++++++++++++++++++++*/

               ++resSTPtr->edDistSL;
               ++winDistUI;

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun06 Sec05 Sub03 Cat05:
               +   - move to next nucleotide
               \++++++++++++++++++++++++++++++++++++++++*/

               nextNt_fun06_sec05_sub03_cat05:;

               if(baseInWinUI >= winSizeUI)
               { /*If: completed a window*/
                  ++resSTPtr->numFullWinUI;

                  resSTPtr->maxWinDistUI =
                     max_genMath(
                        resSTPtr->maxWinDistUI,
                        winDistUI
                     );

                  resSTPtr->minWinDistUI =
                     min_genMath(
                        resSTPtr->minWinDistUI,
                        winDistUI
                     );

                  winDistUI = 0;
                  baseInWinUI = 0;
               } /*If: completed a window*/

               ++baseInWinUI;
               ++uiQry;
               ++uiRef;

               --refValSI;
               --qryValSI;

               --tmpSI;
            } /*Loop: count number of SNPs*/

            break;
         /*Case: known snp or could be snp*/

         /***********************************************\
         * Fun06 Sec05 Sub04:
         *   - check matches or same error type
         \***********************************************/

         /*cases I treat as same*/
         case 122: /*=/=*/
         case 146: /*I/I*/
         case 136: /*D/D*/
         /*Case: treat as same*/
            tmpSI =
               min_genMath(
                  refValSI,
                  qryValSI
               );

            refValSI -= tmpSI;
            qryValSI -= tmpSI;

            if(refSTPtr->cigTypeStr[siRefCig] != 'D')
            { /*If: not deletion event*/
               uiRef += tmpSI;
               uiQry += tmpSI;
               resSTPtr->overlapUI += tmpSI;

               baseInWinUI += tmpSI;

               while(baseInWinUI >= winSizeUI)
               { /*If: complete a window*/
                  ++resSTPtr->numFullWinUI;

                  resSTPtr->maxWinDistUI =
                     max_genMath(
                        resSTPtr->maxWinDistUI,
                        winDistUI
                     );

                  resSTPtr->minWinDistUI =
                     min_genMath(
                        resSTPtr->minWinDistUI,
                        winDistUI
                     );

                  winDistUI = 0;
                  baseInWinUI -= winSizeUI;
               } /*If: complete a window*/
            } /*If: not deletion event*/

            break;
         /*Case: treat as same*/
 
         /***********************************************\
         * Fun06 Sec05 Sub05:
         *   - handel deletion cases
         *   o fun06 sec05 sub05 cat01:
         *     - get number of bases in deletion
         *   o fun06 sec05 sub05 cat02:
         *     - check if keeping deletion
         *   o fun06 sec05 sub05 cat03:
         *     - handle query is insertion for windows
         *   o fun06 sec05 sub05 cat04:
         *     - handle query is deletion for windows
         *   o fun06 sec05 sub05 cat05:
         *     - move to next entries
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun06 Sec05 Sub05 Cat01:
         +   - get number of bases in deletion
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         /*deletion cases*/
         case 145: /*D/M*/
         case 129: /*D/=*/
         case 156: /*D/X*/
         /*Case: deletion (only one sequence)*/
            tmpSI =
               min_genMath(
                  refValSI,
                  qryValSI
               );

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun06 Sec05 Sub05 Cat02:
            +   - check if keeping deletion
            \+++++++++++++++++++++++++++++++++++++++++++*/

            indelBl = 0;

            if(tmpSI >= (sint) indelLenUI)
            { /*If: keeping deletion event*/
               if(
                     resSTPtr->depthAryUI
                  && refSTPtr->cigTypeStr[siRefCig] != 'D'
               ){ /*If: using a depth profile*/
                  if(
                       resSTPtr->depthAryUI[uiRef]
                     < minDepthUI
                  ) goto nextNt_fun06_sec05_sub05_cat03;
               } /*If: using a depth profile*/

               indelBl = 1;
               resSTPtr->edDistSL += tmpSI;
               ++(resSTPtr->indelEventsUI);
               resSTPtr->numIndelUI += tmpSI;
            } /*If: keeping deletion event*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun06 Sec05 Sub05 Cat03:
            +   - handle query is insertion for windows
            \+++++++++++++++++++++++++++++++++++++++++++*/

            nextNt_fun06_sec05_sub05_cat03:;

            if(qrySTPtr->cigTypeStr[siQryCig] == 'D')
            { /*If: query was deletion (ref ins)*/
               uiRef += tmpSI;

               if(indelBl)
               { /*If: keeping indel*/
                  if(baseInWinUI + tmpSI >= winSizeUI)
                  { /*If: indels go beond window size*/
                     indelAddUI = winSizeUI - baseInWinUI;
                     winDistUI += indelAddUI;
                     indelAddUI = tmpSI - indelAddUI;
                  } /*If: indels go beond window size*/

                  else
                  { /*Else: all indels in window*/
                     indelAddUI = 0;
                     winDistUI += tmpSI;
                  } /*Else: all indels in window*/
               } /*If: keeping indel*/

               baseInWinUI += tmpSI;

               while(baseInWinUI >= winSizeUI)
               { /*If: complete a window*/
                  ++resSTPtr->numFullWinUI;

                  resSTPtr->maxWinDistUI =
                     max_genMath(
                        resSTPtr->maxWinDistUI,
                        winDistUI
                     );

                  resSTPtr->minWinDistUI =
                     min_genMath(
                        resSTPtr->minWinDistUI,
                        winDistUI
                     );

                  if(indelBl)
                  { /*If: keeping indels*/
                     if(indelAddUI < winSizeUI)
                        winDistUI = indelAddUI;
                     else
                     { /*Else: changes out of bounds*/
                        winDistUI = winSizeUI;
                        indelAddUI -= winSizeUI;
                     } /*Else: changes out of bounds*/
                  } /*If: keeping indels*/

                  else
                     winDistUI = 0;

                  baseInWinUI -= winSizeUI;
               } /*If: complete a window*/
            } /*If: query was deletion (ref ins)*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun06 Sec05 Sub05 Cat04:
            +   - handle query is deletion for windows
            \+++++++++++++++++++++++++++++++++++++++++++*/

            else
            { /*Else: reference was a deletion*/
               if(indelBl)
                  winDistUI += tmpSI;

               uiQry += tmpSI;
               resSTPtr->overlapUI += tmpSI;
            } /*Else: reference was a deletion*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun06 Sec05 Sub05 Cat05:
            +   - move to next entries
            \+++++++++++++++++++++++++++++++++++++++++++*/

            qryValSI -= tmpSI;
            refValSI -= tmpSI;

            break;
         /*Case: deletion (only one sequence)*/

         /***********************************************\
         * Fun06 Sec05 Sub06:
         *   - handel insertion cases
         *   o fun06 sec05 sub07 cat01:
         *     - check if reference or query is insertion
         *   o fun06 sec05 sub07 cat02:
         *     - check if keeping insertion
         *   o fun06 sec05 sub07 cat03:
         *     - windows when query is insertion
         *   o fun06 sec05 sub07 cat04:
         *     - windows when reference is insertion
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun06 Sec05 Sub07 Cat01:
         +   - check if reference or query is insertion
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 150: /*I/M*/
         case 134: /*I/=*/
         case 161: /*I/X*/
         case 141: /*I/D*/
         /*Case: insertion (only one sequence)*/

            if(refSTPtr->cigTypeStr[siRefCig] == 'I')
            { /*If: reference had insertion*/
               tmpSI = refSTPtr->cigArySI[siRefCig];
               uiRef += tmpSI;
               refValSI = 0;
            } /*If: reference had insertion*/

            else
            { /*Else: query had insertion*/
               tmpSI = qrySTPtr->cigArySI[siQryCig];
               uiQry += tmpSI;
               qryValSI = 0;
            } /*Else: query had insertion*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun06 Sec05 Sub07 Cat02:
            +   - check if keeping insertion
            \+++++++++++++++++++++++++++++++++++++++++++*/

            indelBl = 0;

            if(tmpSI >= (sint) indelLenUI)
            { /*If: keeping query deletion*/
               if(
                     resSTPtr->depthAryUI
                  && refSTPtr->cigTypeStr[siRefCig] == 'I'
               ){ /*If: ref is insertion*/
                if(
                     resSTPtr->depthAryUI[uiRef]
                   < minDepthUI
                 ) goto ignoreIns_fun06_sec05_sub06_cat03;
               } /*If: ref is insertion*/

               resSTPtr->edDistSL += tmpSI;
               indelBl = 1;
               ++(resSTPtr->indelEventsUI);
               resSTPtr->numIndelUI += tmpSI;
            } /*If: keeping query deletion*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun06 Sec05 Sub07 Cat04:
            +   - windows when query is insertion
            \+++++++++++++++++++++++++++++++++++++++++++*/

            ignoreIns_fun06_sec05_sub06_cat03:;

            if(refSTPtr->cigTypeStr[siRefCig] == 'I')
            { /*If: reference is insertion*/
               if(indelBl)
               { /*If: keeping indel*/
                  if(baseInWinUI + tmpSI >= winSizeUI)
                  { /*If: indels go beond window size*/
                     indelAddUI = winSizeUI - baseInWinUI;
                     winDistUI += indelAddUI;
                     indelAddUI = tmpSI - indelAddUI;
                  } /*If: indels go beond window size*/

                  else
                  { /*Else: all indels in window*/
                     indelAddUI = 0;
                     winDistUI += tmpSI;
                  } /*Else: all indels in window*/
               } /*If: keeping indel*/

               baseInWinUI += tmpSI;

               while(baseInWinUI >= winSizeUI)
               { /*If: complete a window*/
                  ++resSTPtr->numFullWinUI;

                  resSTPtr->maxWinDistUI =
                     max_genMath(
                        resSTPtr->maxWinDistUI,
                        winDistUI
                     );

                  resSTPtr->minWinDistUI =
                     min_genMath(
                        resSTPtr->minWinDistUI,
                        winDistUI
                     );

                  if(indelBl)
                  { /*If: keeping indels*/
                     if(indelAddUI < winSizeUI)
                        winDistUI = indelAddUI;
                     else
                     { /*Else: changes out of bounds*/
                        winDistUI = winSizeUI;
                        indelAddUI -= winSizeUI;
                     } /*Else: changes out of bounds*/
                  } /*If: keeping indels*/

                  else
                     winDistUI = 0;

                  baseInWinUI -= winSizeUI;
               } /*If: complete a window*/
            } /*If: reference is insertion*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun06 Sec05 Sub07 Cat05:
            +   - windows when reference is insertion
            \+++++++++++++++++++++++++++++++++++++++++++*/

            else
            { /*Else: query had insertion*/
               if(indelBl)
                  winDistUI += tmpSI;
            } /*Else: query had insertion*/

            break;
         /*Case: insertion (only one sequence)*/
      } /*Switch: find mutation combination*/

      /**************************************************\
      * Fun06 Sec05 Sub07:
      *   - check if move to next cigar entry
      \**************************************************/

      if(refValSI == 0)
      { /*If: on next cigar entry*/
         ++siRefCig;
         refValSI = refSTPtr->cigArySI[siRefCig];
      } /*If: on next cigar entry*/

      if(qryValSI == 0)
      { /*If: on next cigar entry*/
         ++siQryCig;
         qryValSI = qrySTPtr->cigArySI[siQryCig];
      } /*If: on next cigar entry*/
   } /*Loop: get edit distance*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec06:
   ^   - add last window to data
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(baseInWinUI >= winSizeUI)
   { /*If: covered read had no partial windows*/
      ++resSTPtr->numFullWinUI;

      resSTPtr->maxWinDistUI =
         max_genMath(
            resSTPtr->maxWinDistUI,
            winDistUI
         );

      resSTPtr->minWinDistUI =
         min_genMath(
            resSTPtr->minWinDistUI,
            winDistUI
         );

      resSTPtr->extraWinNtUI = 0;
      resSTPtr->lastWinDiffUI = 0;

      resSTPtr->avgWinDistF =
           (float) resSTPtr->edDistSL
         / (float) resSTPtr->numFullWinUI;
   } /*If: covered read had no partial windows*/

   else
   { /*Else: had a partial window*/
      resSTPtr->extraWinNtUI = winSizeUI;
      resSTPtr->lastWinDiffUI = winDistUI;

      resSTPtr->maxWinDistUI =
         max_genMath(
            resSTPtr->maxWinDistUI,
            winDistUI
         );

      resSTPtr->minWinDistUI =
         min_genMath(
            resSTPtr->minWinDistUI,
            winDistUI
         );

      resSTPtr->avgWinDistF =
           (float) resSTPtr->edDistSL
         / ( (float) resSTPtr->numFullWinUI + 1 );
   } /*Else: had a partial window*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec07:
   ^   - return the edit distance
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retSL = resSTPtr->edDistSL;
   goto ret_fun06_sec07;

   noOverlap_fun06_sec07:;
      retSL = def_noOverlap_edDist;
      goto ret_fun06_sec07;

   diffRef_fun06_sec07:;
      retSL = def_diffRef_edDist;
      goto ret_fun06_sec07;

   noSeq_fun06_sec07:;
      retSL = def_noSeq_edDist;
      goto ret_fun06_sec07;

   noMap_fun06_sec07:;
      retSL = def_noMap_edDist;
      goto ret_fun06_sec07;

   ret_fun06_sec07:;
      return retSL;
} /*readCmpDist_edDist*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC:
   '   - gets edit distances between the query & reference
   '   o fun07 sec01:
   '     - variable declerations
   '   o fun07 sec02:
   '     - check if can get edit distance
   '   o fun07 sec03:
   '     - initialize and check if have q-score entreis
   '   o fun07 sec04:
   '     - find edit distance
   '   o fun07 sec06:
   '     - add last window
   '   o fun07 sec07:
   '     - return the edit distance
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint siQryCig = 0;  /*query cigar position*/
   sint qryValSI = 0;
   uint uiQry = 0;     /*query nucleotide on*/
   uint uiRef = 0;

   slong retSL = 0;    /*returned edit distance*/

   schar qryQBl = 0;
   uchar qUC = 0;

   schar indelBl = 0;     /*for adding indesl to windows*/
   uint baseInWinUI = 0;  /*number of bases in window*/
   uint winDistUI = 0;    /*edit distance of window*/
   uint indelAddUI = 0;   /*find offset for indels*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^   - check if can get edit distance
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_res_edDist(resSTPtr);

   if(qrySTPtr->flagUS & 4)
      goto noMap_fun07_sec07;

   if(
         ! qrySTPtr->seqStr
      || *qrySTPtr->seqStr == '*' 
      || *qrySTPtr->seqStr == '\0'
   ) goto noSeq_fun07_sec07;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^   - initialize and check if have q-score entreis
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
         qrySTPtr->qStr[0] == '\0'
      ||
         (
             qrySTPtr->qStr[0] == '*'
          && qrySTPtr->qStr[1] != '\0'
         )
   ) qryQBl = 0;

   else
      qryQBl = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec04:
   ^   - find edit distance
   ^   o fun07 sec04 sub01:
   ^     - get past masking at start
   ^   o fun07 sec04 sub02:
   ^     - start loop and check if at end (hard/soft mask)
   ^   o fun07 sec04 note01:
   ^     - table of cigar entries added together
   ^   o fun07 sec04 sub03:
   ^     - ignore matches (or check if have snps)
   ^   o fun07 sec04 sub04:
   ^     - check if keeping snp
   ^   o fun07 sec04 sub05:
   ^     - handel deletion cases
   ^   o fun07 sec04 sub06:
   ^     - handel insertion cases
   ^   o fun07 sec04 sub07:
   ^     - check if move to next cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun07 Sec04 Sub01:
   *   - get past masking at start
   \*****************************************************/

   if(qrySTPtr->cigTypeStr[siQryCig] == 'S')
   { /*If: have softmasked region*/
      uiQry += qryValSI;
      ++siQryCig;
      qryValSI = qrySTPtr->cigArySI[siQryCig];
   } /*If: have softmasked region*/

   if(qrySTPtr->cigTypeStr[siQryCig] == 'H')
   { /*If: have hardmaksed region*/
      ++siQryCig;
      qryValSI = qrySTPtr->cigArySI[siQryCig];
   } /*If: have hardmasked region*/

   /*****************************************************\
   * Fun07 Sec04 Sub02:
   *   - start loop and check if at end (hard/soft mask)
   \*****************************************************/

   while(siQryCig < (sint) qrySTPtr->lenCigUI)
   { /*Loop: get edit distance*/
      if(qrySTPtr->cigTypeStr[siQryCig] == 'S')
         break; /*soft masking only at ends (finished)*/

      if(qrySTPtr->cigTypeStr[siQryCig] == 'H')
         break; /*hard masking only at ends (finished)*/

      /**************************************************\
      * Fun07 Sec04 Sub03:
      *   - ignore matches (or check if have snps)
      \**************************************************/

      switch(qrySTPtr->cigTypeStr[siQryCig])
      { /*Switch: find mutation combination*/
         /*cases were I have to compare sequences*/
         case 'M':
         /*Case: snp or match*/
            if(refSTPtr)
               goto snp_fun07_sec04_sub04_cat01;
            else
               goto match_fun07_sec04_sub03;

         /*Case: snp or match*/

         case '=':
         /*Case: likley match (or no way to check)*/
            match_fun07_sec04_sub03:;
            uiQry += qryValSI;
            uiRef += qryValSI;
            baseInWinUI += qryValSI;

            while(baseInWinUI >= winSizeUI)
            { /*Loop: save complete a window*/
               ++resSTPtr->numFullWinUI;

               resSTPtr->maxWinDistUI =
                  max_genMath(
                     resSTPtr->maxWinDistUI,
                     winDistUI
                  );

               resSTPtr->minWinDistUI =
                  min_genMath(
                     resSTPtr->minWinDistUI,
                     winDistUI
                  );

               winDistUI = 0;
               baseInWinUI -= winSizeUI;
            } /*Loop: save complete a window*/

            break;
         /*Case: likley match (or no way to check)*/

         /***********************************************\
         * Fun07 Sec04 Sub04:
         *   - check if keeping snp
         *   o fun07 sec04 sub04 cat01:
         *     - check snp quality score + start loop
         *   o fun07 sec04 sub04 cat02:
         *     - ignore transitions (if ignoring)
         *   o fun07 sec04 sub04 cat03:
         *     - count differnce (if kept) & move on
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun07 Sec04 Sub04 Cat01:
         +   - check snp quality score + start loop
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'X':
         /*Case: snp*/
            snp_fun07_sec04_sub04_cat01:;

            while(qryValSI > 0)
            { /*Loop: count number of SNPs*/
               if(qryQBl)
               { /*If: have query q-score entry*/
                  qUC =
                       qrySTPtr->qStr[uiQry]
                     - def_adjQ_samEntry;

                  if(qUC < minQUC)
                     goto nextNt_fun07_sec04_sub03_cat03;
               } /*If: have query q-score entry*/

               if(qrySTPtr->seqStr[uiQry] == ('N' & ~32))
                  goto nextNt_fun07_sec04_sub03_cat03;
                  /*not counting masked bases*/

               if(baseInWinUI >= winSizeUI)
               { /*If: complete a window*/
                  ++resSTPtr->numFullWinUI;

                  resSTPtr->maxWinDistUI =
                     max_genMath(
                        resSTPtr->maxWinDistUI,
                        winDistUI
                     );

                  resSTPtr->minWinDistUI =
                     min_genMath(
                        resSTPtr->minWinDistUI,
                        winDistUI
                     );

                  winDistUI = 0;
                  baseInWinUI = 0;
               } /*If: complete a window*/

               ++resSTPtr->edDistSL;
               ++winDistUI;

              /*+++++++++++++++++++++++++++++++++++++++++\
              + Fun07 Sec04 Sub04 Cat03:
              +   - count differnce (if kept) & move on
              \+++++++++++++++++++++++++++++++++++++++++*/

               nextNt_fun07_sec04_sub03_cat03:;
               ++uiQry;
               ++uiRef;
               --qryValSI;
            } /*Loop: count number of SNPs*/

            break;
         /*Case: snp*/

         /***********************************************\
         * Fun07 Sec04 Sub05:
         *   - handel deletion cases
         *   o fun07 sec04 sub05 cat01:
         *     - check if keeping deletion
         *   o fun07 sec04 sub05 cat02:
         *     - find if deletion overflows window
         *   o fun07 sec04 sub05 cat03:
         *     - update windows until pariatial window
         *   o fun07 sec04 sub05 cat04:
         *     - move to next reference base
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun07 Sec04 Sub05 Cat01:
         +   - check if keeping deletion
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         /*deletion cases*/
         case 'D':
         /*Case: deletion*/
            if(qryValSI >= (sint) indelLenUI)
            { /*If: keeping deletion event*/
               resSTPtr->edDistSL += qryValSI;
               ++resSTPtr->indelEventsUI;
               resSTPtr->numIndelUI += qryValSI;
               indelBl = 1; 
            } /*If: keeping deletion event*/

            else
               indelBl = 0; 

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun07 Sec04 Sub05 Cat02:
            +   - find if deletion overflows window
            \+++++++++++++++++++++++++++++++++++++++++++*/

            if(indelBl)
            { /*If: deletion was kept*/
               if(baseInWinUI + qryValSI >= winSizeUI)
               { /*If: indels go beond window size*/
                  indelAddUI = winSizeUI - baseInWinUI;
                  winDistUI += indelAddUI;
                  indelAddUI = qryValSI - indelAddUI;
               } /*If: indels go beond window size*/

               else
               { /*Else: all indels in window*/
                  indelAddUI = 0;
                  winDistUI += qryValSI;
               } /*Else: all indels in window*/
            } /*If: deletion was kept*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun07 Sec04 Sub05 Cat03:
            +   - update windows until pariatial window
            \+++++++++++++++++++++++++++++++++++++++++++*/

            baseInWinUI += qryValSI;

            while(baseInWinUI >= winSizeUI)
            { /*Loop: save complete a window*/
               ++resSTPtr->numFullWinUI;

               resSTPtr->maxWinDistUI =
                  max_genMath(
                     resSTPtr->maxWinDistUI,
                     winDistUI
                  );

               resSTPtr->minWinDistUI =
                  min_genMath(
                     resSTPtr->minWinDistUI,
                     winDistUI
                  );

               if(indelBl)
               { /*If: keeping indels*/
                  if(indelAddUI < winSizeUI)
                     winDistUI = indelAddUI;
                  else
                  { /*Else: changes out of bounds*/
                     winDistUI = winSizeUI;
                     indelAddUI -= winSizeUI;
                  } /*Else: changes out of bounds*/
               } /*If: keeping indels*/

               else
                  winDistUI = 0;

               baseInWinUI -= winSizeUI;
            } /*Loop: save complete a window*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun07 Sec04 Sub05 Cat04:
            +   - move to next reference base
            \+++++++++++++++++++++++++++++++++++++++++++*/

            uiRef += qryValSI;
            break;
         /*Case: deletion*/

         /***********************************************\
         * Fun07 Sec04 Sub06:
         *   - handel insertion cases
         \***********************************************/

         case 'I':
         /*Case: insertion*/
            uiQry += qryValSI;

            if(qryValSI >= (sint) indelLenUI)
            { /*If: keeping query deletion*/
               resSTPtr->edDistSL += qryValSI;
               ++resSTPtr->indelEventsUI;
               resSTPtr->numIndelUI += qryValSI;

               winDistUI += qryValSI; /*not in reference*/
            } /*If: keeping query deletion*/

            break;
         /*Case: insertion*/
      } /*Switch: find mutation combination*/

      /**************************************************\
      * Fun07 Sec04 Sub07:
      *   - check if move to next cigar entry
      \**************************************************/

      ++siQryCig;
      qryValSI = qrySTPtr->cigArySI[siQryCig];
   } /*Loop: get edit distance*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec06:
   ^   - add last window
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(baseInWinUI >= winSizeUI)
   { /*If: covered read had no partial windows*/
      ++resSTPtr->numFullWinUI;

      resSTPtr->maxWinDistUI =
         max_genMath(
            resSTPtr->maxWinDistUI,
            winDistUI
         );

      resSTPtr->minWinDistUI =
         min_genMath(
            resSTPtr->minWinDistUI,
            winDistUI
         );

      resSTPtr->extraWinNtUI = 0;
      resSTPtr->lastWinDiffUI = 0;

      resSTPtr->avgWinDistF =
           (float) resSTPtr->edDistSL
         / (float) resSTPtr->numFullWinUI;
   } /*If: covered read had no partial windows*/

   else
   { /*Else: had a partial window*/
      resSTPtr->extraWinNtUI = winSizeUI;
      resSTPtr->lastWinDiffUI = winDistUI;

      resSTPtr->maxWinDistUI =
         max_genMath(
            resSTPtr->maxWinDistUI,
            winDistUI
         );

      resSTPtr->minWinDistUI =
         min_genMath(
            resSTPtr->minWinDistUI,
            winDistUI
         );

      resSTPtr->avgWinDistF =
           (float) resSTPtr->edDistSL
         / ( (float) resSTPtr->numFullWinUI + 1 );
   } /*Else: had a partial window*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec07:
   ^   - return the edit distance
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   resSTPtr->startUI = qrySTPtr->refStartUI;
   resSTPtr->overlapUI = qrySTPtr->alnReadLenUI;

   retSL = resSTPtr->edDistSL;
   goto ret_fun07_sec07;

   noSeq_fun07_sec07:;
      retSL = def_noSeq_edDist;
      goto ret_fun07_sec07;

   noMap_fun07_sec07:;
      retSL = def_noMap_edDist;
      goto ret_fun07_sec07;

   ret_fun07_sec07:;
      return retSL;
} /*dist_edDist*/

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
){
   resSTPtr->probUI =
         100
      * (
            resSTPtr->edDistSL
          / (resSTPtr->overlapUI * percErrF)
        );

   resSTPtr->probMaxWinUI =
         100
      * (
            resSTPtr->maxWinDistUI
          / (winSizeUI * percErrF)
        );

   return resSTPtr->probUI;
} /*percDist_edDist*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun09 TOC:
   '   - adds a read to a depth profile
   '   o fun09 sec01:
   '     - variable declarations
   '   o fun09 sec02:
   '     - filter reads
   '   o fun09 sec03:
   '     - prepare for comparisions
   '   o fun09 sec04:
   '     - add depths
   '   o fun09 sec05:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   float overlapF = 0;  /*finding overlap*/

   sint siQryCig = 0;  /*query cigar position*/
   sint siRefCig = 0;  /*reference cigar position*/
   sint refValSI = 0;  /*number reference bases in cigar*/
   sint qryValSI = 0;  /*number query bases in cigar*/

   sint tmpSI = 0;

   uint uiQry = 0;     /*query nucleotide on*/
   uint uiRef = 0;     /*reference nucleotide on*/
   uchar qUC = 0;      /*for q-score finding*/

   schar qryQBl = 0;    /*tells if query has q-score*/
   schar refQBl = 0;    /*tells if ref has q-score*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec02:
   ^   - filter reads
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(qrySTPtr->extraStr[0] == '@')
      goto discard_fun09_sec05;

   if(qrySTPtr->flagUS & 4)
      goto discard_fun09_sec05;
   
   if(qrySTPtr->refEndUI < refSTPtr->refStartUI)
      goto discard_fun09_sec05;

   if(qrySTPtr->refStartUI > refSTPtr->refEndUI)
      goto discard_fun09_sec05;

   if(qrySTPtr->refEndUI < refSTPtr->refStartUI)
      goto discard_fun09_sec05;

   if(minOverlapF > 0)
   { /*If: user wanted overlap checks*/
      overlapF =
        (float)
        min_genMath(
           qrySTPtr->refEndUI,
           refSTPtr->refEndUI
        );

      overlapF -=
        (float)
        max_genMath(
           qrySTPtr->refStartUI,
           refSTPtr->refStartUI
        );

      if(
                overlapF / (float) refSTPtr->alnReadLenUI
              < minOverlapF
         &&
                overlapF / (float) qrySTPtr->alnReadLenUI
              < minOverlapF
      ) goto discard_fun09_sec05;
   } /*If: user wanted overlap checks*/

   if(
         ! qrySTPtr->seqStr
      || *qrySTPtr->seqStr == '*' 
      || *qrySTPtr->seqStr == '\0'
   ) goto discard_fun09_sec05;

   if(
         ! refSTPtr->seqStr
      || *refSTPtr->seqStr == '*' 
      || *refSTPtr->seqStr == '\0'
   ) goto discard_fun09_sec05;

   if(
      eql_ulCp(
         qrySTPtr->refIdStr,
         refSTPtr->refIdStr,
         0,
         '\0'
      )
   ) goto discard_fun09_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec03:
   ^   - prepare for comparisions
   ^   o fun09 sec03 sub01:
   ^     - check if have q-score entries
   ^   o fun09 sec03 sub02:
   ^     - find start positions
   ^   o fun09 sec03 sub03:
   ^     - move past softmasking
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun09 Sec03 Sub01:
   *   - check if have q-score entries
   \*****************************************************/

   if(
         qrySTPtr->qStr[0] == '\0'
      ||
         (
             qrySTPtr->qStr[0] == '*'
          && qrySTPtr->qStr[1] == '\0'
         )
   ) qryQBl = 0;

   else
      qryQBl = 1;

   if(
         refSTPtr->qStr[0] == '\0'
      ||
         (
             refSTPtr->qStr[0] == '*'
          && refSTPtr->qStr[1] == '\0'
         )
   ) refQBl = 0;

   else
      refQBl = 1;

   /*****************************************************\
   * Fun09 Sec03 Sub02:
   *   - find start positions
   \*****************************************************/

   resSTPtr->startUI =
      max_genMath(
         refSTPtr->refStartUI,
         qrySTPtr->refStartUI
      );

   refValSI = refSTPtr->cigArySI[0];
   qryValSI = qrySTPtr->cigArySI[0];

   if(qrySTPtr->refStartUI > refSTPtr->refStartUI)
   { /*If: move the reference foward*/
      uiQry = refSTPtr->refStartUI;

      findRefPos_samEntry(
         refSTPtr,
         &siRefCig,
         &refValSI,
         qrySTPtr->refStartUI,
         (sint *) &uiQry, /*reference pos (discard)*/
         (sint *) &uiRef
      ); /*set reference to first query base*/

      uiQry = 0;
   } /*If: move reference foward*/

   else if(qrySTPtr->refStartUI < refSTPtr->refStartUI)
   { /*Else If: move query foward*/
      uiRef = qrySTPtr->refStartUI;

      findRefPos_samEntry(
         qrySTPtr,
         &siQryCig,
         &qryValSI,
         refSTPtr->refStartUI, /*end query position*/
         (sint *) &uiRef, /*query pos (discard)*/
         (sint *) &uiQry  /*end query position*/
      ); /*set reference to first query base*/

      uiRef = 0;
   } /*Else If: move query foward*/

   /*****************************************************\
   * Fun09 Sec03 Sub03:
   *   - move past softmasking
   \*****************************************************/

   if(refSTPtr->cigTypeStr[siRefCig] == 'S')
   { /*If: have softmasked region*/
      uiRef += refValSI;
      ++siRefCig;
      refValSI = refSTPtr->cigArySI[siRefCig];
   } /*If: have softmasked region*/

   if(qrySTPtr->cigTypeStr[siQryCig] == 'S')
   { /*If: have softmasked region*/
      uiQry += qryValSI;
      ++siQryCig;
      qryValSI = qrySTPtr->cigArySI[siQryCig];
   } /*If: have softmasked region*/

   if(refSTPtr->cigTypeStr[siRefCig] == 'H')
   { /*If: have hardmasked region*/
      ++siRefCig;
      refValSI = refSTPtr->cigArySI[siRefCig];
   } /*If: have hardmasked region*/

   if(qrySTPtr->cigTypeStr[siQryCig] == 'H')
   { /*If: have hardmaksed region*/
      ++siQryCig;
      qryValSI = qrySTPtr->cigArySI[siQryCig];
   } /*If: have hardmasked region*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec04:
   ^   - add depths
   ^   o fun09 sec04 sub01:
   ^     - start loop + break if softmasking at end
   ^   o fun09 sec04 note01:
   ^     - table of cigar entries added together
   ^   o fun09 sec04 sub02:
   ^     - start switch/check match/possible snp combos
   ^   o fun09 sec04 sub03:
   ^     - check ins/ins or del/del
   ^   o fun09 sec04 sub04:
   ^     - handel deletion cases
   ^   o fun09 sec04 sub05:
   ^     - handel insertion cases
   ^   o fun09 sec04 sub06:
   ^     - check if move to next cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun09 Sec04 Sub01:
   *   - start loop + break if softmasking at end
   \*****************************************************/

   while(
         siQryCig < (sint) qrySTPtr->lenCigUI
      && siRefCig < (sint) refSTPtr->lenCigUI
   ){ /*Loop: add depths*/
      if(refSTPtr->cigTypeStr[siRefCig] == 'S')
         break; /*soft masking only at ends (finished)*/

      if(qrySTPtr->cigTypeStr[siQryCig] == 'S')
         break; /*soft masking only at ends (finished)*/

      if(refSTPtr->cigTypeStr[siRefCig] == 'H')
         break; /*hard masking only at ends (finished)*/

      if(qrySTPtr->cigTypeStr[siQryCig] == 'H')
         break; /*hard masking only at ends (finished)*/

      /*`````````````````````````````````````````````````\
      ` Fun09 Sec04 Note01:
      `   - table of cigar entries added together
      \`````````````````````````````````````````````````*/

      /*Values:
      ` =: 61 (match)
      ` X: 88 (snp)
      ` M: 77 (match/snp)
      ` D: 68 (deletion)
      ` I: 73 (insertion)
      ` M + M: 154
      ` = + =: 122
      ` X + X: 176 (need to check if snps agree/disagree)
      ` I + I: 146 (check size; rm low q-scores)
      ` D + D: 136 (check size)
      `
      ` X + M: 165 (need to check if agree/disagree)
      ` = + M: 138 (need to check if agree/disagree)
      ` = + X: 149 (disagree, check q-scores)
      `
      ` D + M: 145 (disagree, check size)
      ` D + =: 129 (disagree, check size)
      ` D + X: 156 (disagree, check siz)
      `
      ` I + M: 150
      ` I + =: 134
      ` I + X: 161
      ` I + D: 141
      */

      /**************************************************\
      * Fun09 Sec04 Sub02:
      *   - start switch/check match/possible snp combos
      *   o fun09 sec04 sub02 cat01:
      *     - find number shared snps/matchs + check case
      *   o fun09 sec04 sub02 cat02:
      *     - find snps passing min quality scores
      *   o fun09 sec04 sub02 cat03:
      *     - check if SNP (passed) is a transversion
      *   o fun09 sec04 sub02 cat04:
      *     - move to next SNP/match
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun09 Sec04 Sub02 Cat01:
      +   - find number shared snps/matchs + check case
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      switch(
           refSTPtr->cigTypeStr[siRefCig]
         + qrySTPtr->cigTypeStr[siQryCig]
      ){ /*Switch: find mutation combination*/
         /*cases were I have to compare sequences*/
         case 154: /*M/M*/
         case 176: /*X/X*/
         case 165: /*X/M*/
         case 138: /*=/M*/
         case 149: /*=/X*/
         case 122: /*=/=*/
         /*Case: match or snp*/
            tmpSI =
               min_genMath(
                  refValSI,
                  qryValSI
               );

            while(tmpSI > 0)
            { /*Loop: add depth for matches*/
               if(
                     refSTPtr->seqStr[uiRef]
                  != qrySTPtr->seqStr[uiQry]
               ) goto nextNt_fun09_sec04_sub02_cat02;

               else if(refSTPtr->seqStr[uiRef] == 'N')
                  goto nextNt_fun09_sec04_sub02_cat02;

               else if(qrySTPtr->seqStr[uiQry] == 'N')
                  goto nextNt_fun09_sec04_sub02_cat02;

               if(qryQBl)
               { /*If: have query q-score entry*/
                  qUC =
                       qrySTPtr->qStr[uiQry]
                     - def_adjQ_samEntry;

                  if(qUC < minQUC)
                     goto nextNt_fun09_sec04_sub02_cat02;
               } /*If: have query q-score entry*/

               if(refQBl)
               { /*If: have reference q-score entry*/
                  qUC =
                       refSTPtr->qStr[uiRef]
                     - def_adjQ_samEntry;

                  if(qUC < minQUC)
                     goto nextNt_fun09_sec04_sub02_cat02;
               } /*If: have reference q-score entry*/

               ++resSTPtr->depthAryUI[uiRef];

               nextNt_fun09_sec04_sub02_cat02:;

               ++uiQry;
               ++uiRef;

               --refValSI;
               --qryValSI;

               --tmpSI;
            } /*Loop: add depth for matches*/

            break;
         /*Case: match or snp*/

         /***********************************************\
         * Fun09 Sec04 Sub03:
         *   - check ins/ins or del/del
         \***********************************************/

         /*cases I treat as same*/
         case 146: /*I/I*/
         case 136: /*D/D*/
         /*Case: treat as same*/
            tmpSI =
               min_genMath(
                  refValSI,
                  qryValSI
               );

            refValSI -= tmpSI;
            qryValSI -= tmpSI;

            if(refSTPtr->cigTypeStr[siRefCig] != 'D')
            { /*If: not deletion event*/
               uiQry += tmpSI;

               while(tmpSI > 0)
               { /*Loop: add in depths*/
                  ++resSTPtr->depthAryUI[uiRef];
                  --tmpSI;
                  ++uiRef;
               } /*Loop: add in depths*/
            } /*If: not deletion event*/

            break;
         /*Case: treat as same*/
 
         /***********************************************\
         * Fun09 Sec04 Sub04:
         *   - handel deletion cases
         \***********************************************/

         /*deletion cases*/
         case 145: /*D/M*/
         case 129: /*D/=*/
         case 156: /*D/X*/
         /*Case: deletion (only one sequence)*/
            tmpSI =
               min_genMath(
                  refValSI,
                  qryValSI
               );

            if(refSTPtr->cigTypeStr[siRefCig] == 'D')
               uiQry += tmpSI;
            else
               uiRef += tmpSI;

            qryValSI -= tmpSI;
            refValSI -= tmpSI;

            break;
         /*Case: deletion (only one sequence)*/

         /***********************************************\
         * Fun09 Sec04 Sub05:
         *   - handel insertion cases
         \***********************************************/

         case 150: /*I/M*/
         case 134: /*I/=*/
         case 161: /*I/X*/
         case 141: /*I/D*/
         /*Case: insertion (only one sequence)*/

            if(refSTPtr->cigTypeStr[siRefCig] == 'I')
            { /*If: reference had insertion*/
               tmpSI = refSTPtr->cigArySI[siRefCig];
               uiRef += tmpSI;
               refValSI = 0;
            } /*If: reference had insertion*/

            else
            { /*Else: query had insertion*/
               tmpSI = qrySTPtr->cigArySI[siQryCig];
               uiQry += tmpSI;
               qryValSI = 0;
            } /*Else: query had insertion*/

            break;
         /*Case: insertion (only one sequence)*/
      } /*Switch: find mutation combination*/

      /**************************************************\
      * Fun09 Sec04 Sub06:
      *   - check if move to next cigar entry
      \**************************************************/

      if(refValSI == 0)
      { /*If: on next cigar entry*/
         ++siRefCig;
         refValSI = refSTPtr->cigArySI[siRefCig];
      } /*If: on next cigar entry*/

      if(qryValSI == 0)
      { /*If: on next cigar entry*/
         ++siQryCig;
         qryValSI = qrySTPtr->cigArySI[siQryCig];
      } /*If: on next cigar entry*/
   } /*Loop: add depths*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec05:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   discard_fun09_sec05:;
   return 1;
} /*addReadToDepth_edDist*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun10 TOC:
   '   - finds depth of every base in samEntry reference
   '   o fun10 sec01:
   '     - variable declarations
   '   o fun10 sec02:
   '     - allocate memory and get first line
   '   o fun10 sec03:
   '     - get depths
   '   o fun10 sec04:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec02:
   ^   - allocate memory and get first line
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(resSTPtr->sizeDepthUL < refSTPtr->readLenUI)
   { /*If: need more memory*/
      free(resSTPtr->depthAryUI);
      resSTPtr->depthAryUI = 0;
   } /*If: need more memory*/

   if(! resSTPtr->depthAryUI)
   { /*If: need to allocate memory for the depth array*/
      resSTPtr->depthAryUI =
         calloc(
            (refSTPtr->readLenUI + 9),
            sizeof(uint)
         );

      if(! resSTPtr->depthAryUI)
         goto memErr_fun10_sec04_sub02;

      resSTPtr->sizeDepthUL = refSTPtr->readLenUI;
   } /*If: need to allocate memory for the depth array*/

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
         goto memErr_fun10_sec04_sub02;

      else
         goto fileErr_fun10_sec04_sub03;
   } /*If: had error*/


   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec03:
   ^   - get depths
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(! errSC)
   { /*Loop: get depths*/
      addReadToDepth_edDist(
         refSTPtr,
         samSTPtr,
         minQUC,
         minOverlapF,
         resSTPtr
      );

      errSC =
         get_samEntry(
            samSTPtr,
            buffStrPtr,
            lenBuffULPtr,
            samFILE
         );
   } /*Loop: get depths*/

   if(errSC == def_memErr_samEntry)
      goto memErr_fun10_sec04_sub02;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec04:
   ^   - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun10_sec04_sub05;

   memErr_fun10_sec04_sub02:;
      errSC = def_memErr_edDist;
      goto ret_fun10_sec04_sub05;

   fileErr_fun10_sec04_sub03:;
      errSC = def_fileErr_edDist;
      goto ret_fun10_sec04_sub05;

   ret_fun10_sec04_sub05:;
      fseek(
         samFILE,
         0,
         SEEK_SET
      );

      return errSC;
} /*mkDepthProfile_edDist*/

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
){

   fprintf(
      (FILE *) outFILE,
     "qry\tref\tdist\tdist_div_err\tstart\taln_len"
   );

   fprintf(
      (FILE *) outFILE,
     "\tnum_indels\tindel_events\tmin_win\tavg_win"
   );

   fprintf(
      (FILE *) outFILE,
     "\tmax_win\tmax_win_div_err\n"
   );
} /*phead_edDist*/

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
){

   fprintf(
      outFILE,
      "%s\t%s",
      qryIdStr,
      refIdStr
   );

   fprintf(
      outFILE,
      "\t%li\t%i\t%u\t%u\t%u\t%u",
      resSTPtr->edDistSL,
      resSTPtr->probUI,
      resSTPtr->startUI,
      resSTPtr->overlapUI,
      resSTPtr->numIndelUI,
      resSTPtr->indelEventsUI
   );

   fprintf(
      outFILE,
      "\t%u\t%0.3f\t%u\t%u\n",
      resSTPtr->minWinDistUI,
      resSTPtr->avgWinDistF,
      resSTPtr->maxWinDistUI,
      resSTPtr->probMaxWinUI
   );
} /*pdist_edDist*/
   
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
