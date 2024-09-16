/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' edDist SOF: Start Of File
'   - has functions to find edit distances (does some
'     filtering)
'   o header:
'     - included libraries
'   o .c fun01: isTransNt_edDist
'     - checks to see if query base is a transversion
'     - not used (manually inlined), but here for future
'   o fun02: readCmpDist_edDist
'     - gets edit distances between the query & reference
'   o fun03: dist_edDist
'     - gets edit distances for reference
'   o fun04: distProb_edDist
'     - finds probablity of reads matching by chance
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
#include "../genLib/charCp.h"
#include "../genLib/samEntry.h"

/*.h files only*/
#include "../genLib/dataTypeShortHand.h"
#include "../genLib/genMath.h" /*min .h macro only*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .c   #include "../geLib/base10str.h"
!   o .c   #include "../geLib/numToStr.h"
!   o .h   #include "../geLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: isTransNt_edDist
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
| Fun02: readCmpDist_edDist
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
|   - overlapUI:
|     o pointer to unsigned int to hold nubmer of
|       reference base covered by query
|   - numIndelUI:
|     o pointer to unisigned int to hold the number of
|       indels in edit distance
|   - indelEventsUI:
|     o pointer to unsigned int to hold the number of
|       indel events (times a group of indels occured)
| Output:
|   - Modifies:
|     o numIndelUI to have number indels kept
|     o indelEventsUI to have number of indel events
|     o overlapUI to hold number of overlapped bases
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
   unsigned int *overlapUI,   /*overlap length*/
   unsigned int *numIndelUI,  /*number indels kept*/
   unsigned int *indelEventsUI /*nubmer indel events*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - gets edit distances between the query & reference
   '   o fun02 sec01:
   '     - variable declerations
   '   o fun02 sec02:
   '     - check if can get edit distance
   '   o fun02 sec03:
   '     - initialize and check if have q-score entreis
   '   o fun02 sec04:
   '     - find start of overlap between query & reference
   '   o fun02 sec05:
   '     - find edit distance
   '   o fun02 sec06:
   '     - return the edit distance
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint siQryCig = 0;  /*query cigar position*/
   sint siRefCig = 0;  /*reference cigar position*/
   sint refValSI = 0;
   sint qryValSI = 0;

   sint tmpSI = 0;

   uint uiQry = 0;     /*query nucleotide on*/
   uint uiRef = 0;     /*reference nucleotide on*/

   slong distSL = 0;    /*edit distance*/

   schar qryQBl = 0;
   schar refQBl = 0;
   uchar qUC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - check if can get edit distance
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *overlapUI = 0;
   *indelEventsUI = 0;
   *numIndelUI = 0;

   if(qrySTPtr->flagUS & 4)
      goto noMap_fun02_sec06;

   if(refSTPtr->flagUS & 4)
      goto noMap_fun02_sec06;

   if(qrySTPtr->refEndUI < refSTPtr->refStartUI)
      goto noOverlap_fun02_sec06;

   if(qrySTPtr->refStartUI > refSTPtr->refEndUI)
      goto noOverlap_fun02_sec06;

   if(
         ! qrySTPtr->seqStr
      || *qrySTPtr->seqStr == '*' 
      || *qrySTPtr->seqStr == '\0'
   ) goto noSeq_fun02_sec06;

   if(
         ! refSTPtr->seqStr
      || *refSTPtr->seqStr == '*' 
      || *refSTPtr->seqStr == '\0'
   ) goto noSeq_fun02_sec06;

   if(
      eql_charCp(
         qrySTPtr->refIdStr,
         refSTPtr->refIdStr,
         '\0'
      )
   ) goto diffRef_fun02_sec06;
     /*reads are mapped to different references*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
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
   ^ Fun02 Sec04:
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
      uiRef = refSTPtr->refStartUI;

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
   ^ Fun02 Sec05:
   ^   - find edit distance
   ^   o fun02 sec05 sub01:
   ^     - get past masking at start
   ^   o fun02 sec05 sub02:
   ^     - start loop and check if at end (hard/soft mask)
   ^   o fun02 sec05 note01:
   ^     - table of cigar entries added together
   ^   o fun02 sec05 sub03:
   ^     - start switch/check match/possible snp combos
   ^   o fun02 sec05 sub04:
   ^     - check matches or same error type
   ^   o fun02 sec05 sub05:
   ^     - handel deletion cases
   ^   o fun02 sec05 sub06:
   ^     - handel insertion cases
   ^   o fun02 sec05 sub07:
   ^     - check if move to next cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec05 Sub01:
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
   * Fun02 Sec05 Sub02:
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
      ` Fun02 Sec05 Note01:
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
      * Fun02 Sec05 Sub03:
      *   - start switch/check match/possible snp combos
      *   o fun02 sec05 sub03 cat01:
      *     - find number shared snps/matchs + check case
      *   o fun02 sec05 sub03 cat02:
      *     - find snps passing min quality scores
      *   o fun02 sec05 sub03 cat03:
      *     - check if SNP (passed) is a transversion
      *   o fun02 sec05 sub03 cat04:
      *     - move to next SNP/match
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun02 Sec05 Sub03 Cat01:
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
         /*Case: known snp or could be snp*/
            tmpSI =
               min_genMath(
                  refValSI,
                  qryValSI
               );

            *overlapUI += tmpSI;

            while(tmpSI > 0)
            { /*Loop: count number of SNPs*/
               if(
                     refSTPtr->seqStr[uiRef]
                  == qrySTPtr->seqStr[uiQry]
               ) goto nextNt_fun02_sec05_sub03_cat02;

               if(qryQBl)
               { /*If: have query q-score entry*/
                  qUC =
                       qrySTPtr->qStr[uiQry]
                     - def_adjQ_samEntry;

                  if(qUC < minQUC)
                     goto nextNt_fun02_sec05_sub03_cat02;
               } /*If: have query q-score entry*/

               if(refQBl)
               { /*If: have reference q-score entry*/
                  qUC =
                       refSTPtr->qStr[uiRef]
                     - def_adjQ_samEntry;

                  if(qUC < minQUC)
                     goto nextNt_fun02_sec05_sub03_cat02;
               } /*If: have reference q-score entry*/

               ++distSL;

               nextNt_fun02_sec05_sub03_cat02:;

               ++uiQry;
               ++uiRef;

               --refValSI;
               --qryValSI;

               --tmpSI;
            } /*Loop: count number of SNPs*/

            break;
         /*Case: known snp or could be snp*/

         /***********************************************\
         * Fun02 Sec05 Sub04:
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
               *overlapUI += tmpSI;
            } /*If: not deletion event*/

            break;
         /*Case: treat as same*/
 
         /***********************************************\
         * Fun02 Sec05 Sub05:
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

            if(tmpSI >= (sint) indelLenUI)
            { /*If: keeping deletion event*/
               distSL += tmpSI;
               ++(*indelEventsUI);
               *numIndelUI += tmpSI;
            } /*If: keeping deletion event*/

            if(refSTPtr->cigTypeStr[siRefCig] == 'D')
               uiQry += tmpSI;

            else
            { /*Else: query was a deletion*/
               uiRef += tmpSI;
               *overlapUI += tmpSI;
            } /*Else: query was a deletion*/

            qryValSI -= tmpSI;
            refValSI -= tmpSI;

            break;
         /*Case: deletion (only one sequence)*/

         /***********************************************\
         * Fun02 Sec05 Sub06:
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

            if(tmpSI >= (sint) indelLenUI)
            { /*If: keeping query deletion*/
               distSL += tmpSI;
               ++(*indelEventsUI);
               *numIndelUI += tmpSI;
            } /*If: keeping query deletion*/

            if(qrySTPtr->cigTypeStr[siRefCig] == 'D')
               *overlapUI += tmpSI;

            break;
         /*Case: insertion (only one sequence)*/
      } /*Switch: find mutation combination*/

      /**************************************************\
      * Fun02 Sec05 Sub07:
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
   ^ Fun02 Sec06:
   ^   - return the edit distance
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun02_sec06;

   noOverlap_fun02_sec06:;
   distSL = def_noOverlap_edDist;
   goto ret_fun02_sec06;

   diffRef_fun02_sec06:;
   distSL = def_diffRef_edDist;
   goto ret_fun02_sec06;

   noSeq_fun02_sec06:;
   distSL = def_noSeq_edDist;
   goto ret_fun02_sec06;

   noMap_fun02_sec06:;
   distSL = def_noMap_edDist;
   goto ret_fun02_sec06;

   ret_fun02_sec06:;
   return distSL;
} /*readCmpDist_edDist*/

/*-------------------------------------------------------\
| Fun03: dist_edDist
|   - gets edit distances for mapped reference
|   - deletions and insertions are only counted if they
|     execed a minimum length.
| Input:
|   - qrySTPtr:
|     o pointer to samEntry structure with read (query) to
|       find the edit distance for
|   - indelLenUI:
|     o minimum length for a indel to count as an event
|   - minQUC:
|     o minimum q-score to not discard an snp
|   - numIndelUI:
|     o pointer to unisigned int to hold the number of
|       indels in edit distance
|   - indelEventsUI:
|     o pointer to unsigned int to hold the number of
|       indel events (times a group of indels occured)
| Output:
|   - Modifies:
|     o numIndelUI to have number indels kept
|     o indelEventsUI to have number of indel events
|     o overlapUI to hold number of overlapped bases
|   - Returns:
|     o edit distance for query and mapped reference
|     o negative value if could not find distance
|       - def_noSeq_edDist if one read is missing sequence
|       - def_noMap_edDist if one read is unmapped
\-------------------------------------------------------*/
signed long
dist_edDist(
   struct samEntry *qrySTPtr, /*read for edit distance*/
   unsigned int indelLenUI,   /*min indel length*/
   unsigned char minQUC,      /*min Q-score (snp/ins)*/
   unsigned int *numIndelUI,  /*number indels kept*/
   unsigned int *indelEventsUI /*nubmer indel events*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - gets edit distances between the query & reference
   '   o fun03 sec01:
   '     - variable declerations
   '   o fun03 sec02:
   '     - check if can get edit distance
   '   o fun03 sec03:
   '     - initialize and check if have q-score entreis
   '   o fun03 sec04:
   '     - find edit distance
   '   o fun03 sec06:
   '     - return the edit distance
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint siQryCig = 0;  /*query cigar position*/
   sint qryValSI = 0;
   uint uiQry = 0;     /*query nucleotide on*/

   slong distSL = 0;    /*edit distance*/

   schar qryQBl = 0;
   uchar qUC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - check if can get edit distance
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *indelEventsUI = 0;
   *numIndelUI = 0;

   if(qrySTPtr->flagUS & 4)
      goto noMap_fun03_sec06;

   if(
         ! qrySTPtr->seqStr
      || *qrySTPtr->seqStr == '*' 
      || *qrySTPtr->seqStr == '\0'
   ) goto noSeq_fun03_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
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
   ^ Fun03 Sec04:
   ^   - find edit distance
   ^   o fun03 sec04 sub01:
   ^     - get past masking at start
   ^   o fun03 sec04 sub02:
   ^     - start loop and check if at end (hard/soft mask)
   ^   o fun03 sec04 note01:
   ^     - table of cigar entries added together
   ^   o fun03 sec04 sub03:
   ^     - start switch/check match/possible snp combos
   ^   o fun03 sec04 sub04:
   ^     - check matches or same error type
   ^   o fun03 sec04 sub05:
   ^     - handel deletion cases
   ^   o fun03 sec04 sub06:
   ^     - handel insertion cases
   ^   o fun03 sec04 sub07:
   ^     - check if move to next cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec04 Sub01:
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
   * Fun03 Sec04 Sub02:
   *   - start loop and check if at end (hard/soft mask)
   \*****************************************************/

   while(siQryCig < (sint) qrySTPtr->lenCigUI)
   { /*Loop: get edit distance*/
      if(qrySTPtr->cigTypeStr[siQryCig] == 'S')
         break; /*soft masking only at ends (finished)*/

      if(qrySTPtr->cigTypeStr[siQryCig] == 'H')
         break; /*hard masking only at ends (finished)*/

      /**************************************************\
      * Fun03 Sec04 Sub03:
      *   - start switch/check match/possible snp combos
      *   o fun03 sec04 sub03 cat01:
      *     - find number shared snps/matchs + check case
      *   o fun03 sec04 sub03 cat02:
      *     - find snps passing min quality scores
      *   o fun03 sec04 sub03 cat03:
      *     - check if SNP (passed) is a transversion
      *   o fun03 sec04 sub03 cat04:
      *     - move to next SNP/match
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec04 Sub03 Cat01:
      +   - find number shared snps/matchs + check case
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      switch(qrySTPtr->cigTypeStr[siQryCig])
      { /*Switch: find mutation combination*/
         /*cases were I have to compare sequences*/
         case 'M':
         case '=':
         /*Case: likley match (or no way to check)*/
            uiQry += qryValSI;
            break;
         /*Case: likley match (or no way to check)*/

         case 'X':
         /*Case: snp*/
            while(qryValSI > 0)
            { /*Loop: count number of SNPs*/
               if(qryQBl)
               { /*If: have query q-score entry*/
                  qUC =
                       qrySTPtr->qStr[uiQry]
                     - def_adjQ_samEntry;

                  if(qUC >= minQUC)
                     ++distSL;
               } /*If: have query q-score entry*/

               else
                  ++distSL;

               ++uiQry;
               --qryValSI;
            } /*Loop: count number of SNPs*/

            break;
         /*Case: snp*/

         /***********************************************\
         * Fun03 Sec04 Sub05:
         *   - handel deletion cases
         \***********************************************/

         /*deletion cases*/
         case 'D':
         /*Case: deletion*/
            if(qryValSI >= (sint) indelLenUI)
            { /*If: keeping deletion event*/
               distSL += qryValSI;
               ++(*indelEventsUI);
               *numIndelUI += qryValSI;
            } /*If: keeping deletion event*/

            break;
         /*Case: deletion*/

         /***********************************************\
         * Fun03 Sec04 Sub06:
         *   - handel insertion cases
         \***********************************************/

         case 'I':
         /*Case: insertion*/
            uiQry += qryValSI;

            if(qryValSI >= (sint) indelLenUI)
            { /*If: keeping query deletion*/
               distSL += qryValSI;
               ++(*indelEventsUI);
               *numIndelUI += qryValSI;
            } /*If: keeping query deletion*/

            break;
         /*Case: insertion*/
      } /*Switch: find mutation combination*/

      /**************************************************\
      * Fun03 Sec04 Sub07:
      *   - check if move to next cigar entry
      \**************************************************/

      ++siQryCig;
      qryValSI = qrySTPtr->cigArySI[siQryCig];
   } /*Loop: get edit distance*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec06:
   ^   - return the edit distance
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun03_sec06;

   noSeq_fun03_sec06:;
   distSL = def_noSeq_edDist;
   goto ret_fun03_sec06;

   noMap_fun03_sec06:;
   distSL = def_noMap_edDist;
   goto ret_fun03_sec06;

   ret_fun03_sec06:;
   return distSL;
} /*dist_edDist*/

/*-------------------------------------------------------\
| Fun04: percDist_edDist
|   - gives a rough idea on precentage of difference from
|     error
|   - not great, but allows lumping reads together
| Input:
|   - distSL:
|     o edit distance to find probablity for
|   - overlapUI:
|     o length of alignment (overlap between query & ref)
|   - percErrF:
|     o expected percent of errors (0 to 1) in reads
|     o if doing read to read comparsions; double this
| Output:
|   - Returns:
|     o probablity of edit distance due to chance
\-------------------------------------------------------*/
signed int
percDist_edDist(
   signed long distSL,
   unsigned int overlapUI,
   float percErrF
){
   return 100 * (distSL / (overlapUI * percErrF));
} /*percDist_edDist*/

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
