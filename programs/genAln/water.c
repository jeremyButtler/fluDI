/*########################################################
# Name water
# Use:
#  o Holds functions to do a Smith Waterman pairwise
#    alignment
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'   o header:
'     - included libraries
'   o .c fun01: scoreIndel_water
'     - gets the indel score for a water alignment
'   o fun01 water:
'     - run a memory efficent Waterman Smith alignment on
'       input sequences
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

#include "water.h"

#include "../genLib/seqST.h"

#include "alnSet.h"
#include "dirMatrix.h"
#include "indexToCoord.h"

/*.h files only*/
#include "../genLib/dataTypeShortHand.h"
#include "../genLib/genMath.h" /*only using .h commands*/
#include "alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden files
!   o std #include <stdio.h>
!   o .c  #include "../genLib/base10str.h"
!   o .c  #include "../genLib/ulCp.h"
!   o .c  #include "../genLib/charCp.h"
!   o .c  #include "../genLib/numToStr.h"
!   o .c  #include "../genLib/samEntry.h"
!   o .h  #include "../genLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: water
|   - run a Waterman Smith alignment on input sequences
| Input:
|   - qrySTPtr:
|     o pointer to seqST with query sequence 
|       - qrySTPtr->offsetUL; first query base to align
|       - qrySTPtr->endAlnUL; last query base to align
|   - refSTPtr:
|     o pointer to seqST with reference sequence 
|       - refSTPtr->offsetUL; 1st reference base to align
|       - refSTPtr->endAlnUL; last reference base to align
|   - matrixSTPtr:
|     o pointer to dirMatrix to use for the alingment
|   - alnSet:
|     o pointer to alnSet with alignment settings
| Output:
|  - Modifies:
|    o errSC in matrixSTPtr to def_memErr_water for
|      memory errors
|    o allocates memory for dirMatrixSC and scoreAryUL
|      if they are to small
|    o updates lenMatrixUL and lenScoreUL if dirMatrixSC
|      or scoreAryUL are resized
|  - Returns:
|    o score for alignment
\-------------------------------------------------------*/
signed long
water(
   struct seqST *qrySTPtr, /*query sequence and data*/
   struct seqST *refSTPtr, /*ref sequence and data*/
   struct dirMatrix *matrixSTPtr, /*direction matrix*/
   struct alnSet *settings     /*settings for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC: WatermanAln
   '  - Run a Waterman Smith alignment on input sequences
   '  o fun01 sec01:
   '    - variable declerations
   '  o fun01 sec02:
   '    - allocate memory for alignment
   '  o fun01 sec03:
   '    - fill in initial negatives for ref
   '  o fun01 sec04:
   '    - fill the matrix with scores
   '  o fun01 sec05:
   '    - set up for returing matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01: Variable declerations
   ^  o fun01 sec01 sub01:
   ^    - variables dealing with the query and reference
   ^      starting positions
   ^  o fun01 sec01 sub02:
   ^    - variables holding the scores (only two rows)
   ^  o fun01 sec01 sub03:
   ^    - directinol matrix variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec01 Sub01:
   *  - variables dealing with the query and reference
   *    starting positions
   \*****************************************************/

   /*Get start & end of query and reference sequences*/
   schar *refSeqStr =
      refSTPtr->seqStr + refSTPtr->offsetUL;

   schar *qrySeqStr =
      qrySTPtr->seqStr + qrySTPtr->offsetUL;

   /*Find the length of the reference and query*/
   ulong lenQryUL =
      qrySTPtr->endAlnUL - qrySTPtr->offsetUL + 1;

   ulong lenRefUL =
      refSTPtr->endAlnUL - refSTPtr->offsetUL + 1;
     /*The + 1 is to account for index 0 of endAlnUL*/

   ulong lenMatrixUL = (lenRefUL + 1) * (lenQryUL + 1);
     /*+1 for the gap column and row*/

   ulong ulRef = 0;
   ulong ulQry = 0;

   /*Set up counters for the query and reference base
   `  index
   */
   /*****************************************************\
   * Fun01 Sec01 Sub02:
   *  - variables holding the scores (only two rows)
   \*****************************************************/

   slong snpScoreSL = 0;    /*Score for deletion*/
   slong nextSnpScoreSL = 0;/*Score for match/snp*/

   slong insScoreSL = 0;    /*Score for deletion*/
   slong delScoreSL = 0;    /*Score for deletion*/

   /*Marks when to reset score buffer (every second row)*/
   slong *scoreArySL = 0;/*scoring row for alignment*/

   /*****************************************************\
   * Fun01 Sec01 Sub03:
   *  - directinol matrix variables
   \*****************************************************/

   /*Direction matrix (one cell holds one direction)*/
   schar *dirMatrixSC = 0;/*Direction matrix*/
   schar *insDir = 0;    /*Direction above cell*/
   ulong indexUL = 0;
   /*slong keepSL = 0; here for unbranched version*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - allocate memory for alignment
   ^   o fun01 sec02 sub01:
   ^     - set up the directional matrix
   ^   o fun01 sec02 sub02:
   ^     - set up score array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *   - set up the directional matrix
   \*****************************************************/

   if(matrixSTPtr->lenMatrixUL < lenMatrixUL)
   { /*If: need to resize the matrix*/
      free(matrixSTPtr->dirMatrixSC);
      matrixSTPtr->dirMatrixSC = 0;
      matrixSTPtr->lenMatrixUL = 0;

      matrixSTPtr->dirMatrixSC =
         malloc((lenMatrixUL + 1) * sizeof(char));

      if(matrixSTPtr->dirMatrixSC == 0)
         goto memErr_fun01_sec05;

       matrixSTPtr->lenMatrixUL = lenMatrixUL;
   } /*If: need to resize the matrix*/

   dirMatrixSC = matrixSTPtr->dirMatrixSC;

   blank_dirMatrix(matrixSTPtr);

   matrixSTPtr->lenRefUL = lenRefUL;
   matrixSTPtr->refOffsetUL = refSTPtr->offsetUL;
   matrixSTPtr->refEndUL = refSTPtr->endAlnUL;

   matrixSTPtr->lenQryUL = lenQryUL;
   matrixSTPtr->qryOffsetUL = qrySTPtr->offsetUL;
   matrixSTPtr->qryEndUL = qrySTPtr->endAlnUL;

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *   - set up score array
   \*****************************************************/

   if(matrixSTPtr->lenScoreUL < lenRefUL + 1)
   { /*If: need to make a larger score array*/
      free(matrixSTPtr->scoreArySL);
      matrixSTPtr->scoreArySL = 0;
      matrixSTPtr->lenScoreUL = 0;
      
      matrixSTPtr->scoreArySL =
         calloc((lenRefUL + 1), sizeof(long));

      if(! matrixSTPtr->scoreArySL)
         goto memErr_fun01_sec05;

      matrixSTPtr->lenScoreUL = lenRefUL;
   } /*If: need to make a larger score array*/

   scoreArySL = matrixSTPtr->scoreArySL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^  - fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      indexUL = 0;
      indexUL <= lenRefUL;
      ++indexUL
   ){ /*Loop: initialize the first row*/
      dirMatrixSC[indexUL] = def_mvStop_alnDefs;
      scoreArySL[indexUL] = 0;
   } /*Loop: initialize the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^  - fill the matrix with scores
   ^  o fun01 sec04 sub01:
   ^    - get initial scores
   ^  o fun01 sec04 sub02:
   ^    - start loops
   ^  o fun01 sec04 sub03:
   ^    - find bests score for each base pair
   ^  o fun01 sec04 sub04:
   ^    - find high score
   ^  o fun01 sec04 sub05:
   ^    - check if keep score (score > 0)
   ^  o fun01 sec04 sub06:
   ^    - find next deletion score
   ^  o fun01 sec04 sub07:
   ^    - check if have new high score
   ^  o fun01 sec04 sub08:
   ^    - set up for scoring next row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec04 Sub01:
   *  - get initial scores
   \*****************************************************/

   /*set up pointers*/
   insDir = dirMatrixSC;

   /*set up scores*/
   delScoreSL = 0;
   nextSnpScoreSL = 0;
   dirMatrixSC[indexUL] = def_mvStop_alnDefs;

   /*move to first base*/
   ++indexUL;
   refSeqStr = refSTPtr->seqStr + refSTPtr->offsetUL - 1;
   qrySeqStr = qrySTPtr->seqStr + qrySTPtr->offsetUL;

   /*****************************************************\
   * Fun01 Sec04 Sub02:
   *   - start loops
   \*****************************************************/

   /*Starting on the first sequence row*/
   for(
      ulQry = 0;
      ulQry < lenQryUL;
      ++ulQry
   ){ /*loop; compare query base against all ref bases*/

      for(
         ulRef = 1;
         ulRef <= lenRefUL;
         ++ulRef
      ){ /*Loop; compare one query to one reference base*/

         /***********************************************\
         * Fun01 Sec04 Sub03:
         *   - find bests score for each base pair
         \***********************************************/

         snpScoreSL =
            getScore_alnSet(
               qrySeqStr[ulQry],
               refSeqStr[ulRef],
               settings
            ); /*Find the score for the base pairs*/

         snpScoreSL += nextSnpScoreSL;
         nextSnpScoreSL = scoreArySL[ulRef];

         #ifdef NOEXTEND
            insScoreSL =
               scoreArySL[ulRef] + settings->gapSS;
         #else
            insScoreSL = scoreArySL[ulRef];
            insScoreSL +=
               settings->insArySS[insDir[ulRef]];
         #endif

         /***********************************************\
         * Fun01 Sec04 Sub04:
         *   - find high score
         \***********************************************/

         scoreArySL[ulRef] =
            max_genMath(
               insScoreSL,
               snpScoreSL
         ); /*find if ins/snp is best (5 Op)*/
 
         /*find direction (4 Op)*/
         dirMatrixSC[indexUL] =
            scoreArySL[ulRef] > delScoreSL;

         dirMatrixSC[indexUL] <<=
            (snpScoreSL < insScoreSL);

         ++dirMatrixSC[indexUL];

         /*Logic:
         `   - noDel: maxSC > delSc:
         `     o 1 if deletion not max score
         `     o 0 if deletion is max score
         `   - type: noDel << (snpSc < insSc):
         `     o 1 << 1 = 2 if insertion is maximum
         `     o 1 << 0 = 1 if snp is maximum
         `     o 0 << 0 = 0 if del is max, and snp > ins
         `     o 0 << 1 = 0 if del is max, but ins >= snp
         `   - dir: type + 1
         `     o adds 1 to change from stop to direction
         */

         /*finish max*/
         scoreArySL[ulRef] =
            max_genMath(
               delScoreSL,
               scoreArySL[ulRef]
         ); /*find if del is best (5 Op)*/

         /***********************************************\
         * Fun01 Sec04 Sub05:
         *   - check if keep score (score > 0)
         \***********************************************/

         if(scoreArySL[ulRef] <= 0)
         {
            dirMatrixSC[indexUL] = 0;
            scoreArySL[ulRef] = 0;
         }

         /*check if have negative or positive score*/
         /*no branched version is slower here
         keepSL = (slong) -(scoreArySL[ulRef] > 0);
         dirMatrixSC[indexUL] &= keepSL;
         scoreArySL[ulRef] &= keepSL;
         */

         /***********************************************\
         * Fun01 Sec04 Sub06:
         *   - find next deletion score
         \***********************************************/

         #ifdef NOEXTEND
            delScoreSL =
               scoreArySL[ulRef] + settings->gapSS;
         #else
            delScoreSL = scoreArySL[ulRef];
            delScoreSL +=
               settings->delArySS[dirMatrixSC[indexUL]];
         #endif

         /***********************************************\
         * Fun01 Sec04 Sub07:
         *   - check if have new high score
         \***********************************************/

         /*This is faster than the branchless option.
         ` I am guessing to much is done in the if and 
         ` that the if if fired rarely.
         */
         if(
             matrixSTPtr->scoreSL
           < scoreArySL[ulRef]
         ){ /*if have a new best score*/
            matrixSTPtr->scoreSL = scoreArySL[ulRef];
            matrixSTPtr->indexUL = indexUL;
         } /*if have a new best score*/

         ++indexUL;
      } /*Loop; compare one query to one reference base*/

      /**************************************************\
      *  Fun01 Sec04 Sub08:
      *    - set up for scoring next row
      \**************************************************/

      delScoreSL = 0;
      nextSnpScoreSL = 0; /*indel column is always 0*/

      insDir += ulRef;
      dirMatrixSC[indexUL] = def_mvStop_alnDefs;

      ++indexUL;
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec05:
   ^  - set up for returing the matrix (clean up/wrap up)
   ^  o fun01 sec05 sub01:
   ^    - no error clean up and return
   ^  o fun01 sec05 sub02:
   ^    - error clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec05 Sub01:
   *   - no error clean up and return
   \*****************************************************/

   /*Move back to the lower right conor cell
   ` This is not needed, but is nice.
   */

   --indexUL; /*get off last -1*/
   dirMatrixSC[indexUL] = def_mvStop_alnDefs;
   return matrixSTPtr->scoreSL; /*best score*/

   /*****************************************************\
   * Fun01 Sec05 Sub02:
   *   - error clean up and return
   \*****************************************************/

   memErr_fun01_sec05:;
   matrixSTPtr->errSC = def_memErr_water;
   goto errCleanUp_fun01_sec05;

   errCleanUp_fun01_sec05:;

   return 0;
} /*water*/

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
