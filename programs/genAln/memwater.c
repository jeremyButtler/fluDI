/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' memwater SOF: Start Of Functions
'   - has memory efficent Waterman Smith alignment
'   o header:
'     - Included libraries
'   o fun01 memwater:
'     - Run a memory efficent Waterman Smith alignment on
'       input sequences
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "memwater.h"

#include "../genBio/seqST.h"

#include "alnSet.h"
#include "indexToCoord.h"

/*.h files only*/
#include "../genLib/dataTypeShortHand.h"
#include "../genLib/genMath.h" /*only using .h portion*/
#include "alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden files
!   o std #include <stdio.h>
!   o .c  #include "../genLib/base10str.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: memwater
|   - Performs a memory efficent Smith Waterman alignment
|     on a pair of sequences
| Input;
|   - qrySeqSTVoidPtr:
|     o Point to an seqST with the query sequence and
|       index 0 coordinates to start (offsetUL)/end
|       (endAlnUL) the alignment.
|   - refSeqSTVoidPtr:
|     o Point to an seqST with the reference sequence
|       and index 0 coordinates to start (offsetUL)/end
|       (endAlnUL) the alignment.
|   - refStartUL:
|     o Pointer to unsigned long to hold the frist
|       reference base in the alignment
|   - refEndUL:
|     o Pointer to unsigned long to hold the last
|       reference base in the alignment
|   - qryStartUL:
|     o Pointer to unsigned long to hold the frist query 
|       base in the alignment
|   - qryEndUL:
|     o Pointer to unsigned long to hold the last query
|       base in the alignment
|   - alnSetVoidPtr:
|     o Pointer to an alnSet structure with the gap open,
|       gap extend, and scoring matrix for the alingment
| Output:
|  - Modifies:
|    o refStartUL to have 1st reference base in alignment
|    o refEndUL to have last reference base in alignment
|    o qryStartUL to have first query base in alignment
|    o qryEndUL to have last query base in alignment
|  - Returns:
|    o score for aligment
|    o negative number for memory errors
\-------------------------------------------------------*/
signed long
memwater(
   struct seqST *qrySTPtr, /*query sequence and data*/
   struct seqST *refSTPtr, /*ref sequence and data*/
   unsigned long *refStartUL,
   unsigned long *refEndUL,
   unsigned long *qryStartUL,
   unsigned long *qryEndUL,
   struct alnSet *settings
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC: memwaterAln
   '  - Run a memory efficent Waterman Smith alignment on
   '    input sequences
   '  o fun01 sec01:
   '    - Variable declerations
   '  o fun01 sec02:
   '    - Allocate memory for alignment
   '  o fun01 sec03:
   '    - Fill in initial negatives for ref
   '  o fun0 sec04:
   '    - Fill the matrix with scores
   '  o fun01 sec05:
   '    - Set up for returing matrix (clean up/wrap up)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01: Variable declerations
   ^  o fun01 sec01 sub01:
   ^    - Variables dealing with the query and reference
   ^      starting positions
   ^  o fun01 sec01 sub02:
   ^    - Variables holding the scores (only two rows)
   ^  o fun01 sec01 sub03:
   ^    - Directinol matrix variables
   ^  o fun01 sec01 sub04:
   ^    - Variables for building returend alignment array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec01 Sub01:
   *  - Variables dealing with the query and reference
   *    starting positions
   \*****************************************************/

   slong scoreSL = 0;     /*score to return*/
   ulong bestStartUL = 0; /*records best starting index*/
   ulong bestEndUL = 0;   /*records best ending index*/

   /*Get start & end of query and reference sequences*/
   schar *refSeqStr = 0;
   schar *qrySeqStr = 0;

   ulong lenRefUL =
      refSTPtr->endAlnUL - refSTPtr->offsetUL + 1;

   ulong lenQryUL =
      qrySTPtr->endAlnUL - qrySTPtr->offsetUL + 1;

   /*Iterators for loops*/
   ulong ulRef = 0;
   ulong ulQry = 0;

   /*****************************************************\
   * Fun01 Sec01 Sub02:
   *  - Variables holding the scores (only two rows)
   \*****************************************************/

   /*slong keepSL = 0; here for unbranched version*/

   slong snpScoreSL = 0;
   slong insScoreSL = 0;
   slong delScoreSL = 0;   /*Score for doing an deletion*/
   slong nextSnpScoreSL = 0;/*Score for next match/snp*/
   slong *scoreHeapArySL = 0; /*matrix to use in alignment*/

   /*****************************************************\
   * Fun01 Sec01 Sub03:
   *  - Directional matrix variables
   \*****************************************************/

   /*Direction matrix (1 cell holds a single direction)*/
   schar *dirRowHeapSC = 0;  /*Holds directions*/

   /*Keeping track of alignment starting positions*/
   ulong indexUL = 0;      /*Index I am at in the matrix*/
   ulong snpIndexUL = 0;   /*last snp index*/
   ulong tmpIndexUL = 0;   /*for getting snp index*/

   ulong *indexHeapUL=0;   /*Row of starting indexes*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^  - Allocate memory for alignment
   ^  o fun01 sec02 sub01:
   ^    - Allocate memory for the alignment
   ^  o fun01 sec02 sub02:
   ^    - Allocate memory for keeping track of indexes
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *  - Allocate memory for the alignment
   \****************************************************/

   dirRowHeapSC =
      calloc(
         (lenRefUL + 1),
         sizeof(schar)
      );

   if(dirRowHeapSC == 0)
      goto memErr_fun01_sec05_sub03;

   scoreHeapArySL =
      calloc(
         (lenRefUL + 1),
          sizeof(slong)
      );
   /*+ 1 is for the indel column*/

   if(scoreHeapArySL == 0)
      goto memErr_fun01_sec05_sub03;

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *  - Get memory for keeping track of starting indexes
   \*****************************************************/

   /*Set up the first row of starting indexes*/
   indexHeapUL =
      calloc(
         (lenRefUL + 1),
         sizeof(ulong)
      );

   if(indexHeapUL == 0)
      goto memErr_fun01_sec05_sub03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^  - Fill in initial negatives for reference
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      indexUL = 0;
      indexUL <= lenRefUL;
      ++indexUL
   ){ /*loop; till have initalized the first row*/
      dirRowHeapSC[indexUL] = def_mvStop_alnDefs;
      indexHeapUL[indexUL] = indexUL + 1;
      scoreHeapArySL[indexUL] = 0;
   } /*loop; till have initalized the first row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^  - Fill the matrix with scores
   ^  o fun01 sec04 sub01:
   ^    - Final set up before scoring the matrix
   ^  o fun01 sec04 sub02:
   ^    - get snp and ins scores + start loop
   ^  o fun01 sec04 sub03:
   ^    - find high score
   ^  o fun01 sec04 sub04:
   ^    - check if keep score (score > 0)
   ^  o fun01 sec04 sub05:
   ^    - find next deletion score
   ^  o fun01 sec04 sub06:
   ^    - check if have new high score
   ^  o fun01 sec04 sub07:
   ^    - prepare to score the next row in the matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec04 Sub01:
   *  - Final set up before scoring the matrix
   \*****************************************************/

   /*Move the row of starting indexes to the last row*/
   indexHeapUL[0] = indexUL;

   /*These are always negative*/
   delScoreSL = 0;
   nextSnpScoreSL = 0;
   snpIndexUL = indexHeapUL[0];

   /*incurment to frist base*/
   ++indexUL;
   qrySeqStr = qrySTPtr->seqStr + qrySTPtr->offsetUL;
   refSeqStr = refSTPtr->seqStr + refSTPtr->offsetUL - 1;
      /*offseting reference by 1 to account for the gap
      `  column
      */

   /*****************************************************\
   * Fun01 Sec04 Sub02:
   *  - get snp and ins scores + start loop
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
      ){ /*loop; compare one query to one reference base*/

         snpScoreSL =
            getScore_alnSet(
               qrySeqStr[ulQry],
               refSeqStr[ulRef],
               settings
            ); /*find the score for the base pairs*/

         snpScoreSL += nextSnpScoreSL;
         nextSnpScoreSL = scoreHeapArySL[ulRef];

         /*insertion score*/
         #ifdef NOEXTEND
            insScoreSL =
               scoreHeapArySL[ulRef] + settings->gapSS;
         #else
            insScoreSL = scoreHeapArySL[ulRef];
            insScoreSL +=
               settings->insArySS[dirRowHeapSC[ulRef]];
         #endif

         /***********************************************\
         * Fun02 Sec04 Sub03:
         *   - find high score
         \***********************************************/

         scoreHeapArySL[ulRef] =
            max_genMath(
               insScoreSL,
               snpScoreSL
         ); /*find if ins/snp is best (5 Op)*/

         tmpIndexUL = indexHeapUL[ulRef];

         indexHeapUL[ulRef] =
            ifmax_genMath(
               insScoreSL,
               snpScoreSL,
               indexHeapUL[ulRef], /*insertion index*/
               snpIndexUL          /*snp index*/
            ); /*get index of high score*/

         snpIndexUL = tmpIndexUL;

         /*find direction (4 Op)*/
         dirRowHeapSC[ulRef] =
            scoreHeapArySL[ulRef] > delScoreSL;

         dirRowHeapSC[ulRef] <<=
            (snpScoreSL < insScoreSL);

         ++dirRowHeapSC[ulRef];

         /*finish finding max's*/
         indexHeapUL[ulRef] =
            ifmax_genMath(
               delScoreSL,
               scoreHeapArySL[ulRef],
               indexHeapUL[ulRef - 1],    /*del index*/
               indexHeapUL[ulRef]         /*current best*/
            ); /*get index of high score*/

         scoreHeapArySL[ulRef] =
            max_genMath(
               delScoreSL,
               scoreHeapArySL[ulRef]
         ); /*find if del is best (5 Op)*/
            
         /***********************************************\
         * Fun02 Sec04 Sub04:
         *   - check if keep score (score > 0)
         \***********************************************/

         if(scoreHeapArySL[ulRef] <= 0)
         {
            dirRowHeapSC[ulRef] = 0;
            scoreHeapArySL[ulRef] = 0;
            indexHeapUL[ulRef] = indexUL + 1;
            /*always one off for new index's*/
         }

         /* branchless method is slower here
         keepSL = (slong) -(scoreHeapArySL[ulRef] > 0);
         dirRowHeapSC[ulRef] &= keepSL;
         scoreHeapArySL[ulRef] &= keepSL;

         indexHeapUL[ulRef] &= keepSL;
         indexHeapUL[ulRef] |= ( indexUL & (~keepSL) );
         */

         /***********************************************\
         * Fun02 Sec04 Sub05:
         *   - find next deletion score
         \***********************************************/

         #ifdef NOEXTEND
            delScoreSL =
               scoreHeapArySL[ulRef] + settings->gapSS;
         #else
            delScoreSL = scoreHeapArySL[ulRef];
            delScoreSL +=
               settings->delArySS[dirRowHeapSC[ulRef]];
         #endif

         /***********************************************\
         * Fun02 Sec04 Sub06:
         *   - check if have new high score
         \***********************************************/

         if(scoreSL < scoreHeapArySL[ulRef])
         { /*If: this is the best score*/
            scoreSL = scoreHeapArySL[ulRef];
            bestStartUL = indexHeapUL[ulRef];
            bestEndUL = indexUL;
         } /*If: this was an snp or match*/

         ++indexUL;
      } /*loop; compare one query to one reference base*/

     /***************************************************\
     *  Fun01 Sec04 Sub07:
     *   - prepare for the next round
     \***************************************************/

     /*Get scores set up for the gap column*/
	  nextSnpScoreSL = 0;
     delScoreSL = 0;

     indexHeapUL[0] = indexUL; /*next index*/
     snpIndexUL = indexHeapUL[0];

     ++indexUL; /*Set index for the next base pair*/
   } /*loop; compare query base against all ref bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec05:
   ^  - Set up for returing the matrix (clean up/wrap up)
   ^  o fun01 sec05 sub01:
   ^    - get coordinates
   ^  o fun01 sec05 sub02:
   ^    - no error clean up
   ^  o fun01 sec05 sub03:
   ^    - memory error clean up
   ^  o fun01 sec05 sub04:
   ^    - general clean up (no error or error)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec05 Sub01:
   *  - get coordinates
   \*****************************************************/

   *refStartUL =
      refCoord_indexToCoord(
         lenRefUL,
         bestStartUL
      ); /*Convert the starting index to coordinates*/

   *qryStartUL =
      qryCoord_indexToCoord(
         lenRefUL,
         bestStartUL
      ); /*Convert the starting index to coordinates*/

   *refStartUL += refSTPtr->offsetUL;
   *qryStartUL += qrySTPtr->offsetUL;

   *refEndUL =
      refCoord_indexToCoord(
         lenRefUL,
         bestEndUL
      ); /*Convert ending index to coordinates*/

   *qryEndUL =
      qryCoord_indexToCoord(
         lenRefUL,
         bestEndUL
      ); /*Convert ending index to coordinates*/

   *refEndUL += refSTPtr->offsetUL;
   *qryEndUL += qrySTPtr->offsetUL;

   /*****************************************************\
   * Fun01 Sec05 Sub02:
   *  - no error clean up
   \*****************************************************/

   goto cleanUp_fun01_sec05_sub04;

   /*****************************************************\
   * Fun01 Sec05 Sub03:
   *  - memory error clean up
   \*****************************************************/

   memErr_fun01_sec05_sub03:;
   scoreSL = -1;
   goto cleanUp_fun01_sec05_sub04;

   /*****************************************************\
   * Fun01 Sec05 Sub04:
   *  - general clean up (no error or error)
   \*****************************************************/

   cleanUp_fun01_sec05_sub04:;

   free(dirRowHeapSC);
   dirRowHeapSC = 0;

   free(scoreHeapArySL);
   scoreHeapArySL = 0;

   free(indexHeapUL);
   indexHeapUL = 0;

   return scoreSL;
} /*memwater*/

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
