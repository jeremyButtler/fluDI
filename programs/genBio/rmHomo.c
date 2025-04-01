/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' rmHomo SOF: Start Of File
'   - removes indels in homopolymers
'   o header:
'     - included libraries
'   o fun01: indel_rmHomo
'     - remove indels from homopolymers
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

#include "rmHomo.h"

#include "../genLib/ulCp.h"
#include "samEntry.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
!   - std #include <stdio.h>
!   - .c  #include "../genLib/base10str.h"
!   - .c  #include "../genLib/numToStr.h"
!   - .c  #include "../genLib/strAry.h"
!   - .h  #include "../genLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: indel_rmHomo
|   - remove indels from homopolymers
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to sequence to remove
|       indels from homopolymers
|   - refStr:
|     o c-string with reference sequence in samSTPtr is
|       mapped to
|   - minHomoSI:
|     o minimum homopolymer length to remove indel from
|   - maxIndelSI:
|     o maximum indel size to remove
|   - maskSC:
|     o base to mask deletions with
|     o if 0 if(! maskSC), then uses reference base
|   - seqBuffStrPtr:
|     o c-string pointer to buffer to use in copy sequence
|     o size to avoid reallocation of buffer
|       * sequence length + number deletions + 8
|   - sizeSeqUIPtr:
|     o unsigned long pointer to with size of
|       seqBuffStrPtr
|   - qBuffStrPtr:
|     o c-string pointer to buffer to use in copy q-score
|       entry
|     o size to avoid reallocation of buffer
|       * sequence length + number deletions + 8
|   - sizeQUIPtr:
|     o unsigned long pointer to with size of qBuffStrPtr
| Output:
|   - Modifies:
|     o seqStr in samSTPtr to have no homopolymer indels
|     o seqBuffStrPtr and qBuffStrPtr to be resized (or
|       allocated) to sequence length + deletions
|     o sizeQUIPtr and sizeSeqUIPtr to have their buffer
|       size if seqBuffStrPtr or qBuffStrPtr are resized
|   - Returns:
|     o 0 for no errors
|     o 1 for memory errors
\-------------------------------------------------------*/
signed char
indel_rmHomo(
   struct samEntry *samSTPtr,/*sequence with indels*/
   signed char *refStr,      /*reference sequence*/
   signed int minHomoSI,     /*min homoplymer length*/
   signed int maxIndelSI,    /*indel to large to remove*/
   signed char maskSC,       /*base to mask dels with*/
   signed char **seqBuffStrPtr, /*sequence buffer*/
   unsigned int *sizeSeqUIPtr,  /*size of seqBuffStrPtr*/
   signed char **qBuffStrPtr,   /*q-score buffer*/
   unsigned int *sizeQUIPtr     /*size of qBuffStrPtr*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - remove indels from homopolymers
   '   o fun01 sec01:
   '     - variable declarations 
   '   o fun01 sec02:
   '     - memory allocation
   '   o fun01 sec03:
   '     - indel removal
   '   o fun01 sec04:
   '     - add null to end of c-strings (sequence/cigar)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned char qUC = 0;

   unsigned int uiCig = 0;
   unsigned int cpCigUI = 0;

   unsigned int seqPosUI = 0;
   unsigned int dupPosUI = 0;

   unsigned int refPosUI = 0;
   unsigned int tmpPosUI = 0;
   unsigned int lenHomoUI = 0;

   unsigned int lenBuffUI = 0;
   signed char qBl = 0; /*1: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - memory allocation
   ^   o fun01 sec02 sub01:
   ^     - get maximum possible length
   ^   o fun01 sec02 sub02:
   ^     - resize sequence buffers
   ^   o fun01 sec02 sub03:
   ^     - resize q-score buffer (if have)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *   - get maximum possible length
   \*****************************************************/

   lenBuffUI = samSTPtr->readLenUI;
   lenBuffUI += samSTPtr->numDelUI;
      /*deletions mean adding in a base, so need to have
      `  some extra room
      */

   lenBuffUI += 8;
      /*prevents read overflows from ulCp functions*/

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *   - resize sequence buffers
   \*****************************************************/

   if(*sizeSeqUIPtr < lenBuffUI)
   { /*If: need to resize buffer*/
      if(*seqBuffStrPtr)
         free(*seqBuffStrPtr);
      *seqBuffStrPtr = 0;
   } /*If: need to resize buffer*/

   if(! *seqBuffStrPtr)
   { /*If: need new sequence buffer*/
      *sizeSeqUIPtr = 0;

      *seqBuffStrPtr =
         malloc(lenBuffUI * sizeof(signed char));
      if(! *seqBuffStrPtr)
         goto memErr_fun01_sec05;

      *sizeSeqUIPtr = lenBuffUI;
   } /*If: need new sequence buffer*/

   /*****************************************************\
   * Fun01 Sec02 Sub03:
   *   - resize q-score buffer (if have)
   \*****************************************************/

   if(! samSTPtr->qStr)
      qBl = 1;
   else if(samSTPtr->qStr[0] == '\0')
      qBl = 1;
   else if(
         samSTPtr->qStr[0] == '*'
      && samSTPtr->qStr[0] == '\0'
   ) ; /*no q-score entry*/

   else
      qBl = 1;

   if(qBl)
   { /*If: have q-score entry*/
      if(*sizeQUIPtr < lenBuffUI)
      { /*If: need to resize buffer*/
         if(*qBuffStrPtr)
            free(*qBuffStrPtr);
         *qBuffStrPtr = 0;
      } /*If: need to resize buffer*/

      if(! *qBuffStrPtr)
      { /*If: need new quence buffer*/
         *sizeQUIPtr = 0;

         *qBuffStrPtr =
            malloc(lenBuffUI * sizeof(signed char));
         if(! *qBuffStrPtr)
            goto memErr_fun01_sec05;

         *sizeQUIPtr = lenBuffUI;
      } /*If: need new quence buffer*/
   } /*If: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - indel removal
   ^   o fun01 sec03 sub01:
   ^     - copy matches/snps/soft masked bases
   ^   o fun01 sec03 sub02:
   ^     - deal with hard masking cases
   ^   o fun01 sec03 sub03:
   ^     - find homopolymer length
   ^   o fun01 sec03 sub04:
   ^     - check if discarding indel
   ^   o fun01 sec03 sub05:
   ^     - remove indels
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec03 Sub01:
   *  - copy matches/snps/soft masked bases
   \*****************************************************/

   refPosUI = samSTPtr->refStartUI;

   for(
      cpCigUI = 0;
      cpCigUI < samSTPtr->lenCigUI;
      ++cpCigUI
   ){ /*Loop: remove indels*/

      samSTPtr->cigTypeStr[uiCig] =
         samSTPtr->cigTypeStr[cpCigUI];
      samSTPtr->cigArySI[uiCig] =
         samSTPtr->cigArySI[cpCigUI];

      if(
            samSTPtr->cigTypeStr[uiCig] == 'M'
         || samSTPtr->cigTypeStr[uiCig] == '='
         || samSTPtr->cigTypeStr[uiCig] == 'X'
         || samSTPtr->cigTypeStr[uiCig] == 'S'
      ){ /*If: non-indel mutation*/
      
         if(samSTPtr->cigTypeStr[uiCig] != 'S')
            refPosUI += samSTPtr->cigArySI[uiCig];

         cpLen_ulCp(
            &(*seqBuffStrPtr)[dupPosUI],
            &samSTPtr->seqStr[seqPosUI],
            samSTPtr->cigArySI[uiCig]
         );

         if(qBl)
         { /*If: have q-score entry*/
            cpLen_ulCp(
               &(*qBuffStrPtr)[dupPosUI],
               &samSTPtr->qStr[seqPosUI],
               samSTPtr->cigArySI[uiCig]
            );
         } /*If: have q-score entry*/

         dupPosUI += samSTPtr->cigArySI[uiCig];
         seqPosUI += samSTPtr->cigArySI[uiCig];
         ++uiCig;

         continue;

      }  /*If: non-indel mutation*/

      /**************************************************\
      * Fun01 Sec03 Sub02:
      *  - deal with hard masking cases
      \**************************************************/

      if(
            samSTPtr->cigTypeStr[uiCig] != 'D'
         && samSTPtr->cigTypeStr[uiCig] != 'I'
      ){ /*If: not an indel (hard masking)*/
         ++uiCig;
         continue;
      }  /*If: not an indel (hard masking)*/

      /**************************************************\
      * Fun01 Sec03 Sub03:
      *  - find homopolymer length
      \**************************************************/

      lenHomoUI = 1;

      if(refStr[refPosUI] != '\0')
      { /*If: have next bases (potentail homopolymer)*/
         tmpPosUI = refPosUI + 1;

         while(refStr[tmpPosUI] == refStr[refPosUI])
         { /*Loop: find end of homopolymer*/
            ++lenHomoUI;
            ++tmpPosUI;
         } /*Loop: find end of homopolymer*/
      } /*If: have next bases (potentail homopolymer)*/



      if(refPosUI > 0)
      { /*If: have previous bases (homopolymer?)*/
         tmpPosUI = refPosUI - 1;

         while(refStr[tmpPosUI] == refStr[refPosUI])
         { /*Loop: find end of homopolymer*/
            ++lenHomoUI;

            if(! tmpPosUI)
               break;

            --tmpPosUI;
         } /*Loop: find end of homopolymer*/
      } /*If: have previous bases (homopolymer?)*/

      /**************************************************\
      * Fun01 Sec03 Sub04:
      *  - check if discarding indel
      \**************************************************/

      if(
            lenHomoUI < (unsigned int) minHomoSI
         || samSTPtr->cigArySI[uiCig] > maxIndelSI
      ){ /*If: indel is to large or homoplymer to small*/
         if(samSTPtr->cigTypeStr[uiCig] == 'I')
         { /*If: need to copy insertion*/

            tmpPosUI = seqPosUI;
            tmpPosUI += samSTPtr->cigArySI[uiCig];

            while(seqPosUI < tmpPosUI)
            { /*Loop: copy insertion*/
               if(qBl)
                  (*qBuffStrPtr)[dupPosUI] =
                     samSTPtr->qStr[seqPosUI];

               (*seqBuffStrPtr)[dupPosUI++] =
                  samSTPtr->seqStr[seqPosUI++];
            } /*Loop: copy insertion*/
         } /*If: need to copy insertion*/

         ++uiCig; /*keep cigar entry*/
         continue;
      }  /*If: indel is to large or homoplymer to small*/
 
      /**************************************************\
      * Fun01 Sec03 Sub05:
      *  - remove indels
      \**************************************************/

      if(samSTPtr->cigTypeStr[uiCig] == 'D')
      { /*If: removing deletions*/
         tmpPosUI = refPosUI;
         tmpPosUI += samSTPtr->cigArySI[uiCig];

         while(refPosUI < tmpPosUI)
         { /*Loop: add reference bases in (remove del)*/
            if(qBl)
            { /*If: have q-score entries to add to*/
               /*I am using 0, becuase replacing a
               `  deletion with a reference base means no
               `  confidence in the call
               */
               ++samSTPtr->qHistUI[0];
               (*qBuffStrPtr)[dupPosUI] =
                  def_adjQ_samEntry; /*adjQ is q-score 0*/
            } /*If: have q-score entries to add to*/

            if(maskSC)
            { /*If: user wanted base masked*/
               (*seqBuffStrPtr)[dupPosUI++] = maskSC;
               ++refPosUI;
            } /*If: user wanted base masked*/

            else
               (*seqBuffStrPtr)[dupPosUI++] =
                  refStr[refPosUI++];
         } /*Loop: add reference bases in (remove del)*/

         samSTPtr->numDelUI -= samSTPtr->cigArySI[uiCig];
         samSTPtr->numMatchUI +=
            samSTPtr->cigArySI[uiCig];
         samSTPtr->readLenUI += samSTPtr->cigArySI[uiCig];


         tmpPosUI += samSTPtr->cigTypeStr[uiCig] = '=';

         if(samSTPtr->cigTypeStr[uiCig - 1] == '=')
         { /*If: this  cigar entry can be merged*/
            samSTPtr->cigArySI[uiCig - 1] +=
               samSTPtr->cigArySI[uiCig];
         }  /*If: this cigar entry can be merged*/

         else
            ++uiCig; /*can not merge*/
      } /*If: removing deletions*/


      else
      { /*Else: insertion*/
         samSTPtr->numInsUI -= samSTPtr->cigArySI[uiCig];
         samSTPtr->readLenUI -= samSTPtr->cigArySI[uiCig];

         tmpPosUI = seqPosUI;
         tmpPosUI += samSTPtr->cigArySI[uiCig];


         if(qBl)
         { /*If: have quality score entry*/
            while(seqPosUI < tmpPosUI)
            { /*Loop: remove insertion quality scores*/
               qUC =
                  (unsigned char)
                  samSTPtr->qStr[seqPosUI++];
               qUC -= def_adjQ_samEntry;

               --samSTPtr->qHistUI[qUC];
               samSTPtr->sumQUL -= qUC;
            } /*Loop: remove insertion quality scores*/
         } /*If: have quality score entry*/

         else
         { /*Else: no quality score entry*/
            seqPosUI = tmpPosUI;
            continue;
         } /*Else: no quality score entry*/
      } /*Else: insertion*/

   } /*Loop: remove indels*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^   - add null to end of c-strings (sequence/cigar)
   ^   o fun01 sec04 sub01:
   ^     - merge duplicate cigar entries
   ^   o fun01 sec04 sub02:
   ^     - copy sequence and q-score entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec04 Sub01:
   *   - merge duplicate cigar entries
   \*****************************************************/

   samSTPtr->lenCigUI = uiCig;
   uiCig = 0;

   for(
     cpCigUI = 0;
     cpCigUI < samSTPtr->lenCigUI;
     ++cpCigUI
   ){ /*Loop: merge duplicate cigar entries*/
      samSTPtr->cigTypeStr[uiCig] =
         samSTPtr->cigTypeStr[cpCigUI];
      samSTPtr->cigArySI[uiCig] =
         samSTPtr->cigArySI[cpCigUI];

      if(! cpCigUI)
      { /*If: first cigar entry*/
         ++uiCig;
         continue;
      } /*If: first cigar entry*/

      if(
            samSTPtr->cigTypeStr[uiCig - 1]
         == samSTPtr->cigTypeStr[uiCig]
      ){ /*If: can merge cigar entries*/
         samSTPtr->cigArySI[uiCig - 1] +=
            samSTPtr->cigArySI[uiCig];
      }  /*If: can merge cigar entries*/

      else
         ++uiCig;
   }  /*Loop: merge duplicate cigar entries*/

   samSTPtr->lenCigUI = uiCig;
   samSTPtr->cigArySI[uiCig] = '\0';

   /*****************************************************\
   * Fun01 Sec04 Sub02:
   *   - copy sequence and q-score entries
   \*****************************************************/

   (*seqBuffStrPtr)[dupPosUI] = '\0';
   (*qBuffStrPtr)[dupPosUI] = '\0';

   if(samSTPtr->lenSeqBuffUI < lenBuffUI)
   { /*If: new sequence buffer is larger*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         malloc(lenBuffUI * sizeof(signed char));
      if(! samSTPtr->seqStr)
         goto memErr_fun01_sec05;

      samSTPtr->lenSeqBuffUI = lenBuffUI - 8;
   } /*If: new sequence buffer is larger*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      (*seqBuffStrPtr),
      samSTPtr->readLenUI
   );

   if(qBl)
   { /*If: have q-score entry*/
      if(samSTPtr->lenQBuffUI < lenBuffUI)
      { /*If: current q-score buffer is larger*/
         free(samSTPtr->qStr);

         samSTPtr->qStr =
            malloc(lenBuffUI * sizeof(signed char));
         if(! samSTPtr->qStr)
            goto memErr_fun01_sec05;

         samSTPtr->lenQBuffUI = lenBuffUI - 8;
      } /*If: current q-score buffer is larger*/

      cpLen_ulCp(
         samSTPtr->qStr,
         (*qBuffStrPtr),
         samSTPtr->readLenUI
      );

      findQScores_samEntry(samSTPtr);
         /*need to re-adjust q-scores for removed indels*/
   } /*If: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec05:
   ^   - return and clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   qUC = 0;
   goto ret_fun01_sec05;

   memErr_fun01_sec05:;
      qUC = 1;
      goto ret_fun01_sec05;

   ret_fun01_sec05:;
      return (signed char) qUC;
} /*indel_rmHomo*/

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
