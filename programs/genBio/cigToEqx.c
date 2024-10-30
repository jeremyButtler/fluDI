/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' cigToEqx SOF: Start Of File
'   - functions to convert an regular sam cigar to an eqx
'     cigar
'   o header:
'     - included libraries and .c only defined variables
'   o fun01: samCigToEqxCig
'     - convert an sam entry cigar to an eqx cigar
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries and .c only defined varaibles
|   o header sec01:
|     - included libraries
|   o header sec02:
|     - defined variables
\-------------------------------------------------------*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec01:
^   - included libraries
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "cigToEqx.h"

#include <stdio.h>

#include "../genLib/ulCp.h"
#include "samEntry.h"

/*No .c files (.h only)*/
#include "../genLib/dataTypeShortHand.h"
#include "ntTo5Bit.h" /*look up tables*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .c  #include "../generalLib/base10str.h"
!   o .c  #include "../generalLib/numToStr.h"
!   o .c  #include "../generalLib/strAry.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec02:
^   - defined variables
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define def_charToSI_cigToEqx 2  /*char to int by shift*/

/*-------------------------------------------------------\
| Fun01: cigToEqx 
|   - convert an sam entry cigar to an eqx cigar
| Input:
|   - samSTPtr:
|     o pointer to an samEntry struct to convert cigar
|       into an eqx cigar entry
|   - refStr:
|     o pointer to an reference sequence to use in the
|       conversion
|   - strBuffPtr:
|     o pointer to an c-string buffer to hold the cigar
|       entry temporarly
|   - valArySI:
|     o pointer to an int array to hold the cigar
|       entry temporarly
|   - lenBuffSI:
|     o pointer to current length of strBuffPtr and
|       valArySI
|   - nAsSnpBl:
|     o 1: Treat anymous bases as SNPs always
|     o 0: Treat annoymous bases as matches, if they match
| Output:
|   - Modifies:
|     o cigTypeStr in samSTPtr to be an eqx entry
|     o cigArySI in samSTPtr to have counts for the
|       eqx entry
|     o strBuffPtr to be larger if needed
|     o valArySI to be larger if needed
|     o lenBuffSI to hold new buffer size if needed
|   - Returns:
|     o 0 for success
|     o 1 for memory error
|       - does not free strBuffPtr (you need to free)
|       - does not free valArySI   (you need to free)
\-------------------------------------------------------*/
signed char
cigToEqx(
   struct samEntry *samSTPtr,
   signed char *refStr,
   signed char **strBuffPtr,
   signed int **valArySI,
   signed int *lenBuffSI,
   signed char nAsSnpBl
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - convert an sam entry cigar to an eqx cigar
   '   o fun01 sec01:
   '     - variable declerations
   '   o fun01 sec02:
   '     - memory allocation for buffers
   '   o fun01 sec03:
   '     - make the new cigar entry
   '   o fun01 sec04:
   '     - copy new cigar entry to sam entry
   '   o fun01 sec05:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar *tmpStr = 0;

   sint siCig = 0;
   sint siNumNtCig = 0;
   sint siNewCig = 0;

   uint refPosUI = 0;
   uint seqPosUI = 0;

   schar refNtSC = 0;
   schar qryNtSC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - memory allocation for buffers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*check if I have an cigar entry*/
   if(samSTPtr->cigTypeStr[0] == '*')
      goto finishAndExit_fun01_sec02_sub01;

   /*check if I have an sequence entry*/
   if(samSTPtr->seqStr[0] == '*')
      goto finishAndExit_fun01_sec02_sub01;

   if(*strBuffPtr == 0)
   { /*If: I need to allocate memory*/
      *lenBuffSI = (sint) samSTPtr->lenCigUI << 1;

      *strBuffPtr =
         malloc((*lenBuffSI + 1) * sizeof(schar));

      if(! strBuffPtr)
         goto memErr_fun01_sec02_sub02;

      *valArySI = malloc(*lenBuffSI * sizeof(sint));

      if(! valArySI)
         goto memErr_fun01_sec02_sub02;
   } /*If: I need to allocate memory*/

   else if(*lenBuffSI < (sint) samSTPtr->readLenUI)
   { /*Else If: I need to reallocate memory*/
      *lenBuffSI = (sint) samSTPtr->readLenUI << 1;

      tmpStr = malloc((*lenBuffSI +1) * sizeof(schar));

      if(! tmpStr)
         goto memErr_fun01_sec02_sub02;

      free(*strBuffPtr);
      *strBuffPtr = 0;
      *strBuffPtr = tmpStr;

      tmpStr= (schar *) malloc(*lenBuffSI * sizeof(sint));

      if(! valArySI)
         goto memErr_fun01_sec02_sub02;

      free(*valArySI);
      *valArySI = 0;
      *valArySI = (sint *) tmpStr;
   } /*Else If: I need to reallocate memory*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - make the new cigar entry
   ^   o fun01 sec03 sub01:
   ^     - set up, start loop, and start cigar checks
   ^   o fun01 sec03 sub02:
   ^     - handle SNP cases
   ^   o fun01 sec03 sub03:
   ^     - handle match cases
   ^   o fun01 sec03 sub04:
   ^     - handle SNP or match cases
   ^   o fun01 sec03 sub05:
   ^     - handle insertion and sofmask cases
   ^   o fun01 sec03 sub06:
   ^     - handle deletion cases
   ^   o fun01 sec03 sub07:
   ^     - handle any other cases (hard mask likely)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec03 Sub01:
   *   - set up, start loop, and start cigar checks
   \*****************************************************/

   refPosUI = samSTPtr->refStartUI;
   siNewCig = 0; /*So gets incurmented to zero*/
   (*strBuffPtr)[0] = 0;

   for(
      siCig = 0;
      siCig < (sint) samSTPtr->lenCigUI;
      ++siCig
   ){ /*Loop: Check all cigar entries*/
      switch(samSTPtr->cigTypeStr[siCig])
      { /*Switch: Check the cigar entry type*/

         /***********************************************\
         * Fun01 Sec03 Sub02:
         *   - handle SNP cases
         \***********************************************/

         case 'X':
         /*Case: SNPs*/
            if(! nAsSnpBl)
            { /*If: I an converting N's to matches*/
               siNumNtCig = samSTPtr->cigArySI[siCig];

               while(siNumNtCig > 0)
               { /*Loop: break up by anonymous bases*/
                  qryNtSC =
                     ntTo5Bit[
                         (uchar)
                         samSTPtr->seqStr[seqPosUI]
                     ]; /*Convert query to code*/

                   refNtSC =
                     ntTo5Bit[
                        (uchar)
                        refStr[refPosUI]
                     ]; /*Convert reference to code*/

                  if(refNtSC & qryNtSC)
                  { /*If: Is match (aonymous bases)*/
                     if((*strBuffPtr)[siNewCig] != '=')
                     { /*If: Need to make an match entry*/
                        siNewCig+= ((*strBuffPtr)[0] !=0);

                        (*strBuffPtr)[siNewCig] = '=';
                        (*valArySI)[siNewCig] = 1;
                     } /*If: Need to make an match entry*/

                     else
                        ++( (*valArySI)[siNewCig] );
                  } /*If: Is match (aonymous bases)*/

                  else if((*strBuffPtr)[siNewCig] != 'X')
                  { /*Else If: started an new entry*/
                     siNewCig += ((*strBuffPtr)[0] != 0);

                     (*strBuffPtr)[siNewCig] = 'X';
                     (*valArySI)[siNewCig] = 1;
                  } /*Else If: started an new entry*/

                  else
                     ++( (*valArySI)[siNewCig] );

                  --siNumNtCig;
                  ++refPosUI;
                  ++seqPosUI;
               } /*Loop: break up by anonymous bases*/
            } /*If: I an converting N's to matches*/

            else
            { /*Else: I am doing no changes*/
               siNewCig += ((*strBuffPtr)[0] != 0);

               (*strBuffPtr)[siNewCig] = 'X';

               (*valArySI)[siNewCig] =
                  samSTPtr->cigArySI[siCig];

               refPosUI += samSTPtr->cigArySI[siCig];
               seqPosUI += samSTPtr->cigArySI[siCig];
            } /*Else: I am doing no changes*/

            break;
         /*Case: SNPs*/

         /***********************************************\
         * Fun01 Sec03 Sub03:
         *   - handle match cases
         \***********************************************/

         case '=':
         /*Case: Matches*/
            if(nAsSnpBl)
            { /*If: I an converting N's to matches*/
               siNumNtCig = samSTPtr->cigArySI[siCig];

               while(siNumNtCig > 0)
               { /*Loop: break up by anonymous bases*/
                  qryNtSC =
                     ntTo5Bit[
                         (uchar)
                         samSTPtr->seqStr[seqPosUI]
                     ];

                   refNtSC =
                      ntTo5Bit[
                         (uchar)
                         refStr[refPosUI]
                      ];

                  if(
                        (refNtSC & def_n_fithBit_ntTo5Bit)
                     || (qryNtSC & def_n_fithBit_ntTo5Bit)
                  ){ /*If: Is match (aonymous bases)*/
                     if((*strBuffPtr)[siNewCig] != 'X')
                     { /*If: Need to make an SNP entry*/
                        siNewCig+= ((*strBuffPtr)[0] !=0);
                        (*strBuffPtr)[siNewCig] = 'X';
                        (*valArySI)[siNewCig] = 1;
                     } /*If: Need to make an SNP entry*/

                     else
                        ++( (*valArySI)[siNewCig] );
                  } /*If: Is match (aonymous bases)*/

                  else if((*strBuffPtr)[siNewCig] != '=')
                  { /*Else If: neither base is an N*/
                      siNewCig += ((*strBuffPtr)[0] != 0);
                     (*strBuffPtr)[siNewCig] = '=';
                     (*valArySI)[siNewCig] = 1;
                  } /*Else If: neither base is an N*/

                  else
                     ++( (*valArySI)[siNewCig] );

                  --siNumNtCig;
                  ++refPosUI;
                  ++seqPosUI;
               } /*Loop: break up by anonymous bases*/
            } /*If: I an converting N's to matches*/

            else
            { /*Else: I am doing no changes*/
               siNewCig += ((*strBuffPtr)[0] != 0);

               (*strBuffPtr)[siNewCig] = '=';

               (*valArySI)[siNewCig] =
                  samSTPtr->cigArySI[siCig];

               refPosUI += samSTPtr->cigArySI[siCig];
               seqPosUI += samSTPtr->cigArySI[siCig];
            } /*Else: I am doing no changes*/

            break;
         /*Case: Matches*/

         /***********************************************\
         * Fun01 Sec03 Sub04:
         *   - handle SNP or match cases
         *   o fun01 sec03 sub04 cat01:
         *     - Start loop and get sixBit code
         *   o fun01 sec03 sub04 cat02:
         *     - mismatch regardless if anonymous
         *   o fun01 sec03 sub04 cat03:
         *     - match, one base is anonymous
         *   o fun01 sec03 sub04 cat04:
         *     - match, nothing is anonymous
         *   o fun01 sec03 sub04 cat05:
         *     - move to next base
         *   o fun01 sec03 sub04 cat06:
         *     - move to next cigar entry and restart
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun01 Sec03 Sub04 Cat01:
         +   - Start loop and get sixBit code
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'M':
         /*Case: Matches and SNPs*/
            siNumNtCig = samSTPtr->cigArySI[siCig];

            while(siNumNtCig > 0)
            { /*Loop: break up by anonymous bases*/
               qryNtSC =
                  ntTo5Bit[
                      (uchar)
                      samSTPtr->seqStr[seqPosUI]
                  ];

                refNtSC =
                   ntTo5Bit[
                      (uchar)
                      refStr[refPosUI]
                   ];

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun01 Sec03 Sub04 Cat02:
               +   - mismatch regardless if anonymous
               \++++++++++++++++++++++++++++++++++++++++*/

               if(! (refNtSC & qryNtSC) )
               { /*If: This is not an match*/
                  if((*strBuffPtr)[siNewCig] != 'X')
                  { /*If: Need to make an SNP entry*/
                     siNewCig += ((*strBuffPtr)[0] != 0);
                     (*strBuffPtr)[siNewCig] = 'X';
                     (*valArySI)[siNewCig] = 1;
                  } /*If: Need to make an SNP entry*/

                  else
                     ++( (*valArySI)[siNewCig] );
               } /*If: This is not an match*/

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun01 Sec03 Sub04 Cat03:
               +   - match, one base is anonymous
               \++++++++++++++++++++++++++++++++++++++++*/

               else if(
                     (refNtSC & def_n_fithBit_ntTo5Bit)
                  || (qryNtSC & def_n_fithBit_ntTo5Bit)
               ){ /*Else If: Is aonymous:*/
                  if(nAsSnpBl)
                  { /*If: anonymous bases are SNPs*/
                     if((*strBuffPtr)[siNewCig] != 'X')
                     { /*If: Need to make an SNP entry*/
                        siNewCig+= ((*strBuffPtr)[0] !=0);
                        (*strBuffPtr)[siNewCig] = 'X';
                        (*valArySI)[siNewCig] = 1;
                     } /*If: Need to make an SNP entry*/

                     else
                        ++( (*valArySI)[siNewCig] );
                  } /*If: anonymous bases are SNPs*/

                  else
                  { /*Anonymous bases are matches*/
                     if((*strBuffPtr)[siNewCig] != '=')
                     { /*Else If: neither base is an N*/
                        siNewCig+= ((*strBuffPtr)[0] !=0);
                        (*strBuffPtr)[siNewCig] = '=';
                        (*valArySI)[siNewCig] = 1;
                     } /*Else If: neither base is an N*/

                     else
                        ++( (*valArySI)[siNewCig] );
                  } /*Anonymous bases are matches*/
               } /*Else If: is anonymous*/

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun01 Sec03 Sub04 Cat04:
               +   - match, nothing is anonymous
               \++++++++++++++++++++++++++++++++++++++++*/

               else
               { /*Else: Is an non-aonymous match*/
                  if((*strBuffPtr)[siNewCig] != '=')
                  { /*Else If: neither base is an N*/
                     siNewCig += ((*strBuffPtr)[0] != 0);
                     (*strBuffPtr)[siNewCig] = '=';
                     (*valArySI)[siNewCig] = 1;
                  } /*Else If: neither base is an N*/

                  else
                     ++( (*valArySI)[siNewCig] );
               } /*Else: Is an non-aonymous match*/

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun01 Sec03 Sub04 Cat05:
               +   - move to next base
               \++++++++++++++++++++++++++++++++++++++++*/

               --siNumNtCig;
               ++refPosUI;
               ++seqPosUI;
            } /*Loop: break up by anonymous bases*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun01 Sec03 Sub04 Cat06:
            +   - move to next cigar entry and restart
            \+++++++++++++++++++++++++++++++++++++++++++*/

            break;
         /*Case: Matches and SNPs*/

         /***********************************************\
         * Fun01 Sec03 Sub05:
         *   - handle insertion and sofmask cases
         \***********************************************/

         case 'S':
         case 'I':
         /*Case: insertion or softmasking*/
            siNewCig += ((*strBuffPtr)[0] != 0);

            (*strBuffPtr)[siNewCig] =
               samSTPtr->cigTypeStr[siCig];

            (*valArySI)[siNewCig] =
               samSTPtr->cigArySI[siCig];

            seqPosUI += samSTPtr->cigArySI[siCig];
            break;
         /*Case: insertion or softmasking*/

         /***********************************************\
         * Fun01 Sec03 Sub06:
         *   - handle deletion cases
         \***********************************************/

         case 'D':
         /*Case: deletion*/
            siNewCig += ((*strBuffPtr)[0] != 0);

            (*strBuffPtr)[siNewCig] = 'D';

            (*valArySI)[siNewCig] =
               samSTPtr->cigArySI[siCig];

            refPosUI += samSTPtr->cigArySI[siCig];
            break;
         /*Case: deletion*/

         /***********************************************\
         * Fun01 Sec03 Sub07:
         *   - handle any other cases (hard mask likely)
         \***********************************************/

         default:
         /*Case: hard mask of some kind*/
            siNewCig += ((*strBuffPtr)[0] != 0);

            (*strBuffPtr)[siNewCig] =
               samSTPtr->cigTypeStr[siCig];

            (*valArySI)[siNewCig] =
               samSTPtr->cigArySI[siCig];

            break;
         /*Case: hard mask of some kind*/
      } /*Switch: Check the cigar entry type*/
   } /*Loop: Check all cigar entries*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^   - copy new cigar entry to sam entry
   ^   o fun01 sec04 sub01:
   ^     - if needed allocate more memory for new cigar
   ^   o fun01 sec04 sub02:
   ^     - copy the new cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec04 Sub01:
   *   - if needed allocate more memory for new cigar
   \*****************************************************/

   ++siNewCig; /*Convert index 0 to index 1*/

   if(samSTPtr->lenCigBuffUI < (uint) (siNewCig + 1))
   { /*If: I need to resize the cigar arrays*/
      samSTPtr->lenCigBuffUI = (uint) siNewCig << 1;

      tmpStr =
         (schar *)
         malloc(samSTPtr->lenCigBuffUI * sizeof(char));

      if(! tmpStr)
         goto memErr_fun01_sec02_sub02;

      free(samSTPtr->cigTypeStr);
      samSTPtr->cigTypeStr = 0;
      samSTPtr->cigTypeStr = tmpStr;

      tmpStr =
         (schar *)
         malloc(samSTPtr->lenCigBuffUI * sizeof(sint));

      if(! tmpStr)
         goto memErr_fun01_sec02_sub02;

      free(samSTPtr->cigArySI);
      samSTPtr->cigArySI = 0;
      samSTPtr->cigArySI = (sint *) tmpStr;
   } /*If: I need to resize the cigar arrays*/

   /*****************************************************\
   * Fun01 Sec04 Sub02:
   *   - copy the new cigar entry
   \*****************************************************/

   samSTPtr->lenCigUI = (uint) siNewCig;

   for(siCig = 0; siCig < siNewCig; ++siCig)
   { /*Loop: Copy the new cigar over*/
      samSTPtr->cigTypeStr[siCig] = (*strBuffPtr)[siCig];
      samSTPtr->cigArySI[siCig] = (*valArySI)[siCig];
   } /*Loop: Copy the new cigar over*/

   samSTPtr->cigTypeStr[siCig] = '\0';
   samSTPtr->cigArySI[siCig] = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec05:
   ^   - clean up
   ^   o fun01 sec05 sub01:
   ^     - clean up after success
   ^   o fun01 sec05 sub02:
   ^     - error clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec05 Sub01:
   *   - clean up after success
   \*****************************************************/

   finishAndExit_fun01_sec02_sub01:;

   return 0;

   /*****************************************************\
   * Fun01 Sec05 Sub02:
   *   - error clean up
   \*****************************************************/

   memErr_fun01_sec02_sub02:;

   return 1;
} /*samCigToEqxCig*/


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
