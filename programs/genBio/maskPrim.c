/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' maskPrim SOF: Start Of File
'   - holds functions for primer masking or trimming with
'     sam files
'   o header:
'     - defined variables and guards
'   o .h fun01: blankFlag_maskPrim
'     - blanks all flags in an primFlag uint
'   o .h fun02: setDir_maskPrim
'     - sets the direction flag for an primer
'   o .h fun03: setPair_maskPrim
'     - sets the paired flag for an primer
'   o .h fun04: setMateIndexPrimIndex
'     - sets the index of the other primer
'   o .h fun05: getDir_maskPrim
'     - gets the primer direction from an primer flag uint
'   o .h fun06: getPair_maskPrim
'     - gets if the primer is to be treated as an pair
'   o .h fun07: getMateIndex_maskPrim
'     - gets the other primers index in the pair
'   o fun08: maskPrim
'     - masks the primers in the sam sequence
'   o .c fun09: swapPrimCoords
'     - swaps two primer coordinates values around and
'       updates the index's
'   o fun09: sortCoords_maskPrim
'     - Sorts the primer coordiantes of a primer flag
'       array setup. Order is least to greates by first
'       starting positon, then ending position
'   o fun11: getCoords_maskPrim
'     - gets primer coordinates and flags from an file
'   o license:
'     - Licensing for this code (public domain / mit)
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

#include "maskPrim.h"

#include <stdio.h>

#include "../genLib/ulCp.h"
#include "../genLib/base10str.h"
#include "../genLib/shellSort.h"
#include "samEntry.h"

/*No .c files (only .h)*/
#include "../genLib/dataTypeShortHand.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries
!   o .c  #include "../genLib/numToStr.h"
!   o .c  #include "../genLib/strAry.h"
!   o .h  #include "ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun08: maskPrim
|   - masks the primers in the sam sequence
| Input:
|   - samSTPtr:
|     o pointer to an sam structure to mask primers in
|   - primStartAryUI:
|     o uint array of primer starting coordiantes
|     o starting position on the reference sequence
|   - primEndAryUI:
|     o uint array of primer ending coordiantes
|     o ending position on the reference sequence
|   - primFlagAryUI:
|     o int array with bits filled for the flags and
|       coordinates
|       - first two bits have the direction,
|         o def_reverse_maskPrim (1) if is an reverse
|           primer
|         o def_pair_maskPrim (2) if you want to only
|           mask primers that have pairs instead off any
|           primer
|       - next 15 bits (3 to 17) index of forward primer
|       - last 15 bits (18 to 32) index of reverse primer
|   - numPrimSI:
|     o number of primers in primStartAryUI and primEndUI
|   - fudgeSI:
|     o sint with the number of bases before (start) or
|       after (end) to 
|       - -1 mask any detected primer
|       - 0 is no bases before start or after end
|       - > 0 Allow 
|   - maskSC:
|     o character to mask with
| Output:
|   - Modifies:
|     o samSTPtr to have any detected primers masked
|   - Returns:
|     o 0 if no problems
|     o -1 for no detected primers
|     o (primer index << 1) | 1 if primer found but not
|        masked. This can only happen for fudge lengths
|        out of bounds or incomplete primer pairs
\-------------------------------------------------------*/
signed int
maskPrim(
   struct samEntry *samSTPtr,
   unsigned int *primStartAryUI,
   unsigned int *primEndAryUI,
   unsigned int *primFlagAryUI,
   signed int numPrimSI,
   signed int fudgeSI, /*max fudge size*/
   signed char maskSC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun08 TOC:
   '   - mask the primers in an sam sequence
   '   o fun08 sec01:
   '     - variable declerations
   '   o fun08 sec02:
   '     - mask regions in the read with primers
   '   o fun08 sec03:
   '     - return back the answer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar revPrimBl = 0;
   schar pairPrimBl = 0;
   sint primIndexSI = 0;
   sint revIndexSI = -1;

   uint forStartUI = 0;
   uint forEndUI = 0;
   uint revStartUI = 0;
   uint revEndUI = 0;

   char firstMaskBl = 1; /*No primers masked yet*/

   /*Values for finding start position*/
   uint uiCig = 0;
   uint cigBaseLeftUI = 0;

   uint refPosUI = 0;
   uint seqPosUI = 0;

   /*Are for finding end of mask region*/
   uint uiEndCig = 0;
   uint cigEndBaseLeftUI = 0;
   uint refEndPosUI = 0;
   uint seqEndPosUI = 0;

   sint retIndexSI = 0; /*For error returing*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec02:
   ^   - mask regions in the read with primers
   ^   o fun08 sec02 sub01:
   ^     - find 1st primers index, start loop, check range
   ^   o fun08 sec02 sub02:
   ^     - check primers paired, get reverse primer
   ^   o fun08 sec02 sub03:
   ^     - check if the primer(s) are in fudge length
   ^   o fun08 sec02 sub04:
   ^     - find start of primer and mask
   ^   o fun08 sec02 sub05:
   ^     - mask the primer position
   ^   o fun08 sec02 sub06:
   ^     - for paired primers, set up reverse primer mask
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   /*****************************************************\
   * Fun08 Sec02 Sub01:
   *   - find first primers index, start loop, check range
   \*****************************************************/

   cigBaseLeftUI = samSTPtr->cigArySI[0];
   refPosUI = samSTPtr->refStartUI;

   /*Account for some primer positions not being
   `  complete. I can always loop to the start later if
   `  needed
   */
   if(samSTPtr->refStartUI > 20)
      forStartUI = samSTPtr->refStartUI - 20;
   else
      forStartUI = 0;

   primIndexSI =
      (sint)
      uiRange_shellSort(
         primStartAryUI,
         forStartUI,
         samSTPtr->refEndUI,
         (signed long) numPrimSI
      );

   if(primIndexSI < 0)
      goto noMask_fun08_sec03_sub02;

   while(primIndexSI < numPrimSI)
   { /*Loop: Mask primers*/
      forStartUI = primStartAryUI[primIndexSI];
      forEndUI = primEndAryUI[primIndexSI];

      /*Only fires after the first case*/
      if(forStartUI > samSTPtr->refEndUI)
      { /*If: the primer is outside of my range*/
         if(firstMaskBl)
            goto noMask_fun08_sec03_sub02;

         break;
      } /*If: the primer is outside of my range*/

      /*If primer is before sequence*/
      if(forEndUI < samSTPtr->refStartUI)
         goto nextPrim_fun08_sec02_sub07;
      
      /**************************************************\
      * Fun08 Sec02 Sub02:
      *   - check primers paired, get reverse primer
      \**************************************************/

      pairPrimBl =
         getPair_maskPrim(primFlagAryUI[primIndexSI]);

      revPrimBl =
         getDir_maskPrim(primFlagAryUI[primIndexSI]);

      if(pairPrimBl)
      { /*If: These primers are paired*/
         /*If this is an pair, then I have already checked
         `  the reverse complement or will check the
         `  the reverse complement primer
         */
         if(revPrimBl)
            goto nextPrim_fun08_sec02_sub07;

         revIndexSI =
            getMateIndex_maskPrim(
               primFlagAryUI[primIndexSI]
            ); /*Get the index of the reverse primer*/

         revStartUI = primStartAryUI[revIndexSI];
         revEndUI = primEndAryUI[revIndexSI];

         if(revStartUI > samSTPtr->refEndUI)
         { /*If: reverse primer is outside of range*/
            if(firstMaskBl)
               goto noMask_fun08_sec03_sub02;

            break;
         } /*If: reverse primer is outside of range*/
      } /*If: These primers are paired*/

      /**************************************************\
      * Fun08 Sec02 Sub03:
      *   - check if the primer(s) are in fudge length
      \**************************************************/

      if(fudgeSI > -1)
      { /*If: The user only wants end primers*/
         /*I need to make sure this primer is at the end
         `  of the read + or - users fudge length
         */
         if(pairPrimBl)
         { /*If: paired primer*/
            if(forStartUI -fudgeSI > samSTPtr->refStartUI)
               goto nextPrim_fun08_sec02_sub07;

            if(forStartUI +fudgeSI < samSTPtr->refStartUI)
               goto nextPrim_fun08_sec02_sub07;

            if(revEndUI + fudgeSI < samSTPtr->refEndUI)
               goto nextPrim_fun08_sec02_sub07;

            if(revEndUI - fudgeSI > samSTPtr->refEndUI)
               goto nextPrim_fun08_sec02_sub07;
         } /*If: paired primer*/

         else if(revPrimBl)
         { /*Else If: unpaired reverse primer*/
            if(forEndUI + fudgeSI < samSTPtr->refEndUI)
               goto nextPrim_fun08_sec02_sub07;

            if(forEndUI - fudgeSI > samSTPtr->refEndUI)
               goto nextPrim_fun08_sec02_sub07;
         } /*Else If: unpaired reverse primer*/

         else
         { /*Else: is unpaired forward primer*/
            if(
               forStartUI - fudgeSI > samSTPtr->refStartUI
            ) goto nextPrim_fun08_sec02_sub07;

            if(
               forStartUI + fudgeSI < samSTPtr->refStartUI
            ) goto nextPrim_fun08_sec02_sub07;
         } /*Else: is unpaired forward primer*/
      } /*If: The user only wants end primers*/

      /**************************************************\
      * Fun08 Sec02 Sub04:
      *   - find the start of the primer
      \**************************************************/

      firstMaskBl = 0; /*At this point masking happens*/

      findRefPos_samEntry(
         samSTPtr, 
         (sint *) &uiCig,
         (sint *) &cigBaseLeftUI,
         forStartUI,
         (sint *) &refPosUI,
         (sint *) &seqPosUI
      );

      uiEndCig = uiCig;
      cigEndBaseLeftUI = cigBaseLeftUI;
      refEndPosUI = refPosUI;
      seqEndPosUI = seqPosUI;

      /**************************************************\
      * Fun08 Sec02 Sub05:
      *   - mask the primer position
      \**************************************************/

      maskLoop_fun08_sec02_sub05:;

      while(refEndPosUI <= forEndUI)
      { /*Loop: Mask bases and adjust cigar*/
         if(uiEndCig >= samSTPtr->lenCigUI)
            break; /*sequence had part of primer only*/

         if(refEndPosUI > samSTPtr->refEndUI)
            break;

         switch(samSTPtr->cigTypeStr[uiEndCig])
         { /*Switch: Find the cigar entry type*/
            case 'M':
            case 'X':
            case '=':
            /*Case: match/snp*/
               while(cigEndBaseLeftUI > 0)
               { /*Loop: Mask the bases*/
                  if(refEndPosUI > forEndUI)
                     goto checkPair_fun08_sec02_sub06;

                  samSTPtr->seqStr[seqEndPosUI] = maskSC;
                  ++seqEndPosUI;
                  ++refEndPosUI;
                  --cigEndBaseLeftUI;
               } /*Loop: Mask the bases*/

               if(cigEndBaseLeftUI < 1)
               { /*If: I used up all the bases*/
                  ++uiEndCig;

                  cigEndBaseLeftUI =
                     samSTPtr->cigArySI[uiEndCig];
               } /*If: I used up all the bases*/

               break;
            /*Case: match/snp*/

            case 'I':
            /*Case: insertions*/
               while(cigEndBaseLeftUI > 0)
               { /*Loop: Mask the bases*/
                  samSTPtr->seqStr[seqEndPosUI] = maskSC;
                  ++seqEndPosUI;
                  --cigEndBaseLeftUI;
               } /*Loop: Mask the bases*/

               ++uiEndCig;

               cigEndBaseLeftUI =
                     samSTPtr->cigArySI[uiEndCig];

               break;
            /*Case: insertions*/

            case 'D':
            /*Case: Deletions*/
               refEndPosUI += cigEndBaseLeftUI;

               if(refEndPosUI > (forEndUI + 1))
               { /*If: I have overshot*/
                   cigEndBaseLeftUI =
                      refEndPosUI - forEndUI - 1;

                   refEndPosUI = forEndUI + 1;
               } /*If: I have overshot*/

               else
               { /*Else: I still have more masking*/
                  ++uiEndCig;

                  cigEndBaseLeftUI =
                     samSTPtr->cigArySI[uiEndCig];
               } /*Else: I still have more masking*/

               break;
            /*Case: Deletions*/

            case 'S':
            /*Case: softmasks*/
               while(cigEndBaseLeftUI > 0)
               { /*Loop: Mask the bases*/
                  samSTPtr->seqStr[seqEndPosUI] = maskSC;
                  ++seqEndPosUI;
                  --cigEndBaseLeftUI;
               } /*Loop: Mask the bases*/

               ++uiEndCig;

               cigEndBaseLeftUI =
                     samSTPtr->cigArySI[uiEndCig];

               break;
            /*Case: softmasks*/

            default:
            /*Case: hard mask of some kind*/
               ++uiEndCig;

               cigEndBaseLeftUI =
                     samSTPtr->cigArySI[uiEndCig];

               break;
            /*Case: hard mask of some kind*/
         } /*Switch: Find the cigar entry type*/
      } /*Loop: Mask bases and adjust cigar*/

      /**************************************************\
      * Fun08 Sec02 Sub06:
      *   - for paired primers, set up reverse primer mask
      \**************************************************/

      checkPair_fun08_sec02_sub06:;

      if(pairPrimBl)
      { /*If: I have primer pairs to check*/
         forStartUI = revStartUI;
         forEndUI = revEndUI;

         findRefPos_samEntry(
            samSTPtr, 
            (sint *) &uiEndCig,
            (sint *) &cigEndBaseLeftUI,
            forStartUI,
            (sint *) &refEndPosUI,
            (sint *) &seqEndPosUI
         );

         pairPrimBl = 0;

         goto maskLoop_fun08_sec02_sub05;
      } /*If: I have primer pairs to check*/

      /**************************************************\
      * Fun08 Sec02 Sub07:
      *   - move to the next set of primers
      \**************************************************/

      nextPrim_fun08_sec02_sub07:;

      ++primIndexSI;
   } /*Loop: Mask primers*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec03:
   ^   - return back the answer
   ^   o fun08 sec03 sub01:
   ^     - return success
   ^   o fun08 sec03 sub02:
   ^     - return no primers masked
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun08 Sec03 Sub01:
   *   - return success
   \*****************************************************/

   return -firstMaskBl; /*Masking was done*/
   /*-1 if firstMaskBl is 1; else 0*/

   /*****************************************************\
   * Fun08 Sec03 Sub02:
   *   - return no primers masked
   \*****************************************************/

   noMask_fun08_sec03_sub02:;

   if(primIndexSI < 0) 
      retIndexSI = (sint) -1;
   else
      retIndexSI = (((sint) primIndexSI) << 1) | 1;

   return retIndexSI;
} /*maskPrim*/

/*-------------------------------------------------------\
| Fun09: sortCoords_maskPrim
|  - sorts the primer coordiantes of a primer flag array
|    setup. Order is least to greates by first starting
|    positon, then ending position
| Input:
|  - startAryUI:
|    o starting coordinates for the primers
|  - endAryUI:
|    o ending coordinates for the primers
|  - flagAryUI:
|    o array full of primFlags with indexs to keep up to
|      date
|  - lenAryUI:
|    o length of the array
| Output:
|  - Modifies:
|    o startAry array to be sorted form least to greatest
|    o endAry to be kept in order with firstAry
|    o flagAry to be kept in order with firstAry and have
|      any index's updated
\-------------------------------------------------------*/
void
sortCoords_maskPrim(
   unsigned int *startAryUI,
   unsigned int *endAryUI,
   unsigned int *flagAryUI,
   unsigned int lenAryUI 
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun09 TOC:
   '  - sorts the primer coordiantes of a primer flag
   '    array setup. Order is least to greates by first
   '    starting positon, then ending position
   '  - Shell sort taken from:
   '    - Adam Drozdek. 2013. Data Structures and
   '      Algorithims in c++. Cengage Leraning. fourth
   '      edition. pages 505-508
   '    - i made some minor changes, but is mostly the
   '      same
   '  o fun09 sec01:
   '    - variable declerations
   '  o fun09 sec02:
   '    - find the number of rounds to sort for
   '  o fun09 sec03:
   '    - sort the coordiantes
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec01:
  ^  - variable declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  
  /*number of sorting rounds*/
  uint subUI = 0;
  uint nextUI = 0;
  uint lastUI = 0;
  uint elmOnUI = 0;
  
  /*Variables to incurment loops*/
  uint uiIndex = 0;
  uint uiElm = 0;
 
  uint swapUI = 0;
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec02:
  ^  - Find the max search value (number rounds to sort)
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  
  /*recursion formula: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
  subUI = 1; /*Initialzie first array*/
  
  while(subUI < lenAryUI)
     subUI = (3 * subUI) + 1;
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec03:
  ^   - sort the coordinates
  ^   o fun09 sec03 sub01:
  ^     - move though sub arrays
  ^   o fun09 sec03 sub02:
  ^     - swap array elments (if needed)
  ^   o fun09 sec03 sub03:
  ^     - do insertion sort on subarray
  ^   o fun09 sec03 sub04:
  ^     - move to next sub array
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  
  /******************************************************\
  * Fun09 Sec03 Sub01:
  *   - move though sub arrays
  \******************************************************/

  while(subUI > 0)
  { /*loop trhough all sub arrays sort the subarrays*/

    for( 
       uiIndex = 0;
       uiIndex <= subUI;
       ++uiIndex
    ){ /*For each element in the subarray*/

      for(uiElm = uiIndex;
          uiElm + subUI < lenAryUI;
          uiElm += subUI
      ){ /*Loop; swap each nth element of the subarray*/
        nextUI = uiElm + subUI;
        
        /************************************************\
        * Fun09 Sec03 Sub02:
        *   - swap array elments (if needed)
        \************************************************/

        if(
              startAryUI[uiElm] > startAryUI[nextUI]
           || endAryUI[uiElm] > endAryUI[nextUI]
        ){ /*If I need to swap an element*/

           startAryUI[uiElm] ^= startAryUI[nextUI];
           startAryUI[nextUI] ^= startAryUI[uiElm];
           startAryUI[uiElm] ^= startAryUI[nextUI];

           endAryUI[uiElm] ^= endAryUI[nextUI];
           endAryUI[nextUI] ^= endAryUI[uiElm];
           endAryUI[uiElm] ^= endAryUI[nextUI];

           /*Changing the paired index. By doing this step
           `  before the swap I can keep self pointers
           `  pointing to self. use self pointers
           `  (index's) for single sequence entries
           */
           swapUI =
              getMateIndex_maskPrim(flagAryUI[uiElm]);

           setMateIndex_maskPrim(
              flagAryUI[swapUI],
              nextUI
           ); /*Update the new mate index*/

           swapUI =
              getMateIndex_maskPrim(flagAryUI[nextUI]);

           setMateIndex_maskPrim(
              flagAryUI[swapUI],
              uiElm
           ); /*Update the new mate index*/

           flagAryUI[uiElm] ^= flagAryUI[nextUI];
           flagAryUI[nextUI] ^= flagAryUI[uiElm];
           flagAryUI[uiElm] ^= flagAryUI[nextUI];

           /************************************************\
           * Fun09 Sec03 Sub03:
           *   - do insertion sort on subarray
           \************************************************/

           lastUI = uiElm;
           elmOnUI = uiElm;
          
           while(lastUI >= subUI)
           { /*loop; move swapped element back*/
              lastUI -= subUI;
              
              if(
                startAryUI[elmOnUI] >startAryUI[lastUI]
               || endAryUI[elmOnUI] > endAryUI[lastUI]
              ) break; /*Positioned the element*/
              

              startAryUI[elmOnUI] ^= startAryUI[lastUI];
              startAryUI[lastUI] ^= startAryUI[elmOnUI];
              startAryUI[elmOnUI] ^= startAryUI[lastUI];

              endAryUI[elmOnUI] ^= endAryUI[lastUI];
              endAryUI[lastUI] ^= endAryUI[elmOnUI];
              endAryUI[elmOnUI] ^= endAryUI[lastUI];

              /*Changing paired index. By doing this step
              `  before the swap I can keep self pointers
              `  pointing to self. use self pointers
              `  (index's) for single sequence entries
              */
              swapUI =
                getMateIndex_maskPrim(flagAryUI[elmOnUI]);

              setMateIndex_maskPrim(
                 flagAryUI[swapUI],
                 lastUI
              ); /*Update the new mate index*/

              swapUI =
                 getMateIndex_maskPrim(flagAryUI[lastUI]);

              setMateIndex_maskPrim(
                 flagAryUI[swapUI],
                 elmOnUI
              ); /*Update the new mate index*/

              flagAryUI[elmOnUI] ^= flagAryUI[lastUI];
              flagAryUI[lastUI] ^= flagAryUI[elmOnUI];
              flagAryUI[elmOnUI] ^= flagAryUI[lastUI];
              
              elmOnUI = lastUI;
           } /*loop; move swapped element back*/
        } /*If I need to swap elements*/
      } /*Loop; swap each nth element of the subarray*/
    } /*For each element in the subarray*/
    
    /****************************************************\
    * Fun09 Sec03 Sub04:
    *   - move to next sub array
    \****************************************************/

    subUI = (subUI - 1) / 3; /*Move to next round*/
  } /*loop through all sub arrays to sort the subarrays*/
} /*sortCoords_maskPrim*/

/*-------------------------------------------------------\
| Fun11: getCoords_maskPrim
|   - gets primer coordinates and flags from an file
| Input:
|   o primFileStr:
|     - c-string with name of file to get primer
|       coordinates from
|       o ignored\tignored\tpair\tstart\tend\tstart\tend
|         - pair: 1, P, T, or + if primers paired
|         - start: start coordinate of foward/reverse
|           primer
|         - end: end coordinate of foward/reverse primer
|   o startAryUI:
|     - pointer to an unsigned int pointer to hold the
|       starting coordinates for each primer
|   o endAryUI:
|     - pointer to an unsigned int pointer to hold the
|       ending coordinates for each primer
|   o flagAryUI:
|     - pointer to an unsigned int pointer to hold the
|       flags and pair index for each primer
|   o errSL:
|     - pointer to unsigned long to hold the error
|       message
| Output:
|   - Modifies:
|     o startAryUI to have the starting coordiantes
|     o endAryUI to have the ending coordiantes
|     o flagAryUI to have the direction and pair flags
|       and the index of the other primer
|     o errSL to hold the error
|       - 0 for no errors
|       - def_memErr_maskPrim for an memory error
|       - (siLine << 8) | def_fileErr_maskPrim for file
|         errors
|       - def_emptyFileErr_maskPrim for empty file
|  - Returns:
|    o 0 for errors
|    o number of primers extracted
\-------------------------------------------------------*/
signed int
getCoords_maskPrim(
   signed char *primFileStr,
   unsigned int **startAryUI,
   unsigned int **endAryUI,
   unsigned int **flagAryUI,
   signed long *errSL
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC:
   '   - gets primer coordinates and flags from an file
   '   o fun11 sec01:
   '     - Variable declerations
   '   o fun11 sec02:
   '     - Get the number of primers (lines * 2) in file
   '   o fun11 sec03:
   '     - Allocate memory for the primers
   '   o fun11 sec04:
   '     - Read in the primers
   '   o fun11 sec05:
   '     - Clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^   - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define lenBuffUS 1024
   signed char buffStr[lenBuffUS];
   signed char *tmpStr = 0;

   signed int numPrimSI = 0;
   signed int siLine = 0;

   schar revBl = 0; /*On the reverse primer*/

   FILE *primFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec02:
   ^   - Get the number of primers (lines * 2) in file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   primFILE =
      fopen(
         (char *) primFileStr,
         "r"
      );

   if(!primFILE)
       goto fileErr_fun11_sec05_sub03;

   /*do not expect line lengths be over 1024 characters*/
   while(fgets((char *) buffStr, lenBuffUS, primFILE))
      numPrimSI += 2;

   numPrimSI -= 2; /*Account for the header*/

   fseek(
      primFILE,
      0,
      SEEK_SET
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec03:
   ^   - Allocate memory for the primers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numPrimSI < 1)
      goto emptyFileErr_fun11_sec05_sub04;
      
   if(*startAryUI)
      free(*startAryUI);

   *startAryUI = 0;

   *startAryUI = malloc(numPrimSI * sizeof(uint));

   if(! *startAryUI)
      goto memErr_fun11_sec05_sub02;

   if(*endAryUI)
      free(*endAryUI);

   *endAryUI = 0;

   *endAryUI = malloc(numPrimSI * sizeof(uint));

   if(! *endAryUI)
       goto memErr_fun11_sec05_sub02;

   if(*flagAryUI)
      free(*flagAryUI);

   *flagAryUI = 0;

   *flagAryUI = malloc(numPrimSI * sizeof(uint));

   if(! *flagAryUI)
      goto memErr_fun11_sec05_sub02;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec04:
   ^   - Read in the primers
   ^   o fun11 sec04 sub01:
   ^     - Get past header and start primer read loop
   ^   o fun11 sec04 sub02:
   ^     - get past reference name
   ^   o fun11 sec04 sub03:
   ^     - get past primer name
   ^   o fun11 sec04 sub04:
   ^     - find if this is an paired primer
   ^   o fun11 sec04 sub05:
   ^     - get primers starting coordinate
   ^   o fun11 sec04 sub06:
   ^     - get primers ending coordinate
   ^   o fun11 sec04 sub07:
   ^     - Move to second primer index and restart loop
   ^   o fun11 sec04 sub08:
   ^     - Get the next line
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun11 Sec04 Sub01:
   *   - Get past header and start primer read loop
   \*****************************************************/

   tmpStr =
      (schar *)
      fgets(
         (char *) buffStr,
         lenBuffUS,
         primFILE
      );

   numPrimSI = 0;
   siLine = 1;

   while(fgets((char *) buffStr, lenBuffUS, primFILE))
   { /*Loop: Read in each primer coordinate*/
      ++siLine;
      tmpStr = buffStr;

      /**************************************************\
      * Fun11 Sec04 Sub02:
      *   - get past reference name
      \**************************************************/

      while(*tmpStr > 32)
        ++tmpStr; /*Move past the reference name*/

      if(*tmpStr == '\n')
         goto fileErr_fun11_sec05_sub03;

      if(*tmpStr == '\0')
         goto fileErr_fun11_sec05_sub03;

      while(*tmpStr < 33)
      { /*Loop: Find the start of the next entry*/
         ++tmpStr;

         if(*tmpStr == '\0')
            goto fileErr_fun11_sec05_sub03;
      } /*Loop: Find the start of the next entry*/

      /**************************************************\
      * Fun11 Sec04 Sub03:
      *   - get past primer name
      \**************************************************/

      while(*tmpStr > 32)
        ++tmpStr; /*Move past the primer name*/

      if(*tmpStr == '\n')
         goto fileErr_fun11_sec05_sub03;

      if(*tmpStr == '\0')
         goto fileErr_fun11_sec05_sub03;

      while(*tmpStr < 33)
      { /*Loop: Find the start of the next entry*/
         ++tmpStr;

         if(*tmpStr == '\0')
            goto fileErr_fun11_sec05_sub03;
      } /*Loop: Find the start of the next entry*/

      /**************************************************\
      * Fun11 Sec04 Sub04:
      *   - find if this is an paired primer
      \**************************************************/

      (*flagAryUI)[numPrimSI] = 0;

      if(*tmpStr == '1')
         {setPair_maskPrim((*flagAryUI)[numPrimSI], 1);}
      else if( ((*tmpStr) & ~32) == 'P')
         {setPair_maskPrim((*flagAryUI)[numPrimSI], 1);}
      else if( ((*tmpStr) & ~32) == 'T')
         {setPair_maskPrim((*flagAryUI)[numPrimSI], 1);}
      else if( *tmpStr == '+')
         {setPair_maskPrim((*flagAryUI)[numPrimSI], 1);}
      else
         {setPair_maskPrim((*flagAryUI)[numPrimSI], 0);}

      /*Should always be an forward primer*/
      setDir_maskPrim((*flagAryUI)[numPrimSI], 0);

      while(*tmpStr > 32)
        ++tmpStr; /*Move to end of entry*/

      if(*tmpStr == '\n')
         goto fileErr_fun11_sec05_sub03;

      if(*tmpStr == '\0')
         goto fileErr_fun11_sec05_sub03;

      while(*tmpStr < 33)
      { /*Loop: Find the start of the next entry*/
         ++tmpStr;

         if(*tmpStr == '\0')
            goto fileErr_fun11_sec05_sub03;
      } /*Loop: Find the start of the next entry*/

      /**************************************************\
      * Fun11 Sec04 Sub05:
      *   - get primers starting coordinate
      \**************************************************/

      getPrimInfo_fun11_sec04_sub05:;

      if(tmpStr[0] != 'N')
      { /*If: have a forward primer*/
         tmpStr +=
            strToUI_base10str(
               tmpStr,
               &(*startAryUI)[numPrimSI]
            );

         if(*tmpStr > 32)
            goto fileErr_fun11_sec05_sub03; /*not entry*/
            
         --(*startAryUI)[numPrimSI]; /*to index 0*/

         if(*tmpStr == '\n')
            goto fileErr_fun11_sec05_sub03;

         if(*tmpStr == '\n')
            goto fileErr_fun11_sec05_sub03;

         while(*tmpStr < 33)
         { /*Loop: Find the start of the next entry*/
            ++tmpStr;

            if(*tmpStr == '\0')
               goto fileErr_fun11_sec05_sub03;
         } /*Loop: Find the start of the next entry*/
      } /*If: have a forward primer*/

      else
      { /*Else: no forward primer*/
         while(*tmpStr++ > 33) ; /*move to next entry*/

         if(*(tmpStr - 1) == '\0')
            goto fileErr_fun11_sec05_sub03;

         while(*tmpStr++ < 33)
         { /*Loop: move to next entry*/
            if(*(tmpStr - 1) == '\0')
               goto fileErr_fun11_sec05_sub03;
         } /*Loop: move to next entry*/

         /*move past ending coordiante*/
         while(*tmpStr++ > 33) ;

         if(*(tmpStr - 1) == '\0')
            goto fileErr_fun11_sec05_sub03;

         while(*tmpStr < 33)
         { /*Loop: Find the start of the next entry*/
            ++tmpStr;

            if(*tmpStr == '\0')
               goto fileErr_fun11_sec05_sub03;
         } /*Loop: Find the start of the next entry*/

         goto skipFor_fun11_sec04_sub07;
      } /*Else: no forward primer*/

      /**************************************************\
      * Fun11 Sec04 Sub06:
      *   - get primers ending coordinate
      \**************************************************/

      tmpStr +=
         strToUI_base10str(
           tmpStr,
           &(*endAryUI)[numPrimSI]
         );

      if(*tmpStr > 32)
         goto fileErr_fun11_sec05_sub03; /*Not entry*/
         
      --(*endAryUI)[numPrimSI]; /*convert to index 0*/

      if(
           (*startAryUI)[numPrimSI]
         > (*endAryUI)[numPrimSI]
      ){ /*If: I need to swap start and end*/
         (*startAryUI)[numPrimSI] ^=
            (*endAryUI)[numPrimSI];

         (*endAryUI)[numPrimSI] ^=
            (*startAryUI)[numPrimSI];

         (*startAryUI)[numPrimSI] ^=
            (*endAryUI)[numPrimSI];
      } /*If: I need to swap start and end*/

      while(*tmpStr < 33)
      { /*Loop: move to next entry*/
         ++tmpStr;

         if(*(tmpStr) == '\0')
         { /*If: at end of line*/
            if(revBl)
               goto nextLine_fun11_sec04_sub08;
               /*have both primers coordinates*/
             else
                goto fileErr_fun11_sec05_sub03;
         } /*If: at end of line*/
      } /*Loop: move to next entry*/

      /**************************************************\
      * Fun11 Sec04 Sub07:
      *   - Move to second primer index and restart loop
      \**************************************************/

      skipFor_fun11_sec04_sub07:;

      if(revBl)
         goto fileErr_fun11_sec05_sub03;
         /*If there was more than two primers*/

      if(*tmpStr == 'N')
         goto nextLine_fun11_sec04_sub08;
         /*only a forward primer*/

      ++numPrimSI;

      /*At this point this is not the end of the line*/
      setMateIndex_maskPrim(
         (*flagAryUI)[numPrimSI - 1],
         numPrimSI
      );

      setPair_maskPrim(
         (*flagAryUI)[numPrimSI],
         getPair_maskPrim((*flagAryUI)[numPrimSI - 1])
      );

      setDir_maskPrim(
         (*flagAryUI)[numPrimSI],
         1
      );

      setMateIndex_maskPrim(
         (*flagAryUI)[numPrimSI],
         numPrimSI - 1
      );

      revBl = 1;

      goto getPrimInfo_fun11_sec04_sub05;

      /**************************************************\
      * Fun11 Sec04 Sub08:
      *   - Get the next line
      \**************************************************/

      nextLine_fun11_sec04_sub08:;

      revBl = 0;
      ++numPrimSI;
   } /*Loop: Read in each primer coordinate*/
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec05:
   ^   - Clean up and return
   ^   o fun11 sec05 sub01:
   ^     - Clean up after no errors
   ^   o fun11 sec05 sub02:
   ^     - Set up the memory error return
   ^   o fun11 sec05 sub03:
   ^     - Set up the file error return
   ^   o fun11 sec05 sub04:
   ^     - Clean up after an error
   \*****************************************************/

   /*****************************************************\
   * Fun11 Sec05 Sub01:
   *   - Clean up after no errors
   \*****************************************************/

   if(primFILE)
      fclose(primFILE);

   primFILE = 0;

   sortCoords_maskPrim(
      *startAryUI,
      *endAryUI,
      *flagAryUI,
      numPrimSI
   ); /*Sort the indexes*/ 

   *errSL = 0;
   goto cleanUp_fun11_sec05_sub05;

   /*****************************************************\
   * Fun11 Sec05 Sub02:
   *   - Set up the memory error return
   \*****************************************************/

   memErr_fun11_sec05_sub02:;
   *errSL = def_memErr_maskPrim;
   numPrimSI = 0;
   goto cleanUp_fun11_sec05_sub05;

   /*****************************************************\
   * Fun11 Sec05 Sub03:
   *   - Set up the file error return
   \*****************************************************/

   fileErr_fun11_sec05_sub03:;
   *errSL = (siLine << 8) | def_fileErr_maskPrim;
   numPrimSI = 0;
   goto cleanUp_fun11_sec05_sub05;

   /*****************************************************\
   * Fun11 Sec05 Sub04:
   *   - Set up to return empty file error
   \*****************************************************/

   emptyFileErr_fun11_sec05_sub04:;
   *errSL = def_emptyFileErr_maskPrim;
   numPrimSI = 0;
   goto cleanUp_fun11_sec05_sub05;

   /*****************************************************\
   * Fun11 Sec05 Sub05:
   *   - Clean up after an error
   \*****************************************************/

   cleanUp_fun11_sec05_sub05:;

   if(primFILE)
      fclose(primFILE);

   primFILE = 0;

   return numPrimSI;
} /*getCoords_maskPrim*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconveint / not possible, this code is under the
:   MIT license
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
