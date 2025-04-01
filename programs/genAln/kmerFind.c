/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' kmerFind SOF: Start Of File
'   - holds functions for kmer detection of primers
'   o header:
'     - included libraries
'   o .h st01: tblST_kmerFind
'     - holds the kmer tables for detecting spoligytpes
'   o .h st02: refST_kmerFind
'     - holds the kmer pattern for the reference
'   o fun01: blank_tblST_kmerFind
'     - blanks all stored values in an tblST_kmerFind
'   o fun02: qckBlank_tblST_kmerFind
'     - does a quick blank on all stored values in a
'       tblST_kmerFind struct (only stored kmers blanked)
'   o fun03: init_tblST_kmerFind
'     - initializes an tblST_kmerFind structure
'   o fun04: setup_tblST_kmerFind
'     - allocates memory for variables in a tblST_kmerFind
'       structure
'   o fun05: freeStack_tblST_kmerFind
'     - frees all variables in an tblST_kmerFind structure
'   o fun06: freeHeap_tblST_kmerFind
'     - frees and tblST_kmerFind structure
'   o fun07: blank_refST_kmerFind
'     - sets the counting and kmer values in an
'       refST_kmerFind to 0 or def_noKmer_kmerFind
'   o fun08: init_refST_kmerFind
'     - initializes an refST_kmerFind structure
'   o fun09: setup_refST_kmerFind
'     - allocates memory for a refST_kmerFind structure
'   o fun10: freeStack_refST_kmerFind
'     - frees the variables in an refST_kmerFind structure
'   o fun11: freeHeap_refST_kmerFind
'     - frees an refST_kmerFind structure
'   o fun12: freeHeapAry_refST_kmerFind
'     - frees an array of refST_kmerFind structure
'   o fun13: uiTwinShort_kmerFind
'     - sorts a unsigned int array from least to greatest
'       and keeps it linked to a second unsigned int array
'   o fun14: addSeqToRefST_kmerFInd
'     - adds a sequence to a refST_kmerFind structure
'   o fun15: prep_tblST_kmerFind
'     - sets up an tblST_kmerFind structure for primer
'       searching
'   o fun16: tsvToAry_refST_kmerFind
'     - makes an array of refST_kmerFind structures from a
'       tsv file
'   o fun17: faToAry_refST_kmerFind
'     - makes an array of refST_kmerFind structures
'   o fun18: nextSeqChunk_tblST_kmerFind
'     - adds a new set of kmers from an sequence to an
'       tblST_kmerFind structure
'   o fun19: forCntMatchs_kmerFind
'     - finds the number of kmers that are in both the
'       kmer table (query) and the pattern (reference)
'   o fun20: revCntMatchs_kmerFind
'     - finds the number of kmers that are shared in the
'       kmer table (query) and the reverse pattern
'       (reference)
'   o fun21: matchCheck_kmerFind
'     - tells if the  match meets the min requirements to
'       do an alignment or not
'   o fun22: findRefInChunk_kmerFind
'     - does an kmer check and alings an single sequence
'       in an refST_kmerFind structure to see if there is
'       an match
'   o fun24: waterFindPrims_kmerFind
'     - finds primers in an sequence (from fastx file)
'       using a slower, but more percise waterman
'   o fun24: fxFindPrims_kmerFind
'     - finds spoligotype spacers in an sequence (from
'       fastx file) using an faster kmer search followed
'       by an slower waterman to finalize alignments
'   o fun25: phit_kmerFind
'     - prints out the primer hits for a sequence
'   o fun26: pHeaderHit_kmerFind
'      - prints header for phit_kmerFind (fun25)
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

#include "kmerFind.h"

#include <stdio.h>

#include "../genLib/ulCp.h"
#include "../genLib/shellSort.h"

#include "../genBio/seqST.h"

#include "alnSet.h"
#include "memwater.h"

/*.h files only (no .c files*/
#include "../genLib/genMath.h" /*only using .h macros*/

#include "alnDefs.h"
#include "../genBio/kmerBit.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../genLib/base10Str.h"
!   - .c  #include "../genAln/indexToCoord.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: blank_tblST_kmerFind
|   - blanks all stored values in an tblST_kmerFind
|     structure
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to blank
|   - blankSeqBl:
|     o 1: blank the seqST (sequence) in tblSTPtr
|     o 0: do not blank the seqST in tblSTPtr
| Output:
|   - Modifies:
|     o tblSI in tblSTPtr to be full of zeros
|     o kmerArySI to be full of -1's and to end in -2
|     o lenLastKmer to be 0
|     o seqPosUL to be 0
|     o if requested seqST with blaknSeqST from
|       ../memwater/seqST.h
\-------------------------------------------------------*/
void
blank_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   signed char blankSeqBl
){
   unsigned int uiKmer = 0;

   if(! tblSTPtr)
      return;

   tblSTPtr->lastKmerUL = 0;
   tblSTPtr->lenLastKmerUL = 0;
   tblSTPtr->seqPosUL = 0;

   if(! tblSTPtr->kmerArySI)
      ;
   else if(tblSTPtr->tblSI)
   { /*Else If: I have an kmer table to blank*/
      for(
         uiKmer = 0;
         uiKmer < tblSTPtr->lenTblUI;
         ++uiKmer
      ) tblSTPtr->tblSI[uiKmer] = 0;
   } /*Else If: I have an kmer table to blank*/

   if(tblSTPtr->kmerArySI)
   { /*If: I have an kmer array to blank*/
      for(
         uiKmer = 0;
         uiKmer < tblSTPtr->numKmerUI;
         ++uiKmer
      ) tblSTPtr->kmerArySI[uiKmer] = def_noKmer_kmerFind;

      tblSTPtr->kmerArySI[uiKmer] = def_endKmers_kmerFind;
   } /*If: I have an kmer array to blank*/

   if(blankSeqBl && tblSTPtr->seqSTPtr)
   { /*If: I am blanking the stored seqeunce*/
      blank_seqST(tblSTPtr->seqSTPtr);
   } /*If: I am blanking the stored seqeunce*/
} /*blank_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun02: qckBlank_tblST_kmerFind
|   - does a quick blank on all stored values in a
|     tblST_kmerFind struct (only stored kmers blanked)
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to blank
|   - blankSeqBl:
|     o 1: blank the seqST (sequence) in tblSTPtr
|     o 0: do not blank the seqST in tblSTPtr
| Output:
|   - Modifies:
|     o tblSI in tblSTPtr to have stored counts set to 0
|     o kmerArySI to have all stored kmers set to -1
|       * ending is set to -2
|     o lenLastKmer to be 0
|     o seqPosUL to be 0
|     o if requested seqST with blaknSeqST from
|       ../memwater/seqST.h
\-------------------------------------------------------*/
void
qckBlank_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   signed char blankSeqBl
){
   unsigned int uiKmer = 0;

   if(! tblSTPtr)
      return;

   tblSTPtr->lastKmerUL = 0;
   tblSTPtr->lenLastKmerUL = 0;
   tblSTPtr->seqPosUL = 0;

   if(! tblSTPtr->kmerArySI)
      ;
   else if(tblSTPtr->tblSI)
   { /*Else If: I have an kmer table to blank*/
      for(
         uiKmer = 0;
         uiKmer < tblSTPtr->numKmerUI;
         ++uiKmer
      ){ /*Loop: remove all known kmers in kmer table*/

         if(tblSTPtr->kmerArySI[uiKmer] >= 0)
            tblSTPtr->tblSI[
               tblSTPtr->kmerArySI[uiKmer]
             ] = 0;

         tblSTPtr->kmerArySI[uiKmer] =def_noKmer_kmerFind;

      } /*Loop: remove all known kmers in kmer table*/
   } /*Else If: I have an kmer table to blank*/

   if(blankSeqBl && tblSTPtr->seqSTPtr)
   { /*If: I am blanking the stored seqeunce*/
      blank_seqST(tblSTPtr->seqSTPtr);
   } /*If: I am blanking the stored seqeunce*/
} /*blank_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun03: init_tblST_kmerFind
|   - initializes a tblST_kmerFind structure
| Input:
|   - tblSTPtr:
|     o pointer to tblST_kmerFind structure to initialize
| Output:
|   - Modifies:
|     o all varaibles in an tblST_kmerFind find structure
|       to be 0
\-------------------------------------------------------*/
void
init_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr
){
   if(! tblSTPtr)
      return;

   tblSTPtr->seqSTPtr = 0; /*just in case have memory error*/
   tblSTPtr->lenTblUI = 0;
   tblSTPtr->tblSI = 0;

   tblSTPtr->kmerArySI = 0;
   tblSTPtr->numKmerUI = 0;
   tblSTPtr->ntInWinUI = 0;
   tblSTPtr->rmNtUI = 0;
   tblSTPtr->seqPosUL = 0;
   tblSTPtr->kmerMaskUL = 0;

   blank_tblST_kmerFind(
      tblSTPtr,
      0           /*no seqST structures*/
   ); /*blank remaing variables*/
} /*init_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun04: setup_tblST_kmerFind
|   - allocates memory for variables in a tblST_kmerFind
|     structure
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to blank
|   - lenKmerUC:
|     o length of one kmer
| Output:
|   - Modifies:
|     o sets lenKmerUC and allocates memory for arrays
|       and seqST structures
|   - Returns:
|     o 0 for no errors
|     o def_memErr_kmerFind for an memory error
\-------------------------------------------------------*/
signed char
setup_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   unsigned char lenKmerUC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun04 TOC:
   '   - allocates memory for variables in a
   '     tblST_kmerFind
   '   o fun04 sec01:
   '     - allocate memory for the table
   '   o fun04 sec02:
   '     - build the kmer mask
   '   o fun04 sec03:
   '     - initialize seqST and blank table
   '   o fun04 sec04:
   '     - return any errors
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec01:
   ^   - allocate memory for the table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned int uiKmer = 0;

   tblSTPtr->lenKmerUC = lenKmerUC;
   tblSTPtr->lenTblUI = 1;

   for(
      uiKmer = 0;
      uiKmer < lenKmerUC;
      ++uiKmer
   ) tblSTPtr->lenTblUI <<= 2; /*multiply by 4*/

   if(tblSTPtr->tblSI)
      free(tblSTPtr->tblSI);

   tblSTPtr->tblSI = 0;

   tblSTPtr->tblSI =
      calloc(
         (tblSTPtr->lenTblUI + 1),
         sizeof(signed int)
      );

   if(! tblSTPtr->tblSI)
      goto memErr_fun04_sec03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec02:
   ^   - build the kmer mask
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   tblSTPtr->kmerMaskUL = mkMask_kmerBit(lenKmerUC);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec03:
   ^   - initialize seqST and blank table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(tblSTPtr->seqSTPtr)
     freeHeap_seqST(tblSTPtr->seqSTPtr);

   tblSTPtr->seqSTPtr = 0;

   tblSTPtr->seqSTPtr = malloc(sizeof(struct seqST));

   if(! tblSTPtr->seqSTPtr)
      goto memErr_fun04_sec03;

   init_seqST(tblSTPtr->seqSTPtr);

   blank_tblST_kmerFind(
      tblSTPtr,
      0 /*sequence structure is already blanked*/
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec04:
   ^   - return any errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun04_sec03:;

   freeStack_tblST_kmerFind(tblSTPtr);

   return def_memErr_kmerFind;
} /*setup_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun05: freeStack_tblST_kmerFind
|   - frees all variables in an tblST_kmerFind structure
|     and sets other values to defaults (calls blank)
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to free
|       variables in
| Output:
|   - Frees:
|     o all allocated varialbes in tblStPtr
\-------------------------------------------------------*/
void
freeStack_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr
){
   if(! tblSTPtr)
      return;

   if(tblSTPtr->kmerArySI)
      free(tblSTPtr->kmerArySI);

   if(tblSTPtr->seqSTPtr)
      freeHeap_seqST(tblSTPtr->seqSTPtr);

   if(tblSTPtr->tblSI)
      free(tblSTPtr->tblSI);

   init_tblST_kmerFind(tblSTPtr);
} /*freeStack_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun06: freeHeap_tblST_kmerFind
|   - frees and tblST_kmerFind structure
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to free
| Output:
|   - Frees:
|     o tblSTPtr (does not set to null)
\-------------------------------------------------------*/
void
freeHeap_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr
){
   if(! tblSTPtr)
      return

   freeStack_tblST_kmerFind(tblSTPtr);
   free(tblSTPtr);
} /*freeHeap_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun07: blank_refST_kmerFind
|   - sets the counting and kmer values in an
|     refST_kmerFind to 0 or def_noKmer_kmerFind
| Input:
|   - refSTPtr:
|     o pointer to an refST_kmerFind structure to blank
|       variables in
| Output:
|   - Modifies:
|     o forKmerArySI to be full of -1's
|     o forRefAryUS to be full of 0's
|     o revKmerArySI to be full of -1's
|     o revRefAryUS to be full of 0's
|     o forSeqST to be blanked
|     o revSeqST to be blanked
\-------------------------------------------------------*/
void
blank_refST_kmerFind(
   struct refST_kmerFind *refSTPtr
){
   signed long slKmer = 0;

   if(! refSTPtr)
      return;

   refSTPtr->mateSI = -1;

   if(refSTPtr->forKmerArySI)
   { /*If: I have kmer arrays allocated*/
      for(
         slKmer = 0;
         slKmer < refSTPtr->lenRepSI;
         ++slKmer
      ){ /*Loop: blank kmer arrays*/
         refSTPtr->forKmerArySI[slKmer] =
            def_noKmer_kmerFind;

         refSTPtr->forRepAryUI[slKmer] = 0;

         refSTPtr->revKmerArySI[slKmer] =
            def_noKmer_kmerFind;

         refSTPtr->revRepAryUI[slKmer] = 0;
      } /*Loop: blank kmer arrays*/

      refSTPtr->forKmerArySI[slKmer] =
         def_endKmers_kmerFind;

      refSTPtr->revKmerArySI[slKmer] =
         def_endKmers_kmerFind;
   } /*If: I have kmer arrays allocated*/

   if(refSTPtr->forSeqST)
      blank_seqST(refSTPtr->forSeqST);

   if(refSTPtr->revSeqST)
      blank_seqST(refSTPtr->revSeqST);
} /*blank_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun08: init_refST_kmerFind
|   - initializes an refST_kmerFind structure
| Input:
|   - refSTPtr:
|     o pointer to refST_kmerFind structure to initialize
| Output:
|   - Modifies:
|     o refSTPtr to be initialized
\-------------------------------------------------------*/
void
init_refST_kmerFind(
   struct refST_kmerFind *refSTPtr
){
   if(! refSTPtr)
      return;

   refSTPtr->lenKmerUC = 0;
   refSTPtr->minKmersUI = 0;
   refSTPtr->maxForScoreF = 0;
   refSTPtr->maxRevScoreF = 0;
   refSTPtr->lenRepSI = 0; /*not changed often*/
   refSTPtr->sizeRepSI = 0;
   refSTPtr->sizeKmerSI = 0;

   refSTPtr->forKmerArySI = 0;
   refSTPtr->forRepAryUI = 0;

   refSTPtr->revKmerArySI = 0;
   refSTPtr->revRepAryUI = 0;

   refSTPtr->forSeqST = 0;
   refSTPtr->revSeqST = 0;

   /*not needed, but for future additions*/
   blank_refST_kmerFind(refSTPtr);
} /*init_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun09: setup_refST_kmerFind
|   - allocates memory for a refST_kmerFind structure
| Input:
|   o refSTPtr,
|     - pointer to a refST_kmerFind structure to setup
|   o lenKmerUC:
|     - length of one kmer
| Output:
|   - Modifies:
|     o refSTPtr to have memory allocated to variables and
|       kmer size
|   - Returns:
|     o 0 if no errors
|     o def_memErr_kmerFind for memory errors
\-------------------------------------------------------*/
signed char
setup_refST_kmerFind(
   struct refST_kmerFind *refSTPtr,
   unsigned char lenKmerUC
){
   refSTPtr->lenKmerUC = lenKmerUC;


   /*set up forward seqST*/

   if(refSTPtr->forSeqST)
      freeHeap_seqST(refSTPtr->forSeqST);
   refSTPtr->forSeqST = 0;

   refSTPtr->forSeqST = malloc(sizeof(struct seqST));

   if(! refSTPtr->forSeqST)
      goto memErr_fun0x;
   init_seqST(refSTPtr->forSeqST);


   /*set up reverse seqST*/

   if(refSTPtr->revSeqST)
      freeHeap_seqST(refSTPtr->revSeqST);
   refSTPtr->revSeqST = 0;

   refSTPtr->revSeqST = malloc(sizeof(struct seqST));

   if(! refSTPtr->revSeqST)
      goto memErr_fun0x;
   init_seqST(refSTPtr->revSeqST);


   /*finished*/

   return 0;

   memErr_fun0x:;
      return def_memErr_kmerFind;
} /*setup_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun10: freeStack_refST_kmerFind
|   - frees the variables in an refST_kmerFind structure
| Input:
|   - refSTPtr:
|     o pointer to refST_kmerFind structure to free
|       varialbes in
| Output:
|   - Frees:
|     o all pointers in refSTPtr
|   - Sets
|     o pointers in refSTPtr to 0
\-------------------------------------------------------*/
void
freeStack_refST_kmerFind(
   struct refST_kmerFind *refSTPtr
){
   if(! refSTPtr)
      return;

   if(refSTPtr->forKmerArySI)
      free(refSTPtr->forKmerArySI);

   if(refSTPtr->forRepAryUI)
      free(refSTPtr->forRepAryUI);

   if(refSTPtr->revKmerArySI)
      free(refSTPtr->revKmerArySI);

   if(refSTPtr->revRepAryUI)
      free(refSTPtr->revRepAryUI);

   if(refSTPtr->forSeqST)
      freeHeap_seqST(refSTPtr->forSeqST);

   if(refSTPtr->revSeqST)
      freeHeap_seqST(refSTPtr->revSeqST);

   init_refST_kmerFind(refSTPtr);
} /*freeStack_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun11: freeHeap_refST_kmerFind
|   - frees an refST_kmerFind structure
| Input:
|   - refSTPtr:
|     o pointer to refST_kmerFind structure to free
| Output:
|   - Frees:
|     o refSTPtr
\-------------------------------------------------------*/
void
freeHeap_refST_kmerFind(
   struct refST_kmerFind *refSTPtr
){
   if(! refSTPtr)
      return;

   freeStack_refST_kmerFind(refSTPtr);
   free(refSTPtr);
} /*freeStack_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun12: freeHeapAry_refST_kmerFind
|   - frees an array of refST_kmerFind structure
| Input:
|   - refSTAry:
|     o pointer to refST_kmerFind structure array to free
|   - lenArySI:
|     o number of refST_kmerFind structures in refSTPtr
| Output:
|   - Frees:
|     o all refST_kmerFind structures in refSTAry
|     o refSTAry (does not set to 0)
\-------------------------------------------------------*/
void
freeHeapAry_refST_kmerFind(
   struct refST_kmerFind *refSTAry,
   signed int lenArySI
){
   if(! refSTAry)
      return;

   --lenArySI; /*convert index 1 to index 0*/

   while(lenArySI > -1)
   { /*Loop: free structres in array*/
      freeStack_refST_kmerFind(&refSTAry[lenArySI]);
      --lenArySI;
   } /*Loop: free structres in array*/

   free(refSTAry);
} /*freeHeapAry_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun13: uiTwinShort_kmerFind
|   - sorts a unsigned int array from least to greatest
|     and keeps it linked to a second unsigned int array
| Input:
|   - uiAry:
|     o array to sort
|   - uiSecAry:
|     o second array to keep in order with uiAry
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o uiAry to be sorted form least to greatest
|     o uiSecAry to be sorted by uiAry
\-------------------------------------------------------*/
void
uiTwinSort_kmerFind(
   unsigned int *uiAry,
   unsigned int *uiSecAry,
   unsigned long startUL,
   unsigned long endUL
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun13 TOC:
   '   - sorts signed int array from least to greatest
   '   o fun13 sec01:
   '     - variable declerations
   '   o fun13 sec02:
   '     - find the number of rounds to sort for
   '   o fun13 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of elements to sort*/
   unsigned long numElmUL = endUL - startUL;

   /*Number of sorting rounds*/
   unsigned long subUL = 0;
   unsigned long nextUL = 0;
   unsigned long lastUL = 0;
   unsigned long onUL = 0;

   /*Variables to incurment loops*/
   unsigned long uiIndex = 0;
   unsigned long uiElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numElmUL <= 1)
      return; /*nothing to sort*/

   /*Recursion formuia: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subUL = 1; /*Initialzie first array*/

   while(subUL < numElmUL - 1)
      subUL = (3 * subUL) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec03:
   ^   - sort arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(subUL > 0)
   { /*Loop: all rounds*/
      for(
         uiIndex = 0;
         uiIndex <= subUL;
         ++uiIndex
      ){ /*Loop: though sub array*/
         uiElm = uiIndex;

         for(
            uiElm = uiIndex;
            uiElm + subUL <= endUL;
            uiElm += subUL
         ){ /*Loop: swap elements in subarray*/
            nextUL = uiElm + subUL;

            if(uiAry[uiElm] > uiAry[nextUL])
            { /*If I need to swap an element*/
               uiAry[uiElm] ^= uiAry[nextUL];
               uiAry[nextUL] ^= uiAry[uiElm];
               uiAry[uiElm] ^= uiAry[nextUL];


               uiSecAry[uiElm] ^= uiSecAry[nextUL];
               uiSecAry[nextUL] ^= uiSecAry[uiElm];
               uiSecAry[uiElm] ^= uiSecAry[nextUL];

               lastUL = uiElm;
               onUL = uiElm;

               while(lastUL >= subUL)
               { /*loop: move swapped element back*/
                  lastUL -= subUL;

                  if(uiAry[onUL] > uiAry[lastUL])
                     break; /*Positioned the element*/
   
                  uiAry[onUL] ^= uiAry[lastUL];
                  uiAry[lastUL] ^= uiAry[onUL];
                  uiAry[onUL] ^= uiAry[lastUL];

                  uiSecAry[onUL] ^= uiSecAry[lastUL];
                  uiSecAry[lastUL] ^= uiSecAry[onUL];
                  uiSecAry[onUL] ^= uiSecAry[lastUL];

                  onUL = lastUL;
               } /*Loop: move swapped element back*/
            } /*If I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subUL = (subUL - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*uiTwinSort_kmerFind*/

/*-------------------------------------------------------\
| Fun14: addSeqToRefST_kmerFInd
|   - adds a sequence to a refST_kmerFind structure
| Input:
|   - tblSTPtr:
|     o pointer to a tblST_kmerFind structure with
|       settings, such as the kmer length, mask, and
|       maximum number of kmers
|   - refSTPtr:
|     o pionter to the refST_kmerFind structure to add the
|       sequence to
|   - seqSTPtr:
|     o pointer to the seqST to get the sequence from
|     o if 0, defaults to refSTPtr->forSeqST
|   - minPercKmersF:
|     o float with minimum percentage of kmers to start
|       considering an window supports an spacer
|   - longestSeqUI:
|     o length of the longest sequence in a refSTPtr
|       structure
|   - alnSetPtr:
|     o pointer to alnSet structure with score matrix
|       (for finding max possible score)
| Output:
|   - Returns:
|     o 0 for memory error
|     o length of longest sequence in an refST_kmerFind
|       structure
\-------------------------------------------------------*/
unsigned int
addSeqToRefST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   struct refST_kmerFind *refSTPtr,
   struct seqST *seqSTPtr,
   float minPercKmersF,
   unsigned int longestSeqUI,
   struct alnSet *alnSetPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun14 TOC:
   '   - add a sequence to a refST_kmerFind structure
   '   o fun14 sec01:
   '     - variable declerations
   '   o fun14 sec02:
   '     - copy sequence and see if is longest sequence
   '   o fun14 sec03:
   '     - find minimum number of kmers (for kmer search)
   '   o fun14 sec04:
   '     - add kmers to the kmer arrays
   '   o fun14 sec05:
   '     - move empty (not in sequence) kmers to the end
   '   o fun14 sec06:
   '     - return the longest sequence or 0 for errors
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned char errUC = 0;
   unsigned char ntUC = 0;
   signed char qckBlankBl = 0; /*1: if did quick blank*/

   signed int siKmer = 0;
   signed int siSeq = 0;

   /*for building kmers*/
   unsigned long forKmerUL = 0;
   unsigned long lenForKmerUL = 0;

   /*minimum length needed for kmer array*/
   signed int minLenSI = 0;

   unsigned long revKmerUL = 0;
   unsigned long lenRevKmerUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec02:
   ^   - copy sequence and see if is longest sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(seqSTPtr)
   { /*If: given a sequence structure*/
      errUC =
         cp_seqST(
            refSTPtr->forSeqST,
            seqSTPtr
         ); /*copy forward sequence*/

      if(errUC)
         goto memErr_fun14_sec06;
   } /*If: given a sequence structure*/

   refSTPtr->forSeqST->endAlnUL =
      refSTPtr->forSeqST->lenSeqUL - 1;

   refSTPtr->forSeqST->offsetUL = 0;

   longestSeqUI =
      max_genMath(
         longestSeqUI,
         (unsigned int) refSTPtr->forSeqST->lenSeqUL
      ); /*find the length of the longest primer*/

   errUC =
      cp_seqST(
         refSTPtr->revSeqST,
         refSTPtr->forSeqST
      ); /*copy the reverse complement sequence*/

   refSTPtr->forSeqST->endAlnUL =
      refSTPtr->forSeqST->lenSeqUL - 1;

   refSTPtr->forSeqST->offsetUL = 0;

   if(errUC)
      goto memErr_fun14_sec06;

   if(refSTPtr->revSeqST->seqStr[0] > 31)
   { /*If: human (a,t,g,c) format*/
      revComp_seqST(refSTPtr->revSeqST);

      /*convert sequences to the correct format*/
      seqToIndex_alnSet(refSTPtr->forSeqST->seqStr);
      seqToIndex_alnSet(refSTPtr->revSeqST->seqStr);
   } /*If: human (a,t,g,c) format*/
   
   else
   { /*Else: sequence is in lookup index format*/
      revCmpIndex_alnSet(
         refSTPtr->revSeqST->seqStr,  /*sequence*/
         refSTPtr->revSeqST->qStr,    /*keep in sync*/
         refSTPtr->revSeqST->lenSeqUL /*sequence length*/
      );
   } /*Else: sequence is in lookup index format*/

   /*I am merging duplicates, so I never expect
   `   more then the maxiumum number of possible kmers
   */
   refSTPtr->lenRepSI = tblSTPtr->lenTblUI;
   refSTPtr->lenKmerUC = tblSTPtr->lenKmerUC;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec03:
   ^   - set up the kmer arrays
   ^   o fun14 sec03 sub01:
   ^     -  find minimum number of kmers (for kmer search)
   ^   o fun14 sec03 sub02:
   ^     -  run quick blank (if possible)
   ^   o fun14 sec03 sub03:
   ^     - allocate memory kmer arrays
   ^   o fun14 sec03 sub04:
   ^     - allocate memory for replicate (count) arrays
   ^   o fun14 sec03 sub05:
   ^     - initialize kmer counts
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun14 Sec03 Sub01:
   *   -  find minimum number of kmers (for kmer search)
   \*****************************************************/

   refSTPtr->minKmersUI = 
      (unsigned int)
      refSTPtr->forSeqST->lenSeqUL;

   refSTPtr->minKmersUI -= refSTPtr->lenKmerUC;
   ++refSTPtr->minKmersUI; /*total kmers*/
   refSTPtr->minKmersUI *= minPercKmersF;

   minLenSI =
      min_genMath(
         refSTPtr->lenRepSI,
         (signed int) refSTPtr->forSeqST->lenSeqUL
      ); /*find shortest length for kmer array*/

   /*****************************************************\
   * Fun14 Sec03 Sub02:
   *   -  run quick blank if possible
   \*****************************************************/

   siSeq =
      (signed int)
      max_genMath(
         refSTPtr->lenForKmerSI,
         refSTPtr->lenRevKmerSI
      ); /*finding max array length*/

   qckBlankBl = 0;

   if(! refSTPtr->forKmerArySI)
      ; /*no forward kmer list*/
   else if(! refSTPtr->revKmerArySI)
      ; /*no reverse kmer list*/
   else if(! refSTPtr->forRepAryUI)
      ; /*no forward table*/
   else if(! refSTPtr->revRepAryUI)
      ; /*no reverse table*/
   else if(refSTPtr->sizeRepSI <= refSTPtr->lenRepSI)
      ; /*table will be resized*/
   else if(refSTPtr->lenRepSI < siSeq)
      ; /*kmer list larger than table*/
   else
   { /*Else: no table reallocation, quickly blank*/
      qckBlankBl = 1;

      for(
         siKmer = 0;
         siKmer < siSeq;
         ++siKmer
      ){ /*Loop: inititalize kmer counts*/
         if(siKmer >= refSTPtr->lenForKmerSI)
            ;
         else if(refSTPtr->forKmerArySI[siKmer] >= 0)
            refSTPtr->forRepAryUI[
               refSTPtr->forKmerArySI[siKmer]
            ] = 0;

         if(siKmer >= refSTPtr->lenRevKmerSI)
            ;
         else if(refSTPtr->revKmerArySI[siKmer] >= 0)
            refSTPtr->revRepAryUI[
               refSTPtr->revKmerArySI[siKmer]
            ] = 0;
      } /*Loop: inititalize kmer counts*/
   
      siKmer = 0;
   } /*Else: no table reallocation, quickly blank*/

   siSeq = 0;

   /*************************************************\
   * Fun14 Sec03 Sub03:
   *   - allocate memory kmer arrays
   \*************************************************/
   
   if(
         refSTPtr->forKmerArySI
      && refSTPtr->sizeKmerSI > minLenSI
   ) ; /*nothing needs to be done*/

   else
   { /*Else: need more memory*/

      if(refSTPtr->forKmerArySI)
         free(refSTPtr->forKmerArySI);
      refSTPtr->forKmerArySI = 0;

      refSTPtr->forKmerArySI =
         malloc(
              (refSTPtr->lenRepSI + 9)
            * sizeof(unsigned int)
         );

      if(! refSTPtr->forKmerArySI)
         goto memErr_fun14_sec06;


      if(refSTPtr->revKmerArySI)
         free(refSTPtr->revKmerArySI);
      refSTPtr->revKmerArySI = 0;

      refSTPtr->revKmerArySI =
         malloc(
              (refSTPtr->lenRepSI + 9)
            * sizeof(unsigned int)
         );

      if(! refSTPtr->revKmerArySI)
         goto memErr_fun14_sec06;

      refSTPtr->sizeKmerSI = minLenSI + 1;
   } /*Else: need more memory*/

   /*************************************************\
   * Fun14 Sec03 Sub04:
   *   - allocate memory for replicate (count) arrays
   \*************************************************/
   
   if(
         refSTPtr->revRepAryUI
      && refSTPtr->lenRepSI < refSTPtr->sizeRepSI
   ) ; /*nothing needs to be done*/

   else
   { /*Else: reverse kmer table need more memory*/
      if(refSTPtr->forRepAryUI)
         free(refSTPtr->forRepAryUI);
      refSTPtr->forRepAryUI = 0;

      refSTPtr->forRepAryUI =
         malloc(
              (refSTPtr->lenRepSI + 9)
            * sizeof(unsigned int)
         );

      if(! refSTPtr->forRepAryUI)
         goto memErr_fun14_sec06;


      if(refSTPtr->revRepAryUI)
         free(refSTPtr->revRepAryUI);
      refSTPtr->revRepAryUI = 0;

      refSTPtr->revRepAryUI =
         malloc(
              (refSTPtr->lenRepSI + 9)
            * sizeof(unsigned int)
         );

      if(! refSTPtr->revRepAryUI)
         goto memErr_fun14_sec06;

      refSTPtr->sizeRepSI = refSTPtr->lenRepSI + 1;
   } /*Else: reverse kmer table need more memory*/
   
   /*************************************************\
   * Fun14 Sec03 Sub04:
   *   - initialize kmer counts
   \*************************************************/

   if(! qckBlankBl)
   { /*If: could not do a quick blank*/
      for(
         siKmer = 0;
         siKmer < (signed int) refSTPtr->lenRepSI;
         ++siKmer
      ){ /*Loop: inititalize kmer counts*/
         refSTPtr->forRepAryUI[siKmer] = 0;
         refSTPtr->revRepAryUI[siKmer] = 0;
      } /*Loop: inititalize kmer counts*/
   } /*If: could not do a quick blank*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec04:
   ^   - add kmers to the kmer arrays
   ^   o fun14 sec04 sub01:
   ^     - get kmer counts
   ^   o fun14 sec04 sub02:
   ^     - get foward counts and max score
   ^   o fun14 sec04 sub03:
   ^     - get reverse counts and max score
   ^   o fun14 sec04 sub04:
   ^     - add end kmer marker to kmer list arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*************************************************\
   * Fun14 Sec04 Sub01:
   *   - get kmer counts
   \*************************************************/

   refSTPtr->lenForKmerSI = 0;
   refSTPtr->lenRevKmerSI = 0;

   for(
      siSeq = 0;
      siSeq < (signed int) refSTPtr->forSeqST->lenSeqUL;
      ++siSeq
   ){ /*Loop: copy the kmers*/

      /**********************************************\
      * Fun14 Sec04 Sub02:
      *   - get foward counts and max score
      \**********************************************/

      ntUC =
         (unsigned char)
         refSTPtr->forSeqST->seqStr[siSeq];

      refSTPtr->maxForScoreF +=
         getScore_alnSet(
            ntUC & def_ntToCode_alnSet,
            ntUC & def_ntToCode_alnSet,
            alnSetPtr
         );

      ntUC = alnNtTo_kmerBit[ ntUC ];

      forKmerUL =
         ntBitToKmer_kmerBit(
            ntUC,
            forKmerUL,
            tblSTPtr->kmerMaskUL
         );

      if(ntUC < def_anonNt_kmerBit)
      { /*If: no anymous bases or errors*/
         ++lenForKmerUL;

         if(lenForKmerUL < refSTPtr->lenKmerUC)
            ;
         else
         { /*Else: have a complete kmer*/
            if(! refSTPtr->forRepAryUI[forKmerUL])
               refSTPtr->forKmerArySI[
                   refSTPtr->lenForKmerSI++
               ] = (signed int) forKmerUL;
                  /*adding new kmer in*/

            ++refSTPtr->forRepAryUI[forKmerUL];
         } /*Else: have a complete kmer*/
      } /*If: no anymous bases or errors*/

      else
         lenForKmerUL = 0;

      /**********************************************\
      * Fun14 Sec04 Sub03:
      *   - get reverse counts and max score
      \**********************************************/

      ntUC =
         (unsigned char)
         refSTPtr->revSeqST->seqStr[siSeq];

      refSTPtr->maxRevScoreF +=
         getScore_alnSet(
            ntUC & def_ntToCode_alnSet,
            ntUC & def_ntToCode_alnSet,
            alnSetPtr
         );

      ntUC = alnNtTo_kmerBit[ ntUC ];

      revKmerUL =
         ntBitToKmer_kmerBit(
            ntUC,
            revKmerUL,
            tblSTPtr->kmerMaskUL
         );

      if(ntUC < def_anonNt_kmerBit)
      { /*If: no anymous bases or errors*/
         ++lenRevKmerUL;

         if(lenForKmerUL < refSTPtr->lenKmerUC)
            ;
         else
         { /*Else: have a complete kmer*/
            if(! refSTPtr->revRepAryUI[revKmerUL])
               refSTPtr->revKmerArySI[
                   refSTPtr->lenRevKmerSI++
               ] = (signed int) revKmerUL;
                  /*adding new kmer in*/

            ++refSTPtr->revRepAryUI[revKmerUL];
         } /*Else: have a complete kmer*/
      } /*If: no anymous bases or errors*/

      else
         lenRevKmerUL = 0;
   } /*Loop: copy the kmers*/

   /*****************************************************\
   * Fun14 Sec04 Sub04:
   *  - add end kmer marker to kmer list arrays
   \*****************************************************/

   refSTPtr->forKmerArySI[refSTPtr->lenForKmerSI] =
      def_endKmers_kmerFind;
   refSTPtr->revKmerArySI[refSTPtr->lenRevKmerSI] =
      def_endKmers_kmerFind;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec05:
   ^   - move empty (not in sequence) kmers to the end
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*this is likley not needed with the new kmer addition
   *  system, were I am adding kmers in one by one
   */

   /*I am converting the kmer numbers to unsigned ints
   `  so that -1's and -2's will be at the ends
   */

   /*
   ui_shellSort(
      (unsigned int *) refSTPtr->forKmerArySI,
      0,
      refSTPtr->lenRepSI - 1
   );

   ui_shellSort(
      (unsigned int *) refSTPtr->revKmerArySI,
      0,
      refSTPtr->lenRepSI - 1
   );
   */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec06:
   ^   - return the longest sequence or 0 for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return longestSeqUI;

   memErr_fun14_sec06:;
      return 0;
} /*addSeqToRefST_kmerFind*/

/*-------------------------------------------------------\
| Fun15: prep_tblST_kmerFind
|   - sets up an tblST_kmerFind structure for primer
|     searching
| Input:
|   - tblSTPtr:
|     o pointer to a tblST_kmerFind structure to set up
|   - percExtraNtInWinF:
|     o float with percentage of extra nucleotides to
|       store in one window (beyond reference length)
|   - percWinShiftF:
|     o float with percentage of bases to shift for each
|       new window in tblSTPtr
|   - longestSeqUI:
|     o longest sequence to map against. Used to find the
|       maximum window size
| Output:
|   - Returns:
|     o 0 for no errors
|     o def_memErr_kmerFind for memory errors
\-------------------------------------------------------*/
signed char
prep_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   float percExtraNtInWinF,
   float percWinShiftF,
   unsigned long longestSeqUI
){
   unsigned int lenKmerAryUI = 0;
   signed int siSeq = 0;

   /*find window size*/
   tblSTPtr->ntInWinUI = longestSeqUI;
   tblSTPtr->ntInWinUI *= percExtraNtInWinF;
   tblSTPtr->ntInWinUI += longestSeqUI;

   /*find the number of kmers to remove per window shift*/
   tblSTPtr->rmNtUI = longestSeqUI;
   tblSTPtr->rmNtUI *= percWinShiftF;

   /*find number of kmers in one window*/
   lenKmerAryUI = tblSTPtr->ntInWinUI;
   lenKmerAryUI -= tblSTPtr->lenKmerUC;
   ++lenKmerAryUI;

   /*set up the kmer window array*/
   if(tblSTPtr->numKmerUI < lenKmerAryUI)
   { /*If: the array is to small*/
      if(tblSTPtr->kmerArySI)
         free(tblSTPtr->kmerArySI); /*have an old table*/

      tblSTPtr->kmerArySI = 0;

      tblSTPtr->kmerArySI =
         malloc((lenKmerAryUI + 1) * sizeof(signed int));

      if(! tblSTPtr->kmerArySI)
         goto memErr_fun17; /*memory error*/

      tblSTPtr->numKmerUI = lenKmerAryUI;
   } /*If: the array is to small*/


   /*make sure have enough memory for table*/
   siSeq = 1;
   for(
      lenKmerAryUI = 0;
      lenKmerAryUI < tblSTPtr->lenKmerUC;
      ++lenKmerAryUI
   ) siSeq *= 4;

   if(tblSTPtr->lenTblUI < (unsigned int) siSeq)
   { /*If: need more memory*/
      if(tblSTPtr->tblSI)
         free(tblSTPtr->tblSI);
      tblSTPtr->tblSI = 0;
      tblSTPtr->tblSI =
          calloc(
             (siSeq + 1),
             sizeof(signed int)
          );
      if(! tblSTPtr->tblSI)
         goto memErr_fun17; /*memory error*/
      tblSTPtr->lenTblUI = (unsigned int) siSeq;
   } /*If: need more memory*/


   /*initialize/blank the array*/
   for(
      siSeq = 0;
      siSeq < (signed int) tblSTPtr->numKmerUI;
      ++siSeq
   ) tblSTPtr->kmerArySI[siSeq] = def_noKmer_kmerFind;
   
   tblSTPtr->kmerArySI[siSeq] = def_endKmers_kmerFind;

   return 0;

   memErr_fun17:;
      return def_memErr_kmerFind;
} /*prep_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun16: tsvToAry_refST_kmerFind
|   - makes an array of refST_kmerFind structures from a
|     tsv file
| Input:
|   - tsvFileStr:
|     o c-string with path to tsv file with sequences
|       to build an refST_kmerFind array from
|       - column 1: primer id
|       - column 2: paired (True, Yes, 1) or unpaired
|       - column 3: foward primer
|       - column 4: reverse primer
|   - lenKmerUC:
|     o length of one kmer
|   - lenArySIPtr:
|     o will hold the number of refST_kmerFind structures
|       made
|   - minPercKmersF:
|     o float with minimum percentage of kmers to start
|       considering an window supports an spacer
|   - maxNtUIPtr:
|     o pointer to unsigned int to hold the largest
|       reference sequence
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to link
|       with this reference set
|   - percExtraNtInWinF:
|     o float with percentage of extra nucleotides to
|       store in one window (beyond reference length)
|   - percWinShiftF:
|     o float with percentage of bases to shift for each
|       new window in tblSTPtr
|   - alnSetPtr:
|     o pointer to an alnSet structure with the alignment
|       settings (for finding max possible score)
|   - errSC:
|     o holds error type
| Output:
|   - Modifies:
|     o lenArySIPtr to hold the number of sequences in
|       refST_kmerFind array
|     o errSC to hold the error type
|       - 0 for no errors
|       - def_memErr_kmerFind for memory errors
|       - def_fileErr_kmerFind for file errors
|     o ntInWinUI in tblSTPtr to have the number of bases
|       in one window (2 * longest primer length)
|     o rmNtUI in tblSTPtr to have the number of bases
|       to remove with each window shift
|       (percWinShiftF * tblSTPtr->ntINWinUI)
|   - Returns:
|     o an array of refST_kmerFind structures with
|       sequences to check
|     o 0 for an error
\-------------------------------------------------------*/
refST_kmerFind *
tsvToAry_refST_kmerFind(
   signed char *tsvFileStr,
   unsigned char lenKmerUC,
   signed int *lenArySIPtr,
   float minPercKmersF,
   struct tblST_kmerFind *tblSTPtr,
   float percExtraNtInWinF,
   float percWinShiftF,
   struct alnSet *alnSetPtr,
   signed char *errSC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun16 TOC:
   '   - makes an array of refST_kmerFind structures
   '   o fun16 sec01:
   '     - variable declerations 
   '   o fun16 sec02:
   '     - check if can open input file
   '   o fun16 sec03:
   '     - get number of sequences in file
   '   o fun16 sec04:
   '     - allocate memory
   '   o fun16 sec05:
   '     - read in sequences
   '   o fun16 sec06:
   '     - update the table window size values
   '   o fun16 sec07:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec01:
   ^   - variable declerations 
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define def_lenBuff_fun16 4096
   signed char buffStr[def_lenBuff_fun16];
   signed char *tmpStr = 0;

   unsigned char pairBl = 0;
   signed int  numSeqSI = 0;
   unsigned int longestPrimUI = 0;
      /*length of longest primer*/

   /*structures to work with*/
   struct refST_kmerFind *retRefHeapAryST = 0;
   struct seqST seqStackST;

   FILE *tsvFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec02:
   ^   - initialize and check if can open input file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_seqST(&seqStackST);
   *lenArySIPtr = 0;

   tsvFILE =
      fopen(
         (char *) tsvFileStr,
         "r"
      );

   if(! tsvFILE)
      goto fileErr_fun16_sec07_sub03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec03:
   ^   - find number of lines in the file
   ^   o fun16 sec03 sub01:
   ^     - find number of lines in the file
   ^   o fun16 sec03 sub02:
   ^     - check for errors and move to start of file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun16 Sec03 Sub01:
   *   - find number of lines in the file
   \*****************************************************/

   while(fgets((char *) buffStr, def_lenBuff_fun16, tsvFILE))
      ++numSeqSI;

   --numSeqSI;     /*account for the header*/
   numSeqSI <<= 1; /*acount for two sequences per line*/

   /*****************************************************\
   * Fun16 Sec03 Sub04:
   *   - check for errors and move to start of file
   \*****************************************************/

   if(numSeqSI < 1)
      goto fileErr_fun16_sec07_sub03; /*no sequences*/

   fseek(
      tsvFILE,
      0,
      SEEK_SET
   ); /*find start of file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec04:
   ^   - allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retRefHeapAryST =
      malloc(numSeqSI * sizeof(struct refST_kmerFind));

   if(! retRefHeapAryST)
      goto memErr_fun16_sec07_sub02;

   for(
      *lenArySIPtr = 0;
      *lenArySIPtr < numSeqSI;
      ++(*lenArySIPtr)
   ){ /*Loop: initialize the new structures*/
      init_refST_kmerFind(&retRefHeapAryST[*lenArySIPtr]);

      *errSC |=
         setup_refST_kmerFind(
            &retRefHeapAryST[*lenArySIPtr],
            lenKmerUC
         );
   } /*Loop: initialize the new structures*/

   if(*errSC)
   { /*If: I had an memroy error*/
      *lenArySIPtr = numSeqSI;
      goto memErr_fun16_sec07_sub02;
   } /*If: I had an memroy error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec05:
   ^   - read in sequences from tsv file
   ^   o fun16 sec05 sub01:
   ^     - set sequence counter to zero and start loop
   ^   o fun16 sec05 sub02:
   ^     - get primer id
   ^   o fun16 sec05 sub11:
   ^     - update the table window size values
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun16 Sec05 Sub01:
   *   - set sequence counter to zero and start loop
   \*****************************************************/

   *lenArySIPtr = 0;

   /*get past header*/
   tmpStr =
      (signed char *)
      fgets(
         (char *) buffStr,
         def_lenBuff_fun16,
         tsvFILE
      );

   while(fgets((char *) buffStr, def_lenBuff_fun16, tsvFILE))
   { /*Loop: read in each sequence in the file*/

       /*************************************************\
       * Fun15 Sec05 Sub02:
       *   - get primer id
       \*************************************************/

       tmpStr = buffStr;
       seqStackST.idStr = buffStr;

       while(*tmpStr > 32)
          ++tmpStr;

       if(*tmpStr != '\t' && *tmpStr != ' ')
          goto fileErr_fun16_sec07_sub03; /*invalid line*/
       
       seqStackST.lenIdUL = tmpStr - buffStr;

       *tmpStr = '\0'; /*c-string for seqST*/
       ++tmpStr;

       /*************************************************\
       * Fun15 Sec05 Sub03:
       *   - check if primers are paired
       \*************************************************/

       while(*tmpStr < 33)
       { /*Loop: get off white space*/
          if(*tmpStr != '\t' && *tmpStr != ' ')
             goto fileErr_fun16_sec07_sub03; /*invalid*/

          ++tmpStr;
       } /*Loop: get off white space*/

       if(*tmpStr < 33)
          goto fileErr_fun16_sec07_sub03; /*invalid line*/

       pairBl = 0;
       pairBl |= (( *tmpStr & (~32) ) == 'T');
       pairBl |= (( *tmpStr & (~32) ) == 'Y');
       pairBl |= (*tmpStr == '1');

       if(pairBl)
       { /*If: the primers are paired*/
          retRefHeapAryST[*lenArySIPtr].mateSI =
             *lenArySIPtr + 1;

          retRefHeapAryST[*lenArySIPtr + 1].mateSI =
             *lenArySIPtr;
       } /*If: the primers are paired*/

       while(*tmpStr++ > 32) ;

       --tmpStr;

       if(*tmpStr != '\t' && *tmpStr != ' ')
          goto fileErr_fun16_sec07_sub03; /*invalid line*/

       /*************************************************\
       * Fun15 Sec05 Sub04:
       *   - get the foward sequence
       \*************************************************/

       while(*tmpStr < 33)
       { /*Loop: get off white space*/
          if(*tmpStr != '\t' && *tmpStr != ' ')
             goto fileErr_fun16_sec07_sub03; /*invalid*/

          ++tmpStr;
       } /*Loop: get off white space*/

       if(*tmpStr < 33)
          goto fileErr_fun16_sec07_sub03; /*invalid line*/

       seqStackST.seqStr = tmpStr;

       while(*tmpStr > 32)
          ++tmpStr;

       if(*tmpStr != '\t' && *tmpStr != ' ')
          goto fileErr_fun16_sec07_sub03; /*invalid line*/

       *tmpStr = '\0'; /*convert sequence to c-string*/
       ++tmpStr;
       
       /*check if this is a blank entry*/
       if(seqStackST.seqStr[0] == '0')
       { /*If: no foward sequence*/
          retRefHeapAryST[*lenArySIPtr].mateSI = -1;
          retRefHeapAryST[*lenArySIPtr + 1].mateSI = -1;

          goto noForSeq_fun16_sec05_sub0x;
       } /*If: no foward sequence*/

       if(
             (seqStackST.seqStr[0] & ~ 32) == 'N'
          && (seqStackST.seqStr[1] & ~ 32) == 'A'
       ){ /*If: no foward sequence (NA)*/
          retRefHeapAryST[*lenArySIPtr].mateSI = -1;
          retRefHeapAryST[*lenArySIPtr +1].mateSI = -1;

          goto noForSeq_fun16_sec05_sub0x;  
       } /*If: no foward sequence (NA)*/

       /*find sequence length*/
       seqStackST.lenSeqUL =
          tmpStr - seqStackST.seqStr - 1;
          /*-1 to account for tmpStr being one off*/

       longestPrimUI =
          addSeqToRefST_kmerFind(
             tblSTPtr, 
             &(retRefHeapAryST[*lenArySIPtr]),
             &seqStackST,
             minPercKmersF,
             longestPrimUI,
             alnSetPtr
           ); /*add the sequence to the reference array*/

       if(! longestPrimUI)
          goto memErr_fun16_sec07_sub02;

       /*************************************************\
       * Fun15 Sec05 Sub06:
       *   - get the reverse sequence
       \*************************************************/

       ++(*lenArySIPtr);

       noForSeq_fun16_sec05_sub0x:;

       while(*tmpStr < 33)
       { /*Loop: get off white space*/
          if(*tmpStr != '\t' && *tmpStr != ' ')
             goto fileErr_fun16_sec07_sub03; /*invalid*/

          ++tmpStr;
       } /*Loop: get off white space*/

       if(*tmpStr < 33)
          goto fileErr_fun16_sec07_sub03; /*invalid line*/

       seqStackST.seqStr = tmpStr;

       while(*tmpStr > 32)
          ++tmpStr;

       *tmpStr = '\0'; /*convert sequence to c-string*/
       
       /*check if this is a blank entry*/
       if(seqStackST.seqStr[0] == '0')
       { /*If: no reverse sequence*/
          retRefHeapAryST[*lenArySIPtr].mateSI = -1;
          retRefHeapAryST[*lenArySIPtr +1].mateSI = -1;

          continue; /*no reverse sequence*/
       } /*If: no reverse sequence*/

       if(
             (seqStackST.seqStr[0] & ~ 32) == 'N'
          && (seqStackST.seqStr[1] & ~ 32) == 'A'
       ){ /*If: no foward sequence (NA)*/
          retRefHeapAryST[*lenArySIPtr].mateSI = -1;
          retRefHeapAryST[*lenArySIPtr +1].mateSI = -1;

          continue; /*no reverse sequence*/
       } /*If: no reverse sequence (NA)*/

       /*find sequence length*/
       seqStackST.lenSeqUL = tmpStr - seqStackST.seqStr;
          /*I did not move tmpStr past the null*/

       longestPrimUI =
          addSeqToRefST_kmerFind(
             tblSTPtr,
             &(retRefHeapAryST[*lenArySIPtr]),
             &seqStackST,
             minPercKmersF,
             longestPrimUI,
             alnSetPtr
           ); /*add the sequence to the reference array*/

       if(! longestPrimUI)
          goto memErr_fun16_sec07_sub02;

       ++(*lenArySIPtr);
   } /*Loop: read in each sequence in the file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec06:
   ^   - update the table window size values
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *errSC =
      prep_tblST_kmerFind(
         tblSTPtr,
         percExtraNtInWinF,
         percWinShiftF,
         longestPrimUI
      );

   if(*errSC)
      goto memErr_fun16_sec07_sub02; /*memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec07:
   ^   - clean up
   ^   o fun16 sec06 sub01:
   ^     - success clean up
   ^   o fun16 sec06 sub02:
   ^     - memory error report
   ^   o fun16 sec06 sub03:
   ^     - file error report
   ^   o fun16 sec06 sub04:
   ^     - error clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun16 Sec07 Sub01:
   *   - success clean up
   \*****************************************************/

   *errSC = 0;

   fclose(tsvFILE);
   tsvFILE = 0;

   return retRefHeapAryST;

   /*****************************************************\
   * Fun16 Sec07 Sub02:
   *   - memory error report
   \*****************************************************/

   memErr_fun16_sec07_sub02:;

   *errSC = def_memErr_kmerFind;

   goto errCleanUp_fun16_sec07_sub04;

   /*****************************************************\
   * Fun16 Sec07 Sub03:
   *   - file error report
   \*****************************************************/

   fileErr_fun16_sec07_sub03:;

   *errSC = def_fileErr_kmerFind;

   goto errCleanUp_fun16_sec07_sub04;

   /*****************************************************\
   * Fun16 Sec07 Sub04:
   *   - error clean up
   \*****************************************************/

   errCleanUp_fun16_sec07_sub04:;

   if(tsvFILE)
      fclose(tsvFILE);

   tsvFILE = 0;

   freeHeapAry_refST_kmerFind(
      retRefHeapAryST,
      *lenArySIPtr
   );

   retRefHeapAryST = 0;
   *lenArySIPtr = 0;

   seqStackST.idStr = 0;
   seqStackST.seqStr = 0;

   freeStack_seqST(&seqStackST);

   return 0;
} /*tsvToAry_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun17: faToAry_refST_kmerFind
|   - makes an array of refST_kmerFind structures from a
|     fasta file
| Input:
|   - faFileStr:
|     o c-string with path to fasta file with sequences
|       to build an refST_kmerFind array from
|   - lenKmerUC:
|     o length of one kmer
|   - lenArySIPtr:
|     o will hold the number of refST_kmerFind structures
|       made
|   - minPercKmersF:
|     o float with minimum percentage of kmers to start
|       considering an window supports an spacer
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to link
|       with this reference set
|   - percExtraNtInWinF:
|     o float with percentage of extra nucleotides to
|       store in one window (beyond reference length)
|   - percWinShiftF:
|     o float with percentage of bases to shift for each
|       new window in tblSTPtr
|   - alnSetPtr:
|     o pointer to an alnSet structure with the alignment
|       settings (for finding max possible score)
|   - errSC:
|     o holds error type
| Output:
|   - Modifies:
|     o lenArySIPtr to hold the number of sequences in
|       th refST_kmerFind array
|     o errSC to hold the error type
|       - 0 for no errors
|       - def_memErr_kmerFind for memory errors
|       - def_fileErr_kmerFind for file errors
|     o ntInWinUI in tblSTPtr to have the number of bases
|       in one window (2 * longest primer length)
|     o rmNtUI in tblSTPtr to have the number of bases
|       to remove with each window shift
|       (percWinShiftF * tblSTPtr->ntINWinUI)
|   - Returns:
|     o an array of refST_kmerFind structures with
|       sequences to check
|     o 0 for an error
\-------------------------------------------------------*/
refST_kmerFind *
faToAry_refST_kmerFind(
   signed char *faFileStr,
   unsigned char lenKmerUC,
   signed int *lenArySIPtr,
   float minPercKmersF,
   struct tblST_kmerFind *tblSTPtr,
   float percExtraNtInWinF,
   float percWinShiftF,
   struct alnSet *alnSetPtr,
   signed char *errSC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun17 TOC:
   '   - makes an array of refST_kmerFind structures
   '   o fun17 sec01:
   '     - variable declerations 
   '   o fun17 sec02:
   '     - check if can open input file
   '   o fun17 sec03:
   '     - get number of sequences in file
   '   o fun17 sec04:
   '     - allocate memory
   '   o fun17 sec05:
   '     - read in sequences
   '   o fun17 sec06:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec01:
   ^   - variable declerations 
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define def_lenBuff_fun17 4096
   signed char buffStr[def_lenBuff_fun17];
   unsigned long bytesUL = 0;
   unsigned long posUL = 0;
   signed char newLineBl = 0;

   signed int  numSeqSI = 0;

   unsigned int longestPrimUI = 0;
      /*length of longest primer*/

   /*structures to work with*/
   struct refST_kmerFind *retRefHeapAryST = 0;

   FILE *faFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec02:
   ^   - initialize and check if can open input file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *lenArySIPtr = 0;

   faFILE = fopen((char *) faFileStr, "r");

   if(! faFILE)
      goto fileErr_fun17_sec06_sub03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec03:
   ^   - get number of sequences in file
   ^   o fun17 sec03 sub01:
   ^     - read in first part of file and start loop
   ^   o fun17 sec03 sub02:
   ^     - see if '>' for new line cases
   ^   o fun17 sec03 sub03:
   ^     - find next newline or end of buffer
   ^   o fun17 sec03 sub04:
   ^     - check for errors and move to start of file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun17 Sec03 Sub01:
   *   - read in first part of file and start loop
   \*****************************************************/

   bytesUL =
      fread(
         (char *) buffStr,
         sizeof(char),
         def_lenBuff_fun17 - 1,
         faFILE
      );

   buffStr[bytesUL] = '\0';
   posUL = 0;
   newLineBl = 1; /*assume I start on an new line*/

   while(bytesUL)
   { /*Loop: find number of sequences in file*/

      /**************************************************\
      * Fun17 Sec03 Sub02:
      *   - see if '>' for new line cases
      \**************************************************/

      if(newLineBl)
      { /*If: I just finished an newline*/
         while(buffStr[posUL] < 33)
         { /*Loop: Get past white space*/
            if(posUL >= bytesUL)
            { /*If: I need to read in more buffer*/
               bytesUL =
                  fread(
                     (char *) buffStr,
                     sizeof(char),
                     def_lenBuff_fun17 - 1,
                     faFILE
                  );

               buffStr[bytesUL] = '\0';
               posUL = 0;

               if(! bytesUL)
                  goto EOF_fun17_sec03_sub04;
            } /*If: I need to read in more buffer*/

            else
               ++posUL;
         } /*Loop: Get past white space*/

         if(buffStr[posUL] == '>')
            ++numSeqSI;

         newLineBl = 0;
      } /*If: I just finished an newline*/

      /**************************************************\
      * Fun17 Sec03 Sub03:
      *   - find next newline or end of buffer
      \**************************************************/

      posUL += endLine_ulCp(&buffStr[posUL]);

      if(buffStr[posUL] == '\0')
      { /*If: I need to get more file*/
         bytesUL =
            fread(
               (char *) buffStr,
               sizeof(char),
               def_lenBuff_fun17,
               faFILE
            );

         posUL = 0;
      } /*If: I need to get more file*/

      else
         newLineBl = 1;
   } /*Loop: find number of sequences in file*/

   /*****************************************************\
   * Fun17 Sec03 Sub04:
   *   - check for errors and move to start of file
   \*****************************************************/

   EOF_fun17_sec03_sub04:;

   if(*errSC)
   { /*If: I had an error of some kind*/
      if(*errSC & 64)
         goto memErr_fun17_sec06_sub02; /*memory error*/

      goto fileErr_fun17_sec06_sub03;
   } /*If: I had an error of some kind*/

   if(numSeqSI < 1)
      goto fileErr_fun17_sec06_sub03; /*no seq in file*/

   fseek(faFILE, 0, SEEK_SET);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec04:
   ^   - allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retRefHeapAryST =
      malloc(numSeqSI * sizeof(struct refST_kmerFind));

   if(! retRefHeapAryST)
      goto memErr_fun17_sec06_sub02;

   for(
      *lenArySIPtr = 0;
      *lenArySIPtr < numSeqSI;
      ++(*lenArySIPtr)
   ){ /*Loop: initialize the new structures*/
      init_refST_kmerFind(&retRefHeapAryST[*lenArySIPtr]);

      *errSC |=
         setup_refST_kmerFind(
            &retRefHeapAryST[*lenArySIPtr],
            lenKmerUC
         );
   } /*Loop: initialize the new structures*/

   if(*errSC)
   { /*If: I had an memroy error*/
      *lenArySIPtr = numSeqSI;
      goto memErr_fun17_sec06_sub02;
   } /*If: I had an memroy error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec05:
   ^   - read in sequences
   ^   o fun17 sec05 sub01:
   ^     - add sequences to the reference array
   ^   o fun17 sec05 sub02:
   ^     - deal with errors or no sequences
   ^   o fun17 sec05 sub11:
   ^     - update the table window size values
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun17 Sec05 Sub01:
   *   - add sequences to the reference array
   \*****************************************************/

   for(
      numSeqSI = 0;
      numSeqSI < *lenArySIPtr;
      ++numSeqSI
   ){ /*Loop: read in each sequence in the file*/
       *errSC = 
          getFaSeq_seqST(
            faFILE,
            retRefHeapAryST[numSeqSI].forSeqST
          );

       if(*errSC > 1)
          goto memErr_fun17_sec06_sub02; /*memory error*/

       longestPrimUI =
          addSeqToRefST_kmerFind(
             tblSTPtr,
             &(retRefHeapAryST[numSeqSI]),
             0,   /*foward sequence has sequence*/
             minPercKmersF,
             longestPrimUI,
             alnSetPtr
           ); /*add the sequence to the reference array*/

       if(! longestPrimUI)
          goto memErr_fun17_sec06_sub02;
   } /*Loop: read in each sequence in the file*/

   /*****************************************************\
   * Fun17 Sec05 Sub02:
   *   - deal with errors or no sequences
   \*****************************************************/

   if(*errSC)
   { /*If: I had an error of some kind*/
      if(*errSC & 64)
         goto memErr_fun17_sec06_sub02; /*memory error*/

      goto fileErr_fun17_sec06_sub03;
   } /*If: I had an error of some kind*/

   if(numSeqSI < 1)
      goto fileErr_fun17_sec06_sub03; /*no seq in file*/

   /*****************************************************\
   * Fun17 Sec05 Sub11:
   *   - update the table window size values
   \*****************************************************/

   *errSC =
      prep_tblST_kmerFind(
         tblSTPtr,
         percExtraNtInWinF,
         percWinShiftF,
         longestPrimUI
      );

   if(*errSC)
      goto memErr_fun17_sec06_sub02; /*memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec06:
   ^   - clean up
   ^   o fun17 sec06 sub01:
   ^     - success clean up
   ^   o fun17 sec06 sub02:
   ^     - memory error report
   ^   o fun17 sec06 sub03:
   ^     - file error report
   ^   o fun17 sec06 sub04:
   ^     - error clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun17 Sec06 Sub01:
   *   - success clean up
   \*****************************************************/

   *errSC = 0;

   fclose(faFILE);
   faFILE = 0;

   return retRefHeapAryST;

   /*****************************************************\
   * Fun17 Sec06 Sub02:
   *   - memory error report
   \*****************************************************/

   memErr_fun17_sec06_sub02:;

   *errSC = def_memErr_kmerFind;

   goto errCleanUp_fun17_sec06_sub04;

   /*****************************************************\
   * Fun17 Sec06 Sub03:
   *   - file error report
   \*****************************************************/

   fileErr_fun17_sec06_sub03:;

   *errSC = def_fileErr_kmerFind;

   goto errCleanUp_fun17_sec06_sub04;

   /*****************************************************\
   * Fun17 Sec06 Sub04:
   *   - error clean up
   \*****************************************************/

   errCleanUp_fun17_sec06_sub04:;

   if(faFILE)
      fclose(faFILE);

   faFILE = 0;

   freeHeapAry_refST_kmerFind(
      retRefHeapAryST,
      *lenArySIPtr
   );

   retRefHeapAryST = 0;
   *lenArySIPtr = 0;

   return 0;
} /*faToAry_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun18: nextSeqChunk_tblST_kmerFind
|   - adds a new set of kmers from an sequence to an
|     tblST_kmerFind structure
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to add
|       kmers to
|   - firstTimeBl:
|     o 1: first time adding sequence (blank kmer array)
|     o 0: updating the kmer window
| Output:
|   - Modifies:
|     o tblSI and seqAryUS in tblSTPtr to have the old
|       kmers (number specified by rmNtUI in tblSI)
|       remove and the new kmers added in
|       - for end of sequence it sets an index to
|         def_endKmers_kmerFind
|    o firstTimeBl:
|      o to be 0 if it is 1
|   - Returns:
|     o 0 for not end of sequence
|     o 1 for end of sequence
\-------------------------------------------------------*/
signed char
nextSeqChunk_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,/*table to add seq to*/
   signed char *firstTimeBl /*1: first set of kmers*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun18 TOC:
   '   o fun18 sec01:
   '     - variable declerations
   '   o fun18 sec02:
   '     - handle first time adding in sequence cases
   '   o fun18 sec03:
   '     - remove the old kmers from the table
   '   o fun18 sec04:
   '     - move keep kmers to start of kmer array
   '   o fun18 sec05:
   '     - add new kmers to table
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned int uiNt = 0;
   unsigned int endWindowUI = 0;
   unsigned int dupNtUI = 0;
   signed char lastWinBl = 0;

   unsigned char ntUC = 0; /*holds one nucleotide*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec02:
   ^   - handle first time adding in sequence cases
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(*firstTimeBl)
   { /*If: this is the first window*/
      *firstTimeBl = 0; /*no longer first time*/
      tblSTPtr->lastKmerUL = 0;    /*no kmers in window*/
      tblSTPtr->lenLastKmerUL = 0; /*no kmers in window*/

      endWindowUI = tblSTPtr->seqPosUL;
      endWindowUI += tblSTPtr->lenKmerUC;
      --endWindowUI; /*I want to be 1 base off the kmer*/

      while(
         tblSTPtr->seqPosUL < (unsigned long) endWindowUI
      ){ /*Loop: build the first kmer*/
         ntUC =
            (unsigned char) 
            tblSTPtr->seqSTPtr->seqStr[
               tblSTPtr->seqPosUL
            ];

         ntUC = alnNtTo_kmerBit[ ntUC ];
 
         tblSTPtr->lastKmerUL =
            ntBitToKmer_kmerBit(
               ntUC,
               tblSTPtr->lastKmerUL,
               tblSTPtr->kmerMaskUL
            );

         ++tblSTPtr->seqPosUL;

         if(ntUC < def_anonNt_kmerBit)
            ++tblSTPtr->lenLastKmerUL;
          else
            tblSTPtr->lenLastKmerUL = 0;
      } /*Loop: build the first kmer*/

      endWindowUI += tblSTPtr->ntInWinUI;
      endWindowUI -= tblSTPtr->lenKmerUC;
      ++endWindowUI; /*index 0 to index 1*/
      dupNtUI = 0;   /*first base in window*/

      goto firstKmers_fun18_sec04;
   } /*If: this is the first window*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec03:
   ^   - remove the old kmers from the table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      uiNt = 0;
      uiNt < (tblSTPtr)->rmNtUI;
      ++uiNt
   ){ /*Loop: remove discarded kmers from table*/
      tblSTPtr->lastKmerUL =
         (unsigned long) tblSTPtr->kmerArySI[uiNt];

      if(tblSTPtr->kmerArySI[uiNt] < 0)
      { /*If: I have an invalid kmer*/
         if(
               tblSTPtr->kmerArySI[uiNt]
            == def_endKmers_kmerFind
         ) break; /*was an end of kmers*/

         continue; /*was an no kmer*/
      } /*If: I have an invalid kmer*/

      --tblSTPtr->tblSI[tblSTPtr->lastKmerUL];
   } /*Loop: remove discarded kmers from table*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec04:
   ^   - move keep kmers to start of kmer array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   dupNtUI = 0;

   for(
      uiNt = (tblSTPtr)->rmNtUI;
      uiNt < (tblSTPtr)->numKmerUI;
      ++uiNt
   ){ /*Loop: copy keep kmers over*/
      tblSTPtr->kmerArySI[ dupNtUI ] =
         tblSTPtr->kmerArySI[ uiNt ];

      if(
            tblSTPtr->kmerArySI[ uiNt ]
         == def_endKmers_kmerFind
      ) break; /*finished copying*/

      ++dupNtUI;
   } /*Loop: copy keep kmers over*/

   endWindowUI = tblSTPtr->ntInWinUI;
   endWindowUI -= tblSTPtr->rmNtUI;
   endWindowUI += tblSTPtr->seqPosUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec05:
   ^   - add new kmers to table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   firstKmers_fun18_sec04:; /*first time for sequence*/

   if(endWindowUI > tblSTPtr->seqSTPtr->lenSeqUL)
   { /*If: this is the last window*/
      endWindowUI = tblSTPtr->seqSTPtr->lenSeqUL;
      lastWinBl = 1;
   } /*If: this is the last window*/

   while(tblSTPtr->seqPosUL < endWindowUI)
   { /*Loop: add kmers in*/
      ntUC =
         (unsigned char) 
         tblSTPtr->seqSTPtr->seqStr[tblSTPtr->seqPosUL];

      ntUC = alnNtTo_kmerBit[ ntUC ];

      tblSTPtr->lastKmerUL =
         ntBitToKmer_kmerBit(
            ntUC,
            tblSTPtr->lastKmerUL,
            tblSTPtr->kmerMaskUL
         );

      /*add kmer to kmer table counts if is complete*/

      if(ntUC < def_anonNt_kmerBit)
      { /*If: no anymous bases or errors*/
         ++tblSTPtr->lenLastKmerUL;

         if(
               (signed long) tblSTPtr->lenLastKmerUL
            >= tblSTPtr->lenKmerUC
         ){ /*If: had a forward kmer*/
            ++tblSTPtr->tblSI[ tblSTPtr->lastKmerUL ];

            tblSTPtr->kmerArySI[ dupNtUI ] =
               (signed int) tblSTPtr->lastKmerUL;
         } /*If: had a forward kmer*/

         else
            tblSTPtr->kmerArySI[ dupNtUI ] =
               def_noKmer_kmerFind;
      } /*If: no anymous bases or errors*/

      else
      { /*Else: kmer incomplete*/
         tblSTPtr->lenLastKmerUL = 0;

         tblSTPtr->kmerArySI[ dupNtUI ] =
            def_noKmer_kmerFind;
      } /*Else: kmer incomplete*/

      ++tblSTPtr->seqPosUL; /*move to next base in seq*/
      ++dupNtUI;              /*add the new kmer in*/
   } /*Loop: add kmers in*/

   if(lastWinBl)
   { /*If: this is the last window*/
      /*mark end of sequence*/
      tblSTPtr->kmerArySI[dupNtUI] =
         def_endKmers_kmerFind;

      return 1;
   } /*If: this is the last window*/

   return 0;
} /*nextSeqChunk_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun19: forCntMatchs_kmerFind
|   - finds the number of kmers that are in both the
|     kmer table (query) and the pattern (reference)
| Input:
|   - tblST_kmerFindPtr:
|     o pointer to an tblST_kmerFind structure with an
|       kmer table to get the number of matches from
|   - refST_kmerFindPtr:
|     o pointer to an refST_kmerFind structure with an
|       array of kmers to get the number of matches for
| Output:
|   - Returns:
|     o number of matching kmers
\-------------------------------------------------------*/
signed long
forCntMatchs_kmerFind(
   struct tblST_kmerFind *tblST_kmerFindPtr,
   struct refST_kmerFind *refST_kmerFindPtr
){
   signed long tmpSL = 0;
   signed long slKmer = 0;
   signed long numMatchSL = 0;

   for(
      slKmer = 0;
      slKmer <
         (signed long) refST_kmerFindPtr->lenForKmerSI;
      ++slKmer
   ){ /*Loop: count number of matching kmers*/
      tmpSL = refST_kmerFindPtr->forKmerArySI[ slKmer ];

      tmpSL =
         min_genMath(
            (signed int)
               refST_kmerFindPtr->forRepAryUI[tmpSL],
            tblST_kmerFindPtr->tblSI[ tmpSL ]
         );

      numMatchSL += tmpSL;
   } /*Loop: count number of matching kmers*/

   return numMatchSL;
} /*forCntMatchs_kmerFind*/ 

/*-------------------------------------------------------\
| Fun20: revCntMatchs_kmerFind
|   - finds the number of kmers that are shared in the
|     kmer table (query) and the reverse pattern
|     (reference)
| Input:
|   - tblST_kmerFindPtr:
|     o pointer to an tblST_kmerFind structure with an
|       kmer table to get the number of matches from
|   - refST_kmerFindPtr:
|     o pointer to an refST_kmerFind structure with an
|       array of kmers to get the number of matches rev
| Output:
|   - Returns:
|     o number of matching kmers
\-------------------------------------------------------*/
signed long
revCntMatchs_kmerFind(
   struct tblST_kmerFind *tblST_kmerFindPtr,
   struct refST_kmerFind *refST_kmerFindPtr
){
   signed long tmpSL = 0;
   signed long slKmer = 0;
   signed long numMatchSL = 0;

   for(
      slKmer = 0;
      slKmer <
         (signed long) refST_kmerFindPtr->lenRevKmerSI;
      ++slKmer
   ){ /*Loop: count number of matching kmers*/
      tmpSL =
         refST_kmerFindPtr->revKmerArySI[ slKmer ];

      tmpSL =
         min_genMath(
            (signed int)
               refST_kmerFindPtr->revRepAryUI[tmpSL],
            tblST_kmerFindPtr->tblSI[ tmpSL ]
         );

      numMatchSL += tmpSL;
   } /*Loop: count number of matching kmers*/

   return numMatchSL;
} /*revCntMatchs_kmerFind*/

/*-------------------------------------------------------\
| Fun21: matchCheck_kmerFind
|   - tells if the  match meets the min requirements to
|     do an alignment or not
| Input:
|   - tblST_kmerFindPtr:
|     o pointer to an tblST_kmerFind structure with the
|       chunk of query (kmer table) to check
|   - refST_kmerFindPtr:
|     o pointer to an refST_kmerFind structure with the
|       reference (primers) kmers to check
| Output:
|   - Returns:
|     o 0 (foward) or 1 (best) if there is no match
|     o 0 | 1 if match and the forward direction is best
|     o 2 | 1 if match and the reverse direction is best
\-------------------------------------------------------*/
signed char
matchCheck_kmerFind(  
   struct tblST_kmerFind *tblST_kmerFindPtr, 
   struct refST_kmerFind *refST_kmerFindPtr  
){ 
   signed char retBl = 0; 
   signed char revBestBl = 0; 
   unsigned long forKmerCntUL = 0; 
   unsigned long revKmerCntUL = 0; 
   unsigned long bestKmerCntUL = 0; 
   
   forKmerCntUL = 
      forCntMatchs_kmerFind( 
         tblST_kmerFindPtr, 
         refST_kmerFindPtr 
      ); 
   
   revKmerCntUL = 
      revCntMatchs_kmerFind( 
         tblST_kmerFindPtr, 
         refST_kmerFindPtr 
      ); 
   
   bestKmerCntUL = 
      max_genMath( 
         forKmerCntUL, 
         revKmerCntUL 
      ); 
   
   revBestBl = (revKmerCntUL == bestKmerCntUL); 
   revBestBl <<= 1; 
   
   retBl=( 
      bestKmerCntUL >= (refST_kmerFindPtr)->minKmersUI
   ); 
   
   return retBl | revBestBl;
} /*matchCheck_kmerFind*/

/*-------------------------------------------------------\
| Fun22: findRefInChunk_kmerFind
|   - does an kmer check and alings an single sequence
|     in an refST_kmerFind structure to see if there is
|     an match
| Input:
|   - tblST_kmerFindPtr:
|     o pointer to an tblST_kmerFind structure with the
|       chunk of query (kmer table) to check
|     o the stored sequence must be converted with
|       seqToIndex_alnSet from alnSetStruct.h
|   - refST_kmerFindPtr:
|     o pointer to an refST_kmerFind structure with the
|       reference (primers) kmers to check
|   - alnSetPtr:
|     o pointer to an alnSet structure with the alignment
|       settings
|   - minPerScoreF:
|     o float with minimum percent score to keep an
|       alingment
|   - scoreSL:
|     o pointer to an signed long to hold the alingment
|       score
|   - qryStartUL:
|     o pointer to an unsigned long to hold the first base
|       in the query alingment
|   - qryEndUL:
|     o pointer to an unsigned long to hold the last base
|       in the query alingment
|   - refStartUL:
|     o pointer to an unsigned long to hold the first base
|       in the reference alingment
|   - refEndUL:
|     o pointer to an unsigned long to hold the last base
|       in the reference alingment
| Output:
|   - Modifies:
|     o scoreSL
|       - 0 if no alignment done
|       - score if an alignment was done
|     o qryStartUL
|       - 0 if no alignment done
|       - first aligned query base if alignment done
|     o qryEndtUL
|       - 0 if no alignment done
|       - last aligned query base if alignment done
|     o refStartUL
|       - 0 if no alignment done
|       - first aligned reference base if alignment done
|     o refEndtUL
|       - 0 if no alignment done
|       - last aligned reference base if alignment done
|   - Returns:
|     o 1 if the reference sequence was found in the
|       kmer table (query) sequence
|     o 2 if the reverse alignment was best (may not have
|       been found)
|     o 0 if reference sequence not found
\-------------------------------------------------------*/
signed char
findRefInChunk_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   struct refST_kmerFind *refSTPtr,
   struct alnSet *alnSetPtr,
   float minPercScoreF,
   signed long *scoreSL,
   unsigned long *qryStartUL,
   unsigned long *qryEndUL,
   unsigned long *refStartUL,
   unsigned long *refEndUL
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun22 TOC:
   '   - finds an sequence pattern in an sam entry
   '   o fun22 sec01:
   '     - variable declerations
   '   o fun22 sec02:
   '     - initialize & see if enough kmers for alignment
   '   o fun22 sec03:
   '     - prepare for alignemnt (if passed kmer check)
   '   o fun22 sec04:
   '     - do alignment and check if passes min score
   '   o fun22 sec05:
   '     - return the answer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char matchBl = 0;
   float percScoreF = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec02:
   ^   - initialize and see if enough kmers for alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *qryStartUL = 0;
   *qryEndUL = 0;

   *refStartUL = 0;
   *refEndUL = 0;

   *scoreSL = 0;

   matchBl =
      matchCheck_kmerFind(
         tblSTPtr,
         refSTPtr
      );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec03:
   ^   - prepare for alignemnt (if passed kmer check)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(matchBl & 1)
   { /*If: I had enough kmers to do an alignment*/
      /*find the max score possible*/
      /*start of alignment region*/
      tblSTPtr->seqSTPtr->offsetUL = tblSTPtr->seqPosUL;
      tblSTPtr->seqSTPtr->offsetUL -= 
         (tblSTPtr->ntInWinUI -1);
         /*-1 to account for seqPosUL being index 1*/

      /*find end of alignment region*/
      tblSTPtr->seqSTPtr->endAlnUL = tblSTPtr->seqPosUL;

      /*convert index 1 to index 0*/
      --tblSTPtr->seqSTPtr->endAlnUL;

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun22 Sec04:
      ^   - do alignment and check if passes min score
      ^   o fun22 sec04 sub01:
      ^     - do the alignment
      ^   o fun22 sec04 sub02:
      ^     - check if it passes the alignment
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      /**************************************************\
      * Fun22 Sec04 Sub01:
      *   - do the alignment
      \**************************************************/

      if(matchBl & 2)
      { /*If: this was an reverse alignment*/
         *(scoreSL) =
            memwater(
              tblSTPtr->seqSTPtr,
              refSTPtr->revSeqST,
              refStartUL,
              refEndUL,
              qryStartUL,
              qryEndUL,
              alnSetPtr
            ); /*align primer to region*/

         percScoreF = (float) *scoreSL;
         percScoreF /= refSTPtr->maxRevScoreF;
      } /*If: this was an reverse alignment*/

      else
      { /*Else: this is an foward alignment*/
         *(scoreSL) =
            memwater(
              tblSTPtr->seqSTPtr,
              refSTPtr->forSeqST,
              refStartUL,
              refEndUL,
              qryStartUL,
              qryEndUL,
              alnSetPtr
            ); /*align primer to region*/

         percScoreF = (float) *scoreSL;
         percScoreF /= refSTPtr->maxForScoreF;
      } /*Else: this is an foward alignment*/

      /**************************************************\
      * Fun22 Sec04 Sub02:
      *   - check if it passes the alignment
      \**************************************************/

      matchBl &= ( -(percScoreF >= minPercScoreF) );
      ++tblSTPtr->seqSTPtr->endAlnUL;
   } /*If: I had enough kmers to do an alignment*/

   else
     matchBl = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec05:
   ^   - return the answer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return matchBl;
} /*findRefInChunk_kmerFind*/

/*-------------------------------------------------------\
| Fun24: waterFindPrims_kmerFind
|   - finds primers in an sequence (from fastx file) using
|     a slower, but more percise waterman
| Input:
|   - refAryST
|     o array of refST_kmerFind structures with reference
|       (primer) sequences to search for
|   - lenRefAryUI:
|     o number of refST_kmerFind structures in refSTAry
|   - seqSTPtr:
|     o pointer to an seqST structer with the query
|       sequence
|   - minPerScoreF:
|     o float with minimum percent score to keep an
|       alingment
|   - codeAryUI:
|     o unsigned integer array to hold the number of times
|       each spacer was detected
|   - dirArySC:
|     o pointer to a signed char array to hold mapped
|       primer directions
|   - scoreArySL:
|     o array of signed longs with the best score for each
|       matched primer
|   - seqStartAryUL:
|     o array of unsigned longs with the starting position
|       one the sequence for each score in scoreArySL
|   - endAyUL:
|     o array of unsigned longs with the ending position
|       on the sequence for each score in scoreArySL
|   - primStartAryUL:
|     o array of unsigned longs with the starting position
|       on the primer for each score in scoreArySL
|   - primEndAyUL:
|     o array of unsigned longs with the ending position
|       on the primer for score in scoreArySL
|   - alnSetPtr:
|     o pointer to an alnSet structure with the alignment
|       settings
| Output:
|   - Modifies:
|     o each position codeStr
|       - 1 if detected at one primer
|       - 2 if detected at both primers
|       - 0 if no primers were found
|     o dirAryStr to have direction of best alignment
|       - F for foward
|       - R for reverse
|     o scoreArySL score of best alignment for each primer
|     o seqStartAryUL starting sequence position of best
|       alignment for each primer
|     o seqEndAryUL ending sequence position of best
|       alignment for each primer
|     o primgStartAryUL starting primer position of best
|       alignment for each primer
|     o primgEndAryUL ending primer position of best
|       alignment for each primer
|   - Returns:
|     o 0 if the reference sequence was found in the kmer
|       table (query) sequence
|     o def_noPrim_kmerFind if no primers were found
\-------------------------------------------------------*/
signed char
waterFindPrims_kmerFind(
   struct refST_kmerFind *refSTAry,
   unsigned int lenRefAryUI,
   struct seqST *seqSTPtr,
   float minPercScoreF,
   unsigned int codeAryUI[],
   signed char dirArySC[],
   signed long scoreArySL[],
   unsigned long seqStartAryUL[],
   unsigned long seqEndAryUL[],
   unsigned long primStartAryUL[],
   unsigned long primEndAryUL[],
   struct alnSet *alnSetPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ^ Fun24 TOC:
   '   o fun24 sec01:
   '     - varaible declerations
   '   o fun24 sec02:
   '     - assign sequence to table
   '   o fun24 sec03:
   '     - check sequence for spacers
   '   o fun24 sec04:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun24 Sec01:
   ^   - varaible declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned char retUC = def_noPrim_kmerFind; /*return*/

   unsigned int uiPrim = 0;
   float percScoreF = 0;

   /*alignemnt variables; I keep*/
   signed long scoreSL = 0;
   unsigned long qryStartUL = 0;
   unsigned long qryEndUL = 0;

   /*alignemnt variables; I ingore*/
   unsigned long refStartUL = 0;
   unsigned long refEndUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun24 Sec02:
   ^   - convert to sequence to index
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   seqToIndex_alnSet(seqSTPtr->seqStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun24 Sec03:
   ^   - check sequence for primers
   ^   o fun24 sec03 sub01:
   ^     - start primer loop
   ^   o fun24 sec03 sub02:
   ^     - foward alignment of primer
   ^   o fun24 sec03 sub03:
   ^     - reverse complement alignment of primer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun24 Sec03 Sub01:
   *   - start primer loop
   \*****************************************************/

   for(
      uiPrim = 0;
      uiPrim < lenRefAryUI;
      ++uiPrim
   ){ /*Loop: detect primers in each chunk*/

      /**************************************************\
      * Fun24 Sec03 Sub02:
      *   - foward alignment of primer
      \**************************************************/

      scoreSL =
         memwater(
            seqSTPtr,
            refSTAry[uiPrim].forSeqST,
            &refStartUL,
            &refEndUL,
            &qryStartUL,
            &qryEndUL,
            alnSetPtr
         );

      percScoreF = (float) scoreSL;
      percScoreF /= refSTAry[uiPrim].maxForScoreF;

      if(percScoreF >= minPercScoreF)
      { /*If: foward primer found*/
         retUC = 0;
         codeAryUI[uiPrim] = 1;

         scoreArySL[uiPrim] = scoreSL;

         seqStartAryUL[uiPrim] = qryStartUL;
         seqEndAryUL[uiPrim] = qryEndUL;

         primStartAryUL[uiPrim] = refStartUL;
         primEndAryUL[uiPrim] = refEndUL;
         dirArySC[uiPrim] = 'F'; /*forward*/

         /*I need to make sure there is not a better
         `  score possible. Otherwise I will miss some
         `  primer pairs
         */
      } /*If: foward primer found*/

      else
      { /*Else: no foward match*/
         codeAryUI[uiPrim] = 0;
         scoreArySL[uiPrim] = 0;

         seqStartAryUL[uiPrim] = 0;
         seqEndAryUL[uiPrim] = 0;

         primStartAryUL[uiPrim] = 0;
         primEndAryUL[uiPrim] = 0;
         dirArySC[uiPrim] = 'N'; /*forward*/
      } /*Else: no foward match*/

      /**************************************************\
      * Fun24 Sec03 Sub03:
      *   - reverse complement alignment of primer
      \**************************************************/

      scoreSL =
         memwater(
            seqSTPtr,
            refSTAry[uiPrim].revSeqST,
            &refStartUL,
            &refEndUL,
            &qryStartUL,
            &qryEndUL,
            alnSetPtr
         );

      percScoreF = (float) scoreSL;
      percScoreF /= refSTAry[uiPrim].maxRevScoreF;

      if(percScoreF >= minPercScoreF)
      { /*If: i found this primer*/
         ++codeAryUI[uiPrim];

         if(scoreSL > scoreArySL[uiPrim])
         { /*If: I found a better score*/
             retUC = 0;
             scoreArySL[uiPrim] = scoreSL;

             seqStartAryUL[uiPrim] = qryStartUL;
             seqEndAryUL[uiPrim] = qryEndUL;

             primStartAryUL[uiPrim] = refStartUL;
             primEndAryUL[uiPrim] = refEndUL;

             dirArySC[uiPrim] = 'R'; /*reverse*/
         } /*If: I found a better score*/
      } /*If: i found this primer*/
   } /*Loop: detect primers in each chunk*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun24 Sec04:
   ^   - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   indexToSeq_alnSet(seqSTPtr->seqStr);
   return (signed char) retUC;
} /*waterFindPrims_kmerFind*/

/*-------------------------------------------------------\
| Fun24: fxFindPrims_kmerFind
|   - finds primers in an sequence (from fastx file) using
|     an faster kmer search followed by an slower waterman
|     to finalize alignments
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure with
|       settings and to hold the query sequence to check
|   - refAryST
|     o array of refST_kmerFind structures with reference
|       (primer) sequences to search for
|   - lenRefAryUI:
|     o number of refST_kmerFind structures in refSTAry
|   - seqSTPtr:
|     o pointer to an seqST structure with the
|       sequence to check for primers in
|   - minPerScoreF:
|     o float with minimum percent score to keep an
|       alingment
|   - codeAryUI:
|     o unsigned integer array to hold the number of times
|       each spacer was detected
|   - dirArySC:
|     o pointer to a signed char array to hold mapped
|       primer directions
|   - scoreArySL:
|     o array of signed longs with the best score for each
|       matched primer
|   - seqStartAryUL:
|     o array of unsigned longs with the starting position
|       one the sequence for each score in scoreArySL
|   - endAyUL:
|     o array of unsigned longs with the ending position
|       on the sequence for each score in scoreArySL
|   - primStartAryUL:
|     o array of unsigned longs with the starting position
|       on the primer for each score in scoreArySL
|   - primEndAyUL:
|     o array of unsigned longs with the ending position
|       on the primer for score in scoreArySL
|   - alnSetPtr:
|     o pointer to an alnSet structure with the alignment
|       settings
| Output:
|   - Modifies:
|     o each position codeStr to have number of times I
|       detected each primer
|     o dirAryStr to have direction of best alignment
|       - F for foward
|       - R for reverse
|     o scoreArySL score of best alignment for each primer
|     o seqStartAryUL starting sequence position of best
|       alignment for each primer
|     o seqEndAryUL ending sequence position of best
|       alignment for each primer
|     o primgStartAryUL starting primer position of best
|       alignment for each primer
|     o primgEndAryUL ending primer position of best
|       alignment for each primer
|   - Returns:
|     o 0 if the reference sequence was found in the kmer
|       table (query) sequence
|     o def_noPrim_kmerFind if no primers were found
\-------------------------------------------------------*/
signed char
fxFindPrims_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   struct refST_kmerFind *refSTAry,
   unsigned int lenRefAryUI,
   struct seqST *seqSTPtr,
   float minPercScoreF,
   unsigned int codeAryUI[],
   signed char dirArySC[],
   signed long scoreArySL[],
   unsigned long seqStartAryUL[],
   unsigned long seqEndAryUL[],
   unsigned long primStartAryUL[],
   unsigned long primEndAryUL[],
   struct alnSet *alnSetPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ^ Fun24 TOC:
   '   o fun24 sec01:
   '     - varaible declerations
   '   o fun24 sec02:
   '     - assign sequence to table
   '   o fun24 sec03:
   '     - check sequence for spacers
   '   o fun24 sec04:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun24 Sec01:
   ^   - varaible declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0; /*error messages*/
   signed char retSC = def_noPrim_kmerFind; /*return*/

   unsigned char matchBl = 0;
   signed char firstTimeBl = 1;

   unsigned int uiPrim = 0;

   /*alignemnt variables; I keep*/
   signed long scoreSL = 0;
   unsigned long qryStartUL = 0;
   unsigned long qryEndUL = 0;

   /*alignemnt variables; I ingore*/
   unsigned long refStartUL = 0;
   unsigned long refEndUL = 0;

   /*to keep the old assigned sequence*/
   struct seqST *oldSeqST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun24 Sec02:
   ^   - check positions and assign sequence to table
   ^   o fun24 sec02 sub01:
   ^     - see if i have an direct repeat region
   ^   o fun24 sec02 sub02:
   ^     - add the sequence to the kmer table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun24 Sec02 Sub01:
   *   - see if i have an direct repeat region
   \*****************************************************/

   for(
      uiPrim = 0;
      uiPrim < lenRefAryUI;
      ++uiPrim
   ){ /*Loop: blank my arrays*/
      codeAryUI[uiPrim] = 0; /*blank for new seq*/
      scoreArySL[uiPrim] = 0;
      dirArySC[uiPrim] = 0;
      seqStartAryUL[uiPrim] = 0;
      seqEndAryUL[uiPrim] = 0;
      primStartAryUL[uiPrim] = 0;
      primEndAryUL[uiPrim] = 0;
   } /*Loop: blank my arrays*/

   /*****************************************************\
   * Fun24 Sec02 Sub02:
   *   - add the sequence to the kmer table
   \*****************************************************/

   retSC = def_noPrim_kmerFind; /*allows detection*/

   oldSeqST = tblSTPtr->seqSTPtr;
   tblSTPtr->seqSTPtr = seqSTPtr;

   qckBlank_tblST_kmerFind(
      tblSTPtr,
      0          /*I do not want to blank the sequence*/
   ); /*using quick blank here, since all non-filled
		` kmers will already be blanked
		*/

   /*prepare for waterman alignment*/
   seqToIndex_alnSet(tblSTPtr->seqSTPtr->seqStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun24 Sec03:
   ^   - check sequence for spacers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   do{
      errSC =
         nextSeqChunk_tblST_kmerFind(
            tblSTPtr,
            &firstTimeBl
         );

      for(
         uiPrim = 0;
         uiPrim < lenRefAryUI;
         ++uiPrim
      ){ /*Loop: detect primers in each chunk*/
         matchBl =
            findRefInChunk_kmerFind(
               tblSTPtr,
               &(refSTAry[uiPrim]),
               alnSetPtr,
               minPercScoreF,
               &scoreSL,
               &qryStartUL,
               &qryEndUL,
               &refStartUL,
               &refEndUL
            );

         codeAryUI[uiPrim] += matchBl;

         if(scoreSL > scoreArySL[uiPrim])
         { /*If: I found a better score*/
             scoreArySL[uiPrim] = scoreSL;

             seqStartAryUL[uiPrim] = qryStartUL;
             seqEndAryUL[uiPrim] = qryEndUL;

             primStartAryUL[uiPrim] = refStartUL;
             primEndAryUL[uiPrim] = refEndUL;

             if(matchBl & 2)
                dirArySC[uiPrim] = 'R'; /*reverse*/
             else
                dirArySC[uiPrim] = 'F'; /*forward*/
         } /*If: I found a better score*/

         retSC &= (!matchBl);
      } /*Loop: detect primers in each chunk*/
   } while(! errSC);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun24 Sec04:
   ^   - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   qckBlank_tblST_kmerFind(
      tblSTPtr,
      0          /*I do not want to blank the sequence*/
   ); /*using quick blank here, since all non-filled
		` kmers will already be blanked
		*/

   /*reset table to original sequences*/
   indexToSeq_alnSet(tblSTPtr->seqSTPtr->seqStr);
   tblSTPtr->seqSTPtr = oldSeqST;

   return (signed char) retSC;
} /*fxFindPrims_kmerFind*/

/*-------------------------------------------------------\
| Fun25: phit_kmerFind
|   - prints out the primer hits for a sequence
| Input:
|   - refAryST:
|     o Array of refST_kmerFind structures with the primer
|       sequences
|   - numRefsSI:
|     o number of refST_kmerFind (primers) structures in
|       refAryST
|   - seqSTPtr:
|     o pointer to an seqST structer with the read id
|   - codeAryUI:
|     o pointer to a unsigned int array with the number of
|       times the sequence had each primer
|   - dirArySC:
|     o pointer to a signed char array with the mapped
|       direction for each primer
|   - scoreArySL:
|     o pointer to a signed long array with the score bes
|       score for each primer
|   - seqStartArySL:
|     o pointer to a unsigned long array with the primers
|       starting position on the sequence for the best
|       score
|   - seqEndArySL:
|     o pointer to a unsigned long array with the primers
|       ending position on the sequence for the best
|       score
|   - primStartArySL:
|     o pointer to a unsigned long array with the
|       sequences starting position on the primer for the
|       best score
|   - primEndArySL:
|     o pointer to a unsigned long array with the
|       sequences ending position on the primer for the
|       best score
|   - outFILE:
|     o file to print the stats to
\-------------------------------------------------------*/
void
phit_kmerFind(
   struct refST_kmerFind *refAryST,
   signed int numRefsSI,
   struct seqST *seqSTPtr,
   unsigned int *codeAryUI,
   signed char *dirArySC,
   signed long *scoreArySL,
   unsigned long *seqStartAryUL,
   unsigned long *seqEndAryUL,
   unsigned long *primStartAryUL,
   unsigned long *primEndAryUL,
   void *outFILE
){
   signed int siRef = 0;
   signed int siMate = 0;
   signed char *oldRefStr = 0;
   signed char *oldSeqStr = 0;
   signed char oldRefBreakSC = 0;
   signed char oldSeqBreakSC = 0;

   unsigned char skip1stCharUC = 0;

   oldSeqStr = seqSTPtr->idStr;
   while(*oldSeqStr++ > 32);
   --oldSeqStr;
   oldSeqBreakSC = *oldSeqStr;
   *oldSeqStr = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun25 Sec02:
   ^   - print out mapped primers
   ^   o fun25 sec02 sub01:
   ^     - start loop and check primers that mapped
   ^   o fun25 sec02 sub02:
   ^     - primer mate (pairing) stats printout
   ^   o fun25 sec02 sub03:
   ^     - non-primer mate (no pairing) stats printout
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun25 Sec02 Sub01:
   *   - start loop and check primers that mapped
   \*****************************************************/

   siRef = 0;
   skip1stCharUC = *refAryST[siRef].forSeqST->idStr =='>';

   while(siRef < numRefsSI)
   { /*Loop: print out hits*/
      if(codeAryUI[siRef] == 0)
         goto nextRef_fun25_sec02_sub04;

      if(dirArySC[siRef] == 'N')
         goto nextRef_fun25_sec02_sub04;

      /**************************************************\
      * Fun25 Sec02 Sub02:
      *   - primer mate (pairing) stats printout
      \**************************************************/

      if(refAryST[siRef].mateSI >= 0)
      { /*If: I have a mate primer*/
         siMate = refAryST[siRef].mateSI;

         if(codeAryUI[siMate] == 0)
         { /*If: the mate did not map*/
            ++siRef; /*move past mate*/
            goto nextRef_fun25_sec02_sub04;
         } /*If: the mate did not map*/

         if(dirArySC[siMate] == 'N')
         { /*If: the mate did not map*/
            ++siRef; /*move past mate*/
            goto nextRef_fun25_sec02_sub04;
         } /*If: the mate did not map*/

         if(dirArySC[siRef] == dirArySC[siMate])
         { /*If: the mate was in the same direction*/
            ++siRef; /*move past mate*/
            goto nextRef_fun25_sec02_sub04;
         } /*If: the mate was in the same direction*/

         fprintf(
             (FILE *) outFILE,
             "%s\t%s",
             seqSTPtr->idStr,
             refAryST[siRef].forSeqST->idStr
                + skip1stCharUC
         ); /*ids*/

         /*print out aligned length*/
         if(seqStartAryUL[siRef] < seqEndAryUL[siMate])
         { /*If: the mate comes last (forward read)*/
            fprintf(
               (FILE *) outFILE,
               "\t%lu",
               seqEndAryUL[siMate] - seqStartAryUL[siRef]
            ); /*sequence length including primers*/
         } /*If: the mate comes last (foward read)*/

         else
         { /*Else: the mate comes first (reverse read)*/
            fprintf(
               (FILE *) outFILE,
               "\t%lu",
               seqEndAryUL[siRef] - seqStartAryUL[siMate]
            ); /*sequence length including primers*/
         } /*Else: the mate comes first (reverse read)*/

         fprintf(
             (FILE *) outFILE,
             "\t%c\t%lu\t%lu",
             dirArySC[siRef],
             seqStartAryUL[siRef] + 1,
             seqEndAryUL[siRef] + 1
         ); /*forward primer mapping coordinates*/

         fprintf(
             (FILE *) outFILE,
             "\t%c\t%lu\t%lu",
             dirArySC[siMate],
             seqStartAryUL[siMate] + 1,
             seqEndAryUL[siMate] + 1
         ); /*reverse primer mapping coordinates*/

         oldRefStr = refAryST[siRef].forSeqST->idStr;
         oldRefStr += (*oldRefStr == '>');

         while(*oldRefStr > 32)
            ++oldRefStr;

         oldRefBreakSC = *oldRefStr;
         *oldRefStr = '\0';

         /*print primer ids*/
         fprintf(
             (FILE *) outFILE,
             "\t%0.2f\t%0.2f\t%i\t%lu\t%lu",
             (float)
                scoreArySL[siRef] / def_scoreAdj_alnDefs,
             (float)
                  refAryST[siRef].maxForScoreF
                / def_scoreAdj_alnDefs,
             codeAryUI[siRef],
             primStartAryUL[siRef] + 1,
             primEndAryUL[siRef] + 1
         ); /*forward primer mapping stats*/

         oldRefStr = refAryST[siMate].forSeqST->idStr;
         *oldRefStr = oldRefBreakSC;
         oldRefStr += (*oldRefStr == '>');

         while(*oldRefStr > 32)
            ++oldRefStr;

         oldRefBreakSC = *oldRefStr;
         *oldRefStr = '\0';

         fprintf(
             (FILE *) outFILE,
             "\t%0.2f\t%0.2f\t%i\t%lu\t%lu\n",
             (float)
                scoreArySL[siMate] / def_scoreAdj_alnDefs,
             (float)
                 refAryST[siRef].maxRevScoreF
               / def_scoreAdj_alnDefs,
             codeAryUI[siMate],
             primStartAryUL[siMate] + 1,
             primEndAryUL[siMate] + 1
         ); /*reverse primer mapping stats*/

         *oldRefStr = oldRefBreakSC;
         ++siRef; /*move past mate*/
      } /*If: I have a mate primer*/

      /**************************************************\
      * Fun25 Sec02 Sub03:
      *   - non-primer mate (no pairing) stats printout
      \**************************************************/

      else
      { /*Else: I have no mate primers*/
         fprintf(
             (FILE *) outFILE,
             "%s\t%s\tNA\t%c\t%lu\t%lu",
             seqSTPtr->idStr,
             refAryST[siRef].forSeqST->idStr
                + skip1stCharUC,
             dirArySC[siRef],
             seqStartAryUL[siRef] + 1,
             seqEndAryUL[siRef] + 1
         ); /*sequence alignment stats*/

         fprintf(
             (FILE *) outFILE,
             "\tNA\tNA\tNA"
         ); /*mate coordinates and mapped length*/

         oldRefStr = refAryST[siRef].forSeqST->idStr;
         oldRefStr += (*oldRefStr == '>');

         while(*oldRefStr > 32)
            ++oldRefStr;

         oldRefBreakSC = *oldRefStr;
         *oldRefStr = '\0';

         fprintf(
             (FILE *) outFILE,
             "\t%0.2f\t%0.2f\t%i\t%lu\t%lu",
             (float)
                scoreArySL[siRef] / def_scoreAdj_alnDefs,
             (float)
                  refAryST[siRef].maxForScoreF
                / def_scoreAdj_alnDefs,
             codeAryUI[siRef],        /*total mapps*/
             primStartAryUL[siRef] + 1,/*start of map*/
             primEndAryUL[siRef] + 1   /*end of map*/
         ); /*foward primer mapping stats*/

         *oldRefStr = oldRefBreakSC;

         fprintf(
             (FILE *) outFILE,
             "\tNA\tNA\tNA\tNA\tNA\n"
         ); /*reverse primer mapping stats*/
      } /*Else: I have no mate primers*/

      /**************************************************\
      * Fun25 Sec02 Sub04:
      *   - move to the next reference
      \**************************************************/

      nextRef_fun25_sec02_sub04:;

      ++siRef;
   } /*Loop: print out hits*/

   *oldSeqStr = oldSeqBreakSC;
} /*phit_kmerFInd*/

/*-------------------------------------------------------\
| Fun26: pHeaderHit_kmerFind
|    - prints header for phit_kmerFind (fun25)
| Input:
|   - outFILE:
|     o file to print header to
| Output:
|   - Prints:
|     o header to outFILE
\-------------------------------------------------------*/
void
pHeaderHit_kmerFind(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "read_id\tprim_id\taln_len\tfor_dir\tfor_seq_start"
   );

   fprintf(
      (FILE *) outFILE,
      "\tfor_seq_end\trev_dir\trev_start\trev_end"
   );

   fprintf(
      (FILE *) outFILE,
      "\tfor_score\tfor_max_score"
   );

   fprintf(
      (FILE *) outFILE,
      "\tfor_num_maps\tfor_prim_start\tfor_prim_end"
   );

   fprintf(
      (FILE *) outFILE,
      "\trev_score\trev_max_score\trev_num_maps"
   );

   fprintf(
      (FILE *) outFILE,
      "\trev_prim_start\trev_prim_end\n"
   );
} /*pHeader_kmerFind*/

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
