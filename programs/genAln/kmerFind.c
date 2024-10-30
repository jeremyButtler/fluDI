/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' kmerFind SOF: Start Of File
'   - holds functions for kmer detection of primers
'   o header:
'     - included libraries
'   o tbl01: alnNtToBit_kmerFind
'     - converts an nucleotide alignment code form
'       alnSetStruct.h to an two bit value, with an 3rd
'       bit being used for anonymous bases and the 4th bit
'       for errors
'   o .h st01: tblST_kmerFind
'     - holds the kmer tables for detecting spoligytpes
'   o .h st02: refST_kmerFind
'     - holds the kmer pattern for the reference
'   o fun01: blank_tblST_kmerFind
'     - blanks all stored values in an tblST_kmerFind
'   o fun02: init_tblST_kmerFind
'     - initializes an tblST_kmerFind structure
'   o fun03: setup_tblST_kmerFind
'     - allocates memory for variables in a tblST_kmerFind
'       structure
'   o fun04: freeStack_tblST_kmerFind
'     - frees all variables in an tblST_kmerFind structure
'   o fun05: freeHeap_tblST_kmerFind
'     - frees and tblST_kmerFind structure
'   o fun06: blank_refST_kmerFind
'     - sets the counting and kmer values in an
'       refST_kmerFind to 0 or def_noKmer_kmerFind
'   o fun07: init_refST_kmerFind
'     - initializes an refST_kmerFind structure
'   o fun08: setup_refST_kmerFind
'     - allocates memory for a refST_kmerFind structure
'   o fun09: freeStack_refST_kmerFind
'     - frees the variables in an refST_kmerFind structure
'   o fun10: freeHeap_refST_kmerFind
'     - frees an refST_kmerFind structure
'   o fun11: freeHeapAry_refST_kmerFind
'     - frees an array of refST_kmerFind structure
'   o fun12: uiTwinShort_kmerFind
'     - sorts a unsigned int array from least to greatest
'       and keeps it linked to a second unsigned int array
'   o fun13: addSeqToRefST_kmerFInd
'     - adds a sequence to a refST_kmerFind structure
'   o fun14: prep_tblST_kmerFind
'     - sets up an tblST_kmerFind structure for primer
'       searching
'   o fun15: tsvToAry_refST_kmerFind
'     - makes an array of refST_kmerFind structures from a
'       tsv file
'   o fun16: faToAry_refST_kmerFind
'     - makes an array of refST_kmerFind structures
'   o fun17: nextSeqChunk_tblST_kmerFind
'     - adds a new set of kmers from an sequence to an
'       tblST_kmerFind structure
'   o fun18: forCntMatchs_kmerFind
'     - finds the number of kmers that are in both the
'       kmer table (query) and the pattern (reference)
'   o fun19: revCntMatchs_kmerFind
'     - finds the number of kmers that are shared in the
'       kmer table (query) and the reverse pattern
'       (reference)
'   o fun20: matchCheck_kmerFind
'     - tells if the  match meets the min requirements to
'       do an alignment or not
'   o fun21: findRefInChunk_kmerFind
'     - does an kmer check and alings an single sequence
'       in an refST_kmerFind structure to see if there is
'       an match
'   o fun22: waterFindPrims_kmerFind
'     - finds primers in an sequence (from fastx file)
'       using a slower, but more percise waterman
'   o fun23: fxFindPrims_kmerFind
'     - finds spoligotype spacers in an sequence (from
'       fastx file) using an faster kmer search followed
'       by an slower waterman to finalize alignments
'   o fun24: phit_kmerFind
'     - prints out the primer hits for a sequence
'   o fun25: pHeaderHit_kmerFind
'      - prints header for phit_kmerFind (fun24)
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
#include "../genLib/dataTypeShortHand.h"
#include "../genLib/genMath.h" /*only using .h macros*/

#include "alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../genLib/base10Str.h"
!   - .c  #include "../genAln/indexToCoord.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_a_kmerFind 0
#define def_c_kmerFind 1
#define def_g_kmerFind 2
#define def_t_kmerFind 3
#define def_anonNt_kmerFind 4 /*anonmyous nucleotide*/
#define def_noNt_kmerFind 8

/*-------------------------------------------------------\
| Tbl01: alnNtToBit_kmerFind
|   - converts an nucleotide alignment code form
|     alnSetStruct.h to an two bit value, with an 3rd
|     bit being used for anonymous bases and the 4th bit
|     for errors
\-------------------------------------------------------*/
unsigned char alnNtToBit_kmerFind[] =
{
   def_noNt_kmerFind,    /*00 = 64 = ,*/  
   def_a_kmerFind,       /*01 = 65 = A/a*/  
   def_anonNt_kmerFind,  /*02 = 66 = B/b*/  
   def_c_kmerFind,       /*03 = 67 = C/c*/  
   def_anonNt_kmerFind,  /*04 = 68 = D/d*/  
   def_noNt_kmerFind,    /*05 = 69 = E/e*/  
   def_noNt_kmerFind,    /*06 = 70 = F/f*/  
   def_g_kmerFind,       /*07 = 71 = G/g*/  
   def_anonNt_kmerFind,  /*08 = 72 = H/h*/  
   def_noNt_kmerFind,    /*09 = 73 = I/i*/  
   def_noNt_kmerFind,    /*10 = 74 = J/j*/  
   def_anonNt_kmerFind,  /*11 = 75 = K/k*/  
   def_noNt_kmerFind,    /*12 = 76 = L/l*/  
   def_anonNt_kmerFind,  /*13 = 77 = M/m*/  
   def_anonNt_kmerFind,  /*14 = 78 = N/n*/  
   def_noNt_kmerFind,    /*15 = 79 = O/o*/  
   def_noNt_kmerFind,    /*16 = 80 = P/p*/  
   def_noNt_kmerFind,    /*17 = 81 = Q/q*/  
   def_anonNt_kmerFind,  /*18 = 82 = R/r*/  
   def_anonNt_kmerFind,  /*19 = 83 = S/s*/  
   def_t_kmerFind,       /*20 = 84 = T/t*/  
   def_t_kmerFind,       /*21 = 85 = U/u*/  
   def_anonNt_kmerFind,  /*22 = 86 = V/v*/  
   def_anonNt_kmerFind,  /*23 = 87 = W/w*/  
   def_anonNt_kmerFind,  /*24 = 88 = X/x (amino acids)*/  
   def_anonNt_kmerFind,  /*25 = 89 = Y/y*/  
   def_noNt_kmerFind     /*26 = 90 = Z/z*/  
}; /*ntAlnCodeToBitTbl*/

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

   tblSTPtr->lastKmerUL = 0;
   tblSTPtr->lenLastKmerUL = 0;
   tblSTPtr->seqPosUL = 0;

   if(tblSTPtr->tblSI)
   { /*If: I have an kmer table to blank*/
      for(
         uiKmer = 0;
         uiKmer < tblSTPtr->lenTblUI;
         ++uiKmer
      ) tblSTPtr->tblSI[uiKmer] = 0;
   } /*If: I have an kmer table to blank*/

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
| Fun02: init_tblST_kmerFind
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
| Fun03: setup_tblST_kmerFind
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
   ' Fun03 TOC:
   '   - allocates memory for variables in a
   '     tblST_kmerFind
   '   o fun03 sec01:
   '     - allocate memory for the table
   '   o fun03 sec02:
   '     - build the kmer mask
   '   o fun03 sec03:
   '     - initialize seqST and blank table
   '   o fun03 sec04:
   '     - return any errors
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - allocate memory for the table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint uiKmer = 0;

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
      malloc((tblSTPtr->lenTblUI + 1) * sizeof(sint));

   if(! tblSTPtr->tblSI)
      goto memErr_fun03_sec03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - build the kmer mask
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   tblSTPtr->kmerMaskUL =
      (lenKmerUC * def_bitsPerKmer_kmerFind);

   tblSTPtr->kmerMaskUL =
      (sizeof(ulong) << 3) - tblSTPtr->kmerMaskUL;

   tblSTPtr->kmerMaskUL =
      ((ulong) -1) >> tblSTPtr->kmerMaskUL;

   /*logic
   `   - mask = lenKmerUC << bitsPerKmer:
   `     o gives me the number of bits an single kmer
   `       takes up
   `   - mask = (sizeof(ulong) << 3) - mask
   `     o sizeof(ulong) << 3 gives the number of bits
   `       in an unsigned long
   `     o bits in long - mask gives the number of bits
   `       that need to be masked (not in kmer)
   `   -  mask = -1 >> mask
   `      o removes all the leading ones that should be
   `        masked out (111...11111 to 000...011111)
   */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - initialize seqST and blank table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(tblSTPtr->seqSTPtr)
     freeHeap_seqST(tblSTPtr->seqSTPtr);

   tblSTPtr->seqSTPtr = 0;

   tblSTPtr->seqSTPtr = malloc(sizeof(struct seqST));

   if(! tblSTPtr->seqSTPtr)
      goto memErr_fun03_sec03;

   init_seqST(tblSTPtr->seqSTPtr);

   blank_tblST_kmerFind(
      tblSTPtr,
      0 /*sequence structure is already blanked*/
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - return any errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun03_sec03:;

   freeStack_tblST_kmerFind(tblSTPtr);

   return def_memErr_kmerFind;
} /*setup_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun04: freeStack_tblST_kmerFind
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
   if(tblSTPtr)
   { /*If: I have an struct to free*/
      if(tblSTPtr->kmerArySI)
         free(tblSTPtr->kmerArySI);

      tblSTPtr->kmerArySI = 0;

      if(tblSTPtr->seqSTPtr)
         freeHeap_seqST(tblSTPtr->seqSTPtr);

      tblSTPtr->seqSTPtr = 0;

      if(tblSTPtr->tblSI)
         free(tblSTPtr->tblSI);

      tblSTPtr->tblSI = 0;

      blank_tblST_kmerFind(
         tblSTPtr,
         0 /*no seqST to blank*/
      );
   } /*If: I have an struct to free*/
} /*freeStack_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun05: freeHeap_tblST_kmerFind
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
   if(tblSTPtr)
   { /*If: have table to free*/
      freeStack_tblST_kmerFind(tblSTPtr);
      free(tblSTPtr);
   } /*If: have table to free*/
} /*freeHeap_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun06: blank_refST_kmerFind
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
   ulong ulKmer = 0;
   refSTPtr->mateSI = -1;

   if(refSTPtr->forKmerArySI)
   { /*If: I have kmer arrays allocated*/
      for(
         ulKmer = 0;
         ulKmer < refSTPtr->lenAryUI;
         ++ulKmer
      ){ /*Loop: blank kmer arrays*/
         refSTPtr->forKmerArySI[ulKmer] =
            def_noKmer_kmerFind;

         refSTPtr->forRepAryUI[ulKmer] = 0;

         refSTPtr->revKmerArySI[ulKmer] =
            def_noKmer_kmerFind;

         refSTPtr->revRepAryUI[ulKmer] = 0;
      } /*Loop: blank kmer arrays*/

      refSTPtr->forKmerArySI[ulKmer] =
         def_endKmers_kmerFind;

      refSTPtr->revKmerArySI[ulKmer] =
         def_endKmers_kmerFind;
   } /*If: I have kmer arrays allocated*/

   if(refSTPtr->forSeqST)
      blank_seqST(refSTPtr->forSeqST);

   if(refSTPtr->revSeqST)
      blank_seqST(refSTPtr->revSeqST);
} /*blank_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun07: init_refST_kmerFind
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
   refSTPtr->lenKmerUC = 0;
   refSTPtr->minKmersUI = 0;
   refSTPtr->maxForScoreF = 0;
   refSTPtr->maxRevScoreF = 0;
   refSTPtr->lenAryUI = 0;
   refSTPtr->lenAryUI = 0;

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
| Fun08: setup_refST_kmerFind
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
| Fun09: freeStack_refST_kmerFind
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
   if(refSTPtr)
   { /*If: I have an structure to free*/
      if(refSTPtr->forKmerArySI)
         free(refSTPtr->forKmerArySI);
      refSTPtr->forKmerArySI = 0;

      if(refSTPtr->forRepAryUI)
         free(refSTPtr->forRepAryUI);
      refSTPtr->forRepAryUI = 0;

      if(refSTPtr->revKmerArySI)
         free(refSTPtr->revKmerArySI);
      refSTPtr->revKmerArySI = 0;

      if(refSTPtr->revRepAryUI)
         free(refSTPtr->revRepAryUI);
      refSTPtr->revRepAryUI = 0;

      if(refSTPtr->forSeqST)
         freeHeap_seqST(refSTPtr->forSeqST);
      refSTPtr->forSeqST = 0;

      if(refSTPtr->revSeqST)
         freeHeap_seqST(refSTPtr->revSeqST);
      refSTPtr->revSeqST = 0;
   } /*If: I have an structure to free*/

   /*not needed, but for future additions*/
   blank_refST_kmerFind(refSTPtr);
} /*freeStack_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun10: freeHeap_refST_kmerFind
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
   if(refSTPtr)
   { /*If: have structure to free*/
      freeStack_refST_kmerFind(refSTPtr);
      free(refSTPtr);
   } /*If: have structure to free*/
} /*freeStack_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun11: freeHeapAry_refST_kmerFind
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
   --lenArySI; /*convert index 1 to index 0*/

   while(lenArySI > -1)
   { /*Loop: free structres in array*/
      freeStack_refST_kmerFind(&refSTAry[lenArySI]);
      --lenArySI;
   } /*Loop: free structres in array*/

   free(refSTAry);
} /*freeHeapAry_refST_kmerFind*/

/*-------------------------------------------------------\
| Fun12: uiTwinShort_kmerFind
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
   ' Fun12 TOC:
   '   - sorts signed int array from least to greatest
   '   o fun12 sec01:
   '     - variable declerations
   '   o fun12 sec02:
   '     - find the number of rounds to sort for
   '   o fun12 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec01:
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
   ^ Fun12 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Recursion formuia: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subUL = 1; /*Initialzie first array*/

   while(subUL < numElmUL - 1)
      subUL = (3 * subUL) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec03:
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

                  uiSecAry[uiElm] ^= uiSecAry[nextUL];
                  uiSecAry[nextUL] ^= uiSecAry[uiElm];
                  uiSecAry[uiElm] ^= uiSecAry[nextUL];

                  onUL = lastUL;
               } /*Loop: move swapped element back*/
            } /*If I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subUL = (subUL - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*uiTwinSort_kmerFind*/

/*-------------------------------------------------------\
| Fun13: addSeqToRefST_kmerFInd
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
   ' Fun13 TOC:
   '   - add a sequence to a refST_kmerFind structure
   '   o fun13 sec01:
   '     - variable declerations
   '   o fun13 sec02:
   '     - copy sequence and see if is longest sequence
   '   o fun13 sec03:
   '     - find minimum number of kmers (for kmer search)
   '   o fun13 sec04:
   '     - add kmers to the kmer arrays
   '   o fun13 sec05:
   '     - move empty (not in sequence) kmers to the end
   '   o fun13 sec06:
   '     - return the longest sequence or 0 for errors
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uchar errUC = 0;
   uchar ntUC = 0;

   sint siKmer = 0;
   sint siSeq = 0;

   /*for building kmers*/
   ulong forKmerUL = 0;
   ulong lenForKmerUL = 0;

   ulong revKmerUL = 0;
   ulong lenRevKmerUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec02:
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
         goto memErr_fun13_sec06;
   } /*If: given a sequence structure*/

   refSTPtr->forSeqST->endAlnUL =
      refSTPtr->forSeqST->lenSeqUL - 1;

   refSTPtr->forSeqST->offsetUL = 0;

   longestSeqUI =
      max_genMath(
         longestSeqUI,
         (uint) refSTPtr->forSeqST->lenSeqUL
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
      goto memErr_fun13_sec06;

   revComp_seqST(refSTPtr->revSeqST);

   /*convert sequences to the correct format*/
   seqToIndex_alnSet(refSTPtr->forSeqST->seqStr);
   seqToIndex_alnSet(refSTPtr->revSeqST->seqStr);
   
   /*I am merging duplicates, so I never expect
   `   more then the maxiumum number of possible kmers
   */
   refSTPtr->lenAryUI = tblSTPtr->lenTblUI;
   refSTPtr->lenKmerUC = tblSTPtr->lenKmerUC;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec03:
   ^   - set up the kmer arrays
   ^   o fun13 sec03 sub01:
   ^     -  find minimum number of kmers (for kmer search)
   ^   o fun13 sec03 sub02:
   ^     - allocate memory for forward kmer arrays
   ^   o fun13 sec03 sub03:
   ^     - allocate reverse kmer array memory
   ^   o fun13 sec03 sub04:
   ^     - initialize kmer counts
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun13 Sec03 Sub01:
   *   -  find minimum number of kmers (for kmer search)
   \*****************************************************/

   refSTPtr->minKmersUI = 
      (uint)
      refSTPtr->forSeqST->lenSeqUL;

   refSTPtr->minKmersUI -= refSTPtr->lenKmerUC;
   ++refSTPtr->minKmersUI; /*total kmers*/
   refSTPtr->minKmersUI *= minPercKmersF;

   /*************************************************\
   * Fun13 Sec03 Sub02:
   *   - allocate memory for forward kmer arrays
   \*************************************************/
   
   refSTPtr->forKmerArySI =
      malloc((refSTPtr->lenAryUI + 1) * sizeof(uint));

   if(! refSTPtr->forKmerArySI)
      goto memErr_fun13_sec06;

   refSTPtr->forRepAryUI =
      malloc((refSTPtr->lenAryUI + 1) * sizeof(uint));

   if(! refSTPtr->forRepAryUI)
      goto memErr_fun13_sec06;

   /*************************************************\
   * Fun13 Sec03 Sub03:
   *   - allocate reverse kmer array memory
   \*************************************************/
   
   refSTPtr->revKmerArySI =
      malloc((refSTPtr->lenAryUI + 1) * sizeof(uint));

   if(! refSTPtr->revKmerArySI)
      goto memErr_fun13_sec06;

   refSTPtr->revRepAryUI =
      malloc((refSTPtr->lenAryUI + 1) * sizeof(uint));

   if(! refSTPtr->revRepAryUI)
      goto memErr_fun13_sec06;

   /*************************************************\
   * Fun13 Sec03 Sub04:
   *   - initialize kmer counts
   \*************************************************/

   for(
      siKmer = 0;
      siKmer < (sint) refSTPtr->lenAryUI;
      ++siKmer
   ){ /*Loop: inititalize kmer counts*/
      refSTPtr->forRepAryUI[siKmer] = 0;
      refSTPtr->revRepAryUI[siKmer] = 0;

      refSTPtr->forKmerArySI[siKmer]= def_noKmer_kmerFind;
      refSTPtr->revKmerArySI[siKmer]= def_noKmer_kmerFind;
   } /*Loop: inititalize kmer counts*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec04:
   ^   - add kmers to the kmer arrays
   ^   o fun13 sec04 sub01:
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*************************************************\
   * Fun13 Sec05 Sub06:
   *   - get kmer counts
   \*************************************************/

   for(
      siSeq = 0;
       siSeq < (sint) refSTPtr->forSeqST->lenSeqUL;
      ++siSeq
   ){ /*Loop: copy the kmers*/

      /**********************************************\
      * Fun13 Sec05 Sub07:
      *   - get foward counts and max score
      \**********************************************/

      ntUC =
         (uchar)
         refSTPtr->forSeqST->seqStr[siSeq];

      refSTPtr->maxForScoreF +=
         getScore_alnSet(
            ntUC & def_ntToCode_alnSet,
            ntUC & def_ntToCode_alnSet,
            alnSetPtr
         );

      ntUC = alnNtToBit_kmerFind[ ntUC ];
      forKmerUL <<= def_bitsPerKmer_kmerFind;
      forKmerUL |= ntUC;
      forKmerUL &= tblSTPtr->kmerMaskUL;

      if(ntUC < def_anonNt_kmerFind)
      { /*If: no anymous bases or errors*/
         ++lenForKmerUL;

         if(lenForKmerUL >= refSTPtr->lenKmerUC)
         { /*If: had a forward kmer*/
            ++refSTPtr->forRepAryUI[forKmerUL];
            
            refSTPtr->forKmerArySI[forKmerUL] =
               (sint)
               forKmerUL;
         } /*If: had a forward kmer*/
      } /*If: no anymous bases or errors*/

      else
         lenForKmerUL = 0;

      /**********************************************\
      * Fun13 Sec05 Sub08:
      *   - get reverse counts and max score
      \**********************************************/

      ntUC =
         (uchar)
         refSTPtr->revSeqST->seqStr[siSeq];

      refSTPtr->maxRevScoreF +=
         getScore_alnSet(
            ntUC & def_ntToCode_alnSet,
            ntUC & def_ntToCode_alnSet,
            alnSetPtr
         );

      ntUC = alnNtToBit_kmerFind[ ntUC ];
      revKmerUL <<= def_bitsPerKmer_kmerFind;
      revKmerUL |= ntUC;
      revKmerUL &= tblSTPtr->kmerMaskUL;

      if(ntUC < def_anonNt_kmerFind)
      { /*If: no anymous bases or errors*/
         ++lenRevKmerUL;

         if(lenRevKmerUL >= refSTPtr->lenKmerUC)
         { /*If: had a revward kmer*/
            ++refSTPtr->revRepAryUI[revKmerUL];
            
            refSTPtr->revKmerArySI[revKmerUL] =
               (sint)
               revKmerUL;
         } /*If: had a revward kmer*/
      } /*If: no anymous bases or errors*/

      else
         lenRevKmerUL = 0;
   } /*Loop: copy the kmers*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec05:
   ^   - move empty (not in sequence) kmers to the end
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*I am converting the kmer numbers to unsigned ints
   `  so that -1's and -2's will be at the ends
   */

   ui_shellSort(
      (uint *) refSTPtr->forKmerArySI,
      0,
      refSTPtr->lenAryUI - 1
   ); /*sort the forward kmers by kmer id*/

   ui_shellSort(
      (uint *) refSTPtr->revKmerArySI,
      0,
      refSTPtr->lenAryUI - 1
   ); /*sort the revers kmers by kmer id*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec06:
   ^   - return the longest sequence or 0 for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return longestSeqUI;

   memErr_fun13_sec06:;
   return 0;
} /*addSeqToRefST_kmerFind*/

/*-------------------------------------------------------\
| Fun14: prep_tblST_kmerFind
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
   uint lenKmerAryUI = 0;
   sint siSeq = 0;

   /*find window size*/
   tblSTPtr->ntInWinUI = longestSeqUI;
   tblSTPtr->ntInWinUI *= percExtraNtInWinF;
   tblSTPtr->ntInWinUI += longestSeqUI;

   /*find the number of kmers to remove per window shift*/
   tblSTPtr->rmNtUI = tblSTPtr->ntInWinUI;
   tblSTPtr->rmNtUI *= percWinShiftF;

   /*find number of kmers in one window*/
   lenKmerAryUI = tblSTPtr->ntInWinUI;
   lenKmerAryUI -= tblSTPtr->lenKmerUC;
   ++lenKmerAryUI;

   /*set up the kmer window array*/
   if(tblSTPtr->numKmerUI < lenKmerAryUI)
   { /*If: the array is to small*/
      if(tblSTPtr->tblSI)
         free(tblSTPtr->kmerArySI); /*have an old table*/

      tblSTPtr->kmerArySI = 0;

      tblSTPtr->kmerArySI =
         malloc((lenKmerAryUI + 1) * sizeof(sint));

      if(! tblSTPtr->kmerArySI)
         goto memErr_fun16; /*memory error*/

      tblSTPtr->numKmerUI = lenKmerAryUI;
   } /*If: the array is to small*/

   /*initialize/blank the array*/
   for(
      siSeq = 0;
      siSeq < (sint) tblSTPtr->numKmerUI;
      ++siSeq
   ) tblSTPtr->kmerArySI[siSeq] = def_noKmer_kmerFind;
   
   tblSTPtr->kmerArySI[siSeq] = def_endKmers_kmerFind;

   return 0;

   memErr_fun16:;
   return def_memErr_kmerFind;
} /*prep_tblST_kmerFind*/

/*-------------------------------------------------------\
| Fun15: tsvToAry_refST_kmerFind
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
   ' Fun15 TOC:
   '   - makes an array of refST_kmerFind structures
   '   o fun15 sec01:
   '     - variable declerations 
   '   o fun15 sec02:
   '     - check if can open input file
   '   o fun15 sec03:
   '     - get number of sequences in file
   '   o fun15 sec04:
   '     - allocate memory
   '   o fun15 sec05:
   '     - read in sequences
   '   o fun15 sec06:
   '     - update the table window size values
   '   o fun15 sec07:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec01:
   ^   - variable declerations 
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define def_lenBuff_fun15 4096
   signed char buffStr[def_lenBuff_fun15];
   signed char *tmpStr = 0;

   uchar pairBl = 0;
   sint  numSeqSI = 0;
   uint longestPrimUI = 0; /*length of longest primer*/

   /*structures to work with*/
   struct refST_kmerFind *retRefHeapAryST = 0;
   struct seqST seqStackST;

   FILE *tsvFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec02:
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
      goto fileErr_fun15_sec07_sub03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec03:
   ^   - find number of lines in the file
   ^   o fun15 sec03 sub01:
   ^     - find number of lines in the file
   ^   o fun15 sec03 sub02:
   ^     - check for errors and move to start of file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun15 Sec03 Sub01:
   *   - find number of lines in the file
   \*****************************************************/

   while(fgets((char *) buffStr, def_lenBuff_fun15, tsvFILE))
      ++numSeqSI;

   --numSeqSI;     /*account for the header*/
   numSeqSI <<= 1; /*acount for two sequences per line*/

   /*****************************************************\
   * Fun15 Sec03 Sub04:
   *   - check for errors and move to start of file
   \*****************************************************/

   if(numSeqSI < 1)
      goto fileErr_fun15_sec07_sub03; /*no sequences*/

   fseek(
      tsvFILE,
      0,
      SEEK_SET
   ); /*find start of file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec04:
   ^   - allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retRefHeapAryST =
      malloc(numSeqSI * sizeof(struct refST_kmerFind));

   if(! retRefHeapAryST)
      goto memErr_fun15_sec07_sub02;

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
      goto memErr_fun15_sec07_sub02;
   } /*If: I had an memroy error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec05:
   ^   - read in sequences from tsv file
   ^   o fun15 sec05 sub01:
   ^     - set sequence counter to zero and start loop
   ^   o fun15 sec05 sub02:
   ^     - get primer id
   ^   o fun15 sec05 sub11:
   ^     - update the table window size values
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun15 Sec05 Sub01:
   *   - set sequence counter to zero and start loop
   \*****************************************************/

   *lenArySIPtr = 0;

   /*get past header*/
   tmpStr =
      (schar *)
      fgets(
         (char *) buffStr,
         def_lenBuff_fun15,
         tsvFILE
      );

   while(fgets((char *) buffStr, def_lenBuff_fun15, tsvFILE))
   { /*Loop: read in each sequence in the file*/

       /*************************************************\
       * Fun14 Sec05 Sub02:
       *   - get primer id
       \*************************************************/

       tmpStr = buffStr;
       seqStackST.idStr = buffStr;

       while(*tmpStr > 32)
          ++tmpStr;

       if(*tmpStr != '\t' && *tmpStr != ' ')
          goto fileErr_fun15_sec07_sub03; /*invalid line*/
       
       seqStackST.lenIdUL = tmpStr - buffStr;

       *tmpStr = '\0'; /*c-string for seqST*/
       ++tmpStr;

       /*************************************************\
       * Fun14 Sec05 Sub03:
       *   - check if primers are paired
       \*************************************************/

       while(*tmpStr < 33)
       { /*Loop: get off white space*/
          if(*tmpStr != '\t' && *tmpStr != ' ')
             goto fileErr_fun15_sec07_sub03; /*invalid*/

          ++tmpStr;
       } /*Loop: get off white space*/

       if(*tmpStr < 33)
          goto fileErr_fun15_sec07_sub03; /*invalid line*/

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
          goto fileErr_fun15_sec07_sub03; /*invalid line*/

       /*************************************************\
       * Fun14 Sec05 Sub04:
       *   - get the foward sequence
       \*************************************************/

       while(*tmpStr < 33)
       { /*Loop: get off white space*/
          if(*tmpStr != '\t' && *tmpStr != ' ')
             goto fileErr_fun15_sec07_sub03; /*invalid*/

          ++tmpStr;
       } /*Loop: get off white space*/

       if(*tmpStr < 33)
          goto fileErr_fun15_sec07_sub03; /*invalid line*/

       seqStackST.seqStr = tmpStr;

       while(*tmpStr > 32)
          ++tmpStr;

       if(*tmpStr != '\t' && *tmpStr != ' ')
          goto fileErr_fun15_sec07_sub03; /*invalid line*/

       *tmpStr = '\0'; /*convert sequence to c-string*/
       ++tmpStr;
       
       /*check if this is a blank entry*/
       if(seqStackST.seqStr[0] == '0')
       { /*If: no foward sequence*/
          retRefHeapAryST[*lenArySIPtr].mateSI = -1;
          retRefHeapAryST[*lenArySIPtr + 1].mateSI = -1;

          goto noForSeq_fun15_sec05_sub0x;
       } /*If: no foward sequence*/

       if(
             (seqStackST.seqStr[0] & ~ 32) == 'N'
          && (seqStackST.seqStr[1] & ~ 32) == 'A'
       ){ /*If: no foward sequence (NA)*/
          retRefHeapAryST[*lenArySIPtr].mateSI = -1;
          retRefHeapAryST[*lenArySIPtr +1].mateSI = -1;

          goto noForSeq_fun15_sec05_sub0x;  
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
          goto memErr_fun15_sec07_sub02;

       /*************************************************\
       * Fun14 Sec05 Sub06:
       *   - get the reverse sequence
       \*************************************************/

       ++(*lenArySIPtr);

       noForSeq_fun15_sec05_sub0x:;

       while(*tmpStr < 33)
       { /*Loop: get off white space*/
          if(*tmpStr != '\t' && *tmpStr != ' ')
             goto fileErr_fun15_sec07_sub03; /*invalid*/

          ++tmpStr;
       } /*Loop: get off white space*/

       if(*tmpStr < 33)
          goto fileErr_fun15_sec07_sub03; /*invalid line*/

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
          goto memErr_fun15_sec07_sub02;

       ++(*lenArySIPtr);
   } /*Loop: read in each sequence in the file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec06:
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
      goto memErr_fun15_sec07_sub02; /*memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec07:
   ^   - clean up
   ^   o fun15 sec06 sub01:
   ^     - success clean up
   ^   o fun15 sec06 sub02:
   ^     - memory error report
   ^   o fun15 sec06 sub03:
   ^     - file error report
   ^   o fun15 sec06 sub04:
   ^     - error clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun15 Sec07 Sub01:
   *   - success clean up
   \*****************************************************/

   *errSC = 0;

   fclose(tsvFILE);
   tsvFILE = 0;

   return retRefHeapAryST;

   /*****************************************************\
   * Fun15 Sec07 Sub02:
   *   - memory error report
   \*****************************************************/

   memErr_fun15_sec07_sub02:;

   *errSC = def_memErr_kmerFind;

   goto errCleanUp_fun15_sec07_sub04;

   /*****************************************************\
   * Fun15 Sec07 Sub03:
   *   - file error report
   \*****************************************************/

   fileErr_fun15_sec07_sub03:;

   *errSC = def_fileErr_kmerFind;

   goto errCleanUp_fun15_sec07_sub04;

   /*****************************************************\
   * Fun15 Sec07 Sub04:
   *   - error clean up
   \*****************************************************/

   errCleanUp_fun15_sec07_sub04:;

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
| Fun16: faToAry_refST_kmerFind
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
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec01:
   ^   - variable declerations 
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define def_lenBuff_fun16 4096
   signed char buffStr[def_lenBuff_fun16];
   ulong bytesUL = 0;
   ulong posUL = 0;
   schar newLineBl = 0;

   sint  numSeqSI = 0;

   uint longestPrimUI = 0; /*length of longest primer*/

   /*structures to work with*/
   struct refST_kmerFind *retRefHeapAryST = 0;

   FILE *faFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec02:
   ^   - initialize and check if can open input file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *lenArySIPtr = 0;

   faFILE = fopen((char *) faFileStr, "r");

   if(! faFILE)
      goto fileErr_fun16_sec06_sub03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec03:
   ^   - get number of sequences in file
   ^   o fun16 sec03 sub01:
   ^     - read in first part of file and start loop
   ^   o fun16 sec03 sub02:
   ^     - see if '>' for new line cases
   ^   o fun16 sec03 sub03:
   ^     - find next newline or end of buffer
   ^   o fun16 sec03 sub04:
   ^     - check for errors and move to start of file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun16 Sec03 Sub01:
   *   - read in first part of file and start loop
   \*****************************************************/

   bytesUL =
      fread(
         (char *) buffStr,
         sizeof(char),
         def_lenBuff_fun16 - 1,
         faFILE
      );

   buffStr[bytesUL] = '\0';
   posUL = 0;
   newLineBl = 1; /*assume I start on an new line*/

   while(bytesUL)
   { /*Loop: find number of sequences in file*/

      /**************************************************\
      * Fun16 Sec03 Sub02:
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
                     def_lenBuff_fun16 - 1,
                     faFILE
                  );

               buffStr[bytesUL] = '\0';
               posUL = 0;

               if(! bytesUL)
                  goto EOF_fun16_sec03_sub04;
            } /*If: I need to read in more buffer*/

            else
               ++posUL;
         } /*Loop: Get past white space*/

         if(buffStr[posUL] == '>')
            ++numSeqSI;

         newLineBl = 0;
      } /*If: I just finished an newline*/

      /**************************************************\
      * Fun16 Sec03 Sub03:
      *   - find next newline or end of buffer
      \**************************************************/

      posUL += endLine_ulCp(&buffStr[posUL]);

      if(buffStr[posUL] == '\0')
      { /*If: I need to get more file*/
         bytesUL =
            fread(
               (char *) buffStr,
               sizeof(char),
               def_lenBuff_fun16,
               faFILE
            );

         posUL = 0;
      } /*If: I need to get more file*/

      else
         newLineBl = 1;
   } /*Loop: find number of sequences in file*/

   /*****************************************************\
   * Fun16 Sec03 Sub04:
   *   - check for errors and move to start of file
   \*****************************************************/

   EOF_fun16_sec03_sub04:;

   if(*errSC)
   { /*If: I had an error of some kind*/
      if(*errSC & 64)
         goto memErr_fun16_sec06_sub02; /*memory error*/

      goto fileErr_fun16_sec06_sub03;
   } /*If: I had an error of some kind*/

   if(numSeqSI < 1)
      goto fileErr_fun16_sec06_sub03; /*no seq in file*/

   fseek(faFILE, 0, SEEK_SET);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec04:
   ^   - allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retRefHeapAryST =
      malloc(numSeqSI * sizeof(struct refST_kmerFind));

   if(! retRefHeapAryST)
      goto memErr_fun16_sec06_sub02;

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
      goto memErr_fun16_sec06_sub02;
   } /*If: I had an memroy error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec05:
   ^   - read in sequences
   ^   o fun16 sec05 sub01:
   ^     - add sequences to the reference array
   ^   o fun16 sec05 sub02:
   ^     - deal with errors or no sequences
   ^   o fun16 sec05 sub11:
   ^     - update the table window size values
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun16 Sec05 Sub01:
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
          goto memErr_fun16_sec06_sub02; /*memory error*/

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
          goto memErr_fun16_sec06_sub02;
   } /*Loop: read in each sequence in the file*/

   /*****************************************************\
   * Fun16 Sec05 Sub02:
   *   - deal with errors or no sequences
   \*****************************************************/

   if(*errSC)
   { /*If: I had an error of some kind*/
      if(*errSC & 64)
         goto memErr_fun16_sec06_sub02; /*memory error*/

      goto fileErr_fun16_sec06_sub03;
   } /*If: I had an error of some kind*/

   if(numSeqSI < 1)
      goto fileErr_fun16_sec06_sub03; /*no seq in file*/

   /*****************************************************\
   * Fun16 Sec05 Sub11:
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
      goto memErr_fun16_sec06_sub02; /*memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec06:
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
   * Fun16 Sec06 Sub01:
   *   - success clean up
   \*****************************************************/

   *errSC = 0;

   fclose(faFILE);
   faFILE = 0;

   return retRefHeapAryST;

   /*****************************************************\
   * Fun16 Sec06 Sub02:
   *   - memory error report
   \*****************************************************/

   memErr_fun16_sec06_sub02:;

   *errSC = def_memErr_kmerFind;

   goto errCleanUp_fun16_sec06_sub04;

   /*****************************************************\
   * Fun16 Sec06 Sub03:
   *   - file error report
   \*****************************************************/

   fileErr_fun16_sec06_sub03:;

   *errSC = def_fileErr_kmerFind;

   goto errCleanUp_fun16_sec06_sub04;

   /*****************************************************\
   * Fun16 Sec06 Sub04:
   *   - error clean up
   \*****************************************************/

   errCleanUp_fun16_sec06_sub04:;

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
| Fun17: nextSeqChunk_tblST_kmerFind
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
   ' Fun17 TOC:
   '   o fun17 sec01:
   '     - variable declerations
   '   o fun17 sec02:
   '     - handle first time adding in sequence cases
   '   o fun17 sec03:
   '     - remove the old kmers from the table
   '   o fun17 sec04:
   '     - move keep kmers to start of kmer array
   '   o fun17 sec05:
   '     - add new kmers to table
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint uiNt = 0;
   uint endWindowUI = 0;
   uint dupNtUI = 0;
   schar lastWinBl = 0;

   uchar ntUC = 0; /*holds one nucleotide*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec02:
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

      while(tblSTPtr->seqPosUL < (ulong) endWindowUI)
      { /*Loop: build the first kmer*/
         ntUC =
            (uchar) 
            tblSTPtr->seqSTPtr->seqStr[
               tblSTPtr->seqPosUL
            ];

         ntUC = alnNtToBit_kmerFind[ ntUC ];
         tblSTPtr->lastKmerUL <<=def_bitsPerKmer_kmerFind;
         tblSTPtr->lastKmerUL |= ntUC;
         tblSTPtr->lastKmerUL &= tblSTPtr->kmerMaskUL;

         ++tblSTPtr->seqPosUL;

         if(ntUC < def_anonNt_kmerFind)
            ++tblSTPtr->lenLastKmerUL;
          else
            tblSTPtr->lenLastKmerUL = 0;
      } /*Loop: build the first kmer*/

      endWindowUI += tblSTPtr->ntInWinUI;
      endWindowUI -= tblSTPtr->lenKmerUC;
      dupNtUI = 0;   /*first base in window*/
      goto firstKmers_fun17_sec04;
   } /*If: this is the first window*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec03:
   ^   - remove the old kmers from the table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      uiNt = 0;
      uiNt < (tblSTPtr)->rmNtUI;
      ++uiNt
   ){ /*Loop: remove discarded kmers from table*/
      tblSTPtr->lastKmerUL =
         (ulong) tblSTPtr->kmerArySI[uiNt];

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
   ^ Fun17 Sec04:
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
   ^ Fun17 Sec05:
   ^   - add new kmers to table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   firstKmers_fun17_sec04:; /*first time for sequence*/

   if(endWindowUI > tblSTPtr->seqSTPtr->lenSeqUL)
   { /*If: this is the last window*/
      endWindowUI = tblSTPtr->seqSTPtr->lenSeqUL;
      lastWinBl = 1;
   } /*If: this is the last window*/

   while(tblSTPtr->seqPosUL < endWindowUI)
   { /*Loop: add kmers in*/
      ntUC =
         (uchar) 
         tblSTPtr->seqSTPtr->seqStr[tblSTPtr->seqPosUL];

      ntUC = alnNtToBit_kmerFind[ ntUC ];
      tblSTPtr->lastKmerUL <<= def_bitsPerKmer_kmerFind;
      tblSTPtr->lastKmerUL |= ntUC;
      tblSTPtr->lastKmerUL &= tblSTPtr->kmerMaskUL;

      /*add kmer to kmer table counts if is complete*/

      if(ntUC < def_anonNt_kmerFind)
      { /*If: no anymous bases or errors*/
         ++tblSTPtr->lenLastKmerUL;

         if(
               (slong) tblSTPtr->lenLastKmerUL
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
| Fun18: forCntMatchs_kmerFind
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
   slong tmpSL = 0;
   slong slKmer = 0;
   slong numMatchSL = 0;

   for(
      slKmer = 0;
      slKmer < (slong) refST_kmerFindPtr->lenAryUI;
      ++slKmer
   ){ /*Loop: count number of matching kmers*/
      tmpSL = refST_kmerFindPtr->forKmerArySI[ slKmer ];

      if(tmpSL < 0)
         break; /*if at end of kmers*/

      tmpSL =
         min_genMath(
            (sint) refST_kmerFindPtr->forRepAryUI[tmpSL],
            tblST_kmerFindPtr->tblSI[ tmpSL ]
         );

      numMatchSL += tmpSL;
   } /*Loop: count number of matching kmers*/

   return numMatchSL;
} /*forCntMatchs_kmerFind*/ 

/*-------------------------------------------------------\
| Fun19: revCntMatchs_kmerFind
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
   slong tmpSL = 0;
   slong slKmer = 0;
   slong numMatchSL = 0;

   for(
      slKmer = 0;
      slKmer < (slong) refST_kmerFindPtr->lenAryUI;
      ++slKmer
   ){ /*Loop: count number of matching kmers*/
      tmpSL =
         refST_kmerFindPtr->revKmerArySI[ slKmer ];

      if(tmpSL < 0)
         break; /*if at end of kmers*/

      tmpSL =
         min_genMath(
            (sint) refST_kmerFindPtr->revRepAryUI[tmpSL],
            tblST_kmerFindPtr->tblSI[ tmpSL ]
         );

      numMatchSL += tmpSL;
   } /*Loop: count number of matching kmers*/

   return numMatchSL;
} /*revCntMatchs_kmerFind*/

/*-------------------------------------------------------\
| Fun20: matchCheck_kmerFind
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
   schar retBl = 0; 
   schar revBestBl = 0; 
   ulong forKmerCntUL = 0; 
   ulong revKmerCntUL = 0; 
   ulong bestKmerCntUL = 0; 
   
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
| Fun21: findRefInChunk_kmerFind
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
   ' Fun21 TOC:
   '   - finds an sequence pattern in an sam entry
   '   o fun21 sec01:
   '     - variable declerations
   '   o fun21 sec02:
   '     - initialize & see if enough kmers for alignment
   '   o fun21 sec03:
   '     - prepare for alignemnt (if passed kmer check)
   '   o fun21 sec04:
   '     - do alignment and check if passes min score
   '   o fun21 sec05:
   '     - return the answer
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar matchBl = 0;
   float percScoreF = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec02:
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
   ^ Fun21 Sec03:
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
      ^ Fun21 Sec04:
      ^   - do alignment and check if passes min score
      ^   o fun21 sec04 sub01:
      ^     - do the alignment
      ^   o fun21 sec04 sub02:
      ^     - check if it passes the alignment
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      /**************************************************\
      * Fun21 Sec04 Sub01:
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
      * Fun21 Sec04 Sub02:
      *   - check if it passes the alignment
      \**************************************************/

      matchBl &= ( -(percScoreF >= minPercScoreF) );
      ++tblSTPtr->seqSTPtr->endAlnUL;
   } /*If: I had enough kmers to do an alignment*/

   else
     matchBl = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec05:
   ^   - return the answer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return matchBl;
} /*findRefInChunk_kmerFind*/

/*-------------------------------------------------------\
| Fun22: waterFindPrims_kmerFind
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
   ^ Fun22 TOC:
   '   o fun22 sec01:
   '     - varaible declerations
   '   o fun22 sec02:
   '     - assign sequence to table
   '   o fun22 sec03:
   '     - check sequence for spacers
   '   o fun22 sec04:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec01:
   ^   - varaible declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uchar retUC = def_noPrim_kmerFind; /*return*/

   uint uiPrim = 0;
   float percScoreF = 0;

   /*alignemnt variables; I keep*/
   slong scoreSL = 0;
   ulong qryStartUL = 0;
   ulong qryEndUL = 0;

   /*alignemnt variables; I ingore*/
   ulong refStartUL = 0;
   ulong refEndUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec02:
   ^   - convert to sequence to index
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   seqToIndex_alnSet(seqSTPtr->seqStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec03:
   ^   - check sequence for primers
   ^   o fun22 sec03 sub01:
   ^     - start primer loop
   ^   o fun22 sec03 sub02:
   ^     - foward alignment of primer
   ^   o fun22 sec03 sub03:
   ^     - reverse complement alignment of primer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun22 Sec03 Sub01:
   *   - start primer loop
   \*****************************************************/

   for(
      uiPrim = 0;
      uiPrim < lenRefAryUI;
      ++uiPrim
   ){ /*Loop: detect primers in each chunk*/

      /**************************************************\
      * Fun22 Sec03 Sub02:
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
      * Fun22 Sec03 Sub03:
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
   ^ Fun22 Sec04:
   ^   - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   indexToSeq_alnSet(seqSTPtr->seqStr);
   return (schar) retUC;
} /*waterFindPrims_kmerFind*/

/*-------------------------------------------------------\
| Fun23: fxFindPrims_kmerFind
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
   ^ Fun23 TOC:
   '   o fun23 sec01:
   '     - varaible declerations
   '   o fun23 sec02:
   '     - assign sequence to table
   '   o fun23 sec03:
   '     - check sequence for spacers
   '   o fun23 sec04:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec01:
   ^   - varaible declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0; /*error messages*/
   schar retSC = def_noPrim_kmerFind; /*return*/

   uchar matchBl = 0;
   signed char firstTimeBl = 1;

   uint uiPrim = 0;

   /*alignemnt variables; I keep*/
   slong scoreSL = 0;
   ulong qryStartUL = 0;
   ulong qryEndUL = 0;

   /*alignemnt variables; I ingore*/
   ulong refStartUL = 0;
   ulong refEndUL = 0;

   /*to keep the old assigned sequence*/
   struct seqST *oldSeqST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec02:
   ^   - check positions and assign sequence to table
   ^   o fun23 sec02 sub01:
   ^     - see if i have an direct repeat region
   ^   o fun23 sec02 sub02:
   ^     - add the sequence to the kmer table
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun23 Sec02 Sub01:
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
   * Fun23 Sec02 Sub02:
   *   - add the sequence to the kmer table
   \*****************************************************/

   retSC = def_noPrim_kmerFind; /*allows detection*/

   oldSeqST = tblSTPtr->seqSTPtr;
   tblSTPtr->seqSTPtr = seqSTPtr;

   blank_tblST_kmerFind(
      tblSTPtr,
      0          /*I do not want to blank the sequence*/
   );

   /*prepare for waterman alignment*/
   seqToIndex_alnSet(tblSTPtr->seqSTPtr->seqStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec03:
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
   ^ Fun23 Sec04:
   ^   - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*reset table to original sequences*/
   indexToSeq_alnSet(tblSTPtr->seqSTPtr->seqStr);
   tblSTPtr->seqSTPtr = oldSeqST;

   return (schar) retSC;
} /*fxFindPrims_kmerFind*/

/*-------------------------------------------------------\
| Fun24: phit_kmerFind
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
   sint siRef = 0;
   sint siMate = 0;
   schar *oldRefStr = 0;
   schar *oldSeqStr = 0;
   schar oldRefBreakSC = 0;
   schar oldSeqBreakSC = 0;

   uchar skip1stCharUC = 0;

   oldSeqStr = seqSTPtr->idStr;
   while(*oldSeqStr++ > 32);
   --oldSeqStr;
   oldSeqBreakSC = *oldSeqStr;
   *oldSeqStr = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun24 Sec02:
   ^   - print out mapped primers
   ^   o fun24 sec02 sub01:
   ^     - start loop and check primers that mapped
   ^   o fun24 sec02 sub02:
   ^     - primer mate (pairing) stats printout
   ^   o fun24 sec02 sub03:
   ^     - non-primer mate (no pairing) stats printout
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun24 Sec02 Sub01:
   *   - start loop and check primers that mapped
   \*****************************************************/

   siRef = 0;
   skip1stCharUC = *refAryST[siRef].forSeqST->idStr =='>';

   while(siRef < numRefsSI)
   { /*Loop: print out hits*/
      if(codeAryUI[siRef] == 0)
         goto nextRef_fun24_sec02_sub04;

      if(dirArySC[siRef] == 'N')
         goto nextRef_fun24_sec02_sub04;

      /**************************************************\
      * Fun24 Sec02 Sub02:
      *   - primer mate (pairing) stats printout
      \**************************************************/

      if(refAryST[siRef].mateSI >= 0)
      { /*If: I have a mate primer*/
         siMate = refAryST[siRef].mateSI;

         if(codeAryUI[siMate] == 0)
         { /*If: the mate did not map*/
            ++siRef; /*move past mate*/
            goto nextRef_fun24_sec02_sub04;
         } /*If: the mate did not map*/

         if(dirArySC[siMate] == 'N')
         { /*If: the mate did not map*/
            ++siRef; /*move past mate*/
            goto nextRef_fun24_sec02_sub04;
         } /*If: the mate did not map*/

         if(dirArySC[siRef] == dirArySC[siMate])
         { /*If: the mate was in the same direction*/
            ++siRef; /*move past mate*/
            goto nextRef_fun24_sec02_sub04;
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
             seqStartAryUL[siRef],
             seqEndAryUL[siRef]
         ); /*forward primer mapping coordinates*/

         fprintf(
             (FILE *) outFILE,
             "\t%c\t%lu\t%lu",
             dirArySC[siMate],
             seqStartAryUL[siMate],
             seqEndAryUL[siMate]
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
             primStartAryUL[siRef],
             primEndAryUL[siRef]
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
             primStartAryUL[siMate],
             primEndAryUL[siMate]
         ); /*reverse primer mapping stats*/

         *oldRefStr = oldRefBreakSC;
         ++siRef; /*move past mate*/
      } /*If: I have a mate primer*/

      /**************************************************\
      * Fun24 Sec02 Sub03:
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
             seqStartAryUL[siRef],
             seqEndAryUL[siRef]
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
             primStartAryUL[siRef],   /*start of map*/
             primEndAryUL[siRef]      /*end of map*/
         ); /*foward primer mapping stats*/

         *oldRefStr = oldRefBreakSC;

         fprintf(
             (FILE *) outFILE,
             "\tNA\tNA\tNA\tNA\tNA\n"
         ); /*reverse primer mapping stats*/
      } /*Else: I have no mate primers*/

      /**************************************************\
      * Fun24 Sec02 Sub04:
      *   - move to the next reference
      \**************************************************/

      nextRef_fun24_sec02_sub04:;

      ++siRef;
   } /*Loop: print out hits*/

   *oldSeqStr = oldSeqBreakSC;
} /*phit_kmerFInd*/

/*-------------------------------------------------------\
| Fun25: pHeaderHit_kmerFind
|    - prints header for phit_kmerFind (fun24)
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
