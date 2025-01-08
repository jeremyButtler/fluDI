/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' tbCon SOF: Start Of File
'  - holds functions to make a quick consensus, more like
'    ivar. this is a bit memory intensive,
'    (40 + # insertions * (17) + bases in insertions)
'    bytes 
'   o header:
'     - Included libraries
'   o .h st01: insBase
'     - structure to hold an list of insetions for a
'       single position
'   o .h st02: conNt_tbCon
'     - structure to hold an single base in a sequence
'   o .h st03: set_tbCon
'     - holds the settings for tbCon
'   o fun01: initInsBase
'     - initializes an initBase structure
'   o fun02: mkIns_tbCon
'     - makes and initializes an ins_tbCon structure
'   o fun03: freeHeap_ins_tbCon
'     - Frees an heap allocated insBase structure
'   o fun04: blank_conNt_tbCon
'     - sets all variables in conNt_tbCon structure to 0
'   o fun05: init_conNt_tbCon
'     - initializes a conNt_tbCon structure
'   o fun06: mk_conBase_tbCon
'     - Makes an initialized, heap allocated conNt_tbCon
'       structure
'   o fun07: freeStack_conNt_tbCon
'     - frees heap allocated variables in a conNt_tbCon
'       structure (also initializes structure)
'   o fun08: freeHeap_conNt_tbCon
'     - frees a conNt_tbCon structure
'   o fun09: freeHeapAry_conNt_tbCon
'     - frees array of heap allocated conNt_tbCon structs
'   o fun10: init_set_tbCon
'     - initialize a set_tbCon struct to default settings
'   o fun11: freeStack_set_tbCon
'     - frees variables inside a set_tbCon struct
'   o fun12: addRead_tbCon
'     - adds read to a conNt_tbCon array
'   o fun13: collapse_tbCon
'     - Collapses a conNt_tbCon array into an array of
'       samEntry structures
'   o fun14: noFragCollapse_tbCon
'     - collapses a conNt_tbCon struct array into a single
'       samEntry struct (low read depth is masked)
'   o fun15: pvar_tbCon
'     - print entries in a conNt_tbCon structure array
'       that are above the minimum read depth
'   o license:
'     - licensing for this code (public domain / mit)
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

#include <stdio.h>

#include "tbCon.h"

#include "../genLib/ulCp.h"
#include "../genLib/numToStr.h"
#include "samEntry.h"

/*.h files only (no .c file or not using .c file*/
#include "../genLib/dataTypeShortHand.h"
#include "../genLib/genMath.h" /*only using .h max macro*/
#include "tbConDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries
!   o .c  #include "../genLib/base10str.h"
!   o .c  #include "../genLib/strAry.h"
!   o .h  #include "ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: initInsBase
|   - Initializes an initBase structure
| Input:
|   - insSTPtr:
|     o A Pointer to an insBase structure to initialize
| Output:
|   - Modifies:
|     o All values in insBase to be 0
\-------------------------------------------------------*/
void
init_insBase_tbCon(
   struct ins_tbCon *insSTPtr
){
   insSTPtr->insStr = 0;
   insSTPtr->lenInsSI = 0;
   insSTPtr->numHitsSI = 0;
   insSTPtr->nextIns = 0;
} /*initInsBase*/

/*-------------------------------------------------------\
| Fun02: mkIns_tbCon
|   - makes and initializes an ins_tbCon structure
| Input:
| Output:
|   - Returns:
|     o pointer to initialized ins_tbCon structure
|     o 0 for memory errors
\-------------------------------------------------------*/
struct ins_tbCon *
mkIns_tbCon(
){
   struct ins_tbCon *retMacST =
      malloc(sizeof(struct ins_tbCon));

   if(retMacST)
      init_insBase_tbCon(retMacST);

   return retMacST;
} /*mkIns_tbCon*/

/*-------------------------------------------------------\
| Fun03: freeHeap_ins_tbCon
|   - frees heap allocated insBase structure
| Input:
|   - insSTPtr:
|     o pointer to insBase structure to free
| Output:
|   - Frees:
|     o insSTPtr
|   - Returns:
|     o next insBase structure (nextIns) in insBase list
\-------------------------------------------------------*/
struct ins_tbCon *
freeHeap_ins_tbCon(
   struct ins_tbCon *insSTPtr
){
   struct ins_tbCon *nextIns = 0;
 
   if(insSTPtr)
   { /*If: have structure to free*/
      nextIns = insSTPtr->nextIns;
 
      if(insSTPtr->insStr)
         free(insSTPtr->insStr);
 
      free(insSTPtr);
   } /*If: have structure to free*/
 
   return nextIns;
} /*freeHeap_ins_tbCon*/

/*-------------------------------------------------------\
| Fun04: blank_conNt_tbCon
|   - sets all variables in conNt_tbCon structure to 0
| Input:
|   - conNtSTPtr:
|     o pointer to conNt_tbCon structure to blank
| Output:
|   - Modifies:
|     o all non-pointer values in conNtSTPtr to be 0
\-------------------------------------------------------*/
void
blank_conNt_tbCon(
   struct conNt_tbCon *conNtSTPtr
){
   conNtSTPtr->numASI = 0;
   conNtSTPtr->numTSI = 0;
   conNtSTPtr->numGSI = 0;
   conNtSTPtr->numCSI = 0;
   conNtSTPtr->numDelSI = 0;
   
   conNtSTPtr->numInsSI = 0;

   conNtSTPtr->totalNtSI = 0;
   conNtSTPtr->ntKeptSI = 0;
} /*blank_conNt_tbCon*/

/*-------------------------------------------------------\
| Fun05: init_conNt_tbCon
|   - initializes a conNt_tbCon structure
| Input:
|   - conNtSTPtr:
|     o pointer to an conNt_tbCon structure to initialize
| Output:
|   - Modifies:
|     o all values (including pointers) in conNtSTPtr to
|       be 0
\-------------------------------------------------------*/
void
init_conNt_tbCon(
   struct conNt_tbCon *conNtSTPtr
){
   blank_conNt_tbCon(conNtSTPtr);
   conNtSTPtr->insList = 0;
} /*init_conNt_tbCon*/

/*-------------------------------------------------------\
| Fun06: mk_conBase_tbCon
|   - makes initialized, heap allocated conNt_tbCon struct
| Input:
| Output:
|   - Returns:
|     o pointer to an initialized conNt_tbCon structure
|     o 0 for memory errors
\-------------------------------------------------------*/
struct conNt_tbCon *
mk_conBase_tbCon(
){
   struct conNt_tbCon *retST =
      malloc(sizeof(struct conNt_tbCon));
   
   if(retST)
      init_conNt_tbCon(retST);
   
   return retST;
} /*mk_conBase_tbCon*/

/*-------------------------------------------------------\
| Fun07: freeStack_conNt_tbCon
|   - frees heap allocated variables in a conNt_tbCon
|     structure (also initializes structure)
| Input:
|   - conNtSTPtr:
|     o pointer to conNt_tbCon struct to free variables
| Output:
|   - Frees:
|     o frees the conNt_tbConST->insList 
|   - Sets:
|     o all values in conNt_tbConSTPTr to 0
\-------------------------------------------------------*/
void
freeStack_conNt_tbCon(
   struct conNt_tbCon *conNtSTPtr
){
  if(conNtSTPtr)
  { /*If: have structure to free*/
     while(conNtSTPtr->insList)
     { /*Loop: Free insertions*/
        conNtSTPtr->insList =
            freeHeap_ins_tbCon(conNtSTPtr->insList);
     } /*Loop: Free insertions*/
     
     blank_conNt_tbCon(conNtSTPtr);
  } /*If: have structure to free*/
} /*freeStack_conNt_tbCon*/

/*-------------------------------------------------------\
| Fun08: freeHeap_conNt_tbCon
|   - frees a conNt_tbCon structure
| Input:
|   - conNtSTPtr:
|     o pointer toa conNt_tbCon structure to free
| Output:
|   - Frees:
|     o frees conNt_tbConST
\-------------------------------------------------------*/
void
freeHeap_conNt_tbCon(
   struct conNt_tbCon *conNtSTPtr
){
  if(conNtSTPtr)
  { /*If: have structure to free*/
     freeStack_conNt_tbCon(conNtSTPtr);
     free(conNtSTPtr);
  } /*If: have structure to free*/
} /*freeHeap_conNt_tbCon*/

/*-------------------------------------------------------\
| Fun09: freeHeapAry_conNt_tbCon
|   - frees array of heap allocated conNt_tbCon structs
| Input:
|   - conNtAryST:
|     o pointer conNt_tbCon struct array to free
|   - lenArySI:
|     o number of conNt_tbCon structurs in array
| Output:
|   - Frees:
|     o array in conNtSTPtr
\-------------------------------------------------------*/
void
freeHeapAry_conNt_tbCon(
   struct conNt_tbCon *conNtAryST,
   signed int lenArySI
){
  int siIndex = 0;

  if(conNtAryST)
  { /*If: have structures to free*/
    for(
       siIndex = 0;
       siIndex < lenArySI;
       ++siIndex
    ) freeStack_conNt_tbCon(&conNtAryST[siIndex]);
    
    free(conNtAryST);
  } /*If: have structures to free*/
} /*freeHeap_conNt_tbCon*/

/*-------------------------------------------------------\
| Fun10: init_set_tbCon
|   - initialize a set_tbCon struct to default settings
| Input:
|   - setSTPtr:
|     o Pointer to the set_tbCon structure with settings
|       for tbCon to set to defaults
| Output:
|   - Modifies:
|     o set_tbConSTPtr to default settings
\-------------------------------------------------------*/
void
init_set_tbCon(
   struct set_tbCon *setSTPtr
){
   setSTPtr->minQSI = def_minNtQ_tbConDefs;
   setSTPtr->minInsQSI = def_minInsQ_tbConDefs;
   setSTPtr->minMapqUC = def_minMapq_tbConDefs;

   /*Min percent/depth/length to keep an match/mutation*/
   setSTPtr->minLenSI = def_minLen_tbConDefs;
   setSTPtr->minDepthSI = def_minDepth_tbConDefs;
   setSTPtr->minPercSnpF = def_minSnpPerc_tbConDefs;
   setSTPtr->minPercDelF = def_minDelPerc_tbConDefs;
   setSTPtr->minPercInsF = def_minInsPerc_tbConDefs;

   /*To mask low quality/support cosnesus bases with*/
   setSTPtr->maskSC = def_mask_tbConDefs;

   setSTPtr->minPrintDepthSI = def_minPrintDepth_tbConDefs;

   setSTPtr->printMinSupSnpF =
      def_minSnpPrintPerc_tbConDefs;

   setSTPtr->printMinSupInsF =
      def_minInsPrintPerc_tbConDefs;

   setSTPtr->printMinSupDelF =
      def_minDelPrintPerc_tbConDefs;
} /*init_set_tbCon*/

/*-------------------------------------------------------\
| Fun11: freeStack_set_tbCon
|   - frees variables inside a set_tbCon struct
| Input:
|   - setSTPtr:
|     o pointer to the set_tbCon struct
| Output:
|   - Nothing:
|     o here for future use
\-------------------------------------------------------*/
void
freeStack_set_tbCon(
   struct set_tbCon *setSTPtr
){
   setSTPtr = setSTPtr;
} /*freeStack_set_tbCon*/

/*-------------------------------------------------------\
| Fun12: addRead_tbCon
|   - adds read to a conNt_tbCon array
| Input:
|   - sameEntrySTPtr:
|     o pointer to a samEntry struct with read to add
|   - conNtAryST:
|     o pointer to conNt_tbCon struct array to add read to
|     o the array is resize if to small
|   - lenRefUI:
|     o length of the reference sequence; use 0 for unkown
|     o is updated as conNt_tbConSTAry expands
|   - settings:
|     o pointer to set_tbCon struct with settings
| Output:
|   - Modifies:
|     o conNtAryST to have new read
|     o lenRefUI to have new conNtAryST size if conNtAryST
|       is resized
|   - Returns:
|     o for no errors
|     o def_header_tbConDefs for sam file header entry
|     o def_noMap_tbConDefs for an unmapped read
|     o def_noSeq_tbConDefs for a read missing a sequence
|     o def_memErr_tbConDefs for memory errors
\-------------------------------------------------------*/
signed char
addRead_tbCon(
   struct samEntry *samSTPtr, /*read to add to consensus*/
   struct conNt_tbCon *conNtAryST[], /*consensus array*/
   unsigned int *lenRefUI,       /*length of reference*/
   struct set_tbCon *settings
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun12 TOC: addRead_tbCon
   '   - Adds the read into a conNt_tbCon array
   '   o fun12 sec01:
   '     - variable declerations
   '   o fun12 sec02:
   '     - check read and allocate/reallocate memory
   '   o fun12 sec03:
   '     - Add bases to the consensus
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;        /*err reporting at end*/

   /*For loops*/
   uint uiCig = 0;         /*for iterating throug cigars*/
   uint uiBase = 0;        /*for iterating through bases*/
   uint uiRef = 0;         /*position at in the reference*/
   uint endCigUI = 0;      /*end of cigar entry*/

   /*For adding insertions to the consensus*/
   sint siIns = 0;         /*position at in insertion*/
   sint lenInsSI = 0;      /*length of found ins*/
   schar *insHeapStr  = 0; /*memory to hold insertion*/

   /*for adding insertions or looking for insertions*/
   struct ins_tbCon *insST = 0;
   struct ins_tbCon *lastInsST = 0;

   struct conNt_tbCon *baseST = 0; /*for reallocs*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec02:
   ^   - check read and allocate/reallocate memory
   ^   o fun12 sec02 sub01:
   ^     - check if read can be added to consensus
   ^   o fun12 sec02 sub02:
   ^     - set up consnesus array (first time only)
   ^   o fun12 sec02 sub03:
   ^     - increase consnesus array size (if needed)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun12 Sec02 Sub01:
   *   - check if read can be added to consensus
   \*****************************************************/

   if(samSTPtr->extraStr[0] =='\0')
      goto headerErr_fun12_sec04_sub03;

   if(samSTPtr->flagUS & 4)
      goto noMapErr_fun12_sec04_sub04;

   if(samSTPtr->seqStr == 0)
      goto noSeqErr_fun12_sec04_sub05;

   if(samSTPtr->seqStr[0] == '*')
      goto noSeqErr_fun12_sec04_sub05;

   if(samSTPtr->mapqUC < settings->minMapqUC)
      goto lowMapErr_fun12_sec04_sub06;

   /*****************************************************\
   * Fun12 Sec02 Sub02:
   *   - set up consnesus array (first time only)
   \*****************************************************/

   /*First time*/
   if(! *conNtAryST)
   { /*If: First round*/
      if(*lenRefUI < samSTPtr->refEndUI)
         *lenRefUI = samSTPtr->refEndUI + 128;

      *conNtAryST =
         malloc(*lenRefUI * sizeof(struct conNt_tbCon));

      if(! *conNtAryST)
         goto memErr_fun12_sec04_sub02;

      /*Initialize all the strutures*/
      for(
         uiBase = 0;
         uiBase < *lenRefUI;
         ++uiBase
      ) init_conNt_tbCon(&(*conNtAryST)[uiBase]);
   } /*If: First round*/

   /*****************************************************\
   * Fun12 Sec02 Sub03:
   *   - increase consnesus array size (if needed)
   \*****************************************************/

   else if(*lenRefUI < samSTPtr->refEndUI)
   { /*If: need to add more bases*/
      uiBase = *lenRefUI;

      /*add some extra memory to reduce future reallocs*/
      *lenRefUI = samSTPtr->refEndUI + 128;

      /*TODO: having issue here for repeated consensuses*/
      baseST =
         realloc(
           *conNtAryST,
           *lenRefUI * sizeof(struct conNt_tbCon)
         ); /*Add more memory for the bases*/

      if(! baseST)
         goto memErr_fun12_sec04_sub02;

      *conNtAryST = baseST; 

      while(uiBase < *lenRefUI)
      { /*Loop: Initialize the new structures*/
          init_conNt_tbCon(&(*conNtAryST)[uiBase]);
          ++uiBase;
      } /*Loop: Initialize the new structures*/
   } /*If: need to add more bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec03:
   ^   - add bases to consensus
   ^   o fun12 sec03 sub01:
   ^     - setup and start loop
   ^   o fun12 sec03 sub02:
   ^     - check matches and snp cases
   ^   o fun12 sec03 sub03:
   ^     - check deletion cases
   ^   o fun12 sec03 sub04:
   ^     - check the insertion cases
   ^   o fun12 sec03 sub05:
   ^     - handle softmasking
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun12 Sec03 Sub01:
   *   - setup and start loop
   \*****************************************************/

   uiBase = 0;
   uiRef = samSTPtr->refStartUI;

   for(
      uiCig = 0;
      uiCig < samSTPtr->lenCigUI;
      ++uiCig
   ){ /*Loop: add bases to the consensus*/

      switch(samSTPtr->cigTypeStr[uiCig])
      { /*Switch: Check the cigar entry type*/

         /***********************************************\
         * Fun12 Sec03 Sub02:
         *   - check matches and snp cases
         *   o fun12 sec03 sub02 cat01:
         *     - add snps/matches to total & check q-score
         *   o fun12 sec03 sub02 cat02:
         *     - passed qc; add to nucleotide count
         *   o fun12 sec03 sub02 cat03:
         *     - move to next base in match/snp
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun12 Sec03 Sub02 Cat01:
         +   - add snps/matches to total and check q-score
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'M':
         case 'X':
         case '=':
         /*Case: match or snp*/

            endCigUI = uiBase + samSTPtr->cigArySI[uiCig];

            while(uiBase < endCigUI)
            { /*Loop: add matchs/snps to consensus*/
               ++(*conNtAryST)[uiRef].totalNtSI;

               if(
                     samSTPtr->qStr[0] != '*'
                  && samSTPtr->qStr[1] != '\0'
               ){ /*If: I have a q-score entry*/
                  if(
                       samSTPtr->qStr[uiBase]
                         - def_adjQ_samEntry
                     < settings->minQSI
                  ){ /*If: base is beneath min q*/
                     ++uiBase;
                     ++uiRef;
                     continue;
                  } /*If: base is beneath min q*/
               } /*If: I have a q-score entry*/

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun12 Sec03 Sub02 Cat02:
               +   - passed qc; add to nucleotide count
               \++++++++++++++++++++++++++++++++++++++++*/

               /*Incurmenting total base count here so
               `  masked positions are still kept
               */
               ++(*conNtAryST)[uiRef].ntKeptSI;

               switch(samSTPtr->seqStr[uiBase] & ~32)
               { /*Switch: Check the base type*/
                  case 'A':
                     ++(*conNtAryST)[uiRef].numASI;
                     break;

                  case 'T':
                     ++(*conNtAryST)[uiRef].numTSI;
                     break;

                  case 'C':
                     ++(*conNtAryST)[uiRef].numCSI;
                     break;

                  case 'G':
                     ++(*conNtAryST)[uiRef].numGSI;
                     break;
               } /*Switch: Check the base type*/

               /*++++++++++++++++++++++++++++++++++++++++\
               + Fun12 Sec03 Sub02 Cat03:
               +   - move to next base in match/snp
               \++++++++++++++++++++++++++++++++++++++++*/

               ++uiBase;
               ++uiRef;
            } /*Loop: Add matchs/snps to consensus*/

            break;
         /*Case: match or snp*/

         /***********************************************\
         * Fun12 Sec03 Sub03:
         *   - check deletion cases
         \***********************************************/

         case 'D':
         /*Case: Deletion*/ 
            endCigUI = uiRef + samSTPtr->cigArySI[uiCig];

            while(uiRef < endCigUI)
            { /*Loop: add in deletions*/
               ++(*conNtAryST)[uiRef].totalNtSI;
               ++(*conNtAryST)[uiRef].ntKeptSI;
               ++(*conNtAryST)[uiRef].numDelSI;
               ++uiRef;
            } /*Loop: add in deletions*/

            break;
         /*Case: Deletion*/ 

         /***********************************************\
         * Fun12 Sec03 Sub04:
         *   - check insertion cases
         *   o fun12 sec02 sub04 cat01:
         *     - set up and allocate memory for ins
         *   o fun12 sec02 sub04 cat02:
         *     - copy insertion bases with min q-score
         *   o fun12 sec02 sub04 cat03:
         *     - check if kept any insertion bases
         *   o fun12 sec02 sub04 cat04:
         *     - see if insertion sequence already added
         *   o fun12 sec02 sub04 cat05:
         *     - insertion is new; add to insertion list
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun12 Sec02 Sub04 Cat01:
         +   - set up and allocate memory for ins
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'I':
         /*Case: Insertions*/
            lenInsSI = samSTPtr->cigArySI[uiCig];

            endCigUI =
               uiBase + samSTPtr->cigArySI[uiCig];

            insHeapStr =
               malloc((lenInsSI + 1) * sizeof(char));

            if(! insHeapStr)
               goto memErr_fun12_sec04_sub02;

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun12 Sec02 Sub04 Cat02:
            +   - copy insertion bases with min q-score
            \+++++++++++++++++++++++++++++++++++++++++++*/

            siIns = 0;

            /*Get the total Q-score for the insertion*/
            while(uiBase < endCigUI)
            { /*Loop: Copy the insertion*/

              if(
                     samSTPtr->qStr[0] != '*'
                  && samSTPtr->qStr[1] != '\0'
              ){ /*If: have q-score entry*/

                 if(
                      samSTPtr->qStr[uiBase]
                         - def_adjQ_samEntry
                    < settings->minInsQSI
                 ){ /*If: This insertion is low quality*/
                    ++uiBase;
                    continue;
                 } /*If: This insertion is low quality*/
              } /*If: have q-score entry*/

              insHeapStr[siIns] =
                 samSTPtr->seqStr[uiBase];

              ++uiBase;
              ++siIns;
            } /*Loop: Copy the insertion*/
 
            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun12 Sec02 Sub04 Cat03:
            +   - check if kept any insertion bases
            \+++++++++++++++++++++++++++++++++++++++++++*/

            if(siIns == 0)
            { /*If: I discarded the insertion*/
               free(insHeapStr);
               insHeapStr = 0;
               break;
            } /*If: I discarded the insertion*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun12 Sec02 Sub04 Cat04:
            +   - see if insertion sequence already added
            \+++++++++++++++++++++++++++++++++++++++++++*/

            lenInsSI = siIns;
            insHeapStr[lenInsSI] = '\0';

            insST = (*conNtAryST)[uiRef].insList;
            lastInsST = insST;

            while(insST)
            { /*Loop: check if already have insertion*/
               if(lenInsSI != insST->lenInsSI)
               { /*If: lengths are differnt (not match)*/
                  lastInsST = insST;
                  insST = insST->nextIns;
                  continue;
               } /*If: lengths are differnt (not match)*/
 
               /*Check if I have an match*/
               for(
                  siIns = 0;
                  insHeapStr[siIns]
                    == insST->insStr[siIns];
                  ++siIns
               ) if(insHeapStr[siIns] =='\0') break;

               if(
                  insHeapStr[siIns] ==insST->insStr[siIns]
               ){ /*If: I found a match*/
                  free(insHeapStr);
                  insHeapStr = 0;
                  ++insST->numHitsSI;
                  break;
               } /*If: I found a match*/

               lastInsST = insST;
               insST = insST->nextIns;
            } /*Loop: check if already have insertion*/
 
            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun12 Sec02 Sub04 Cat05:
            +   - insertion is new; add to insertion list
            \+++++++++++++++++++++++++++++++++++++++++++*/

            if(! insST)
            { /*If: I did not find a match*/

               if(lastInsST)
               { /*If: there are previous insertions*/
                  lastInsST->nextIns = mkIns_tbCon();

                  if(! lastInsST->nextIns)
                     goto memErr_fun12_sec04_sub02;

                  lastInsST = lastInsST->nextIns;
               } /*If: there are previous insertions*/

               else
               { /*Else: is first insertion*/
                  (*conNtAryST)[uiRef].insList =
                     mkIns_tbCon();

                  if(! (*conNtAryST)[uiRef].insList)
                     goto memErr_fun12_sec04_sub02;

                  lastInsST= (*conNtAryST)[uiRef].insList;
               } /*Else: is first insertion*/
 
               lastInsST->insStr = insHeapStr;
               lastInsST->lenInsSI = lenInsSI;
               lastInsST->numHitsSI = 1;

               insHeapStr = 0;
               lenInsSI = 0;
            } /*If: I did not find a match*/

            break;
         /*Case: Insertions*/

         /***********************************************\
         * Fun12 Sec03 Sub05:
         *   - Handle softmasking
         \***********************************************/

         case 'S':
         /*Case: soft masking*/
            endCigUI = uiBase + samSTPtr->cigArySI[uiCig];

            while(uiBase < endCigUI)
               ++uiBase;

            break;
         /*Case: soft masking*/

         /*Other cases invovle hard masking*/
      } /*Switch: Check the cigar entry type*/
   } /*Loop: Add bases to the consensus*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec04:
   ^   - clean up and return
   ^   o fun12 sec04 sub01:
   ^     - no error clean up
   ^   o fun12 sec04 sub02:
   ^     - memory error clean up
   ^   o fun12 sec04 sub03:
   ^     - header error clean up
   ^   o fun12 sec04 sub04:
   ^     - no map error clean up
   ^   o fun12 sec04 sub05:
   ^     - no sequence error clean up
   ^   o fun12 sec04 sub06:
   ^     - low mapq error clean up
   ^   o fun12 sec04 sub07:
   ^     - general clean up (both no error or error)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * fun12 sec04 sub01:
   *   - no error clean up
   \*****************************************************/

   errSC = 0;
   goto cleanUp_fun12_sec04_sub07;

   /*****************************************************\
   * fun12 sec04 sub02:
   *   - memory error clean up
   \*****************************************************/

   memErr_fun12_sec04_sub02:;
   errSC = def_memErr_tbConDefs;
   goto cleanUp_fun12_sec04_sub07;

   /*****************************************************\
   * fun12 sec04 sub03:
   *   - header error clean up
   \*****************************************************/

   headerErr_fun12_sec04_sub03:;
   errSC = def_header_tbConDefs;
   goto cleanUp_fun12_sec04_sub07;

   /*****************************************************\
   * fun12 sec04 sub04:
   *   - no map error clean up
   \*****************************************************/

   noMapErr_fun12_sec04_sub04:;
   errSC = def_noMap_tbConDefs;
   goto cleanUp_fun12_sec04_sub07;

   /*****************************************************\
   * fun12 sec04 sub05:
   *   - no sequence error clean up
   \*****************************************************/

   noSeqErr_fun12_sec04_sub05:;
   errSC = def_noSeq_tbConDefs;
   goto cleanUp_fun12_sec04_sub07;

   /*****************************************************\
   * fun12 sec04 sub06:
   *   - low mapq error clean up
   \*****************************************************/

   lowMapErr_fun12_sec04_sub06:;
   errSC = def_lowMapq_tbConDefs;
   goto cleanUp_fun12_sec04_sub07;

   /*****************************************************\
   * fun12 sec04 sub07:
   *   - general clean up (both no error or error)
   \*****************************************************/

   cleanUp_fun12_sec04_sub07:;

   if(insHeapStr)
      free(insHeapStr);

   insHeapStr = 0;

   return errSC;
}/*addRead_tbCon*/

/*-------------------------------------------------------\
| Fun13: collapse_tbCon
|   - collapses a conNt_tbCon strrucvt array into an array
|     of samEntry structs
| Input:
|   - conNtAryST:
|     o pointer to a conNt_tbCon struct array to collapse
|   - lenConAryUI:
|     o length of the conNt_tbCon array
|   - lenSamArySI:
|     o pointer to integer to have number of samEntry
|       structures made
|   - refIdStr:
|     o c-string with referernce sequence name
|   - settings:
|     o pointer to set_tbCon struct with consensus
|       building settings
|   - errSC:
|     o pointer to character to hold error type
| Output:
|   - Modifies:
|     o lenSamST to hold the returned samEntry arrray
|       length
|     o errSC:
|       - 0 for no errors
|       - def_noSeq_tbConDefs if could not build consensus
|       - def_memErr_tbConDefs if had memory error
|   - Returns:
|     o array of sam entry structures with consensus 
|       fragments.
|       - idealy this would be one structure, but for
|         amplicons it will be an array
|     o 0 for memory erors
\-------------------------------------------------------*/
struct samEntry *
collapse_tbCon(
   struct conNt_tbCon conNtAryST[], /*to collapse*/
   unsigned int lenConAryUI,   /*length of consensus*/
   signed int *lenSamArySI,    /*set to out array length*/
   signed char *refIdStr,      /*name of reference seq*/
   struct set_tbCon *settings, /*settings for collapsing*/
   signed char *errSC          /*error reports*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun13 TOC: collapse_tbCon
   '   - collapses a conNt_tbCon strrucvt array into an
   '     array of samEntry structs
   '   o fun13 sec01:
   '     - variable declerations
   '   o fun13 sec02:
   '     - count the number of fragments I have
   '   o fun13 sec03:
   '     - allocate memory for each fragment
   '   o fun13 sec04:
   '     - collapse each fragment into a samEntry struct
   '   o fun13 sec05:
   '     - collapse the consensus fragment
   '   o fun13 sec06:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint siFrag = 0;
   uint uiRef = 0;     /*ref base on*/
   uint uiEndRef = 0;  /*last base in amplicon/fragment*/
   uint lenFragUI = 0; /*length of fragment/amplicon*/
   uint uiBase = 0;    /*position on fragment (base on)*/

   uint uiCig = 0;     /*cigar entry on*/
   uint lenCigUI = 0;
   schar lastCigSC = '0';
   uint extraInsUI = 0; /*max extra bases from ins's*/

   /*General support*/
   float snpPerSupF = 0;
   float insPerSupF = 0;
   float delPerSupF = 0;

   /*Support for a particler snp/match*/
   uint nonMaskBaseUI = 0; /*number of non-anonymous*/
   float aPercSupF = 0;
   float tPercSupF = 0;
   float cPercSupF = 0;
   float gPercSupF = 0;

   /*Find the most supported insertion*/
   int keptReadsSI = 0;
   struct ins_tbCon *curInsST = 0;
   struct ins_tbCon *bestInsST = 0;

   struct samEntry *retSamST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec02:
   ^   - count number of fragments I have
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *lenSamArySI = 0;
   *errSC = 0;
   uiRef = 0;

   while(uiRef < lenConAryUI)
   { /*Loop: Find the number of fragments*/

      /*Find the starting position for this fragment*/
      while(
           conNtAryST[uiRef].ntKeptSI
         < settings->minDepthSI
      ){ /*Loop: Find the first base with enough depth*/
         ++uiRef;

         if(uiRef >= lenConAryUI)
            goto collapse_fun13_sec03;
      } /*Loop: Find the first base with enough depth*/

      /*Find the length of this fragment*/

      for(
         uiEndRef = uiRef;
         uiEndRef < lenConAryUI;
         ++uiEndRef
      ){ /*Loop: Find the memory for the fragment*/
         if(
               conNtAryST[uiEndRef].ntKeptSI
            < settings->minDepthSI
         ) break;
      } /*Loop: Find the memory for the fragment*/

      /*Find the length*/
      lenFragUI = uiRef - uiEndRef;

      if(lenFragUI < (uint) settings->minLenSI)
      { /*If: This is beneath the min length to keep*/
         uiRef = uiEndRef;
         continue;
      } /*If: This is beneath the min length to keep*/

      ++(*lenSamArySI);

      uiRef = uiEndRef;
   } /*Loop: Find the number of fragments*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec03:
   ^   - allocate memory for each fragment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   collapse_fun13_sec03:;

   if(*lenSamArySI == 0)
      goto noConErr_fun13_sec06_sub03;
      /*could not build the consensus*/

   retSamST =
      malloc(*lenSamArySI * sizeof(struct samEntry));

   if(!retSamST)
   { /*If: failed to allocate memory*/
      *lenSamArySI = 0;
      goto memErr_fun13_sec06_sub02;
   } /*If: failed to allocate memory*/

   /*Initialize the structures*/
   for(
      siFrag = 0;
      siFrag < *lenSamArySI;
      ++siFrag
   ) init_samEntry(&retSamST[siFrag]);

   /*set up memory*/
   for(
      siFrag = 0;
      siFrag < *lenSamArySI;
      ++siFrag
   ){ /*Loop: set up (allocate memory) to samEntry*/
      if(setup_samEntry(&retSamST[siFrag]))
         goto memErr_fun13_sec06_sub02;
   } /*Loop: set up (allocate memory) to samEntry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec04:
   ^   - collapse each fragment into a samEntry struct
   ^   o fun13 sec01 sub01:
   ^     - get length & check if fragment meets min length
   ^   o fun13 sec01 sub02:
   ^     - find the length and memory I need to allocate
   ^   o fun13 sec01 sub03:
   ^     - set up the sam entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun13 Sec01 Sub01:
   *   - get length & check if fragment meets min length
   \*****************************************************/

   uiRef = 0;
   siFrag = 0;

   while(uiRef < lenConAryUI)
   { /*Loop: Find the number of fragments*/

      /*Find the starting position for this fragment*/
      while(
           conNtAryST[uiRef].ntKeptSI
         < settings->minDepthSI
      ){ /*Loop: Find the first base with enough depth*/
         ++uiRef;

         if(uiRef >= lenConAryUI)
            goto noErr_fun13_sec06_sub01;
            /*at end of consensusq*/
      } /*Loop: Find the first base with enough depth*/

      /**************************************************\
      * Fun13 Sec01 Sub02:
      *   - find the length and memory I need to allocate
      *   o fun13 sec04 sub02 cat01:
      *     - find the length/start memory count loop
      *   o fun13 sec04 sub02 cat02:
      *     - find the percent support for each option
      *   o fun13 sec04 sub02 cat03:
      *     - check if have ins (if so get extra bases)
      *   o fun13 sec04 sub02 cat04:
      *     - check if an snp/match, del, or mask will
      *       be choosen
      *   o fun13 sec04 sub02 cat05:
      *     - find fragement length/check if long enough
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub02 Cat01:
      +   - find the length/start memory count loop
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      extraInsUI = 0;

      for(
         uiEndRef = uiRef;
         uiEndRef < lenConAryUI;
         ++uiEndRef
      ){ /*Loop: Find the memory for the fragment*/
         if(
               conNtAryST[uiEndRef].ntKeptSI
            < settings->minDepthSI
         ) break;

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun13 Sec04 Sub02 Cat02:
         +   - find the percent support for each option
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         /*The insertion is inbetween bases, which makes
         ` it hard to gauge. So, I am taking the worst
         ` case
         */
         keptReadsSI =
           max_genMath(
             conNtAryST[uiEndRef].ntKeptSI,
             conNtAryST[
                uiEndRef - (uiEndRef > 0)
             ].ntKeptSI
           );
            /*uiEndRef - (uiEndRef > 0) turns into
            `  uiEndRef - 1, when uiEndRef is > 0; else is
            `  uiEndRef
            */

         if(keptReadsSI > 0)
            insPerSupF =
                 (float) conNtAryST[uiEndRef].numInsSI
               / (float) keptReadsSI;
         else
            insPerSupF = 0;

         if(conNtAryST[uiEndRef].ntKeptSI > 0)
            delPerSupF =
                 (float) conNtAryST[uiEndRef].numDelSI
               / (float) conNtAryST[uiEndRef].ntKeptSI;
         else
            delPerSupF = 0;

         nonMaskBaseUI = conNtAryST[uiRef].numASI;
         nonMaskBaseUI += conNtAryST[uiRef].numTSI;
         nonMaskBaseUI += conNtAryST[uiRef].numGSI;
         nonMaskBaseUI += conNtAryST[uiRef].numCSI;

         snpPerSupF = nonMaskBaseUI;
         nonMaskBaseUI += conNtAryST[uiRef].numDelSI;


         if(nonMaskBaseUI > 0)
            snpPerSupF /= (float) nonMaskBaseUI;
         else
            snpPerSupF = 0;

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun13 Sec04 Sub02 Cat03:
         +   - check if have ins (if so get extra bases)
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(insPerSupF >= settings->minPercInsF)
         { /*If: I have an insertion*/
            lenCigUI += (lastCigSC != 'I');
            lastCigSC = 'I';

            bestInsST = conNtAryST[uiEndRef].insList;
            curInsST = bestInsST;

            while(curInsST)
            { /*Loop: Find the most supported insertion*/
              if(
                   bestInsST->numHitsSI
                 < curInsST->numHitsSI
              ) bestInsST = curInsST;

              curInsST = curInsST->nextIns;
            } /*Loop: Find the most supported insertion*/

            /*How many bases the insertions could add*/
            extraInsUI += bestInsST->lenInsSI;
         } /*If: I have an insertion*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun13 Sec04 Sub02 Cat04:
         +   - check if an snp/match, del, or mask will
         +     be choosen
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(snpPerSupF >= settings->minPercSnpF)
         { /*If: I am keeping an mathc/snp*/
            lenCigUI += (lastCigSC != 'M');
            lastCigSC = 'M';
         } /*If: I am keeping an mathc/snp*/

         else if(delPerSupF >= settings->minPercDelF)
         { /*Else If: I am keeping an deletion*/
            lenCigUI += (lastCigSC != 'D');
            lastCigSC = 'D';
         } /*Else If: I am keeping an deletion*/

         /*For the cigar a mask is same as an snp/match*/
         else
         { /*Else: I am masking (snp/match case)*/
            lenCigUI += (lastCigSC != 'M');
            lastCigSC = 'M';
         } /*Else: I am masking (snp/match case)*/
      } /*Loop: Find the memory for the fragment*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub02 Cat05:
      +   - find the fragement length/check if long enough
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*Find the length*/
      lenFragUI = uiEndRef - uiRef;
      ++lenFragUI; /*convert to index 1*/

      if(lenFragUI < (uint) settings->minLenSI)
      { /*If: This is beneath the min length to keep*/
         uiRef = uiEndRef;
         continue;
      } /*If: This is beneath the min length to keep*/

      /**************************************************\
      * Fun13 Sec04 Sub03:
      *   - set up the sam entry
      *   o fun13 sec04 sub03 cat01:
      *     - make the query id
      *   o fun13 sec04 sub03 cat02:
      *     - make the reference id
      *   o fun13 sec04 sub03 cat03:
      *     - set up the RNEXT entry
      *   o fun13 sec04 sub03 cat04:
      *     - set up reference positions/alinged length
      *   o fun13 sec04 sub03 cat05:
      *     - allocate memory for the sequence
      *   o fun13 sec04 sub03 cat06:
      *     - set up the Q-score entry
      *   o fun13 sec04 sub03 cat07:
      *     - allocate memory for the cigar types buffer
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub03 Cat01:
      +   - make the query id
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      retSamST[siFrag].lenQryIdUC =
         cpStr_ulCp(
             retSamST[siFrag].qryIdStr,
             refIdStr
         );

      retSamST[siFrag].qryIdStr[
         retSamST[siFrag].lenQryIdUC
      ] = '_';

      ++retSamST[siFrag].lenQryIdUC;

      retSamST[siFrag].lenQryIdUC +=
         numToStr(
            &retSamST[siFrag].qryIdStr[
               retSamST[siFrag].lenQryIdUC
            ],
         uiRef
      );

      retSamST[siFrag].qryIdStr[
         retSamST[siFrag].lenQryIdUC
      ] = '-';

      ++retSamST[siFrag].lenQryIdUC;

      retSamST[siFrag].lenQryIdUC +=
         (uchar)
         numToStr(
            &retSamST[siFrag].qryIdStr[
               retSamST[siFrag].lenQryIdUC
            ],
            uiEndRef
         ); /*Copy the ending position*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub03 Cat02:
      +   - make the reference id
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*Copy the reference id*/
      retSamST[siFrag].lenRefIdUC =
         cpDelim_ulCp(
            retSamST[siFrag].refIdStr,
            refIdStr,
            0,    /*'\0' = 0*/
            '\0'
         ); /*Copy the reference name*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub03 Cat03:
      +   - set up the RNEXT entry
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      retSamST[siFrag].rNextStr[0] = '*';
      retSamST[siFrag].rNextStr[1] = '\0';
      retSamST[siFrag].lenRNextUC = 1;

      /*The flag is already set up*/
      /*Not sure what to put for mapping quality, so using
      ` 0 (default)
      */

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub03 Cat04:
      +   - set up the reference positions/alinged length
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*Get the reference lengths*/
      retSamST[siFrag].refStartUI = uiRef;
      retSamST[siFrag].alnReadLenUI = lenFragUI;
      retSamST[siFrag].refEndUI = uiEndRef;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub03 Cat05:
      +   - allocate memory for the sequence
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      lenFragUI += extraInsUI;

      /*Make the sequence buffer*/
      if(retSamST[siFrag].lenSeqBuffUI < lenFragUI + 1)
      { /*If: I need to expand memory*/
         free(retSamST[siFrag].seqStr);

         retSamST[siFrag].seqStr =
            malloc((lenFragUI + 9) * (sizeof(schar)));

         if(! retSamST[siFrag].seqStr)
            goto memErr_fun13_sec06_sub02;

         retSamST[siFrag].lenSeqBuffUI = lenFragUI;
      } /*If: I need to expand memory*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub03 Cat06:
      +   - set up the Q-score entry
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*initSamEntry allocates over 2 bytes, so I am
      ` safe here
      */
      retSamST[siFrag].qStr[0] = '*';
      retSamST[siFrag].qStr[1] = '\0';
      retSamST[siFrag].lenQBuffUI = 2;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun13 Sec04 Sub03 Cat07:
      +   - allocate memory for the cigar types buffer
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(retSamST[siFrag].lenCigBuffUI < lenCigUI + 1)
      { /*If: I need to allocate memory for cigar buff*/
          /*Make the cigar types buffer*/
          retSamST[siFrag].cigTypeStr =
             malloc((lenCigUI + 9) * sizeof(schar));

          if(!retSamST[siFrag].cigTypeStr)
             goto memErr_fun13_sec06_sub02;

          retSamST[siFrag].lenCigBuffUI = lenCigUI;

          retSamST[siFrag].cigArySI =
             malloc((lenCigUI + 9) * sizeof(sint));

          if(!retSamST[siFrag].cigArySI)
             goto memErr_fun13_sec06_sub02;
      } /*If: I need to allocate memory for cigar buff*/
      
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun13 Sec05:
      ^   - Collapse the consensus fragment
      ^   o fun13 sec05 sub01:
      ^     - get the amount of support for each position
      ^   o fun13 sec05 sub02:
      ^     - check if I have an insertion. 
      ^     - insertion is always one base behind
      ^   o fun13 sec05 sub03:
      ^     - check if I have support for an snp/match
      ^   o fun13 sec05 sub04:
      ^     - check if this was a deletion instead of snp
      ^   o fun13 sec05 sub05:
      ^     - handle not enough support cases; likely
      ^       mixed infections
      ^   o fun13 sec05 sub06:
      ^     - check if I need to resize sequence buffer
      ^   o fun13 sec05 sub07:
      ^     - check if I need to resize the cigar buffer
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      /**************************************************\
      * Fun13 Sec05 Sub01:
      *   - get the amount of support for each position
      \**************************************************/

      uiBase = 0;
      uiCig = 0;

      retSamST[siFrag].cigTypeStr[uiCig] = '\0';

      while(uiRef < uiEndRef)
      { /*Loop: Collapse the fragment*/

         /*The insertion is inbetween bases, which makes
         ` it hard to gauge. So, I am taking the worst
         ` case. Here I want all positoins
         */
         keptReadsSI =
           max_genMath(
             conNtAryST[uiRef].ntKeptSI,
             conNtAryST[uiRef -(uiRef > 0)].ntKeptSI
           );
            /*uiRef - (uiRef > 0) turns into
            `  uiRef - 1, when uiRef is > 0; else is
            `  uiRef
            */

         if(keptReadsSI > 0)
            insPerSupF =
                 (float) conNtAryST[uiRef].numInsSI
               / (float) keptReadsSI;
         else
            insPerSupF = 0;

         /*For deletions an masked base is equivlent to no
         `   support
         */
         if(conNtAryST[uiRef].numDelSI > 0)
            delPerSupF =
                 (float) conNtAryST[uiRef].numDelSI
               / (float) conNtAryST[uiRef].ntKeptSI;
         else
            delPerSupF = 0;

         /*Find the number of non-anonymous bases
         `   For tbCon all anonymous bases are N's (masked)
         */
         nonMaskBaseUI = conNtAryST[uiRef].numASI;
         nonMaskBaseUI += conNtAryST[uiRef].numTSI;
         nonMaskBaseUI += conNtAryST[uiRef].numGSI;
         nonMaskBaseUI += conNtAryST[uiRef].numCSI;

         snpPerSupF = nonMaskBaseUI;
         nonMaskBaseUI += conNtAryST[uiRef].numDelSI;


         if(nonMaskBaseUI > 0)
            snpPerSupF /= (float) nonMaskBaseUI;
         else
            snpPerSupF = 0;

         /***********************************************\
         * Fun13 Sec05 Sub02:
         *   - check if I have an insertion.
         *   - insertion is always one base behind
         *   o fun13 sec05 sub02 cat01:
         *     - check if need new cigar entry + find
         *       insertion with the most support
         *   o fun13 sec05 sub02 cat01:
         *     - check if keeping inse / copy insertion
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun13 Sec05 Sub02 Cat01:
         +   - check if need new cigar entry + find
         +     insertion with the most support
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(insPerSupF >= settings->minPercInsF)
         { /*If: An insertions is supported here*/
            if(retSamST[siFrag].cigTypeStr[uiCig] != 'I')
            { /*If: This is a new cigar entry*/
               ++uiCig;
               retSamST[siFrag].cigTypeStr[uiCig] = 'I';
               retSamST[siFrag].cigArySI[uiCig] = 0;
            } /*If: This is a new cigar entry*/

            bestInsST = conNtAryST[uiRef].insList;
            curInsST = bestInsST;

            while(curInsST)
            { /*Loop: Find the most supported insertion*/
              if(
                   bestInsST->numHitsSI
                 < curInsST->numHitsSI
              ) bestInsST = curInsST;

              curInsST = curInsST->nextIns;
            } /*Loop: Find the most supported insertion*/

            keptReadsSI = curInsST->numHitsSI;

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun13 Sec05 Sub02 Cat02:
            +   - check if keeping ins / copy insertion
            \+++++++++++++++++++++++++++++++++++++++++++*/

            insPerSupF =
                 (float) bestInsST->numHitsSI
               / (float) keptReadsSI;

            if(insPerSupF >= settings->minPercInsF)
            { /*If: keeping insertion*/
               cpLen_ulCp(
                  &retSamST[siFrag].seqStr[uiBase],
                  bestInsST->insStr,
                  bestInsST->lenInsSI
               );

               uiBase += bestInsST->lenInsSI;

               retSamST[siFrag].cigArySI[uiCig] +=
                  bestInsST->lenInsSI;

               retSamST[siFrag].numInsUI +=
                  bestInsST->lenInsSI;
            } /*If: keeping insertion*/
         } /*If: an Insertion is supported here*/

         /***********************************************\
         * Fun13 Sec05 Sub03:
         *   - check if I have support for an snp/match
         *   o fun13 sec05 sub03 cat01:
         *     - check if need new cigar entry + find
         *       support for each base type
         *   o fun13 sec05 sub03 cat01:
         *     - find the base with the most support and
         *       see if it had enough support to keep
         \***********************************************/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun13 Sec05 Sub03 Cat01:
         +   - check if need new cigar entry + find
         +     support for each base type
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(nonMaskBaseUI < (uint) settings->minDepthSI)
            goto maskPos_fun13_sec05_sub05;

         else if(snpPerSupF >= settings->minPercSnpF)
         { /*If: snps were the best choice*/

            if(retSamST[siFrag].cigTypeStr[uiCig] != 'M')
            { /*If: This is a new cigar entry*/
               uiCig +=
                  (
                       retSamST[siFrag].cigTypeStr[uiCig]
                    != '\0'
                  );

               retSamST[siFrag].cigTypeStr[uiCig] = 'M';
               retSamST[siFrag].cigArySI[uiCig] = 0;
            } /*If: This is a new cigar entry*/

             aPercSupF = (float) conNtAryST[uiRef].numASI;
             aPercSupF /= (float) nonMaskBaseUI;

             tPercSupF = (float) conNtAryST[uiRef].numTSI;
             tPercSupF /= (float) nonMaskBaseUI;

             gPercSupF = (float) conNtAryST[uiRef].numGSI;
             gPercSupF /= (float) nonMaskBaseUI;

             cPercSupF = (float) conNtAryST[uiRef].numCSI;
             cPercSupF /= (float) nonMaskBaseUI;

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun13 Sec05 Sub03 Cat02:
            +   - find the base with the most support and
            +     see if it had enough support to keep
            \+++++++++++++++++++++++++++++++++++++++++++*/

            if(   aPercSupF > tPercSupF
               && aPercSupF > cPercSupF
               && aPercSupF > gPercSupF
            ) { /*If: A has the most support*/

               if(aPercSupF >= settings->minPercSnpF)
               { /*If: this was an A*/
                  retSamST[siFrag].seqStr[uiBase] = 'A';
                  ++retSamST[siFrag].numMatchUI;
               } /*If: this was an A*/

               else
               { /*Else: not enough support to call*/
                  retSamST[siFrag].seqStr[uiBase] =
                     settings->maskSC;

                  ++retSamST[siFrag].numMaskUI;
               } /*Else: not enough support to call*/
            } /*If: A has the most support*/

            else if(
                  tPercSupF > cPercSupF
               && tPercSupF > gPercSupF
            ) { /*Else If: T has the most support*/

               if(tPercSupF >= settings->minPercSnpF)
               { /*If: this was an T*/
                  retSamST[siFrag].seqStr[uiBase] = 'T';
                  ++retSamST[siFrag].numMatchUI;
               } /*If: this was an T*/

               else
               { /*Else: not enough support to call*/
                  retSamST[siFrag].seqStr[uiBase] =
                     settings->maskSC;

                  ++retSamST[siFrag].numMaskUI;
               } /*Else: not enough support to call*/
            } /*Else If: T has the most support*/

            else if(cPercSupF > gPercSupF)
            { /*Else If: C has the most support*/

               if(cPercSupF >= settings->minPercSnpF)
               { /*If: this was an C*/
                  retSamST[siFrag].seqStr[uiBase] = 'C';
                  ++retSamST[siFrag].numMatchUI;
               } /*If: this was an C*/

               else
               { /*Else: not enough support to call*/
                  retSamST[siFrag].seqStr[uiBase] =
                     settings->maskSC;

                  ++retSamST[siFrag].numMaskUI;
               } /*Else: not enough support to call*/
            } /*Else If: C has the most support*/

            else
            { /*Else If: G has the most support*/

               if(gPercSupF >= settings->minPercSnpF)
               { /*If: this was an G*/
                  retSamST[siFrag].seqStr[uiBase] = 'G';
                  ++retSamST[siFrag].numMatchUI;
               } /*If: this was an G*/

               else
               { /*Else: not enough support to call*/
                  retSamST[siFrag].seqStr[uiBase] =
                     settings->maskSC;

                  ++retSamST[siFrag].numMaskUI;
               } /*Else: not enough support to call*/
            } /*Else If: G has the most support*/

            ++uiBase;
            ++retSamST[siFrag].cigArySI[uiCig];
         } /*If: snps were the best choice*/

         /***********************************************\
         * Fun13 Sec05 Sub04:
         *   - check if this was a deletion instead of snp
         \***********************************************/

         /*Do nothing if deletion was selected*/
         else if(delPerSupF >= settings->minPercDelF)
         { /*Else if: there was a deletion*/
            if(retSamST[siFrag].cigTypeStr[uiCig] != 'D')
            { /*If: This is a new cigar entry*/
               uiCig +=
                  (
                       retSamST[siFrag].cigTypeStr[uiCig]
                    != '\0'
                  );

               retSamST[siFrag].cigTypeStr[uiCig] = 'D';
               retSamST[siFrag].cigArySI[uiCig] = 0;
            } /*If: This is a new cigar entry*/

            ++(retSamST[siFrag].cigArySI[uiCig]);
            ++retSamST[siFrag].numDelUI;
         } /*Else if: there was a deletion*/

         /***********************************************\
         * Fun13 Sec05 Sub05:
         *   - handle not enough support cases; likely
         *     mixed infections
         \***********************************************/

         else
         { /*Else: I have no support, assume mask snp*/
            maskPos_fun13_sec05_sub05:;

            if(retSamST[siFrag].cigTypeStr[uiCig] != 'M')
            { /*If: new cigar entry*/
               uiCig +=
                  (
                       retSamST[siFrag].cigTypeStr[uiCig]
                    != '\0'
                  );

               retSamST[siFrag].cigTypeStr[uiCig] = 'M';
               retSamST[siFrag].cigArySI[uiCig] = 0;
            } /*If: new cigar entry*/

            retSamST[siFrag].seqStr[uiBase] =
               settings->maskSC;

            ++uiBase;
            ++retSamST[siFrag].cigArySI[uiCig];
            ++retSamST[siFrag].numMaskUI;
         } /*Else: I have no support, assume mask snp*/

         ++uiRef;
      } /*Loop: Collapse the fragment*/

      retSamST[siFrag].readLenUI = uiBase;
      retSamST[siFrag].lenCigUI = uiCig + 1;
         /*The read length is at index 1, but the cigar is
         `   at index 0, so needs a + 1
         */

      ++siFrag;
   } /*Loop: Find the number of fragments*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec06:
   ^   - clean up and return
   ^   o fun13 sec07 sub01:
   ^     - no error clean up
   ^   o fun13 sec07 sub01:
   ^     - memory error clean up
   ^   o fun13 sec07 sub0x:
   ^     - clean up after an error (acutal clean up)
   ^   o fun13 sec07 sub0y:
   ^     - return the result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun13 Sec07 Sub01:
   *   - no error clean up
   \*****************************************************/

   noErr_fun13_sec06_sub01:;
   *errSC = 0;
   goto ret_fun13_sec06_sub06;

   /*****************************************************\
   * Fun13 Sec07 Sub02:
   *   - memory error clean up
   \*****************************************************/

   memErr_fun13_sec06_sub02:
   *errSC = def_memErr_tbConDefs;
   goto errCleanUp_fun13_sec06_sub0x;

   /*****************************************************\
   * Fun13 Sec07 Sub03:
   *   - no consensus error clean up
   \*****************************************************/

   noConErr_fun13_sec06_sub03:;
   *errSC = def_noSeq_tbConDefs;
   goto errCleanUp_fun13_sec06_sub0x;

   /*****************************************************\
   * Fun13 Sec07 Sub0x:
   *   - clean up after an error (actual clean up)
   \*****************************************************/

   errCleanUp_fun13_sec06_sub0x:;

   for(
      siFrag = 0;
      siFrag < *lenSamArySI;
      ++siFrag
   ) freeStack_samEntry(&retSamST[siFrag]);

   free(retSamST);
   retSamST = 0;
   *lenSamArySI = 0;

   goto ret_fun13_sec06_sub06;

   /*****************************************************\
   * Fun13 Sec07 Sub0y:
   *   - return the result
   \*****************************************************/

   ret_fun13_sec06_sub06:;

   return retSamST;
} /*collapse_tbCon*/

/*-------------------------------------------------------\
| Fun14: noFragCollapse_tbCon
|   - collapses a conNt_tbCon struct array into a single
|     samEntry struct (low read depth is masked)
| Input:
|   - conNtAryST:
|     o pointer to a conNt_tbCon struct array to collapse
|   - lenConAryUI:
|     o length of the conNt_tbCon array
|   - refIdStr:
|     o c-string with referernce sequence name
|   - numMaskUIPtr:
|     o pointer to unsigned int to have number of bases
|       masked in consensus
|   - settings:
|     o pointer to set_tbCon struct with consensus
|       building settings
|   - errSC:
|     o pointer to character to hold error type
| Output:
|   - Modifies:
|     o lenSamST to hold the returned samEntry arrray
|       length
|     o numMaskUIPtr to have number of masked bases
|     o errSC:
|       - 0 for no errors
|       - def_noSeq_tbConDefs if could not build consensus
|       - def_memErr_tbConDefs if had memory error
|   - Returns:
|     o array of sam entry structures with consensus 
|       fragments.
|       - idealy this would be one structure, but for
|         amplicons it will be an array
|     o 0 for memory erors
\-------------------------------------------------------*/
struct samEntry *
noFragCollapse_tbCon(
   struct conNt_tbCon conNtAryST[], /*to collapse*/
   unsigned int lenConAryUI,   /*length of consensus*/
   signed char *refIdStr,      /*name of reference seq*/
   unsigned int *numMaskUIPtr, /*# bases masked*/
   struct set_tbCon *settings, /*settings for collapsing*/
   signed char *errSC          /*error reports*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun14 TOC: collapse_tbCon
   '   - collapses a conNt_tbCon struct array into a
   '     single samEntry struct (low read depth is masked)
   '   o fun14 sec01:
   '     - variable declerations
   '   o fun14 sec02:
   '     - count the number of fragments I have
   '   o fun14 sec03:
   '     - allocate memory for each fragment
   '   o fun14 sec04:
   '     - collapse each fragment into a samEntry struct
   '   o fun14 sec05:
   '     - collapse the consensus fragment
   '   o fun14 sec06:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint uiRef = 0;     /*ref base on*/
   uint lenFragUI = 0; /*length of fragment/amplicon*/
   uint uiBase = 0;    /*position on fragment (base on)*/

   uint uiCig = 0;     /*cigar entry on*/
   uint lenCigUI = 0;
   schar lastCigSC = '0';
   uint extraInsUI = 0; /*max extra bases from ins's*/

   /*General support*/
   float snpPerSupF = 0;
   float insPerSupF = 0;
   float delPerSupF = 0;

   /*Support for a particler snp/match*/
   uint nonMaskBaseUI = 0; /*number of non-anonymous*/
   float aPercSupF = 0;
   float tPercSupF = 0;
   float cPercSupF = 0;
   float gPercSupF = 0;

   /*Find the most supported insertion*/
   int keptReadsSI = 0;
   struct ins_tbCon *curInsST = 0;
   struct ins_tbCon *bestInsST = 0;

   struct samEntry *retSamST = 0;

   sint startSI = -1;
   sint endSI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec02:
   ^   - allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *numMaskUIPtr = 0;

   retSamST = malloc(sizeof(struct samEntry));

   if(!retSamST)
      goto memErr_fun14_sec06_sub02;

   init_samEntry(retSamST);

   if( setup_samEntry(retSamST) )
      goto memErr_fun14_sec06_sub02;

   /*****************************************************\
   * Fun14 Sec04 Sub01:
   *   - find length and memory I need to allocate
   *   o fun14 sec04 sub01 cat01:
   *     - find start/end & check read depth + start loop
   *   o fun14 sec04 sub01 cat02:
   *     - find the percent support for each option
   *   o fun14 sec04 sub01 cat03:
   *     - check if have ins (if so get extra bases)
   *   o fun14 sec04 sub01 cat04:
   *     - check if an snp/match, del, or mask will
   *       be choosen
   *   o fun14 sec04 sub01 cat05:
   *     - find fragement length
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun14 Sec04 Sub01 Cat01:
   +   - find start/end and check read depth + start loop
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   extraInsUI = 0;

   for(
      uiRef = 0;
      uiRef < lenConAryUI;
      ++uiRef
   ){ /*Loop: find memory usage*/

      if(startSI < 0)
      { /*If: have not found the start yet*/
         if(
              conNtAryST[uiRef].ntKeptSI
            < settings->minDepthSI
         ) continue;

         else
            startSI = (sint) uiRef;
      } /*If: have not found the start yet*/

      if(
           conNtAryST[uiRef].ntKeptSI
         >= settings->minDepthSI
      ) endSI = (sint) uiRef; /*find last base*/

      else
      { /*Else: low read depth*/
         lenCigUI += (lastCigSC != 'M');
         lastCigSC = 'M';
      } /*Else: low read depth*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun14 Sec04 Sub01 Cat02:
      +   - find the percent support for each option
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*The insertion is inbetween bases, which makes
      ` it hard to gauge. So, I am taking the worst
      ` case
      */
      keptReadsSI =
        max_genMath(
          conNtAryST[uiRef].ntKeptSI,
          conNtAryST[
             uiRef - (uiRef > 0)
          ].ntKeptSI
        );
         /*uiRef - (uiRef > 0) turns into
         `  uiRef - 1, when uiRef is > 0; else is
         `  uiRef
         */

      if(keptReadsSI > 0)
         insPerSupF =
              (float) conNtAryST[uiRef].numInsSI
            / (float) keptReadsSI;
      else
         insPerSupF = 0;

      if(conNtAryST[uiRef].ntKeptSI > 0)
         delPerSupF =
              (float) conNtAryST[uiRef].numDelSI
            / (float) conNtAryST[uiRef].ntKeptSI;
      else
         delPerSupF = 0;

      nonMaskBaseUI = conNtAryST[uiRef].numASI;
      nonMaskBaseUI += conNtAryST[uiRef].numTSI;
      nonMaskBaseUI += conNtAryST[uiRef].numGSI;
      nonMaskBaseUI += conNtAryST[uiRef].numCSI;

      snpPerSupF = nonMaskBaseUI;
      nonMaskBaseUI += conNtAryST[uiRef].numDelSI;


      if(nonMaskBaseUI > 0)
         snpPerSupF /= (float) nonMaskBaseUI;
      else
         snpPerSupF = 0;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun14 Sec04 Sub01 Cat03:
      +   - check if have ins (if so get extra bases)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(insPerSupF >= settings->minPercInsF)
      { /*If: I have an insertion*/
         lenCigUI += (lastCigSC != 'I');
         lastCigSC = 'I';

         bestInsST = conNtAryST[uiRef].insList;
         curInsST = bestInsST;

         while(curInsST)
         { /*Loop: Find the most supported insertion*/
           if(
                bestInsST->numHitsSI
              < curInsST->numHitsSI
           ) bestInsST = curInsST;

           curInsST = curInsST->nextIns;
         } /*Loop: Find the most supported insertion*/

         /*How many bases the insertions could add*/
         extraInsUI += bestInsST->lenInsSI;
      } /*If: I have an insertion*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun14 Sec04 Sub01 Cat04:
      +   - check if an snp/match, del, or mask will
      +     be choosen
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(snpPerSupF >= settings->minPercSnpF)
      { /*If: I am keeping an mathc/snp*/
         lenCigUI += (lastCigSC != 'M');
         lastCigSC = 'M';
      } /*If: I am keeping an mathc/snp*/

      else if(delPerSupF >= settings->minPercDelF)
      { /*Else If: I am keeping an deletion*/
         lenCigUI += (lastCigSC != 'D');
         lastCigSC = 'D';
      } /*Else If: I am keeping an deletion*/

      /*For the cigar a mask is same as an snp/match*/
      else
      { /*Else: I am masking (snp/match case)*/
         lenCigUI += (lastCigSC != 'M');
         lastCigSC = 'M';
      } /*Else: I am masking (snp/match case)*/
   } /*Loop: find memory usage*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun14 Sec04 Sub01 Cat05:
   +   - find the fragement length
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(startSI < 0)
      goto noConErr_fun14_sec06_sub03;

   lenFragUI = endSI - startSI;
   /*will add in extra insertion nucleotides later*/
   ++lenFragUI; /*convert to index 1*/

   if(lenFragUI < (uint) settings->minLenSI)
      goto noConErr_fun14_sec06_sub03; /*to small*/

   /*****************************************************\
   * Fun14 Sec04 Sub02:
   *   - set up the sam entry
   *   o fun14 sec04 sub02 cat01:
   *     - make the query id
   *   o fun14 sec04 sub02 cat02:
   *     - make the reference id
   *   o fun14 sec04 sub02 cat03:
   *     - set up the RNEXT entry
   *   o fun14 sec04 sub02 cat04:
   *     - set up reference positions/alinged length
   *   o fun14 sec04 sub02 cat05:
   *     - allocate memory for the sequence
   *   o fun14 sec04 sub02 cat06:
   *     - set up the Q-score entry
   *   o fun14 sec04 sub02 cat07:
   *     - allocate memory for the cigar types buffer
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun14 Sec04 Sub02 Cat01:
   +   - make the query id
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   retSamST->lenQryIdUC =
      cpStr_ulCp(
          retSamST->qryIdStr,
          refIdStr
      ); /*get reference id*/

   retSamST->qryIdStr[retSamST->lenQryIdUC] = '_';
   ++retSamST->lenQryIdUC;

   retSamST->lenQryIdUC +=
      (uchar)
      numToStr(
         &retSamST->qryIdStr[retSamST->lenQryIdUC],
         (uint) startSI
      ); /*add starting coordinate*/

   retSamST->qryIdStr[retSamST->lenQryIdUC] = '-';
   ++retSamST->lenQryIdUC;

   retSamST->lenQryIdUC +=
      (uchar)
      numToStr(
         &retSamST->qryIdStr[retSamST->lenQryIdUC],
         (uint) endSI
      ); /*Copy the ending position*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun14 Sec04 Sub02 Cat02:
   +   - make the reference id
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /*Copy the reference id*/
   retSamST->lenRefIdUC =
      cpDelim_ulCp(
         retSamST->refIdStr,
         refIdStr,
         0,    /*'\0' = 0*/
         '\0'
      ); /*Copy the reference name*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun14 Sec04 Sub02 Cat03:
   +   - set up the RNEXT entry
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   retSamST->rNextStr[0] = '*';
   retSamST->rNextStr[1] = '\0';
   retSamST->lenRNextUC = 1;

   /*The flag is already set up*/
   /*Not sure what to put for mapping quality, so using
   ` 0 (default)
   */

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun14 Sec04 Sub02 Cat04:
   +   - set up the reference positions/alinged length
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /*Get the reference lengths*/
   retSamST->refStartUI = startSI;
   retSamST->alnReadLenUI = endSI - startSI;
   retSamST->refEndUI = endSI;

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun14 Sec04 Sub02 Cat05:
   +   - allocate memory for the sequence
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   lenFragUI += extraInsUI;

   /*Make the sequence buffer*/
   if(retSamST->lenSeqBuffUI < lenFragUI + 1)
   { /*If: I need to expand memory*/
      free(retSamST->seqStr);

      retSamST->seqStr =
         malloc((lenFragUI + 9) * (sizeof(schar)));

      if(! retSamST->seqStr)
         goto memErr_fun14_sec06_sub02;

      retSamST->lenSeqBuffUI = lenFragUI;
   } /*If: I need to expand memory*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun14 Sec04 Sub02 Cat06:
   +   - set up the q-score entry
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   /*initSamEntry allocates over 2 bytes, so I am
   ` safe here
   */
   retSamST->qStr[0] = '*';
   retSamST->qStr[1] = '\0';
   retSamST->lenQBuffUI = 2;

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun14 Sec04 Sub02 Cat07:
   +   - allocate memory for the cigar types buffer
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(retSamST->lenCigBuffUI < lenCigUI + 1)
   { /*If: I need to allocate memory for cigar buff*/
       /*Make the cigar types buffer*/
       retSamST->cigTypeStr =
          malloc((lenCigUI + 9) * sizeof(schar));

       if(!retSamST->cigTypeStr)
          goto memErr_fun14_sec06_sub02;

       retSamST->lenCigBuffUI = lenCigUI;

       retSamST->cigArySI =
          malloc((lenCigUI + 9) * sizeof(sint));

       if(!retSamST->cigArySI)
          goto memErr_fun14_sec06_sub02;
   } /*If: I need to allocate memory for cigar buff*/
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec05:
   ^   - collapse consensus
   ^   o fun14 sec05 sub01:
   ^     - check read depth + start loop
   ^   o fun14 sec05 sub02:
   ^     - get the amount of support for each position
   ^   o fun14 sec05 sub03:
   ^     - check if I have an insertion. 
   ^     - insertion is always one base behind
   ^   o fun14 sec05 sub04:
   ^     - check if I have support for an snp/match
   ^   o fun14 sec05 sub05:
   ^     - check if this was a deletion instead of snp
   ^   o fun14 sec05 sub06:
   ^     - handle not enough support cases; likely
   ^       mixed infections
   ^   o fun14 sec05 sub07:
   ^     - check if I need to resize sequence buffer
   ^   o fun14 sec05 sub08:
   ^     - check if I need to resize the cigar buffer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/


   /*****************************************************\
   * Fun14 Sec05 Sub01:
   *   - check read depth + start loop
   \*****************************************************/

   uiBase = 0;
   uiCig = 0;

   retSamST->cigTypeStr[uiCig] = '\0';

   while(startSI <= endSI)
   { /*Loop: collapse consensus*/

      if(
           conNtAryST[startSI].ntKeptSI
         < settings->minDepthSI
      ){ /*If: low read depth*/

         retSamST->seqStr[uiBase] = settings->maskSC;
         ++(*numMaskUIPtr);

         ++uiBase;
         ++startSI;

         if(retSamST->cigTypeStr[uiCig] != 'M')
         { /*If: need to add a new cigar entry*/
            uiCig +=
               (
                    retSamST->cigTypeStr[uiCig]
                 != '\0'
               );

            retSamST->cigTypeStr[uiCig] = 'M';
            retSamST->cigArySI[uiCig] = 0;
         } /*If: need to add a new cigar entry*/

         ++retSamST->cigArySI[uiCig];
         ++retSamST->numMaskUI;

         continue;
      } /*If: low read depth*/

      /**************************************************\
      * Fun14 Sec05 Sub02:
      *   - get the amount of support for each position
      \**************************************************/

      /*The insertion is inbetween bases, which makes
      ` it hard to gauge. So, I am taking the worst
      ` case. Here I want all positoins
      */
      keptReadsSI =
        max_genMath(
          conNtAryST[startSI].ntKeptSI,
          conNtAryST[startSI -(startSI > 0)].ntKeptSI
        );
         /*startSI - (startSI > 0) turns into
         `  startSI - 1, when startSI is > 0; else is
         `  startSI
         */

      if(keptReadsSI > 0)
         insPerSupF =
              (float) conNtAryST[startSI].numInsSI
            / (float) keptReadsSI;
      else
         insPerSupF = 0;

      /*For deletions an masked base is equivlent to no
      `   support
      */
      if(conNtAryST[startSI].ntKeptSI > 0)
         delPerSupF =
              (float) conNtAryST[startSI].numDelSI
            / (float) conNtAryST[startSI].ntKeptSI;
      else
         delPerSupF = 0;

      /*Find the number of non-anonymous bases
      `   For tbCon all anonymous bases are N's (masked)
      */
      nonMaskBaseUI = conNtAryST[startSI].numASI;
      nonMaskBaseUI += conNtAryST[startSI].numTSI;
      nonMaskBaseUI += conNtAryST[startSI].numGSI;
      nonMaskBaseUI += conNtAryST[startSI].numCSI;

      snpPerSupF = nonMaskBaseUI;
      nonMaskBaseUI += conNtAryST[startSI].numDelSI;

      if(nonMaskBaseUI > 0)
         snpPerSupF /= (float) nonMaskBaseUI;
      else
         snpPerSupF = 0;

      /**************************************************\
      * Fun14 Sec05 Sub03:
      *   - check if I have an insertion.
      *   - insertion is always one base behind
      *   o fun14 sec05 sub03 cat01:
      *     - check if need new cigar entry + find
      *       insertion with the most support
      *   o fun14 sec05 sub03 cat01:
      *     - check if keeping inse / copy insertion
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun14 Sec05 Sub03 Cat01:
      +   - check if need new cigar entry + find
      +     insertion with the most support
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(insPerSupF >= settings->minPercInsF)
      { /*If: An insertions is supported here*/
         if(retSamST->cigTypeStr[uiCig] != 'I')
         { /*If: This is a new cigar entry*/
            ++uiCig;
            retSamST->cigTypeStr[uiCig] = 'I';
            retSamST->cigArySI[uiCig] = 0;
         } /*If: This is a new cigar entry*/

         bestInsST = conNtAryST[startSI].insList;
         curInsST = bestInsST;

         while(curInsST)
         { /*Loop: Find the most supported insertion*/
           if(
                bestInsST->numHitsSI
              < curInsST->numHitsSI
           ) bestInsST = curInsST;

           curInsST = curInsST->nextIns;
         } /*Loop: Find the most supported insertion*/

         keptReadsSI = curInsST->numHitsSI;

         /*+++++++++++++++++++++++++++++++++++++++++++\
         + Fun14 Sec05 Sub03 Cat02:
         +   - check if keeping ins / copy insertion
         \+++++++++++++++++++++++++++++++++++++++++++*/

         insPerSupF =
              (float) bestInsST->numHitsSI
            / (float) keptReadsSI;

         if(insPerSupF >= settings->minPercInsF)
         { /*If: keeping insertion*/
            cpLen_ulCp(
               &retSamST->seqStr[uiBase],
               bestInsST->insStr,
               bestInsST->lenInsSI
            );

            uiBase += bestInsST->lenInsSI;

            retSamST->cigArySI[uiCig] +=
               bestInsST->lenInsSI;

            retSamST->numInsUI +=
               bestInsST->lenInsSI;
         } /*If: keeping insertion*/
      } /*If: an Insertion is supported here*/

      /**************************************************\
      * Fun14 Sec05 Sub04:
      *   - check if I have support for an snp/match
      *   o fun14 sec05 sub04 cat01:
      *     - check if need new cigar entry + find
      *       support for each base type
      *   o fun14 sec05 sub04 cat01:
      *     - find the base with the most support and
      *       see if it had enough support to keep
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun14 Sec05 Sub04 Cat01:
      +   - check if need new cigar entry + find
      +     support for each base type
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(nonMaskBaseUI < (uint) settings->minDepthSI)
         goto maskPos_fun14_sec05_sub05;

      else if(snpPerSupF >= settings->minPercSnpF)
      { /*If: snps were the best choice*/

         if(retSamST->cigTypeStr[uiCig] != 'M')
         { /*If: This is a new cigar entry*/
            uiCig +=
               (
                    retSamST->cigTypeStr[uiCig]
                 != '\0'
               );

            retSamST->cigTypeStr[uiCig] = 'M';
            retSamST->cigArySI[uiCig] = 0;
         } /*If: This is a new cigar entry*/

          aPercSupF = (float) conNtAryST[startSI].numASI;
          aPercSupF /= (float) nonMaskBaseUI;

          tPercSupF = (float) conNtAryST[startSI].numTSI;
          tPercSupF /= (float) nonMaskBaseUI;

          gPercSupF = (float) conNtAryST[startSI].numGSI;
          gPercSupF /= (float) nonMaskBaseUI;

          cPercSupF = (float) conNtAryST[startSI].numCSI;
          cPercSupF /= (float) nonMaskBaseUI;

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun14 Sec05 Sub04 Cat02:
         +   - find the base with the most support and
         +     see if it had enough support to keep
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(   aPercSupF > tPercSupF
            && aPercSupF > cPercSupF
            && aPercSupF > gPercSupF
         ) { /*If: A has the most support*/

            if(aPercSupF >= settings->minPercSnpF)
            { /*If: this was an A*/
               retSamST->seqStr[uiBase] = 'A';
               ++retSamST->numMatchUI;
            } /*If: this was an A*/

            else
            { /*Else: not enough support to call*/
               retSamST->seqStr[uiBase] =settings->maskSC;
               ++(*numMaskUIPtr);
               ++retSamST->numMaskUI;
            } /*Else: not enough support to call*/
         } /*If: A has the most support*/

         else if(
               tPercSupF > cPercSupF
            && tPercSupF > gPercSupF
         ) { /*Else If: T has the most support*/

            if(tPercSupF >= settings->minPercSnpF)
            { /*If: this was an T*/
               retSamST->seqStr[uiBase] = 'T';
               ++retSamST->numMatchUI;
            } /*If: this was an T*/

            else
            { /*Else: not enough support to call*/
               retSamST->seqStr[uiBase] =settings->maskSC;
               ++(*numMaskUIPtr);
               ++retSamST->numMaskUI;
            } /*Else: not enough support to call*/
         } /*Else If: T has the most support*/

         else if(cPercSupF > gPercSupF)
         { /*Else If: C has the most support*/

            if(cPercSupF >= settings->minPercSnpF)
            { /*If: this was an C*/
               retSamST->seqStr[uiBase] = 'C';
               ++retSamST->numMatchUI;
            } /*If: this was an C*/

            else
            { /*Else: not enough support to call*/
               retSamST->seqStr[uiBase] =settings->maskSC;
               ++(*numMaskUIPtr);
               ++retSamST->numMaskUI;
            } /*Else: not enough support to call*/
         } /*Else If: C has the most support*/

         else
         { /*Else If: G has the most support*/

            if(gPercSupF >= settings->minPercSnpF)
            { /*If: this was an G*/
               retSamST->seqStr[uiBase] = 'G';
               ++retSamST->numMatchUI;
            } /*If: this was an G*/

            else
            { /*Else: not enough support to call*/
               retSamST->seqStr[uiBase] =settings->maskSC;
               ++(*numMaskUIPtr);
               ++retSamST->numMaskUI;
            } /*Else: not enough support to call*/
         } /*Else If: G has the most support*/

         ++uiBase;
         ++retSamST->cigArySI[uiCig];
      } /*If: snps were the best choice*/

      /***********************************************\
      * Fun14 Sec05 Sub05:
      *   - check if this was a deletion instead of snp
      \***********************************************/

      /*Do nothing if deletion was selected*/
      else if(delPerSupF >= settings->minPercDelF)
      { /*Else if: there was a deletion*/
         if(retSamST->cigTypeStr[uiCig] != 'D')
         { /*If: This is a new cigar entry*/
            uiCig +=
               (
                    retSamST->cigTypeStr[uiCig]
                 != '\0'
               );

            retSamST->cigTypeStr[uiCig] = 'D';
            retSamST->cigArySI[uiCig] = 0;
         } /*If: This is a new cigar entry*/

         ++(retSamST->cigArySI[uiCig]);
         ++retSamST->numDelUI;
      } /*Else if: there was a deletion*/

      /***********************************************\
      * Fun14 Sec05 Sub06:
      *   - handle not enough support cases; likely
      *     mixed infections
      \***********************************************/

      else
      { /*Else: I have no support, assume mask snp*/
         maskPos_fun14_sec05_sub05:;

         if(retSamST->cigTypeStr[uiCig] != 'M')
         { /*If: new cigar entry*/
            uiCig +=
               (
                    retSamST->cigTypeStr[uiCig]
                 != '\0'
               );

            retSamST->cigTypeStr[uiCig] = 'M';
            retSamST->cigArySI[uiCig] = 0;
         } /*If: new cigar entry*/

         retSamST->seqStr[uiBase] = settings->maskSC;
         ++(*numMaskUIPtr);

         ++uiBase;
         ++retSamST->cigArySI[uiCig];
         ++retSamST->numMaskUI;
      } /*Else: I have no support, assume mask snp*/

      ++startSI;
   } /*Loop: Collapse the fragment*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec06:
   ^   - clean up and return
   ^   o fun14 sec07 sub01:
   ^     - no error clean up
   ^   o fun14 sec07 sub01:
   ^     - memory error clean up
   ^   o fun14 sec07 sub0x:
   ^     - clean up after an error (acutal clean up)
   ^   o fun14 sec07 sub0y:
   ^     - return the result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun14 Sec07 Sub01:
   *   - no error clean up
   \*****************************************************/

   retSamST->readLenUI = uiBase;
   retSamST->lenCigUI = uiCig + 1;
   *errSC = 0;
   goto ret_fun14_sec06_sub06;

   /*****************************************************\
   * Fun14 Sec07 Sub02:
   *   - memory error clean up
   \*****************************************************/

   memErr_fun14_sec06_sub02:
   *errSC = def_memErr_tbConDefs;
   goto errCleanUp_fun14_sec06_sub0x;

   /*****************************************************\
   * Fun14 Sec07 Sub03:
   *   - no consensus error clean up
   \*****************************************************/

   noConErr_fun14_sec06_sub03:;
   *errSC = def_noSeq_tbConDefs;
   goto errCleanUp_fun14_sec06_sub0x;

   /*****************************************************\
   * Fun14 Sec07 Sub0x:
   *   - clean up after an error (actual clean up)
   \*****************************************************/

   errCleanUp_fun14_sec06_sub0x:;
   freeHeap_samEntry(retSamST);
   retSamST = 0;
   goto ret_fun14_sec06_sub06;

   /*****************************************************\
   * Fun14 Sec07 Sub0y:
   *   - return the result
   \*****************************************************/

   ret_fun14_sec06_sub06:;
   return retSamST;
} /*noFragcollapse_tbCon*/

/*-------------------------------------------------------\
| Fun15: pvar_tbCon
|   - print entries in a conNt_tbCon structure array that
|     are above the minimum read depth
| Input:
|   - conNtAryST:
|     o pointer to a conNt_tbCon structure array to print
|   - lenConAryUI:
|     o length of the conNt_tbCon array
|   - minDepthSI:
|     o minimum read depth to print out a position
|   - refIdStr:
|     o c-string with referernce sequence name
|   - outFILE:
|     o c-string with name of file to print everything to
| Output:
|   - Prints:
|     o entries in conNtAryST to outFILE
|   - Returns:
|     o 0 for success
|     o def_fileErr_tbConDefs for file errors
\-------------------------------------------------------*/
char
pvar_tbCon(
   struct conNt_tbCon conNtAryST[], /*consensus array*/
   unsigned int lenConAryUI,     /*length of conNtAryST*/
   signed char *refIdStr,        /*name of reference seq*/
   struct set_tbCon *settings,   /*settings for printing*/
   signed char *outStr           /*file to print to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun15 TOC: pvar_tbCon
   '   o fun15 sec01:
   '     - Variable declerations
   '   o fun15 sec02:
   '     - Print out the header
   '   o fun15 sec03:
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec01:
   ^   - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0; /*for error reporting*/
   uint uiBase = 0;
   struct ins_tbCon *insST = 0;

   uint nonMaskBaseUI = 0; /*number of non-anonymous*/
   uint maskedBasesUI = 0; /*number of non-anonymous*/
   float percSupF = 0;

   sint leastSupInsSI = 0;
   sint insTotalSI = 0;

   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec02:
   ^   - Print out the header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
         ! outStr
      || *outStr == '*'
   ) outFILE = stdout;

   else
   { /*Else: I need to open a file*/
      outFILE =
         fopen(
            (char *) outStr,
            "w"
         );

      if(! outFILE)
         goto fileErr_fun15_sec04;
   } /*Else: I need to open a file*/

   /*Print the header*/
   fprintf(
      outFILE,
      "refId\tposition\ttype\tsequence\tsupport"
   );

   fprintf(
     outFILE,
     "\tpercSupport\tmasked\tkeptBases\ttotalBases\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec03:
   ^   - Print out each base
   ^   o fun15 sec03 sub01:
   ^     - Start loop and print insertion entries
   ^   o fun15 sec03 sub02:
   ^     - Print out the snp/match entry for A
   ^   o fun15 sec03 sub03:
   ^     - Print out the snp/match entry for T
   ^   o fun15 sec03 sub04:
   ^     - Print out the snp/match entry for C
   ^   o fun15 sec03 sub05:
   ^     - Print out the snp/match entry for G
   ^   o fun15 sec03 sub06:
   ^     - Print out the deletion entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun15 Sec03 Sub01:
   *   - Start loop and print insertion entries
   \*****************************************************/

   for(uiBase = 0; uiBase < lenConAryUI; ++uiBase)
   { /*Loop: Print out the consensus*/

      if(!uiBase)
         goto firstBaseNoIns;

      insST = conNtAryST[uiBase].insList;

      leastSupInsSI =
        max_genMath(
           conNtAryST[uiBase].ntKeptSI,
           conNtAryST[uiBase - 1].ntKeptSI
        );

     insTotalSI =
        max_genMath(
           conNtAryST[uiBase].totalNtSI,
           conNtAryST[uiBase - 1].totalNtSI
        );

      while(insST)
      { /*Loop: Check insertions*/
         if(insST->numHitsSI < settings->minPrintDepthSI)
         { /*If: There is not enough read depth*/
            insST = insST->nextIns;
            continue;
         } /*If: There is not enough read depth*/

         percSupF =
           (float) insST->numHitsSI/(float) leastSupInsSI;

         if(percSupF < settings->printMinSupInsF)
         { /*If: I am not printing this ins*/
            insST = insST->nextIns;
            continue;
         } /*If: I am not printing this ins*/

         fprintf(
            outFILE,
            "%s\t%i\tins\t%s\t%i\t%f\tNA\t%i\t%i\n",
            refIdStr,
            uiBase,
            insST->insStr,
            insST->numHitsSI,
            percSupF,
            leastSupInsSI,
            insTotalSI
         );

         insST = insST->nextIns;
      } /*Loop: Check insertions*/

      /**************************************************\
      * Fun15 Sec03 Sub02:
      *   - Print out the snp/match entry for A
      \**************************************************/

      firstBaseNoIns:;

      if(
           conNtAryST[uiBase].ntKeptSI
         < settings->minPrintDepthSI
      ) continue; /*Not enough support to print out*/

      nonMaskBaseUI = conNtAryST[uiBase].numASI;
      nonMaskBaseUI += conNtAryST[uiBase].numTSI;
      nonMaskBaseUI += conNtAryST[uiBase].numGSI;
      nonMaskBaseUI += conNtAryST[uiBase].numCSI;
      nonMaskBaseUI += conNtAryST[uiBase].numDelSI;

      maskedBasesUI = conNtAryST[uiBase].ntKeptSI;
      maskedBasesUI -= nonMaskBaseUI;

      percSupF = (float) conNtAryST[uiBase].numASI;
      percSupF /= (float) nonMaskBaseUI;

      if(
           conNtAryST[uiBase].numASI
         < settings->minPrintDepthSI
      ) ;

      else if(percSupF >= settings->printMinSupSnpF)
      { /*If: I had enough support to print out A*/
         fprintf(
            outFILE,
            "%s\t%i\tbase\tA\t%i\t%f\t%u\t%i\t%i\n",
            refIdStr,
            uiBase + 1,
            conNtAryST[uiBase].numASI,
            percSupF,
            maskedBasesUI,
            conNtAryST[uiBase].ntKeptSI,
            conNtAryST[uiBase].totalNtSI
         );
      } /*If: I had enough support to print out A*/

      /**************************************************\
      * Fun15 Sec03 Sub03:
      *   - Print out the snp/match entry for T
      \**************************************************/

      percSupF = (float) conNtAryST[uiBase].numTSI;
      percSupF /= (float) nonMaskBaseUI;

      if(
           conNtAryST[uiBase].numTSI
         < settings->minPrintDepthSI
      ) ;

      else if(percSupF >= settings->printMinSupSnpF)
      { /*If: I had enough support to print an T*/
         fprintf(
            outFILE,
            "%s\t%i\tbase\tT\t%i\t%f\t%u\t%i\t%i\n",
            refIdStr,
            uiBase + 1,
            conNtAryST[uiBase].numTSI,
            percSupF,
            maskedBasesUI,
            conNtAryST[uiBase].ntKeptSI,
            conNtAryST[uiBase].totalNtSI
         );
      } /*If: I had enough support to print an T*/

      /**************************************************\
      * Fun15 Sec03 Sub04:
      *   - Print out the snp/match entry for C
      \**************************************************/

      percSupF = (float) conNtAryST[uiBase].numCSI;
      percSupF /= (float) nonMaskBaseUI;

      if(
           conNtAryST[uiBase].numCSI
         < settings->minPrintDepthSI
      ) ;

      else if(percSupF >= settings->printMinSupSnpF)
      { /*If: I had enough support to print an C*/
         fprintf(
            outFILE,
            "%s\t%i\tbase\tC\t%i\t%f\t%u\t%i\t%u\n",
            refIdStr,
            uiBase + 1,
            conNtAryST[uiBase].numCSI,
            percSupF,
            maskedBasesUI,
            conNtAryST[uiBase].ntKeptSI,
            conNtAryST[uiBase].totalNtSI
         );
      } /*If: I had enough support to print an C*/

      /**************************************************\
      * Fun15 Sec03 Sub05:
      *   - Print out the snp/match entry for G
      \**************************************************/

      percSupF = (float) conNtAryST[uiBase].numGSI;
      percSupF /= (float) nonMaskBaseUI;

      if(
           conNtAryST[uiBase].numGSI
         < settings->minPrintDepthSI
      ) ;

      else if(percSupF >= settings->printMinSupSnpF)
      { /*If: I had enough support to print an G*/
         fprintf(
            outFILE,
            "%s\t%i\tbase\tG\t%i\t%f\t%u\t%i\t%u\n",
            refIdStr,
            uiBase + 1,
            conNtAryST[uiBase].numGSI,
            percSupF,
            maskedBasesUI,
            conNtAryST[uiBase].ntKeptSI,
            conNtAryST[uiBase].totalNtSI
         );
      } /*If: I had enough support to print an G*/

      /**************************************************\
      * Fun15 Sec03 Sub06:
      *   - Print out the deletion entry
      \**************************************************/

      if(
           conNtAryST[uiBase].numDelSI
         > settings->minPrintDepthSI
      ){ /*If: The a base has enough depth*/
         percSupF = (float) conNtAryST[uiBase].numDelSI;

         percSupF /=
            (float) conNtAryST[uiBase].ntKeptSI;

         if(percSupF >= settings->printMinSupDelF)
            fprintf(
               outFILE,
               "%s\t%i\tdel\tdel\t%i\t%f\t%u\t%i\t%u\n",
               refIdStr,
               uiBase + 1,
               conNtAryST[uiBase].numDelSI,
               percSupF,
               maskedBasesUI,
               conNtAryST[uiBase].ntKeptSI,
               conNtAryST[uiBase].totalNtSI
            );
      } /*If: The a base has enough depth*/
   } /*Loop: Print out the consensus*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec04:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto cleanUp_fun15_sec04;

   fileErr_fun15_sec04:;
   errSC = def_fileErr_tbConDefs;
   goto cleanUp_fun15_sec04;

   cleanUp_fun15_sec04:;

   if(outFILE != stdout)
      fclose(outFILE);

   outFILE = 0;

   return errSC;
} /*pConBasAry*/

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
