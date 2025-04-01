/*-------------------------------------------------------\
' fluST SOF: Start Of File
'   - holds fluST structure and the supporting structures
'     for finding flu segment once primer location is
'     found
'   o header:
'     - included libraries and defined variables
'   o .h st01: segIdSearch_fluST
'     - is a single link in a flu segment search
'   o .h st02: fluST
'     - has the segement id list to search and segments
'   o fun01: blank_segIdSearch_fluST
'     - sets a segIdSearch_fluST structure to defaults
'   o fun02: init_segIdSearch_fluST
'     - initializes variables in a segIdSearch_fluST
'       structure
'   o fun03: freeStack_segIdSearch_fluST
'     - free variables in a segIdSearch_fluST structure
'     - This will free all nodes, except the previous node
'       in the linked list
'   o fun04: freeHeap_segIdSearch_fluST
'     - free variables in a segIdSearch_fluST structure
'     - This will free all nodes, except previous nodes
'       in the linked list
'   o .c fun05: addSeqTo_segIdSeach_fluST
'     - adds a new pattern (sequence) to a
'       segIdSearch_fluST linked list
'   o .c fun06: findSeq_segIdSeach_fluST
'     - finds a sequence in segIdSearch_fluST linked list
'   o .c fun07: cpCompSeq_fluST
'     - copies the complement sequence to a new string.
'       anonymous bases and '-'s are convereted to 'N's
'   o fun08: addSegTo_fluST
'     - adds a segment to a fluST structure
'   o fun09: rmSegSeqFrom_fluST
'     - removes a segment sequence from a fluST structure
'   o fun10: findSeg_fluST
'     - find the segement a sequence belongs to
'   o fun11: blank_fluST
'     - sets a fluST stucture to defaults
'   o fun12: init_fluST
'     - initializes all values in a fluST struct to 0
'   o fun13: setup_fluST
'     - allocates and sets up a fluST struture to default
'       values (call init_fluST first)
'   o fun14: freeStack_fluST
'     - frees all variables inside a fluST structure
'   o fun15: freeHeap_fluST
'     - frees a fluST structure
'   o fun16: detectDI_fluST
'     - finds segment number and then detects if read is
'       a diRNA or a vRNA (full)
'   o fun17: pidHeader_fluST
'     - prints out the header for pid_fluST
'   o fun18: pid_fluST
'     - prints out read id and stats for a flu segment
'   o license:
'     - licensing for this code (public domain / mit)
\-------------------------------------------------------*/

/*-------------------------------------------------------\
| Header:
|   - included libraries and defined variables
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

#include "fluST.h"

#include <stdio.h>
#include "../genLib/shellSort.h"

/*only .h files*/
#include "../genLib/dataTypeShortHand.h"
#include "fluSeg.h"

#include "../genAln/alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: blank_segIdSearch_fluST
|   - sets a segIdSearch_fluST structure to defaults
| Input:
|   - segIdSearch_fluST_ptr:
|     o pointer to a segIdSearch_fluST to set to defaults
| Output:
|   - Modifies:
|     o values in segIdSearch_fluST_ptr to be defaults
\-------------------------------------------------------*/
void
blank_segIdSearch_fluST(
   struct segIdSearch_fluST *segIdSearch_fluST_ptr
){
   unsigned char ucSegMac = 0;

   for(
      ucSegMac = 0;
      ucSegMac < def_numSeg_fluST;
      ++ucSegMac
   )
      segIdSearch_fluST_ptr->segArySC[ucSegMac] =
         def_noSeg_fluST;

   segIdSearch_fluST_ptr->numSegsUC = 0;
} /*blank_segIdSearch_fluST*/

/*-------------------------------------------------------\
| Fun02: init_segIdSearch_fluST
|   - initializes variables in a segIdSearch_fluST
|     structure
| Input:
|   - segIdSearch_fluST_ptr:
|     o pointer to a segIdSearch_fluST to initialize
| Output:
|   - Modifies:
|     o values in segIdSearch_fluST_ptr to be defaults or
|       null
\-------------------------------------------------------*/
void
init_segIdSearch_fluST(
   struct segIdSearch_fluST *segIdSearch_fluST_ptr
){
   segIdSearch_fluST_ptr->aNtST = 0;
   segIdSearch_fluST_ptr->cNtST = 0;
   segIdSearch_fluST_ptr->gNtST = 0;
   segIdSearch_fluST_ptr->tNtST = 0;

   segIdSearch_fluST_ptr->lastNtST = 0;

   blank_segIdSearch_fluST(segIdSearch_fluST_ptr);
} /*init_segIdSearch_fluST*/

/*-------------------------------------------------------\
| Fun03: freeStack_segIdSearch_fluST
|   - free variables in a segIdSearch_fluST structure
|   - This will free all nodes, except the previous node
|     in the linked list
| Input:
|   - segSTPtr:
|     o pointer to a segIdSearch_fluST to free
| Output:
|   - Frees:
|     o the A, C, G, and T linked lists in segSTPtr
|   - Modifies:
|     o sets segArySC to 0
|     o sets the A, C, G, and T pointers to 0
\-------------------------------------------------------*/
void
freeStack_segIdSearch_fluST( 
   struct segIdSearch_fluST *segSTPtr
){
   if(! segSTPtr)
      return;

   if(segSTPtr->aNtST)
   { /*If: have A nucleotide*/
      freeStack_segIdSearch_fluST(segSTPtr->aNtST);
      free(segSTPtr->aNtST);
      segSTPtr->aNtST = 0;
   } /*If: have A nucleotide*/

   if(segSTPtr->cNtST)
   { /*If: have C nucleotide*/
      freeStack_segIdSearch_fluST(segSTPtr->cNtST);
      free(segSTPtr->cNtST);
      segSTPtr->cNtST = 0;
   } /*If: have C nucleotide*/

   if(segSTPtr->gNtST)
   { /*If: have G nucleotide*/
      freeStack_segIdSearch_fluST(segSTPtr->gNtST);
      free(segSTPtr->gNtST);
      segSTPtr->gNtST = 0;
   } /*If: have G nucleotide*/

   if(segSTPtr->tNtST)
   { /*If: have T nucleotide*/
      freeStack_segIdSearch_fluST(segSTPtr->tNtST);
      free(segSTPtr->tNtST);
      segSTPtr->tNtST = 0;
   } /*If: have T nucleotide*/

   blank_segIdSearch_fluST(segSTPtr);
   return;
} /*freeStack_segIdSearch_fluST*/

/*-------------------------------------------------------\
| Fun04: freeHeap_segIdSearch_fluST
|   - free variables in a segIdSearch_fluST structure
|   - This will free all nodes, except previous nodes
|     in the linked list
| Input:
|   - segSTPtr:
|     o pointer to a segIdSearch_fluST to free
| Output:
|   - Frees:
|     o segSTPtr and all values in it
|       - you must set segSTPtr to null
|   - Modifies:
|     o segSTPtr->lastNtST to have the link to segSTPtr
|       set to null
\-------------------------------------------------------*/
void
freeHeap_segIdSearch_fluST( 
   struct segIdSearch_fluST *segSTPtr
){
   if(segSTPtr->aNtST)
   { /*If: have A nucleotide*/
      freeStack_segIdSearch_fluST(segSTPtr->aNtST);
      free(segSTPtr->aNtST);
      segSTPtr->aNtST = 0;
   } /*If: have A nucleotide*/

   if(segSTPtr->cNtST)
   { /*If: have C nucleotide*/
      freeStack_segIdSearch_fluST(segSTPtr->cNtST);
      free(segSTPtr->cNtST);
      segSTPtr->cNtST = 0;
   } /*If: have C nucleotide*/

   if(segSTPtr->gNtST)
   { /*If: have G nucleotide*/
      freeStack_segIdSearch_fluST(segSTPtr->gNtST);
      free(segSTPtr->gNtST);
      segSTPtr->gNtST = 0;
   } /*If: have G nucleotide*/

   if(segSTPtr->tNtST)
   { /*If: have T nucleotide*/
      freeStack_segIdSearch_fluST(segSTPtr->tNtST);
      free(segSTPtr->tNtST);
      segSTPtr->tNtST = 0;
   } /*If: have T nucleotide*/

   blank_segIdSearch_fluST(segSTPtr);

   /*remove any links to this structure*/
   if(segSTPtr->lastNtST)
   { /*If: had previous list*/
      if(segSTPtr->lastNtST->aNtST == segSTPtr)
         segSTPtr->lastNtST->aNtST = 0; /*remove link*/

      else if(segSTPtr->lastNtST->cNtST == segSTPtr)
         segSTPtr->lastNtST->cNtST = 0; /*remove link*/

      else if(segSTPtr->lastNtST->gNtST == segSTPtr)
         segSTPtr->lastNtST->gNtST = 0; /*remove link*/

      else if(segSTPtr->lastNtST->tNtST == segSTPtr)
         segSTPtr->lastNtST->tNtST = 0; /*remove link*/
   } /*If: had previous list*/

   free(segSTPtr); /*free this structure*/

   return;
} /*freeStack_segIdSearch_fluST*/

/*-------------------------------------------------------\
| Fun05: addSeqTo_segIdSeach_fluST
|   - adds a new pattern (sequence) to a segIdSearch_fluST
|     linked list
| Input:
|   - segSTPtr:
|     o pointer to segIdSearch_fluST to add new sequence 
|   - seqStr:
|     o c-string with sequence/pattern to add
|   - segArySC:
|     o segment id to add; use def_XXXNum_fluSeg variables
| Output:
|   - Modifies:
|     o segSTPtr linked list to have the new sequence
|       - last link in seqStr has segArySC 
|   - Returns:
|     o 0 for no errors
|     o 1 for memory error
|     o 2 for unkown base
\-------------------------------------------------------*/
signed char
addSeqTo_segIdSearch_fluST(
   struct segIdSearch_fluST *segSTPtr,
   signed char *seqStr,
   signed char segArySC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun05 TOC:
   '   - adds a new pattern (sequence) to a
   '     segIdSearch_fluST linked list
   '   o fun05 sec01:
   '     - check if base in sequence is 'A'
   '   o fun05 sec02:
   '     - check if base in sequence is 'C'
   '   o fun05 sec03:
   '     - check if base in sequence is 'G'
   '   o fun05 sec04:
   '     - check if base in sequence is 'T'
   '   o fun05 sec05:
   '     - add segment number if at end sequence or return
   '       unkown base error
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec01:
   ^   - check if base in sequence is 'A'
   ^   o fun05 sec01 sub01:
   ^     - check if is a 'A' and if 'A' already exists
   ^   o fun05 sec01 sub02:
   ^     - 'A' does not exists, make a 'A' structure
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec01 Sub01:
   *   - check if is a 'A' and if 'A' already exists
   \*****************************************************/

   if((*seqStr & ~32) == 'A')
   { /*If: this was A*/
      if(segSTPtr->aNtST)
      { /*If: A already exists*/
         return
            addSeqTo_segIdSearch_fluST(
               segSTPtr->aNtST,
               seqStr + 1,
               segArySC
            );
      } /*If: A already exists*/

      /**************************************************\
      * Fun05 Sec01 Sub02:
      *   - 'A' does not exists, make a 'A' structure
      \**************************************************/

      segSTPtr->aNtST = malloc(sizeof(segIdSearch_fluST));

      if(! segSTPtr->aNtST)
         return 1; /*memory error*/

      init_segIdSearch_fluST(segSTPtr->aNtST);
      segSTPtr->aNtST->lastNtST = segSTPtr;

      return
         addSeqTo_segIdSearch_fluST(
            segSTPtr->aNtST,
            seqStr + 1,
            segArySC
         );
   } /*If: this was A*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec02:
   ^   - check if base in sequence is 'C'
   ^   o fun05 sec02 sub01:
   ^     - check if is a 'C' and if 'C' already exists
   ^   o fun05 sec02 sub02:
   ^     - 'C' does not exists, make a 'C' structure
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec02 Sub01:
   *   - check if is a 'A' and if 'A' already exists
   \*****************************************************/

   if((*seqStr & ~32) == 'C')
   { /*If: this was C*/
      if(segSTPtr->cNtST)
      { /*If: C already exists*/
         return
            addSeqTo_segIdSearch_fluST(
               segSTPtr->cNtST,
               seqStr + 1,
               segArySC
            );
      } /*If: C already exists*/

      /**************************************************\
      * Fun05 Sec02 Sub02:
      *   - 'C' does not exists, make a 'C' structure
      \**************************************************/

      segSTPtr->cNtST = malloc(sizeof(segIdSearch_fluST));

      if(! segSTPtr->cNtST)
         return 1; /*memory error*/

      init_segIdSearch_fluST(segSTPtr->cNtST);
      segSTPtr->cNtST->lastNtST = segSTPtr;

      return
         addSeqTo_segIdSearch_fluST(
            segSTPtr->cNtST,
            seqStr + 1,
            segArySC
         );
   } /*If: this was C*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec03:
   ^   - check if base in sequence is 'G'
   ^   o fun05 sec03 sub01:
   ^     - check if is a 'G' and if 'G' already exists
   ^   o fun05 sec03 sub02:
   ^     - 'G' does not exists, make a 'G' structure
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec03 Sub01:
   *   - check if is a 'G' and if 'G' already exists
   \*****************************************************/

   if((*seqStr & ~32) == 'G')
   { /*If: this was G*/
      if(segSTPtr->gNtST)
      { /*If: G already exists*/
         return
            addSeqTo_segIdSearch_fluST(
               segSTPtr->gNtST,
               seqStr + 1,
               segArySC
            );
      } /*If: G already exists*/

      /**************************************************\
      * Fun05 Sec03 Sub02:
      *   - 'G' does not exists, make a 'G' structure
      \**************************************************/

      segSTPtr->gNtST = malloc(sizeof(segIdSearch_fluST));

      if(! segSTPtr->gNtST)
         return 1; /*memory error*/

      init_segIdSearch_fluST(segSTPtr->gNtST);
      segSTPtr->gNtST->lastNtST = segSTPtr;

      return
         addSeqTo_segIdSearch_fluST(
            segSTPtr->gNtST,
            seqStr + 1,
            segArySC
         );
   } /*If: this was G*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec04:
   ^   - check if base in sequence is 'T'
   ^   o fun05 sec04 sub01:
   ^     - check if is a 'T' and if 'T' already exists
   ^   o fun05 sec04 sub02:
   ^     - 'T' does not exists, make a 'T' structure
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec04 Sub01:
   *   - check if is a 'T' and if 'T' already exists
   \*****************************************************/

   if((*seqStr & ~32) == 'T')
   { /*If: this was T*/
      if(segSTPtr->tNtST)
      { /*If: T already exists*/
         return
            addSeqTo_segIdSearch_fluST(
               segSTPtr->tNtST,
               seqStr + 1,
               segArySC
            );
      } /*If: T already exists*/

      /**************************************************\
      * Fun05 Sec04 Sub02:
      *   - 'T' does not exists, make a 'T' structure
      \**************************************************/

      segSTPtr->tNtST = malloc(sizeof(segIdSearch_fluST));

      if(! segSTPtr->tNtST)
         return 1; /*memory error*/

      init_segIdSearch_fluST(segSTPtr->tNtST);
      segSTPtr->tNtST->lastNtST = segSTPtr;

      return
         addSeqTo_segIdSearch_fluST(
            segSTPtr->tNtST,
            seqStr + 1,
            segArySC
         );
   } /*If: this was T*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec05:
   ^   - add segment number if at end sequence or return
   ^     unkown base error
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(*seqStr == '\0')
   { /*If: this was a null*/
      segSTPtr->segArySC[segSTPtr->numSegsUC] = segArySC;

      ++segSTPtr->numSegsUC;

      if(segSTPtr->numSegsUC == 1)
      { /*If: have two segments*/
         if(
              (uchar) segSTPtr->segArySC[0]
            > (uchar) segSTPtr->segArySC[1])
         { /*If: need to swap segment numbers*/
            segSTPtr->segArySC[0] ^=segSTPtr->segArySC[1];
            segSTPtr->segArySC[1] ^=segSTPtr->segArySC[0];
            segSTPtr->segArySC[0] ^=segSTPtr->segArySC[1];
         } /*If: need to swap segment numbers*/
      } /*If: have two segments*/

      else if(segSTPtr->numSegsUC > 1)
      { /*Else If: need to sort segment numbers*/
         uc_shellSort(
            (uchar *) segSTPtr->segArySC,
            0,
            (ulong) segSTPtr->numSegsUC - 1
         );
      } /*Else If: need to sort segment numbers*/

      return 0;
   } /*If: this was a null*/

   return 2; /*no idea what this was*/
} /*addSeqTo_segIdSearch_fluST*/

/*-------------------------------------------------------\
| Fun06: findSeq_segIdSeach_fluST
|   - finds a sequence in a segIdSearch_fluST linked list
| Input:
|   - segSTPtr:
|     o pointer to segIdSearch_fluST to add new sequence 
|   - seqStr:
|     o c-string with sequence to search for
|   - revBl:
|     o 1: search backwards (reverse sequences)
|     o 0: search forwards (forward sequences)
| Output:
|   - Returns:
|     o 0 for error
|     o pointer to last node with the segment
\-------------------------------------------------------*/
segIdSearch_fluST *
findSeq_segIdSeach_fluST(
   struct segIdSearch_fluST *segSTPtr,
   signed char *seqStr,
   signed char revBl
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun06 TOC:
   '   - finds a sequence in segIdSearch_fluST linked list
   '   o fun06 sec01:
   '     - variable declarations
   '   o fun06 sec02:
   '     - find longest match to sequence
   '   o fun06 sec03:
   '     - see if any segments were detected in sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   struct segIdSearch_fluST *nodeST = segSTPtr;

   if(revBl)
      revBl = -1;
   else
      revBl = 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec02:
   ^   - find longest match to sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(*seqStr)
   { /*Loop: find last node mathcing sequence*/
       switch(*seqStr & ~32)
       { /*Switch: check base*/
          case 'A':
          /*Case: sequence has 'A'*/
             if(nodeST->aNtST)
                nodeST = nodeST->aNtST;
             else
                goto foundSeq_fun06_sec03;

             break;
          /*Case: sequence has 'A'*/

          case 'C':
          /*Case: sequence has 'C'*/
             if(nodeST->cNtST)
                nodeST = nodeST->cNtST;
             else
                goto foundSeq_fun06_sec03;

             break;
          /*Case: sequence has 'C'*/

          case 'G':
          /*Case: sequence has 'G'*/
             if(nodeST->gNtST)
                nodeST = nodeST->gNtST;
             else
                goto foundSeq_fun06_sec03;

             break;
          /*Case: sequence has 'G'*/

          case 'T':
          /*Case: sequence has 'T'*/
             if(nodeST->tNtST)
                nodeST = nodeST->tNtST;
             else
                goto foundSeq_fun06_sec03;

             break;
          /*Case: sequence has 'T'*/

          case 'U':
          /*Case: sequence has 'T' (RNA 'U')*/
             if(nodeST->tNtST)
                nodeST = nodeST->tNtST;
             else
                goto foundSeq_fun06_sec03;

             break;
          /*Case: sequence has 'T' (RNA 'U')*/

          default:
             goto foundSeq_fun06_sec03;
             /*no idea what base is*/
       } /*Switch: check base*/

      seqStr += revBl;
   } /*Loop: find last node mathcing sequence*/
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun06 Sec03:
   ^   - see if any segments were detected in sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   foundSeq_fun06_sec03:;

   while(nodeST->lastNtST)
   { /*Loop: see if any matches in the path*/
      if(nodeST->numSegsUC > 0)
         break; /*found a match*/
      else
         nodeST = nodeST->lastNtST;
   } /*Loop: see if any matches in the path*/

   return nodeST;
} /*findSeq_segIdSeach_fluST*/

/*-------------------------------------------------------\
| Fun07: cpCompSeq_fluST
|   - copies the complement sequence to a new string.
|     anonymous bases and '-'s are convereted to 'N's
| Input:
|   - dupStr:
|     o c-string to copy the string to
|   - seqStr:
|     o c-string with sequence to complement and copy
| Output:
|   - Modifies:
|     o dupStr to have complement sequence of seqStr.
\-------------------------------------------------------*/
void
cpCompSeq_fluST(
   signed char *dupStr,
   signed char *seqStr
){
   signed short posMacSS = 0;
   signed short dupMacSS = 0;

   while(seqStr[posMacSS] != '\0' )
   { /*Loop: copy complement sequence*/
      switch(seqStr[posMacSS] & ~32)
      { /*Switch: copy complement base*/
         case 'A':
            dupStr[dupMacSS] = 'T';
            break;

         case 'C':
            dupStr[dupMacSS] = 'G';
            break;

         case 'G':
            dupStr[dupMacSS] = 'C';
            break;

         case 'T':
            dupStr[dupMacSS] = 'A';
            break;

         case 'U':
            dupStr[dupMacSS] = 'A';
            break;

         default:
            dupStr[dupMacSS] = 'N'; /*no idea*/
            break;
      } /*Switch: copy complement base*/

      ++posMacSS;
      ++dupMacSS;
   } /*Loop: copy complement sequence*/

   dupStr[dupMacSS] = '\0';
} /*cpCompSeq_fluST*/

/*-------------------------------------------------------\
| Fun08: addSegTo_fluST
|   - adds a segment to a fluST structure
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to add segment to
|   - seqStr:
|     o c-string with the sequence to add
|   - lenSegSI:
|     o expected length of the segment
|   - revPrimBl:
|     o direction of primer (1 = reverse; 0 = foward)
|   - segSC:
|     o segment id to add (use def_XXXNum_fluSeg variables)
| Output:
|   - Modifies:
|     o all variables in fluST to support the new pattern
|   - Returns:
|     o 0 for no errors
|     o 1 for memory error
|     o 2 for invalid base in a seqStr
\-------------------------------------------------------*/
signed char
addSegTo_fluST(
   struct fluST *fluSTPtr,
   signed char *seqStr,
   signed int lenSegSI,
   signed char revPrimBl,
   signed char segSC
){
   signed char errSC = 0;         /*for error checking*/
   signed char compSeqStr[4096];  /*compementing seqStr*/

   fluSTPtr->lenSegArySI[segSC] = lenSegSI;

   /*add in complement sequence*/
   cpCompSeq_fluST(
      compSeqStr,
      seqStr
   );

   /*add in the provided sequence*/
   if(revPrimBl)
   { /*If: is a reverse primer*/
      errSC =
         addSeqTo_segIdSearch_fluST(
            fluSTPtr->revSearchST,
            seqStr,
            segSC
         );

      if(errSC)
         return errSC;

      errSC =
         addSeqTo_segIdSearch_fluST(
            fluSTPtr->revCompSearchST,
            compSeqStr,
            segSC
         );
   } /*If: is a reverse primer*/

   else
   { /*Else: is a foward primer*/
      errSC =
         addSeqTo_segIdSearch_fluST(
            fluSTPtr->forSearchST,
            seqStr,
            segSC
         );

      if(errSC)
         return errSC;

      errSC =
         addSeqTo_segIdSearch_fluST(
            fluSTPtr->forCompSearchST,
            compSeqStr,
            segSC
         );
   } /*Else: is a foward primer*/

   if(errSC)
      return errSC;

   return 0;
} /*addSegTo_fluST*/

/*-------------------------------------------------------\
| Fun09: rmSegSeqFrom_fluST
|   - removes a segment sequence from a fluST structure
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to remove segment
|       sequence
|   - seqStr:
|     o c-string with the sequence to remove
|   - revPrimBl:
|     o direction of primer (1 = reverse; 0 = foward)
|   - segArySC:
|     o segment id to add (use def_XXXNum_fluSeg variables)
| Output:
|   - Modifies:
|     o removes the sequence from the segment linked list
|   - Returns:
|     o 0 for no errors
|     o 1 for sequence not in list
\-------------------------------------------------------*/
signed char
rmSegSeqFrom_fluST(
   struct fluST *fluSTPtr,
   signed char *seqStr,
   signed char revPrimBl,
   signed char rmSegNumSC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun09 TOC:
   '   - removes a segment sequence from a fluST structure
   '   o fun09 sec01:
   '     - variable decleraions
   '   o fun09 sec02:
   '     - find segment
   '   o fun09 sec03:
   '     - delete (lazy) segment id
   '   o fun09 sec04:
   '     - find complement sequence segment id
   '   o fun09 sec05:
   '     - remove complement sequence segment id
   '   o fun09 sec06:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec01:
   ^   - variable decleraions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar noSegBl = 0;
   uchar ucSeg = 0;
   schar compSeqStr[4096];  /*compementing seqStr*/
   struct segIdSearch_fluST *segIdST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec02:
   ^   - find segment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(revPrimBl)
      segIdST = fluSTPtr->revSearchST;
   else
      segIdST = fluSTPtr->forSearchST;

   segIdST =
      findSeq_segIdSeach_fluST(
         segIdST,
         seqStr,
         0       /*is foward direction*/
      ); /*find the pattern*/
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec03:
   ^   - delete (lazy) segment id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(segIdST)
   { /*If: had a sequence match*/
      for(
         ucSeg = 0;
         ucSeg < segIdST->numSegsUC;
         ++ucSeg
      ){ /*Loop: find and remove the id*/
         if(rmSegNumSC == segIdST->segArySC[ucSeg])
         { /*If: found the segment*/
            segIdST->segArySC[ucSeg] = def_noSeg_fluST;
            break;
         } /*If: found the segment*/
      } /*Loop: find and remove the id*/

      if(ucSeg == segIdST->numSegsUC)
         noSegBl = 1;

      while(ucSeg < segIdST->numSegsUC - 1)
      { /*Loop: move keep segments back*/
            segIdST->segArySC[ucSeg] ^=
            segIdST->segArySC[ucSeg + 1];

            segIdST->segArySC[ucSeg + 1] ^=
            segIdST->segArySC[ucSeg];

            segIdST->segArySC[ucSeg] ^=
            segIdST->segArySC[ucSeg + 1];
      } /*Loop: move keep segments back*/
   } /*If: had a sequence match*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec04:
   ^   - find complement sequence segment id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*add in complement sequence*/
    cpCompSeq_fluST(
      compSeqStr,
      seqStr
    );

   if(revPrimBl)
      segIdST = fluSTPtr->revCompSearchST;
   else
      segIdST = fluSTPtr->forCompSearchST;

   segIdST =
      findSeq_segIdSeach_fluST(
         segIdST,
         compSeqStr,
         0  /*should be forward*/
      ); /*find the pattern*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec05:
   ^   - remove complement sequence segment id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(segIdST)
   { /*If: had a sequence match*/
      for(
         ucSeg = 0;
         ucSeg < segIdST->numSegsUC;
         ++ucSeg
      ){ /*Loop: find and remove the id*/
         if(rmSegNumSC == segIdST->segArySC[ucSeg])
         { /*If: found the segment*/
            segIdST->segArySC[ucSeg] = def_noSeg_fluST;
            break;
         } /*If: found the segment*/
      } /*Loop: find and remove the id*/

      while(ucSeg < segIdST->numSegsUC - 1)
      { /*Loop: move keep segments back*/
         segIdST->segArySC[ucSeg] ^=
            segIdST->segArySC[ucSeg + 1];

         segIdST->segArySC[ucSeg + 1] ^=
            segIdST->segArySC[ucSeg];

         segIdST->segArySC[ucSeg] ^=
            segIdST->segArySC[ucSeg + 1];
      } /*Loop: move keep segments back*/
   } /*If: had a sequence match*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun09 Sec06:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return noSegBl;
} /*rmSegSeqFrom_fluST*/

/*-------------------------------------------------------\
| Fun10: findSeg_fluST
|   - find the segement a sequence belongs to
| Input:
|   - fluSTPtr:
|     o pointer to fluST structure to search for segment
|       in
|   - seqStr:
|     o c-sting pointing to first base that is part of
|       segment pattern in the sequence
|   - cmpBl:
|     o 1: sequence is a complement sequence
|     o 0: sequence is normal
|   - revBl:
|     o 1: is a reverse primer
|     o 0: forward primer
|   - numSegsUCPtr:
|     o pointer to unsigned char to hold number of
|       segments with the pattern
| Output:
|   - Modifies:
|     o numSegUCPtr to have the number of segments in the
|       returned array
|   - Returns:
|     o pointer to array with the segment number
|     o 0 if no segment found
\-------------------------------------------------------*/
signed char *
findSeg_fluST(
   struct fluST *fluSTPtr,
   signed char *seqStr,
   signed char cmpBl,
   signed char revBl,
   unsigned char *numSegsUCPtr
){
   struct segIdSearch_fluST *nodeST = 0;

   if(cmpBl)
   { /*If: sequence is reverse complement*/
      if(revBl)
         nodeST = fluSTPtr->revCompSearchST;
      else
         nodeST = fluSTPtr->forCompSearchST;

      nodeST =
         findSeq_segIdSeach_fluST(
            nodeST,
            seqStr,
            1
         );
   } /*If: sequence is reverse complement*/

   else
   { /*Else: sequence is normal*/
      if(revBl)
         nodeST = fluSTPtr->revSearchST;
      else
         nodeST = fluSTPtr->forSearchST;

      nodeST =
         findSeq_segIdSeach_fluST(
            nodeST,
            seqStr,
            0
         );
   } /*Else: sequence is normal*/

   if(! nodeST)
   { /*If: could not find the segment*/
      *numSegsUCPtr = 0;
      return 0;
   } /*If: could not find the segment*/

   *numSegsUCPtr = nodeST->numSegsUC;
   return nodeST->segArySC;
} /*findSeg_fluST*/

/*-------------------------------------------------------\
| Fun11: blank_fluST
|   - sets a fluST stucture to defaults
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to blank
| Output:
|   - Modifies:
|     o id/mvRNA detection thresholds
\-------------------------------------------------------*/
void
blank_fluST(
   struct fluST *fluSTPtr
){
   fluSTPtr->minDiDelSI = def_diMinDel_fluST;
   fluSTPtr->maxMvLenSI = def_maxMvRNALen_fluST;
} /*blank_fluST*/

/*-------------------------------------------------------\
| Fun12: init_fluST
|   - initializes all values in a fluST struct to 0
| Input:
|   - fluSTPtr:
|     o pointer to fluST structure to initialize
| Output:
|   - Modifies:
|     o all values in fluSTPtr to be defaults or 0
\-------------------------------------------------------*/
void
init_fluST(
  struct fluST *fluSTPtr
){
   sshort ssSeg = 0;

   fluSTPtr->forSearchST = 0;
   fluSTPtr->forCompSearchST = 0;
   fluSTPtr->revSearchST = 0;
   fluSTPtr->revCompSearchST = 0;

   for(
      ssSeg = 0;
      ssSeg < def_numSeg_fluST;
      ++ssSeg
   ){ /*Loop: blank segments*/
      fluSTPtr->segIdAryStr[ssSeg][0] = '\0';
      fluSTPtr->lenSegArySI[ssSeg] = 0;
   } /*Loop: blank segments*/

   blank_fluST(fluSTPtr);
} /*init_fluST*/

/*-------------------------------------------------------\
| Fun13: setup_fluST
|   - allocates and sets up a fluST struture to default
|     values (call init_fluST first)
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to initialize
| Output:
|   - Modifies:
|     o all values in fluSTPtr to be defaults.
|   - Returns:
|     o 0 for no errors
|     o 1 for memory error
|     o 2 for invalid base in a seqStr
\-------------------------------------------------------*/
signed char
setup_fluST(
  struct fluST *fluSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   '   o fun13 sec01:
   '     - variable declerations
   '   o fun13 sec02:
   '     - memory assignment for search arrays
   '   o fun13 sec03:
   '     - set up flu segments (search, id, and numbers)
   '   o fun13 sec04:
   '     - call blank function (for future)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   schar *tmpStr = 0;
   sshort ssSeg = 0;
   sshort ssId = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec02:
   ^   - memory assignment for search arrays
   ^   o fun13 sec02 sub01:
   ^     - forward search list head memory assignment
   ^   o fun13 sec02 sub02:
   ^     - reverse search list head memory assignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun13 Sec02 Sub01:
   *   - forwad search list head memory assignment
   \*****************************************************/

   if(fluSTPtr->forSearchST)
      free(fluSTPtr->forSearchST);

   fluSTPtr->forSearchST = 0;

   if(fluSTPtr->forCompSearchST)
      free(fluSTPtr->forCompSearchST);

   fluSTPtr->forCompSearchST = 0;

   fluSTPtr->forSearchST =
      malloc(sizeof(struct segIdSearch_fluST));

   if(! fluSTPtr->forSearchST)
      return 1;

   init_segIdSearch_fluST(fluSTPtr->forSearchST);

   fluSTPtr->forCompSearchST =
      malloc(sizeof(struct segIdSearch_fluST));

   if(! fluSTPtr->forCompSearchST)
      return 1;

   init_segIdSearch_fluST(fluSTPtr->forCompSearchST);

   /*****************************************************\
   * Fun13 Sec02 Sub02:
   *   - reverse search list head memory assignment
   \*****************************************************/

   if(fluSTPtr->revSearchST)
      free(fluSTPtr->revSearchST);

   fluSTPtr->revSearchST = 0;

   if(fluSTPtr->revCompSearchST)
      free(fluSTPtr->revCompSearchST);

   fluSTPtr->revCompSearchST = 0;

   fluSTPtr->revSearchST =
      malloc(sizeof(struct segIdSearch_fluST));

   if(! fluSTPtr->revSearchST)
      return 1;

   init_segIdSearch_fluST(fluSTPtr->revSearchST);

   fluSTPtr->revCompSearchST =
      malloc(sizeof(struct segIdSearch_fluST));

   if(! fluSTPtr->revCompSearchST)
      return 1;

   init_segIdSearch_fluST(fluSTPtr->revCompSearchST);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec03:
   ^   - set up flu segments (search, id, and numbers)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      ssSeg = 0;
      ssSeg < def_numSeg_fluST;
      ++ssSeg
   ){ /*Loop: add segments to search list*/
      errSC =
         addSegTo_fluST(
            fluSTPtr,
            forSeqAryStr_fluSeg[ssSeg],
            segLenArySS_fluSeg[ssSeg],
            0,
            ssSeg
         ); /*add forward sequence to search*/

      if(errSC)
         return errSC;

      errSC =
         addSegTo_fluST(
            fluSTPtr,
            revSeqAryStr_fluSeg[ssSeg],
            segLenArySS_fluSeg[ssSeg],
            1,
            ssSeg
         ); /*add reverse sequence to search*/

      /*copy the segment id*/
      tmpStr = fluSTPtr->segIdAryStr[ssSeg];
      ssId = 0;

      while(segIdAryStr_fluSeg[ssSeg][ssId] != '\0')
      { /*Loop: copy segment id*/
         tmpStr[ssId] = segIdAryStr_fluSeg[ssSeg][ssId];
         ++ssId;
      } /*Loop: copy segment id*/

      tmpStr[ssId] = '\0';
   } /*Loop: add segments to search list*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec04:
   ^   - call blank function (for future)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_fluST(fluSTPtr);

   return 0;
} /*setup_fluST*/

/*-------------------------------------------------------\
| Fun14: freeStack_fluST
|   - frees all variables inside a fluST structure
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to free variables in
| Output:
|   - Frees:
|     o all variables in fluSTPtr and sets pointers to 0
\-------------------------------------------------------*/
void
freeStack_fluST(
   struct fluST *fluSTPtr
){
   if(! fluSTPtr)
      return;

   /*free forward search structures*/
   if(fluSTPtr->forSearchST)
     freeHeap_segIdSearch_fluST(fluSTPtr->forSearchST);
   fluSTPtr->forSearchST = 0;

   if(fluSTPtr->forCompSearchST)
      freeHeap_segIdSearch_fluST(
         fluSTPtr->forCompSearchST
      );
   fluSTPtr->forCompSearchST = 0;

   /*free reverse search structures*/
   if(fluSTPtr->revSearchST)
     freeHeap_segIdSearch_fluST(fluSTPtr->revSearchST);
   fluSTPtr->revSearchST = 0;

   if(fluSTPtr->revCompSearchST)
      freeHeap_segIdSearch_fluST(
         fluSTPtr->revCompSearchST
      );

   fluSTPtr->revCompSearchST = 0;
} /*freeStack_fluST*/

/*-------------------------------------------------------\
| Fun15: freeHeap_fluST
|   - frees a fluST structure
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to free
| Output:
|   - Frees:
|     o fluSTPtr (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_fluST(
   struct fluST *fluSTPtr
){
   freeStack_fluST(fluSTPtr);
} /*freeStack_fluST*/

/*-------------------------------------------------------\
| Fun16: detectDI_fluST
|   - finds segment number and then detects if read is
|     a diRNA or a vRNA (full)
| Input:
|   - fluSTPtr:
|     o pointer to fluST structure with thresholds and
|       segment number
|   - seqStr:
|     o c-string with sequence to see if is diRNA
|   - hitAryUI:
|     o unsigned in array with number of hits for both
|       primers
|   - dirArySC:
|     o signed char array with the direction "F" for
|       foward of the best mapped primer
|   - seqStartAryUL:
|     o unsigned long array with the primers starting
|       coordinate on the sequence
|   - seqEndAryUL:
|     o unsigned long array with the primers ending
|       coordinate on the sequence
|   - segSC:
|     o pointer to signed char to hold the segment number
|       decided on
|   - mappedLenUL:
|     o pionter to unsigned long to hold the length of
|       region from start to end of primers
| Output:
|   - Modifies:
|     o segArySC to have the segment number
|   - Returns:
|     o 0 if both primers were not detected, if primers
|       are in same direction, or if no segments were
|       found or if have multiple segments for one primer
|     o segement stats | genome type
|       - segment stats:
|         o 0 for no segment found
|         o def_segFound_fluST | def_revSegSup_fluST if
|           both primers support the same segment
|         o def_partSeg_fluST if only forward primer
|           supported a segment
|         o def_partSeg_fluST | def_revSegSup_fluST if
|           only reverse primer supported a segment
|         o def_noSeg_fluST if neither primer supports a
|           segment
|         o def_multiSeg_fluST if more than one segment
|           was supported
|       - genome type:
|         o def_diFound_fluST if meets diRNA length
|         o def_mvFound_fluST if is mvRNA length
|         o def_fullFound_fluST if is larger then diRNA
|           length (full genome)
\-------------------------------------------------------*/
signed char
detectDI_fluST(
   struct fluST *fluSTPtr,
   signed char *seqStr,
   unsigned int hitAryUI[],       /*primer hits*/
   signed char dirArySC[],        /*primer direction*/
   unsigned long seqStartAryUL[], /*seq map coordiantes*/
   unsigned long seqEndAryUL[],   /*seq map coordiantes*/
   signed char *segSC,
   unsigned long *mappedLenUL
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun16 TOC:
   '   - find segment number and if is di sequence
   '   o fun16 sec01:
   '     - variable declarations
   '   o fun16 sec02:
   '     - check if have primers in oppisite directions
   '   o fun16 sec03:
   '     - find segment type
   '   o fun16 sec04:
   '     - check if have full segment
   '   o fun16 sec05:
   '     - check if di/mv/fullRNA
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar retSC = 0;

   /*primer segment ids*/
   schar *forArySC = 0;
   schar *revArySC = 0;
   uchar numForSegUC = 0;
   uchar numRevSegUC = 0;

   /*narrowing down segment id*/
   schar keepArySC[def_numSeg_fluST];
   uchar numSegsUC = 0;
   uchar ucFor = 0;
   uchar ucRev = 0;

   unsigned long lenSeqUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec02:
   ^   - check if have both primers in oppisite directions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *segSC = def_noSeg_fluST;/*chagned if find segment*/

   if(hitAryUI[0] == 0)
      return 0; /*no support for foward primer*/

   if(hitAryUI[1] == 0)
      return 0; /*no support for ending primer*/

   if(dirArySC[0] == dirArySC[1])
      return 0; /*primers are same direction*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec03:
   ^   - find segment type
   ^   o fun16 sec03 sub01:
   ^     - find the forward primers segment type
   ^   o fun16 sec03 sub02:
   ^     - find the reverse primers segment type
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun16 Sec03 Sub01:
   *   - find the forward primers segment type
   \*****************************************************/

   if(dirArySC[0] == 'F')
   { /*If: foward primer mapped forward*/
      forArySC =
         findSeg_fluST(
            fluSTPtr,
            &seqStr[seqEndAryUL[0] + 1],
            0, /*not complement*/
            0, /*not reverse primer*/
            &numForSegUC
         );
   } /*If: foward primer mapped forward*/

   else
   { /*Else: foward primer mapped reverse*/
      forArySC =
         findSeg_fluST(
            fluSTPtr,
            &seqStr[seqStartAryUL[0] - 1],
            1, /*is complement*/
            0, /*not reverse primer*/
            &numForSegUC
         );
   } /*Else: foward primer mapped reverse*/

   if(seqEndAryUL[1] > seqStartAryUL[0])
      *mappedLenUL = seqEndAryUL[1] - seqStartAryUL[0];

   else
      *mappedLenUL = seqEndAryUL[0] - seqStartAryUL[1];
    
   /*****************************************************\
   * Fun16 Sec03 Sub02:
   *   - find the reverse primers segment type
   \*****************************************************/

   if(dirArySC[1] == 'F')
   { /*If: reverse primer mapped forward*/
      revArySC =
         findSeg_fluST(
            fluSTPtr,
            &seqStr[seqEndAryUL[1] + 1],
            0, /*not complement*/
            1, /*is reverse primer*/
            &numRevSegUC
         );
   } /*If: reverse primer mapped forward*/

   else
   { /*Else: foward primer mapped reverse*/
      revArySC =
         findSeg_fluST(
            fluSTPtr,
            &seqStr[seqStartAryUL[1] - 1],
            1, /*is complement*/
            1, /*is reverse primer*/
            &numRevSegUC
         );
   } /*Else: foward primer mapped reverse*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec04:
   ^   - check for one primer supported segmet ids
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numRevSegUC < 1)
   { /*If: only forward supports id*/
      if(numForSegUC != 1)
         return 0; /*can not get id*/

      *segSC = *forArySC;
      retSC = def_partSeg_fluST;
   } /*If: only forward supports id*/

   else if(numForSegUC < 1)
   { /*If: only reverse supports id*/
      if(numRevSegUC != 1)
         return 0; /*can not get id*/

      *segSC = *revArySC;
      retSC = def_revSegSup_fluST | def_partSeg_fluST;
   } /*If: only reverse supports id*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec05:
   ^   - reduce segment ids to one possible id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   else
   { /*Else: have two primers supporting a segment*/
      retSC = 0;
      ucFor = 0;
      ucRev = 0;
      numSegsUC = 0;

      while(
            ucFor < numForSegUC
         && ucRev < numRevSegUC
      ){ /*Loop: narrow down segment ids*/
         
         if(forArySC[ucFor] < 0)
            break; /*no more ids*/

         if(revArySC[ucRev] < 0)
            break; /*no more ids*/

         if(forArySC[ucFor] < revArySC[ucRev])
            ++ucFor;

         else if(forArySC[ucFor] > revArySC[ucRev])
            ++ucRev;

         else
         { /*Else: same segment*/
            keepArySC[numSegsUC] = forArySC[ucFor];
            ++numSegsUC;

            ++ucFor;
            ++ucRev;
         } /*Else: same segment*/
      } /*Loop: narrow down segment ids*/

      if(numSegsUC > 1)
         return def_multiSeg_fluST; /*to many segments*/

      else if(numSegsUC)
      { /*Else If: have a segment*/
         *segSC = keepArySC[0];
         retSC = def_segFound_fluST | def_revSegSup_fluST;
      } /*Else If: have a segment*/

      else
         return def_multiSeg_fluST; /*primers disagree*/
   } /*Else: have two primers supporting a segment*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec06:
   ^   - check if di/mv/fullRNA
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   lenSeqUL =
        fluSTPtr->lenSegArySI[*segSC]
      - fluSTPtr->minDiDelSI;

   if(*mappedLenUL < lenSeqUL)
      return retSC | def_diFound_fluST;

   return retSC | def_fullFound_fluST; /*at least full*/
} /*dectectDI_fluST*/

/*-------------------------------------------------------\
| Fun17: pidHeader_fluST
|   - prints out the header for pid_fluST
| Input:
|   - outFILE:
|     o file to print header to
| Output:
|   - Prints:
|     o header to outFILE
\-------------------------------------------------------*/
void
pidHeader_fluST(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "id\tsegment\trnaType\tmap_len\tread_len\tseg_len"
   );

   fprintf(
      (FILE *) outFILE,
      "\tseg_id_prim\tfor_score\tfor_max_score\tfor_dir"
   );

   fprintf(
      (FILE *) outFILE,
      "\tfor_seq_start\tfor_seq_end"
   );

   fprintf(
      (FILE *) outFILE,
      "\tfor_prim_start\tfor_prim_end"
   );

   fprintf(
      (FILE *) outFILE,
      "\trev_score\trev_max_score\tref_dir"
   );

   fprintf(
      (FILE *) outFILE,
      "\trev_seq_start\trev_seq_end\trev_prim_start"
   );
 
   fprintf(
      (FILE *) outFILE,
      "\trev_prim_end\n"
   );
} /*pidHeader_fluST*/

/*-------------------------------------------------------\
| Fun18: pid_fluST
|   - prints out read id and stats for a flu segment
| Input:
|   - fluSTPtr:
|     o pointer to a fluSTPtr structure with lengths and
|       segment names
|   - idStr:
|     o c-string with read id
|   - segArySC:
|     o segment number found
|   - forSegSC:
|     o segment number of the foward primer
|   - revSegSC:
|     o segment number of the reverse primer
|   - diFlagSC:
|     o return from detect_DI_fluST (fun16)
|   - dirArySC:
|     o array of primer directions
|   - scoreArySL:
|     o array of alignment scores
|   - maxForScoreF:
|     o maximum score for foward primer
|   - maxRevScoreF:
|     o maximum score for reverse primer
|   - seqStartAryUL:
|     o array with primer starting positions on sequence
|   - seqEndAryUL:
|     o array with primer ending positions on sequence
|   - primStartAryUL:
|     o array with sequence starting positions on primer
|   - primEndAryUL:
|     o array with sequence ending positions on primer
|   - seqLenUL:
|     o length of sequence
|   - mappedLenUL:
|     o length of region between primers (start to end)
|   - outFILE:
|     o file to print id and stats to
| Output:
|   - Prints:
|     o id and other stats to outFILE
\-------------------------------------------------------*/
void
pid_fluST(
   struct fluST *fluSTPtr,
   signed char *idStr,
   signed char segArySC,
   signed char diFlagSC,
   signed char dirArySC[],
   signed long scoreArySL[],
   float maxForScoreF,
   float maxRevScoreF,
   unsigned long seqStartAryUL[],
   unsigned long seqEndAryUL[],
   unsigned long primStartAryUL[],
   unsigned long primEndAryUL[],
   unsigned long seqLenUL,
   unsigned long mappedLenUL,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun18 TOC:
   '   - print out ids for flu sequences
   '   o fun18 sec01:
   '     - check if this was a flu segment
   '   o fun18 sec02:
   '     - print read id and segment type
   '   o fun18 sec03:
   '     - print out the RNA type
   '   o fun18 sec04:
   '     - print out lengths and if both primers supported
   '   o fun18 sec05:
   '     - print foward primer stats
   '   o fun18 sec06:
   '     - print reverse primer stats
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec01:
   ^   - check if this was a flu segment and set up id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   
   signed char *nullSCPtr = idStr;
   signed char *endStr = 0;
   signed char *tmpStr = idStr;
   signed char nullSC = 0;

   /*get truncate id at first white space*/
   while(*nullSCPtr++ > 32) ;

   --nullSCPtr;
   nullSC = *nullSCPtr;
   *nullSCPtr = '\0';
   endStr = nullSCPtr;
   nullSCPtr = idStr;

   if(
         *nullSCPtr == '>'
      || *nullSCPtr == '@'
   ) ++nullSCPtr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec02:
   ^   - print read id and segment type
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "%s\t%s",
      nullSCPtr, /*read id; after header (> or @) symbol*/
      fluSTPtr->segIdAryStr[segArySC]
   ); /*print out read id and segment id*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec03:
   ^   - print out the RNA type
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(diFlagSC & def_diFound_fluST)
      tmpStr = (schar *) "diRna";

   else
      tmpStr = (schar *) "vRna";

   fprintf(
      (FILE *) outFILE,
      "\t%s",
      tmpStr
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec04:
   ^   - print out lengths and if both primers supported
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "\t%lu\t%lu\t%i",
      mappedLenUL,
      seqLenUL,
      fluSTPtr->lenSegArySI[segArySC]
   ); /*print expected length for full length segment*/

   if(diFlagSC & def_partSeg_fluST)
   { /*If: only primer supported a segment*/
      if(diFlagSC & def_revSegSup_fluST)
         fprintf(
            (FILE *) outFILE,
            "\trev"
         ); /*only reverse primer supported segment*/
      else
         fprintf(
            (FILE *) outFILE,
            "\tfor"
         ); /*only foward primer supported segment*/
   } /*If: only primer supported a segment*/

   else
      fprintf(
         (FILE *) outFILE,
         "\tboth"
      ); /*both primers supported the segment*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec05:
   ^   - print foward primer stats
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "\t%0.2f\t%0.2f\t%c",
      (float) scoreArySL[0] / def_scoreAdj_alnDefs,
      (float) maxForScoreF / def_scoreAdj_alnDefs,
      dirArySC[0]
   );

   fprintf(
      (FILE *) outFILE,
      "\t%lu\t%lu\t%lu\t%lu",
      seqStartAryUL[0] + 1,
      seqEndAryUL[0] + 1,
      primStartAryUL[0] + 1,
      primEndAryUL[0] + 1
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec06:
   ^   - print reverse primer stats
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "\t%0.2f\t%0.2f\t%c",
      (float) scoreArySL[1] / def_scoreAdj_alnDefs,
      (float) maxRevScoreF / def_scoreAdj_alnDefs,
      dirArySC[1]
   );

   fprintf(
      (FILE *) outFILE,
      "\t%lu\t%lu\t%lu\t%lu\n",
      seqStartAryUL[1] + 1,
      seqEndAryUL[1] + 1,
      primStartAryUL[1] + 1,
      primEndAryUL[1] + 1
   );

   *endStr = nullSC; /*revert character back*/
} /*pid_fluST*/

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
