/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' varFilt SOF:
'   - has functions for filtering reads with variants
'   o header:
'     - included libraries
'   o fun02: getProf_varFilt
'     - builds a profile for filtering
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PLAN9
   #include <u.h>
   #include <libc.h
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../genLib/strAry.h"
#include "../genBio/samEntry.h"
#include "../genBio/tbCon.h"

#include "../genBio/tbConDefs.h"
#include "../genLib/genMath.h"
#include "../genLib/dataTypeShortHand.h"

#define def_memErr_varFilt 1
#define def_fileErr_varFilt 2

/*for filtering variants*/
#define def_rmA_varFilt 1
#define def_rmC_varFilt 2
#define def_rmG_varFilt 4
#define def_rmT_varFilt 8

/*-------------------------------------------------------\
| Fun0x: conFilt_varFilt
|   - filters reads by consensus
\-------------------------------------------------------*/
void
conFilt_varFilt(
   float minSupF,                 /*min % sup to filter*/
   float minVarF,            /*min % var to discard read*/
   struct set_tbCon *conSetSTPtr, /*consensus settings*/
   struct samEntry *samSTPtr,     /*for reading sam file*/
   signed char **headStrPtr,      /*holds header*/
   unsigned long *lenHeadULPtr,   /*length of headStrPtr*/
   signed char **buffStrPtr,      /*for reading sam file*/
   unsigned long *lenBuffULPtr,   /*for reading sam file*/
   void *samFILE,                 /*sam file with reads*/
   void *outFILE,                 /*file to save reads*/
   signed char *errSCPtr          /*holds errors*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun0x TOC:
   '   o fun0x sec01:
   '     - variable declarations
   '   o fun0x sec02:
   '     - get reference lengths
   '   o fun0x sec03:
   '     - allocate memory
   '   o fun0x sec04:
   '     - build consensuses
   '   o fun0x sec05:
   '     - collapse consensuses
   '   o fun0x sec06:
   '     - detect low scoring variants
   '   o fun0x sec07:
   '     - filter out reads with low scoring variants
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun0x Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   struct refs_samEntry refStackST;

   struct conNt_tbCon **conHeapAryST = 0;
   struct conNt_tbCon *conOnST = 0; /*sanity*/
   unsigned int *lenConAryUI = 0; /*consensus lengths*/
   unsigned int uiCon = 0;        /*position in con/ref*/
   unsigned int uiRead = 0;       /*position in read*/

   signed char qBl = 0;

   unsigned int cigIndexUI = 0;
   unsigned int totalVarUI = 0;/*number filter variants*/
   unsigned int readVarUI = 0; /*# read filter variants*/
   unsigned int mapVarUI = 0;  /*number mapped variants*/
   signed int cigValSI = 0;
   
   signed long refIndexSL = 0;    /*index of reference*/

   /*for finding variant percent support*/
   float aPercSupF = 0;
   float cPercSupF = 0;
   float gPercSupF = 0;
   float tPercSupF = 0;
   float bestPercSupF = 0; /*best support*/

   float percCoverF = 0;
   float percVarF = 0;

   schar **rmArySCPtr = 0; /*array of variant arrays*/
   schar *rmArySC = 0;     /*sanity for rmArySCPtr*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun0x Sec02:
   ^   - get reference lengths
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_refs_samEntry(&refStackST);

   *errSCPtr =
      setup_refs_samEntry(
         &refStackST,
         16
      ); /*start out with 16 refs*/

   if(*errSCPtr)
     goto memErr_fun0x_sec0x_sub02;

   *errSCPtr =
     getRefLen_samEntry(
        &refStackST,
        samSTPtr,
        buffStrPtr,
        lenBuffULPtr,
        samFILE,
        0,
        headStrPtr,
        lenHeadULPtr
     );

   if(*errSCPtr)
   { /*If: had error*/
      if(*errSCPtr == def_memErr_samEntry)
         goto memErr_fun0x_sec0x_sub02;
      else
         goto fileErr_fun0x_sec0x_sub03;
   } /*If: had error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun0x Sec03:
   ^   - allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   conHeapAryST =
      malloc(
           refStackST.numRefUI
         * sizeof(struct conNt_tbCon *)
      ); /*allocate memory for uncollapsed consensus*/

   if(! conHeapAryST)
      goto memErr_fun0x_sec0x_sub02;

   for(
      refIndexSL = 0;
      refIndexSL < (slong) refStackST.numRefUI;
      ++refIndexSL
   ) conHeapAryST[refIndexSL] = 0;


   lenConAryUI =
      malloc(
           refStackST.numRefUI
         * sizeof(uint)
      ); /*allocate memory for consensus array lengths*/

   if(! lenConAryUI)
      goto memErr_fun0x_sec0x_sub02;

   for(
      refIndexSL = 0;
      refIndexSL < (slong) refStackST.numRefUI;
      ++refIndexSL
   ) lenConAryUI[refIndexSL] =
        refStackST.lenAryUI[refIndexSL];


   rmArySCPtr =
      malloc(
           refStackST.numRefUI
         * sizeof(schar *)
      ); /*allocate memory for filter variants array*/

   if(! rmArySCPtr)
      goto memErr_fun0x_sec0x_sub02;

   for(
      refIndexSL = 0;
      refIndexSL < (slong) refStackST.numRefUI;
      ++refIndexSL
   ) rmArySCPtr[refIndexSL] = 0; /*initialize*/

   for(
      refIndexSL = 0;
      refIndexSL < (slong) refStackST.numRefUI;
      ++refIndexSL
   ){ /*Loop: allocate sub arrays*/
      rmArySCPtr[refIndexSL] = 
         malloc(
              (refStackST.lenAryUI[refIndexSL] + 1)
            * sizeof(schar)
         ); /*allocate memory for sub filter arrays*/

      if(! rmArySCPtr[refIndexSL])
         goto memErr_fun0x_sec0x_sub02;
   } /*Loop: allocate sub arrays*/
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun0x Sec04:
   ^   - build consensuses
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(! *errSCPtr)
   { /*Loop: build consensus*/
      refIndexSL =
         find_strAry(
            refStackST.idAryStr,
            samSTPtr->refIdStr,
            refStackST.numRefUI
         ); /*find reference id*/

      if(refIndexSL < 0)
         goto nextRead_fun0x_sec0y_sub0z;

      *errSCPtr =
         addRead_tbCon(
            samSTPtr,
            &(conHeapAryST[refIndexSL]),
            &lenConAryUI[refIndexSL],
            conSetSTPtr
         ); /*build consensus for reference*/

      nextRead_fun0x_sec0y_sub0z:;

      *errSCPtr =
         get_samEntry(
            samSTPtr,
            buffStrPtr,
            lenBuffULPtr,
            samFILE
         );
   } /*Loop: build consensus*/

   if(*errSCPtr == def_memErr_samEntry)
     goto memErr_fun0x_sec0x_sub02;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun0x Sec06:
   ^   - detect low scoring variants
   ^   o fun0X sec06 sub01:
   ^     - select consensus + start loop
   ^   o fun0X sec06 sub02:
   ^     - find consensus position support + start loop
   ^   o fun0X sec06 sub03:
   ^     - find most supported base
   ^   o fun0X sec06 sub04:
   ^     - find bases to remove
   ^   o fun0X sec06 sub05:
   ^     - free consensus found variants for
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun0X Sec06 Sub01:
   *   - select consensus + start loop
   \*****************************************************/

   for(
      refIndexSL = 0;
      refIndexSL < refStackST.numRefUI;
      ++refIndexSL
   ){ /*Loop: find low scoring variants*/
      conOnST = conHeapAryST[refIndexSL];

      rmArySC = rmArySCPtr[refIndexSL];

      /**************************************************\
      * Fun0X Sec06 Sub02:
      *   - find consensus position support + start loop
      \**************************************************/

      for(
         uiCon = 0;
         uiCon < lenConAryUI[refIndexSL];
         ++uiCon
      ){ /*Loop: find  number of variants*/
          rmArySC[uiCon] = 0;

          if(conOnST[uiCon].ntKeptSI < 0)
             continue;

          /*find percent probabilities*/

          aPercSupF = (float) conOnST[uiCon].numASI;
          aPercSupF /= (float) conOnST[uiCon].ntKeptSI;

          cPercSupF = (float) conOnST[uiCon].numCSI;
          cPercSupF /= (float) conOnST[uiCon].ntKeptSI;

          gPercSupF = (float) conOnST[uiCon].numGSI;
          gPercSupF /= (float) conOnST[uiCon].ntKeptSI;

          tPercSupF = (float) conOnST[uiCon].numTSI;
          tPercSupF /= (float) conOnST[uiCon].ntKeptSI;

          /**********************************************\
          * Fun0X Sec06 Sub03:
          *   - find most supported base
          \**********************************************/

          if(aPercSupF > cPercSupF)
             bestPercSupF = aPercSupF;
          else
             bestPercSupF = cPercSupF;

          if(bestPercSupF < gPercSupF)
             bestPercSupF = gPercSupF;

          if(bestPercSupF < tPercSupF)
             bestPercSupF = tPercSupF;

          /**********************************************\
          * Fun0X Sec06 Sub04:
          *   - find bases to remove
          \**********************************************/

          if(
                aPercSupF < bestPercSupF
             && aPercSupF >= minSupF
          ) rmArySC[uiCon] |= def_rmA_varFilt;

          if(
                cPercSupF < bestPercSupF
             && cPercSupF >= minSupF
          ) rmArySC[uiCon] |= def_rmC_varFilt;

          if(
                gPercSupF < bestPercSupF
             && gPercSupF >= minSupF
          ) rmArySC[uiCon] |= def_rmG_varFilt;

          if(
                tPercSupF < bestPercSupF
             && tPercSupF >= minSupF
          ) rmArySC[uiCon] |= def_rmT_varFilt;

          if(rmArySC[uiCon])
             ++totalVarUI;
      } /*Loop: find  number of variants*/

      /**************************************************\
      * Fun0X Sec06 Sub05:
      *   - free consensus found variants for
      \**************************************************/

      freeHeapAry_conNt_tbCon(
         conHeapAryST[refIndexSL],
         lenConAryUI[refIndexSL]
      );

      conHeapAryST[refIndexSL] = 0;
   } /*Loop: find low scoring variants*/

   free(conHeapAryST);
   conHeapAryST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun0x Sec07:
   ^   - filter out reads with low scoring variants
   ^   o fun0X sec07 sub01:
   ^     - find start of file, start loop, print headers
   ^   o fun0X sec07 sub02:
   ^     - find reference, check q-score entry
   ^   o fun0X sec07 sub03:
   ^     - move past start masking
   ^   o fun0X sec07 sub04:
   ^     - check cigar for ignored entries + start loop
   ^   o fun0X sec07 sub05:
   ^     - see if snp is a read to remove
   ^   o fun0X sec07 sub06:
   ^     - if pass filters; print sam entry
   ^   o fun0X sec07 sub07:
   ^     - move to next read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun0X Sec07 Sub01:
   *   - find start of file, start loop, print headers
   \*****************************************************/

   fseek(
      samFILE,
      0,
      SEEK_SET
   ); /*move to start of file*/

   *errSCPtr =
      get_samEntry(
         samSTPtr,
         buffStrPtr,
         lenBuffULPtr,
         samFILE
      ); /*already read once, so no worry of errors*/

   while(! *errSCPtr)
   { /*Loop: find reads*/
      if(samSTPtr->extraStr[0] == '@')
      { /*If: header*/
         fprintf(
            outFILE,
            "%s\n",
            samSTPtr->extraStr
         );

         goto nextRead_fun0x_sec07_sub0z;
      } /*If: header*/

      /**************************************************\
      * Fun0X Sec07 Sub02:
      *   - find reference, check q-score entry
      \**************************************************/

      refIndexSL =
         find_strAry(
            refStackST.idAryStr,
            samSTPtr->refIdStr,
            refStackST.numRefUI
         ); /*find reference id*/

      if(refIndexSL < 0)
         goto nextRead_fun0x_sec07_sub0z;

      uiRead = 0;
      rmArySC = rmArySCPtr[refIndexSL];

      /*check if have q-score entry*/

      if(
            ! samSTPtr->qStr
         || (
                 samSTPtr->qStr[0] == '*'
              && samSTPtr->qStr[1] < 32
            )
      ) qBl = 0;

      else
         qBl = 1;

      /**************************************************\
      * Fun0X Sec07 Sub03:
      *   - move past start masking
      \**************************************************/

      if(
            samSTPtr->cigTypeStr[0] == 'S'
         || samSTPtr->cigTypeStr[0] == 'H'
      ){ /*If: masked*/
         cigIndexUI = 1;
         cigValSI = samSTPtr->cigArySI[1];
      } /*If: masked*/

      else
      { /*Else: no masking*/
         cigIndexUI = 0;
         cigValSI = samSTPtr->cigArySI[0];
      } /*Else: no masking*/

      /**************************************************\
      * Fun0X Sec07 Sub04:
      *   - check cigar for ignored entries + start loop
      \**************************************************/

      uiCon = 0;
      readVarUI = 0;
      mapVarUI = 0;

      for(
         uiRead = 0;
         uiRead < samSTPtr->readLenUI;
         ++uiRead
      ){ /*Loop: detect removed variants*/
         if(
               samSTPtr->cigTypeStr[cigIndexUI] == 'S'
            || samSTPtr->cigTypeStr[cigIndexUI] == 'H'
         ) break; /*end of alignment*/

         --cigValSI;

         if(samSTPtr->cigTypeStr[cigIndexUI] == 'D')
            ++uiCon;
            
         else if(samSTPtr->cigTypeStr[cigIndexUI] =='I') ;

         else if(
               qBl
            &&    samSTPtr->qStr[uiRead]
               < conSetSTPtr->minQSI
         ) ++uiCon; /*snp to low quality*/

         /***********************************************\
         * Fun0X Sec07 Sub05:
         *   - see if snp is a read to remove
         \***********************************************/

         else if(! rmArySC[uiCon])
            ++uiCon; /*no filter variants at position*/

         else
         { /*Else: is snp*/
            switch(samSTPtr->seqStr[uiRead] & ~32)
            { /*Switch: check if filter base*/
               case 'A':
                  if(rmArySC[uiCon] & def_rmA_varFilt)
                     ++readVarUI;
                  ++mapVarUI;
                  break;

               case 'C':
                  if(rmArySC[uiCon] & def_rmC_varFilt)
                     ++readVarUI;
                  ++mapVarUI;
                  break;

               case 'G':
                  if(rmArySC[uiCon] & def_rmG_varFilt)
                     ++readVarUI;
                  ++mapVarUI;
                  break;

               case 'T':
                  if(rmArySC[uiCon] & def_rmT_varFilt)
                     ++readVarUI;
                  ++mapVarUI;
                  break;
            } /*Switch: check if filter base*/

            ++uiCon;
         } /*Else: is snp*/

         if(cigValSI == 0)
         { /*If: next cigar entry*/
            ++cigIndexUI;
            cigValSI = samSTPtr->cigArySI[cigIndexUI];
         } /*If: next cigar entry*/
      } /*Loop: detect removed variants*/
     
      percCoverF = (float) lenConAryUI[refIndexSL];
      percCoverF /= (float) samSTPtr->alnReadLenUI;

      percVarF = (float) readVarUI;
      percVarF /= (float) mapVarUI;
      percVarF *= percCoverF;

      if(percVarF > minVarF)
         goto nextRead_fun0x_sec07_sub0z; /*to many var*/

      /**************************************************\
      * Fun0X Sec07 Sub06:
      *   - if pass filters; print sam entry
      \**************************************************/

      p_samEntry(
         samSTPtr,
         buffStrPtr,
         lenBuffULPtr,
         0,
         outFILE
      ); /*print out sam entry keeping*/

      /**************************************************\
      * Fun0X Sec07 Sub07:
      *   - move to next read
      \**************************************************/

      nextRead_fun0x_sec07_sub0z:;

      *errSCPtr =
         get_samEntry(
            samSTPtr,
            buffStrPtr,
            lenBuffULPtr,
            samFILE
         );

   } /*Loop: find reads*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun0x Sec0x:
   ^   - return result + clean up
   ^   o fun0x sec0x sub01:
   ^     - no errors
   ^   o fun0x sec0x sub02:
   ^     - memory errors
   ^   o fun0x sec0x sub03:
   ^     - file errors
   ^   o fun0x sec0x sub04:
   ^     - error clean up (all errors call)
   ^   o fun0x sec0x sub05:
   ^     - return and clean up (everything calls)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun0x Sec0x Sub01:
   *   - no errors
   \*****************************************************/

   *errSCPtr = 0;
   goto ret_fun0x_sec0x_sec04;

   /*****************************************************\
   * Fun0x Sec0x Sub02:
   *   - memory errors
   \*****************************************************/

   memErr_fun0x_sec0x_sub02:;
   *errSCPtr = def_memErr_varFilt;
   goto errCleanUp_fun0x_sec0x_sub03;

   /*****************************************************\
   * Fun0x Sec0x Sub03:
   *   - file errors
   \*****************************************************/

   fileErr_fun0x_sec0x_sub03:;
   *errSCPtr = def_fileErr_varFilt;
   goto errCleanUp_fun0x_sec0x_sub03;

   /*****************************************************\
   * Fun0x Sec0x Sub04:
   *   - error clean up (all errors call)
   \*****************************************************/

   errCleanUp_fun0x_sec0x_sub03:;
   goto ret_fun0x_sec0x_sec04;

   /*****************************************************\
   * Fun0x Sec0x Sub05:
   *   - return and clean up (everything calls)
   \*****************************************************/

   ret_fun0x_sec0x_sec04:;

   freeStack_refs_samEntry(&refStackST);

   if(conHeapAryST)
   { /*If: need to free consensus arrays*/
      for(
         refIndexSL = 0;
         refIndexSL < (slong) refStackST.numRefUI;
         ++refIndexSL
      ){ /*Loop: free consensus arrays*/
         freeHeapAry_conNt_tbCon(
            conHeapAryST[refIndexSL],
            lenConAryUI[refIndexSL]
         );

         conHeapAryST[refIndexSL] = 0;
         lenConAryUI[refIndexSL] = 0;
      } /*Loop: free consensus arrays*/
   } /*If: need to free consensus arrays*/

   free(conHeapAryST);
   conHeapAryST = 0;

   if(lenConAryUI)
      free(lenConAryUI);
   lenConAryUI = 0;

   if(rmArySCPtr)
   { /*If: need to free variant filter arrays*/
      for(
         refIndexSL = 0;
         refIndexSL < (slong) refStackST.numRefUI;
         ++refIndexSL
      ){ /*Loop: free variant filter arrays*/
         free(rmArySCPtr[refIndexSL]);
         rmArySCPtr[refIndexSL] = 0;
      } /*Loop: free variant filter arrays*/
   } /*If: need to free variant filter arrays*/
} /*conFilt_varFilt*/

/*-------------------------------------------------------\
| Fun02: getProf_varFilt
|   - builds a profile for filtering
| Input:
|   - minSupF:
|     o minimum percent support to keep variant
|   - varUIPtr:
|     o pointer to unsigned int to hold number of variants
|       kept
|   - conAryST:
|     o conNt_tbCon array with the uncollapsed consensus
|       to use in variant filtering
|   - lenConUI:
|     o length (number of bases) in conAryST
| Output:
|   - Modifies:
|     o varUIPtr to have number of found variants
|   - Returns:
|     o pointer to signed char array with variants to
|       filter out
\-------------------------------------------------------*/
signed char *
getProf_varFilt(
   float minSupF,                /*min % sup to filter*/
   unsigned int *varUIPtr,       /*number variants found*/
   struct conNt_tbCon *conAryST, /*uncollapsed con*/
   unsigned int lenConUI         /*consensus length*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   o fun02 sec01:
   '     - variable declarations
   '   o fun02 sec02:
   '     - initialize and allocate memory
   '   o fun02 sec03:
   '     - find filter profile
   '   o fun02 sec04:
   '     - return results
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned int uiCon = 0;        /*position in con/ref*/

   /*for finding variant percent support*/
   float aPercSupF = 0;
   float cPercSupF = 0;
   float gPercSupF = 0;
   float tPercSupF = 0;
   float bestPercSupF = 0; /*best support*/

   schar *varHeapArySC = 0;     /*sanity for rmArySCPtr*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - initialize and allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *varUIPtr = 0;
   varHeapArySC = malloc((lenConUI + 1) * sizeof(schar));

   if(! varHeapArySC)
      goto memErr_fun02_sec0x;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - find filter profile
   ^   o fun02 sec03 sub01:
   ^     - find support for each base + start loop
   ^   o fun02 sec03 sub02:
   ^     - find most supported base
   ^   o fun02 sec03 sub03:
   ^     - check if add base to profile
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec03 Sub01:
   *   - find support for each base + start loop
   \*****************************************************/

   for(
      uiCon = 0;
      uiCon < lenConUI;
      ++uiCon
   ){ /*Loop: build variant array*/
      varHeapArySC[uiCon] = 0;

      if(conAryST[uiCon].ntKeptSI < 0)
         continue;

      /*find percent probabilities*/

      aPercSupF = (float) conAryST[uiCon].numASI;
      aPercSupF /= (float) conAryST[uiCon].ntKeptSI;

      cPercSupF = (float) conAryST[uiCon].numCSI;
      cPercSupF /= (float) conAryST[uiCon].ntKeptSI;

      gPercSupF = (float) conAryST[uiCon].numGSI;
      gPercSupF /= (float) conAryST[uiCon].ntKeptSI;

      tPercSupF = (float) conAryST[uiCon].numTSI;
      tPercSupF /= (float) conAryST[uiCon].ntKeptSI;

      /**************************************************\
      * Fun02 Sec03 Sub02:
      *   - find most supported base
      \**************************************************/

      if(aPercSupF > cPercSupF)
         bestPercSupF = aPercSupF;
      else
         bestPercSupF = cPercSupF;

      if(bestPercSupF < gPercSupF)
         bestPercSupF = gPercSupF;

      if(bestPercSupF < tPercSupF)
         bestPercSupF = tPercSupF;

      /**************************************************\
      * Fun02 Sec03 Sub03:
      *   - check if add base to profile
      \**************************************************/

      if(
            aPercSupF < bestPercSupF
         && aPercSupF >= minSupF
      ) varHeapArySC[uiCon] |= def_rmA_varFilt;

      if(
            cPercSupF < bestPercSupF
         && cPercSupF >= minSupF
      ) varHeapArySC[uiCon] |= def_rmC_varFilt;

      if(
            gPercSupF < bestPercSupF
         && gPercSupF >= minSupF
      ) varHeapArySC[uiCon] |= def_rmG_varFilt;

      if(
            tPercSupF < bestPercSupF
         && tPercSupF >= minSupF
      ) varHeapArySC[uiCon] |= def_rmT_varFilt;

      if(varHeapArySC[uiCon])
         ++(*varUIPtr);
   } /*Loop: build variant array*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^   - return results
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun02_sec0x;

   memErr_fun02_sec0x:;
      goto errCleanUp_fun02_sec0x;

   errCleanUp_fun02_sec0x:;
      if(varHeapArySC)
         free(varHeapArySC);
      varHeapArySC = 0;
      goto ret_fun02_sec0x;

   ret_fun02_sec0x:;
   return varHeapArySC;
} /*getProf_varFilt*/

/*-------------------------------------------------------\
| Fun03: filtReads_varFilt
|   - filters reads by variants in profile
| Input:
|   - minPercVarF:
|     o min percent of variants needed to keep a read
|     o adjusted by read length (% variants * % coverage)
|   - clustSI:
|     o cluster working on
|   - varArySC:
|     o c-string with variants to remove reads by
|   - lenVarUI:
|     o length of varArySC (used to find % coverage)
|   - indexSTPtr:
|     o pointer to index_clustST struct with read indexs
|       and assigned clusters
|   - samSTPtr:
|     o pointer to samEntry struct to get reads with
|   - buffStr:
|     o c-string to use in reading sam file
|     o needs to be as long as longest entry in sam file
|       (scoreReads_clustST does this for you)
|   - samFILE:
|     o FILE pointer to sam file to scan
\-------------------------------------------------------*/
signed char *
filtReads_varFilt(
   float minPercVarF,            /*min % sup to filter*/
   unsigned int clustSI,         /*cluster working on*/
   signed char *varArySC,        /*array of variants*/
   unsigned int lenVarUI,        /*length of varArySC*/
   struct index_clustST *indexSTPtr, /*read index's*/
   struct samEntry *samSTPtr,    /*for geting sam entrys*/
   signed char *buffStr,         /*for reading sam file*/
   void *samFILE                 /*sam file with reads*/
){
} /*filtReads_varFilt*/
   
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   o fun03 sec01:
   '     - variable declarations
   '   o fun03 sec02:
   '     - initialize and allocate memory
   '   o fun03 sec03:
   '     - find filter profile
   '   o fun03 sec04:
   '     - return results
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned int uiCon = 0;        /*position in con/ref*/

   /*for finding variant percent support*/
   float aPercSupF = 0;
   float cPercSupF = 0;
   float gPercSupF = 0;
   float tPercSupF = 0;
   float bestPercSupF = 0; /*best support*/

   schar *varHeapArySC = 0;     /*sanity for rmArySCPtr*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - initialize and allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *varUI = 0;
   varHeapArySC = malloc((lenConUI + 1) * sizeof(schar));

   if(! varHeapArySC)
      goto memErr_fun03_sec0x;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - find filter profile
   ^   o fun03 sec03 sub01:
   ^     - find support for each base + start loop
   ^   o fun03 sec03 sub02:
   ^     - find most supported base
   ^   o fun03 sec03 sub03:
   ^     - check if add base to profile
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - find support for each base + start loop
   \*****************************************************/

   for(
      uiCon = 0;
      uiCon < lenConUI;
      ++uiCon
   ){ /*Loop: build variant array*/
      varHeapArySC[uiCon] = 0;

      if(conAryST[uiCon].ntKeptSI < 0)
         continue;

      /*find percent probabilities*/

      aPercSupF = (float) conAryST[uiCon].numASI;
      aPercSupF /= (float) conAryST[uiCon].ntKeptSI;

      cPercSupF = (float) conAryST[uiCon].numCSI;
      cPercSupF /= (float) conAryST[uiCon].ntKeptSI;

      gPercSupF = (float) conAryST[uiCon].numGSI;
      gPercSupF /= (float) conAryST[uiCon].ntKeptSI;

      tPercSupF = (float) conAryST[uiCon].numTSI;
      tPercSupF /= (float) conAryST[uiCon].ntKeptSI;

      /**************************************************\
      * Fun03 Sec03 Sub02:
      *   - find most supported base
      \**************************************************/

      if(aPercSupF > cPercSupF)
         bestPercSupF = aPercSupF;
      else
         bestPercSupF = cPercSupF;

      if(bestPercSupF < gPercSupF)
         bestPercSupF = gPercSupF;

      if(bestPercSupF < tPercSupF)
         bestPercSupF = tPercSupF;

      /**************************************************\
      * Fun03 Sec03 Sub03:
      *   - check if add base to profile
      \**************************************************/

      if(
            aPercSupF < bestPercSupF
         && aPercSupF >= minSupF
      ) varHeapArySC[uiCon] |= def_rmA_varFilt;

      if(
            cPercSupF < bestPercSupF
         && cPercSupF >= minSupF
      ) varHeapArySC[uiCon] |= def_rmC_varFilt;

      if(
            gPercSupF < bestPercSupF
         && gPercSupF >= minSupF
      ) varHeapArySC[uiCon] |= def_rmG_varFilt;

      if(
            tPercSupF < bestPercSupF
         && tPercSupF >= minSupF
      ) varHeapArySC[uiCon] |= def_rmT_varFilt;

      if(varHeapArySC[uiCon])
         ++(*varUI);
   } /*Loop: build variant array*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - return results
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun03_sec0x;

   memErr_fun03_sec0x:;
      goto errCleanUp_fun03_sec0x;

   errCleanUp_fun03_sec0x:;
      if(varHeapArySC)
         free(varHeapArySC);
      varHeapArySC = 0;
      goto ret_fun03_sec0x;

   ret_fun03_sec0x:;
   return varHeapArySC;
} /*getProf_varFilt*/


/*-------------------------------------------------------\
\-------------------------------------------------------*/
int
main(
   int numArgsSI,
   char *argAryStr[]
){
   /*schar *samFileStr = (schar *) "/home/reason/files/testing/mixedInfect/SRR13951165-map.sam";*/
   schar *samFileStr = (schar *) argAryStr[numArgsSI - 1];

   schar errSC = 0;

   float minPercSupF = 0.1; /*discard at 10% support*/
   float minVarF = 0.2;     /*min % var to discard read*/
   struct samEntry samStackST;
   struct set_tbCon conSetStackST;

   schar *headHeapStr = 0;
   ulong sizeHeadUL = 0;

   schar *buffHeapStr = 0;
   ulong lenBuffUL = 0;

   FILE *samFILE = fopen((char *) samFileStr, "r");
   FILE *outFILE = stdout;

   init_set_tbCon(&conSetStackST);

   init_samEntry(&samStackST);
   errSC = setup_samEntry(&samStackST);

   if(errSC)
      goto err_main_sec0x;

   conFilt_varFilt(
      minPercSupF,
      minVarF,
      &conSetStackST,
      &samStackST,
      &headHeapStr,
      &sizeHeadUL,
      &buffHeapStr,
      &lenBuffUL,
      samFILE,
      outFILE,
      &errSC
   );

   errSC = 0;
   goto cleanUp_main_sec0x;

   err_main_sec0x:;
   errSC = 1;
   goto cleanUp_main_sec0x;


   cleanUp_main_sec0x:;

   if(headHeapStr)
      free(headHeapStr);
   headHeapStr = 0;
   sizeHeadUL = 0;

   if(buffHeapStr)
      free(buffHeapStr);
   buffHeapStr = 0;
   lenBuffUL = 0;

   freeStack_samEntry(&samStackST);
   freeStack_set_tbCon(&conSetStackST);

   if(
         samFILE
      && samFILE != stdin
      && samFILE != stdout
   ) fclose(samFILE);

   samFILE = 0;

   if(
         outFILE
      && outFILE != stdin
      && outFILE != stdout
   ) fclose(outFILE);

   outFILE = 0;

   return errSC;
} /*main*/
