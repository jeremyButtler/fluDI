/*########################################################
# Name: trimSam
#   - Holds functions for trimming softmasked regoins at
#     the start and end of sequences in an sam file
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - Included libraries and defined variables
'   o fun01 trimSeq_trimSam:
'     o Trim soft mask regions off end of sam entry
'   o fun02: trimCoords_trimSam
'     - Trim an sam file entry by coordinates
'   o fun03 trimReads_trimSam:
'     o Trims soft mask regions for all reads with a
'       sequence in a sam file
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Included libraries and defined variables
|   o header Sec01:
|     - Included libraries
|   o header Sec02:
\-------------------------------------------------------*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec01:
^   - Included libraries
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "trimSam.h"

#include <stdio.h>

#include "samEntry.h"
#include "ulCp.h"

/*Has no .c files*/
#include "dataTypeShortHand.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden depenencies
!   o .c #include "../generalLib/base10StrToNum.h"
!   o .c #include "../generalLib/numToStr.h"
!   o .h #include "../generalLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec02:
^   - Defined varialbes
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*This .c version is slightly slower, but it works*/
/*This will work so long as the integer is beneath
`   128 bits, which is non-realistic value.
*/
#define def_uiToUC_trimSam \
 ((sizeof(unsigned int) >>1) - (sizeof(unsigned int) >>3))

/*Logic:
`   - sizeof(unsigned int) >> 1:
`     o This gets the number of bytes to shift to get
`       an 32 bit or 16 bit number to one bytes. Basically
`       it is divide by 2
`   - (sizeof(unsigned int) >> 3)
`     o This adds in 1 if I have eight bytes. It fails
`       for ints and shorts
`   - (sizeof(int) >> 1) - (sizeof(int) >> 3)
`     o Gives me the number of bits I need to shift to
`       get the number of characters in an int
`   - (sizeof(int) >> 3) << (sizeof(int) >> 4)
`     o This is not used, but would correct for a 128 bit
`       number. Each time I increase the bit size I need
`       to add a correction in
*/

/*-------------------------------------------------------\
| Fun01: trimSeq_trimSam
| Use:
|  - Trims off the soft masked regions of a sam entry
| Input:
|  - samVoidST:
|    o samEntry struct with sam entry to trim. Named
|      void so I can cast it as a samEntry to samST
| Output:
|  - Returns:
|    o 0 if suceeded
|    o def_header_trimSam if header (invalid and ignored)
|    o def_noMap_trimSam if unmapped read (no reference)
|    o def_noSeq_trimSam if no sequence line
|  - Modifies:
|    o Trims cigar, sequence, & q-score entries in samST.
\-------------------------------------------------------*/
signed char
trimSeq_trimSam(
    struct samEntry *samST   /*has sam line to trim*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC: trimSeq_trimSam
   '  - Trims soft masked regions at start & end of sam
   '    entry
   '  o fun01 sec01:
   '    - Variable declerations
   '  o fun01 sec02:
   '    - Check if I have something to trim
   '  o fun01 sec03:
   '    - Finding trimming postions and trim the cigar
   '  o fun01 sec04:
   '    - Trim the sequence and Q-score entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec01:
    ^  - Variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    schar errSC = 0;

    /*Number of bases soft masked at start or end*/
    uchar qUC = 0;

    sint startPosSI = 0;
    sint endTrimSI = 0;
    uint uiCig = 0;
    uint uiSeq = 0;

    schar *seqCpStr = 0;
    schar *seqDupStr = 0;
    ulong *seqCpUL = 0;
    ulong *seqDupUL = 0;

    schar *qCpStr = 0;
    schar *qDupStr = 0;
    ulong *qCpUL = 0;
    ulong *qDupUL = 0;
        
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec02:
    ^   - Check if I have something to trim
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if is an header is an header entry*/
    if(samST->extraStr[0] == '@')
       goto header_fun01_sec05;

    /*Check if is an unmapped read*/
    if(samST->flagUS & 4)
       goto noMap_fun01_sec05;

    /*Check if is an unmapped read*/
    if(samST->cigTypeStr[0] == '*')
       goto noMap_fun01_sec05;

    /*Check if has a sequence to trim*/
    if(samST->seqStr[0] == '*')
       goto noSeq_fun01_sec05;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec03:
    ^   - Finding trimming postions and trim the cigar
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Position my self at the start of the sequence*/
    startPosSI = -((int) samST->cigTypeStr[0] == 'S');
    uiCig -= startPosSI; /*Marks if skipping the start*/

    endTrimSI =
       -((int) samST->cigTypeStr[samST->lenCigUI-1]=='S');

    endTrimSI &= samST->cigArySI[samST->lenCigUI - 1];

    samST->lenCigUI += startPosSI; /*+0 or + -1*/
    samST->lenCigUI -= (endTrimSI > 0);

    startPosSI &= samST->cigArySI[0];
    samST->readLenUI -= (startPosSI + endTrimSI);

    cpLen_ulCp(
       samST->cigTypeStr,
       (samST->cigTypeStr + uiCig),
       samST->lenCigUI
    ); /*Trim the cigar types array*/

    cpLen_ulCp(
       (schar *) samST->cigArySI,
       (schar *) (samST->cigArySI + uiCig),
       (samST->lenCigUI << def_uiToUC_trimSam)
    ); /*Trim the cigar value (bases per type) array*/

    samST->numMaskUI = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec04:
    ^   - Trim the sequence and Q-score entry
    ^   o fun01 sec04 sub01:
    ^     - Set up the pionters for trimming (copying)
    ^   o fun01 sec04 sub02:
    ^     - Adjust Q-scores for removing bases from start
    ^   o fun01 sec04 sub03:
    ^     - Copy kept parts of Q-score & sequence entries
    ^   o fun01 sec04 sub04:
    ^     - Adjust for the Q-scores removed at the end
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /****************************************************\
    * Fun01 Sec04 Sub01:
    *   - Set up the pionters for trimming (copying)
    \****************************************************/

    /*The longs are for 8 byte copying, while the c-string
    ' pionters are for the end (when sequence was not a
    ' multiple of 8
    */

    seqCpUL = (ulong *) (samST->seqStr + startPosSI);
    seqDupUL = (ulong *) samST->seqStr;
    seqCpStr = samST->seqStr + startPosSI;
    seqDupStr = samST->seqStr;

    qCpUL = (ulong *) (samST->qStr + startPosSI);
    qDupUL = (ulong *) samST->qStr;
    qCpStr = samST->qStr + startPosSI;
    qDupStr = samST->qStr;

    /****************************************************\
    * Fun01 Sec04 Sub02:
    *   - Adjust Q-scores for removing bases from start
    \****************************************************/

    for(
       uiSeq = 0;
       uiSeq < (uint) startPosSI;
      ++uiSeq
    ){ /*Loop: adjust for the timmed starting Q-scores*/
       qUC =
          (uchar) samST->qStr[uiSeq] - def_adjQ_samEntry;

       --(samST->qHistUI[qUC]);
       samST->sumQUL -= qUC;
    } /*Loop: adjust for the timmed starting Q-scores*/

    /****************************************************\
    * Fun01 Sec04 Sub03:
    *   - Copy kept parts of Q-score & sequence entries
    \****************************************************/

    /*I am not using my copy functions here to leverage dual
    ` dual accumulators, which the cpu should hanlde.
    ` easier.
    */

    /*Copy With longs*/
    for(
       uiSeq = 0;
       uiSeq < (samST->readLenUI >> def_shiftULBy_ulCp);
       ++uiSeq
    ) { /*Loop: Copy sequences and q-score entries*/
       seqDupUL[uiSeq] = seqCpUL[uiSeq];
       qDupUL[uiSeq] = qCpUL[uiSeq];
    } /*Loop: Copy sequences and q-score entries*/

    /*Finish coping the kept region using characters*/
    for(
       uiSeq =
           samST->readLenUI
         - (samST->readLenUI & def_modUL_ulCp);
       uiSeq < samST->readLenUI;
       ++uiSeq
    ) { /*Loop: Copy sequences and q-score entries*/
       seqDupStr[uiSeq] = seqCpStr[uiSeq];
       qDupStr[uiSeq] = qCpStr[uiSeq];
    } /*Loop: Copy sequences and q-score entries*/
   
    /****************************************************\
    * Fun01 Sec04 Sub04:
    *   - Adjust for the Q-scores removed at the end
    \****************************************************/

    /*Adjust for having removed q-scores at the end*/
    while(qCpStr[uiSeq] != '\0')
    { /*Loop: adjust for the timmed ending Q-scores*/
       qUC =
          (uchar) samST->qStr[uiSeq] - def_adjQ_samEntry;

       --(samST->qHistUI[qUC]);
       samST->sumQUL -= qUC;
       ++uiSeq;
    } /*Loop: adjust for the timmed ending Q-scores*/

    /*Make sure they are c-strings*/
    seqDupStr[samST->readLenUI] = '\0';
    qDupStr[samST->readLenUI] = '\0';

    samST->meanQF=(float) samST->sumQUL /samST->readLenUI;
    qhistToMed_samEntry(samST);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec05:
    ^   - clean up
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    errSC = 0;
    goto cleanUp_fun01_sec05;

    header_fun01_sec05:;
    errSC = def_header_trimSam;
    goto cleanUp_fun01_sec05;

    noMap_fun01_sec05:;
    errSC = def_noMap_trimSam;
    goto cleanUp_fun01_sec05;

    noSeq_fun01_sec05:;
    errSC = def_noSeq_trimSam;
    goto cleanUp_fun01_sec05;

    cleanUp_fun01_sec05:;
    return errSC;
} /*trimSeq_trimSam*/

/*-------------------------------------------------------\
| Fun02: trimCoords_trimSam
|   - Trim an sam file entry by coordinates
| Input:
|   - samSTPtr:
|     o Pointer to an sam entry structure with an read to
|       trim
|   - startSI:
|     o Singed integer with the starting coordinate
|   - endSI:
|     o Singed integer with the ending coordinate
|     o < 1 is treated as do not trim
| Output:
|   - Modifies:
|     o samSTPtr to be trimmed
|   - Returns:
|     o 0 for no errors
|     o def_rangeErr_trimSam for coordinates out of range
|     o def_header_trimSam if header (invalid and ignored)
|     o def_noMap_trimSam if unmapped read (no reference)
|     o def_noSeq_trimSam if no sequence line
\-------------------------------------------------------*/
signed char
trimCoords_trimSam(
   struct samEntry *samSTPtr,
   signed int startSI,
   signed int endSI
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - Trim an sam file entry by coordinates
   '   o fun02 sec01:
   '     - variable declerations
   '   o fun02 sec02:
   '     - check ranges and set up variables
   '   o fun02 sec03:
   '     - find the starting position
   '   o fun02 sec04:
   '     - find the ending position
   '   o fun02 sec05:
   '     - set up the lengths and ending coordinates
   '   o fun02 sec06:
   '     - trim the sequence, Q-score, and cigar entries
   '   o fun02 sec07:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;

   uint seqStartUI = 0;
   sint firstCigSI = 0;

   sint siCig = 0;
   sint cigBaseOnSI = 0;
   sint refPosSI = 0;
   sint seqPosSI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - check ranges and set up variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if is an header is an header entry*/
    if(samSTPtr->extraStr[0] == '@')
       goto header_fun02_sec07;

    /*Check if is an unmapped read*/
    if(samSTPtr->flagUS & 4)
       goto noMap_fun02_sec07;

    if(samSTPtr->cigTypeStr[0] == '*')
       goto noMap_fun02_sec07;

    /*Check if has a sequence to trim*/
    if(samSTPtr->seqStr[0] == '*')
       goto noSeq_fun02_sec07;

   if(samSTPtr->refStartUI > (uint) endSI)
      goto outOfRange_fun02_sec07;

   if(samSTPtr->refEndUI < (uint) startSI)
      goto outOfRange_fun02_sec07;

   cigBaseOnSI = samSTPtr->cigArySI[0];
   refPosSI = (sint) samSTPtr->refStartUI;

   /*Setting to zero becuase it is easier to recount*/
   samSTPtr->numMaskUI = 0;
   samSTPtr->numInsUI = 0;
   samSTPtr->numDelUI = 0;
   samSTPtr->numSnpUI = 0;
   samSTPtr->numMatchUI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - find the starting position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   findRefPos_samEntry(
      samSTPtr,
      &siCig,
      &cigBaseOnSI,
      startSI,
      &refPosSI,
      &seqPosSI
   );

   samSTPtr->cigArySI[siCig] = cigBaseOnSI;
   samSTPtr->refStartUI = (uint) refPosSI;

   seqStartUI = (uint) seqPosSI;
   firstCigSI = siCig;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^   - find the ending position
   ^   o fun02 sec04 sub01:
   ^     - start loop and handle soft masking entires
   ^   o fun02 sec04 sub02:
   ^     - handle insertion entries
   ^   o fun02 sec04 sub03:
   ^     - handle deletion entries
   ^   o fun02 sec04 sub04:
   ^     - handle snp entries
   ^   o fun02 sec04 sub05:
   ^     - move to next cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec04 Sub01:
   *   - start loop and handle soft masking entires
   \*****************************************************/

   while(refPosSI < endSI)
   { /*Loop: till I am on the target base*/
      switch((samSTPtr)->cigTypeStr[(siCig)])
      { /*Switch: check what the next entry is*/
         case 'S':
         /*Case: Softmasking*/
            seqStartUI += cigBaseOnSI;
            samSTPtr->numMaskUI +=  cigBaseOnSI;

            ++siCig;
            cigBaseOnSI = 0;
            break;
         /*Case: Softmasking*/

         /***********************************************\
         * Fun02 Sec04 Sub02:
         *   - handle insertion entries
         \***********************************************/

         case 'I':
         /*Case: insertions*/
            seqStartUI += cigBaseOnSI;
            samSTPtr->numInsUI +=  cigBaseOnSI;

            ++siCig;
            cigBaseOnSI = 0;
            break;
         /*Case: insertions*/

         /***********************************************\
         * Fun02 Sec04 Sub03:
         *   - handle deletion entries
         \***********************************************/

         case 'D':
         /*Case: Deletion*/
            refPosSI += cigBaseOnSI;
            samSTPtr->numDelUI += cigBaseOnSI;

            if(refPosSI <= endSI)
            { /*If: I have not found target position*/
               ++siCig;
               cigBaseOnSI = 0;
            } /*If: I have not found target position*/

            else 
            { /*Else: I overshot the target*/
               /*Find how many bases overshot by*/
               cigBaseOnSI = (refPosSI - endSI);

               /*Make corrections for overshooting*/
               refPosSI -= cigBaseOnSI;
               samSTPtr->numDelUI -= cigBaseOnSI;
            } /*Else: I overshot the target*/

            break;
         /*Case: Deletion*/

         /***********************************************\
         * Fun02 Sec04 Sub03:
         *   - handle match entries
         \***********************************************/

         case 'M':
         case '=':
         /*Case: match (M or =)*/
            refPosSI += cigBaseOnSI;
            seqPosSI += cigBaseOnSI;
            samSTPtr->numMatchUI += cigBaseOnSI;

            if(refPosSI <= endSI)
            { /*If: I have not found target position*/
               cigBaseOnSI = 0;
               ++siCig;
            } /*If: I have not found target position*/

            else
            { /*Else: I overshot the target*/
               /*Find how many bases overshot by*/
               cigBaseOnSI = (int) (refPosSI - endSI);

               /*Make corrections for overshooting*/
               refPosSI -= cigBaseOnSI;
               seqPosSI -= cigBaseOnSI;
               samSTPtr->numMatchUI -= cigBaseOnSI;
            } /*Else: I overshot the target*/

            break;
         /*Case: match (M or =)*/

         /***********************************************\
         * Fun02 Sec04 Sub04:
         *   - handle snp entries
         \***********************************************/

         case 'X':
         /*Case: SNP (X)*/
            refPosSI += cigBaseOnSI;
            seqPosSI += cigBaseOnSI;
            samSTPtr->numSnpUI += cigBaseOnSI;

            if(refPosSI <= endSI)
            { /*If: I have not found target position*/
               cigBaseOnSI = 0;
               ++siCig;
            } /*If: I have not found target position*/

            else
            { /*Else: I overshot the target*/
               /*Find how many bases overshot by*/
               cigBaseOnSI = (int) (refPosSI - endSI);

               /*Make corrections for overshooting*/
               refPosSI -= cigBaseOnSI;
               seqPosSI -= cigBaseOnSI;

               samSTPtr->numSnpUI = cigBaseOnSI;
            } /*Else: I overshot the target*/

            break;
         /*Case: SNP (X)*/
      } /*Switch: check what the next entry is*/

      /**************************************************\
      * Fun02 Sec04 Sub05:
      *   - move to next cigar entry
      \**************************************************/

      if(siCig >= (sint) samSTPtr->lenCigUI)
         break; /*End of the sequence*/

      /*This case will be true most of the time, unless
      `   the start has already been found
      */
      if(cigBaseOnSI == 0)
         cigBaseOnSI = samSTPtr->cigArySI[siCig];
   } /*Loop: till I am on the target base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec05:
   ^   - set up the lengths and ending coordinates
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   samSTPtr->cigArySI[siCig] = cigBaseOnSI;
   samSTPtr->lenCigUI = siCig - firstCigSI;

   samSTPtr->readLenUI = (uint) (seqPosSI - seqStartUI);
   samSTPtr->refEndUI = refPosSI;
   samSTPtr->alnReadLenUI = samSTPtr->refEndUI;
   samSTPtr->alnReadLenUI -= samSTPtr->refStartUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec06:
   ^   - trim the sequence, Q-score, and cigar entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      &(samSTPtr->seqStr[seqStartUI]),
      samSTPtr->readLenUI
   );

   if(samSTPtr->qStr[0] !='*' && samSTPtr->qStr[1] !='\0')
   { /*If: There is no Q-score entry*/
      cpQEntry_samEntry(
         samSTPtr,
         &(samSTPtr->qStr[seqStartUI]),
         1
      ); /*Uses samSTPtr->readLenUI to get length*/
         /*Also finds median and mean Q-scores*/
   } /*If: There is no Q-score entry*/

   cpLen_ulCp(
      samSTPtr->cigTypeStr,
      &(samSTPtr->cigTypeStr[firstCigSI]),
      samSTPtr->lenCigUI
   );

   cpLen_ulCp(
      (schar *) samSTPtr->cigArySI,
      (schar *) &(samSTPtr->cigArySI[firstCigSI]),
      (samSTPtr->lenCigUI << 2)
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec07:
   ^   - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto cleanUp_fun02_sec07;

   outOfRange_fun02_sec07:;
   errSC = def_rangeErr_trimSam;
   goto cleanUp_fun02_sec07;

   header_fun02_sec07:;
   errSC = def_header_trimSam;
   goto cleanUp_fun02_sec07;

   noMap_fun02_sec07:;
   errSC = def_noMap_trimSam;
   goto cleanUp_fun02_sec07;

   noSeq_fun02_sec07:;
   errSC = def_noSeq_trimSam;
   goto cleanUp_fun02_sec07;

   cleanUp_fun02_sec07:;
   return errSC;
} /*trimCoords_trimSam*/

/*-------------------------------------------------------\
| Fun03: trimReads_trimSam
| Use:
|  - Goes though sam file and calls trimSeq_trimSam for each
|    entry
| Input:
|  - samFILE:
|    o Pointer to FILE. Sam file with entries to trim
|  - outFILE:
|    o File to write the trimmed sam file entries to
|  - keepUmapBl:
|    o Also print out entries with unmapped reads
| Output:
|  - Returns:
|    o 0 if suceeded
|    o def_memErr_trimSam for memory errors
|  - Prints:
|    o trimmed sam entries with sequences to outFILE, but
|      ignores sam entries without sequences
\-------------------------------------------------------*/
signed char
trimReads_trimSam(
    void *samFILE,           /*sam file to trim*/
    void *outFILE,           /*file to store output*/
    signed char keepUmapBl   /*1: keep unmapped reads*/
){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ' Fun03 TOC: trimReads_trimSam
   '  - Goes though sam file and calls trimSeq_trimSam for
   '    each entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    
    schar *buffHeapStr = 0;
    ulong lenBuffUL = 0;
    schar errSC = 0;  /*Tells if memory error*/
    struct samEntry samStackST;
    
    init_samEntry(&samStackST);
    errSC = setup_samEntry(&samStackST);
    
    if(errSC)
       goto memErr_fun03;
   
    errSC =
       get_samEntry(
          &samStackST,
          &buffHeapStr,
          &lenBuffUL,
          (FILE *) samFILE
      );
    
    while(! errSC)
    { /*While there are lines in sam file to convert*/
        if(*(samStackST.extraStr) == '@')
        { /*If was a header*/
           errSC =
              p_samEntry(
                 &samStackST,
                 &buffHeapStr,
                 &lenBuffUL,
                 0,          /*want full line output*/
                 (FILE *) outFILE
              );

           if(errSC)
              goto memErr_fun03;
            
           errSC =
              get_samEntry(
                 &samStackST,
                 &buffHeapStr,
                 &lenBuffUL,
                 (FILE *) samFILE
              );
           
           if(errSC)
              goto memErr_fun03;

           continue; /*header line, move to next line*/
        } /*If was a header*/
        
        /*Convert & print out sam file entry*/
        errSC = trimSeq_trimSam(&samStackST);

        if(errSC)
        { /*If: had problematic input*/
           if(errSC & def_noMap_trimSam && keepUmapBl)
           { /*If: printing umapped reads*/
              p_samEntry(
                  &samStackST,
                  &buffHeapStr,
                  &lenBuffUL,
                  0,          /*want full line output*/
                  (FILE *) outFILE
              );
           } /*If: printing umapped reads*/

           /*remainder is no sequence or header*/
        } /*If: had problematic input*/

        else
        { /*Else: read was trimmed*/
           p_samEntry(
               &samStackST,
               &buffHeapStr,
               &lenBuffUL,
               0,          /*want full line output*/
               (FILE *) outFILE
           ); /*print the trimed entry*/
        } /*Else: read was trimmed*/
        
        errSC =
           get_samEntry(
              &samStackST,
              &buffHeapStr,
              &lenBuffUL,
              (FILE *) samFILE
           );
    } /*While there are lines in sam file to convert*/
    
    if(errSC == def_memErr_samEntry)
       goto memErr_fun03;
     
    errSC = 0;
    goto cleanUp_fun03;

    memErr_fun03:;
    errSC = def_memErr_trimSam;
    goto cleanUp_fun03;

    cleanUp_fun03:;

    freeStack_samEntry(&samStackST);

    free(buffHeapStr);
    buffHeapStr = 0;

    return errSC;
} /*trimReads_trimSam*/

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
