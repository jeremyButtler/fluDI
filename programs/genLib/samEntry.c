/*########################################################
# Name: samEntryStruct
#  - Holds structer to hold a sam file entry. This also
#    includes the functions needed to support this
#    structure.
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - Included libraries
'   o .h st01 samEntry:
'     - Holds a single samfile entry
'   o fun01 blank_samEntry:
'     - Sets all non-alloacted variables in samEntryST to 0
'   o fun02 init_samEntry:
'     - Initalize a samEntry struct to 0's
'   o fun03: setup_samEntry
'     - allocates memory for a samEntry structure (call
'       after init_samEntry (fun02))
'   o fun04 freeStack_samEntry:
'     - Frees heap allocations in a stack allocated
'       samEntry struct
'   o fun05 freeHeap_samEntry:
'     - Frees a samEntry structer (and sets to null)
'   o fun06: mk_samEntry
'     - Makes an heap allocated samEntry structure
'   o fun07: qhistToMed_samEntry
'     - Gets the median q-score for an histogram of
'       q-scores in a samStruct
'   o fun08: findQScores_samEntry
'      - Gets the median and mean q-scores from a samEntry
'        Structure.
'   o fun09: cpQEntry_samEntry
'     - Copies q-scores from a string into a samEntry
'       structure
'   o fun10: get_samEntry
'     - Reads in a single line from a sam file
'   o fun11: findRefPos_samEntry
'     - Find an reference coordinate in an sequence in
'       an sam entry structure
'   o fun12: p_samEntry
'     - Prints the sam file entry to a file. This does not
'       print any extra stats that were found.
'   o fun13: pfq_samEntry
'     - Prints the sam entry as a fastq entry to a fastq
'       file
'   o fun14: pfa_samEntry
'     - Prints the sam entry as a fasta entry to a fasta
'       file
'   o fun15: pstats_samEntry
'     - Prints out the stats in a samEntry struct to a file
'   o fun16: revCmp_samEntry
'     - reverse complements a sam file sequence entry
'   o .h note01:
'      - Notes about the sam file format from the sam file
'        pdf
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

#include "samEntry.h"

#include <stdio.h>

#include "base10str.h"
#include "numToStr.h"
#include "ulCp.h"

/*These have no .c files*/
#include "dataTypeShortHand.h"
#include "ntTo5Bit.h"
   /*look up table to see if have anonymous bases*/

/*-------------------------------------------------------\
| Fun01: blank_samEntry
|   - Sets all values to 0, or for c-strings to '\0'
| Input:
|   - samSTPtr:
|     o Pointer to samEntry structure to blank
| Output:
|   - Modifies:
|     o Sets every variable to 0
|     o the first character in every c-string to be '\0'
\-------------------------------------------------------*/
void
blank_samEntry(
   struct samEntry *samSTPtr
){
    uint uiIter = 0;

    samSTPtr->qryIdStr[0] = '\0'; /*query id/name*/
    samSTPtr->lenQryIdUC = 0;     /*Length of query id*/

    samSTPtr->refIdStr[0] = '\0'; /*reference id/name*/
    samSTPtr->lenRefIdUC = 0;     /*Length; reference id*/

    samSTPtr->lenCigUI = 0;       /*Length; cigar entry*/

    for(
       uiIter=0;
       uiIter < samSTPtr->lenCigBuffUI;
       ++uiIter
    ){ /*Loop: Clear all cigar entries*/
       samSTPtr->cigTypeStr[uiIter] = '\0';
       samSTPtr->cigArySI[uiIter] = 0;
    } /*Loop: Clear all cigar entries*/

    samSTPtr->rNextStr[0] = '\0'; /*rNext entry*/
    samSTPtr->lenRNextUC = 0;     /*Length of rNext*/

    if(samSTPtr->seqStr)
       samSTPtr->seqStr[0] = '\0';

    if(samSTPtr->qStr)
       samSTPtr->qStr[0] = '\0';

    if(samSTPtr->extraStr)
       samSTPtr->extraStr[0] = '\0';

    samSTPtr->lenExtraUI = 0;

    /*Flags/single numeric values in the sam entry*/
    samSTPtr->mapqUC = 0;         /*mapping quality*/
    samSTPtr->flagUS = 0;         /*flag in samEntryCStr*/
    samSTPtr->pNextSI = 0;         /*PNEXT postion*/
    samSTPtr->tLenSI = 0;          /*TLEN*/

    /*Stats*/

    /*The q-scores excluded the softmasked regions*/
    samSTPtr->medianQF = 0;       /*median read q-score*/
    samSTPtr->meanQF = 0;         /*mean read q-score*/

    samSTPtr->refStartUI = 0;     /*First reference base*/
    samSTPtr->refEndUI = 0;       /*Last reference base*/
    samSTPtr->readLenUI = 0;      /*Holds read length*/
    samSTPtr->alnReadLenUI = 0;   /*Aligned length*/

    samSTPtr->numMatchUI = 0;     /*number of matches*/
    samSTPtr->numSnpUI = 0;       /*number of mismatches*/
    samSTPtr->numInsUI = 0;       /*number of insertions*/
    samSTPtr->numDelUI = 0;       /*number of deletions*/
    samSTPtr->numMaskUI = 0;      /*number soft masked*/

    /*These variables are used in finding the q-scores*/
    for(
       uiIter = 0;
       uiIter < def_maxQ_samEntry;
       ++uiIter
    ) samSTPtr->qHistUI[uiIter] = 0;

    samSTPtr->sumQUL = 0;
} /*blank_samEntry*/

/*-------------------------------------------------------\
| Fun02: init_samEntry
|   - Initializes a samEntry structure for use. This 
|     function should only ever be called once per
|     structure or after freeStack_samEntry has been used.
|     Use blank_samEntry to reset the values in a samEntry
|     struct after initialization.
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to initialize
| Output:
|   - Modifies:
|     o all variables in samEntry to be 0 or null
|   - Returns:
|     o 0 for success
|     o 1 for memory error
\-------------------------------------------------------*/
void
init_samEntry(
   struct samEntry *samSTPtr
){
    (samSTPtr)->seqStr = 0;
    (samSTPtr)->lenSeqBuffUI = 0;

    (samSTPtr)->qStr = 0;
    (samSTPtr)->lenQBuffUI = 0;

    (samSTPtr)->cigTypeStr = 0;
    (samSTPtr)->lenCigBuffUI = 0;
    (samSTPtr)->cigArySI = 0;

    (samSTPtr)->extraStr = 0;
    (samSTPtr)->lenExtraBuffUI = 0;

    blank_samEntry(samSTPtr);
} /*init_samEntry*/

/*-------------------------------------------------------\
| Fun03: setup_samEntry
|   - allocates memory for a samEntry structure (call
|     after init_samEntry (fun02))
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to allocate memory to
| Output:
|   - Allocates:
|     o memory for seqStr, qStr, cigTypStr, cigArySI,
|       and extraStr
|   - Returns:
|     o 0 for success
|     o def_memErr_samEntry for memory error
\-------------------------------------------------------*/
signed char
setup_samEntry(
   struct samEntry *samSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03: setup_samEntry
   '   o fun03 sec01:
   '     - allocate sequence memory
   '   o fun03 sec02:
   '     - allocate q-score memory
   '   o fun03 sec03:
   '     - allocate cigar type (snp/match/ins/del) memory
   '   o fun03 sec04:
   '     - allocate cigar count array
   '   o fun03 sec05:
   '     - allocate extra array
   '   o fun03 sec06:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec01:
    ^   - allocate sequence memory
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samSTPtr->seqStr)
       free(samSTPtr->seqStr);
  
    samSTPtr->seqStr = 0;

    samSTPtr->seqStr =
       calloc(
          1024,
          sizeof(char)
       ); /*allocate some memory for sequence*/

    if(! samSTPtr->seqStr)
       goto memErr_fun03_sec06;

    samSTPtr->lenSeqBuffUI = 1024;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec02:
    ^   - allocate q-score memory
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    if(samSTPtr->qStr)
       free(samSTPtr->qStr);
  
    samSTPtr->qStr = 0;

    samSTPtr->qStr =
       calloc(
          1024,
          sizeof(char)
       ); /*allocate some memory for q-score*/

    if(! samSTPtr->qStr)
       goto memErr_fun03_sec06;

    samSTPtr->lenQBuffUI = 1024;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec03:
    ^   - allocate cigar type (snp/match/ins/del) memory
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    if(samSTPtr->cigTypeStr)
       free(samSTPtr->cigTypeStr);
  
    samSTPtr->cigTypeStr = 0;

    samSTPtr->cigTypeStr =
       calloc(
          256,
          sizeof(char)
       ); /*allocate memory for cigary type*/

    if(! samSTPtr->cigTypeStr)
       goto memErr_fun03_sec06;

    samSTPtr->lenCigBuffUI = 256;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec04:
    ^   - allocate cigar count array
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    if(samSTPtr->cigArySI)
       free(samSTPtr->cigArySI);
  
    samSTPtr->cigArySI = 0;

    samSTPtr->cigArySI =
       calloc(
          256,
          sizeof(sint)
       ); /*allocate memory for cigar counts*/

    if(! samSTPtr->cigArySI)
       goto memErr_fun03_sec06;
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec05:
    ^   - allocate extra array
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(samSTPtr->extraStr)
       free(samSTPtr->extraStr);
  
    samSTPtr->extraStr = 0;

    samSTPtr->extraStr =
       calloc(
          1024,
          sizeof(char)
       ); /*allocate memory for extra entry*/

    if(! samSTPtr->extraStr)
       goto memErr_fun03_sec06;

    samSTPtr->lenExtraBuffUI = 1024;
    
    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun03 Sec06:
    ^   - return
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    return 0;

    memErr_fun03_sec06:;
    return def_memErr_samEntry;
} /*setup_samEntry*/

/*-------------------------------------------------------\
| Fun04: freeStack_samEntry
| Use:
|   - Frees all variables in samEntry, but not samEntry
| Input:
|   - samSTPtr
|     o Pointer to samEntry struct to free the interal
|       memory of
| Output:
|   - Frees:
|     o allocated memory in samSTPtr
\-------------------------------------------------------*/
void
freeStack_samEntry(
   struct samEntry *samSTPtr
){
    blank_samEntry((samSTPtr));
    
    free((samSTPtr)->seqStr);
    (samSTPtr)->seqStr = 0;
    (samSTPtr)->lenSeqBuffUI = 0;
    
    free((samSTPtr)->qStr);
    (samSTPtr)->qStr = 0;
    (samSTPtr)->lenQBuffUI = 0;
    
    free((samSTPtr)->cigTypeStr);
    (samSTPtr)->cigTypeStr = 0;
    (samSTPtr)->lenCigBuffUI = 0;
    
    free((samSTPtr)->cigArySI);
    (samSTPtr)->cigArySI = 0;
    
    free((samSTPtr)->extraStr);
    (samSTPtr)->extraStr = 0;
    (samSTPtr)->lenExtraBuffUI = 0;
} /*freeStack_samEntry*/

/*-------------------------------------------------------\
| Fun05: freeHeap_samEntry
|   - Frees a samEntry struct
| Input:
|   - samSTPtr
|     o Pointer to Sam entry struct to free
| Output:
|   - Frees:
|     o samSTPtr and its variables with allocated memory
\-------------------------------------------------------*/
void freeHeap_samEntry(
   struct samEntry **samSTPtr
){
    freeStack_samEntry((*samSTPtr));
    free((*samSTPtr));
    (*samSTPtr) = 0;
} /*freeHeap_samEntry*/

/*-------------------------------------------------------\
| Fun06: mk_samEntry
|   - Makes an heap allocated samEntry structure
| Input:
| Output:
|   - Returns:
|     o Pointer to a samEntry structure
|     o 0 if had an memory error
\-------------------------------------------------------*/
struct samEntry *
mk_samEntry(
){
  struct samEntry *samST=malloc(sizeof(struct samEntry));
  uchar errUC = 0;
  
  if(samST)
  { /*If: I did not have a memory error*/
     init_samEntry(samST);
     errUC = setup_samEntry(samST);

     if(errUC)
        freeHeap_samEntry(&samST);
  } /*If: I did not have a memory error*/
  
  return samST;
} /*makeSamEntry*/

/*-------------------------------------------------------\
| Fun07: qhistToMed_samEntry
|   - Gets the median q-score for an histogram of q-scores
|     in a samStruct
| Input:
|   - samSTPTr:
|     o Pointer to samEntry struct with a q-score
|       histogram to find the median q-score for
| Output:
|   - Modifies:
|     o samSTPtr->medianQF to have the median q-score
\-------------------------------------------------------*/
void
qhistToMed_samEntry(
   struct samEntry *samSTPtr 
){
    unsigned int numBasesUI = 0;
    unsigned int midPointUL = samSTPtr->readLenUI >> 1;
    unsigned int uiQ = 0;

    for(uiQ = 0; uiQ < def_maxQ_samEntry; uiQ++)
    { /*Loop: through q-score histogram for midpoint*/
        numBasesUI += samSTPtr->qHistUI[uiQ];

        if(numBasesUI >= midPointUL)
        { /*if found the midpoint, then find the median*/
            if(numBasesUI > midPointUL || numBasesUI & 1)
                samSTPtr->medianQF = (float) uiQ;

            else
            { /*Else: even & 2 differnt Q's at midpoint*/
                numBasesUI = uiQ;
                ++uiQ;

                while(samSTPtr->qHistUI[uiQ++] == 0) ;

                --uiQ; /*account for being one off*/

                samSTPtr->medianQF =
                    (numBasesUI + uiQ) / ((float) 2);
            } /*Else: even & 2 differnt Q's at midpoint*/

            break;
        } /*if found the midpoint, then find the median*/
    } /*Loop: through q-score histogram for midpoint*/
} /*qhistToMed_samEntry*/


/*-------------------------------------------------------\
| Fun08: findQScores_samEntry
|   - Gets the median and mean q-scores from a samEntry
|     Structure.
| Input:
|   - samSTPTr:
|     o Pointer to samEntry struct to find the median and
|       mean q-scores for
| Output:
|   - Modifies:
|     o samSTPtr->medianQF to have the median q-score
|     o samSTPtr->meanQF to have the mean q-score
|     o samSTPtr->qHistUI to have a histogram of q-scores
|     o samSTPtr->sumQUL to have the sum of all q-scores
\-------------------------------------------------------*/
void
findQScores_samEntry(
   struct samEntry *samSTPtr
){
    ulong qScoresUL = 0;
    ulong *qPtrUL = (ulong *) (samSTPtr)->qStr;

    uchar *scoreAryUC = 0;
    uint uiQScore = 0;
    uint uiChar = 0;
    
    ulong qAdjustUL =
       mkDelim_ulCp((schar) def_adjQ_samEntry);

    /*Find the number of q-score characters in buffer*/
    for(
       uiQScore = 0;
       uiQScore <
          ((samSTPtr)->readLenUI >> def_shiftULBy_ulCp);
       ++uiQScore
    ) { /*Loop: Update the q-score historgram and sum*/
       qScoresUL = qPtrUL[uiQScore] - qAdjustUL;

       scoreAryUC =
          (uchar *)
          &qScoresUL;
       
       for(
          uiChar = 0;
          uiChar < def_charInUL_ulCp;
          ++uiChar
       ){ /*Loop: Get the q-score entries*/
         ++(samSTPtr)->qHistUI[scoreAryUC[uiChar]];
         (samSTPtr)->sumQUL += (uchar) scoreAryUC[uiChar];
       } /*Loop: Get the q-score entries*/
    } /*Loop: Update the q-score historgram and sum*/
    
    uiQScore = (samSTPtr)->readLenUI;
    scoreAryUC = (uchar *) (samSTPtr)->qStr;
    
    for(
       uiQScore -=
          ((samSTPtr)->readLenUI & def_modUL_ulCp);
       uiQScore < (samSTPtr)->readLenUI;
       ++uiQScore
    ) { /*Loop: Copy the q-score entries*/
       (samSTPtr)->qStr[uiQScore] = scoreAryUC[uiQScore];

        uiChar = scoreAryUC[uiQScore] - def_adjQ_samEntry;
       
       ++(samSTPtr)->qHistUI[uiChar];
       (samSTPtr)->sumQUL += uiChar;
    } /*Loop: Copy the q-score entries*/
    
    (samSTPtr)->meanQF =
       (float) (samSTPtr)->sumQUL /(samSTPtr)->readLenUI;
    
    qhistToMed_samEntry((samSTPtr));
} /*findQScores_samEntry*/

/*-------------------------------------------------------\
| Fun09: cpQEntry_samEntry
|   - Copies q-scores from a string into a samEntry
|     structure
| Input:
|   - samSTPtr:
|     o Pionter to sam entry struct to copys q-scores to
|   _ cpQStr:
|     o C-string with q-scores to copy to samSTPtr
|   - blankQHistBl:
|     o 1: Blank q-score vars (histogram/sum/mean/median)
|     o 0: do not blank the q-score variables
| Output:
|   - Mofidies:
|     o qStr in samSTPtry to have the q-scores
|     o medianQF in samSTPtr to have the median q-score
|     o meanQF in samSTPtr to have the mean q-score
|     o qHistUI in samSTPtr to have histogram of q-scores
|     o samQUL in samSTPtr to have sum off all q-scores
|   - Returns
|     o The value in samSTPtr->readLenUI
\-------------------------------------------------------*/
unsigned int
cpQEntry_samEntry(
   struct samEntry *samSTPtr, /*Copy q-scores to*/
   signed char *cpQStr,       /*q-scores to copy*/
   signed char blankQHistBl   /*1: to blank q-score hist*/
){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  ' Fun09 TOC: cpQEntry_samEntry
  '   o fun09 sec01:
  '     - Variable declerations
  '   o fun09 sec02:
  '     - Check and if asked blank the q-score values
  '   o fun09 sec03:
  '     - Copy q-scores using unsigned longs
  '   o fun09 sec04:
  '     - Copy last q-scores I could not copy with
  '       unsigned longs
  '   o fun09 sec05:
  '     - Find the median and mean q-scores
  \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec01:
  ^   - Variable declerations
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  uchar *tmpStr = 0;
  uint uiQ = 0;
  uint uiChar = 0;
  ulong qAdjustUL=mkDelim_ulCp((schar) def_adjQ_samEntry);
  ulong *cpPtrUL = (ulong *) (cpQStr);
  ulong *dupPtrUL = (ulong *) samSTPtr->qStr;
  ulong qScoreUL = 0;
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec02:
  ^   - Check and if asked blank the q-score values
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  if(blankQHistBl)
  { /*If: I need to blank the q-score histograms*/
    for(uiChar = 0; uiChar < def_maxQ_samEntry; ++uiChar)
       samSTPtr->qHistUI[uiChar] = 0;

    samSTPtr->sumQUL = 0;
    samSTPtr->medianQF = 0;
    samSTPtr->meanQF = 0;
  } /*If: I need to blank the q-score histograms*/

  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec03:
  ^   - Copy q-scores using unsigned longs
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

  for(
     uiQ = 0;
     uiQ < (samSTPtr->readLenUI >> def_shiftULBy_ulCp);
     ++uiQ
  ) { /*Loop: Copy the q-score entries*/
     dupPtrUL[uiQ] = cpPtrUL[uiQ];
     qScoreUL = dupPtrUL[uiQ] - qAdjustUL;
     tmpStr = (uchar *) &qScoreUL;
     
     for(
        uiChar = 0;
        uiChar < def_charInUL_ulCp;
        ++uiChar
     ) { /*Loop: Get the q-score entries*/
        ++samSTPtr->qHistUI[tmpStr[uiChar]];
        samSTPtr->sumQUL += tmpStr[uiChar];
     } /*Loop: Get the q-score entries*/
  } /*Loop: Copy the q-score entries*/
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec04:
  ^   - Copy last q-scores I could not copy with unsgined
  ^     longs
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  
  uiQ = (samSTPtr)->readLenUI;
  
  for(
     uiQ -= ((samSTPtr)->readLenUI & def_modUL_ulCp);
     uiQ < (samSTPtr)->readLenUI;
     ++uiQ
  ) { /*Loop: Copy the q-score entries*/
     samSTPtr->qStr[uiQ] = cpQStr[uiQ];
     qScoreUL = (uchar) cpQStr[uiQ] - def_adjQ_samEntry;

     ++samSTPtr->qHistUI[qScoreUL];
     samSTPtr->sumQUL += qScoreUL;
  } /*Loop: Copy the q-score entries*/
  
  samSTPtr->qStr[uiQ] = '\0';
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun09 Sec05:
  ^   - Find the median and mean q-scores
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
  
  samSTPtr->meanQF = (float) samSTPtr->sumQUL;
  samSTPtr->meanQF /= (float) samSTPtr->readLenUI;
  
  qhistToMed_samEntry(samSTPtr);
  return uiQ;
} /*cpQEntry_samEntry*/

/*-------------------------------------------------------\
| Fun10: get_samEntry
|  - Reads in a single line from a sam file
| Input:
|  - samSTPtr:
|    o pionter to samEntry structure to store the sam file
|      line in. This structure should be initialized
|  - buffStr:
|    o buffer to read the sam file line temporarly into.
|      this is resized if needed. You can input NULL to
|      create a new buffer.
|  - lenBuffUL:
|    o length of buffStr (updated if buffStr is resized)
|  - samVoidFILE:
|    o sam file to read a line from. The void is so that
|      I can use samFILE in the function.
| Output:
|  - Modifies:
|    o samSTPtr to have the next line
|      - Comments are in extraStr
|    o samFILE to be on the next line
|    o buffStr to hold a sam file line (resized if needed)
|    o lenBuffUL to hold the resized length of buffStr
|  - Returns:
|    o 0 for success
|    o def_EOF_samEntry for EOF (End Of File)
|    o def_memErr_samEntry for memory errors
\-------------------------------------------------------*/
signed char
get_samEntry(
   struct samEntry *samSTPtr,
   signed char **buffStr,
   ulong *lenBuffUL,
   void *samVoidFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun10 TOC: get_samEntry
   '   - Reads a single line from a sam file into samSTPtr
   '   o fun10 sec01:
   '     - Variable declerations
   '   o fun10 sec02:
   '     - Read in one line from the sam file
   '   o fun10 sec03:
   '     - Get the query id from the buffer
   '   o fun10 sec04:
   '     - Get the flag
   '   o fun10 sec05:
   '     - REad in the reference name/id
   '   o fun10 sec06:
   '     - Get the reference position
   '   o fun10 sec07:
   '     - Get the mapping quality
   '   o fun10 sec08:
   '     - Get the cigar entry
   '   o fun10 sec09:
   '     - Get the RNEXT entry
   '   o fun10 sec10:
   '     - Get the PNEXT entry
   '   o fun10 sec11:
   '     - Get the TLEN entry
   '   o fun10 sec12:
   '     - Get the sequence entry
   '   o fun10 sec13:
   '     - Get the q-score entry
   '   o fun10 sec14:
   '     - Copy the extra entry; after strict sam entries
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec01:
   ^   - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ushort extraBuffUS = 4096;

   ulong oldLenUL = 0;
   ulong tabUL = mkDelim_ulCp((schar) '\t');
   signed char *tmpStr = 0;
   signed char *iterStr = 0;

   FILE *samFILE = (FILE *) samVoidFILE;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec02:
   ^   - Read in one line from the sam file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_samEntry(samSTPtr);

   if(*buffStr == 0)
   { /*If: I need to create a buffer*/
      *lenBuffUL = 0;
      *buffStr = malloc(extraBuffUS * sizeof(char));

      if(*buffStr == 0)
         goto memErr_fun10_sec14;

      *lenBuffUL = extraBuffUS;
   } /*If: I need to create a buffer*/

   tmpStr =
      (schar *)
      fgets(
         (char *) *buffStr,
         *lenBuffUL,
         samFILE
      );

   if(! tmpStr )
      goto eof_fun10_sec14; /*EOF*/

   goto checkEOL_fun10_sec02_samEntry;

   oldLenUL += endLine_ulCp(tmpStr);
   tmpStr = *buffStr + oldLenUL;

   while(*tmpStr != '\n')
   { /*Loop: Find the length of  the line*/
      *lenBuffUL <<= 1;
        /*this is a little agressive for memory usage, but
        `   it is fast*/

      tmpStr =
         realloc(
            *buffStr,
            (*lenBuffUL + 1) * sizeof(schar)
         );

      /*check for memory errors; let user handle
      `   freeing buffStr when have memory errors
      */
      if(! tmpStr)
         goto memErr_fun10_sec14;

      *buffStr = tmpStr;
      tmpStr = *buffStr + oldLenUL;

      tmpStr =
         (schar *)
         fgets(
            (char *) tmpStr,
            *lenBuffUL >> 1,
            samFILE
         );

      if(! tmpStr)
         break; /*End of file*/

      checkEOL_fun10_sec02_samEntry:;

      oldLenUL += endLine_ulCp(tmpStr);
      tmpStr = *buffStr + oldLenUL;
   } /*Loop: Find the length of  the line*/
   
   /*Add null in for line end (not end of file)*/
   *tmpStr = '\0';
   iterStr = *buffStr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec03:
   ^   - Get the query id from the buffer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*This is a comment*/
   if(iterStr[0] == '@') goto extraEntry;
 
   /*Query id, single byte copy may be better here,
    ` hard to tell
   */
   samSTPtr->lenQryIdUC =
     (uchar)
      cpDelim_ulCp(
         samSTPtr->qryIdStr,
         iterStr,
         tabUL,
         '\t'
      ); /*Copy the reference id/name*/

   iterStr += samSTPtr->lenQryIdUC + 1; /*+1 get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec04:
   ^   - Get the flag for the alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr +=
      strToUS_base10str(
         iterStr,
         &samSTPtr->flagUS
      ); /*Get the flag for the alignment*/

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec05:
   ^   - REad in the reference name/id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Not clear wich version is faster here*/
   samSTPtr->lenRefIdUC =
      (uchar)
      cpDelim_ulCp(
         samSTPtr->refIdStr,
         iterStr,
         tabUL,
         '\t'
      ); /*Copy the reference id/name*/
   
   iterStr += samSTPtr->lenRefIdUC + 1; /*+1 get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec06:
   ^   - Get the reference position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr +=
      strToUI_base10str(
         iterStr,
         &samSTPtr->refStartUI
      ); /*Get the starting base in the reference*/

    /*Convert the starting positionto index 0*/
   samSTPtr->refStartUI -= (samSTPtr->refStartUI > 0);

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec07:
   ^   - Get the mapping quality
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr +=
      strToUC_base10str(
         iterStr,
         &samSTPtr->mapqUC
      ); /*Get the mapping quality of the alignment*/

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec08:
   ^   - Get the cigar entry
   ^   o fun10 sec08 sub01:
   ^     - Check if there is a cigar entry
   ^   o fun10 sec08 sub02:
   ^     - Read in the cigar entry
   ^   o fun10 sec08 sub03:
   ^     - Count number of matchs/snps/dels/inss/masks in
   ^       the cigar entry
   ^   o fun10 sec08 sub04:
   ^     - Get the read lengths from the cigar entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun10 Sec08 Sub01:
   *   - Check if there is a cigar entry
   \*****************************************************/

   if(iterStr[0] == '*')
   { /*If: the cigar entry was not present*/
      samSTPtr->cigTypeStr[0] = '*';
      samSTPtr->cigArySI[0] = 0;
      iterStr += 2;
      (samSTPtr)->lenCigUI = 1;

      goto rNextEntry;
   } /*If: the cigar entry was not present*/

   /*****************************************************\
   * Fun10 Sec08 Sub02:
   *   - Read in the cigar entry
   \*****************************************************/
   
   while(iterStr[0] > 32)
   { /*Loop: Read in the cigar entry*/

      /*Using -1 to account for the null I will add at the
      `end
      */
      if((samSTPtr)->lenCigUI >= samSTPtr->lenCigBuffUI-1)
      { /*If: I need to increase the cigar buff size*/
         samSTPtr->lenCigBuffUI <<= 1;

         tmpStr =
            realloc(
               samSTPtr->cigTypeStr,
               samSTPtr->lenCigBuffUI * sizeof(char)
            ); /*Resize the type cigar buffer*/

          if(tmpStr == 0)
             goto memErr_fun10_sec14;

          samSTPtr->cigTypeStr = tmpStr;

         tmpStr = (schar *)
            realloc(
               samSTPtr->cigArySI,
               samSTPtr->lenCigBuffUI * sizeof(int)
            ); /*Resize the value cigar buffer*/

          if(tmpStr == 0)
             goto memErr_fun10_sec14;

          samSTPtr->cigArySI = (int *) tmpStr;
      } /*If: I need to increase the cigar buff size*/

      /*Record the cigar entry*/
      iterStr +=
         strToSI_base10str(
            iterStr,
            &samSTPtr->cigArySI[samSTPtr->lenCigUI]
          ); /*Get the number of bases for this type*/

      samSTPtr->cigTypeStr[samSTPtr->lenCigUI]=iterStr[0];

      /**************************************************\
      * Fun10 Sec08 Sub03:
      *   - Count number of matchs/snps/dels/inss/masks in
      *     the cigar entry
      \**************************************************/

      switch(iterStr[0])
      { /*Switch: Check the cigar entry type*/
         case '=':
         /*Case: This was a match*/
             samSTPtr->numMatchUI +=
                samSTPtr->cigArySI[samSTPtr->lenCigUI];

             break;
         /*Case: This was a match*/

         case 'M':
         /*Case: This is an snp or match*/
             samSTPtr->numMatchUI +=
                samSTPtr->cigArySI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an snp or match*/

         case 'X':
         /*Case: This is an snp*/
             samSTPtr->numSnpUI +=
                samSTPtr->cigArySI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an snp*/

         case 'I':
         /*Case: This is an insertion*/
             samSTPtr->numInsUI +=
                samSTPtr->cigArySI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an insertion*/

         case 'D':
         /*Case: This is an deletion*/
             samSTPtr->numDelUI +=
                samSTPtr->cigArySI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an deletion*/

         case 'S':
         /*Case: This is an softmasked region*/
             samSTPtr->numMaskUI +=
                samSTPtr->cigArySI[samSTPtr->lenCigUI];

             break;
         /*Case: This is an softmasked region*/
      } /*Switch: Check the cigar entry type*/

      ++iterStr; /*Get off the tab*/
      ++(samSTPtr)->lenCigUI;
   } /*Loop: Read in the cigar entry*/

   /*****************************************************\
   * Fun10 Sec08 Sub04:
   *   - Get the read lengths from the cigar entries
   \*****************************************************/

    samSTPtr->cigArySI[samSTPtr->lenCigUI] = 0;
    samSTPtr->cigTypeStr[samSTPtr->lenCigUI] = '\0';

   samSTPtr->readLenUI =
        samSTPtr->numMatchUI
      + samSTPtr->numSnpUI
      + samSTPtr->numInsUI
      + samSTPtr->numMaskUI;

   samSTPtr->alnReadLenUI =
        samSTPtr->numMatchUI
      + samSTPtr->numSnpUI
      + samSTPtr->numDelUI;

   samSTPtr->refEndUI = samSTPtr->refStartUI;
   samSTPtr->refEndUI += samSTPtr->alnReadLenUI;

   samSTPtr->refEndUI -= (samSTPtr->alnReadLenUI > 0);
      /*-1 from (alnReadLen > 0) converts to index 0*/

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec09:
   ^   - Get the RNEXT entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   rNextEntry:

   /*Not sure which is better here*/
   samSTPtr->lenRNextUC =
      (uchar)
      cpDelim_ulCp(
         samSTPtr->rNextStr,
         iterStr,
         tabUL,
         '\t'
      ); /*Copy the query id/name*/

   iterStr += samSTPtr->lenRNextUC + 1; /*+1 get off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec10:
   ^   - Get the PNEXT entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr +=
      strToSI_base10str(
         iterStr,
         &samSTPtr->pNextSI
       ); /*Get the number of bases for this type*/

   samSTPtr->pNextSI -= (samSTPtr->pNextSI > 0);

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec11:
   ^   - Get the TLEN entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   iterStr +=
      strToSI_base10str(
         iterStr,
         &samSTPtr->tLenSI
       ); /*Get the number of bases for this type*/

   ++iterStr; /*Get off the tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec12:
   ^   - Get the sequence entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(samSTPtr->readLenUI == 0 && iterStr[0] != '*')
      samSTPtr->readLenUI =
         lenStr_ulCp(iterStr, tabUL, '\t');

   else if(iterStr[0] == '*')
   { /*Else If: There  is no sequence entry*/
      samSTPtr->seqStr[0] = '*';
      samSTPtr->seqStr[1] = '\0';
      iterStr += 2;

      goto noQEntry;
   } /*Else If: There  is no sequence entry*/

   if(samSTPtr->readLenUI + 8 > samSTPtr->lenSeqBuffUI)
   { /*If: I need to resize sequence & q-score buffers*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         malloc((samSTPtr->readLenUI + 9) * sizeof(char));

      if(!samSTPtr->seqStr)
         goto memErr_fun10_sec14;

      free(samSTPtr->qStr);
      samSTPtr->qStr = 0;

      samSTPtr->qStr =
         malloc((samSTPtr->readLenUI + 9) * sizeof(char));

      if(!samSTPtr->qStr)
         goto memErr_fun10_sec14;

      samSTPtr->lenSeqBuffUI = samSTPtr->readLenUI;
      samSTPtr->lenQBuffUI = samSTPtr->readLenUI;
   } /*If: I need to resize sequence & q-score buffers*/

   cpLen_ulCp(
      samSTPtr->seqStr,
      iterStr,
      samSTPtr->readLenUI
   );

   iterStr += samSTPtr->readLenUI + 1; /*+1 gets off tab*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec13:
   ^   - Get the q-score entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(iterStr[0] == '*' && iterStr[1] == '\t')
   { /*If: there is no q-score entry*/
      noQEntry:

      samSTPtr->qStr[0] = '*';
      samSTPtr->qStr[1] = '\0';
      iterStr += 2;
   } /*If: there is no q-score entry*/

   else
   { /*Else: is a q-score entry*/
      iterStr +=
         cpQEntry_samEntry(
            samSTPtr,
            iterStr,
            0
         ) + 1;
   } /*Else: is a q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec14:
   ^   - Copy the extra entry; after strict sam entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   extraEntry:

   if(! iterStr || *iterStr == '\0')
   { /*If: no extra entry*/
      if(samSTPtr->lenExtraUI > 1)
      { /*If: I need to resize the buffer*/
         if(samSTPtr->extraStr)
            free(samSTPtr->extraStr);

         samSTPtr->extraStr = 0;

         samSTPtr->extraStr = malloc(10 * sizeof(schar));
         
         if(samSTPtr->extraStr == 0)
            goto memErr_fun10_sec14;

         samSTPtr->extraStr[0] = '\0';

         samSTPtr->lenExtraUI = 0;
      } /*If: I need to resize the buffer*/
   } /*If: no extra entry*/

   else
   { /*Else: have extra entry*/
      /*not sure if char or ul copy better here*/
      samSTPtr->lenExtraUI =
         lenStr_ulCp(
            iterStr,
            0,
            '\0'
         );

      if(samSTPtr->lenExtraUI > samSTPtr->lenExtraBuffUI)
      { /*If: I need to resize the buffer*/
         free(samSTPtr->extraStr);
         samSTPtr->extraStr = 0;

         samSTPtr->extraStr =
            calloc(
               samSTPtr->lenExtraUI + 9,
               sizeof(schar)
            );
         
         if(samSTPtr->extraStr == 0)
            goto memErr_fun10_sec14;
      } /*If: I need to resize the buffer*/

      cpLen_ulCp(
         samSTPtr->extraStr,
         iterStr,
         samSTPtr->lenExtraUI
      ); /*Copy the extra entry*/
   } /*Else: have extra entry*/

   return 0;

   memErr_fun10_sec14:;
   return def_memErr_samEntry;

   eof_fun10_sec14:;
   return def_EOF_samEntry;
} /*get_samEntry*/

/*-------------------------------------------------------\
| Fun11: findRefPos_samEntry
|   - Find an reference coordinate in an sequence in
|     an sam entry structure
| Input:
|   - samSTPtr:
|     o Pointer to samEntry structure with sequence
|       (query) to find target reference position in
|   - siCig:
|     o pionter to signed int with index of the cigar
|       entry I am on
|   - cigNtSI:
|     o pionter to a singed int with the number of bases
|       left in the cigar entry
|   - targPosSI:
|     o Reference position I wish to find
|   - refPosSI:
|     o pointer to signed int with the current position
|       in reference sequence
|   - seqPosSI:
|     o pointer to signed int with the current position
|       in query sequence
| Output:
|   - Modifies:
|     o siCig to point to the next open cigar entry
|       - will be > samSTPtr->lenCigUI when the sequence
|         does not end at at targPosSI
|     o cigNtSI to have the number of bases remianing
|       in the current siCig entry
|     o refPosSI to be the index to the position in targSI
|     o seqPosSI to be the index of the sequence base at
|       refPosSI (check cigar to see if deletion)
|   - Returns:
|     o Index of last cigar entry I incurmented
\-------------------------------------------------------*/
signed int
findRefPos_samEntry(
   struct samEntry *samSTPtr,/*has cigar and sequence*/
   signed int *siCig,        /*current cigar entry on*/
   signed int *cigNtSI,      /*# of bases left in cigar*/
   signed int targPosSI,     /*position looking for*/
   signed int *refPosSI,     /*current ref position*/
   signed int *seqPosSI      /*current sequence position*/
){/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
  ` Fun11 TOC: findRefPos_samEntry
  `   - Find an reference coordinate in an sequence in
  `     an sam entry structure
  `   o fun11 sec01:
  `     - Start loop and check insertions/soft masking
  `   o fun11 sec02:
  `     - Move position in deletion cases
  `   o fun11 sec03:
  `     - Move position for snp/match cases
  `   o fun11 sec04:
  `     - Move to the next cigar entry
  \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^   - Start loop and check insertions/soft masking
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int lastCigMacSI = *siCig;

   if(samSTPtr->cigTypeStr[*siCig] == 'S')
      goto softMask_fun11_sec01;

   if(samSTPtr->cigTypeStr[*siCig] == 'H')
      goto hardMask_fun11_sec03;

   /*Check if I did a full cigar entry move*/
   while(*refPosSI < targPosSI)
   { /*Loop: till I am on the target base*/
      lastCigMacSI = *siCig;

      switch(samSTPtr->cigTypeStr[*siCig])
      { /*Switch: check what the next entry is*/
         case 'S':
         case 'I':
         /*Case: Softmasking or insertions*/
            softMask_fun11_sec01:;
            *seqPosSI += *cigNtSI;
            ++(*siCig);
            *cigNtSI = 0;
            break;
         /*Case: Softmasking or insertions*/

         /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
         ^ Fun11 Sec02:
         ^   - Move position in deletion cases
         \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

         case 'D':
         /*Case: Deletion*/
            *refPosSI += *cigNtSI;

            if(*refPosSI <= (targPosSI))
            { /*If: I have not found target position*/
               ++(*siCig);
               *cigNtSI = 0;
            } /*If: I have not found target position*/

            else
            { /*Else: I overshot the target*/
               *cigNtSI = *refPosSI - targPosSI;
               *refPosSI -= *cigNtSI; /*correct overshot*/
            } /*Else: I overshot the target*/

            break;
         /*Case: Deletion*/

         /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
         ^ Fun11 Sec03:
         ^   - Move position for snp/match cases\
         \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

         case 'M':
         case 'X':
         case '=':
         /*Case: match (M or =) or snp (M or X)*/
            *refPosSI += *cigNtSI;
            *seqPosSI += *cigNtSI;

            if(*refPosSI <= targPosSI)
            { /*If: I have not found target position*/
               *cigNtSI = 0;
               ++(*siCig);
            } /*If: I have not found target position*/

            else 
            { /*Else: I overshot the target*/
               *cigNtSI = *refPosSI - targPosSI;

               *refPosSI -= *cigNtSI; /*correct overshot*/
               *seqPosSI -= *cigNtSI; /*correct overshot*/
            } /*Else: I overshot the target*/

            break;
         /*Case: match (M or =) or snp (M or X)*/

         default:
         /*Case: hard mask of some kind*/
            hardMask_fun11_sec03:;
            ++(*siCig);
            *cigNtSI = 0;
         /*Case: hard mask of some kind*/
      } /*Switch: check what the next entry is*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun11 Sec04:
      ^   - Move to the next cigar entry
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      if(*siCig >= (signed int) samSTPtr->lenCigUI)
         break; /*End of the sequence*/

      if(*cigNtSI == 0)
         *cigNtSI = samSTPtr->cigArySI[*siCig];
   } /*Loop: till I am on the target base*/

   return lastCigMacSI;
} /*findRefPos_samEntry*/

/*-------------------------------------------------------\
| Fun12: p_samEntry
|   - Prints the sam file entry to a file. This does not
|     print any extra stats that were found.
| Input:
|   - samST
|     o Pointer to samEntry struct with sam entry to print
|   - buffStr:
|     o Pointer to c-string buffer to temporarly hold the
|       cigar entry (speeds things up)
|   - lenBuffUL:
|     o Current length of buffer, adjusted if buffStr is
|       expanded
|   - pNoNewLineBl:
|     o 1: do not print a new line after; you will do this
|     o 0: end of sam entry, print a new line
|   - outFILE:
|     o File to print the sam entry to
| Output:
|  - Prints:
|    o Sam file entry in samST to outFILE.
|  - Modifies:
|    o inceases the size of buffStr, if is not 8x the
|      cigar length
|    o Sets lenBuffUL to the new buffStr size when buffStr
|      is resized
|  - Returns:
|    o 0 for no problems
|    o def_memErr_samEntry for memory errors
\-------------------------------------------------------*/
signed char
p_samEntry(
   struct samEntry *samSTPtr,
   signed char **buffStr,
   ulong *lenBuffUL,
   signed char pNoNewLineBl,
   void *outFILE
){
   uint uiCig = 0;
   schar *tmpStr = *buffStr;
   ulong maxLenUL = 0; 

   maxLenUL =
        (samSTPtr)->lenQryIdUC
      + (samSTPtr)->lenRefIdUC
      + ((samSTPtr)->lenCigUI << 3)
      + (samSTPtr)->lenRNextUC
      + ((samSTPtr)->readLenUI << 1)
      + (samSTPtr)->lenExtraBuffUI
      + 128; /*128 to account for numbers/tabs*/

   if(*(lenBuffUL) < maxLenUL)
   { /*If: I need to add more buffer*/
      free(*buffStr);
      *(lenBuffUL) = maxLenUL;
      *buffStr = malloc((*lenBuffUL + 1) * sizeof(char));
   } /*If: I need to add more buffer*/

   tmpStr = *buffStr;
   
   if(! (*buffStr))
      goto memErr_fun12;
   
   else if(
         (samSTPtr)->extraStr[0] == '@'
      && (samSTPtr)->qryIdStr[0] == '\0'
   ){
      fprintf(
         (FILE *) (outFILE),
         "%s\n",
         (samSTPtr)->extraStr
      ); /*Print out the header (in extra since commnet)*/

      goto noErr_fun12;
   } /*Else If: this was a header*/
   
   /*Copy the query id to the buffer*/
   cpLen_ulCp(
      tmpStr,
      (samSTPtr)->qryIdStr,
      (samSTPtr)->lenQryIdUC
   );

   tmpStr += (samSTPtr)->lenQryIdUC;
   *tmpStr++ = '\t';

   /*Copy the flag over*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->flagUS);
   *tmpStr++ = '\t';

   /*Copy the referenced id to the buffer*/
   cpLen_ulCp(
      tmpStr,
      (samSTPtr)->refIdStr,
      (samSTPtr)->lenRefIdUC
   );

   tmpStr += (samSTPtr)->lenRefIdUC;
   *tmpStr++ = '\t';

   /*Reference position*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->refStartUI + 1);
   *tmpStr++ = '\t';

   /*mapq*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->mapqUC);
   *tmpStr++ = '\t';

   /*Check if there is a cigar entry*/
   if((samSTPtr)->cigTypeStr[0]=='*') *tmpStr++ ='*';

   else
   { /*Else: convert the cigar entry*/
      for(
         uiCig=0;
         uiCig < (samSTPtr)->lenCigUI;
         ++uiCig
      ){ /*Loop: Convert cigar to string*/
         tmpStr +=
            numToStr(
               tmpStr,
               (samSTPtr)->cigArySI[uiCig]
            );
         *tmpStr++ = (samSTPtr)->cigTypeStr[uiCig];
      } /*Loop: Convert cigar to string*/
   } /*Else: convert the cigar entry*/

   *tmpStr++ = '\t';

   /*RNEXT*/
   cpLen_ulCp(
      tmpStr,
      (samSTPtr)->rNextStr,
      (samSTPtr)->lenRNextUC
   );

   tmpStr += (samSTPtr)->lenRNextUC;
   *tmpStr++ = '\t';

   /*PNEXT*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->mapqUC);
   *tmpStr++ = '\t';

   /*TLEN*/
   tmpStr += numToStr(tmpStr, (samSTPtr)->mapqUC);
   *tmpStr++ = '\t';

   /*Copy the sequence to the buffer*/
   cpLen_ulCp(
      tmpStr,
      (samSTPtr)->seqStr,
      (samSTPtr)->readLenUI
   );

   tmpStr += (samSTPtr)->readLenUI;
   *tmpStr++ = '\t';

   /*Copy the q-score entry to the buffer*/
   if((samSTPtr)->qStr[1] == '\0') *tmpStr++ = '*';

   else
   { /*Else: there is an q-score entry*/
      cpLen_ulCp(
         tmpStr,
         (samSTPtr)->qStr,
         (samSTPtr)->readLenUI
      );

      tmpStr += (samSTPtr)->readLenUI;
   } /*Else: there is an q-score entry*/

   if((samSTPtr)->lenExtraUI)
   { /*If: have extra items*/
      *tmpStr++ = '\t';

      /*Copy the extra entry*/
      cpLen_ulCp(
         tmpStr,
         (samSTPtr)->extraStr,
         (samSTPtr)->lenExtraUI
      );

      tmpStr += (samSTPtr)->lenExtraUI;
   } /*If: have extra items*/

   if(! pNoNewLineBl)
      *tmpStr++ = '\n';

   fwrite(
      *buffStr,
      sizeof(char),
      (tmpStr - *buffStr),
      (FILE *) (outFILE)
   );
   
   noErr_fun12:;
   return 0;

   memErr_fun12:;
   return def_memErr_samEntry;
} /*p_samEntry*/

/*-------------------------------------------------------\
| Fun13: pfq_samEntry
|   - Prints the sam entry as fastq entry to a fastq file
| Input:
|   - samST:
|     o Pointer to samEntry structure with fastq entry to
|       print out
|   - outFILE:
|     o Fastq file to print the new fastq entry to
| Output:
|   - Prints:
|     o fastq entry from samST to outFILE
\-------------------------------------------------------*/
void
pfq_samEntry(
   struct samEntry *samSTPtr,
   void *outFILE
){
   if(   (samSTPtr)->seqStr[1] != '\0'
      && (samSTPtr)->qStr[1] != '\0'
   ){ /*If: This entry can be printed out*/
      fprintf(
        (FILE *) (outFILE),
        "@%s ref=%s start=%u len=%u refAlnLen=%u",
        (samSTPtr)->qryIdStr,
        (samSTPtr)->refIdStr,
        (samSTPtr)->refStartUI + 1,
        (samSTPtr)->readLenUI,
        (samSTPtr)->alnReadLenUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  flag=%u mapq=%u match=%u snp=%u ins=%u",
        (samSTPtr)->flagUS,
        (samSTPtr)->mapqUC,
        (samSTPtr)->numMatchUI,
        (samSTPtr)->numSnpUI,
        (samSTPtr)->numInsUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  del=%u softMasked=%u meanQ=%f medianQ=%f\n",
        (samSTPtr)->numDelUI,
        (samSTPtr)->numMaskUI,
        (samSTPtr)->meanQF,
        (samSTPtr)->medianQF
      );
      fprintf(
        (FILE *) (outFILE),
        "%s\n+\n%s\n",
        (samSTPtr)->seqStr,
        (samSTPtr)->qStr
      );
   } /*If: This entry can be printed out*/
} /*p_samEntryAsFq*/

/*-------------------------------------------------------\
| Fun14: pfa_samEntry
|   - Prints the sam entry as fasta entry to a fasta file
| Input:
|   - samST:
|     o Pointer to samEntry structure with fastq entry to
|       print out
|   - outFILE:
|     o Fasta file to print the new fasta entry to
| Output:
|   - Prints:
|     o fasta entry from samST to outFILE
\-------------------------------------------------------*/
void
pfa_samEntry(
   struct samEntry *samSTPtr,
   void *outFILE
){
   if((samSTPtr)->seqStr[1] != '\0')
   { /*If: This entry can be printed out*/
      fprintf(
        (FILE *) (outFILE),
        ">%s ref=%s start=%u len=%u refAlnLen=%u",
        (samSTPtr)->qryIdStr,
        (samSTPtr)->refIdStr,
        (samSTPtr)->refStartUI + 1,
        (samSTPtr)->readLenUI,
        (samSTPtr)->alnReadLenUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  flag=%u mapq=%u match=%u snp=%u ins=%u",
        (samSTPtr)->flagUS,
        (samSTPtr)->mapqUC,
        (samSTPtr)->numMatchUI,
        (samSTPtr)->numSnpUI,
        (samSTPtr)->numInsUI
      );
      
      fprintf(
        (FILE *) (outFILE),
        "  del=%u softMasked=%u meanQ=%f medianQ=%f\n",
        (samSTPtr)->numDelUI,
        (samSTPtr)->numMaskUI,
        (samSTPtr)->meanQF,
        (samSTPtr)->medianQF
      );
      fprintf((outFILE), "%s\n", (samSTPtr)->seqStr);
   } /*If: This entry can be printed out*/
} /*pfa_samEntry*/

/*-------------------------------------------------------\
| Fun15: pstats_samEntry
|   - Prints out the stats in a samEntry struct to a file
| Input:
|   - samEntryST:
|     o Pointer to samEntry struct to print stats for
|   - pHeadBl:
|     o 1: Print the header for the stats tsv file
|     o 0: Do not print the header
|   - pNsBl:
|     o 1: find and print out the anonymous base counts
|     o 0: do not print out anonymous base counts
|   - outFILE:
|     o TSV (tab deliminated) file to print stats to
| Output:
|   - Prints:
|     o line with stats from samEntryStruct to file (tsv)
|   - Modifies:
|     o printHeaderChar to be 0 if set to 1
\-------------------------------------------------------*/
void
pstats_samEntry(
   struct samEntry *samSTPtr,
   signed char *pHeadBl,
   signed char pNsBl,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun15 TOC:
   '   - Prints stats in a samEntry struct to file
   '   o fun15 sec01:
   '     - variable declerations
   '   o fun15 sec02:
   '     - check if comment, if not check if print header
   '   o fun15 sec03:
   '     - print out general stats
   '   o fun15 sec04:
   '     - print matches, snps, ins, dels, and masking.
   '     - if asked also include anonymous bases
   '   o fun15 sec05:
   '     - print out the accuracy stats
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec01:
   ^   - variable declerations 
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint matchNCntUI = 0;
   uint snpNCntUI = 0;
   uint insNCntUI = 0;
   uint maskNCntUI = 0;
   uint *cntUI = 0;

   sint siCig = 0;
   sint numNtSI = 0;
   schar *tmpStr = 0;
   uchar ntUC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec02:
   ^   - check if comment, if not check if print header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(   (samSTPtr)->extraStr[0] != '@'
      && (samSTPtr)->qryIdStr[0] != '\0'
   ){ /*If: This is not a comment*/

      if(*(pHeadBl))
      { /*If: I need to print the header*/
        fprintf((FILE *) (outFILE), "Read\tRef\tFlag");
        fprintf((FILE *) (outFILE), "\tMapQ\tRefPos"); 
        fprintf((FILE *) (outFILE), "\tReadLength");
        fprintf((FILE *) (outFILE), "\tRefAlnLength");

        if((pNsBl))
        { /*If: I am printing the anonymous base counts*/
           fprintf((FILE *) (outFILE), "\tmatch_total");
           fprintf((FILE *) (outFILE), "\tmatch");
           fprintf((FILE *)(outFILE),"\tmatch_anonymous");

           fprintf((FILE *) (outFILE), "\tsnp_total");
           fprintf((FILE *) (outFILE), "\tsnp");
           fprintf((FILE *) (outFILE), "\tsnp_anonymous");

           fprintf((FILE *) (outFILE), "\tins_total");
           fprintf((FILE *) (outFILE), "\tins");
           fprintf((FILE *) (outFILE), "\tins_anonymous");

           fprintf((FILE *) (outFILE), "\tdel");

           fprintf((FILE *) (outFILE), "\tmask_total");
           fprintf((FILE *) (outFILE), "\tmask");
           fprintf((FILE *) (outFILE),"\tmask_anonymous");
        } /*If: I am printing the anonymous base counts*/

        else
        { /*Else: Print out the normal header*/
           fprintf((FILE *) (outFILE), "\tMatches\tSnps");
           fprintf((FILE *) (outFILE), "\tInss\tDels");
           fprintf((FILE *) (outFILE), "\tmask");
        } /*Else: Print out the normal header*/

        fprintf((FILE *) (outFILE), "\tMedianQ\tMeanQ\n");

        *(pHeadBl) = 0;
      } /*If: I need to print the header*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun15 Sec03:
      ^   - print out general stats
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      fprintf(
        (FILE *) (outFILE),
        "%s\t%s\t%u\t%u\t%u\t%u\t%u",
        (samSTPtr)->qryIdStr,
        (samSTPtr)->refIdStr,
        (samSTPtr)->flagUS,
        (samSTPtr)->mapqUC,
        (samSTPtr)->refStartUI + 1,
        (samSTPtr)->readLenUI,
        (samSTPtr)->alnReadLenUI
      );

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun15 Sec04:
      ^   - print matches, snps, ins, dels, and masking.
      ^   - if asked also include anonymous bases
      ^   o fun15 sec04 sub01:
      ^     - print anonymous bases; have sequence?
      ^   o fun15 sec04 sub02:
      ^     - find number anonymous bases
      ^   o fun15 sec04 sub03:
      ^     - print out anonymous base and other counts
      ^   o fun15 sec04 sub04:
      ^     - not printing anoynous bases do regular print
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      /**************************************************\
      * Fun15 Sec04 Sub01:
      *   - (print anonymous bases) check if have sequence
      \**************************************************/

      if( pNsBl )
      { /*If: I am finding the anonymous base counts*/
         tmpStr = (schar *) (samSTPtr)->seqStr;

         if(*tmpStr == '*')
         { /*If: no sequence to check for anonymous*/
            fprintf(
              (FILE *) (outFILE),
              "\t%u\tNA\tNA",
              (samSTPtr)->numMatchUI
            );

            fprintf(
              (FILE *) (outFILE),
              "\t%u\tNA\tNA",
              (samSTPtr)->numSnpUI
            );

            fprintf(
              (FILE *) (outFILE),
              "\t%u\tNA\tNA\t%u",
              (samSTPtr)->numInsUI,
              (samSTPtr)->numDelUI
            );

            fprintf(
              (FILE *) (outFILE),
              "\t%u\tNA\tNA",
              (samSTPtr)->numMaskUI
            );

            goto pqscore_fun15_sec05;
         } /*If: no sequence to check for anonymous*/

         /***********************************************\
         * Fun15 Sec04 Sub02:
         *   - find number anonymous bases
         \***********************************************/

         siCig = 0;

         while(siCig < (sint) (samSTPtr)->lenCigUI)
         { /*Loop: count number of anonymous bases*/
            numNtSI = (samSTPtr)->cigArySI[siCig];

            switch((samSTPtr)->cigTypeStr[siCig])
            { /*Switch: check what type of cigar entry*/
               case 'I':
               /*Case: Insertions*/
                  cntUI = &insNCntUI;
                  break;
               /*Case: Insertions*/

               case 'X':
               /*Case: Mismatches*/
                  cntUI = &snpNCntUI;
                  break;
               /*Case: Mismatches*/

               case '=':
               case 'M':
               /*Case: Matches or non-eqx entry*/
                  cntUI = &matchNCntUI;
                  break;
               /*Case: Matches or non-eqx entry*/

               case 'S':
               /*Case: softmasking*/
                  cntUI = &maskNCntUI;
                  break;
               /*Case: softmasking*/

               /*Handle missing information cases*/
               default:
               /*Case: deletions or hard masks*/
                  cntUI = 0;
                  numNtSI = 0;
               /*Case: deletions or hard masks*/
            } /*Switch: check what type of cigar entry*/

            while(numNtSI > 0)
            { /*Loop: Check each nucleotide*/
               ntUC = ntTo5Bit[(uchar) *tmpStr];
               (*cntUI) +=
                  (!!(ntUC & def_n_fithBit_ntTo5Bit));

               ++tmpStr;
               --numNtSI;
            } /*Loop: Check each nucleotide*/

            ++siCig;
         } /*Loop: count number of anonymous bases*/

         /***********************************************\
         * Fun15 Sec04 Sub03:
         *   - print anonymous base and other counts
         \***********************************************/

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u",
           (samSTPtr)->numMatchUI,
           (samSTPtr)->numMatchUI - matchNCntUI,
           matchNCntUI
         );

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u",
           (samSTPtr)->numSnpUI,
           (samSTPtr)->numSnpUI - snpNCntUI,
           snpNCntUI
         );

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u\t%u",
           (samSTPtr)->numInsUI,
           (samSTPtr)->numInsUI - insNCntUI,
           insNCntUI,
           (samSTPtr)->numDelUI
         );

         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u",
           (samSTPtr)->numMaskUI,
           (samSTPtr)->numMaskUI - maskNCntUI,
           maskNCntUI
         );
      } /*If: I am finding the anonymous base counts*/

      /**************************************************\
      * Fun15 Sec04 Sub04:
      *   - not printing anoynous bases do regular print
      \**************************************************/

      else
      { /*Else: I am not counting anonymous bases*/
         fprintf(
           (FILE *) (outFILE),
           "\t%u\t%u\t%u\t%u\t%u",
           (samSTPtr)->numMatchUI,
           (samSTPtr)->numSnpUI,
           (samSTPtr)->numInsUI,
           (samSTPtr)->numDelUI,
           (samSTPtr)->numMaskUI
         );
      } /*Else: I am not counting anonymous bases*/
      
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun15 Sec05:
      ^   - print out the q-scores
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      pqscore_fun15_sec05:;

      fprintf(
        (FILE *) (outFILE),
        "\t%f\t%f\n",
        (samSTPtr)->meanQF,
        (samSTPtr)->medianQF
      );

   } /*If: This is not a comment*/
} /*pstats_samEntry*/

/*-------------------------------------------------------\
| Fun16: revCmp_samEntry
|   - reverse complements a sam file sequence entry
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to reverse complement
| Output:
|   - Modifies:
|     o samSTPtr (includes flag) to be reverse complment
\-------------------------------------------------------*/
void
revCmp_samEntry(
   struct samEntry *samSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun16 TOC:
   '   - reverse complements a sam file sequence entry
   '   o fun16 sec01:
   '     - variable declarations
   '   o fun16 sec02:
   '     - check if have sequence and qscore entries
   '   o fun16 sec03:
   '     - reverse complement sequence and q-scores
   '   o fun16 sec04:
   '     - reverse the cigar and flag
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar *startStr = 0;
   schar *endStr = 0;

   schar *qStartStr = 0;
   schar *qEndStr = 0;

   sint *startSIPtr;
   sint *endSIPtr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec02:
   ^   - check if have sequence and qscore entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
         ! samSTPtr->seqStr
      || (
            samSTPtr->seqStr[0] == '*'
            && samSTPtr->seqStr[1] < 32
         )
   ) ; /*no sequence entry*/

   else
   { /*Else: have sequence entry*/
      startStr = samSTPtr->seqStr;
      endStr = startStr + samSTPtr->readLenUI;
   } /*Else: have sequence entry*/

   if(
         ! samSTPtr->qStr
      || (
            samSTPtr->qStr[0] == '*'
            && samSTPtr->qStr[1] < 32
         )
   ) ;

   else
   { /*Else: have q-score entry*/
      qStartStr = samSTPtr->qStr;
      qEndStr = qStartStr + samSTPtr->readLenUI;
   } /*Else: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec03:
   ^   - reverse complement sequence and q-scores
   ^   o fun16 sec03 sub01:
   ^     - reverse q-score entries + start loop
   ^   o fun16 sec03 sub02:
   ^     - reverse complement sequence entry
   ^   o fun16 sec03 sub03:
   ^     - make sure last sequence base is complement
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun16 Sec03 Sub01:
   *   - reverse q-score entries + start loop
   \*****************************************************/

   while(startStr < endStr)
   { /*Loop: reverse complement*/
      if(qStartStr)
      { /*If: have qscore entry*/
         *qStartStr ^= *qEndStr;
         *qEndStr ^= *qStartStr;
         *qStartStr ^= *qEndStr;

         ++qStartStr;
         --qEndStr;
      } /*If: have qscore entry*/

      /**************************************************\
      * Fun16 Sec03 Sub02:
      *   - reverse complement sequence entry
      \**************************************************/

      switch(*startStr & (~32))
      { /*Switch: complement*/
          case 'A': *startStr = 'T'; break;
          case 'C': *startStr = 'G'; break;
          case 'G': *startStr = 'C'; break;
          case 'T': *startStr = 'A'; break;
          case 'U': *startStr = 'A'; break;
          case 'W': *startStr = 'W'; break;
          case 'S': *startStr = 'S'; break;
          case 'M': *startStr = 'K'; break; /*A/C*/
          case 'K': *startStr = 'M'; break; /*T/G*/
          case 'R': *startStr = 'Y'; break; /*A/G*/
          case 'Y': *startStr = 'R'; break; /*T/C*/
          case 'B': *startStr = 'V'; break; /*C/G/T*/
          case 'D': *startStr = 'H'; break; /*G/T/A*/
          case 'H': *startStr = 'D'; break; /*C/A/T*/
          case 'V': *startStr = 'B'; break; /*A/C/G*/
          case 'N': *startStr = 'N'; break; /*A/C/G/T*/
          default: *startStr = 0;
      } /*Switch: complement*/

      switch(*endStr & (~32))
      { /*Switch: complement*/
          case 'A': *endStr = 'T'; break;
          case 'C': *endStr = 'G'; break;
          case 'G': *endStr = 'C'; break;
          case 'T': *endStr = 'A'; break;
          case 'U': *endStr = 'A'; break;
          case 'W': *endStr = 'W'; break;
          case 'S': *endStr = 'S'; break;
          case 'M': *endStr = 'K'; break; /*A/C*/
          case 'K': *endStr = 'M'; break; /*T/G*/
          case 'R': *endStr = 'Y'; break; /*A/G*/
          case 'Y': *endStr = 'R'; break; /*T/C*/
          case 'B': *endStr = 'V'; break; /*C/G/T*/
          case 'D': *endStr = 'H'; break; /*G/T/A*/
          case 'H': *endStr = 'D'; break; /*C/A/T*/
          case 'V': *endStr = 'B'; break; /*A/C/G*/
          case 'N': *endStr = 'N'; break; /*A/C/G/T*/
          default: *endStr = 0;
      } /*Switch: complement*/

      *startStr ^= *endStr;
      *endStr ^= *startStr;
      *startStr ^= *endStr;

      ++startStr;
      --endStr;
   } /*Loop: reverse complement*/

   /*****************************************************\
   * Fun16 Sec03 Sub03:
   *   - make sure last sequence base is complement
   \*****************************************************/

   if(
         startStr
      && startStr == endStr
   ){ /*If: have middle base*/
      switch(*startStr & (~32))
      { /*Switch: complement*/
          case 'A': *startStr = 'T'; break;
          case 'C': *startStr = 'G'; break;
          case 'G': *startStr = 'C'; break;
          case 'T': *startStr = 'A'; break;
          case 'U': *startStr = 'A'; break;
          case 'W': *startStr = 'W'; break;
          case 'S': *startStr = 'S'; break;
          case 'M': *startStr = 'K'; break; /*A/C*/
          case 'K': *startStr = 'M'; break; /*T/G*/
          case 'R': *startStr = 'Y'; break; /*A/G*/
          case 'Y': *startStr = 'R'; break; /*T/C*/
          case 'B': *startStr = 'V'; break; /*C/G/T*/
          case 'D': *startStr = 'H'; break; /*G/T/A*/
          case 'H': *startStr = 'D'; break; /*C/A/T*/
          case 'V': *startStr = 'B'; break; /*A/C/G*/
          case 'N': *startStr = 'N'; break; /*A/C/G/T*/
          default: *startStr = 0;
      } /*Switch: complement*/
   } /*If: have middle base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun16 Sec04:
   ^   - reverse the cigar and flag
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(samSTPtr->flagUS & 4) ;
      /*unmapped read; do nothing*/
   else if(samSTPtr->flagUS & 16)
      samSTPtr->flagUS &= ~16;
   else
      samSTPtr->flagUS |= 16;

   if(
         samSTPtr->cigTypeStr
      && samSTPtr->cigTypeStr[0] != '*'
   ){ /*If: have cigar to flip*/
      startStr = samSTPtr->cigTypeStr;
      endStr = startStr + samSTPtr->lenCigUI;

      startSIPtr = samSTPtr->cigArySI;
      endSIPtr = startSIPtr + samSTPtr->lenCigUI;

      while(startStr > endStr)
      { /*Loop: reverse cigar*/
         *startStr ^= *endStr;
         *endStr ^= *startStr;
         *startStr ^= *endStr;

         ++startStr;
         --endStr;

         *startSIPtr ^= *endSIPtr;
         *endSIPtr ^= *startSIPtr;
         *startSIPtr ^= *endSIPtr;

         ++startSIPtr;
         --endSIPtr;
      } /*Loop: reverse cigar*/
   } /*If: have cigar to flip*/
} /*revCmp_samEntry*/

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
