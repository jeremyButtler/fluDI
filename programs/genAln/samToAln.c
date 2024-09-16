/*#######################################################\
# Name: samToAln
#   - has functions to convert a sam entry and the
#     reference sequence to an alignment
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o fun01: psamPg_samToAln
'     - prints sam file program entries in the header
'   o fun02: phead_samToAln
'     - prints expanded cigar format header for alignment
'   o fun03: startLine_samToAln
'     - setup start of line in buffer
'   o fun04: paln_samToAln
'     - prints an alignment in expanded cigar format
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

#include <stdio.h>

#include "../genLib/samEntry.h"
#include "../genLib/seqST.h"
#include "../genLib/numToStr.h"
#include "../genLib/base10str.h"

#include "alnSet.h"

/*no .h files*/
#include "../genLib/dataTypeShortHand.h"
#include "alnDefs.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .c  #include "../generalLib/ulCp.h"
!   o .h  #include "../generalLib/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_maxDig_samToAln 9
   /*maximum digits in coordinates*/

/*-------------------------------------------------------\
| Fun01: psamPg_samToAln
|   - prints sam file program entries in the header
| Input:
|   - samSTPtr:
|     o pointer to a sam entry struct with first header
|   - buffStrPtr:
|     o pointer to buffer to use in read sam file headers
|   - lenBuffULPtr:
|     o pointer to unsigned long with length of buffStrPtr
|   - samFILE:
|     o sam file to get header from
|   - outFILE:
|     o file to print header to
| Output:
|   - Modifies:
|     o samSTPtr to have first entry after header
|     o buffStrPtr to be larger (if needed)
|     o lenBuffULPtr to have current buffStrPtr size
|   - Prints:
|     o aln master header to outFILe
|   - Returns:
|     o 0 for no errors
|     o def_EOF_samEntry for eof
|     o def_memErr_samEntry for memory errors
\-------------------------------------------------------*/
signed char
psamPg_samToAln(
   struct samEntry *samSTPtr,
   signed char **buffStrPtr,
   unsigned long *lenBuffULPtr,
   void *samFILE,
   void *outFILE
){
   schar errSC = 0;

   while(! errSC)
   { /*Loop: read in sam file entries*/
      if(samSTPtr->extraStr[0] != '@')
         break; /*finished*/

      if(
            samSTPtr->extraStr[1] == 'P'
         && samSTPtr->extraStr[2] == 'G'
      ){ /*If: is the program entry*/
         samSTPtr->extraStr[0] = '#';

         fprintf(
            (FILE *) outFILE,
            "%s\n",
            samSTPtr->extraStr
         );
      } /*If: is the program entry*/

      errSC =
         get_samEntry(
            samSTPtr,
            buffStrPtr,
            lenBuffULPtr,
            samFILE
         ); /*get next samfile entry*/
   } /*Loop: read in sam file entries*/

   if(errSC == def_EOF_samEntry)
      goto eof_fun02;

   if(errSC == def_memErr_samEntry)
      goto memErr_fun02;

   errSC = 0;
   goto cleanUp_fun02;

   memErr_fun02:;
   errSC = def_memErr_samEntry;
   goto cleanUp_fun02;

   eof_fun02:;
   errSC = def_EOF_samEntry;
   goto cleanUp_fun02;

   cleanUp_fun02:;
   return errSC;
} /*psamPg_samToAln*/

/*-------------------------------------------------------\
| Fun02: phead_samToAln
|   - prints expanded cigar format header for alignment
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct with aligned sequence
|   - refIdStr:
|     o c-string with reference name
|   - scoreSL:
|     o score for the alignment
|   - numAnonSL:
|     o number of anonymous matches
|   - outFILE:
|     o pointer to file to print header to
| Outputs:
|   - Prints:
|     o header for expanded cigar format to outFILE
\-------------------------------------------------------*/
void
phead_samToAln(
   struct samEntry *samSTPtr,
   signed char *refIdStr,
   signed long scoreSL,
   signed long numAnonSL,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints expanded cigar format header for alignment
   '   o fun02 sec01:
   '     - print out start of header
   '   o fun02 sec02:
   '     - print query information (id/coordinates)
   '   o fun02 sec03:
   '     - print refernece information (id/coordinates)
   '   o fun02 sec04:
   '     - print stats
   '   o fun02 sec05:
   '     - print legend
   '   o fun02 sec06:
   '     - print end of header
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - print out start of header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "#############################"
   ); /*first half of top header*/

   fprintf(
      (FILE *) outFILE,
      "#############################\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - print query information (id/coordinates)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "# qry: %s\n",
      samSTPtr->qryIdStr
   );


   if(
         samSTPtr->cigTypeStr[0] == 'S'
      || samSTPtr->cigTypeStr[0] == 'H'
   ){ /*If; masking at start*/
      fprintf(
         (FILE *) outFILE,
         "#    aligned: %u to",
         samSTPtr->cigArySI[0] + 1
      );
   } /*If: masking at start*/

   else
   { /*Else: no masking at start*/
      fprintf(
         (FILE *) outFILE,
         "#    aligned: 1 to"
      );
   } /*Else: no masking at start*/

   if(
         samSTPtr->cigTypeStr[samSTPtr->lenCigUI - 1]=='S'
      || samSTPtr->cigTypeStr[samSTPtr->lenCigUI - 1]=='H'
   ){ /*If; masking at end*/
      fprintf(
         (FILE *) outFILE,
         " %u\n",
         samSTPtr->readLenUI
          - samSTPtr->cigArySI[samSTPtr->lenCigUI - 1]
      );
   } /*If; masking at end*/

   else
   { /*Else: no masking at end*/
      fprintf(
         (FILE *) outFILE,
         " %u\n",
         samSTPtr->readLenUI
      );
   } /*Else: no masking at end*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - print refernece information (id/coordinates)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "# ref: %s\n",
      refIdStr
   );

   fprintf(
      (FILE *) outFILE,
      "#    aligned: %u to %u\n#\n",
      samSTPtr->refStartUI + 1,
      samSTPtr->refEndUI + 1
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^   - print stats
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "# matches: %u\n",
      samSTPtr->numMatchUI
   );

   fprintf(
      (FILE *) outFILE,
      "# snps: %u\n",
      samSTPtr->numSnpUI
   );

   fprintf(
      (FILE *) outFILE,
      "# insertions: %u\n",
      samSTPtr->numInsUI
   );

   fprintf(
      (FILE *) outFILE,
      "# deletions: %u\n",
      samSTPtr->numDelUI
   );

   fprintf(
      (FILE *) outFILE,
      "# masked: %u\n",
      samSTPtr->numMaskUI
   );

   fprintf(
      (FILE *) outFILE,
      "# anonymous: %li\n",
      numAnonSL
   );

   fprintf(
      (FILE *) outFILE,
      "# edit distance: %u\n",
       samSTPtr->numSnpUI
         + samSTPtr->numInsUI
         + samSTPtr->numDelUI
   );

   fprintf(
      (FILE *) outFILE,
      "# score: %li\n",
      scoreSL
   );

   fprintf(
      (FILE *) outFILE,
      "# mean q-score: %f\n",
      samSTPtr->meanQF
   );

   fprintf(
      (FILE *) outFILE,
      "# median q-score: %f\n#\n",
      samSTPtr->meanQF
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec05:
   ^   - print legend
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "# Eqx legend:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "#    = is match\n"
   );

   fprintf(
      (FILE *) outFILE,
      "#    X is snp\n"
   );

   fprintf(
      (FILE *) outFILE,
      "#    I is insertion\n"
   );

   fprintf(
      (FILE *) outFILE,
      "#    D is deletion\n"
   );

   fprintf(
      (FILE *) outFILE,
      "#    S is soft mask\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec06:
   ^   - print end of header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "#############################"
   ); /*first half of bottom header*/

   fprintf(
      (FILE *) outFILE,
      "#############################\n\n"
   );
} /*phead_samToAln*/

/*-------------------------------------------------------\
| Fun03: startLine_samToAln
|   - setup start of line in buffer
| Input:
|   - pposBl:
|     o 1: add teh positoin to the buffer
|   - refBuffStr:
|     o buffer to set up for reference sequence
|   - qryBuffStr:
|     o buffer to set up for query sequence
|   - cigBuffStr:
|     o buffer to set up for expanded cigar entry
|   - refPosUL:
|     o reference position to add to buffer (i pposBl = 1)
|   - qryPosUL:
|     o query position to add to buffer (i pposBl = 1)
|   - lenLineUS:
|     o pointer to unsigned short to have new buff length
| Output:
|   - Modifies:
|     o refBuffStr to have the reference entry
|     o qryBuffStr to have the query entry
|     o cigHeapStr to have the expand cigar entry
|     o lenbuffUS to length of start of alignment
|     o lenLineUS to length of start of alignment
\-------------------------------------------------------*/
void
startLine_samToAln(
   signed char pposBl,        /*1: print positions*/
   signed char *refBuffStr,
   signed char *qryBuffStr,
   signed char *cigBuffStr,
   unsigned long qryPosUL,    /*qury positoin to add*/
   unsigned long refPosUL,    /*ref position to add*/
   unsigned long *lenLineUL   /*characters in buffer*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - add line start ids and numbers to buffer
   '   o fun03 sec01:
   '     - variable declarations
   '   o fun03 sec02:
   '     - copy entry name
   '   o fun03 sec03:
   '     - add in coordinates
   '   o fun03 sec04:
   '     - add in ending space
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned short numDigUS = 0;
   unsigned short digOnUS = 0;
   signed char refDigStr[8 + (def_maxDig_samToAln << 2)];
   signed char qryDigStr[8 + (def_maxDig_samToAln << 2)];

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - copy entry name
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *lenLineUL = 0;

   refBuffStr[*lenLineUL] = 'r';
   qryBuffStr[*lenLineUL] = 'q';
   cigBuffStr[*lenLineUL] = 'c';
   ++(*lenLineUL);

   refBuffStr[*lenLineUL] = 'e';
   qryBuffStr[*lenLineUL] = 'r';
   cigBuffStr[*lenLineUL] = 'i';
   ++(*lenLineUL);

   refBuffStr[*lenLineUL] = 'f';
   qryBuffStr[*lenLineUL] = 'y';
   cigBuffStr[*lenLineUL] = 'g';
   ++(*lenLineUL);

   refBuffStr[*lenLineUL] = ':';
   qryBuffStr[*lenLineUL] = ':';
   cigBuffStr[*lenLineUL] = ':';
   ++(*lenLineUL);

   refBuffStr[*lenLineUL] = ' ';
   qryBuffStr[*lenLineUL] = ' ';
   cigBuffStr[*lenLineUL] = ' ';
   ++(*lenLineUL);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - add in coordinates
   ^   o fun03 sec03 sub01:
   ^     - add blank spaces to cigar
   ^   o fun03 sec03 sub02:
   ^     - get reference start coordinate
   ^   o fun03 sec03 sub03:
   ^     - get query start coordinate
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - add blank spaces to cigar
   \*****************************************************/

   if(pposBl)
   { /*If: printing coordinates*/
      for(
          numDigUS = 0;
          numDigUS < def_maxDig_samToAln;
          ++numDigUS
      ) cigBuffStr[*lenLineUL + numDigUS] = ' ';

      /**************************************************\
      * Fun03 Sec03 Sub02:
      *   - get reference start coordinate
      \**************************************************/

      numDigUS =
         numToStr(
            refDigStr,
            refPosUL
         );

      for(
         digOnUS = 0;
         digOnUS + numDigUS < def_maxDig_samToAln;
         ++digOnUS
      ) refBuffStr[*lenLineUL + digOnUS] = ' ';

      numDigUS = digOnUS;
      digOnUS = 0;

      while(numDigUS < def_maxDig_samToAln)
      { /*Loop: add white space*/
         refBuffStr[*lenLineUL + numDigUS++] =
            refDigStr[digOnUS++];
      } /*Loop: add white space*/

      /***********************************************\
      * Fun03 Sec02 Sub02:
      *   - get query start coordinate
      \***********************************************/

      numDigUS =
         numToStr(
            qryDigStr,
            qryPosUL
         );

      for(
         digOnUS = 0;
         digOnUS + numDigUS < def_maxDig_samToAln;
         ++digOnUS
      ) qryBuffStr[*lenLineUL + digOnUS] = ' ';

      numDigUS = digOnUS;
      digOnUS = 0;

      while(numDigUS < def_maxDig_samToAln)
      { /*Loop: add white space*/
         qryBuffStr[*lenLineUL + numDigUS++] =
            qryDigStr[digOnUS++];
      } /*Loop: add white space*/

      *lenLineUL += def_maxDig_samToAln;
   } /*If: printing coordinates*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - add in ending space
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   refBuffStr[*lenLineUL] = ' ';
   qryBuffStr[*lenLineUL] = ' ';
   cigBuffStr[*lenLineUL] = ' ';

   ++(*lenLineUL);
} /*startLine_samToAln*/

/*-------------------------------------------------------\
| Fun04: paln_samToAln
|   - prints an alignment in expanded cigar format
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct with aligned sequence
|   - refSTPtr:
|     o pointer to seqST struct with reference sequence
|   - scoreSL:
|     o score for the alignment
|   - setSTPtr:
|     o pointer to alnSet struct with alignment setSTPtr
|   - outFILE:
|     o pointer to file to print alignment to
| Outputs:
|   - Prints:
|     o header for expanded cigar format to outFILE
\-------------------------------------------------------*/
signed char
paln_samToAln(
   struct samEntry *samSTPtr,
   struct seqST *refSTPtr,
   struct alnSet *setSTPtr,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun04 TOC:
   '   - prints an alignment in expanded cigar format
   '   o fun04 sec01:
   '     - variable declarations
   '   o fun04 sec02:
   '     - memory allocation
   '   o fun04 sec03:
   '     - print header
   '   o fun04 sec04:
   '     - print soft masked bases at start
   '   o fun04 sec05:
   '     - print alignment
   '   o fun04 sec06:
   '     - print soft masked bases at end
   '   o fun04 sec07:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;

   ulong lenBuffUL = 0;
   schar buffBl = 0; /*something is in the buffer*/

   /*buffers to rpint out*/
   schar *qryHeapStr = 0;
   schar *refHeapStr = 0;
   schar *cigHeapStr = 0;

   ulong ulNt = 0;        /*number nucleotides in buffer*/

   ulong ulRef = 0;       /*reference base on*/
   ulong ulQry = 0;       /*query base on*/
   ulong qryOffsetUL = 0; /*hard masking offset*/
   ulong ulCig = 0;       /*cigar entry on*/

   slong scoreSL = 0;     /*for getting alingment score*/
   slong numAnonSL = 0;   /*holds number anonymous bases*/

   /*for dealing with starting softmasking*/
   schar *tmpSeqStr = 0;
   schar *tmpMaskStr = 0;
   schar *tmpNoSeqStr = 0;
   ulong *tmpPosUL = 0;

   ulong diffUL = 0;      /*also for cigar entry length*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec02:
   ^   - memory allocation
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(setSTPtr->lineWrapUI > 0)
      lenBuffUL = setSTPtr->lineWrapUI;

   else
      lenBuffUL= samSTPtr->readLenUI + refSTPtr->lenSeqUL;

   qryHeapStr =
      calloc(
         lenBuffUL + 9,
         sizeof(schar)
      );

   if(! qryHeapStr)
      goto memErr_fun04_sec07;


   refHeapStr =
      calloc(
         lenBuffUL + 9,
         sizeof(schar)
      );

   if(! refHeapStr)
      goto memErr_fun04_sec07;


   cigHeapStr =
      calloc(
         lenBuffUL + 9,
         sizeof(schar)
      );

   if(! cigHeapStr)
      goto memErr_fun04_sec07;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec03:
   ^   - print header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ulQry = 0;

   while(samSTPtr->extraStr[ulQry] != '\0')
   { /*Loop: get alignment score*/
      if(samSTPtr->extraStr[ulQry] == 'A')
      { /*If: could be alignment score*/
         ++ulQry;

         if(samSTPtr->extraStr[ulQry++] != 'S')
            goto nextEntry_fun04_sec03;

         if(samSTPtr->extraStr[ulQry++] != ':')
            goto nextEntry_fun04_sec03;

         if(samSTPtr->extraStr[ulQry++] != 'i')
            goto nextEntry_fun04_sec03;

         if(samSTPtr->extraStr[ulQry++] != ':')
            goto nextEntry_fun04_sec03;

         ulQry +=
            strToSL_base10str(
               &samSTPtr->extraStr[ulQry],
               &scoreSL
            );

         if(samSTPtr->extraStr[ulQry] > 32)
            scoreSL = -1; /*invalid score*/
      } /*If: could be alignment score*/

      
      if(samSTPtr->extraStr[ulQry] == 'n')
      { /*If: could be anonymous matching base count*/
         ++ulQry;

         if(samSTPtr->extraStr[ulQry++] != 'n')
            goto nextEntry_fun04_sec03;

         if(samSTPtr->extraStr[ulQry++] != ':')
            goto nextEntry_fun04_sec03;

         if(samSTPtr->extraStr[ulQry++] != 'i')
            goto nextEntry_fun04_sec03;

         if(samSTPtr->extraStr[ulQry++] != ':')
            goto nextEntry_fun04_sec03;

         ulQry +=
            strToSL_base10str(
               &samSTPtr->extraStr[ulQry],
               &numAnonSL
            );

         if(samSTPtr->extraStr[ulQry] > 32)
            numAnonSL = -1; /*invalid score*/
      } /*If: could be anonymous matching base count*/

      nextEntry_fun04_sec03:;
      while(samSTPtr->extraStr[ulQry++] > 32) ;
   } /*Loop: get alignment score*/

   phead_samToAln(
      samSTPtr,
      refSTPtr->idStr,
      scoreSL,
      numAnonSL,
      (FILE *) outFILE
   ); /*print out header for aligment*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec04:
   ^   - print soft masked bases at start
   ^   o fun04 sec04 sub01:
   ^     - find number of softmasked bases
   ^   o fun04 sec04 sub02:
   ^     - find if query/reference has more masked bases
   ^   o fun04 sec04 sub03:
   ^     - print extra soft masked bases
   ^   o fun04 sec04 sub04:
   ^     - print remaining masked bases for query/ref
   ^   o fun04 sec04 sub05:
   ^     - not printing masked bases; move to aligned part
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun04 Sec04 Sub01:
   *   - find number of softmasked bases
   \*****************************************************/

   ulQry = 0;
   ulRef = samSTPtr->refStartUI;

   if(samSTPtr->cigTypeStr[0] == 'S')
      ulQry = samSTPtr->cigArySI[0];

   if(samSTPtr->cigTypeStr[0] == 'H')
      qryOffsetUL = samSTPtr->cigArySI[0];

   /*****************************************************\
   * Fun04 Sec04 Sub02:
   *   - find if query or reference has more masked bases
   \*****************************************************/

   if(setSTPtr->pFullAlnBl)
   { /*If: printing masked bases*/
      startLine_samToAln(
         setSTPtr->pBasePosBl,
         refHeapStr,
         qryHeapStr,
         cigHeapStr,
         ulQry + qryOffsetUL, /*query*/
         1,               /*reference*/
         &ulNt
      ); /*set up buffers*/

      if(ulRef > ulQry)
      { /*If: reference has more soft masked bases*/
         diffUL = ulRef - ulQry;

         tmpSeqStr = refSTPtr->seqStr;
         tmpMaskStr = refHeapStr;
         tmpNoSeqStr = qryHeapStr;
         tmpPosUL = &ulRef;
      } /*If: reference has more soft masked bases*/

      else
      { /*Else: query has more soft masked bases*/
         diffUL = ulQry - ulRef;

         tmpSeqStr = samSTPtr->seqStr;
         tmpMaskStr = qryHeapStr;
         tmpNoSeqStr = refHeapStr;
         tmpPosUL = &ulQry;
      } /*Else: query has more soft masked bases*/

      ulQry = 0;
      ulRef = 0;

      /**************************************************\
      * Fun04 Sec04 Sub03:
      *   - print extra soft masked bases
      \**************************************************/

      while(*tmpPosUL < diffUL)
      { /*Loop: add soft masked bases in*/
         tmpMaskStr[ulNt] = tmpSeqStr[*tmpPosUL];

         tmpNoSeqStr[ulNt] = '-';
         cigHeapStr[ulNt] = 'S';

         ++(*tmpPosUL);
         ++ulNt;

         if(ulNt >= lenBuffUL)
         { /*If: buffer is full (print line)*/
            refHeapStr[ulNt] = '\0';
            qryHeapStr[ulNt] = '\0';
            cigHeapStr[ulNt] = '\0';

            fprintf(
               (FILE *) outFILE,
               "%s\n%s\n%s\n\n",
               refHeapStr,
               qryHeapStr,
               cigHeapStr
            );

            startLine_samToAln(
               setSTPtr->pBasePosBl,
               refHeapStr,
               qryHeapStr,
               cigHeapStr,
               ulQry + qryOffsetUL + 1,
               ulRef + 1,
               &ulNt
            ); /*set up buffers*/
         } /*If: buffer is full (print line)*/
      } /*Loop: add soft masked bases in*/

      /**************************************************\
      * Fun04 Sec04 Sub04:
      *   - print remaining masked bases for query/ref
      \**************************************************/

      if(samSTPtr->cigTypeStr[0] == 'S')
      { /*If: query has softmasking*/
         while(ulQry < (ulong) samSTPtr->cigArySI[0])
         { /*Loop: print out final softmasked bases*/
            qryHeapStr[ulNt] = samSTPtr->seqStr[ulQry];
            refHeapStr[ulNt] = refSTPtr->seqStr[ulRef];
            cigHeapStr[ulNt] = 'S';

            ++ulNt;
            ++ulQry;
            ++ulRef;

            if(ulNt >= lenBuffUL)
            { /*If: buffer is full (print line)*/
               refHeapStr[ulNt] = '\0';
               qryHeapStr[ulNt] = '\0';
               cigHeapStr[ulNt] = '\0';

               fprintf(
                  (FILE *) outFILE,
                  "%s\n%s\n%s\n\n",
                  refHeapStr,
                  qryHeapStr,
                  cigHeapStr
               );

               startLine_samToAln(
                  setSTPtr->pBasePosBl,
                  refHeapStr,
                  qryHeapStr,
                  cigHeapStr,
                  ulQry + qryOffsetUL + 1,
                  ulRef + 1,
                  &ulNt
               ); /*set up buffers*/
            } /*If: buffer is full (print line)*/
         } /*Loop: print out final softmasked bases*/
      } /*If: query has softmasking*/
   } /*If: printing masked bases*/

   /*****************************************************\
   * Fun04 Sec04 Sub05:
   *   - not printing masked bases; move to aligned part
   \*****************************************************/

   else
   { /*Else: skip soft masked bases*/
      startLine_samToAln(
         setSTPtr->pBasePosBl,
         refHeapStr,
         qryHeapStr,
         cigHeapStr,
         ulQry + qryOffsetUL + 1,
         ulRef + 1,
         &ulNt
      ); /*set up buffers*/
   } /*Else: skip soft masked bases*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec05:
   ^   - print alignment
   ^   o fun04 sec05 sub01:
   ^     - move past mask and start loop
   ^   o fun04 sec05 sub02:
   ^     - add nucleotide to alignment
   ^   o fun04 sec05 sub03:
   ^     - check if need to print entry out
   ^   o fun04 sec05 sub04:
   ^     - move to next nucleotide/cigar entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun04 Sec05 Sub01:
   *   - move past mask and start loop
   \*****************************************************/

   ulCig = (samSTPtr->cigTypeStr[0] == 'S');
   diffUL = (ulong) samSTPtr->cigArySI[ulCig];

   while(ulCig < samSTPtr->lenCigUI)
   { /*Loop: print out alignment*/

      /**************************************************\
      * Fun04 Sec05 Sub02:
      *   - add nucleotide to alignment
      *   o fun04 sec05 sub02 cat01:
      *     - assign base with math/mismatch flag 'M'
      *   o fun04 sec05 sub02 cat02:
      *     - assign base with match flag '='
      *   o fun04 sec05 sub02 cat03:
      *     - assign base with SNP flag 'X'
      *   o fun04 sec05 sub02 cat04:
      *     - assign base with deletion flag 'D'
      *   o fun04 sec05 sub02 cat05:
      *     - assign base with insertion flag 'I'
      *   o fun04 sec05 sub02 cat06:
      *     - stop at soft mask flag
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun04 Sec05 Sub02 Cat01:
      +   - assign base with math/mismatch flag 'M'
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      switch(samSTPtr->cigTypeStr[ulCig])
      { /*Switch: find cigar entry*/
         case 'M':
         /*Case: match or mismatch*/
            if(
               getMatch_alnSet(
                  refSTPtr->seqStr[ulRef]
                     & def_ntToCode_alnSet,
                  samSTPtr->seqStr[ulQry]
                     & def_ntToCode_alnSet,
                  setSTPtr
               ) /*get if bases are a match*/
            ) cigHeapStr[ulNt] = '=';

            else
               cigHeapStr[ulNt] = 'X';

            refHeapStr[ulNt] = refSTPtr->seqStr[ulRef];
            qryHeapStr[ulNt] = samSTPtr->seqStr[ulQry];

            ++ulRef;
            ++ulQry;
            --diffUL;

            break;
         /*Case: match or mismatch*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun04 Sec05 Sub02 Cat02:
         +   - assign base with match flag '='
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case '=':
         /*Case: match*/
            cigHeapStr[ulNt] = '=';
            refHeapStr[ulNt] = refSTPtr->seqStr[ulRef];
            qryHeapStr[ulNt] = samSTPtr->seqStr[ulQry];

            ++ulRef;
            ++ulQry;

            --diffUL;
            break;
         /*Case: match*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun04 Sec05 Sub02 Cat03:
         +   - assign base with SNP flag 'X'
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'X':
         /*Case: mismatch*/
            cigHeapStr[ulNt] = 'X';
            refHeapStr[ulNt] = refSTPtr->seqStr[ulRef];
            qryHeapStr[ulNt] = samSTPtr->seqStr[ulQry];

            ++ulRef;
            ++ulQry;
            --diffUL;
            break;
         /*Case: mismatch*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun04 Sec05 Sub02 Cat04:
         +   - assign base with deletion flag 'D'
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'D':
         /*Case: deletion*/
            cigHeapStr[ulNt] = 'D';
            refHeapStr[ulNt] = refSTPtr->seqStr[ulRef];
            qryHeapStr[ulNt] = '-';

            ++ulRef;
            --diffUL;
            break;
         /*Case: deletion*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun04 Sec05 Sub02 Cat05:
         +   - assign base with insertion flag 'I'
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'I':
         /*Case: insertion*/
            cigHeapStr[ulNt] = 'I';
            refHeapStr[ulNt] = '-';
            qryHeapStr[ulNt] = samSTPtr->seqStr[ulQry];

            ++ulQry;
            --diffUL;
            break;
         /*Case: insertion*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun04 Sec05 Sub02 Cat06:
         +   - stop at soft mask flag
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         case 'S': /*soft mask*/
            goto pendMask_fun04_sec06_sub01; /*finished*/

         case 'H': /*hard mask*/
            goto finalPrint_fun04_sec06_sub04;
            /*finished*/
      } /*Switch: find cigar entry*/

      /**************************************************\
      * Fun04 Sec05 Sub03:
      *   - check if need to print entry out
      \**************************************************/

      ++ulNt;
      buffBl = 1;

      if(ulNt >= lenBuffUL)
      { /*If: buffer is full (print line)*/
         buffBl = 0;

         refHeapStr[ulNt] = '\0';
         qryHeapStr[ulNt] = '\0';
         cigHeapStr[ulNt] = '\0';

         fprintf(
            (FILE *) outFILE,
            "%s\n%s\n%s\n\n",
            refHeapStr,
            qryHeapStr,
            cigHeapStr
         );

         startLine_samToAln(
            setSTPtr->pBasePosBl,
            refHeapStr,
            qryHeapStr,
            cigHeapStr,
            ulQry + qryOffsetUL + 1,
            ulRef + 1,
            &ulNt
         ); /*set up buffers*/
      } /*If: buffer is full (print line)*/

      /**************************************************\
      * Fun04 Sec05 Sub04:
      *   - move to next nucleotide/cigar entry
      \**************************************************/

      if(diffUL == 0)
      { /*If: moving to next cigar entry*/
         ++ulCig;
         diffUL = (ulong) samSTPtr->cigArySI[ulCig];
      } /*If: moving to next cigar entry*/
   } /*Loop: print out alignment*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec06:
   ^   - print soft masked bases at end
   ^   o fun04 sec06 sub01:
   ^     - check if printing ending masked bases
   ^   o fun04 sec06 sub02:
   ^     - print out the alignment ends
   ^   o fun04 sec06 sub03:
   ^     - print reference ending masked bases
   ^   o fun04 sec06 sub04:
   ^     - print out the alignment ends
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun04 Sec06 Sub01:
   *   - check if printing ending masked bases
   \*****************************************************/

   pendMask_fun04_sec06_sub01:;

   if(! setSTPtr->pFullAlnBl)
      goto finalPrint_fun04_sec06_sub04;
      /*skipping masked*/

   if(samSTPtr->cigTypeStr[ulCig] != 'S')
      goto prefMask_fun04_sec06_sub03;
      /*hard mask or no mask*/

   /*****************************************************\
   * Fun04 Sec06 Sub02:
   *   - print query ending masked bases
   \*****************************************************/

   diffUL = (ulong) samSTPtr->cigArySI[ulCig];

   while(diffUL > 0)
   { /*Loop: print masked bases*/
      if(refSTPtr->seqStr[ulRef] != '\0')
      { /*If: have reference bases to print*/
         refHeapStr[ulNt] = refSTPtr->seqStr[ulRef];
         ++ulRef;
      } /*If: have reference bases to print*/

      else
         refHeapStr[ulNt] = '-';

      qryHeapStr[ulNt] = samSTPtr->seqStr[ulQry];
      cigHeapStr[ulNt] = 'S';
      ++ulQry;
      ++ulNt;
      --diffUL;
      buffBl = 1;

      if(ulNt >= lenBuffUL)
      { /*If: buffer is full (print line)*/
         buffBl = 0;

         refHeapStr[ulNt] = '\0';
         qryHeapStr[ulNt] = '\0';
         cigHeapStr[ulNt] = '\0';

         fprintf(
            (FILE *) outFILE,
            "%s\n%s\n%s\n\n",
            refHeapStr,
            qryHeapStr,
            cigHeapStr
         );

         startLine_samToAln(
            setSTPtr->pBasePosBl,
            refHeapStr,
            qryHeapStr,
            cigHeapStr,
            ulQry + qryOffsetUL + 1,
            ulRef + 1,
            &ulNt
         ); /*set up buffers*/
      } /*If: buffer is full (print line)*/
   } /*Loop: print masked bases*/

   /*****************************************************\
   * Fun04 Sec06 Sub03:
   *   - print reference ending masked bases
   \*****************************************************/

   prefMask_fun04_sec06_sub03:;

   while(refSTPtr->seqStr[ulRef] != '\0')
   { /*Loop: get remaing reference masked bases*/
      refHeapStr[ulNt] = refSTPtr->seqStr[ulRef];
      ++ulRef;

      qryHeapStr[ulNt] = '-';
      cigHeapStr[ulNt] = 'S';

      ++ulNt;
      buffBl = 1;

      if(ulNt >= lenBuffUL)
      { /*If: buffer is full (print line)*/
         buffBl = 0;

         refHeapStr[ulNt] = '\0';
         qryHeapStr[ulNt] = '\0';
         cigHeapStr[ulNt] = '\0';

         fprintf(
            (FILE *) outFILE,
            "%s\n%s\n%s\n\n",
            refHeapStr,
            qryHeapStr,
            cigHeapStr
         );

         startLine_samToAln(
            setSTPtr->pBasePosBl,
            refHeapStr,
            qryHeapStr,
            cigHeapStr,
            ulQry + qryOffsetUL + 1,
            ulRef + 1,
            &ulNt
         ); /*set up buffers*/
      } /*If: buffer is full (print line)*/
   } /*Loop: get remaing reference masked bases*/

   /*****************************************************\
   * Fun04 Sec06 Sub04:
   *   - print out the alignment ends
   \*****************************************************/

   finalPrint_fun04_sec06_sub04:;

   if(refHeapStr[ulNt] == ' ')
      goto cleanUp_fun04_sec07; /*only starting entry*/

   if(! buffBl)
      goto cleanUp_fun04_sec07; /*nothing in buffer*/

   refHeapStr[ulNt] = '\0';
   qryHeapStr[ulNt] = '\0';
   cigHeapStr[ulNt] = '\0';

   fprintf(
      (FILE *) outFILE,
      "%s\n%s\n%s\n\n",
      refHeapStr,
      qryHeapStr,
      cigHeapStr
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec07:
   ^   - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "\n"
   ); /*print final separator*/

   errSC = 0;
   goto cleanUp_fun04_sec07;
   
   memErr_fun04_sec07:;
   errSC = def_memErr_alnDefs;
   goto cleanUp_fun04_sec07;

   cleanUp_fun04_sec07:;

   free(qryHeapStr);
   qryHeapStr = 0;

   free(refHeapStr);
   refHeapStr = 0;

   free(cigHeapStr);
   cigHeapStr = 0;

   return errSC;
} /*paln_samToAln*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However,for cases were the public domain is not
:   suitable,such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconvient / not possible,this code is under the
:   MIT license.
: 
: Public domain:
: 
: This is free and unencumbered software released into the
:   public domain.
: 
: Anyone is free to copy,modify,publish,use,compile,
:   sell,or distribute this software,either in source
:   code form or as a compiled binary,for any purpose,
:   commercial or non-commercial,and by any means.
: 
: In jurisdictions that recognize copyright laws,the
:   author or authors of this software dedicate any and
:   all copyright interest in the software to the public
:   domain. We make this dedication for the benefit of the
:   public at large and to the detriment of our heirs and
:   successors. We intend this dedication to be an overt
:   act of relinquishment in perpetuity of all present and
:   future rights to this software under copyright law.
: 
: THE SOFTWARE IS PROVIDED "AS IS\n",WITHOUT WARRANTY OF
:   ANY KIND,EXPRESS OR IMPLIED,INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO
:   EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM,
:   DAMAGES OR OTHER LIABILITY,WHETHER IN AN ACTION OF
:   CONTRACT,TORT OR OTHERWISE,ARISING FROM,OUT OF OR
:   IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
:   DEALINGS IN THE SOFTWARE.
: 
: For more information,please refer to
:   <https://unlicense.org>
: 
: MIT License:
: 
: Copyright (c) 2024 jeremyButtler
: 
: Permission is hereby granted,free of charge,to any
:   person obtaining a copy of this software and
:   associated documentation files (the "Software"),to
:   deal in the Software without restriction,including
:   without limitation the rights to use,copy,modify,
:   merge,publish,distribute,sublicense,and/or sell
:   copies of the Software,and to permit persons to whom
:   the Software is furnished to do so,subject to the
:   following conditions:
: 
: The above copyright notice and this permission notice
:   shall be included in all copies or substantial
:   portions of the Software.
: 
: THE SOFTWARE IS PROVIDED "AS IS\n",WITHOUT WARRANTY OF
:   ANY KIND,EXPRESS OR IMPLIED,INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM,DAMAGES OR OTHER LIABILITY,WHETHER IN
:   AN ACTION OF CONTRACT,TORT OR OTHERWISE,ARISING
:   FROM,OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
