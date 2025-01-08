/*########################################################
# Name: getDiDis
#   - finds diRNA and mvRNA sequences in reads
########################################################*/

/*-------------------------------------------------------\
' SOF: Start Of File
'   o header:
'     - included libraries and defined variables
'   o fun01: pversion_diIDs
'     - prints the version for primFind to the input file
'   o fun02: phelp_diIDs
'     - prints the help message for primFind to the input
'       file
'   o fun03: getSeq_diIDs
'     - gets a sequence from a sequence user input flag
'   o fun04: getInput_diIDs
'     - gets user input for primFind
'   o main:
'     - driver function
'   o license:
'     - licensing for this code (public domain / mit)
\-------------------------------------------------------*/

/*-------------------------------------------------------\
| Header:
|   - included libraries and defined variables
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../genLib/base10str.h"
#include "../genLib/charCp.h"
#include "../genLib/ulCp.h"
#include "../genBio/seqST.h"
#include "../genAln/alnSet.h"
#include "../genAln/kmerFind.h"

#include "fluST.h"

/*.h files only*/
#include "../genLib/dataTypeShortHand.h"
#include "../genBio/ntTo2Bit.h"

#include "../genAln/alnDefs.h"

#include "fluSeg.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../genLib/shellSort.h"
!   - .c  #include "../genAln/memwater.h"
!   - .c  #include "../genAln/indexToCoord.h"
!   - .h  #include "../genLib/genMath.h" (using .h macros)
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_year_diIDs 2024
#define def_month_diIDs 8
#define def_day_diIDs 29

#define def_faPrimFile_diIDs 1
#define def_tsvPrimFile_diIDs 2

#define def_minLen_diIDs 40
#define def_maxLen_diIDs 3000
#define def_minPercLen_diIDs 0.85f
#define def_maxPercLen_diIDs 1.1f /*110%*/

#define def_fastSearch_diIDs 1 /*1 = fast search*/

#define def_pDIRna_diIDs 1
#define def_pMVRna_diIDs 1
#define def_pVRna_diIDs 1
#define def_partSeq_diIDs 1

#define def_noAgree_diIDs 0
   /*primers in read not agreeing for segment*/

static signed char
   *forPrimStr_diIds =
      (signed char *) "AGCGAAAGCAGG";
static signed char
   *revPrimStr_diIds =
      (signed char *) "AGTAGAAACAAGG";

/*-------------------------------------------------------\
| Fun01: pversion_diIDs
|   - prints the version for primFind to the input file
| Input:
|   - outFILE:
|     o file to print version number to
| Output:
|   - Prints;
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_diIDs(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "diIDs version: %i-%02i-%02i\n",
       def_year_diIDs,
       def_month_diIDs,
       def_day_diIDs
   );
} /*pverson_diIDs*/

/*-------------------------------------------------------\
| Fun02: phelp_diIDs
|   - prints the help message for primFind to the input
|     file
| Input:
|   - pSegHelpBl:
|     o 1 print the segment help message
|     o 0 print the shorter help message
|   - outFILE:
|     o file to print the help message to
| Output:
|   - Prints;
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_diIDs(
   signed char pSegHelpBl,
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints help message
   '   o fun02 sec01:
   '     - print out usage
   '   o fun02 sec02:
   '     - print out input
   '   o fun02 sec03:
   '     - print output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - print out usage
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint siSeg = 0;

   fprintf(
      (FILE *) outFILE,
      "getDIIds -fq reads.fasta\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - Finds full length flu segments and checks if\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    segment is diRNA or vRNA (full length)\n"
   );

   if(pSegHelpBl)
       goto segInput_fun02_sec02_sub03_cat02;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - print out input
   ^   o fun02 sec02 sub01:
   ^      - input block header
   ^   o fun02 sec02 sub02:
   ^     - file IO
   ^   o fun02 sec02 sub03:
   ^     - changing sequence input
   ^   o fun02 sec02 sub04:
   ^     - filtering input
   ^   o fun02 sec02 sub05:
   ^     - kmer search options
   ^   o fun02 sec02 sub06:
   ^     - help message and version number options
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec02 Sub01:
   *   - input block header
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "Input:\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub02:
   *   - file IO
   *   o fun02 sec02 sub02 cat01:
   *     - read file (fastx)
   *   o fun02 sec02 sub02 cat02:
   *     - print di reads
   *   o fun02 sec02 sub02 cat03:
   *     - print full length reads
   *   o fun02 sec02 sub02 cat04:
   *     - print genomes with one primer id for segment
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat01:
   +   - read file (fastx)
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    fprintf(
      (FILE *) outFILE,
       "  -fq reads.fastq: [Required]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o fastq file with reads to search\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -fa reads.fasta: [Replaces -fq]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o fasta file with reads to search\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat02:
   +   - print di reads
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    if(def_pDIRna_diIDs)
       fprintf(
         (FILE *) outFILE,
          "  -di: [Yes]\n"
       );

    else
       fprintf(
         (FILE *) outFILE,
          "  -di: [No]\n"
       );

    fprintf(
      (FILE *) outFILE,
       "    o print out diRNA read ids\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o disable with \"-no-di\"\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat03:
   +   - print full length reads
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    if(def_pVRna_diIDs)
       fprintf(
         (FILE *) outFILE,
          "  -vrna: [Yes]\n"
       );

    else
       fprintf(
         (FILE *) outFILE,
          "  -vrna: [No]\n"
       );

    fprintf(
      (FILE *) outFILE,
       "    o print out read ids for vRNA (full length)\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o disable with \"-no-vrna\"\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub02 Cat04:
   +   - print genoems with on solid primer
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    if(def_partSeq_diIDs)
       fprintf(
         (FILE *) outFILE,
          "  -part: [Yes]\n"
       );

    else
       fprintf(
         (FILE *) outFILE,
          "  -part: [No]\n"
       );

    fprintf(
      (FILE *) outFILE,
       "    o print out read ids with only one primer\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      having a segment id\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o disable with \"-no-part\"\n"
    );

   /*****************************************************\
   * Fun02 Sec02 Sub03:
   *   - changing sequence input
   *   o fun02 sec02 sub03 cat01:
   *     - primer sequence input
   *   o fun02 sec02 sub03 cat02:
   *     - segment input
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub03 Cat01:
   +   - primer sequence input
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    fprintf(
      (FILE *) outFILE,
       "  -prims %s,%s: [Optional]\n",
       forPrimStr_diIds,
       revPrimStr_diIds
    );

    fprintf(
      (FILE *) outFILE,
       "    o primer pair (forward,reverse) to find\n"
    );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub03 Cat02:
   +   - segment input
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    fprintf(
      (FILE *) outFILE,
       "  -segment forward-seq,reverse-seq: [Optional]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o sequences used to id a segment\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      HA would be -HA GG,GTGT\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -segment-len length: [Optional]\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o expected length of a segment\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      HA would be -HA 1778\n"
    );

   if(pSegHelpBl)
   { /*If: printing segment help message*/
      segInput_fun02_sec02_sub03_cat02:;

      for(
         siSeg = 0;
         siSeg < def_NSNum_fluSeg + 1;
         ++siSeg
      ){ /*Loop: print out segment sequences and length*/
         fprintf(
            (FILE *) outFILE,
            "    o -%s-len %i\t-%s %s,%s\n",
            segIdAryStr_fluSeg[siSeg],
            segLenArySS_fluSeg[siSeg],
            segIdAryStr_fluSeg[siSeg],
            forSeqAryStr_fluSeg[siSeg],
            revSeqAryStr_fluSeg[siSeg]
         );
      } /*Loop: print out segment sequences and length*/

      goto helpInput_fun02_sec02_sub06;
   } /*If: printing segment help message*/

   /*****************************************************\
   * Fun02 Sec02 Sub04:
   *   - filtering input
   \*****************************************************/

    fprintf(
      (FILE *) outFILE,
       "  -min-len %i: [Optional]\n",
       def_minLen_diIDs
    );

    fprintf(
      (FILE *) outFILE,
       "    o minimum read length to keep a read\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -max-len %i: [Optional]\n",
       def_maxLen_diIDs
    );

    fprintf(
      (FILE *) outFILE,
       "    o maximum read length to keep a read\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -min-perc-len %0.2f: [Optional]\n",
       def_minPercLen_diIDs
    );

    fprintf(
      (FILE *) outFILE,
       "    o minimum percent length to keep a read id\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o (length between primers) / (read length)\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -max-perc-len %0.2f: [Optional]\n",
       def_maxPercLen_diIDs
    );

    fprintf(
      (FILE *) outFILE,
       "    o maximum percent length to keep a read id\n"
    );

    fprintf(
      (FILE *) outFILE,
       "    o (mapped read length) / (expected length)\n"
    );

    fprintf(
      (FILE *) outFILE,
       "  -min-perc-score %0.2f: [Optional]\n",
       def_minPercScore_kmerFind
    );

    fprintf(
      (FILE *) outFILE,
       "    o minimum percent score (alignment) to\n"
    );

    fprintf(
      (FILE *) outFILE,
       "      count a primer as mapped\n"
    );

   /*****************************************************\
   * Fun02 Sec02 Sub05:
   *   - kmer search options
   *   o fun02 sec02 sub05 cat01:
   *     - fast/slow options
   *   o fun02 sec02 sub05 cat02:
   *     - kmer length options
   *   o fun02 sec02 sub05 cat03:
   *     - min percent kmers setting
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat01:
   +   - fast/slow options
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_fastSearch_diIDs)
      fprintf(
        (FILE *) outFILE,
         "  -fast: [True]\n"
      );

   else
      fprintf(
        (FILE *) outFILE,
         "  -fast: [False]\n"
      );

    fprintf(
      (FILE *) outFILE,
       "    o do faster, but less percise kmer search\n"
    );

   fprintf(
     (FILE *) outFILE,
      "    o use \"-slow\" to search by Waterman\n"
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat02:
   +   - kmer length options
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
     (FILE *) outFILE,
      "  -len-kmer %i: [Optional]\n",
      def_lenKmer_kmerFind
   );

   fprintf(
     (FILE *) outFILE,
      "    o kmer size for -fast; bigger is faster, but\n"
   );

   fprintf(
     (FILE *) outFILE,
      "      also less sensitive (misses more)\n"
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat03:
   +   - min percent kmers setting
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
     (FILE *) outFILE,
      "  -min-perc-kmer %0.2f: [Optional]\n",
      def_minKmerPerc_kmerFind
   );

   fprintf(
     (FILE *) outFILE,
      "    o minimum percent of total kmers needed to\n"
   );

   fprintf(
     (FILE *) outFILE,
      "      do a Waterman on a window\n"
   );

   fprintf(
     (FILE *) outFILE,
      "    o higher is faster, but also less sensitive\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub06:
   *   - help message and version number options
   \*****************************************************/

   helpInput_fun02_sec02_sub06:;

   if(! pSegHelpBl)
   { /*If: printing the normal shorter help message*/
      fprintf(
        (FILE *) outFILE,
         "  -h: print this help message and exit\n"
      );

      fprintf(
        (FILE *) outFILE,
         "  -h-seg: print segment defaults and exit\n"
      );
   } /*If: printing the normal shorter help message*/

   else
   { /*Else: this is the segment help message*/
      fprintf(
        (FILE *) outFILE,
         "  -h: print non-segmet help message and exit \n"
      );

      fprintf(
        (FILE *) outFILE,
         "  -h-seg: print this help message and exit\n"
      );
   } /*Else: this is the segment help message*/

   fprintf(
     (FILE *) outFILE,
      "  -v: print version number and exit\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - print output block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "Output:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - prints reads with detected primers to stdout\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o column 1 is read ID; see output for rest\n"
   );
} /*phelp_diIDs*/

/*-------------------------------------------------------\
| Fun03: getSeq_diIDs
|   - gets a sequence from a sequence user input flag
| Input:
|   - firstSeqStr:
|     o c-string to hold the first sequence read in
|   - secSeqStr:
|     o c-string to hold the second sequence read in
|   - inStr:
|     o c-string with user input
| Output:
|   - Modifies:
|     o firstSeqStr to point to frist sequence in user
|       input
|     o secSeqStr to point to second sequence in user
|       input
|     o inStr to have a null separating the two sequences
|   - Returns:
|     o 0: for no problems
|     o 1: for only one sequence
|     o 2: for more than two sequences
\-------------------------------------------------------*/
signed char
getSeq_diIDs(
   signed char **firstSeqStr,
   signed char **secSeqStr,
   signed char *inStr
){
   *firstSeqStr = inStr; 

   while(
      ntTo2Bit[(uchar) *inStr++] != def_err3rdBit_ntTo2Bit
   ) ;

   --inStr;

   if(*inStr == '\0')
      return 1; /*only one sequence*/

   *inStr = '\0';
   ++inStr;

   *secSeqStr = inStr;

   while(
      ntTo2Bit[(uchar) *inStr++] != def_err3rdBit_ntTo2Bit
   ) ;

   --inStr;

   if(*inStr != '\0')
      return 2; /*more than two sequences*/

   return 0;
} /*getSeq_diIDs*/

/*-------------------------------------------------------\
| Fun04: getInput_diIDs
|   - gets user input for primFind
| Input:
|   - numArgsSI:
|     o number of arguments/parameters in argAryStr
|   - argAryStr:
|     o array of c-strings with user input
|   - fluSTPtr:
|     o pointer to a fluST structure with settings to
|       change
|   - readFileStr:
|     o pointer to c-string to point to file with reads
|   - fqBl:
|     o pointer to char to be set to
|       - 1 = readFileStr is a fastq file
|       - 0 = readFileStr is a fasta file
|   - forPrimStr:
|     o pointer to c-string to hold forward primer seq
|   - revPrimStr:
|     o pointer to c-string to hold reverse primer seq
|   - minLenUI:
|     o pointer to unsigned long to hold minimum primer
|       length
|   - maxlenUI:
|     o pointer to unsigned long to hold maximum primer
|       length
|   - minPercLenF:
|     o pointer to float to hold the minimum percent
|       length to keep a read
|   - maxPercLenF:
|     o pointer to float to hold the maximum percent
|       length (of expected segment length) to keep a read
|   - minPerScoreF:
|     o pointer to float to hold the minimum percent score
|       to count a primer as mapped (found by waterman)
|   - pDiRnaBl:
|     o pointer to set to 1 (print diRNA) or 0 (no print)
|   - pVRnaBl:
|     o pointer to set to 1 (print vRNA) or 0 (no print)
|   - pPartBl:
|     o pointer to set to 1 to print out segments with
|       only one primer having segment id
|     o pointer to set to 0 to not print out 
|   - fastBl:
|     o ponter to char to be set to
|       - 1 = use faster, but less senistive kmer method
|       - 0 = use slower, but more percise waterman method
|   - lenKmerUC:
|     o pointer to unsigned char to hold kmer length for
|       the faster kmer search
|   - minPercKmerF:
|     o pointer to float to hold min percentage of kmers
|       needed to do waterman alignment (kmer search only)
| Output:
|   - Modifies:
|     o all input, except numArgsSI and argAryStr to hold
|       the users input
|   - Prints:
|     o help message and version number to stdout if
|       requested
|     o prints any errors to stderr
|   - Returns:
|     o 0 for no errors
|     o 1 for printed help message or version number
|     o 2 for an error
\-------------------------------------------------------*/
signed char
getInput_diIDs(
   int numArgsSI,
   char *argAryStr[],
   struct fluST *fluSTPtr,
   signed char **readFileStr,
   signed char *fqBl,
   signed char **forPrimStr,
   signed char **revPrimStr,
   unsigned int *minLenUI,
   unsigned int *maxLenUI,
   float *minPercLenF,
   float *maxPercLenF,
   float *minPercScoreF,
   signed char *pDiRnaBl,
   signed char *pVRnaBl,
   signed char *pPartBl,
   signed char *fastBl,
   unsigned char *lenKmerUC,
   float *minPercKmerF
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun04 TOC:
   '   - gets user input for diIDs
   '   o fun04 sec01:
   '     - variable declerations
   '   o fun04 sec01:
   '     - get input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec01:
   ^   - variable declerations
   ^   fun04 sec01 sub01:
   ^     - declare variables
   ^   fun04 sec01 sub02:
   ^     - sert up flags for flu segments (pattern/length)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun04 Sec01 Sub01:
   *   - declare variables
   \*****************************************************/

   schar *tmpStr = 0;
   sint siArg = 1;
   schar errSC = 0;
   sint eqlSI = 0;


   schar segFlagsStr[(def_NSNum_fluSeg + 1) << 1][10];
   sint siSeg = 0;

   schar *forSeqStr = 0;
   schar *revSeqStr = 0;

   /*****************************************************\
   * Fun04 Sec01 Sub02:
   *   - sert up flags for flu segments (pattern/length)
   \*****************************************************/

   for(
      siSeg = 0;
      siSeg < def_NSNum_fluSeg + 1;
      siSeg += 2
   ){ /*Loop: build segment id flags*/
      /*set up sequence flag*/
      tmpStr = segFlagsStr[siSeg];

      *tmpStr++ = '-';

      tmpStr +=
         eql_charCp(
            tmpStr,
            segIdAryStr_fluSeg[siSeg >> 1],
            '\0'
         );

      /*set up length flag*/
      tmpStr = segFlagsStr[siSeg + 1];

      *tmpStr++ = '-';

      tmpStr +=
         cpDelim_charCp(
            tmpStr,
            segIdAryStr_fluSeg[siSeg >> 1],
            '\0'
         );

      tmpStr +=
         cpDelim_charCp(
            tmpStr,
            (schar *) "-len",
            '\0'
         );
   } /*Loop: build segment id flags*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun04 Sec02:
   ^   - get input
   ^   o fun04 sec02 sub01:
   ^     - check if have input and start loop
   ^   o fun04 sec02 sub02:
   ^     - file IO
   ^   o fun04 sec02 sub03:
   ^     - check if primer sequences
   ^   o fun04 sec02 sub04:
   ^     - filtering input (length/score)
   ^   o fun04 sec02 sub05:
   ^     - fast alignment input
   ^   o fun04 sec02 sub06:
   ^     - help messages
   ^   o fun04 sec02 sub07:
   ^     - check segment input or if invalid
   ^   o fun04 sec02 sub08:
   ^     - move to next argument
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun04 Sec02 Sub01:
   *   - check if have input and start loop
   \*****************************************************/

   if(numArgsSI < 2)
   { /*If: nothing was input; assume help message wanted*/
      phelp_diIDs(
         0,      /*print short help message*/
         stdout
      );

      return 1;
   } /*If: nothing was input; assume help message wanted*/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/

      /**************************************************\
      * Fun04 Sec02 Sub02:
      *   - file IO
      \**************************************************/

      if(
         ! eql_charCp(
            (schar *) "-fa",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*If: a fasta sequence file was input*/
         ++siArg; /*move to argument*/
         *fqBl = 0;
         *readFileStr = (schar *) argAryStr[siArg];
      } /*If: a fasta sequence file was input*/

      else if(
         ! eql_charCp(
            (schar *) "-fq",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: a fastq sequence file was input*/
         ++siArg; /*move to argument*/
         *fqBl = 1;
         *readFileStr = (schar *) argAryStr[siArg];
      } /*Else If: a fastq sequence file was input*/

      else if(
         ! eql_charCp(
            (schar *) "-di",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pDiRnaBl = 1;

      else if(
         ! eql_charCp(
            (schar *) "-no-di",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pDiRnaBl = 0;

      else if(
         ! eql_charCp(
            (schar *) "-vrna",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pVRnaBl = 1;

      else if(
         ! eql_charCp(
            (schar *) "-no-vrna",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pVRnaBl = 0;

      else if(
         ! eql_charCp(
            (schar *) "-part",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pPartBl = 1;

      else if(
         ! eql_charCp(
            (schar *) "-no-part",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pPartBl = 0;

      /**************************************************\
      * Fun04 Sec02 Sub03:
      *   - check if primer sequences
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-prims",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: primer sequences input*/
         ++siArg; /*move to argument*/

         errSC = 
            getSeq_diIDs(
               forPrimStr,
               revPrimStr,
               (schar *) argAryStr[siArg]
            ); /*set up pointers for primer sequence*/

         if(errSC)
         { /*If: there was a error*/
            if(errSC == 1)
               fprintf(
                  stderr,
                  "-prims %s has only one sequence\n",
                  argAryStr[siArg]
               );
            else
               fprintf(
                  stderr,
                  "-prims %s has 3 or more sequences\n",
                  argAryStr[siArg]
               );

            return 2;
         } /*If: there was a error*/
      } /*Else If: primer sequences input*/

      /**************************************************\
      * Fun04 Sec02 Sub04:
      *   - filtering input (length/score)
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-min-len",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: the min length was input*/
         ++siArg; /*move to argument*/
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToUI_base10str(
               tmpStr,
               minLenUI
            );

         if(*tmpStr != '\0')
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-min-len %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: the min length was input*/

      else if(
         ! eql_charCp(
            (schar *) "-max-len",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: the max length was input*/
         ++siArg; /*move to argument*/
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToUI_base10str(
               tmpStr,
               maxLenUI
            );

         if(*tmpStr != '\0')
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-max-len %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: the max length was input*/

      else if(
         ! eql_charCp(
            (schar *) "-min-perc-len",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: the min percent length was input*/
         ++siArg; /*move to argument*/

          *minPercLenF = atof(argAryStr[siArg]);

         if(minPercLenF == 0)
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-min-perc-len %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: the min percent length was input*/

      else if(
         ! eql_charCp(
            (schar *) "-max-perc-len",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: the max percent length was input*/
         ++siArg; /*move to argument*/

          *maxPercLenF = atof(argAryStr[siArg]);

         if(minPercLenF == 0)
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-max-perc-len %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: the max percent length was input*/

      else if(
          ! eql_charCp(
               (schar *) "-min-perc-score",
               (schar *) argAryStr[siArg],
               '\0'
            )
      ){ /*Else If: the min percent score was input*/
         ++siArg; /*move to argument*/
         *minPercScoreF = atof(argAryStr[siArg]);

         if(*minPercScoreF == 0)
         { /*If: invalid input*/
            fprintf(
               stderr,
               "-min-perc-score %s is non-numeric or 0\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: invalid input*/
      } /*Else If: the min percent score was input*/

      /**************************************************\
      * Fun04 Sec02 Sub05:
      *   - fast alignment input
      \**************************************************/

      /*check if using a fast (kmer) or slow alignment*/
      else if(
         ! eql_charCp(
            (schar *)  "-fast",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *fastBl = 1;

      else if(
         ! eql_charCp(
            (schar *)  "-slow",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *fastBl = 0;

      else if(
        ! eql_charCp(
           (schar *) "-len-kmer",
           (schar *) argAryStr[siArg],
           (schar) '\0'
         )
      ){ /*Else If: kmer length was input*/
         ++siArg; /*move to argument*/
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToUC_base10str(
               tmpStr,
               lenKmerUC
            );

         if(*tmpStr != '\0')
         { /*If: this is non-numeric*/
            fprintf(
               stderr,
               "-len-kmer %s is non-numeric\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: this is non-numeric*/
      } /*Else If: kmer length was input*/

      else if(
          ! eql_charCp(
               (schar *) "-min-perc-kmer",
               (schar *) argAryStr[siArg],
               (schar) '\0'
          )
      ){ /*Else If: min percent kmers was input*/
         ++siArg; /*move to argument*/
         *minPercKmerF = atof(argAryStr[siArg]);

         if(*minPercScoreF == 0)
         { /*If: invalid input*/
            fprintf(
               stderr,
               "-min-perc-kmer %s is non-numeric or 0\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: invalid input*/
      } /*Else If: min percent kmers was input*/

      /**************************************************\
      * Fun04 Sec02 Sub06:
      *   - help and version messages
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-h",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         phelp_diIDs(
            0,     /*print shorter help message*/
            stdout
         );

         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "--h",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         phelp_diIDs(
            0,     /*print shorter help message*/
            stdout
         );

         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "help",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         phelp_diIDs(
            0,     /*print shorter help message*/
            stdout
         );

         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "-help",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         phelp_diIDs(
            0,     /*print shorter help message*/
            stdout
         );

         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "--help",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         phelp_diIDs(
            0,     /*print shorter help message*/
            stdout
         );

         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "-h-seg",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         phelp_diIDs(
            1,     /*print shorter help message*/
            stdout
         );

         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "-v",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         pversion_diIDs(stdout);
         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "--v",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         pversion_diIDs(stdout);
         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "version",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         pversion_diIDs(stdout);
         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "-version",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         pversion_diIDs(stdout);
         return 1;
      } /*Else If: help message*/

      else if(
         ! eql_charCp(
            (schar *) "--version",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message*/
         pversion_diIDs(stdout);
         return 1;
      } /*Else If: help message*/

      /**************************************************\
      * Fun04 Sec02 Sub07:
      *   - check segment input or if invalid
      *   o fun04 sec02 sub07 cat01:
      *     - check if segment length (and start loop)
      *   o fun04 sec02 sub07 cat02:
      *     - check if segment sequence
      *   o fun04 sec02 sub07 cat03:
      *     - check if valid input
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun04 Sec02 Sub07 Cat01:
      +   - check if segment length (and start loop)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else
      { /*Else: segment input or not recognized*/

         for(
            siSeg = 0;
            siSeg < def_NSNum_fluSeg + 1;
            siSeg += 2
         ){ /*Loop: check if segment input*/
            eqlSI =
               eql_charCp(
                  segFlagsStr[siSeg + 1],
                  (schar *) argAryStr[siArg],
                  '\0'
               ); /*check if is segment length*/

            if(! eqlSI)
            { /*If: input is length of segment*/
               ++siArg;
               tmpStr = (schar *) argAryStr[siArg];

               tmpStr +=
                  strToSI_base10str(
                     tmpStr,
                     &fluSTPtr->lenSegArySI[siSeg >> 1]
                  );

               if(*tmpStr != '\0')
               { /*If: I had a non-numeric entry*/
                  fprintf(
                     stderr,
                     "%s %s is non-numeric\n",
                     segFlagsStr[siSeg + 1],
                     argAryStr[siArg]
                  );

                  return 2;
               } /*If: I had a non-numeric entry*/

               break;
            } /*If: input is length of segment*/

            eqlSI =
               eql_charCp(
                  segFlagsStr[siSeg],
                  (schar *) argAryStr[siArg],
                  '\0'
               ); /*check if is segement sequence*/

            /*+++++++++++++++++++++++++++++++++++++++++++\
            + Fun04 Sec02 Sub07 Cat02:
            +   - check if segment sequence
            \+++++++++++++++++++++++++++++++++++++++++++*/

            if(! eqlSI)
            { /*Else If: segment sequences input*/
               ++siArg; /*move to argument*/

               errSC = 
                  getSeq_diIDs(
                     &forSeqStr,
                     &revSeqStr,
                     (schar *) argAryStr[siArg]
                  ); /*set up pointers for sequences*/


               errSC =
                  rmSegSeqFrom_fluST(
                     fluSTPtr,
                     forSeqAryStr_fluSeg[siSeg >> 1],
                     0,
                     segLenArySS_fluSeg[siSeg >> 1]
                  );

               errSC =
                  addSegTo_fluST(
                     fluSTPtr,
                     forSeqStr,
                     segLenArySS_fluSeg[siSeg >> 1],
                     0,
                     siSeg >> 1
                  );

               if(errSC)
                  goto segErr_fun04_sec0x_sub0y;

               errSC =
                  rmSegSeqFrom_fluST(
                     fluSTPtr,
                     revSeqAryStr_fluSeg[siSeg >> 1],
                     1,
                     segLenArySS_fluSeg[siSeg >> 1]
                  );

               errSC =
                  addSegTo_fluST(
                     fluSTPtr,
                     revSeqStr,
                     segLenArySS_fluSeg[siSeg >> 1],
                     1,
                     siSeg >> 1
                  );

               if(errSC)
               { /*If: there was a error*/
                  segErr_fun04_sec0x_sub0y:;

                  if(errSC == 1)
                     fprintf(
                        stderr,
                        "%s %s has only one sequence\n",
                        segFlagsStr[siSeg],
                        argAryStr[siArg]
                     );
                  else
                     fprintf(
                        stderr,
                        "%s %s has 3 or more sequences\n",
                        segFlagsStr[siSeg],
                        argAryStr[siArg]
                     );

                  return 2;
               } /*If: there was a error*/
            } /*Else If: segment sequences input*/
         } /*Loop: check if segment input*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun04 Sec02 Sub07 Cat03:
         +   - check if valid input
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(siSeg >= def_NSNum_fluSeg + 1)
         { /*If: input no recongized*/
            fprintf(
               stderr,
               "%s is not recognized\n",
               argAryStr[siArg]
            );

            return 2;
         } /*If: input no recongized*/
      } /*Else: segment input or not recognized*/

      /**************************************************\
      * Fun04 Sec02 Sub08:
      *   - move to next argument
      \**************************************************/

      ++siArg;
   } /*Loop: get user input*/

   return 0;
} /*getInput_diIDs*/

/*-------------------------------------------------------\
| Main:
|   - main driver function for primFind
| Input:
|   - numArgsSI:
|     o number of arguments user input in argAryStr
|   - argAryStr:
|     o array of c-strings with user arguments
| Output:
|   - Prints:
|     o
\-------------------------------------------------------*/
int
main(
   int numArgsSI,
   char *argAryStr[]
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '   o main sec01:
   '     - variable declerations
   '   o main sec02:
   '     - get user input and initialize structures
   '   o main sec03:
   '     - set up primer sequences and kmer tables
   '   o main sec04:
   '     - find and print primer positions
   '   o main sec05:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^   - variable declerations
   ^   o main sec01 sub01:
   ^     - output/input variables (one tempory)
   ^   o main sec01 sub02:
   ^     - variables for filtering or storing results
   ^   o main sec01 sub03:
   ^     - variables unique to kmer search
   ^   o main sec01 sub04:
   ^     - structure and file variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec01 Sub01:
   *   - output/input variables (one tempory)
   \*****************************************************/

   schar errSC = 0;      /*returned error*/

   schar *seqFileStr = 0; /*fasta/q (x) with sequences*/
   schar fqBl = 1;       /*1 is a fastq file; 0 is fasta*/

   schar *forSeqStr = forPrimStr_diIds;
   schar *revSeqStr = revPrimStr_diIds;

   /*settings for printing out ids*/
   schar diRnaBl = def_pDIRna_diIDs;
   schar vRnaBl = def_pVRna_diIDs;
   schar partBl = def_partSeq_diIDs;

   schar keepBl = 0;  /*holds if keep di/v RNA read*/
   schar diFlagSC = 0;/*holds result from detectDI_fluST*/

   schar *outFileStr = 0; /*output file*/

   ulong entryOnUL = 0;
   schar *tmpStr = 0;

   /*****************************************************\
   * Main Sec01 Sub02:
   *   - variables for filtering or storing results
   \*****************************************************/

   /*variables for filtering*/
   uint minLenUI = def_minLen_diIDs;
   uint maxLenUI = def_maxLen_diIDs;

   /*variables for holding search output*/
   /*usign 3 to avoid overflow errors*/
   uint codeAryUI[3]; /*matchs/primer*/

   slong scoreArySL[3];
   schar dirArySC[3];

   ulong seqStartAryUL[3];
   ulong seqEndAryUL[3];

   ulong primStartAryUL[3];
   ulong primEndAryUL[3];

   float minPercLenF = def_minPercLen_diIDs;
      /*minimum % of read length between primers*/

   float maxPercLenF = def_maxPercLen_diIDs;
      /*length between primers / expected segment length*/

   /*variables holding return values*/
   float percLenDiffF = 0; /*primers map len versus read*/

   schar segNumSC = 0; /*holds segment number of read*/

   ulong mapLenUL = 0; /*length between primers*/

   /*****************************************************\
   * Main Sec01 Sub03:
   *   - variables unique to kmer search
   \*****************************************************/

   schar fastBl = def_fastSearch_diIDs;
      /*1 = kmer search, 0 = waterman*/

   float minPercScoreF = def_minPercScore_kmerFind;
      /*min score to keep mapping*/

   uchar lenKmerUC = def_lenKmer_kmerFind; /*kmer length*/

   float frameShiftF = def_percShift_kmerFind;
      /*percentage of bases to move when shifting a win*/

   float minPercKmerF = def_minKmerPerc_kmerFind;
      /*minimum percentage of kmers to do a waterman*/

   float extraNtInWinF = def_extraNtInWin_kmerFind;
     /*number of extra (beyond primer) bases in window*/

   /*****************************************************\
   * Main Sec01 Sub04:
   *   - structure and file variables
   \*****************************************************/

   /*structures for primer search*/
   struct refST_kmerFind refStackAryST[2];
   unsigned int longestSeqUI = 0;

   struct tblST_kmerFind kmerTblStackST;
   struct seqST seqStackST;
   struct alnSet alnSetStackST;
   struct fluST fluStackST;

   FILE *seqFILE = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - get user input and initialize structures
   ^   o main sec02 sub01:
   ^     - initialize structures
   ^   o main sec02 sub02:
   ^     - get input
   ^   o main sec02 sub03:
   ^     - initialize structures needing user input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize structures
   \*****************************************************/

   init_seqST(&seqStackST);
   init_alnSet(&alnSetStackST);

   init_tblST_kmerFind(&kmerTblStackST);
   init_refST_kmerFind(&refStackAryST[0]);
   init_refST_kmerFind(&refStackAryST[1]);
   init_fluST(&fluStackST);

   alnSetStackST.gapSS = -4 * def_scoreAdj_alnDefs;
   alnSetStackST.extendSS = -1 * def_scoreAdj_alnDefs;

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get input
   \*****************************************************/

   errSC = setup_fluST(&fluStackST);

   if(errSC)
   { /*If: had a memory error*/
      fprintf(
         stderr,
         "memory error setting up fluST\n"
      );

      goto memErr_main_sec05_sub02;
   } /*If: had a memory error*/

   errSC =
      getInput_diIDs(
         numArgsSI,
         argAryStr,
         &fluStackST,
         &seqFileStr,
         &fqBl,
         &forSeqStr,
         &revSeqStr,
         &minLenUI,
         &maxLenUI,
         &minPercLenF,
         &maxPercLenF,
         &minPercScoreF,
         &diRnaBl,
         &vRnaBl,
         &partBl,
         &fastBl,
         &lenKmerUC,
         &minPercKmerF
      );

   if(errSC)
   { /*If: had error*/
      --errSC; /*remove help message error*/
      errSC <<= 5; /*avoid conflicts with other errors*/
      goto cleanUp_main_sec05_sub04;
   } /*If: had error*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - initialize structures needing user input
   \*****************************************************/

   setup_tblST_kmerFind(
      &kmerTblStackST,
      lenKmerUC
   );

   setup_refST_kmerFind(
      &refStackAryST[0],
      lenKmerUC
   ); /*foward primer*/

   setup_refST_kmerFind(
      &refStackAryST[1],
      lenKmerUC
   ); /*reverse primer*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - open files
   *   o main sec02 sub03 cat01:
   *     - open output file
   *   o main sec02 sub03 cat02:
   *     - open read file
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat01:
   +   - open output file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(!outFileStr || *outFileStr == '-')
      outFILE = stdout;

   else
   { /*Else: user supplied an output file*/
      outFILE =
         fopen(
            (char *) outFileStr,
            "w"
         );

      if(! outFILE)
      { /*If: could not open output file*/
         fprintf(
            stderr,
            "could not open -out %s\n",
            outFileStr
         );
      } /*If: could not open output file*/

      goto fileErr_main_sec05_sub03;
   } /*Else: user supplied an output file*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat02:
   +   - open read file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(! seqFileStr || *seqFileStr == '-')
      seqFILE = stdin;

   else
   { /*Else: user supplied an read file*/
      seqFILE =
         fopen(
            (char *) seqFileStr,
            "r"
         );

      if(! seqFILE)
      { /*If: could not open readFxput file*/
         if(fqBl)
            tmpStr = (schar *) "-fq";
         else
            tmpStr = (schar *) "-fa";

         fprintf(
            stderr,
            "could not open %s %s\n",
            tmpStr,
            seqFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*If: could not open readFxput file*/
   } /*Else: user supplied an read file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^   - set up primer sequences and kmer tables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   seqStackST.seqStr = forSeqStr;

   seqStackST.lenSeqUL =
     lenStr_ulCp(
        forSeqStr,
        0,
        0
      );
  
   seqStackST.idStr = (schar *) "for";

   longestSeqUI = 
      addSeqToRefST_kmerFind(
         &kmerTblStackST,
         &refStackAryST[0],
         &seqStackST,
         minPercKmerF,
         longestSeqUI,
         &alnSetStackST
      ); /*set up the foward primer sequence*/

   seqStackST.seqStr = 0;
   seqStackST.lenSeqUL = 0;
   seqStackST.idStr = 0;

   if(longestSeqUI == 0)
   { /*If: I had an error*/
      fprintf(
         stderr,
         "Memory error when reading primer sequences\n"
      );

      goto memErr_main_sec05_sub02;
   } /*If: I had an error*/


   seqStackST.seqStr = revSeqStr;

   seqStackST.lenSeqUL =
     lenStr_ulCp(
        revSeqStr,
        0,
        0
      );
  
   seqStackST.idStr = (schar *) "rev";

   longestSeqUI = 
      addSeqToRefST_kmerFind(
         &kmerTblStackST,
         &refStackAryST[1],
         &seqStackST,
         minPercKmerF,
         longestSeqUI,
         &alnSetStackST
      ); /*set up the reverse primer sequence*/

   seqStackST.seqStr = 0;
   seqStackST.lenSeqUL = 0;
   seqStackST.idStr = 0;

   if(longestSeqUI == 0)
   { /*If: I had an error*/
      fprintf(
         stderr,
         "Memory error when reading primer sequences\n"
      );

      goto memErr_main_sec05_sub02;
   } /*If: I had an error*/

   errSC =
      (schar)
      prep_tblST_kmerFind(
         &kmerTblStackST,
         extraNtInWinF,
         frameShiftF,
         longestSeqUI
      ); /*set up the kmer table*/

   if(errSC)
   { /*If: I had an error*/
      fprintf(
         stderr,
         "Memory error when reading primer sequences\n"
      );

      goto memErr_main_sec05_sub02;
   } /*If: I had an error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - find and print primer positions
   ^   o main sec04 sub01:
   ^     - read in first sequence and start loop
   ^   o main sec04 sub02:
   ^     - check if sequence passes filters
   ^   o main sec04 sub03:
   ^     - find primer positions
   ^   o main sec04 sub04:
   ^     - check if has primers and is di/mvRNA
   ^   o main sec04 sub05:
   ^     - check if read is within print settings
   ^   o main sec04 sub06:
   ^     - print read (passed checks)
   ^   o main sec04 sub07:
   ^     - get the next sequence
   ^   o main sec04 sub08:
   ^     - check for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - read in first sequence and start loop
   \*****************************************************/

   pidHeader_fluST(outFILE);

   if(fqBl)
   { /*If: i have a fastq file*/
      errSC =
         getFqSeq_seqST(
            seqFILE,
            &seqStackST
         );
   } /*If: i have a fastq file*/

   else
   { /*Else: i have a fasta file*/
      errSC =
         getFaSeq_seqST(
            seqFILE,
            &seqStackST
         );
   } /*Else: i have a fasta file*/

   while(! errSC)
   { /*Loop: find primers in sequences*/

      /**************************************************\
      * Main Sec04 Sub02:
      *   - check if sequence passes filters
      \**************************************************/

      ++entryOnUL;

      if(seqStackST.lenSeqUL < minLenUI)
         goto nextSeq_main_sec04_sub05;

      if(seqStackST.lenSeqUL > maxLenUI)
         goto nextSeq_main_sec04_sub05;

      /**************************************************\
      * Main Sec04 Sub03:
      *   - find primer positions
      \**************************************************/

      if(fastBl)
      { /*If: i am searching with kmers*/
         errSC =
            fxFindPrims_kmerFind(
               &kmerTblStackST,    /*table structure*/
               refStackAryST,      /*primers looking for*/
               2,                  /*number of primers*/
               &seqStackST,        /*sequence to look at*/
               minPercScoreF,      /*min % for waterman*/
               codeAryUI,      /*# times found prim*/
               dirArySC,       /*direction of primer*/
               scoreArySL,     /*scores for prims*/
               seqStartAryUL,  /*1st align seq base*/
               seqEndAryUL,    /*last align seq base*/
               primStartAryUL, /*1st align prim base*/
               primEndAryUL,   /*last align prim bas*/
               &alnSetStackST      /*settings*/
            ); 
      } /*If: i am searching with kmers*/

      else
      { /*Else: slower waterman search*/
         errSC =
            waterFindPrims_kmerFind(
               refStackAryST,      /*primers looking for*/
               2,                  /*number of primers*/
               &seqStackST,        /*sequence to look at*/
               minPercScoreF,      /*min % for waterman*/
               codeAryUI,      /*# times found prim*/
               dirArySC,       /*direction of primer*/
               scoreArySL,     /*scores for prims*/
               seqStartAryUL,  /*1st align seq base*/
               seqEndAryUL,    /*last align seq base*/
               primStartAryUL, /*1st align prim base*/
               primEndAryUL,   /*last align prim bas*/
               &alnSetStackST      /*settings*/
            ); 
      } /*Else: slower waterman search*/

      /**************************************************\
      * Main Sec04 Sub04:
      *   - check if has primers and is di/mvRNA
      \**************************************************/

      if(! errSC)
      { /*If: I found at least on primer*/
         diFlagSC =
            detectDI_fluST(
               &fluStackST,
               seqStackST.seqStr,  /*sequence*/
               codeAryUI,      /*# times found a primer*/
               dirArySC,       /*direction of primer*/
               seqStartAryUL,  /*first aligned seq base*/
               seqEndAryUL,    /*last aligned seq base*/
               &segNumSC,      /*gets segment number*/
               &mapLenUL
         );

         /***********************************************\
         * Main Sec04 Sub05:
         *   - check if read is within print settings
         \***********************************************/

         /*check to see if printing out*/
         keepBl |=
            ((!!(diFlagSC & def_fullFound_fluST))&vRnaBl);

         keepBl |=
            ((!!(diFlagSC & def_diFound_fluST)) &diRnaBl);

         keepBl &=
              ((!!(diFlagSC & def_partSeg_fluST)) &partBl)
            | (!!(diFlagSC & def_segFound_fluST));

         /*find percent of read is mapped region*/
         percLenDiffF = (float) mapLenUL;
         percLenDiffF /= (float) seqStackST.lenSeqUL;
         keepBl &= (percLenDiffF >= minPercLenF);

         /*see if meets maximum length for segment*/
         percLenDiffF = (float) mapLenUL;
         percLenDiffF /=
            (float) fluStackST.lenSegArySI[segNumSC];

         keepBl &= (percLenDiffF <= maxPercLenF);

         /*check if could id the segment (found result)*/
         if(keepBl)
         { /*If: print segment (passes filters)*/

            /********************************************\
            * Main Sec04 Sub06:
            *   - print read (passed checks)
            \********************************************/

            pid_fluST(
                &fluStackST,
                (schar *) seqStackST.idStr, /*read id*/
                segNumSC,      /*segment number*/
                diFlagSC,
                dirArySC,      /*primer direction*/
                scoreArySL,    /*primer scores*/
                refStackAryST[0].maxForScoreF, /*for*/
                refStackAryST[1].maxForScoreF, /*rev*/
                seqStartAryUL, /*1st seq base*/
                seqEndAryUL,   /*last seq base*/
                primStartAryUL,/*1st primer base*/
                primEndAryUL,  /*last primer base*/
                seqStackST.lenSeqUL,
                mapLenUL,
                outFILE
             ); /*print out read ids*/
         } /*If: print segment (passes filters)*/
      } /*If: I found at least on primer*/

      /**************************************************\
      * Main Sec04 Sub07:
      *   - get the next sequence
      \**************************************************/

      nextSeq_main_sec04_sub05:;

      if(fqBl)
      { /*If: i have a fastq file*/
         errSC =
            (schar)
            getFqSeq_seqST(
               seqFILE,
               &seqStackST
            );
      } /*If: i have a fastq file*/

      else
      { /*Else: i have a fasta file*/
         errSC =
            (schar)
            getFaSeq_seqST(
               seqFILE,
               &seqStackST
            );
      } /*Else: i have a fasta file*/
   } /*Loop: find primers in sequences*/

   /*****************************************************\
   * Main Sec04 Sub08:
   *   - check for errors
   \*****************************************************/

   if(errSC != def_EOF_seqST)
   { /*If: i had an error of some kind*/
      if(fqBl)
         tmpStr = (schar *) "-fq";
      else
         tmpStr = (schar *) "-fa";

      if(errSC & def_memErr_seqST)
      { /*If: i had a memory error*/
         fprintf(
            stderr,
            "Memory error when %s %s\n",
            tmpStr,
            seqFileStr
         );

         goto memErr_main_sec05_sub02;
      } /*If: i had a memory error*/

      fprintf(
         stderr,
         "Error on entry %lu of %s %s\n",
         entryOnUL,
         tmpStr,
         seqFileStr
      );

      goto fileErr_main_sec05_sub03;
   } /*If: i had an error of some kind*/

   else
      errSC = 0; /*no errors*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05:
   ^   - clean up
   ^   o main sec05 sub01:
   ^     - clean up after no errors
   ^   o main sec05 sub02:
   ^     - deal with memory errors
   ^   o main sec05 sub03:
   ^     - deal with file errors
   ^   o main sec05 sub04:
   ^     - clean up and exit
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec05 Sub01:
   *   - clean up after no errors
   \*****************************************************/

   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub02:
   *   - deal with memory errors
   \*****************************************************/

   memErr_main_sec05_sub02:;

   errSC = def_memErr_kmerFind;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub03:
   *   - deal with file errors
   \*****************************************************/

   fileErr_main_sec05_sub03:;

   errSC = def_fileErr_kmerFind;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub04:
   *   - clean up and exit
   \*****************************************************/

   cleanUp_main_sec05_sub04:;

   /*fluST for segment searching*/
   freeStack_fluST(&fluStackST);

   /*alignment settings*/
   freeStack_alnSet(&alnSetStackST);

   /*kmer find structures*/
   freeStack_seqST(&seqStackST);
   freeStack_tblST_kmerFind(&kmerTblStackST);

   freeStack_refST_kmerFind(&refStackAryST[0]);
   freeStack_refST_kmerFind(&refStackAryST[1]);

   if(outFILE && outFILE != stdout)
      fclose(outFILE);

   outFILE = 0;

   if(seqFILE && seqFILE != stdin)
      fclose(seqFILE);

   seqFILE = 0;

   return errSC;
} /*main*/

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
