/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' fluDI SOF: Start Of File
'   - driver function for fluDI
'   o header:
'     - included librares and defined variables for fluDI
'   o fun01: pversion_fluDI
'     - prints version number for fluDI
'   o fun02: phelp_fluDI
'     - prints help message for fluDI
'   o fun03: input_fluDI
'     - get user input
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included librares and defined variables for fluDI
\-------------------------------------------------------*/

#define def_repInterveral_fluDI 5000

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../genLib/base10str.h"
#include "../genLib/numToStr.h"
#include "../genLib/ulCp.h"

#include "../genBio/seqST.h"
#include "../genBio/kmerCnt.h"
#include "../genBio/samEntry.h"
#include "../genBio/tbCon.h"

#include "../genAln/alnSet.h"
#include "../genAln/dirMatrix.h"
#include "../genAln/water.h"
#include "../genAln/kmerFind.h"

#include "../genClust/clustST.h"
#include "../genClust/edClust.h"

#include "../diFragSrc/diScan.h"

#include "../diIdsSrc/fluST.h"

/*no .c files*/
#include "../fluDI.h" /*version numbers*/
#include "../genLib/dataTypeShortHand.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../genLib/charCp.h"
!   - .c  #include "../genLib/shellSort.h"
!   - .c  #include "../genLib/genMath.h"
!   - .c  #include "../genLib/strAry.h"
!   - .c  #include "../genBio/edDist.h"
!   - .c  #include "../genAln/indexToCoord.h"
!   - .c  #include "../genAln/memwater.h"
!   - .h  #include "../genBio/ntTo2Bit.h"
!   - .h  #include "../genBio/ntTo5Bit.h"
!   - .h  #include "../genBio/tbConDefs.h"
!   - .h  #include "../genAln/alnDefs.h"
!   - .h  #include "../diIdsSrc/fluSeg.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_fq_fluDI 1
#define def_fa_fluDI 2

#define def_lenKmer_fluDI 7 /*length of one kmer*/
#define def_minPercScore_fluDI 0.5f /*90% min socre*/
#define def_minKmerPerc_fluDI 0.4f  /*40% min kmers*/

#define def_extend_fluDI -10    /*gap extend; -10 = -0.1*/
#define def_gap_fluDI -1000     /*gap open; -1000 = -10*/

#define def_minDels_fluDI 20   /*min del size for DI*/
#define def_startTrim_fluDI 13 /*first x bases to trim*/
#define def_endTrim_fluDI 13   /*last x bases to trim*/

#define def_minPercLen_fluDI 0.85f
#define def_maxPercLen_fluDI 1.1f /*110%*/

signed char *glob_prefixStr = (schar *) "out";

#define def_diPrimDI_fluDI 1
#define def_diPrimVRna_fluDI 2
#define def_diPrimNoCall_fluDI 4
#define def_diFragDI_fluDI 8
#define def_diFragVRna_fluDI 16
#define def_segMatch_fluDI 32
#define def_segMismatch_fluDI 64

static signed char
   *forPrimStr_fluDI =
      (signed char *) "AGCGAAAGCAGG";
static signed char
   *revPrimStr_fluDI =
      (signed char *) "AGTAGAAACAAGG";

/*this is copied from fluSeg. Here to avoid unused
`  variables error
*/
#define def_PB2Num_fluDI 0
#define def_PB1Num_fluDI 1
#define def_PANum_fluDI 2
#define def_HANum_fluDI 3
#define def_NPNum_fluDI 4
#define def_NANum_fluDI 5
#define def_MNum_fluDI 6
#define def_NSNum_fluDI 7 /*always last segment*/

static signed char
   segIdAryStr_fluDI[def_NSNum_fluDI + 1][4] =
   {
      "PB2",
      "PB1",
      "PA", 
      "HA", 
      "NP", 
      "NA", 
      "M",  
      "NS"  
   }; /*forIdAryStr_fluSeg*/


/*-------------------------------------------------------\
| Fun01: pversion_fluDI
|   - prints version number for fluDI
| Input:
|   - outFILE:
|     o pionter to FILE to print version number to
| Output:
|   - Prints:
|     o fluDI version number to outFILE
\-------------------------------------------------------*/
void
pversion_fluDI(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "fluDI from fluDI version:%i-%02i-%02i\n",
      def_year_fluDI,
      def_month_fluDI,
      def_day_fluDI
   );
} /*pversion_fluDI*/

/*-------------------------------------------------------\
| Fun02: phelp_fluDI
|   - prints help message for fluDI
| Input:
|   - outFILE:
|     o pionter to FILE to print help message to
| Output:
|   - Prints:
|     o fluDI help message to outFILE
\-------------------------------------------------------*/
void
phelp_fluDI(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints help message for fluDI
   '   o fun02 sec01:
   '     - print usage entry 
   '   o fun02 sec02:
   '     - print input
   '   o fun02 sec03:
   '     - print output
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - print usage entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
     (FILE *) outFILE,
     "fluDI -fq reads.fastq -ref ref.fasta -prefix name\n"
   );

   fprintf(
     (FILE *) outFILE,
     "  find DI reads and build consensues for DI reads\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - print input
   ^   o fun02 sec02 sub01:
   ^     - input block header
   ^   o fun02 sec02 sub02:
   ^     - file IO
   ^   o fun02 sec02 sub03:
   ^     - DI filtering
   ^   o fun02 sec02 sub04:
   ^     - score and kmer variables
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
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  File IO:\n"
   );


   fprintf(
      (FILE *) outFILE,
      "    -fq reads.fastq: [Required]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      o fastq with DI reads to scan\n"
   );


   fprintf(
      (FILE *) outFILE,
      "    -fa reads.fasta: [Replaces -fq]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      o fasta with DI reads to scan\n"
   );


   fprintf(
      (FILE *) outFILE,
      "    -ref ref.fasta: [Required]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      o fasta with reference sequences, each\n"
   );

   fprintf(
      (FILE *) outFILE,
      "        sequence header should start with\n"
   );

   fprintf(
      (FILE *) outFILE,
      "        segment/segment number (>segment_name)\n"
   );


   fprintf(
      (FILE *) outFILE,
      "    -prefix name: [Optional; %s]\n",
      glob_prefixStr
   );

   fprintf(
      (FILE *) outFILE,
      "      o prefix for output file names\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub03:
   *   - DI filtering
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  DI Filtering:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    -min-del %u: [Optinal; %u]\n",
      def_minDels_fluDI,
      def_minDels_fluDI
   );

   fprintf(
      (FILE *) outFILE,
      "      o minimum number of deletions to have a DI\n"
   );


   fprintf(
      (FILE *) outFILE,
      "    -start-trim %u: [Optinal; %u]\n",
      def_startTrim_fluDI,
      def_startTrim_fluDI
   );

   fprintf(
      (FILE *) outFILE,
      "      o remove large deltions before -start-trim\n"
   );


   fprintf(
      (FILE *) outFILE,
      "    -end-trim %u: [Optinal; %i]\n",
      def_endTrim_fluDI,
      def_endTrim_fluDI
   );

   fprintf(
      (FILE *) outFILE,
      "      o remove large deltions after -end-trim\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub04:
   *   - score and kmer variables
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  Reference Detection:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    -score-min-perc %0.2f: [Optinal; %0.2f]\n",
      def_minPercScore_fluDI,
      def_minPercScore_fluDI
   );

   fprintf(
      (FILE *) outFILE,
      "      o minimum percent score to keep alignment\n"
   );


   fprintf(
      (FILE *) outFILE,
      "    -kmer-min-perc %0.2f: [Optinal; %0.2f]\n",
      def_minKmerPerc_fluDI,
      def_minKmerPerc_fluDI
   );

   fprintf(
     (FILE *) outFILE,
     "      o minimum percent shared kmers to keep read\n"
   );


   fprintf(
      (FILE *) outFILE,
      "    -len-kmer %i: [Optinal; %i]\n",
      def_lenKmer_fluDI,
      def_lenKmer_fluDI
   );

   fprintf(
      (FILE *) outFILE,
      "      o length of one kmer (kmer size)\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub0x:
   *   - score and kmer variables
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  Misc:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    -h: print this help message and exit\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    -v: print version number and exit\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - print output
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "Output:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - DI reads to prefix-di.sam\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - DI consensuses to prefix-di-cons.sam\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - vRNA reads to prefix-vRNA.sam\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - vRNA consensuses to prefix-vRNA-cons.sam\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - DI and vRNA ids to prefix.tsv\n"
   );
} /*phelp_fluDI*/

/*-------------------------------------------------------\
| Fun03: input_fluDI
|   - get user input
| Input:
|   - argAryStr:
|     o c-string array with user input
|   - numArgsSI:
|     o number arguments user input
|   - fxFileStrPtr:
|     o pointer to c-string to piont to fastx file with
|       reads to scan for DI's
|   - fxFlagBlPtr:
|     o set to def_fq_fluDI if fastq file input
|     o set to def_fa_fluDI if fasta file
|   - refFileStrPtr:
|     o pointer to c-string to point to reference fasta
|   - preifixFileStrPtr:
|     o pointer to c-string to point to prefix
|   - minDelUIPtr:
|     o unsigned int pointer to hold min DI deletion size
|   - startTrimSIPtr:
|     o signed int pointer to hold first base to have DI
|   - endTrimSIPtr:
|     o signed int pointer to hold last base to have DI
|   - lenKmerUCPtr:
|     o unsigned char pointer to hold kmer length
|   - minScorePercFPtr:
|     o unsigned float pointer to hold minimum percent
|       alignment score
|   - minScorePercFPtr:
|     o unsigned float pointer to hold minimum percent
|       number of kmers needed to keep reference
| Output:
|   - Prints:
|     - if requested; version number to stdout
|     - if requested; help messsage to stdout
|     - errors to stderr
|   - Modifies:
|     o all input variables to have user input
|   - Returns:
|     o 0 for no errors
|     o 1 for version number
|     o 1 for help message
|     o 2 for errors
\-------------------------------------------------------*/
signed char
input_fluDI(
   char *argAryStr[],
   int numArgsSI,
   signed char **fxFileStrPtr,
   signed char *fxFlagBlPtr,
   signed char **refFileStrPtr,
   signed char **prefixStrPtr,
   unsigned int *minDelUIPtr, /*minimum deletion size*/
   signed int *startTrimSIPtr,   /*trim DI until*/
   signed int *endTrimSIPtr,     /*trim all DI's after*/
   unsigned char *lenKmerUCPtr,  /*length of one kmer*/
   float *minScorePercFPtr,      /*min % score for aln*/
   float *minKmerPercFPtr        /*min % kmer score*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - get user input
   '   o fun03 sec01:
   '     - variable declarations
   '   o fun03 sec02:
   '     - check number input arguments
   '   o fun03 sec03:
   '     - get user input
   '   o fun03 sec04:
   '     - return results
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   int siArg = 1;
   schar *tmpStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - check number input arguments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numArgsSI < 1)
   { /*If: printing help message*/
      phelp_fluDI(stdout);
      goto phelp_fun03_sec04;
   } /*If: printing help message*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - get user input
   ^   o fun03 sec03 sub01:
   ^     - file IO + start loop
   ^   o fun03 sec03 sub02:
   ^     - DI filtering
   ^   o fun03 sec03 sub0x:
   ^     - help message checks
   ^   o fun03 sec03 sub0y:
   ^     - version number checks
   ^   o fun03 sec03 sub0z:
   ^     - unkown input
   ^   o fun03 sec03 sub0w:
   ^     - move to next argument
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - file IO + start loop
   \*****************************************************/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/

      if(
         ! eqlNull_ulCp(
            (schar *) "-fq",
            (schar *) argAryStr[siArg]
         )
      ){ /*If: fastq file input*/
         ++siArg;
         *fxFileStrPtr = (schar *) argAryStr[siArg];
         *fxFlagBlPtr = def_fq_fluDI;
      } /*If: fastq file input*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-fa",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: fasta file input*/
         ++siArg;
         *fxFileStrPtr = (schar *) argAryStr[siArg];
         *fxFlagBlPtr = def_fa_fluDI;
      } /*Else If: fasta file input*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-ref",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: refernce file input*/
         ++siArg;
         *refFileStrPtr = (schar *) argAryStr[siArg];
      } /*Else If: reference file input*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-prefix",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: prefix input*/
         ++siArg;
         *prefixStrPtr = (schar *) argAryStr[siArg];
      } /*Else If: prefix Input*/

      /**************************************************\
      * Fun03 Sec03 Sub02:
      *   - DI filtering
      \**************************************************/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-min-del",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: DI deltion size*/
         ++siArg;
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToUI_base10str(
               tmpStr,
               minDelUIPtr
            );

         if(*tmpStr != '\0')
         { /*If: error*/
            fprintf(
               stderr,
               "-min-del %s is to large or non-numeric\n",
               argAryStr[siArg]
            );

            goto err_fun03_sec04;
         } /*If: error*/
      } /*Else If: DI deltion size*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-start-trim",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: last base in start to trim at*/
         ++siArg;
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToSI_base10str(
               tmpStr,
               startTrimSIPtr
            );

         if(*tmpStr != '\0')
         { /*If: error*/
            fprintf(
              stderr,
              "-start-trim %s; to large or non-numeric\n",
              argAryStr[siArg]
            );

            goto err_fun03_sec04;
         } /*If: error*/
      } /*Else If: last base in start to trim at*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-end-trim",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: last base to support DI*/
         ++siArg;
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToSI_base10str(
               tmpStr,
               endTrimSIPtr
            );

         if(*tmpStr != '\0')
         { /*If: error*/
            fprintf(
              stderr,
              "-end-trim %s is to large or non-numeric\n",
              argAryStr[siArg]
            );

            goto err_fun03_sec04;
         } /*If: error*/
      } /*Else If: last base to support DI*/

      /**************************************************\
      * Fun03 Sec03 Sub03:
      *   - scoring/kmer variables
      \**************************************************/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-len-kmer",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: kmer length*/
         ++siArg;
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToUC_base10str(
               tmpStr,
               lenKmerUCPtr
            );

         if(*tmpStr != '\0')
         { /*If: error*/
            fprintf(
              stderr,
              "-len-kmer %s is to large or non-numeric\n",
              argAryStr[siArg]
            );

            goto err_fun03_sec04;
         } /*If: error*/
      } /*Else If: kmer length*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-score-min-perc",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: min percent aligment score*/
         ++siArg;
         *minScorePercFPtr = atof(argAryStr[siArg]);
      } /*Else If: min percent aligment score*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-kmer-min-perc",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: min percent kmers*/
         ++siArg;
         *minKmerPercFPtr = atof(argAryStr[siArg]);
      } /*Else If: min percent kmers*/

      /**************************************************\
      * Fun03 Sec03 Sub0x:
      *   - help message checks
      \**************************************************/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-h",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: help message request*/
         phelp_fluDI(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message request*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "--h",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: help message request*/
         phelp_fluDI(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message request*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "help",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: help message request*/
         phelp_fluDI(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message request*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-help",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: help message request*/
         phelp_fluDI(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message request*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "--help",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: help message request*/
         phelp_fluDI(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message request*/

      /**************************************************\
      * Fun03 Sec03 Sub0y:
      *   - version number checks
      \**************************************************/
 
      else if(
         ! eqlNull_ulCp(
            (schar *) "-v",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: version number request*/
         pversion_fluDI(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version request*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "--v",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: version number request*/
         pversion_fluDI(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version request*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "version",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: version number request*/
         pversion_fluDI(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version request*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "-version",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: version number request*/
         pversion_fluDI(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version request*/

      else if(
         ! eqlNull_ulCp(
            (schar *) "--version",
            (schar *) argAryStr[siArg]
         )
      ){ /*Else If: version number request*/
         pversion_fluDI(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version request*/

      /**************************************************\
      * Fun03 Sec03 Sub0z:
      *   - unkown input
      \**************************************************/

      else
      { /*Else: no idea what is*/
         fprintf(
            stderr,
            "%s is not recongnized\n",
            argAryStr[siArg]
         );

         goto err_fun03_sec04;
      } /*Else: no idea what is*/

      /**************************************************\
      * Fun03 Sec03 Sub0w:
      *   - move to next argument
      \**************************************************/

      ++siArg;
   } /*Loop: get user input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - return results
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   phelp_fun03_sec04:;
   pversion_fun03_sec04:;
      return 1;

   err_fun03_sec04:;
      return 2;
} /*input_fluDI*/

/*-------------------------------------------------------\
| Main:
|   - driver function for fluDI
| Input:
|   - argAryStr:
|     o c-string array with user input
|   - numArgsSI:
|     o number of input arguments
| Output:
|   - Prints:
|     o file output:
|       * DI consensus to "prefix-di-con.fa"
|       * non-DI consensus to "prefix-vRNA-con.fa"
|       * DI ids to "prefix-di.ids"
|       * non-DI ids to "prefix-vRNA.ids"
|       * DI reads to "prefix-di-reference-cluster.sam"
|       * non-DI reads to "prefix-vRNA-ref-cluster.sam"
|     o non-analysis output:
|       * errors to stderr
|       * help message (if requested) to stdout
|       * version number (if requested) to stdout
|   - Returns:
|     o 0 for no errors
|     o 1 for an error
\-------------------------------------------------------*/
int
main(
   int numArgsSI,
   char *argAryStr[]
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '   - driver function for fluDI
   '   o main sec01:
   '     - variable declarations
   '   o main sec02:
   '     - initialize, get input, and check input
   '   o main sec03:
   '     - set up primer sequences for kmer scan
   '   o main sec04:
   '     - print headers
   '   o main Sec05:
   '     - find DI reads
   '   o main sec06:
   '     - cluster reads
   '   o main sec07:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^   - variable declarations
   ^   o main sec01 sub01:
   ^     - variable declarations
   ^   o main sec01 sub02:
   ^     - reference sequence detection (alignment scan)
   ^   o main sec01 sub03:
   ^     - DI waterman alignment variables
   ^   o main sec01 sub04:
   ^     - di primer detection variables
   ^   o main sec01 sub05:
   ^     - clustering
   ^   o main sec01 sub06:
   ^     - file IO
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec01 Sub01:
   *   - variable declarations
   \*****************************************************/

   schar errSC = 0;
   schar *tmpStr = 0; /*temporary variable*/

   schar *fxFileStr = 0; /*fastx file with reads*/
   schar fxFlagBl = def_fq_fluDI; /*fastq or fasta input*/
   schar *refFileStr = 0;/*fasta file with references*/
   schar *prefixStr = 0; /*prefix to name output*/

   /*****************************************************\
   * Main Sec01 Sub02:
   *   - reference sequence detection (alignment scan)
   \*****************************************************/

   /*for reference sequences*/
   struct kmerCnt *refHeapAryST = 0; /*sequences + kmers*/
   uint numRefUI = 0;                /*number references*/
   sint *kmerHeapTblSI = 0;          /*kmer table*/
   sint *cntHeapArySI = 0;           /*kmer counts*/
   sint fragSegSI = 0;

   ulong numSeqUL = 0;
   uint tmpUI = 0;

   /*****************************************************\
   * Main Sec01 Sub03:
   *   - DI waterman alignment variables
   \*****************************************************/

   /*for alignment*/
   uchar lenKmerUC = def_lenKmer_fluDI;
   float scorePercF = def_minPercScore_fluDI;
   float kmerPercF = def_minKmerPerc_fluDI;
   sint numKmersSI = 0;

   /*last/first base to start trimming ends at*/
   uint minDelUI = def_minDels_fluDI;
   sint startTrimSI = def_startTrim_fluDI;
   sint endTrimSI = def_endTrim_fluDI;

   struct alnSet alnSetStackST; /*alingment settings*/
   struct dirMatrix matrixStackST; /*get alignment*/
   struct seqST seqStackST;     /*holds reads*/
   sint numDIEventsSI = 0;      /*# DI events in read*/
   
   struct samEntry samStackST; /*for printing/ other*/
   schar *buffHeapStr = 0;     /*sam printing/reading*/
   ulong lenBuffUL = 0;        /*sam printing/reading*/

   schar alnSegSC = 0; /*for comparing segments*/

   /*****************************************************\
   * Main Sec01 Sub04:
   *   - di primer detection variables
   \*****************************************************/

   schar diFlagSC = 0; /*flag for how classified*/
   schar diResultSC = 0;

   /*omni primer sequences*/
   schar *forSeqStr = forPrimStr_fluDI;
   schar *revSeqStr = revPrimStr_fluDI;

   /*settings*/
   float minPercLenF = def_minPercLen_fluDI;
   float maxPercLenF = def_maxPercLen_fluDI;

   uchar primKmerUC = def_lenKmer_kmerFind;
   float minPercScoreF = def_minPercScore_kmerFind;
   float frameShiftF = def_percShift_kmerFind;
   float minPrimPercKmerF = def_minKmerPerc_kmerFind;
   float extraNtInWinF = def_extraNtInWin_kmerFind;

   struct refST_kmerFind refKmerStackST[2];
   uint maxSeqLenUI = 0;

   struct alnSet alnPrimSetStackST;
   struct fluST fluStackST;
   struct tblST_kmerFind kmerTblStackST;

   /*hold primer identification stats*/
   schar diIdSegSC = 0; /*holds segment number of read*/
   ulong mapLenUL = 0; /*length between primers*/
   float percLenDiffF = 0; /*primers map len versus read*/

   /*holds coordinates (3 is for safety)*/
   uint codeAryUI[3]; /*matchs/primer*/

   slong scoreArySL[3];
   schar dirArySC[3];

   ulong seqStartAryUL[3];
   ulong seqEndAryUL[3];

   ulong primStartAryUL[3];
   ulong primEndAryUL[3];

   /*****************************************************\
   * Main Sec01 Sub05:
   *   - clustering
   \*****************************************************/

   /*for clustering*/
   struct con_clustST *conHeapST = 0; /*has consensuses*/
   struct index_clustST *indexHeapST = 0; /*clusters*/
   struct set_clustST clustSetStackST;/*settings*/
   struct set_tbCon conSetStackST;    /*settings*/

   #define sizeHead_main 1 << 17
   schar headStr[sizeHead_main]; /*header for clusters*/
   schar *tmpHeadStr = 0;

   /*****************************************************\
   * Main Sec01 Sub06:
   *   - file IO
   \*****************************************************/

   FILE *inFILE = 0;           /*input file*/

   FILE *diFILE = 0;        /*output file*/
   schar diStr[1024];        /*output file name*/

   FILE *vRnaFILE = 0;      /*output file*/
   schar vRnaStr[1024];    /*output file name*/

   FILE *diPrimScanFILE = 0;    /*primer scan tsv file*/
   schar diPrimScanStr[1024];

   FILE *tsvFILE = 0;       /*tsv output file*/
   schar tsvStr[1024];     /*output tsv file name*/

   FILE *reportFILE = 0;       /*report tsv file*/
   schar reportStr[1024];      /*report tsv file name*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - initialize, get input, and check input
   ^   o main sec02 sub01:
   ^     - initialize variables
   ^   o main sec02 sub02:
   ^     - get user input
   ^   o main sec02 sub03:
   ^     - setup structures
   ^   o main sec02 sub04:
   ^     - get reference sequences
   ^   o main sec02 sub05:
   ^     - allocate memory for kmer arrays
   ^   o main sec02 sub06:
   ^     - open read fastx file
   ^   o main sec02 sub07:
   ^     - open DI output file
   ^   o main sec02 sub08:
   ^     - open vRNA output file
   ^   o main sec02 sub09:
   ^     - open tsv output file
   ^   o main sec02 sub10:
   ^     - open di primer scan output file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize variables
   \*****************************************************/

   init_samEntry(&samStackST);
   init_alnSet(&alnSetStackST);
   init_dirMatrix(&matrixStackST);
   init_seqST(&seqStackST);

   init_set_clustST(&clustSetStackST);
   clustSetStackST.minMapqUC = 0; /*I do not assing mapq*/
   clustSetStackST.repIntervalSL=def_repInterveral_fluDI;

   init_set_tbCon(&conSetStackST);

   init_alnSet(&alnPrimSetStackST);
   init_tblST_kmerFind(&kmerTblStackST);
   init_refST_kmerFind(&refKmerStackST[0]);
   init_refST_kmerFind(&refKmerStackST[1]);
   init_fluST(&fluStackST);

   alnSetStackST.extendSS = def_extend_fluDI;
   alnSetStackST.gapSS = def_gap_fluDI;

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get user input
   \*****************************************************/

   errSC =
      input_fluDI(
         argAryStr,
         numArgsSI,
         &fxFileStr,
         &fxFlagBl,
         &refFileStr,
         &prefixStr,
         &minDelUI,
         &startTrimSI,
         &endTrimSI,
         &lenKmerUC,
         &scorePercF,
         &kmerPercF
      ); /*get user input*/

   if(errSC)
   { /*If: had error*/
      --errSC; /*help/version error goes to no error (0)*/
      goto cleanUp_main_sec07_sub03;
   } /*If: had error*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - setup structures
   \*****************************************************/

   errSC = setup_samEntry(&samStackST);

   if(errSC)
   { /*If: had error*/
      fprintf(
         stderr,
         "error allocating memory for first samStruct\n"
       );

      goto err_main_sec07_sub02;
   } /*If: had error*/


   errSC = setup_fluST(&fluStackST);

   if(errSC)
   { /*If: had error*/
      fprintf(
         stderr,
         "error allocating memory for fluST struct\n"
       );

      goto err_main_sec07_sub02;
   } /*If: had error*/


   errSC =
      setup_tblST_kmerFind(
         &kmerTblStackST,
         primKmerUC
      );

   if(errSC)
   { /*If: had error*/
      fprintf(
         stderr,
         "error allocating memory for kmer table\n"
       );

      goto err_main_sec07_sub02;
   } /*If: had error*/


   errSC =
      setup_refST_kmerFind(
         &refKmerStackST[0],
         lenKmerUC
      );

   if(errSC)
   { /*If: had error*/
      fprintf(
         stderr,
         "error allocating memory for first kmer ref\n"
       );

      goto err_main_sec07_sub02;
   } /*If: had error*/


   errSC =
      setup_refST_kmerFind(
         &refKmerStackST[1],
         lenKmerUC
      );

   if(errSC)
   { /*If: had error*/
      fprintf(
         stderr,
         "error allocating memory for second kmer ref\n"
       );

      goto err_main_sec07_sub02;
   } /*If: had error*/

   /*****************************************************\
   * Main Sec02 Sub04:
   *   - get reference sequences
   \*****************************************************/

   refHeapAryST =
      faToKmerCnt_kmerCnt(
         refFileStr,
         lenKmerUC,
         &numRefUI,
         &errSC
   ); /*get reference sequences*/

   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_kmerCnt)
         fprintf(
            stderr,
            "MEMORY error reading -ref %s\n",
            refFileStr
         );

      else
         fprintf(
            stderr,
            "coult not open/invalid entry in -ref %s\n",
            refFileStr
         );

      goto err_main_sec07_sub02;
   } /*If: had error*/

   /*****************************************************\
   * Main Sec02 Sub05:
   *   - allocate memory for kmer arrays
   \*****************************************************/

   /*these variables are set to null, so safe to free
   `  when no memory was allocated
   */
   tmpUI = 1;

   for(
     fragSegSI = 0;
     fragSegSI < lenKmerUC;
     ++fragSegSI
   ) tmpUI <<= 2;

   kmerHeapTblSI =
      malloc((tmpUI + 1) * sizeof(sint));

   if(! kmerHeapTblSI)
      goto err_main_sec07_sub02;

   cntHeapArySI = malloc((tmpUI + 1) * sizeof(sint));

   if(! cntHeapArySI)
      goto err_main_sec07_sub02;

   /*****************************************************\
   * Main Sec02 Sub06:
   *   - open read fastx file
   \*****************************************************/

   if(
         ! fxFileStr
      || *fxFileStr == '-'
   ) inFILE = stdin;
   
   else
   { /*Else: user input file*/
      inFILE =
         fopen(
            (char *) fxFileStr,
            "r"
         );

      if(! inFILE)
      { /*If: unable to open input*/
         if(fxFlagBl == def_fq_fluDI)
            fprintf(
               stderr,
               "unable to open -fq %s\n",
               fxFileStr
            );

         else
            fprintf(
               stderr,
               "unable to open -fa %s\n",
               fxFileStr
            );

         goto err_main_sec07_sub02;
      } /*If: unable to open input*/
   } /*Else: user input file*/

   /*****************************************************\
   * Main Sec02 Sub07:
   *   - open DI output file
   \*****************************************************/

   tmpStr = diStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   *tmpStr++ = '-';
   *tmpStr++ = 'd';
   *tmpStr++ = 'i';
   *tmpStr++ = 'R';
   *tmpStr++ = 'n';
   *tmpStr++ = 'a';
   *tmpStr++ = '.';
   *tmpStr++ = 's';
   *tmpStr++ = 'a';
   *tmpStr++ = 'm';
   *tmpStr = '\0';

   diFILE =
      fopen(
         (char *) diStr,
         "w"
      );

   if(! diFILE)
   { /*If: could not open DI reads file*/
      fprintf(
         stderr,
         "unable to open %s\n",
         diStr
      );

      goto err_main_sec07_sub02;
   } /*If: could not open DI reads file*/

   /*****************************************************\
   * Main Sec02 Sub08:
   *   - open vRNA output file
   \*****************************************************/

   tmpStr = vRnaStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   *tmpStr++ = '-';
   *tmpStr++ = 'v';
   *tmpStr++ = 'R';
   *tmpStr++ = 'N';
   *tmpStr++ = 'A';
   *tmpStr++ = '.';
   *tmpStr++ = 's';
   *tmpStr++ = 'a';
   *tmpStr++ = 'm';
   *tmpStr = '\0';

   vRnaFILE =
      fopen(
         (char *) vRnaStr,
         "w"
      );

   if(! vRnaFILE)
   { /*If: could not open vRNA reads file*/
      fprintf(
         stderr,
         "unable to open %s\n",
         vRnaStr
      );

      goto err_main_sec07_sub02;
   } /*If: could not open vRNA reads file*/

   /*****************************************************\
   * Main Sec02 Sub09:
   *   - open tsv output file
   \*****************************************************/

   tmpStr = tsvStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   *tmpStr++ = '-';
   *tmpStr++ = 'f';
   *tmpStr++ = 'r';
   *tmpStr++ = 'a';
   *tmpStr++ = 'g';
   *tmpStr++ = '.';
   *tmpStr++ = 't';
   *tmpStr++ = 's';
   *tmpStr++ = 'v';
   *tmpStr = '\0';

   tsvFILE =
      fopen(
         (char *) tsvStr,
         "w"
      );

   if(! tsvFILE)
   { /*If: could not open tsv file*/
      fprintf(
         stderr,
         "unable to open %s\n",
         tsvStr
      );

      goto err_main_sec07_sub02;
   } /*If: could not open tsv file*/

   /*****************************************************\
   * Main Sec02 Sub10:
   *   - open di primer scan output file
   \*****************************************************/

   tmpStr = diPrimScanStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   *tmpStr++ = '-';
   *tmpStr++ = 'd';
   *tmpStr++ = 'i';
   *tmpStr++ = '-';
   *tmpStr++ = 'I';
   *tmpStr++ = 'D';
   *tmpStr++ = 's';
   *tmpStr++ = '.';
   *tmpStr++ = 't';
   *tmpStr++ = 's';
   *tmpStr++ = 'v';
   *tmpStr = '\0';

   diPrimScanFILE =
      fopen(
         (char *) diPrimScanStr,
         "w"
      );

   if(! diPrimScanFILE)
   { /*If: could not open tsv file*/
      fprintf(
         stderr,
         "unable to open %s\n",
         diPrimScanStr
      );

      goto err_main_sec07_sub02;
   } /*If: could not open tsv file*/

   /*****************************************************\
   * Main Sec02 Sub11:
   *   - open report output file
   \*****************************************************/

   tmpStr = reportStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   *tmpStr++ = '-';
   *tmpStr++ = 'r';
   *tmpStr++ = 'e';
   *tmpStr++ = 'p';
   *tmpStr++ = 'o';
   *tmpStr++ = 'r';
   *tmpStr++ = 't';
   *tmpStr++ = '.';
   *tmpStr++ = 't';
   *tmpStr++ = 's';
   *tmpStr++ = 'v';
   *tmpStr = '\0';

   reportFILE =
      fopen(
         (char *) reportStr,
         "w"
      );

   if(! reportFILE)
   { /*If: could not open tsv file*/
      fprintf(
         stderr,
         "unable to open %s\n",
         reportStr
      );

      goto err_main_sec07_sub02;
   } /*If: could not open tsv file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^   - set up primer sequences for kmer scan
   ^   o main sec03 sub01:
   ^     - setup forward primer sequence for kmer scan
   ^   o main sec03 sub02:
   ^     - setup reverse primer sequence for kmer scan
   ^   o main sec03 sub03:
   ^     - setup table for kmer scan
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec03 Sub01:
   *   - setup forward primer sequence for kmer scan
   \*****************************************************/

   seqStackST.seqStr = forSeqStr;

   seqStackST.lenSeqUL =
     lenStr_ulCp(
        forSeqStr,
        0,
        0
      );

   seqStackST.idStr = (schar *) "for";

   maxSeqLenUI = 
      addSeqToRefST_kmerFind(
         &kmerTblStackST,
         &refKmerStackST[0],
         &seqStackST,
         minPrimPercKmerF,
         maxSeqLenUI,
         &alnPrimSetStackST
      ); /*set up the foward primer sequence*/

   seqStackST.seqStr = 0;
   seqStackST.lenSeqUL = 0;
   seqStackST.idStr = 0;

   if(! maxSeqLenUI)
   { /*If: I had an error*/
      fprintf(
         stderr,
         "Memory error assiging forward kmer scan prim\n"
      );

      goto err_main_sec07_sub02;
   } /*If: I had an error*/

   /*****************************************************\
   * Main Sec03 Sub02:
   *   - setup reverse primer sequence for kmer scan
   \*****************************************************/

   seqStackST.seqStr = revSeqStr;

   seqStackST.lenSeqUL =
     lenStr_ulCp(
        revSeqStr,
        0,
        0
      );
  
   seqStackST.idStr = (schar *) "rev";

   maxSeqLenUI = 
      addSeqToRefST_kmerFind(
         &kmerTblStackST,
         &refKmerStackST[1],
         &seqStackST,
         minPrimPercKmerF,
         maxSeqLenUI,
         &alnPrimSetStackST
      ); /*set up the reverse primer sequence*/

   seqStackST.seqStr = 0;
   seqStackST.lenSeqUL = 0;
   seqStackST.lenSeqBuffUL = 0;
   seqStackST.idStr = 0;
   seqStackST.lenIdBuffUL = 0;

   if(! maxSeqLenUI)
   { /*If: I had an error*/
      fprintf(
         stderr,
         "Memory error assiging reverse kmer scan prim\n"
      );

      goto err_main_sec07_sub02;
   } /*If: I had an error*/

   /*****************************************************\
   * Main Sec03 Sub03:
   *   - setup table for kmer scan
   \*****************************************************/

   errSC =
      (schar)
      prep_tblST_kmerFind(
         &kmerTblStackST,
         extraNtInWinF,
         frameShiftF,
         maxSeqLenUI
      ); /*set up the kmer table*/

   if(errSC)
   { /*If: I had an error*/
      fprintf(
         stderr,
         "Memory error preparing table for kmer scan\n"
      );

      goto err_main_sec07_sub02;
   } /*If: I had an error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - print headers
   ^   o main sec04 sub01:
   ^     - print tsv header
   ^   o main sec04 sub02:
   ^     - print sam file header
   ^   o main sec04 sub03:
   ^     - print fluDI program header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - print tsv header
   \*****************************************************/

   phead_diScan(
      tsvFILE,
      lenKmerUC
   );

   pidHeader_fluST(diPrimScanFILE);

   /*header for report*/
   fprintf(
      reportFILE,
      "id\tseg_agree\tcall_agree\tlen\taln_len\tprim_len"
   );

   fprintf(
      reportFILE,
      "\taln_seg\tdi_seg\taln_call\tdi_call\n"
   );

   /*****************************************************\
   * Main Sec04 Sub02:
   *   - print sam file header
   *   o main sec04 sub02 cat01:
   *     - print out main header
   *   o main sec04 sub02 cat02:
   *     - print out sequence entries of header
   *   o main sec04 sub02 cat03:
   *     - print fluDI version (program header)
   *   o main sec04 sub02 cat04:
   *     - print fluDI input file
   *   o main sec04 sub02 cat05:
   *     - print fluDI reads and prefix
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec04 Sub02 Cat01:
   +   - print out main header
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
      diFILE,
      "@HD\tVN:1.6\tSO:unsorted\tGO:none\n"
   );

   fprintf(
      vRnaFILE,
      "@HD\tVN:1.6\tSO:unsorted\tGO:none\n"
   );

   tmpHeadStr = headStr;

   tmpHeadStr +=
      cpLine_ulCp(
         tmpHeadStr,
         (schar *) "@HD\tVN:1.6\tSO:unsorted\tGO:none\n"
      );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec04 Sub02 Cat02:
   +   - print out sequence entries of header
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   tmpStr =
      (schar *)
      refHeapAryST[fragSegSI].forSeqST->idStr;

   while(*tmpStr++ > 32) ;
   *(tmpStr - 1) = '\0';

   for(
      fragSegSI = 0;
      fragSegSI < (sint) numRefUI;
      ++fragSegSI
   ){ /*Loop: print out sequence headers*/

      tmpStr =
         (schar *)
         refHeapAryST[fragSegSI].forSeqST->idStr;

      while(*tmpStr++ > 32) ;
      *(tmpStr - 1) = '\0';

      fprintf(
         diFILE,
         "@SQ\tSN:%s\tLN:%lu\n",
         refHeapAryST[fragSegSI].forSeqST->idStr,
         refHeapAryST[fragSegSI].forSeqST->lenSeqUL
      );

      fprintf(
         vRnaFILE,
         "@SQ\tSN:%s\tLN:%lu\n",
         refHeapAryST[fragSegSI].forSeqST->idStr,
         refHeapAryST[fragSegSI].forSeqST->lenSeqUL
      );
   } /*Loop: print out sequence headers*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec04 Sub02 Cat03:
   +   - print fluDI version (program header)
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
     diFILE,
     "@PG\tID:fluDI\tVN:fluDI_%i-%02i-%02i\tfluDI",
     def_year_fluDI,
     def_month_fluDI,
     def_day_fluDI
   ); /*print out first part of program id tag*/

   fprintf(
     vRnaFILE,
     "@PG\tID:fluDI\tVN:fluDI_%i-%02i-%02i\tfluDI",
     def_year_fluDI,
     def_month_fluDI,
     def_day_fluDI
   ); /*print out first part of program id tag*/

   tmpHeadStr +=
      cpLine_ulCp(
         tmpHeadStr,
         (schar *) "@PG\tID:fluDI\tVN:fluDI_"
      );

   numSeqUL = def_year_fluDI;
   *tmpHeadStr++ = 48 + (numSeqUL / 1000);
   numSeqUL %= 1000;
   *tmpHeadStr++ = 48 + (numSeqUL / 100);
   numSeqUL %= 100;
   *tmpHeadStr++ = 48 + (numSeqUL / 10);
   numSeqUL %= 10;
   *tmpHeadStr++ = 48 + numSeqUL;
   numSeqUL = 0;

   *tmpHeadStr++ = '-';
   *tmpHeadStr++ = 48 + (def_month_fluDI / 10);
   *tmpHeadStr++ = 48 + (def_month_fluDI % 10);

   *tmpHeadStr++ = '-';
   *tmpHeadStr++ = 48 + (def_day_fluDI / 10);
   *tmpHeadStr++ = 48 + (def_day_fluDI % 10);

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec04 Sub02 Cat04:
   +   - print fluDI input file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(fxFlagBl == def_fq_fluDI)
   { /*If: fastq file input*/
      fprintf(diFILE, " -fq");
      fprintf(vRnaFILE, " -fq");

      tmpHeadStr +=
         cpLine_ulCp(
            tmpHeadStr,
            (schar *) "\tfluDI -fq "
         );
   } /*If: fastq file input*/

   else
   { /*Else: fasta file input*/
      fprintf(diFILE, " -fa");
      fprintf(vRnaFILE, " -fa");

      tmpHeadStr +=
         cpLine_ulCp(
            tmpHeadStr,
            (schar *) "\tfluDI -fa "
         );
   } /*Else: fasta file input*/

   if(! fxFileStr)
   { /*If: stdin input*/
         fprintf(diFILE, " -");
         fprintf(vRnaFILE, " -");
         *tmpHeadStr++ = '-';
   } /*If: stdin input*/

   else
   { /*Else: file input*/
         fprintf(diFILE, " %s", fxFileStr);
         fprintf(vRnaFILE, " %s", fxFileStr);

         tmpHeadStr +=
            cpLine_ulCp(
               tmpHeadStr,
               fxFileStr
            );
   } /*Else: file input*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec04 Sub02 Cat05:
   +   - print fluDI reads and prefix
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(diFILE, " -ref %s", refFileStr);
   fprintf(vRnaFILE, " -ref %s", refFileStr);

   tmpHeadStr +=
      cpLine_ulCp(
         tmpHeadStr,
         (schar *) " -ref "
      );

   tmpHeadStr +=
      cpLine_ulCp(
         tmpHeadStr,
         refFileStr
      );


   fprintf(diFILE, " -prefix %s\n", prefixStr);
   fprintf(vRnaFILE, " -prefix %s\n", prefixStr);

   tmpHeadStr +=
      cpLine_ulCp(
         tmpHeadStr,
         (schar *) " -prefix "
      );

   tmpHeadStr +=
      cpLine_ulCp(
         tmpHeadStr,
         prefixStr
      ); /*puts null at end*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05:
   ^   - find DI reads
   ^   o main sec05 sub01:
   ^     - read in first entry in fastq or fasta file
   ^   o main sec05 sub02:
   ^     - scan for DI sequences + start loop
   ^   o main sec05 sub03:
   ^     - do primer scan for DI events
   ^   o main sec05 sub04:
   ^     - check if alignment result matches primer scan
   ^   o main sec05 sub05:
   ^     - print out diFrag results
   ^   o main sec05 sub06:
   ^     - move to next sequence
   ^   o main sec05 sub07:
   ^     - check for errors
   ^   o main sec05 sub08:
   ^     - clean up (minor)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec05 Sub01:
   *   - read in first entry in fastq or fasta file
   \*****************************************************/

   if(fxFlagBl == def_fq_fluDI)
      errSC =
        (schar)
        getFqSeq_seqST(
           inFILE,
           &seqStackST
        );
   else
      errSC =
        (schar)
        getFaSeq_seqST(
           inFILE,
           &seqStackST
        );

   while(! errSC)
   { /*Loop: find DI sequences*/

      if(clustSetStackST.repIntervalSL > 0)
      { /*If: reporting status*/

         if(! (numSeqUL % clustSetStackST.repIntervalSL) )
         { /*If: reporting status (multiple of interval)*/
            fprintf(
               stderr,
               "%lu sequences mapped\n",
               numSeqUL
            ); /*allows some progress reporting*/

            fflush(stderr); 
         } /*If: reporting status (multiple of interval)*/

      } /*If: reporting status*/

      ++numSeqUL;

      /*
      fprintf(
         stderr,
         "%09lu ",
         numSeqUL
      );*/ /*DEBUG for debugging*/

      /**************************************************\
      * Main Sec05 Sub02:
      *   - scan for DI sequences + start loop
      \**************************************************/

      diFlagSC = 0;

      numDIEventsSI =
         waterScan_diScan(
            &seqStackST,    /*sequence*/
            refHeapAryST,   /*references*/
            numRefUI,       /*number of refs to map*/
            kmerHeapTblSI,  /*gets kmers in sequence*/
            cntHeapArySI,   /*gets sequenced kmer counts*/
            scorePercF,     /*min align score to keep*/
            kmerPercF,      /*min % of kmers needed*/
            lenKmerUC,      /*length of one kmer*/
            minDelUI,       /*min deletion size to be DI*/
            endTrimSI,      /*min nt at ends before DI*/
            &samStackST,    /*will have aligned sequence*/
            &fragSegSI,         /*returned segment number*/
            &numKmersSI,    /*holds number kmers shared*/
            &alnSetStackST, /*alignment settings*/
            &matrixStackST  /*matrix for waterman to use*/
         );

      if(numDIEventsSI < 0)
      { /*If: had an error; check if memory error*/
         if(numDIEventsSI == def_memErr_diScan)
            goto err_main_sec07_sub02;

         fragSegSI =
             findSeg_diScan(
                &seqStackST,    /*sequence*/
                refHeapAryST,   /*references*/
                numRefUI,       /*number of refs to map*/
                kmerHeapTblSI,  /*gets kmers in sequence*/
                cntHeapArySI, /*gets sequence kmer count*/
                lenKmerUC,      /*length of one kmer*/
                0,
                &numKmersSI     /*number kmers found*/
             ); /*find the best segment*/
      } /*If: had an error; check if memory error*/

      /*
      else
         rmEndDels_diScan(
            &samStackST,
            (sint) minDelUI,
            startTrimSI,
            endTrimSI
         );*/ /*remove DI events at end of reads*/

      /**************************************************\
      * Main Sec05 Sub03:
      *   - do primer scan for DI events
      *   o main sec05 sub03 cat01:
      *     - do primer scan for DI events
      *   o main sec05 sub03 cat02:
      *     - detect if is a DI event
      *   o main sec05 sub03 cat03:
      *     - detect if had enough coverage to make call
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec05 Sub03 Cat01:
      +   - do primer scan for DI events
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      errSC =
         fxFindPrims_kmerFind(
            &kmerTblStackST,    /*table structure*/
            refKmerStackST,     /*primers looking for*/
            2,                  /*number of primers*/
            &seqStackST,        /*sequence to look at*/
            minPercScoreF,      /*min % for waterman*/
            codeAryUI,          /*# times found prim*/
            dirArySC,           /*direction of primer*/
            scoreArySL,         /*scores for prims*/
            seqStartAryUL,      /*1st align seq base*/
            seqEndAryUL,        /*last align seq base*/
            primStartAryUL,     /*1st align prim base*/
            primEndAryUL,       /*last align prim bas*/
            &alnPrimSetStackST  /*settings*/
         ); 

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec05 Sub03 Cat02:
      +   - detect if is a DI event
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(! errSC)
      { /*If: found one (mabye two) primers*/
         diResultSC =
            detectDI_fluST(
               &fluStackST,
               seqStackST.seqStr,/*sequence*/
               codeAryUI,       /*# times found a primer*/
               dirArySC,        /*direction of primer*/
               seqStartAryUL,   /*first aligned seq base*/
               seqEndAryUL,     /*last aligned seq base*/
               &diIdSegSC,       /*gets segment number*/
               &mapLenUL        /*has mapped length*/
         );

         if(diIdSegSC < 0)
         { /*If: unable to classify read*/
            diFlagSC |= def_diPrimNoCall_fluDI;
            goto findAlnSeg_main_sec05_sub04_cat01;
         } /*If: unable to classify read*/
        
         if(diResultSC & def_fullFound_fluST)
         { /*If: was vRNA*/
            diFlagSC |= def_diPrimVRna_fluDI;

            pid_fluST(
                &fluStackST,
                (schar *) seqStackST.idStr, /*read id*/
                diIdSegSC,      /*segment number*/
                diResultSC,
                dirArySC,      /*primer direction*/
                scoreArySL,    /*primer scores*/
                refKmerStackST[0].maxForScoreF, /*for*/
                refKmerStackST[1].maxForScoreF, /*rev*/
                seqStartAryUL, /*1st seq base*/
                seqEndAryUL,   /*last seq base*/
                primStartAryUL,/*1st primer base*/
                primEndAryUL,  /*last primer base*/
                seqStackST.lenSeqUL,
                mapLenUL,
                diPrimScanFILE
             ); /*print out read ids*/

            goto findAlnSeg_main_sec05_sub04_cat01;
         } /*If: was vRNA*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Main Sec05 Sub03 Cat03:
         +   - detect if had enough coverage to make call
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         /*see if have the minimum read coverage*/
         percLenDiffF = (float) mapLenUL;
         percLenDiffF /= (float) seqStackST.lenSeqUL;

         if(percLenDiffF < minPercLenF)
         { /*If: does not cover enough of read*/
            diFlagSC |= def_diPrimNoCall_fluDI;
            goto findAlnSeg_main_sec05_sub04_cat01;
         } /*If: does not cover enough of read*/

         /*see if meets maximum length for segment*/
         percLenDiffF = (float) mapLenUL;
         percLenDiffF /=
            (float) fluStackST.lenSegArySI[diIdSegSC];

         if(percLenDiffF > maxPercLenF)
         { /*If: read primers are past consensus length*/
            diFlagSC |= def_diPrimNoCall_fluDI;
            goto findAlnSeg_main_sec05_sub04_cat01;
         } /*If: read primers are past consensus length*/

         diFlagSC |= def_diPrimDI_fluDI;

         pid_fluST(
             &fluStackST,
             (schar *) seqStackST.idStr, /*read id*/
             diIdSegSC,      /*segment number*/
             diResultSC,
             dirArySC,      /*primer direction*/
             scoreArySL,    /*primer scores*/
             refKmerStackST[0].maxForScoreF, /*for*/
             refKmerStackST[1].maxForScoreF, /*rev*/
             seqStartAryUL, /*1st seq base*/
             seqEndAryUL,   /*last seq base*/
             primStartAryUL,/*1st primer base*/
             primEndAryUL,  /*last primer base*/
             seqStackST.lenSeqUL,
             mapLenUL,
             diPrimScanFILE
         ); /*print out read ids*/
      } /*If: found one (mabye two) primers*/

      else
         diFlagSC |= def_diPrimNoCall_fluDI;

      /**************************************************\
      * Main Sec05 Sub04:
      *   - check if alignment result matches primer scan
      *   o main sec05 sub04 cat01:
      *     - find aligment flu segment
      *   o main sec05 sub04 cat02:
      *     - check diIds only classified or nothing did
      *   o main sec05 sub04 cat03:
      *     - check: ID and kmer count aggree on segment
      *   o main sec05 sub04 cat04:
      *     - check if diFrag only classified (no diIds)
      *   o main sec05 sub04 cat05:
      *     - print diID result
      *   o main sec05 sub04 cat06:
      *     - check both diFrag and diIds got an answer
      *   o main sec05 sub04 cat07:
      *     - print diFrag and diIds results
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec05 Sub04 Cat01:
      +   - find aligment flu segment
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      findAlnSeg_main_sec05_sub04_cat01:;

      alnSegSC = 0;

      if(
            fragSegSI >= 0
         && refHeapAryST[fragSegSI].forSeqST->idStr[0] > '9'
      ){ /*If: non-numeric segment id*/

         switch(
            refHeapAryST[fragSegSI].forSeqST->idStr[0]
         ){ /*Switch: find flu segment id*/
            case 'P':
            /*Case: polymerase: PB2, PB1, PA*/
               if(
                    refHeapAryST[
                       fragSegSI
                    ].forSeqST->idStr[1]
                 == 'A'
               ) alnSegSC = def_PANum_fluDI;

               else if(
                    refHeapAryST[
                       fragSegSI
                    ].forSeqST->idStr[1]
                 != 'B'
               ) alnSegSC = -1; /*no idea what is*/

               else if(
                    refHeapAryST[
                       fragSegSI
                    ].forSeqST->idStr[2]
                 == '1'
               ) alnSegSC = def_PB1Num_fluDI;

               else if(
                    refHeapAryST[
                       fragSegSI
                    ].forSeqST->idStr[2]
                 == '2'
               ) alnSegSC = def_PB2Num_fluDI;

               else
                  alnSegSC = -1; /*no idea what is*/
            
               break;
            /*Case: polymerase: PB2, PB1, PA*/

            case 'H':
            /*Case: HA*/
               if(
                    refHeapAryST[
                       fragSegSI
                    ].forSeqST->idStr[1]
                 == 'A'
               ) alnSegSC = def_HANum_fluDI;

               else
                  alnSegSC = -1; /*no idea what is*/

               break;
            /*Case: HA*/

            case 'N':
            /*Case: NP, NA, and NS*/
               if(
                    refHeapAryST[
                       fragSegSI
                    ].forSeqST->idStr[1]
                 == 'P'
               ) alnSegSC = def_NPNum_fluDI;

               else if(
                    refHeapAryST[
                       fragSegSI
                    ].forSeqST->idStr[1]
                 == 'A'
               ) alnSegSC = def_NANum_fluDI;

               else if(
                    refHeapAryST[
                       fragSegSI
                    ].forSeqST->idStr[1]
                 == 'S'
               ) alnSegSC = def_NSNum_fluDI;

               else
                  alnSegSC = -1; /*no idea what is*/

               break;
            /*Case: NP, NA, and NS*/

            case 'M':
            /*Case: M*/
               alnSegSC = def_MNum_fluDI;
               break;
            /*Case: M*/

            default:
               alnSegSC = -1; /*no idea what is*/

               break;
         } /*Switch: find flu segment id*/
      } /*If: non-numeric segment id*/

      else if(fragSegSI >= 0)
      { /*Else If: numeric (segment was found)*/
         alnSegSC =
              refHeapAryST[fragSegSI].forSeqST->idStr[0]
            - 48  /*48 is ascii 0*/
            - 1;  /*-1 for index one conversion*/
      } /*Else If: numeric (segment was found)*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec05 Sub04 Cat02:
      +   - check if diIds  only classified or nothing did
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(
            numDIEventsSI < 0
         || alnSegSC < 0
      ){ /*If: alignment could not find DI events*/
          if(diFlagSC & def_diPrimNoCall_fluDI)
             goto getNextSeq_main_sec05_sub04;
             /*read not classified by either system*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Main Sec05 Sub04 Cat03:
         +   - check: ID and kmer count aggree on segment
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(alnSegSC == diIdSegSC)
         { /*If: both systems agreeded on same segment*/
            /*read is already aligned*/
            percLenDiffF =(float) samStackST.alnReadLenUI;
            percLenDiffF /=(float) samStackST.readLenUI;

            if(percLenDiffF < def_minPercLen_fluDI)
               alnSegSC = -1;
         } /*If: both systems agreeded on same segment*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Main Sec05 Sub04 Cat04:
         +   - print diID result
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         tmpStr = seqStackST.idStr;

         while(*tmpStr > 32)
            ++tmpStr;

         *tmpStr = '\0';

         fprintf(
            reportFILE,
            "%s",
            seqStackST.idStr + 1
         );


         if(alnSegSC)
         { /*If: segment found by kmer count*/
            if(alnSegSC == diIdSegSC)
            { /*If: same segments found*/
               fprintf(
                  reportFILE,
                  "\tTRUE\t*\t%u\t%lu\t%u",
                  samStackST.readLenUI,
                  mapLenUL,
                  samStackST.alnReadLenUI
               );

               diFlagSC |= def_segMatch_fluDI;
            } /*If: same segments found*/

            else
               fprintf(
                  reportFILE,
                  "\tFALSE\t*\t*\t%lu\t*",
                  mapLenUL
               );

            fprintf(
               reportFILE,
               "\t%s\t%s\t*",
               segIdAryStr_fluDI[alnSegSC],
               segIdAryStr_fluDI[diIdSegSC]
            ); /*If: failed alignment, but segment found*/
         } /*If: segment found by kmer count*/

         else
         { /*Else: only diIDs found a segment*/
            fprintf(
               reportFILE,
               "\t*\t*\t*\t%lu\t*",
               mapLenUL
            );

            fprintf(
               reportFILE,
               "\t*\t%s\t*",
               segIdAryStr_fluDI[diIdSegSC]
            );
         } /*Else: only diIDs found a segment*/

         if(diFlagSC & def_diPrimDI_fluDI)
            fprintf(
               reportFILE,
               "\tdiRna"
            );
         else if(diFlagSC & def_diPrimVRna_fluDI)
            fprintf(
               reportFILE,
               "\tvRna"
            );
         else
            fprintf(
               reportFILE,
               "\t*"
            );

         fprintf(
            reportFILE,
            "\n"
         );
      } /*If: alignment could not find DI events*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec05 Sub04 Cat03:
      +   - check if diFrag only classified (no diIds)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(diFlagSC & def_diPrimNoCall_fluDI)
      { /*Else If: DI ids could not make a call*/
         tmpStr = seqStackST.idStr;

         while(*tmpStr > 32)
            ++tmpStr;

         *tmpStr = '\0';

         fprintf(
            reportFILE,
            "%s",
            seqStackST.idStr + 1
         );


         fprintf(
            reportFILE,
            "\t*\t*\t%u\t*\t%u\t%s\t*",
            samStackST.readLenUI,
            samStackST.alnReadLenUI,
            segIdAryStr_fluDI[alnSegSC]
         );

         if(numDIEventsSI > 0)
            fprintf(
               reportFILE,
               "\tdiRna\t*"
            );
         else
            fprintf(
               reportFILE,
               "\tvRna\t*"
            );

         fprintf(
            reportFILE,
            "\n"
         );
      } /*Else If: DI ids could not make a call*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec05 Sub04 Cat05:
      +   - check if diFrag and diIds got same answers
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else
      { /*Else: diFrag and diIds found an anwser*/
         tmpStr = seqStackST.idStr;

         while(*tmpStr > 32)
            ++tmpStr;

         *tmpStr = '\0';

         fprintf(
            reportFILE,
            "%s",
            seqStackST.idStr + 1
         );

         if(diIdSegSC == alnSegSC)
         { /*If: diFrag and diIds agree on segement*/
            fprintf(
               reportFILE,
               "\tTRUE"
            );

            diFlagSC |= def_segMatch_fluDI;
         } /*If: diFrag and diIds agree on segement*/

         else
         { /*Else: diFrag and diIds disagree on segment*/
            fprintf(
               reportFILE,
               "\tFALSE"
            );

            diFlagSC |= def_segMismatch_fluDI;
            /*likely due to primers having wrong tags*/
         } /*Else: diFrag and diIds disagree on segment*/


         if(
               (
                     numDIEventsSI > 0
                  && (diFlagSC & def_diPrimDI_fluDI)
               )
            ||
               (
                     numDIEventsSI == 0
                  && (diFlagSC & def_diPrimVRna_fluDI)
               )
         ){ /*If: diFrag and diIds have same call*/
            fprintf(
               reportFILE,
               "\tTRUE"
            );
         } /*If: diFrag and diIds have same call*/


         else
         { /*Else: diFrag and diIds disagree about DI*/
            fprintf(
               reportFILE,
               "\tFALSE"
            );
         } /*Else: diFrag and diIds disagree about DI*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Main Sec05 Sub04 Cat06:
         +   - print diFrag and diIds results
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         fprintf(
            reportFILE,
            "\t%u\t%lu\t%u",
            samStackST.readLenUI,
            mapLenUL,
            samStackST.alnReadLenUI
         );

         fprintf(
            reportFILE,
            "\t%s\t%s",
            segIdAryStr_fluDI[alnSegSC],
            segIdAryStr_fluDI[diIdSegSC]
         );

         if(numDIEventsSI == 0)
            fprintf(
               reportFILE,
               "\tvRna"
            );
         else
            fprintf(
               reportFILE,
               "\tdiRna"
            );

         if(diFlagSC & def_diPrimDI_fluDI)
            fprintf(
               reportFILE,
               "\tdiRna"
            );
         else if(diFlagSC & def_diPrimVRna_fluDI)
            fprintf(
               reportFILE,
               "\tvRna"
            );
         else
            fprintf(
               reportFILE,
               "\t*"
            );

         fprintf(
            reportFILE,
            "\n"
         );
      } /*Else: diFrag and diIds found an anwser*/
      
      /**************************************************\
      * Main Sec05 Sub07:
      *   - print out diFrag results
      *   o main sec05 sub07 cat01:
      *     - check if diID unique result (diRna only)
      *   o main sec05 sub07 cat02:
      *     - check if diFrag found DI event
      *   o main sec05 sub07 cat03:
      *     - DI: frag supports ID segment (or no call)
      *   o main sec05 sub07 cat04:
      *     - vRna: frag support ID segment (or no call)
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec05 Sub07 Cat01:
      +   - check if diID unique result (diRna only)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(
            numDIEventsSI < 0
         &&
            diFlagSC & def_segMatch_fluDI
         &&
            diFlagSC & def_diPrimDI_fluDI
      ){ /*If: only diIDs could make a call*/
         /*In this case both agree on the segment, but
         `   only diIDs could call diRna or vRna
         */

         p_samEntry(
            &samStackST,
            &buffHeapStr,
            &lenBuffUL,
            1,           /*want to add a tag*/
            diFILE
         ); /*print alignment*/

         fprintf(
            diFILE,
            "\tcall:s:diID\n"
         );

         goto getNextSeq_main_sec05_sub04;
      } /*If: only diIDs could make a call*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec05 Sub07 Cat02:
      +   - check if diFrag found DI event
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(numDIEventsSI < 0)
         goto getNextSeq_main_sec05_sub04;
         /*likely not a flu segment*/

      samStackST.lenRefIdUC =
         cpDelim_ulCp(
            samStackST.refIdStr,
            refHeapAryST[fragSegSI].forSeqST->idStr,
            0,
            '\0'
         ); /*copy the mapped segment id*/

      if(numDIEventsSI >= 0)
      { /*If: diFrag could make a call*/
            pfrag_diScan(
               &samStackST,
               numDIEventsSI,
               refHeapAryST[fragSegSI].forSeqST->lenSeqUL,
               matrixStackST.scoreSL,
               numKmersSI,
               tsvFILE
            );
      } /*If: diFrag could make a call*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec05 Sub07 Cat03:
      +   - frag supports ID segment (or no call) + is DI
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(
             (diFlagSC & def_segMatch_fluDI)
          || (def_diPrimNoCall_fluDI)
      ){ /*If: segments matched or diIDs could not call*/

         if(numDIEventsSI)
         { /*If: diFrag found a DI sequence*/
            p_samEntry(
               &samStackST,
               &buffHeapStr,
               &lenBuffUL,
               1,           /*want to add flag*/
               diFILE
            ); /*print alignment*/

            if(diFlagSC & def_diPrimDI_fluDI)
               fprintf(
                  diFILE,
                  "\tcall:s:both\n"
               );
            else
               fprintf(
                  diFILE,
                  "\tcall:s:frag\n"
               );
         } /*If: if diFrag found a DI sequence*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Main Sec05 Sub07 Cat04:
         +   - vRna: frag support ID segment (or no call)
         \+++++++++++++++++++++++++++++++++++++++++++++*/

         else
         { /*Else: diFrag does not support DI sequence*/
            if(diFlagSC & def_diPrimDI_fluDI)
            { /*If: diIds classifed as DI (assume DI)*/
               p_samEntry(
                  &samStackST,
                  &buffHeapStr,
                  &lenBuffUL,
                  1,
                  diFILE
               ); /*print alignment*/

               fprintf(
                  diFILE,
                  "\tcall:s:diID\n"
               );
            } /*If: diIds classifed as DI (assume DI)*/
            
            else
            { /*Else: assuming vRNA*/
               p_samEntry(
                  &samStackST,
                  &buffHeapStr,
                  &lenBuffUL,
                  1,
                  vRnaFILE
               ); /*print alignment*/

            if(diFlagSC & def_diPrimVRna_fluDI)
               fprintf(
                  vRnaFILE,
                  "\tcall:s:both\n"
               );
            else
               fprintf(
                  vRnaFILE,
                  "\tcall:s:frag\n"
               );
            } /*Else: assuming vRNA*/
         } /*Else: diFrag does not support DI sequence*/
      } /*If: segments matched or diIDs could not call*/

      /**************************************************\
      * Main Sec05 Sub06:
      *   - move to next sequence
      \**************************************************/

      getNextSeq_main_sec05_sub04:;

      if(fxFlagBl == def_fq_fluDI)
         errSC =
           (schar)
           getFqSeq_seqST(
              inFILE,
              &seqStackST
           );
      else
         errSC =
           (schar)
           getFaSeq_seqST(
              inFILE,
              &seqStackST
           );
   } /*Loop: find DI sequences*/

   /*****************************************************\
   * Main Sec05 Sub07:
   *   - check for errors
   \*****************************************************/

   if(errSC != def_EOF_seqST)
   { /*If: had an error*/
      if(errSC == def_memErr_seqST)
         fprintf(
            stderr,
            "memory error\n"
         );
      
      else
      { /*Else: file error*/
         if(fxFlagBl == def_fq_fluDI)
            tmpStr = (schar *) "fq";
         else
            tmpStr = (schar *) "fa";

         fprintf(
            stderr,
            "problem with entry %lu in -%s %s\n",
             numSeqUL,
             tmpStr,
             fxFileStr
         );
      } /*Else: file error*/

      goto err_main_sec07_sub02;
   } /*If: had an error*/

   if(! numSeqUL)
   { /*If: no sequences in file*/
      if(fxFlagBl == def_fq_fluDI)
         tmpStr = (schar *) "fq";
      else
         tmpStr = (schar *) "fa";
      fprintf(
         stderr,
         "-%s %s has no reads\n",
         tmpStr,
         fxFileStr
      );

      goto err_main_sec07_sub02;
   } /*If: no sequences in file*/

   /*****************************************************\
   * Main Sec05 Sub08:
   *   - clean up (minor)
   \*****************************************************/

   if(inFILE != stdin)
      fclose(inFILE);
   inFILE = 0;

   if(tsvFILE != stdout)
      fclose(tsvFILE);
   tsvFILE = 0;

   if(diPrimScanFILE != stdout)
      fclose(diPrimScanFILE);
   diPrimScanFILE = 0;

   if(reportFILE != stdout)
      fclose(reportFILE);
   reportFILE = 0;

   freeStack_alnSet(&alnSetStackST);
   freeStack_dirMatrix(&matrixStackST);
   freeStack_seqST(&seqStackST);

   freeStack_fluST(&fluStackST);
   freeStack_alnSet(&alnPrimSetStackST);
   freeStack_tblST_kmerFind(&kmerTblStackST);
   freeStack_refST_kmerFind(&refKmerStackST[0]);
   freeStack_refST_kmerFind(&refKmerStackST[1]);

   free(cntHeapArySI);
   cntHeapArySI = 0;

   free(kmerHeapTblSI);
   kmerHeapTblSI = 0;

   if(refHeapAryST)
      freeHeapAry_kmerCnt(
         refHeapAryST,
         numRefUI
      );
   refHeapAryST = 0;
   numRefUI = 0;

   /*will reopen as read later*/
   fclose(diFILE);
   diFILE = 0;

   fclose(vRnaFILE);
   vRnaFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec06:
   ^   - cluster reads
   ^   o main sec06 sub01:
   ^     - cluster DI reads
   ^   o main sec05 sub02:
   ^     - cluster vRNA reads
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec06 Sub01:
   *   - cluster DI reads
   *   o main sec06 sub01 cat01:
   *     - cluster DI reads
   *   o main sec06 sub01 cat02:
   *     - print DI read consensuses
   *   o main sec06 sub01 cat03:
   *     - print DI read clusters
   *   o main sec06 sub01 cat04:
   *     - clean up (di)
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec06 Sub01 Cat01:
   +   - cluster DI reads
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   diFILE =
      fopen(
         (char *) diStr,
         "r"
      );

   tmpStr = diStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         (schar *) "-DI-log.txt",
         0,
         '\0'
      );

   inFILE =
      fopen(
         (char *) diStr,
         "w"
      ); /*if errors out, clusters prints to stderr*/
   
   conHeapST =
      cluster_edClust(
         &indexHeapST,     /*holds clusters*/
         &clustSetStackST, /*cluster settings*/
         &conSetStackST,   /*consensus building settings*/
         &samStackST,      /*for sam file reading*/
         &buffHeapStr,     /*for sam file reading*/
         &lenBuffUL,       /*size of buffHeapStr*/
         diFILE,           /*clustering DI hits*/
         inFILE,           /*log for clustering*/
         &errSC
      );
   
   /*check if failed to open file*/
   if(inFILE)
      fclose(inFILE);
   inFILE = 0;

   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_edClust)
      { /*If: memory error*/
         fprintf(
            stderr,
            "memory error during DI cluster step\n"
         );

         goto err_main_sec07_sub02;
      } /*If: memory error*/

      else if(errSC == def_fileErr_edClust)
      { /*If: file error (should never happen)*/
         fprintf(
            stderr,
            "file error during DI cluster step\n"
         );

         goto err_main_sec07_sub02;
      } /*If: file error (should never happen)*/

      fprintf(
         stderr,
         "unable to cluster di reads\n"
      );

      goto clustVRna_main_sec05_sub02_cat01;
   } /*If: had error*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec06 Sub01 Cat02:
   +   - print DI read consensuses
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   tmpStr = diStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         (schar *) "-DI-cons.sam",
         0,
         '\0'
      );

   inFILE =
      fopen(
         (char *) diStr,
         "w"
      ); /*if errors out, clusters prints to stderr*/

   if(! inFILE)
   { /*If: could not open file*/
      fprintf(
         stderr,
         "unable to open output file for DI consensuses\n"
      );

      goto err_main_sec07_sub02;
   } /*If: could not open file*/

   errSC =
      plist_con_clustST(
         conHeapST,
         headStr,
         0,        /*program header is in header*/
         &buffHeapStr,
         &lenBuffUL,
         inFILE    /*file to save consensuses to*/
      ); /*print the consensus to the output file*/
      
   if(errSC)
   { /*If: had error*/
      fprintf(
         stderr,
         "memory error printing DI consensuses\n"
      );

      goto err_main_sec07_sub02;
   } /*If: had memory error*/

   fclose(inFILE);
   inFILE = 0;

   freeHeapList_con_clustST(conHeapST);
   conHeapST = 0;

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec06 Sub01 Cat03:
   +   - print DI read clusters
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   tmpStr = diStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         (schar *) "-DI",
         0,
         '\0'
      );

   errSC =
      pbins_clustST(
         diStr,
         clustSetStackST.clustSI,
         indexHeapST,
         0,          /*no program header*/
         &samStackST,
         &buffHeapStr,
         &lenBuffUL,
         diFILE
      ); /*print DI clusters*/

   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_clustST)
      { /*If: had memory error*/
         fprintf(
            stderr,
            "memory error printing DI clusters\n"
         );

         goto err_main_sec07_sub02;
      } /*If: had memory error*/

      else
      { /*Else: file error of some kind (to few reads)*/
         fprintf(
            stderr,
            "file error, maybe no DI clusters found\n"
         );
      } /*Else: file error of some kind (to few reads)*/
   } /*If: had error*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec06 Sub01 Cat04:
   +   - clean up (di)
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   freeHeap_index_clustST(indexHeapST);
   indexHeapST = 0;

   /*****************************************************\
   * Main Sec06 Sub02:
   *   - cluster vRNA reads
   *   o main sec06 sub02 cat01:
   *     - cluster vRNA reads
   *   o main sec06 sub02 cat02:
   *     - print vRNA read consensuses
   *   o main sec06 sub02 cat03:
   *     - print vRNA read clusters
   *   o main sec06 sub02 cat04:
   *     - clean up (vRNA)
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec06 Sub02 Cat01:
   +   - cluster vRNA reads
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   clustVRna_main_sec05_sub02_cat01:;

   blank_set_clustST(&clustSetStackST);
   /*remove cluster number*/

   vRnaFILE =
      fopen(
         (char *) vRnaStr,
         "r"
      );

   tmpStr = vRnaStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         (schar *) "-vRNA-log.txt",
         0,
         '\0'
      );

   inFILE =
      fopen(
         (char *) vRnaStr,
         "w"
      ); /*if errors out, clusters prints to stderr*/
   
   conHeapST =
      cluster_edClust(
         &indexHeapST,     /*holds clusters*/
         &clustSetStackST, /*cluster settings*/
         &conSetStackST,   /*consensus building settings*/
         &samStackST,      /*for sam file reading*/
         &buffHeapStr,     /*for sam file reading*/
         &lenBuffUL,       /*size of buffHeapStr*/
         vRnaFILE,         /*clustering DI hits*/
         inFILE,           /*log for vRNA clustering*/
         &errSC
     );
   
   /*check if failed to open file*/
   if(inFILE)
      fclose(inFILE);
   inFILE = 0;

   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_edClust)
      { /*If: memory error*/
         fprintf(
            stderr,
            "memory error during vRNA cluster step\n"
         );

         goto err_main_sec07_sub02;
      } /*If: memory error*/

      else if(errSC == def_fileErr_edClust)
      { /*If: file error (should never happen)*/
         fprintf(
            stderr,
            "file error during vRNA cluster step\n"
         );

         goto err_main_sec07_sub02;
      } /*If: file error (should never happen)*/

      fprintf(
         stderr,
         "unable to cluster vRNA reads\n"
      );

      goto cleanUp_main_sec07_sub03;
   } /*If: had error*/


   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec06 Sub02 Cat02:
   +   - print vRNA read consensuses
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   tmpStr = vRnaStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         (schar *) "-vRNA-cons.sam",
         0,
         '\0'
      );

   inFILE =
      fopen(
         (char *) vRnaStr,
         "w"
      ); /*if errors out, clusters prints to stderr*/

   if(! inFILE)
   { /*If: could not open file*/
      fprintf(
         stderr,
         "unable to open output file; vRNA consensuses\n"
      );

      goto err_main_sec07_sub02;
   } /*If: could not open file*/

   errSC =
      plist_con_clustST(
         conHeapST,
         headStr,
         0,        /*program header is in header*/
         &buffHeapStr,
         &lenBuffUL,
         inFILE    /*file to save consensuses to*/
      ); /*print the consensus to the output file*/
      
   if(errSC)
   { /*If: had error*/
      fprintf(
         stderr,
         "memory error printing vRNA consensuses\n"
      );

      goto err_main_sec07_sub02;
   } /*If: had memory error*/

   fclose(inFILE);
   inFILE = 0;

   freeHeapList_con_clustST(conHeapST);
   conHeapST = 0;

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec06 Sub02 Cat03:
   +   - print vRNA read clusters
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   tmpStr = vRnaStr;

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         prefixStr,
         0,
         '\0'
      );

   tmpStr +=
      cpDelim_ulCp(
         tmpStr,
         (schar *) "-vRNA",
         0,
         '\0'
      );

   errSC =
      pbins_clustST(
         vRnaStr,
         clustSetStackST.clustSI,
         indexHeapST,
         0,          /*no program header*/
         &samStackST,
         &buffHeapStr,
         &lenBuffUL,
         vRnaFILE
      ); /*print DI clusters*/

   fclose(vRnaFILE);
   vRnaFILE = 0;

   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_clustST)
      { /*If: had memory error*/
         fprintf(
            stderr,
            "memory error printing vRNA clusters\n"
         );

         goto err_main_sec07_sub02;
      } /*If: had memory error*/

      else
      { /*Else: file error of some kind (to few reads)*/
         fprintf(
            stderr,
            "file error, maybe no vRNA clusters found\n"
         );
      } /*Else: file error of some kind (to few reads)*/
   } /*If: had error*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec06 Sub02 Cat04:
   +   - clean up (vRna)
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   freeHeap_index_clustST(indexHeapST);
   indexHeapST = 0;

   freeStack_set_clustST(&clustSetStackST);
   freeStack_set_tbCon(&conSetStackST);
   freeStack_samEntry(&samStackST);

   free(buffHeapStr);
   buffHeapStr = 0;
   lenBuffUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec07:
   ^   - clean up and return
   ^   o main sec07 sub01:
   ^     - no error clean up
   ^   o main sec07 sub02:
   ^     - error clean up
   ^   o main sec07 sub03:
   ^     - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec07 Sub01:
   *   - no error clean up
   \*****************************************************/

   errSC = 0;
   goto cleanUp_main_sec07_sub03;

   /*****************************************************\
   * Main Sec07 Sub02:
   *   - error clean up
   \*****************************************************/

   err_main_sec07_sub02:;
      errSC = 1;
      goto cleanUp_main_sec07_sub03;

   /*****************************************************\
   * Main Sec07 Sub03:
   *   - clean up and return
   \*****************************************************/

   cleanUp_main_sec07_sub03:;
      freeStack_alnSet(&alnSetStackST);
      freeStack_dirMatrix(&matrixStackST);
      freeStack_seqST(&seqStackST);
      freeStack_set_clustST(&clustSetStackST);
      freeStack_set_tbCon(&conSetStackST);

      freeStack_fluST(&fluStackST);
      freeStack_alnSet(&alnPrimSetStackST);
      freeStack_tblST_kmerFind(&kmerTblStackST);
      freeStack_refST_kmerFind(&refKmerStackST[0]);
      freeStack_refST_kmerFind(&refKmerStackST[1]);

      if(refHeapAryST)
         freeHeapAry_kmerCnt(
            refHeapAryST,
            numRefUI
         );
      refHeapAryST = 0;
      numRefUI = 0;

      if(kmerHeapTblSI)
         free(kmerHeapTblSI);
      kmerHeapTblSI = 0;

      if(cntHeapArySI)
         free(cntHeapArySI);
      cntHeapArySI = 0;

      freeStack_samEntry(&samStackST);

      if(buffHeapStr)
         free(buffHeapStr);
      buffHeapStr = 0;
      lenBuffUL = 0;

      if(conHeapST)
         freeHeapList_con_clustST(conHeapST);
      conHeapST = 0;

      if(indexHeapST)
         freeHeap_index_clustST(indexHeapST);
      indexHeapST = 0;

      if(
            inFILE
         && inFILE != stdin
         && inFILE != stdout
         && inFILE != stderr
      ) fclose(inFILE);
      inFILE = 0;

      if(
            diFILE
         && diFILE != stdin
         && diFILE != stdout
         && diFILE != stderr
      ) fclose(diFILE);
      diFILE = 0;

      if(
            vRnaFILE
         && vRnaFILE != stdin
         && vRnaFILE != stdout
         && vRnaFILE != stderr
      ) fclose(vRnaFILE);
      vRnaFILE = 0;

      if(
            tsvFILE
         && tsvFILE != stdin
         && tsvFILE != stdout
         && tsvFILE != stderr
      ) fclose(tsvFILE);
      tsvFILE = 0;

      if(
            diPrimScanFILE
         && diPrimScanFILE != stdin
         && diPrimScanFILE != stdout
         && diPrimScanFILE != stderr
      ) fclose(diPrimScanFILE);
      diPrimScanFILE = 0;

      if(
            reportFILE
         && reportFILE != stdin
         && reportFILE != stdout
         && reportFILE != stderr
      ) fclose(reportFILE);
      reportFILE = 0;

      return errSC;
} /*main*/
