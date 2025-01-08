/*########################################################
# Name: diFrag
#   - find defective influenza segments using a waterman
#     alignment
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries and defined variables
'   o fun01: pversion_diFrag
'     - prints the version number for diFrag
'   o fun02: phelp_diFrag
'     - prints the help message for diFrag
'   o fun03: input_diFrag
'     - gets user input
'   o main:
'     - driver function for diFrag
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

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
#include "../genLib/ulCp.h"
#include "../genLib/charCp.h"

#include "../genBio/seqST.h"
#include "../genBio/kmerCnt.h"
#include "../genBio/samEntry.h"

#include "../genAln/alnSet.h"
#include "../genAln/dirMatrix.h"
#include "../genAln/water.h"

#include "diScan.h"

/*.h only*/
#include "../genLib/dataTypeShortHand.h"

#include "../fluDI.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../genLib/genMath.h"
!   - .c  #include "../genLib/strAry.h"
!   - .c  #include "../genAln/indexToCoord.h"
!   - .c  #include "../diCoordsSrc/diCoords.h"
!   - .h  #include "../geBio/ntTo2Bit.h"
!   - .h  #include "../geBio/ntTo5Bit.h"
!   - .h  #include "../geAln/alnDefs.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_lenKmer_diFrag 7 /*length of one kmer*/
#define def_minPercScore_diFrag 0.5f /*90% min socre*/
#define def_minKmerPerc_diFrag 0.4f  /*40% min kmers*/

#define def_minDels_diFrag 20 /*min del size for DI*/
#define def_minPadNt_diFrag 12

#define def_pDI_diFrag 1 /*1 = print DI reads*/
#define def_pVRna_diFrag 1 /*1 = print vRNA reads*/

/*-------------------------------------------------------\
| Fun01: pversion_diFrag
|   - prints the version number for diFrag
| Input:
|   - outFILE:
|     o file to print version number to
| Output:
|   - Prints:
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_diFrag(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "diFrag from fluDI version: %i-%02i-%02i\n",
      def_year_fluDI,
      def_month_fluDI,
      def_day_fluDI
   );
} /*pversion_diFrag*/

/*-------------------------------------------------------\
| Fun02: phelp_diFrag
|   - prints the help message for diFrag
| Input:
|   - outFILE:
|     o file to print help message to
| Output:
|   - Prints:
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_diFrag(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints the help message for diFrag
   '   o fun02 sec01:
   '     - usage
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - usage
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
    (FILE *) outFILE,
    "diFrag -fq reads.fq -ref ref.fa"
   );

   fprintf(
      (FILE *) outFILE,
      " -out-sam out.sam > out.tsv\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - finds the number of DI events in reads\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - print the input options
   ^   o fun02 sec02 sub01:
   ^     - print input block header
   ^   o fun02 sec02 sub02:
   ^     - print file io
   ^   o fun02 sec02 sub03:
   ^     - print filter for DI events
   ^   o fun02 sec02 sub04:
   ^     - score and kmer variables
   ^   o fun02 sec02 sub05:
   ^     - filters for printing
   ^   o fun02 sec02 sub06:
   ^     - help and version flag
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec02 Sub01:
   *   - print input block header
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "Input:\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub02:
   *   - print file io
   \*****************************************************/

   /*input fastq*/
   fprintf(
      (FILE *) outFILE,
      "  -fq reads.fq: [Required]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o fastq file with reads to find DI events in\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o use \"-\" for stdin\n"
   );

   /*input references (fasta file)*/
   fprintf(
      (FILE *) outFILE,
      "  -ref ref.fa: [Required]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o fasta file with references to map reads to\n"
   );

   /*output tsv*/
   fprintf(
      (FILE *) outFILE,
      "  -out out.tsv: [Optinal; stdout]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o file to save tsv to \n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o use \"-\" for stdout\n"
   );

   /*output alignments (sam file)*/
   fprintf(
      (FILE *) outFILE,
      "  -out-sam out.sam: [Optinal; No]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o sam file to save aligned reads to\n"
   );

   fprintf(
     (FILE *) outFILE,
     "    o use \"-sam-out - -out file.tsv\" for stdout\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub03:
   *   - print filter for DI events
   \*****************************************************/

   /*print min del size*/
   fprintf(
      (FILE *) outFILE,
      "  -min-del %i: [Optinal; %i]\n",
      def_minDels_diFrag,
      def_minDels_diFrag
   );

   fprintf(
      (FILE *) outFILE,
      "    o minimum number of deletions to have a DI\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      event\n"
   );

   /*min nucleotides at start/end of sequence*/
   fprintf(
      (FILE *) outFILE,
      "  -min-nt-pad %u: [Optinal; %u]\n",
      def_minPadNt_diFrag,
      def_minPadNt_diFrag
   );

   fprintf(
      (FILE *) outFILE,
      "    o minimum number of matches, SNPs, and\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      insertions needed at start/end of read\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      before counting large deltions as DI\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub04:
   *   - score and kmer variables
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  -score-min-perc %0.2f: [Optinal; %0.2f]\n",
      def_minPercScore_diFrag,
      def_minPercScore_diFrag
   );

   fprintf(
      (FILE *) outFILE,
      "    o minimum percent score to keep an alignment\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  -kmer-min-perc %0.2f: [Optinal; %0.2f]\n",
      def_minKmerPerc_diFrag,
      def_minKmerPerc_diFrag
   );

   fprintf(
      (FILE *) outFILE,
      "    o minimum percent shared kmers to do an\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      alignment (this is between the best\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      reference and read)\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  -len-kmer %i: [Optinal; %i]\n",
      def_lenKmer_diFrag,
      def_lenKmer_diFrag
   );

   fprintf(
      (FILE *) outFILE,
      "    o length of one kmer (kmer size)\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub05:
   *   - filters for printing
   \*****************************************************/

   /*diRNA flag*/
   if(def_pDI_diFrag)
      fprintf(
         (FILE *) outFILE,
         "  -di: [Optinal; Yes]\n"
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -di: [Optinal; No]\n"
      );

   fprintf(
      (FILE *) outFILE,
      "    o print out reads with DI events\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o disable with \"-no-di\"\n"
   );

   /*vRNA flag*/
   if(def_pVRna_diFrag)
      fprintf(
         (FILE *) outFILE,
         "  -vrna: [Optinal; Yes]\n"
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -vrna: [Optinal; No]\n"
      );

   fprintf(
      (FILE *) outFILE,
      "    o print out reads with no DI events\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o disable with \"-no-vrna\"\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub06:
   *   - help and version flag
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  -h: print this help message and exit\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  -v: print version number and exit\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - output block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "Output:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - prints: tsv with read ids, mapped segment,\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    and number of DI events to -out\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o vRNA means no large deletions were found\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      (could be a fragment or a complete read)\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o diRNA means one or more large deletions\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      were found\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  - prints: alignments to -out-sam if file input\n"
   );
} /*phelp_diFrag*/


/*-------------------------------------------------------\
| Fun03: input_diFrag
|   - gets user input
| Input:
|   - numArgsSI:
|     o number of arguments the user input
|   - argAryStr:
|     o c-string array with user input
|   - fqFileStrPtr:
|     o pointer to c-sting to set to fastq file name
|   - refFileStrPtr:
|     o pointer to c-sting to set to reference
|       (fasta file) name
|   - outFileStrPtr:
|     o pointer to c-sting to set to output tsv file name
|   - samFileStrPtr:
|     o pointer to c-sting to set to output sam file name
|   - scoreMinPercF:
|     o pointer to float to hold the min percent score to
|       keep a read
|   - kmerMinPercF:
|     o pointer to float to hold the min percent of shared
|       kmers need to keep a read
|   - lenKmerUCPtr:
|     o pointer to unsigned char to hold the length of one
|       kmer
|   - minDIDelUIPtr:
|     o pointer to unsigned int to hold min deletion size
|       to call DI
|   - minPadNtUIPtr:
|     o pointer to unsigned int to hold the minimum
|       number of nucleotides before/after a large DI
|       event
|   - pDIBlPtr:
|     o pointer to signed char, set to 1 if printing DI
|       reads
|   - pVRnaBlPtr:
|     o pointer to signed char, set to 1 if printing vRNA
|       reads
| Output:
|   - Modifies:
|     o all input except numArgsSI and argAryStr to hold
|       user input
|   - Returns:
|     o 0 for no errors
|     o 1 if printed help message or version number
|     o 2 for an error
\-------------------------------------------------------*/
signed char
input_diFrag(
   int numArgsSI,
   char *argAryStr[],
   signed char **fqFileStrPtr,
   signed char **refFileStrPtr,
   signed char **outFileStrPtr,
   signed char **samFileStrPtr,
   float *scoreMinPercF,         /*holds min % score*/
   float *kmerMinPercF,          /*holds min % kmers*/
   unsigned char *lenKmerUCPtr,  /*holds kmer size*/
   unsigned int *minDIDelUIPtr,
   unsigned int *minPadNtUIPtr,  /*ends non-del length*/
   signed char *pDIBlPtr,
   signed char *pVRnaBlPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - gets user input
   '   o fun03 sec01:
   '     - variable declerations
   '   o fun03 sec02:
   '     - check if something was input
   '   o fun03 sec03:
   '     - get user input
   '   o fun03 sec04:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint siArg = 1;
   schar *tmpStr = 0;
   schar errSC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - check if something was input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numArgsSI < 1)
      goto phelp_fun03_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - get user input
   ^   o fun03 sec03 sub01:
   ^     - start loop and check for file io
   ^   o fun03 sec03 sub02:
   ^     - get deletion counts
   ^   o fun03 sec03 sub03:
   ^     - read filter (score and kmer) settings
   ^   o fun03 sec03 sub04:
   ^     - kmer size setting
   ^   o fun03 sec03 sub05:
   ^     - print settings
   ^   o fun03 sec03 sub06:
   ^     - check for help message requests
   ^   o fun03 sec03 sub07:
   ^     - check for invalid input
   ^   o fun03 sec03 sub08:
   ^     - move to the next argument
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - start loop and check for file io
   \*****************************************************/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/

      if(
         ! eql_charCp(
            (schar *) "-fq",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) { /*If: fastq file was input*/
         ++siArg;
         *fqFileStrPtr = (schar *) argAryStr[siArg];
      } /*If: fastq file was input*/

      else if(
         ! eql_charCp(
            (schar *) "-ref",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*If: a reference fasta file was input*/
         ++siArg;
         *refFileStrPtr = (schar *) argAryStr[siArg];
      } /*If: a reference fasta file was input*/

      else if(
         ! eql_charCp(
            (schar *) "-out",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: an output file was input*/
         ++siArg;
         *outFileStrPtr = (schar *) argAryStr[siArg];
      } /*Else If: an ouput file was input*/

      else if(
         ! eql_charCp(
            (schar *) "-out-sam",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: an sam output file was input*/
         ++siArg;
         *samFileStrPtr = (schar *) argAryStr[siArg];
      } /*Else If: an sam output file was input*/

      /**************************************************\
      * Fun03 Sec03 Sub02:
      *   - DI detection settings
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-min-del",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: number minimum deletions input*/
         ++siArg;
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToUI_base10str(
               (schar *) argAryStr[siArg],
               minDIDelUIPtr
            );

         if(*tmpStr != '\0')
         { /*If: had non-numeric input or overflow*/
            fprintf(
               stderr,
               "-min-del %s is non-numeric or to large\n",
               argAryStr[siArg]
            );

            goto err_fun03_sec04;
         } /*If: had non-numeric input or overflow*/
      } /*Else If: number minimum deletions input*/

      else if(
         ! eql_charCp(
            (schar *) "-min-nt-pad",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: is the minimum DI deletion length*/
         ++siArg;
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToUI_base10str(
               (schar *) argAryStr[siArg],
               minPadNtUIPtr
            );

         if(*tmpStr != '\0')
         { /*If: had non-numeric entry*/
            fprintf(
               stderr,
               "-min-nt-pad %s; non-numeric or to big\n",
               argAryStr[siArg]
            );
 
            goto err_fun03_sec04;
         } /*If: had non-numeric entry*/
      } /*Else If: is the minimum DI deletion length*/

      /**************************************************\
      * Fun03 Sec03 Sub03:
      *   - read filter (score and kmer) settings
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-score-min-perc",
            (schar *)  argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: input the min score (as %)*/
         ++siArg;
         *scoreMinPercF = atof(argAryStr[siArg]);
         /*error value is 0, which is valid input*/

         if(
               *scoreMinPercF > 1
            || *scoreMinPercF < 1
         ){ /*If: score is to high*/
            fprintf(
               stderr,
               "-score-min-perc %s must be between 0 and",
               argAryStr[siArg]
            );

            fprintf(
               stderr,
               " 1\n"
            );

            goto err_fun03_sec04;
         } /*If: score is to high*/
      } /*Else If: input the min score (as %)*/

      else if(
         ! eql_charCp(
            (schar *) "-kmer-min-perc",
            (schar *)  argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: input the min score (as %)*/
         ++siArg;
         *kmerMinPercF = atof(argAryStr[siArg]);
         /*error value is 0, which is valid input*/

         if(
               *kmerMinPercF > 1
            || *kmerMinPercF < 1
         ){ /*If: score is to high*/
            fprintf(
               stderr,
               "-kmer-min-perc %s must be between 0 and",
               argAryStr[siArg]
            );

            fprintf(
               stderr,
               " 1\n"
            );

            goto err_fun03_sec04;
         } /*If: score is to high*/
      } /*Else If: input the min score (as %)*/

      /**************************************************\
      * Fun03 Sec03 Sub04:
      *   - kmer size setting
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-len-kmer",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: kmer length was input*/
         ++siArg;
         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToUC_base10str(
               (schar *) argAryStr[siArg],
               lenKmerUCPtr
            );

         if(*tmpStr != '\0')
         { /*If: invalid input*/
            fprintf(
               stderr,
               "-len-kmer %s is non-numeric or > 255\n",
               argAryStr[siArg]
            );

            goto err_fun03_sec04;
         } /*If: invalid input*/
      } /*Else If: kmer length was input*/

      /**************************************************\
      * Fun03 Sec03 Sub05:
      *   - print settings
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-di",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pDIBlPtr = 1;

      else if(
         ! eql_charCp(
            (schar *) "-no-di",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pDIBlPtr = 0;

      else if(
         ! eql_charCp(
            (schar *) "-vrna",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pVRnaBlPtr = 1;

      else if(
         ! eql_charCp(
            (schar *) "-no-vrna",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pVRnaBlPtr = 0;

      /**************************************************\
      * Fun03 Sec03 Sub06:
      *   - check for help message requests
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-h",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto phelp_fun03_sec04;

      else if(
         ! eql_charCp(
            (schar *) "--h",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto phelp_fun03_sec04;

      else if(
         ! eql_charCp(
            (schar *) "help",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto phelp_fun03_sec04;

      else if(
         ! eql_charCp(
            (schar *) "-help",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto phelp_fun03_sec04;

      else if(
         ! eql_charCp(
            (schar *) "--help",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto phelp_fun03_sec04;

      /**************************************************\
      * Fun03 Sec03 Sub07:
      *   - check for version number requests
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-v",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto pversion_fun03_sec04;

      else if(
         ! eql_charCp(
            (schar *) "--v",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto pversion_fun03_sec04;

      else if(
         ! eql_charCp(
            (schar *) "version",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto pversion_fun03_sec04;

      else if(
         ! eql_charCp(
            (schar *) "-version",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto pversion_fun03_sec04;

      else if(
         ! eql_charCp(
            (schar *) "--version",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) goto pversion_fun03_sec04;

      /**************************************************\
      * Fun03 Sec03 Sub08:
      *   - check for invalid input
      \**************************************************/

      else
      { /*Else: input not recongnized*/
         fprintf(
            stderr,
            "%s is not recongized\n",
            argAryStr[siArg]
         );

         goto err_fun03_sec04;
      } /*Else: input not recongnized*/

      /**************************************************\
      * Fun03 Sec03 Sub09:
      *   - move to the next argument
      \**************************************************/

      ++siArg;
   } /*Loop: get user input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto cleanUp_fun03_sec04;

   phelp_fun03_sec04:;
   phelp_diFrag(stdout);
   errSC = 1;
   goto cleanUp_fun03_sec04;

   pversion_fun03_sec04:;
   pversion_diFrag(stdout);
   errSC = 1;
   goto cleanUp_fun03_sec04;

   err_fun03_sec04:;
   errSC = 2;
   goto cleanUp_fun03_sec04;

   cleanUp_fun03_sec04:;
   return(errSC);
} /*input_diFrag*/

/*-------------------------------------------------------\
| Main:
|   - driver function for diFrag
| Input:
|   - numArgsSI:
|     o number of arguments the user input
|   - argAryStr:
|     o array of c-strings with user input
| Output:
|   - Prints:
|     o prints tsv with ids, scores, reference mapped to
|       and other data of detected DI fragment to stdout
\-------------------------------------------------------*/
int 
main(
   int numArgsSI,         /*number input arguments*/
   char *argAryStr[]     /*arguments input*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main:
   '   - driver function for diFrag
   '   o main sec01:
   '     - variable declerations
   '   o main sec02:
   '     - initialize and get/check user input
   '   o main sec03:
   '     - open files
   '   o main sec04:
   '     - check for DI sequences
   '   o main sec05:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar *fqFileStr = 0;
   schar *refFileStr = 0;
   schar *outFileStr = 0;
   schar *samFileStr = 0; /*save alignments to*/

   schar pDIRnaBl = def_pDI_diFrag;
   schar pVRnaBl = def_pVRna_diFrag;

   schar errSC = 0;
   schar *tmpStr = 0;

   sint segSI = 0; /*index of segment read mapped to*/
   sint numDIEventsSI = 0; /*number DI events in seq*/
   uint minDIDelUI = def_minDels_diFrag;
   uint minPadNtUI = def_minPadNt_diFrag;

   /*scoring variables for waterman*/
   float minPercScoreF = def_minPercScore_diFrag;
   float minKmerPercF = def_minKmerPerc_diFrag;
      /*min % of shared kmers needed to keep*/

   uchar lenKmerUC = def_lenKmer_diFrag;/*kmer len*/
   sint numKmersSI = 0;     /*how many kmers shared*/
   sint *kmerHeapTblSI = 0; /*for get_kmerCnt*/;
   sint *cntHeapArySI = 0;  /*for get_kmerCnt*/
   uint maxKmerUI = 0;      /*for allocating arrays*/

   struct samEntry samStackST;
   schar *buffHeapStr = 0;
   ulong lenBuffUL = 0;

   struct seqST seqStackST;
   struct alnSet alnSetStackST;
   struct dirMatrix matrixStackST;

   struct kmerCnt *kmerHeapAryST = 0;
   uint numRefUI = 0;

   uint seqEntryUI = 0; /*entry error happend on*/
   FILE *fqFILE = 0;
   FILE *outFILE = 0;
   FILE *samFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - initialize and get/check user input
   ^   o main sec02 sub01:
   ^     - initialize non-error causing structures
   ^   o main sec02 sub02:
   ^     - get user input
   ^   o main sec02 sub03:
   ^     - allocate memory for samEntry structure
   ^   o main sec02 sub04:
   ^     - read in reference sequences
   ^   o main sec02 sub05:
   ^     - allocate memory for kmer arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize non-error causing structures
   \*****************************************************/

   init_seqST(&seqStackST);
   init_alnSet(&alnSetStackST);
   init_dirMatrix(&matrixStackST);
   init_samEntry(&samStackST);

   alnSetStackST.extendSS = -10; /*is -0.1*/

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get user input
   \*****************************************************/

   errSC =
      input_diFrag(
         numArgsSI,
         argAryStr,
         &fqFileStr,
         &refFileStr,
         &outFileStr,
         &samFileStr,
         &minPercScoreF,
         &minKmerPercF,
         &lenKmerUC,
         &minDIDelUI,
         &minPadNtUI,   /*min pad length to count DI*/
         &pDIRnaBl,
         &pVRnaBl
      ); /*get user input*/

   if(errSC)
   { /*If: had input error*/
      --errSC; /*convert help/version error to no error*/
      goto cleanUp_main_sec05_sub04; /*help/version number*/
   } /*If: had input error*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - allocate memory for samEntry structure
   \*****************************************************/

   errSC = setup_samEntry(&samStackST);

   if(errSC)
      goto memErr_main_sec05_sub02;

   /*****************************************************\
   * Main Sec02 Sub04:
   *   - read in reference sequences
   \*****************************************************/

   kmerHeapAryST =
      faToKmerCnt_kmerCnt(
         refFileStr,
         lenKmerUC,
         &numRefUI,
         &errSC
      );

   if(errSC)
   { /*If: had an error*/
      if(errSC == def_memErr_kmerCnt)
      { /*If: had memory error*/
         fprintf(
            stderr,
            "MEMORY error reading -ref %s\n",
            refFileStr
         );

         goto memErr_main_sec05_sub02;
      } /*If: had memory error*/

      else
      { /*Else:was a file error*/
         fprintf(
            stderr,
            "coult not open/invalid entry in -ref %s\n",
            refFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*Else:was a file error*/
   } /*If: had an error*/

   /*****************************************************\
   * Main Sec02 Sub05:
   *   - allocate memory for kmer arrays
   \*****************************************************/

   /*these variables are set to null, so safe to free
   `  when no memory was allocated
   */
   maxKmerUI = 1;

   for(
     segSI = 0;
     segSI < lenKmerUC;
     ++segSI
   ) maxKmerUI <<= 2;

   kmerHeapTblSI = malloc((maxKmerUI + 1) * sizeof(sint));

   if(! kmerHeapTblSI)
      goto memErr_main_sec05_sub02;

   cntHeapArySI = malloc((maxKmerUI + 1) * sizeof(sint));

   if(! cntHeapArySI)
      goto memErr_main_sec05_sub02;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^   - open files
   ^   o main sec03 sub01:
   ^     - open fastq file
   ^   o main sec03 sub02:
   ^     - open output file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec03 Sub01:
   *   - open fastq file
   \*****************************************************/

   if(! fqFileStr || *fqFileStr == '-')
      fqFILE = stdin;

   else
   { /*Else: user input fastq file*/
      fqFILE = 
         fopen(
            (char *) fqFileStr,
            "r"
         ); /*open the fastq file*/

      if(! fqFILE)
      { /*If: could not open file*/
         fprintf(
            stderr,
            "could not open -fq %s\n",
            fqFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*If: could not open file*/
   } /*Else: user input fastq file*/

   /*****************************************************\
   * Main Sec03 Sub02:
   *   - open output file
   \*****************************************************/

   if(! outFileStr || *outFileStr == '-')
      outFILE = stdout;

   else
   { /*Else: user input an output file*/
      outFILE = 
         fopen(
            (char *) outFileStr,
            "w"
         ); /*open output file*/

      if(! outFILE)
      { /*If: could not open file*/
         fprintf(
            stderr,
            "could not open -out %s\n",
            outFileStr
         );

         goto fileErr_main_sec05_sub03;
      } /*If: could not open file*/
   } /*Else: user input an output file*/


   /*****************************************************\
   * Main Sec03 Sub03:
   *   - open sam output file
   \*****************************************************/

   if(samFileStr)
   { /*If: user input an sam output file*/
      if(*samFileStr == '-')
      { /*If: outputing to stdout*/
         if(outFILE == stdout)
         { /*If: have output conflict*/
            fprintf(
               stderr,
               "you can not set both -out and -out-sam\n"
            );

            fprintf(
               stderr,
               "  to stdout (-out is stdout by default)\n"
            );

            goto fileErr_main_sec05_sub03;
         } /*If: have output conflict*/

         samFILE = stdout;
      } /*If: outputing to stdout*/

      else
      { /*Else: sam file is a file*/
          samFILE = 
             fopen(
                (char *) samFileStr,
                "w"
             ); /*open output file*/

          if(! samFILE)
          { /*If: could not open file*/
             fprintf(
                stderr,
                "could not open -out-sam %s\n",
                samFileStr
             );

             goto fileErr_main_sec05_sub03;
          } /*If: could not open file*/
      } /*Else: sam file is a file*/
   } /*If: user input an sam output file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - check for DI sequences
   ^   o main sec04 sub01:
   ^     - print out header(s)
   ^   o main sec04 sub02:
   ^     - read in first sequence and start loop
   ^   o main sec04 sub03:
   ^     - scan for DI sequences
   ^   o main sec04 sub04:
   ^     - print out results
   ^   o main sec04 sub05:
   ^     - move to next sequence
   ^   o main sec04 sub05:
   ^     - check for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - print out header(s)
   \*****************************************************/

   phead_diScan(
      outFILE,
      lenKmerUC
   ); /*print header*/

   if(samFILE)
   { /*If: printing to sam file*/
      fprintf(
         samFILE,
         "@HD\tVN:1.6\tSO:unsorted\tGO:none\n"
      );

      tmpStr =
         (schar *)
         kmerHeapAryST[segSI].forSeqST->idStr;

      while(*tmpStr++ > 32) ;
      *(tmpStr - 1) = '\0';

      for(
         segSI = 0;
         segSI < (sint) numRefUI;
         ++segSI
      ){ /*Loop: print out sequence headers*/

         tmpStr =
            (schar *)
            kmerHeapAryST[segSI].forSeqST->idStr;

         while(*tmpStr++ > 32) ;

         *(tmpStr - 1) = '\0';

         fprintf(
            samFILE,
            "@SQ\tSN:%s\tLN:%lu\n",
            kmerHeapAryST[segSI].forSeqST->idStr + 1,
            kmerHeapAryST[segSI].forSeqST->lenSeqUL
         );
      } /*Loop: print out sequence headers*/

      fprintf(
        samFILE,
        "@PG\tID:%s\tVN:fluDI_%i-%02i-%02i\t%s",
        "diFrag",
        def_year_fluDI,
        def_month_fluDI,
        def_day_fluDI,
        "diFrag"
     ); /*print out first part of program id tag*/

      if(pDIRnaBl)
        fprintf( samFILE, " -di");
      else
        fprintf( samFILE, " -no-di");

      if(pVRnaBl)
        fprintf( samFILE, " -vrna");
      else
        fprintf( samFILE, " -no-vrna");

     if(fqFileStr)
        fprintf( samFILE, " -fq %s", fqFileStr);
      else
        fprintf( samFILE, " -fq -");

      if(outFileStr)
        fprintf( samFILE, " -out %s", outFileStr);
      else
        fprintf( samFILE, " -out -");

     fprintf(samFILE, " -out-sam %s\n", samFileStr);
   } /*If: printing to sam file*/

   else
   { /*Else: remove extra header*/
      for(
         segSI = 0;
         segSI < (sint) numRefUI;
         ++segSI
      ){ /*Loop: remove header metadata*/

         tmpStr =
            (schar *)
            kmerHeapAryST[segSI].forSeqST->idStr;

         while(*tmpStr++ > 32) ;
         *(tmpStr - 1) = '\0';
      } /*Loop: remove header metadata*/
   } /*Else: remove extra header*/

   /*****************************************************\
   * Main Sec04 Sub02:
   *   - read in first sequence and start loop
   \*****************************************************/

   errSC =
     (schar)
     getFqSeq_seqST(
        fqFILE,
        &seqStackST
     );

   while(! errSC)
   { /*Loop: find DI sequences*/
      ++seqEntryUI;

      /**************************************************\
      * Main Sec04 Sub03:
      *   - scan for DI sequences
      \**************************************************/

      /*TODO: FIX MEMORY CRASH WHEN RELLOCATING Q-SCORE
      `   ARRAY FOR SAMSTACKST
      */
      numDIEventsSI =
         waterScan_diScan(
            &seqStackST,    /*sequence*/
            kmerHeapAryST,  /*references*/
            numRefUI,       /*number of refs to map*/
            kmerHeapTblSI,  /*gets kmers in sequence*/
            cntHeapArySI,   /*gets sequenced kmer counts*/
            minPercScoreF,  /*min align score to keep*/
            minKmerPercF,   /*min % of kmers needed*/
            lenKmerUC,      /*length of one kmer*/
            minDIDelUI,     /*min deletion size to be DI*/
            minPadNtUI,     /*min nt at ends before DI*/
            &samStackST,    /*will have aligned sequence*/
            &segSI,         /*returned segment number*/
            &numKmersSI,    /*holds number kmers shared*/
            &alnSetStackST, /*alignment settings*/
            &matrixStackST  /*matrix for waterman to use*/
         );

      if(numDIEventsSI < 0)
      { /*If: had an error*/
         if(numDIEventsSI == def_noMatch_diScan)
            goto getNextSeq_main_sec04_sub04;
            /*likely not a flu segment*/

         if(numDIEventsSI == def_memErr_diScan)
            goto memErr_main_sec05_sub02;
            /*def_memErr_diScan is always negative*/
      } /*If: had an error*/

      /**************************************************\
      * Main Sec04 Sub04:
      *   - print out results
      \**************************************************/

      samStackST.lenRefIdUC =
         cpDelim_charCp(
            samStackST.refIdStr,
            kmerHeapAryST[segSI].forSeqST->idStr,
            '\0'
         ); /*copy the mapped segment id*/

      if(numDIEventsSI)
      { /*If: DI sequence*/
         if(pDIRnaBl)
            pfrag_diScan(
               &samStackST,
               numDIEventsSI,
               kmerHeapAryST[segSI].forSeqST->lenSeqUL,
               matrixStackST.scoreSL,
               numKmersSI,
               outFILE
            );

         if(samFILE)
            p_samEntry(
               &samStackST,
               &buffHeapStr,
               &lenBuffUL,
               0,
               samFILE
            ); /*print alignment*/
      } /*If: DI sequence*/

      else if(pVRnaBl)
      { /*Else: not DI sequence*/
         pfrag_diScan(
            &samStackST,
            numDIEventsSI,
            kmerHeapAryST[segSI].forSeqST->lenSeqUL,
            matrixStackST.scoreSL,
            numKmersSI,
            outFILE
         );

         if(samFILE)
            p_samEntry(
               &samStackST,
               &buffHeapStr,
               &lenBuffUL,
               0,
               samFILE
            ); /*print alignment*/
      } /*Else: not DI sequence*/

      /**************************************************\
      * Main Sec04 Sub05:
      *   - move to next sequence
      \**************************************************/

      getNextSeq_main_sec04_sub04:;

      errSC =
         (schar)
         getFqSeq_seqST(
            fqFILE,
            &seqStackST
         );
   } /*Loop: find DI sequences*/

   /*****************************************************\
   * Main Sec04 Sub06:
   *   - check for errors
   \*****************************************************/

   if(errSC != def_EOF_seqST)
   { /*If: had an error*/
      if(errSC == def_memErr_seqST)
      { /*If: memory error*/
         fprintf(
            stderr,
            "memory error\n"
         );

         goto memErr_main_sec05_sub02;
      } /*If: memory error*/
      
      fprintf(
         stderr,
         "problem with entry %u in -fq %s\n",
          seqEntryUI,
          fqFileStr
      );

      goto memErr_main_sec05_sub02;
   } /*If: had an error*/

   if(! seqEntryUI)
   { /*If: no sequences in file*/
      fprintf(
         stderr,
         "-fq %s has no reads\n",
          fqFileStr
      );

      goto memErr_main_sec05_sub02;
   } /*If: no sequences in file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05:
   ^   - clean up
   ^   o main sec05 sub01:
   ^     - clean up after no errors
   ^   o main sec05 sub02:
   ^     - memory error; set return value
   ^   o main sec05 sub03:
   ^     - file error; set return value
   ^   o main sec05 sub04:
   ^     - clean up after no error or after errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec05 Sub01:
   *   - clean up after no errors
   \*****************************************************/

   goto cleanUp_main_sec05_sub04; /*no error*/

   /*****************************************************\
   * Main Sec05 Sub02:
   *   - memory error; set return value
   \*****************************************************/

   memErr_main_sec05_sub02:;
   errSC = 1;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub03:
   *   - file error; set return value
   \*****************************************************/

   fileErr_main_sec05_sub03:;
   errSC = 2;
   goto cleanUp_main_sec05_sub04;

   /*****************************************************\
   * Main Sec05 Sub04:
   *   - clean up after no error or after errors
   \*****************************************************/

   cleanUp_main_sec05_sub04:;

   if(
         fqFILE
      && fqFILE != stdin
      && fqFILE != stdout
   ) fclose(fqFILE);

   fqFILE = 0;

   if(
         outFILE
      && outFILE != stdin
      && outFILE != stdout
   ) fclose(outFILE);

   outFILE = 0;

   if(
         samFILE
      && samFILE != stdin
      && samFILE != stdout
   ) fclose(samFILE);

   samFILE = 0;

   free(kmerHeapTblSI);
   free(cntHeapArySI);

   freeStack_seqST(&seqStackST);
   freeStack_alnSet(&alnSetStackST);
   freeStack_dirMatrix(&matrixStackST);

   freeHeapAry_kmerCnt(
      kmerHeapAryST,
      numRefUI
   );

   kmerHeapAryST = 0;

   free(buffHeapStr);
   buffHeapStr = 0;

   return errSC;
} /*main*/
