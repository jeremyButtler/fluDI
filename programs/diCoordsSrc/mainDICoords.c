/*#######################################################\
# Name: mainDiCoords
#   - driver function to get di coordinates
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries and definitions
'   o fun01: pversion_mainDiCoords
'     - prints version number
'   o fun02: phelp_mainDiCoords
'     - prints help message
'   o fun03: input_mainDiCoords
'     - gets user input from an array of c-strings
'   o main:
'     - main driver function to get DI coordinates
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries and definitions
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
#include "../genBio/samEntry.h"

#include "diCoords.h"

/*.h files only*/
#include "../genLib/dataTypeShortHand.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o std #include <stdio.h>
!   o .c  #include "../genLib/numToStr.h"
!   o .c  #include "../genLib/ulCp.h"
!   o .c  #include "../genLib/strAry.h"
!   o .h  #include "../genBio/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_year_mainDiCoords 2024
#define def_month_mainDiCoords 8
#define def_day_mainDiCoords 29

#define def_minDIDelLen_mainDiCoords 20
   /*minimum deletion size to count as a DI event*/

#define def_minPadNt_mainDiCoords 12
   /*minimum bases at end/start needed to count a DI
   `   this is here to avoid large DI and then 4 matches
   */

#define def_pNonDI_mainDiCoords 0
   /*print non-DI entries*/

/*-------------------------------------------------------\
| Fun01: pversion_mainDiCoords
|   - prints version number
| Input:
|   - outFILE:
|     o file to print version number to
| Output:
|   - Prints:
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_mainDiCoords(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "mainDiCoords version: %i-%02i-%02i\n",
      def_year_mainDiCoords,
      def_month_mainDiCoords,
      def_day_mainDiCoords
   );
} /*pversion_mainDiCoords*/

/*-------------------------------------------------------\
| Fun02: phelp_mainDiCoords
|   - prints help message
| Input:
|   - outFILE:
|     o file to print help message to
| Output:
|   - Prints:
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_mainDiCoords(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints help message
   '   o fun02 sec01:
   '     - usage block
   '   o fun02 sec02:
   '     - input block
   '   o fun02 sec03:
   '     - output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - usage block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "mainDiCoords -sam reads.sam > out.tsv"
   );

   fprintf(
      (FILE *) outFILE,
      "  - find coordinates of flu DIs in mapped reads\n"
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - input block
   ^   o fun02 sec02 sub01:
   ^     - input header
   ^   o fun02 sec02 sub02:
   ^     - file io input
   ^   o fun02 sec02 sub03:
   ^     - DI filter settings
   ^   o fun02 sec02 sub04:
   ^     - print settings
   ^   o fun02 sec02 sub05:
   ^     - help and version number
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec02 Sub01:
   *   - input header
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "Input:\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub02:
   *   - file io input
   \*****************************************************/

   /*input*/
   fprintf(
      (FILE *) outFILE,
      "  -sam reads.sam: [Required]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o sam file with mapped reads to search\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o use \"-sam -\" for stdin\n"
   );

   /*output*/
   fprintf(
      (FILE *) outFILE,
      "  -out out.tsv: [Optinal: stdout]\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o file to output to\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o use \"-out -\" for stdout\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub03:
   *   - DI filter settings
   \*****************************************************/

   /*print min del size*/
   fprintf(
      (FILE *) outFILE,
      "  -min-del %u: [Optinal; %u]\n",
      def_minDIDelLen_mainDiCoords,
      def_minDIDelLen_mainDiCoords
   );

   fprintf(
      (FILE *) outFILE,
      "    o minimum number of deleted bases to flag a\n"
   );

   fprintf(
      (FILE *) outFILE,
      "      deletion as a DI event\n"
   );

   /*min nucleotides at start/end of sequence*/
   fprintf(
      (FILE *) outFILE,
      "  -min-nt-pad %u: [Optinal; %u]\n",
      def_minPadNt_mainDiCoords,
      def_minPadNt_mainDiCoords
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
   *   - print settings
   \*****************************************************/

   /*print non-DI reads*/
   if(def_pNonDI_mainDiCoords)
      fprintf(
         (FILE *) outFILE,
         "  -non-di: [Optinal; Yes]\n"
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -non-di: [Optinal; No]\n"
      );

   fprintf(
      (FILE *) outFILE,
      "    o print entries with no DI events\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    o disable with \"-di-only\"\n"
   );

   /*****************************************************\
   * Fun02 Sec02 Sub05:
   *   - help and version number
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  -h: print this help message\n"
   );

   fprintf(
      (FILE *) outFILE,
      "  -v: print the version number\n"
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
      "  - tsv with read id, reference name,\n"
   );

   fprintf(
      (FILE *) outFILE,
      "    DI start coordinate, and DI end coordinate\n"
   );
} /*phelp_mainDiCoords*/

/*-------------------------------------------------------\
| Fun03: input_mainDiCoords
|   - gets user input from an array of c-strings
| Input:
|   - numArgsSI:
|     o number of arguments input
|   - argAryStr:
|     o array of c-strings with user input
|   - samFileStrPtr:
|     o pointer to c-sting to point to sam file name
|   - outFileStrPtr:
|     o pointer to c-sting to point to output file name
|   - minDILenUIPtr:
|     o pointer to unsigned int to hold the minimum
|       deletion length to be a DI event
|   - minPadNtUIPtr:
|     o pointer to unsigned int to hold the minimum
|       number of nucleotides before/after a large DI
|       event
|   - pNoDIBlPtr:
|     o pointer to signed char (1: print non-DI read ids)
| Output:
|   - Modifies:
|     o all variables except numArgsSI and argAryStr to
|       hold/point to user input
|   - Returns:
|     o 0 for no errors
|     o 1 for printed version number or help message
|     o 2 for invalid input
\-------------------------------------------------------*/
signed char
input_mainDiCoords(
   int numArgsSI,                /*number of arguments*/
   char *argAryStr[],            /*has input arguments*/
   signed char **samFileStrPtr,  /*will have sam file*/
   signed char **outFileStrPtr,  /*will have out file*/
   unsigned int *minDILenUIPtr,  /*will have min DI len*/
   unsigned int *minPadNtUIPtr,  /*ends non-del length*/
   signed char *pNoDIBlPtr       /*1 if keeping non-DI*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   '   o fun03 sec01:
   '     - variable declarations
   '   o fun03 sec02:
   '     - check if user input something
   '   o fun03 sec03:
   '     - get user input
   '   o fun03 sec04:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar errSC = 0;
   sint siArg = 1;
   schar *tmpStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - check if user input something
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numArgsSI == 0)
   { /*If: nothing input*/
      phelp_mainDiCoords(stdout);
      goto phelp_fun03_sec04;
   } /*If: nothing input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - get user input
   ^   o fun03 sec03 sub01:
   ^     - start loop and general input
   ^   o fun03 sec03 sub02:
   ^     - check for help messages
   ^   o fun03 sec03 sub03:
   ^     - check for version number requests
   ^   o fun03 sec03 sub04:
   ^     - invalid entry and move to next entry
   \\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - start loop and general input
   \*****************************************************/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/
      if(
         ! eql_charCp(
            (schar *) "-sam",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*If: is the sam file name*/
         ++siArg;
         *samFileStrPtr = (schar *) argAryStr[siArg];
      } /*If: is the sam file name*/

      else if(
         ! eql_charCp(
            (schar *) "-out",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: is the output file name*/
         ++siArg;
         *outFileStrPtr = (schar *) argAryStr[siArg];
      } /*Else If: is the output file name*/

      else if(
         ! eql_charCp(
            (schar *) "-non-di",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pNoDIBlPtr = 1;

      else if(
         ! eql_charCp(
            (schar *) "-di-only",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ) *pNoDIBlPtr = 0;

      else if(
         ! eql_charCp(
            (schar *) "-min-del",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: is the minimum DI deletion lenght*/
         ++siArg;

         tmpStr = (schar *) argAryStr[siArg];

         tmpStr +=
            strToUI_base10str(
               (schar *) argAryStr[siArg],
               minDILenUIPtr
            );

         if(*tmpStr != '\0')
         { /*If: had non-numeric entry*/
            fprintf(
               stderr,
               "-min-del %s is non-numeric or to large\n",
               argAryStr[siArg]
            );
 
            goto invalidInput_fun03_sec04;
         } /*If: had non-numeric entry*/
      } /*Else If: is the minimum DI deletion lenght*/

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
 
            goto invalidInput_fun03_sec04;
         } /*If: had non-numeric entry*/
      } /*Else If: is the minimum DI deletion length*/

      /**************************************************\
      * Fun03 Sec03 Sub02:
      *   - check for help messages
      \**************************************************/

      else if(
         ! eql_charCp(
            (schar *) "-h",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainDiCoords(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eql_charCp(
            (schar *) "--h",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainDiCoords(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eql_charCp(
            (schar *) "help",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainDiCoords(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eql_charCp(
            (schar *) "-help",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainDiCoords(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eql_charCp(
            (schar *) "--help",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainDiCoords(stdout);
         goto phelp_fun03_sec04;
      } /*Else If: help message requested*/

      /**************************************************\
      * Fun03 Sec03 Sub03:
      *   - check for version number requests
      \**************************************************/
      
      else if(
         ! eql_charCp(
            (schar *) "-v",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainDiCoords(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eql_charCp(
            (schar *) "--v",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainDiCoords(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eql_charCp(
            (schar *) "version",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainDiCoords(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eql_charCp(
            (schar *) "-version",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainDiCoords(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eql_charCp(
            (schar *) "--version",
            (schar *) argAryStr[siArg],
            (schar) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainDiCoords(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      /**************************************************\
      * Fun03 Sec03 Sub04:
      *   - invalid entry and move to next entry
      \**************************************************/

      else
      { /*Else: uknown input*/
         fprintf(
            stderr,
            "%s is not recongized\n",
            argAryStr[siArg]
         );
 
         goto invalidInput_fun03_sec04;
      } /*Else: uknown input*/

      ++siArg;
   } /*Loop: get user input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto cleanUp_fun03_sec04;

   phelp_fun03_sec04:;
   errSC = 1;
   goto cleanUp_fun03_sec04;

   pversion_fun03_sec04:;
   errSC = 1;
   goto cleanUp_fun03_sec04;

   invalidInput_fun03_sec04:;
   errSC = 2;
   goto cleanUp_fun03_sec04;

   cleanUp_fun03_sec04:;
   return(errSC);
} /*input_mainDiCoords*/

/*-------------------------------------------------------\
| Main:
|   - main driver function to get DI coordinates
| Input:
|   - numArgsSI:
|     o number of arguments the user input
|   - argAryStr:
|     o array of c-strings with user input
| Output:
|   - Prints:
|     o tsv with ids, number DI events & coordiantes
\-------------------------------------------------------*/
int
main(
   int numArgsSI,
   char *argAryStr[]
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '   - main driver function to get DI coordinates
   '   o main sec01:
   '     - variable declerations
   '   o main sec02:
   '     - initialize structs, get input, setup memory
   '   o main sec03:
   '     - get coordinates
   '   o main sec04:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   schar *samFileStr = 0;
   schar *outFileStr = 0;

   uint minDILenUI = def_minDIDelLen_mainDiCoords;
   uint minPadNtUI = def_minPadNt_mainDiCoords;

   schar pNoDIBl = def_pNonDI_mainDiCoords;
   schar errSC = 0;

   uint *startHeapAryUI = 0;
   uint *endHeapAryUI = 0;
   uint *delSizeHeapAryUI = 0;
   uint lenArysUI = 0;
   sint numDIsSI = 0;

   struct samEntry samStackST;
   schar *buffHeapStr = 0; /*buffer for sam file*/
   ulong lenBuffUL = 0;   /*length of buffHeapStr*/
   ulong numEntryUL = 0;   /*number entries in sam file*/

   FILE *samFILE = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - initialize structs, get input, setup memory
   ^   o main sec02 sub01:
   ^     - initialize structures
   ^   o main sec02 sub02:
   ^     - get input
   ^   o main sec02 sub03:
   ^     - allocate memory
   ^   o main sec02 sub04:
   ^     - open sam file
   ^   o main sec02 sub05:
   ^     - open output file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize structures
   \*****************************************************/

   init_samEntry(&samStackST);

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get input
   \*****************************************************/


   errSC =
      input_mainDiCoords(
         numArgsSI,     /*number of arguments*/
         argAryStr,     /*has input arguments*/
         &samFileStr,   /*will have sam file*/
         &outFileStr,   /*will have out file*/
         &minDILenUI,   /*will have min DI len*/
         &minPadNtUI,   /*min pad length to count DI*/
         &pNoDIBl       /*1 if keeping non-DI*/
     ); /*get user input*/

   if(errSC)
   { /*If: had error*/
      --errSC; /*help message/version error to no error*/
      goto cleanUp_sec04_sub04;
   } /*If: had error*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - allocate memory
   \*****************************************************/

   setup_samEntry(&samStackST);

   /*****************************************************\
   * Main Sec02 Sub04:
   *   - open sam file
   \*****************************************************/

   if(
         ! samFileStr
      || *samFileStr == '-'
   ) samFILE = stdin;

   else
   { /*Else: user provided a sam file*/
      samFILE =
         fopen(
            (char *) samFileStr,
            "r"
         );

      if(! samFILE)
      { /*If: no sam file*/
         fprintf(
            stderr,
            "could not open -sam %s\n",
            samFileStr
         );

         goto fileErr_main_sec04_sub03;
      } /*If: no sam file*/
   } /*Else: user provided a sam file*/

   /*****************************************************\
   * Main Sec02 Sub05:
   *   - open output file
   \*****************************************************/

   if(
         ! outFileStr
      || *outFileStr == '-'
   ) outFILE = stdout;

   else
   { /*Else: user provided an output file*/
      outFILE =
         fopen(
            (char *) outFileStr,
            "r"
         );

      if(! outFILE)
      { /*If: no output file*/
         fprintf(
            stderr,
            "could not open -sam %s\n",
            outFileStr
         );

         goto fileErr_main_sec04_sub03;
      } /*If: no output file*/
   } /*Else: user provided an output file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^   - get coordinates
   ^   o main sec03 sub01:
   ^     - get first sam file entry and start loop
   ^   o main sec03 sub02:
   ^     - remove comments, and non-primary aligned reads
   ^   o main sec03 sub04:
   ^     - move to the next entry in the sam file
   ^   o main sec03 sub05:
   ^     - after loop; check for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec03 Sub01:
   *   - get first sam file entry and start loop
   \*****************************************************/

   pDIHead_diCoords(outFILE);

   errSC =
      (schar)
      get_samEntry(
         &samStackST,
         &buffHeapStr,
         &lenBuffUL,
         samFILE
      ); /*get first sam entry*/

   while(! errSC)
   { /*Loop: find DI events and coordinates*/

      /**************************************************\
      * Main Sec03 Sub02:
      *   - remove comments, and non-primary aligned reads
      \**************************************************/

      ++numEntryUL;

      if(samStackST.extraStr[0] == '@')
         goto nextEntry_main_sec03_sub04;

      if(samStackST.flagUS & (4 | 256 | 2048))
         goto nextEntry_main_sec03_sub04;
         /* Explenation:
         `   - 4 = unmapped read
         `   - 256 = secondary alignment (alternative),
         `   - 2048 = supplemental alignment (diff ref)
         */
          
      /**************************************************\
      * Main Sec03 Sub03:
      *   - scan for DI entries
      \**************************************************/

      numDIsSI =
         get_diCoords(
            &samStackST,
            minDILenUI,
            minPadNtUI,
            &startHeapAryUI,
            &endHeapAryUI,
            &delSizeHeapAryUI,
            &lenArysUI
        ); /*find DI events and their coordinates*/
   
      if(numDIsSI < 0)
      { /*If: memory error*/
         fprintf(
            stderr,
            "MEMORY error on -sam %s line number %lu\n",
            samFileStr,
            numEntryUL
         );

         goto memErr_main_sec04_sub02;
      } /*If: memory error*/

      if(numDIsSI > 0 || pNoDIBl)
      { /*If: printing out read*/
         pDI_diCoords(
            (schar *) samStackST.qryIdStr,
            (schar *) samStackST.refIdStr,
            startHeapAryUI,
            endHeapAryUI,
            delSizeHeapAryUI,
            numDIsSI,
            outFILE
        );
      } /*If: printing out read*/

      /**************************************************\
      * Main Sec03 Sub04:
      *   - move to the next entry in the sam file
      \**************************************************/
      nextEntry_main_sec03_sub04:;

      errSC =
         (schar)
         get_samEntry(
            &samStackST,
            &buffHeapStr,
            &lenBuffUL,
            samFILE
         ); /*get first sam entry*/
   } /*Loop: find DI events and coordinates*/

   /*****************************************************\
   * Main Sec03 Sub05:
   *   - after loop; check for errors
   \*****************************************************/

   if(errSC != 1)
   { /*If: had memory error*/
      fprintf(
         stderr,
         "MEMORY ERROR on line %lu of -sam %s\n",
         numEntryUL,
         samFileStr
      );

      goto memErr_main_sec04_sub02;
   } /*If: had memory error*/

   if(! numEntryUL)
   { /*If: had no sam file entrys*/
      fprintf(
         stderr,
         "-sam %s has no mapped reads\n",
         samFileStr
      );

      goto memErr_main_sec04_sub02;
   } /*If: had no sam file entrys*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - clean up
   ^   o main sec04 sub01:
   ^     - no error clean up
   ^   o main sec04 sub02:
   ^     - memory error clean up
   ^   o main sec04 sub03:
   ^     - file error clean up
   ^   o main sec04 sub04:
   ^     - variable clean up (all errors and no errors)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - no error clean up
   \*****************************************************/

   errSC = 0;
   goto cleanUp_sec04_sub04;

   /*****************************************************\
   * Main Sec04 Sub02:
   *   - memory error clean up
   \*****************************************************/

   memErr_main_sec04_sub02:;
   errSC = -1;
   goto cleanUp_sec04_sub04;

   /*****************************************************\
   * Main Sec04 Sub03:
   *   - file error clean up
   \*****************************************************/

   fileErr_main_sec04_sub03:;
   errSC = 2;
   goto cleanUp_sec04_sub04;

   /*****************************************************\
   * Main Sec04 Sub04:
   *   - variable clean up (all errors and no errors)
   \*****************************************************/

   cleanUp_sec04_sub04:;

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

   if(startHeapAryUI)
      free(startHeapAryUI);

   startHeapAryUI = 0;

   if(endHeapAryUI)
      free(endHeapAryUI);

   endHeapAryUI = 0;

   if(delSizeHeapAryUI)
      free(delSizeHeapAryUI);

   delSizeHeapAryUI = 0;

   if(buffHeapStr)
      free(buffHeapStr);

   buffHeapStr = 0;

   freeStack_samEntry(&samStackST);

   exit(errSC);
} /*main*/

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
