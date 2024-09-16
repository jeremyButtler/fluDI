/*########################################################
# Name alnSetStruct
# Use:
#  o Holds the settings structures and supporting
#    functions for setting structures for alnSeq's
#    pairwise aligners.
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - included libraries
'  o .h st01 alnSet:
'     o Holds settings for my alignment program
'  o .h fun01 setScore_alnSet:
'    - Sets the score for a base pair (reference/query)
'  o .h fun02 setMatch_alnSet:
'    - Sets if two bases are a match or mismtach
'  o fun03 freeStack_alnSet:
'    o Frees variables inside alnSet
'  o fun04 freeHeap_alnSet:
'    o Frees an alnSet structure (and sets to 0)
'  o .h fun05 getScore_alnSet:
'    - Get the score for a pair of bases from an alignment
'  o .h fun06 getMatch_alnSet:
'    - Check if two bases were a match or mismatch
'  o fun07 readScoreFile_alnSet
'     - Reads in a file of scores for a scoring matrix
'  o fun08 readMatchFile_alnSet:
'    - Reads in a file of matches
'  o fun09 seqToIndex_alnSet:
'    - Converts a sequence to a look up table index
'  o fun10 indexToSeq_alnSet:
'    - Converts a sequence of lookup indexs back into
'      uppercase characters (a-z)
'  o fun11: changeGap_alnSet
'    - changes the gap penalties in an alnSet structure
'  o fun12: maxScore_alnSet
'    - finds maximum score possible for a sequence
'  o fun13 init_alnSet:
'    - Set all values in altSet (alingment settings)
'      structure to defaults
'  o fun14: pDefMatchMatrix_alnSet
'    - print out the default match matrix
'  o license:
'    - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "alnSet.h"

#include <stdio.h>

#include "../genLib/base10str.h"

/*no .c files*/
#include "../genLib/dataTypeShortHand.h"
#include "alnDefs.h"

/*-------------------------------------------------------\
| Fun03: freeStack_alnSet
|  - Does a stack free of an alnSet structer
| Input:
|  - alnSetSTPtr:
|    o alnSetSTPtr to free internal variables in
| Output:
|  - Free:
|    o Nothing, there are no heap allocated variables.
|      Is here in case I decide to have heap allocated
|      variables on day.
\-------------------------------------------------------*/
void
freeStack_alnSet(
   struct alnSet *alnSetSTPtr
){
   alnSetSTPtr = alnSetSTPtr; /*quites error message*/
} /*freeStack_alnSet*/

/*-------------------------------------------------------\
| Fun04: freeHeap_alnSet
|  - Frees and alnSet (alignment settings) structure
| Input:
|  - alnSetSTPtr:
|    o Pionter to alnSetST to free
| Output:
|  - Free:
|    o alnSetSTPtr
|  - Set:
|    o alnSetSTPtr to 0 (NULL)
\-------------------------------------------------------*/
void
freeHeap_alnSet(
   struct alnSet *alnSetSTPtr
){
   freeStack_alnSet(alnSetSTPtr);
   free(alnSetSTPtr);
} /*freeHeap_alnSet*/

/*-------------------------------------------------------\
| Fun07: readScoreFile_alnSet
|  - Reads in a file of scores for a scoring matrix
| Input:
|  - alnSetSTPtr:
|    o pointer to an alnSetST (settings) structure with
|      the score matrix to modify
|  - scoreFILE:
|    o File to get scores from
| Output:
|  - Modifies:
|    o Score matrix in alngSetST to hold the scores from
|      the file (scoreFILE)
|    o scoreFILE to point to the end of the file
|  - Returns:
|    o 0 for no errors
|    o position of error in file
\-------------------------------------------------------*/
unsigned long
readScoreFile_alnSet(
    struct alnSet *alnSetSTPtr, /*score matrix to change*/
    void *scoreFILE         /*File scoring matrix scores*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC: readScoreFile_alnSet
   '  o fun07 sec01:
   '    - Variable declerations & set up
   '  o fun07 sec02:
   '    - Blank the scoring matrix
   '  o fun07 sec03:
   '    - Read in line and check if comment
   '  o fun07 sec04:
   '    - Convert score & add to matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01:
   ^  - Variable declerations and set up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define lenBuffUS 1024
   schar buffStr[lenBuffUS];
   schar *tmpStr = 0;
   sshort scoreSS = 0;

   uchar colUC = 0;
   uchar rowUC = 0;

   buffStr[lenBuffUS - 1] = '\0';
   buffStr[lenBuffUS - 2] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^  - Blank the scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < def_lenCol_alnSet; ++colUC)
   { /*Loop: blank all values in the scoring matrix*/
       for(rowUC = 0; rowUC < def_lenCol_alnSet; ++rowUC)
           alnSetSTPtr->scoreMatrixSS[colUC][rowUC] = 0;
   } /*Loop: blank all values in the scoring matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^  - Read in line and check if comment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(fgets((char *) buffStr, 1024, scoreFILE))
   { /*While I have scores to read in*/
       
       if(buffStr[0] == '/' && buffStr[1] == '/')
       { /*On a comment, move onto the next line*/
           while(
               buffStr[lenBuffUS - 2] != '\0' &&
               buffStr[lenBuffUS - 2] != '\n'
           ){ /*Loop: read in more buffer*/
               buffStr[lenBuffUS - 2] = '\0';

               tmpStr =
                  (schar *)
                  fgets(
                     (char *) buffStr,
                     1024,
                     (FILE *) scoreFILE
                  );
           } /*Loop: read in more buffer*/

           /*Reset the buffer*/
           buffStr[lenBuffUS - 2] = '\0';

           continue;
       } /*On a comment, move onto the next line*/

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun07 Sec04:
       ^  - Convert score & add to matrix
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       if(buffStr[0] == '\n')
           continue;                        /*Blank line*/

       if(buffStr[0] < 64 && buffStr[2] < 64)
           return ftell(scoreFILE);  /*Invalid character*/
       
       tmpStr +=
           strToSS_base10str(
              &buffStr[4],
              &scoreSS
           );

       setScore_alnSet(
         buffStr[0],
         buffStr[2],
         scoreSS,
         alnSetSTPtr
       ); /*Add the score to the matrix*/

       if(tmpStr == &buffStr[3])
           return ftell(scoreFILE);         /*No score*/

       while(
           buffStr[lenBuffUS - 2] != '\0' &&
           buffStr[lenBuffUS - 2] != '\n'
       ){ /*While have more buffer to read in*/
           buffStr[lenBuffUS - 2] = '\0';

           tmpStr =
              (schar *)
              fgets(
                 (char *) buffStr,
                 1024,
                 (FILE *) scoreFILE
              );
       } /*While have more buffer to read in*/

       /*Reset the buffer*/
       buffStr[lenBuffUS - 2] = '\0';
   } /*While I have scores to read in*/

   return 0;
} /*readScoreFile_alnSet*/

/*-------------------------------------------------------\
| Fun08: readMatchFile_alnSet
|  - Reads in a file of matches
| Input:
|  - alnSetSTPtr:
|    o pointer to an alnSetST (settings) structure with
|      the match matrix to modify
|  - matchFILE:
|    o File to get matchs from
| Output:
|  - Modifies:
|    o Match matrix in alngSetST to hold the matchs from
|      the file (matchFILE)
|    o matchFILE to point to the end of the file
|  - Returns:
|    o 0 for no errors
|    o position of error in file
\-------------------------------------------------------*/
unsigned long
readMatchFile_alnSet(
    struct alnSet *alnSetSTPtr,
    void *matchFILE  /*File of matchs for scoring matrix*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun08 TOC: readMatchFile_alnSet
   '  o fun08 sec01:
   '    - Variable declerations & set up
   '  o fun08 sec02:
   '    - Blank the match matrix
   '  o fun08 sec03:
   '    - Read in line and check if comment
   '  o fun08 sec04:
   '    - Add match/snp (mismatch) to match matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec01:
   ^  - Variable declerations and set up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define lenBuffUS 1024
   schar buffStr[lenBuffUS];

   uchar colUC = 0;
   uchar rowUC = 0;

   buffStr[lenBuffUS - 1] = '\0';
   buffStr[lenBuffUS - 2] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec02:
   ^  - Blank the scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < def_lenCol_alnSet; ++colUC)
   { /*Loop: blank all values in the scoring matrix*/
       for(rowUC = 0; rowUC < def_lenCol_alnSet; ++rowUC)
           alnSetSTPtr->matchMatrixSC[colUC][rowUC] = 0;
   } /*Loop: blank all values in the scoring matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun08 Sec03:
   ^  - Read in line and check if comment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(
      fgets((char *) buffStr, 1024, (FILE *) matchFILE)
   ){ /*While I have matchs to read in*/
       
       if(buffStr[0] == '/' && buffStr[1] == '/')
       { /*On a comment, move onto the next line*/
           while(
               buffStr[lenBuffUS - 2] != '\0' &&
               buffStr[lenBuffUS - 2] != '\n'
           ) { /*While have more buffer to read in*/
               buffStr[lenBuffUS - 2] = '\0';

               fgets(
                  (char *) buffStr,
                  1024,
                  (FILE *) matchFILE
               );
           } /*While have more buffer to read in*/

           /*Reset the buffer*/
           buffStr[lenBuffUS - 2] = '\0';

           continue;
       } /*On a comment, move onto the next line*/

       /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
       ^ Fun08 Sec04:
       ^  - Convert match & add to matrix
       \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

       if(buffStr[0] == '\n')
           continue;                        /*Blank line*/

       if(buffStr[4] != '1' && buffStr[4] != '0')
           return ftell((FILE *) matchFILE);
           /*This error means I do not know if match/snp*/

       setMatch_alnSet(
         buffStr[0],
         buffStr[2],
         buffStr[4],
         alnSetSTPtr
       ); /*Add the match to the matrix*/

       while(
           buffStr[lenBuffUS - 2] != '\0' &&
           buffStr[lenBuffUS - 2] != '\n'
       ){ /*While have more buffer to read in*/
           buffStr[lenBuffUS - 2] = '\0';

           fgets(
              (char *) buffStr,
              1024,
              (FILE *) matchFILE
           );
       } /*While have more buffer to read in*/

       /*Reset the buffer*/
       buffStr[lenBuffUS - 2] = '\0';
   } /*While I have matchs to read in*/

   return 0;
} /*readMatchFile_alnSet*/

/*-------------------------------------------------------\
| Fun09: seqToIndex_alnSet
|  - Converts a sequence to a look up table index
|    (table is in alnSetStruct.c/h)
| Input:
|  - seqStr:
|    o pointer c-string with the sequence to convert
| Output:
|  - Modifies:
|    o seqST->seqStr to have look up table indexs (1-27, 
|      with null as 0) instead of bases
\-------------------------------------------------------*/
void
seqToIndex_alnSet(
   signed char *seqStr
){
   while(*seqStr != '\0')
   { /*Loop: convert bases to lookup table values*/
      *seqStr &= def_ntToCode_alnSet;
      ++seqStr;
   } /*Loop: convert bases to lookup table values*/
} /*seqToIndex_alnSet*/

/*-------------------------------------------------------\
| Fun10: indexToSeq_alnSet
|  - Converts a sequence of lookup indexs back into
|    uppercase characters (a-z)
| Input:
|  - seqStructPtr:
|    o Pointer to sequence structer with converte sequence
|      to deconvert (lookup index to base)
| Output:
|  - Modifies:
|    o seqST->seqStr to have bases instead of look up
|      table indexs
\-------------------------------------------------------*/
void
indexToSeq_alnSet(
   signed char *seqStr
){
   while(*seqStr != '\0')
   { /*Loop: convert bases to lookup table values*/
      *seqStr |= 64;
      ++seqStr;
   } /*Loop: convert bases to lookup table values*/
} /*indexToSeq_alnSet*/

/*-------------------------------------------------------\
| Fun11: changeGap_alnSet
|   - changes the gap penalties in an alnSet structure
| Input:
|   - alnSetSTPtr:
|     o pointer to alnSet structure to update gap scores
|   - gapSS:
|     o new gap opening penalty
|   - extendSS:
|     o new gap extension penalty
| Output:
|   - Modifies:
|     o delArySS to have gap extension/opening penalties
|     o insArySS to have gap extension/opening penalties
|     o gapSS in alnSetSTPtr to have gap opening penalty
|     o extendSS to have gap extension penalty
\-------------------------------------------------------*/
void
changeGap_alnSet(
   struct alnSet *alnSetSTPtr,
   signed short gapSS,
   signed short extendSS
){

   alnSetSTPtr->gapSS = gapSS;
   alnSetSTPtr->extendSS = extendSS;
   
   alnSetSTPtr->delArySS[def_mvStop_alnDefs] = gapSS;
   alnSetSTPtr->delArySS[def_mvDel_alnDefs] = extendSS;
   alnSetSTPtr->delArySS[def_mvIns_alnDefs] = gapSS;
   alnSetSTPtr->delArySS[def_mvSnp_alnDefs] = gapSS;

   /*insetion gap penalty matrix*/
   alnSetSTPtr->insArySS[def_mvStop_alnDefs] = gapSS;
   alnSetSTPtr->insArySS[def_mvDel_alnDefs] = gapSS;
   alnSetSTPtr->insArySS[def_mvIns_alnDefs] = extendSS;
   alnSetSTPtr->insArySS[def_mvSnp_alnDefs] = gapSS;
} /*changeGap_alnSet*/

/*-------------------------------------------------------\
| Fun12: maxScore_alnSet
|  - finds maximum score possible for a sequence
| Input:
|  - seqStr:
|    o c-string with sequence to find max score for
|    o needs to be converted to index with
|      indexToSeq_alnSet (prepared for alignment)
|  - startUL:
|    o index of first base to check (index 0)
|  - lenSeqUL:
|    o length of the sequence to check
|  - alnSetPtr:
|    o pointer to an alnSet (settings) structure to
|      initialize
| Output:
|  o Returns:
|    - maximum possible score for alignment
\-------------------------------------------------------*/
signed long
maxScore_alnSet(
   signed char *seqStr,
   unsigned long startUL,
   unsigned long lenSeqUL,
   struct alnSet *alnSetPtr
){
   slong scoreSL = 0;

   while(startUL < lenSeqUL)
   { /*Loop: align sequences*/
      scoreSL +=
         getScore_alnSet(
            seqStr[startUL],
            seqStr[startUL],
            alnSetPtr
         );

      ++startUL;
   } /*Loop: align sequences*/

   return scoreSL;
} /*maxScore_alnSet*/

/*-------------------------------------------------------\
| Fun13: init_alnSet
|  - Set values in altSet (alingment settings) structure
|    to default values
| Input:
|  - alnSetSTPtr:
|    o poineter to an alnSet (settings) structure to
|      initialize
| Output:
|  o Modifies:
|    - alnSetST to have default alignment settings values
\-------------------------------------------------------*/
void
init_alnSet(
    struct alnSet *alnSetST /*Has settings to initialize*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun13 TOC: init_alnSet
   '  - Set values in altSet (alingment settings)
   '    structure to defaults
   '  o fun13 sec01:
   '    - Set non-matrix variables
   '  o fun13 sec02:
   '    - Initialize scoring matrix
   '  o fun13 sec03:
   '    - Initialize match matrix
   '  o fun13 sec04:
   '    - set up scoring matrix for nucleotides
   '  o fun13 sec05:
   '    - set up matching matrix for nucleotides
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec01:
   ^  - Set non-matrix variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Variables for my for loop*/
   uchar colUC = 0;
   uchar rowUC = 0;


   changeGap_alnSet(
      alnSetST,
      def_gapOpen_alnDefs,
      def_gapExtend_alnDefs
   ); /*add in gap scoring*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec02:
   ^  - Initialize scoring matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(colUC = 0; colUC < def_lenCol_alnSet; ++colUC)
   { /*loop for all columns in the comparison matrix*/
       for(rowUC = 0; rowUC < def_lenCol_alnSet; ++rowUC)
           alnSetST->scoreMatrixSS[colUC][rowUC] = 0;
           /*Most of these cells will never be used*/
           /*But are needed to build the table*/
   } /*loop for all columns in the comparison matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec03:
   ^  - Initialize match matrix
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Both words, DNA, and AA all are the same when both
    ` characters/bases/amino acids are the same.
    */
    for(colUC = 0; colUC < def_lenCol_alnSet; ++colUC)
    { /*Loop: fill in comparison matrix*/
        for(rowUC = 0; rowUC < def_lenCol_alnSet; ++rowUC)
        { /*Loop: Fill in the entire matrix*/ 
           if(colUC == rowUC)
             alnSetST->matchMatrixSC[colUC][rowUC] =
                def_ntEql_alnDefs;
           else
             alnSetST->matchMatrixSC[colUC][rowUC] =
                def_ntNotEql_alnDefs;
        } /*Loop: Fill in the entire matrix*/ 
    } /*Loop: fill in comparison matrix*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec04:
   ^   - set up scoring matrix for nucleotides
   ^   o fun13 sec04 sub01:
   ^     - a as first base
   ^   o fun13 sec04 sub02:
   ^     - t as first base
   ^   o fun13 sec04 sub03:
   ^     - u (t) as first base
   ^   o fun13 sec04 sub04:
   ^     - g as first base
   ^   o fun13 sec04 sub05:
   ^     - c as first base
   ^   o fun13 sec04 sub06:
   ^     - w (anonymous) as first base
   ^   o fun13 sec04 sub07:
   ^     - s (anonymous) as first base
   ^   o fun13 sec04 sub08:
   ^     - m (anonymous) as first base
   ^   o fun13 sec04 sub09:
   ^     - k (anonymous) as first base
   ^   o fun13 sec04 sub10:
   ^     - r (anonymous) as first base
   ^   o fun13 sec04 sub11:
   ^     - y (anonymous) as first base
   ^   o fun13 sec04 sub12:
   ^     - b (anonymous) as first base
   ^   o fun13 sec04 sub13:
   ^     - d (anonymous) as first base
   ^   o fun13 sec04 sub14:
   ^     - h (anonymous) as first base
   ^   o fun13 sec04 sub15:
   ^     - v (anonymous) as first base
   ^   o fun13 sec04 sub16:
   ^     - n (anonymous) as first base
   ^   o fun13 sec04 sub17:
   ^     - x (anonymous) as first base (technically aa)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun13 Sec04 Sub01:
   *   - a as first base
   \*****************************************************/

   setScore_alnSet('a','a',def_AToA_alnDefs,alnSetST);
   setScore_alnSet('a','t',def_AToT_alnDefs,alnSetST);
   setScore_alnSet('a','u',def_AToU_alnDefs,alnSetST);
   setScore_alnSet('a','g',def_AToG_alnDefs,alnSetST);
   setScore_alnSet('a','c',def_AToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('a','w',def_AToW_alnDefs,alnSetST);
   setScore_alnSet('a','s',def_AToS_alnDefs,alnSetST);
   setScore_alnSet('a','m',def_AToM_alnDefs,alnSetST);
   setScore_alnSet('a','k',def_AToK_alnDefs,alnSetST);
   setScore_alnSet('a','r',def_AToR_alnDefs,alnSetST);
   setScore_alnSet('a','y',def_AToY_alnDefs,alnSetST);
   setScore_alnSet('a','b',def_AToB_alnDefs,alnSetST);
   setScore_alnSet('a','d',def_AToD_alnDefs,alnSetST);
   setScore_alnSet('a','h',def_AToH_alnDefs,alnSetST);
   setScore_alnSet('a','v',def_AToV_alnDefs,alnSetST);
   setScore_alnSet('a','n',def_AToN_alnDefs,alnSetST);
   setScore_alnSet('a','x',def_AToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub02:
   *   - t as first base
   \*****************************************************/

   setScore_alnSet('t','a',def_TToA_alnDefs,alnSetST);
   setScore_alnSet('t','t',def_TToT_alnDefs,alnSetST);
   setScore_alnSet('t','u',def_TToU_alnDefs,alnSetST);
   setScore_alnSet('t','g',def_TToG_alnDefs,alnSetST);
   setScore_alnSet('t','c',def_TToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('t','w',def_TToW_alnDefs,alnSetST);
   setScore_alnSet('t','s',def_TToS_alnDefs,alnSetST);
   setScore_alnSet('t','m',def_TToM_alnDefs,alnSetST);
   setScore_alnSet('t','k',def_TToK_alnDefs,alnSetST);
   setScore_alnSet('t','r',def_TToR_alnDefs,alnSetST);
   setScore_alnSet('t','y',def_TToY_alnDefs,alnSetST);
   setScore_alnSet('t','b',def_TToB_alnDefs,alnSetST);
   setScore_alnSet('t','d',def_TToD_alnDefs,alnSetST);
   setScore_alnSet('t','h',def_TToH_alnDefs,alnSetST);
   setScore_alnSet('t','v',def_TToV_alnDefs,alnSetST);
   setScore_alnSet('t','n',def_TToN_alnDefs,alnSetST);
   setScore_alnSet('t','x',def_TToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub03:
   *   - u (t) as first base
   \*****************************************************/

   setScore_alnSet('u','a',def_UToA_alnDefs,alnSetST);
   setScore_alnSet('u','t',def_UToT_alnDefs,alnSetST);
   setScore_alnSet('u','u',def_UToU_alnDefs,alnSetST);
   setScore_alnSet('u','g',def_UToG_alnDefs,alnSetST);
   setScore_alnSet('u','c',def_UToC_alnDefs,alnSetST);

   /*Set u & t to same scores (U is RNA version of T)*/
   setScore_alnSet('u','w',def_UToW_alnDefs,alnSetST);
   setScore_alnSet('u','s',def_UToS_alnDefs,alnSetST);
   setScore_alnSet('u','m',def_UToM_alnDefs,alnSetST);
   setScore_alnSet('u','k',def_UToK_alnDefs,alnSetST);
   setScore_alnSet('u','r',def_UToR_alnDefs,alnSetST);
   setScore_alnSet('u','y',def_UToY_alnDefs,alnSetST);
   setScore_alnSet('u','b',def_UToB_alnDefs,alnSetST);
   setScore_alnSet('u','d',def_UToD_alnDefs,alnSetST);
   setScore_alnSet('u','h',def_UToH_alnDefs,alnSetST);
   setScore_alnSet('u','v',def_UToV_alnDefs,alnSetST);
   setScore_alnSet('u','n',def_UToN_alnDefs,alnSetST);
   setScore_alnSet('u','x',def_UToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub04:
   *   - g as first base
   \*****************************************************/

   setScore_alnSet('g','a',def_GToA_alnDefs,alnSetST);
   setScore_alnSet('g','t',def_GToT_alnDefs,alnSetST);
   setScore_alnSet('g','u',def_GToU_alnDefs,alnSetST);
   setScore_alnSet('g','g',def_GToG_alnDefs,alnSetST);
   setScore_alnSet('g','c',def_GToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('g','w',def_GToW_alnDefs,alnSetST);
   setScore_alnSet('g','s',def_GToS_alnDefs,alnSetST);
   setScore_alnSet('g','m',def_GToM_alnDefs,alnSetST);
   setScore_alnSet('g','k',def_GToK_alnDefs,alnSetST);
   setScore_alnSet('g','r',def_GToR_alnDefs,alnSetST);
   setScore_alnSet('g','y',def_GToY_alnDefs,alnSetST);
   setScore_alnSet('g','b',def_GToB_alnDefs,alnSetST);
   setScore_alnSet('g','d',def_GToD_alnDefs,alnSetST);
   setScore_alnSet('g','h',def_GToH_alnDefs,alnSetST);
   setScore_alnSet('g','v',def_GToV_alnDefs,alnSetST);
   setScore_alnSet('g','n',def_GToN_alnDefs,alnSetST);
   setScore_alnSet('g','x',def_GToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub05:
   *   - c as first base
   \*****************************************************/

   setScore_alnSet('c','a',def_CToA_alnDefs,alnSetST);
   setScore_alnSet('c','t',def_CToT_alnDefs,alnSetST);
   setScore_alnSet('c','u',def_CToT_alnDefs,alnSetST);
   setScore_alnSet('c','g',def_CToG_alnDefs,alnSetST);
   setScore_alnSet('c','c',def_CToC_alnDefs,alnSetST);

   /*anonymous matches*/
   setScore_alnSet('c','w',def_CToW_alnDefs,alnSetST);
   setScore_alnSet('c','s',def_CToS_alnDefs,alnSetST);
   setScore_alnSet('c','m',def_CToM_alnDefs,alnSetST);
   setScore_alnSet('c','k',def_CToK_alnDefs,alnSetST);
   setScore_alnSet('c','r',def_CToR_alnDefs,alnSetST);
   setScore_alnSet('c','y',def_CToY_alnDefs,alnSetST);
   setScore_alnSet('c','b',def_CToB_alnDefs,alnSetST);
   setScore_alnSet('c','d',def_CToD_alnDefs,alnSetST);
   setScore_alnSet('c','h',def_CToH_alnDefs,alnSetST);
   setScore_alnSet('c','v',def_CToV_alnDefs,alnSetST);
   setScore_alnSet('c','n',def_CToN_alnDefs,alnSetST);
   setScore_alnSet('c','x',def_CToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub06:
   *   - w (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('w','a',def_WToA_alnDefs,alnSetST);
   setScore_alnSet('w','c',def_WToC_alnDefs,alnSetST);
   setScore_alnSet('w','g',def_WToG_alnDefs,alnSetST);
   setScore_alnSet('w','t',def_WToT_alnDefs,alnSetST);
   setScore_alnSet('w','u',def_WToT_alnDefs,alnSetST);

   setScore_alnSet('w','w',def_WToW_alnDefs,alnSetST);
   setScore_alnSet('w','s',def_WToS_alnDefs,alnSetST);
   setScore_alnSet('w','m',def_WToM_alnDefs,alnSetST);
   setScore_alnSet('w','k',def_WToK_alnDefs,alnSetST);
   setScore_alnSet('w','r',def_WToR_alnDefs,alnSetST);
   setScore_alnSet('w','y',def_WToY_alnDefs,alnSetST);
   setScore_alnSet('w','b',def_WToB_alnDefs,alnSetST);
   setScore_alnSet('w','d',def_WToD_alnDefs,alnSetST);
   setScore_alnSet('w','h',def_WToH_alnDefs,alnSetST);
   setScore_alnSet('w','v',def_WToV_alnDefs,alnSetST);
   setScore_alnSet('w','n',def_WToN_alnDefs,alnSetST);
   setScore_alnSet('w','x',def_WToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub07:
   *   - s (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('s','a',def_SToA_alnDefs,alnSetST);
   setScore_alnSet('s','c',def_SToC_alnDefs,alnSetST);
   setScore_alnSet('s','g',def_SToG_alnDefs,alnSetST);
   setScore_alnSet('s','t',def_SToT_alnDefs,alnSetST);
   setScore_alnSet('s','u',def_SToT_alnDefs,alnSetST);

   setScore_alnSet('s','w',def_SToW_alnDefs,alnSetST);
   setScore_alnSet('s','s',def_SToS_alnDefs,alnSetST);
   setScore_alnSet('s','m',def_SToM_alnDefs,alnSetST);
   setScore_alnSet('s','k',def_SToK_alnDefs,alnSetST);
   setScore_alnSet('s','r',def_SToR_alnDefs,alnSetST);
   setScore_alnSet('s','y',def_SToY_alnDefs,alnSetST);
   setScore_alnSet('s','b',def_SToB_alnDefs,alnSetST);
   setScore_alnSet('s','d',def_SToD_alnDefs,alnSetST);
   setScore_alnSet('s','h',def_SToH_alnDefs,alnSetST);
   setScore_alnSet('s','v',def_SToV_alnDefs,alnSetST);
   setScore_alnSet('s','n',def_SToN_alnDefs,alnSetST);
   setScore_alnSet('s','x',def_SToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub08:
   *   - m (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('m','a',def_MToA_alnDefs,alnSetST);
   setScore_alnSet('m','c',def_MToC_alnDefs,alnSetST);
   setScore_alnSet('m','g',def_MToG_alnDefs,alnSetST);
   setScore_alnSet('m','t',def_MToT_alnDefs,alnSetST);
   setScore_alnSet('m','u',def_MToT_alnDefs,alnSetST);

   setScore_alnSet('m','w',def_MToW_alnDefs,alnSetST);
   setScore_alnSet('m','s',def_MToS_alnDefs,alnSetST);
   setScore_alnSet('m','m',def_MToM_alnDefs,alnSetST);
   setScore_alnSet('m','k',def_MToK_alnDefs,alnSetST);
   setScore_alnSet('m','r',def_MToR_alnDefs,alnSetST);
   setScore_alnSet('m','y',def_MToY_alnDefs,alnSetST);
   setScore_alnSet('m','b',def_MToB_alnDefs,alnSetST);
   setScore_alnSet('m','d',def_MToD_alnDefs,alnSetST);
   setScore_alnSet('m','h',def_MToH_alnDefs,alnSetST);
   setScore_alnSet('m','v',def_MToV_alnDefs,alnSetST);
   setScore_alnSet('m','n',def_MToN_alnDefs,alnSetST);
   setScore_alnSet('m','x',def_MToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub09:
   *   - k (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('k','a',def_KToA_alnDefs,alnSetST);
   setScore_alnSet('k','c',def_KToC_alnDefs,alnSetST);
   setScore_alnSet('k','g',def_KToG_alnDefs,alnSetST);
   setScore_alnSet('k','t',def_KToT_alnDefs,alnSetST);
   setScore_alnSet('k','u',def_KToT_alnDefs,alnSetST);

   setScore_alnSet('k','w',def_KToW_alnDefs,alnSetST);
   setScore_alnSet('k','s',def_KToS_alnDefs,alnSetST);
   setScore_alnSet('k','m',def_KToM_alnDefs,alnSetST);
   setScore_alnSet('k','k',def_KToK_alnDefs,alnSetST);
   setScore_alnSet('k','r',def_KToR_alnDefs,alnSetST);
   setScore_alnSet('k','y',def_KToY_alnDefs,alnSetST);
   setScore_alnSet('k','b',def_KToB_alnDefs,alnSetST);
   setScore_alnSet('k','d',def_KToD_alnDefs,alnSetST);
   setScore_alnSet('k','h',def_KToH_alnDefs,alnSetST);
   setScore_alnSet('k','v',def_KToV_alnDefs,alnSetST);
   setScore_alnSet('k','n',def_KToN_alnDefs,alnSetST);
   setScore_alnSet('k','x',def_KToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub10:
   *   - r (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('r','a',def_RToA_alnDefs,alnSetST);
   setScore_alnSet('r','c',def_RToC_alnDefs,alnSetST);
   setScore_alnSet('r','g',def_RToG_alnDefs,alnSetST);
   setScore_alnSet('r','t',def_RToT_alnDefs,alnSetST);
   setScore_alnSet('r','u',def_RToT_alnDefs,alnSetST);

   setScore_alnSet('r','w',def_RToW_alnDefs,alnSetST);
   setScore_alnSet('r','s',def_RToS_alnDefs,alnSetST);
   setScore_alnSet('r','m',def_RToM_alnDefs,alnSetST);
   setScore_alnSet('r','k',def_RToK_alnDefs,alnSetST);
   setScore_alnSet('r','r',def_RToR_alnDefs,alnSetST);
   setScore_alnSet('r','y',def_RToY_alnDefs,alnSetST);
   setScore_alnSet('r','b',def_RToB_alnDefs,alnSetST);
   setScore_alnSet('r','d',def_RToD_alnDefs,alnSetST);
   setScore_alnSet('r','h',def_RToH_alnDefs,alnSetST);
   setScore_alnSet('r','v',def_RToV_alnDefs,alnSetST);
   setScore_alnSet('r','n',def_RToN_alnDefs,alnSetST);
   setScore_alnSet('r','x',def_RToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub11:
   *   - y (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('y','a',def_YToA_alnDefs,alnSetST);
   setScore_alnSet('y','c',def_YToC_alnDefs,alnSetST);
   setScore_alnSet('y','g',def_YToG_alnDefs,alnSetST);
   setScore_alnSet('y','t',def_YToT_alnDefs,alnSetST);
   setScore_alnSet('y','u',def_YToT_alnDefs,alnSetST);

   setScore_alnSet('y','w',def_YToW_alnDefs,alnSetST);
   setScore_alnSet('y','s',def_YToS_alnDefs,alnSetST);
   setScore_alnSet('y','m',def_YToM_alnDefs,alnSetST);
   setScore_alnSet('y','k',def_YToK_alnDefs,alnSetST);
   setScore_alnSet('y','r',def_YToR_alnDefs,alnSetST);
   setScore_alnSet('y','y',def_YToY_alnDefs,alnSetST);
   setScore_alnSet('y','b',def_YToB_alnDefs,alnSetST);
   setScore_alnSet('y','d',def_YToD_alnDefs,alnSetST);
   setScore_alnSet('y','h',def_YToH_alnDefs,alnSetST);
   setScore_alnSet('y','v',def_YToV_alnDefs,alnSetST);
   setScore_alnSet('y','n',def_YToN_alnDefs,alnSetST);
   setScore_alnSet('y','x',def_YToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub12:
   *   - b (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('b','a',def_BToA_alnDefs,alnSetST);
   setScore_alnSet('b','c',def_BToC_alnDefs,alnSetST);
   setScore_alnSet('b','g',def_BToG_alnDefs,alnSetST);
   setScore_alnSet('b','t',def_BToT_alnDefs,alnSetST);
   setScore_alnSet('b','u',def_BToT_alnDefs,alnSetST);

   setScore_alnSet('b','w',def_BToW_alnDefs,alnSetST);
   setScore_alnSet('b','s',def_BToS_alnDefs,alnSetST);
   setScore_alnSet('b','m',def_BToM_alnDefs,alnSetST);
   setScore_alnSet('b','k',def_BToK_alnDefs,alnSetST);
   setScore_alnSet('b','r',def_BToR_alnDefs,alnSetST);
   setScore_alnSet('b','y',def_BToY_alnDefs,alnSetST);
   setScore_alnSet('b','b',def_BToB_alnDefs,alnSetST);
   setScore_alnSet('b','d',def_BToD_alnDefs,alnSetST);
   setScore_alnSet('b','h',def_BToH_alnDefs,alnSetST);
   setScore_alnSet('b','v',def_BToV_alnDefs,alnSetST);
   setScore_alnSet('b','n',def_BToN_alnDefs,alnSetST);
   setScore_alnSet('b','x',def_BToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub13:
   *   - d (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('d','a',def_DToA_alnDefs,alnSetST);
   setScore_alnSet('d','c',def_DToC_alnDefs,alnSetST);
   setScore_alnSet('d','g',def_DToG_alnDefs,alnSetST);
   setScore_alnSet('d','t',def_DToT_alnDefs,alnSetST);
   setScore_alnSet('d','u',def_DToT_alnDefs,alnSetST);

   setScore_alnSet('d','w',def_DToW_alnDefs,alnSetST);
   setScore_alnSet('d','s',def_DToS_alnDefs,alnSetST);
   setScore_alnSet('d','m',def_DToM_alnDefs,alnSetST);
   setScore_alnSet('d','k',def_DToK_alnDefs,alnSetST);
   setScore_alnSet('d','r',def_DToR_alnDefs,alnSetST);
   setScore_alnSet('d','y',def_DToY_alnDefs,alnSetST);
   setScore_alnSet('d','b',def_DToB_alnDefs,alnSetST);
   setScore_alnSet('d','d',def_DToD_alnDefs,alnSetST);
   setScore_alnSet('d','h',def_DToH_alnDefs,alnSetST);
   setScore_alnSet('d','v',def_DToV_alnDefs,alnSetST);
   setScore_alnSet('d','n',def_DToN_alnDefs,alnSetST);
   setScore_alnSet('d','x',def_DToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub14:
   *   - h (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('h','a',def_HToA_alnDefs,alnSetST);
   setScore_alnSet('h','c',def_HToC_alnDefs,alnSetST);
   setScore_alnSet('h','g',def_HToG_alnDefs,alnSetST);
   setScore_alnSet('h','t',def_HToT_alnDefs,alnSetST);
   setScore_alnSet('h','u',def_HToT_alnDefs,alnSetST);

   setScore_alnSet('h','w',def_HToW_alnDefs,alnSetST);
   setScore_alnSet('h','s',def_HToS_alnDefs,alnSetST);
   setScore_alnSet('h','m',def_HToM_alnDefs,alnSetST);
   setScore_alnSet('h','k',def_HToK_alnDefs,alnSetST);
   setScore_alnSet('h','r',def_HToR_alnDefs,alnSetST);
   setScore_alnSet('h','y',def_HToY_alnDefs,alnSetST);
   setScore_alnSet('h','b',def_HToB_alnDefs,alnSetST);
   setScore_alnSet('h','d',def_HToD_alnDefs,alnSetST);
   setScore_alnSet('h','h',def_HToH_alnDefs,alnSetST);
   setScore_alnSet('h','v',def_HToV_alnDefs,alnSetST);
   setScore_alnSet('h','n',def_HToN_alnDefs,alnSetST);
   setScore_alnSet('h','x',def_HToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub15:
   *   - v (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('v','a',def_VToA_alnDefs,alnSetST);
   setScore_alnSet('v','c',def_VToC_alnDefs,alnSetST);
   setScore_alnSet('v','g',def_VToG_alnDefs,alnSetST);
   setScore_alnSet('v','t',def_VToT_alnDefs,alnSetST);
   setScore_alnSet('v','u',def_VToT_alnDefs,alnSetST);

   setScore_alnSet('v','w',def_VToW_alnDefs,alnSetST);
   setScore_alnSet('v','s',def_VToS_alnDefs,alnSetST);
   setScore_alnSet('v','m',def_VToM_alnDefs,alnSetST);
   setScore_alnSet('v','k',def_VToK_alnDefs,alnSetST);
   setScore_alnSet('v','r',def_VToR_alnDefs,alnSetST);
   setScore_alnSet('v','y',def_VToY_alnDefs,alnSetST);
   setScore_alnSet('v','b',def_VToB_alnDefs,alnSetST);
   setScore_alnSet('v','d',def_VToD_alnDefs,alnSetST);
   setScore_alnSet('v','h',def_VToH_alnDefs,alnSetST);
   setScore_alnSet('v','v',def_VToV_alnDefs,alnSetST);
   setScore_alnSet('v','n',def_VToN_alnDefs,alnSetST);
   setScore_alnSet('v','x',def_VToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub16:
   *   - n (anonymous) as first base
   \*****************************************************/

   setScore_alnSet('n','a',def_NToA_alnDefs,alnSetST);
   setScore_alnSet('n','c',def_NToC_alnDefs,alnSetST);
   setScore_alnSet('n','g',def_NToG_alnDefs,alnSetST);
   setScore_alnSet('n','t',def_NToT_alnDefs,alnSetST);
   setScore_alnSet('n','u',def_NToT_alnDefs,alnSetST);

   setScore_alnSet('n','w',def_NToW_alnDefs,alnSetST);
   setScore_alnSet('n','s',def_NToS_alnDefs,alnSetST);
   setScore_alnSet('n','m',def_NToM_alnDefs,alnSetST);
   setScore_alnSet('n','k',def_NToK_alnDefs,alnSetST);
   setScore_alnSet('n','r',def_NToR_alnDefs,alnSetST);
   setScore_alnSet('n','y',def_NToY_alnDefs,alnSetST);
   setScore_alnSet('n','b',def_NToB_alnDefs,alnSetST);
   setScore_alnSet('n','d',def_NToD_alnDefs,alnSetST);
   setScore_alnSet('n','h',def_NToH_alnDefs,alnSetST);
   setScore_alnSet('n','v',def_NToV_alnDefs,alnSetST);
   setScore_alnSet('n','n',def_NToN_alnDefs,alnSetST);
   setScore_alnSet('n','x',def_NToX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec04 Sub17:
   *   - x (anonymous) as first base (technically aa)
   \*****************************************************/

   setScore_alnSet('x','a',def_XToA_alnDefs,alnSetST);
   setScore_alnSet('x','c',def_XToC_alnDefs,alnSetST);
   setScore_alnSet('x','g',def_XToG_alnDefs,alnSetST);
   setScore_alnSet('x','t',def_XToT_alnDefs,alnSetST);
   setScore_alnSet('x','u',def_XToT_alnDefs,alnSetST);

   setScore_alnSet('x','w',def_XToW_alnDefs,alnSetST);
   setScore_alnSet('x','s',def_XToS_alnDefs,alnSetST);
   setScore_alnSet('x','m',def_XToM_alnDefs,alnSetST);
   setScore_alnSet('x','k',def_XToK_alnDefs,alnSetST);
   setScore_alnSet('x','r',def_XToR_alnDefs,alnSetST);
   setScore_alnSet('x','y',def_XToY_alnDefs,alnSetST);
   setScore_alnSet('x','b',def_XToB_alnDefs,alnSetST);
   setScore_alnSet('x','d',def_XToD_alnDefs,alnSetST);
   setScore_alnSet('x','h',def_XToH_alnDefs,alnSetST);
   setScore_alnSet('x','v',def_XToV_alnDefs,alnSetST);
   setScore_alnSet('x','n',def_XToN_alnDefs,alnSetST);
   setScore_alnSet('x','x',def_XToX_alnDefs,alnSetST);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec05:
   ^  - set up matching matrix for nucleotides
   ^  o fun13 sec05 sub01:
   ^    - a as first base
   ^  o fun13 sec05 sub02:
   ^    - t as first base
   ^  o fun13 sec05 sub03:
   ^    - u (t) as first base
   ^  o fun13 sec05 sub04:
   ^    - g as first base
   ^  o fun13 sec05 sub05:
   ^    - c as first base
   ^  o fun13 sec05 sub06:
   ^    - w (anonymous) as first base
   ^  o fun13 sec05 sub07:
   ^    - s (anonymous) as first base
   ^  o fun13 sec05 sub08:
   ^    - m (anonymous) as first base
   ^  o fun13 sec05 sub09:
   ^    - k (anonymous) as first base
   ^  o fun13 sec05 sub10:
   ^    - r (anonymous) as first base
   ^  o fun13 sec05 sub11:
   ^    - y (anonymous) as first base
   ^  o fun13 sec05 sub12:
   ^    - b (anonymous) as first base
   ^  o fun13 sec05 sub13:
   ^    - d (anonymous) as first base
   ^  o fun13 sec05 sub14:
   ^    - h (anonymous) as first base
   ^  o fun13 sec05 sub15:
   ^    - v (anonymous) as first base
   ^  o fun13 sec05 sub16:
   ^    - n (anonymous) as first base
   ^  o fun13 sec05 sub17:
   ^    - x (anonymous) as first base (technically aa)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun13 Sec05 Sub01:
   *   - a as first base
   \*****************************************************/

   setMatch_alnSet('a','a',def_AEqlA_alnDefs,alnSetST);
   setMatch_alnSet('a','t',def_AEqlT_alnDefs,alnSetST);
   setMatch_alnSet('a','u',def_AEqlU_alnDefs,alnSetST);
   setMatch_alnSet('a','g',def_AEqlG_alnDefs,alnSetST);
   setMatch_alnSet('a','c',def_AEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('a','w',def_AEqlW_alnDefs,alnSetST);
   setMatch_alnSet('a','s',def_AEqlS_alnDefs,alnSetST);
   setMatch_alnSet('a','m',def_AEqlM_alnDefs,alnSetST);
   setMatch_alnSet('a','k',def_AEqlK_alnDefs,alnSetST);
   setMatch_alnSet('a','r',def_AEqlR_alnDefs,alnSetST);
   setMatch_alnSet('a','y',def_AEqlY_alnDefs,alnSetST);
   setMatch_alnSet('a','b',def_AEqlB_alnDefs,alnSetST);
   setMatch_alnSet('a','d',def_AEqlD_alnDefs,alnSetST);
   setMatch_alnSet('a','h',def_AEqlH_alnDefs,alnSetST);
   setMatch_alnSet('a','v',def_AEqlV_alnDefs,alnSetST);
   setMatch_alnSet('a','n',def_AEqlN_alnDefs,alnSetST);
   setMatch_alnSet('a','x',def_AEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub02:
   *   - t as first base
   \*****************************************************/

   setMatch_alnSet('t','a',def_TEqlA_alnDefs,alnSetST);
   setMatch_alnSet('t','t',def_TEqlT_alnDefs,alnSetST);
   setMatch_alnSet('t','u',def_TEqlU_alnDefs,alnSetST);
   setMatch_alnSet('t','g',def_TEqlG_alnDefs,alnSetST);
   setMatch_alnSet('t','c',def_TEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('t','w',def_TEqlW_alnDefs,alnSetST);
   setMatch_alnSet('t','s',def_TEqlS_alnDefs,alnSetST);
   setMatch_alnSet('t','m',def_TEqlM_alnDefs,alnSetST);
   setMatch_alnSet('t','k',def_TEqlK_alnDefs,alnSetST);
   setMatch_alnSet('t','r',def_TEqlR_alnDefs,alnSetST);
   setMatch_alnSet('t','y',def_TEqlY_alnDefs,alnSetST);
   setMatch_alnSet('t','b',def_TEqlB_alnDefs,alnSetST);
   setMatch_alnSet('t','d',def_TEqlD_alnDefs,alnSetST);
   setMatch_alnSet('t','h',def_TEqlH_alnDefs,alnSetST);
   setMatch_alnSet('t','v',def_TEqlV_alnDefs,alnSetST);
   setMatch_alnSet('t','n',def_TEqlN_alnDefs,alnSetST);
   setMatch_alnSet('t','x',def_TEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub03:
   *   - u (t) as first base
   \*****************************************************/

   setMatch_alnSet('u','a',def_UEqlA_alnDefs,alnSetST);
   setMatch_alnSet('u','g',def_UEqlG_alnDefs,alnSetST);
   setMatch_alnSet('u','c',def_UEqlC_alnDefs,alnSetST);
   setMatch_alnSet('u','t',def_UEqlT_alnDefs,alnSetST);
   setMatch_alnSet('u','u',def_UEqlU_alnDefs,alnSetST);

   /*Set u & t to same scores (U is RNA version of T)*/
   setMatch_alnSet('u','w',def_UEqlW_alnDefs,alnSetST);
   setMatch_alnSet('u','s',def_UEqlS_alnDefs,alnSetST);
   setMatch_alnSet('u','m',def_UEqlM_alnDefs,alnSetST);
   setMatch_alnSet('u','k',def_UEqlK_alnDefs,alnSetST);
   setMatch_alnSet('u','r',def_UEqlR_alnDefs,alnSetST);
   setMatch_alnSet('u','y',def_UEqlY_alnDefs,alnSetST);
   setMatch_alnSet('u','b',def_UEqlB_alnDefs,alnSetST);
   setMatch_alnSet('u','d',def_UEqlD_alnDefs,alnSetST);
   setMatch_alnSet('u','h',def_UEqlH_alnDefs,alnSetST);
   setMatch_alnSet('u','v',def_UEqlV_alnDefs,alnSetST);
   setMatch_alnSet('u','n',def_UEqlN_alnDefs,alnSetST);
   setMatch_alnSet('u','x',def_UEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub04:
   *   - g as first base
   \*****************************************************/

   setMatch_alnSet('g','a',def_GEqlA_alnDefs,alnSetST);
   setMatch_alnSet('g','t',def_GEqlT_alnDefs,alnSetST);
   setMatch_alnSet('g','u',def_GEqlU_alnDefs,alnSetST);
   setMatch_alnSet('g','g',def_GEqlG_alnDefs,alnSetST);
   setMatch_alnSet('g','c',def_GEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('g','w',def_GEqlW_alnDefs,alnSetST);
   setMatch_alnSet('g','s',def_GEqlS_alnDefs,alnSetST);
   setMatch_alnSet('g','m',def_GEqlM_alnDefs,alnSetST);
   setMatch_alnSet('g','k',def_GEqlK_alnDefs,alnSetST);
   setMatch_alnSet('g','r',def_GEqlR_alnDefs,alnSetST);
   setMatch_alnSet('g','y',def_GEqlY_alnDefs,alnSetST);
   setMatch_alnSet('g','b',def_GEqlB_alnDefs,alnSetST);
   setMatch_alnSet('g','d',def_GEqlD_alnDefs,alnSetST);
   setMatch_alnSet('g','h',def_GEqlH_alnDefs,alnSetST);
   setMatch_alnSet('g','v',def_GEqlV_alnDefs,alnSetST);
   setMatch_alnSet('g','n',def_GEqlN_alnDefs,alnSetST);
   setMatch_alnSet('g','x',def_GEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub05:
   *   - c as first base
   \*****************************************************/

   setMatch_alnSet('c','a',def_CEqlA_alnDefs,alnSetST);
   setMatch_alnSet('c','t',def_CEqlT_alnDefs,alnSetST);
   setMatch_alnSet('c','u',def_CEqlT_alnDefs,alnSetST);
   setMatch_alnSet('c','g',def_CEqlG_alnDefs,alnSetST);
   setMatch_alnSet('c','c',def_CEqlC_alnDefs,alnSetST);

   /*anonymous matches*/
   setMatch_alnSet('c','w',def_CEqlW_alnDefs,alnSetST);
   setMatch_alnSet('c','s',def_CEqlS_alnDefs,alnSetST);
   setMatch_alnSet('c','m',def_CEqlM_alnDefs,alnSetST);
   setMatch_alnSet('c','k',def_CEqlK_alnDefs,alnSetST);
   setMatch_alnSet('c','r',def_CEqlR_alnDefs,alnSetST);
   setMatch_alnSet('c','y',def_CEqlY_alnDefs,alnSetST);
   setMatch_alnSet('c','b',def_CEqlB_alnDefs,alnSetST);
   setMatch_alnSet('c','d',def_CEqlD_alnDefs,alnSetST);
   setMatch_alnSet('c','h',def_CEqlH_alnDefs,alnSetST);
   setMatch_alnSet('c','v',def_CEqlV_alnDefs,alnSetST);
   setMatch_alnSet('c','n',def_CEqlN_alnDefs,alnSetST);
   setMatch_alnSet('c','x',def_CEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub06:
   *   - w (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('w','a',def_WEqlA_alnDefs,alnSetST);
   setMatch_alnSet('w','c',def_WEqlC_alnDefs,alnSetST);
   setMatch_alnSet('w','g',def_WEqlG_alnDefs,alnSetST);
   setMatch_alnSet('w','t',def_WEqlT_alnDefs,alnSetST);
   setMatch_alnSet('w','u',def_WEqlT_alnDefs,alnSetST);

   setMatch_alnSet('w','w',def_WEqlW_alnDefs,alnSetST);
   setMatch_alnSet('w','s',def_WEqlS_alnDefs,alnSetST);
   setMatch_alnSet('w','m',def_WEqlM_alnDefs,alnSetST);
   setMatch_alnSet('w','k',def_WEqlK_alnDefs,alnSetST);
   setMatch_alnSet('w','r',def_WEqlR_alnDefs,alnSetST);
   setMatch_alnSet('w','y',def_WEqlY_alnDefs,alnSetST);
   setMatch_alnSet('w','b',def_WEqlB_alnDefs,alnSetST);
   setMatch_alnSet('w','d',def_WEqlD_alnDefs,alnSetST);
   setMatch_alnSet('w','h',def_WEqlH_alnDefs,alnSetST);
   setMatch_alnSet('w','v',def_WEqlV_alnDefs,alnSetST);
   setMatch_alnSet('w','n',def_WEqlN_alnDefs,alnSetST);
   setMatch_alnSet('w','x',def_WEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub07:
   *   - s (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('s','a',def_SEqlA_alnDefs,alnSetST);
   setMatch_alnSet('s','c',def_SEqlC_alnDefs,alnSetST);
   setMatch_alnSet('s','g',def_SEqlG_alnDefs,alnSetST);
   setMatch_alnSet('s','t',def_SEqlT_alnDefs,alnSetST);
   setMatch_alnSet('s','u',def_SEqlT_alnDefs,alnSetST);

   setMatch_alnSet('s','w',def_SEqlW_alnDefs,alnSetST);
   setMatch_alnSet('s','s',def_SEqlS_alnDefs,alnSetST);
   setMatch_alnSet('s','m',def_SEqlM_alnDefs,alnSetST);
   setMatch_alnSet('s','k',def_SEqlK_alnDefs,alnSetST);
   setMatch_alnSet('s','r',def_SEqlR_alnDefs,alnSetST);
   setMatch_alnSet('s','y',def_SEqlY_alnDefs,alnSetST);
   setMatch_alnSet('s','b',def_SEqlB_alnDefs,alnSetST);
   setMatch_alnSet('s','d',def_SEqlD_alnDefs,alnSetST);
   setMatch_alnSet('s','h',def_SEqlH_alnDefs,alnSetST);
   setMatch_alnSet('s','v',def_SEqlV_alnDefs,alnSetST);
   setMatch_alnSet('s','n',def_SEqlN_alnDefs,alnSetST);
   setMatch_alnSet('s','x',def_SEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub08:
   *   - m (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('m','a',def_MEqlA_alnDefs,alnSetST);
   setMatch_alnSet('m','c',def_MEqlC_alnDefs,alnSetST);
   setMatch_alnSet('m','g',def_MEqlG_alnDefs,alnSetST);
   setMatch_alnSet('m','t',def_MEqlT_alnDefs,alnSetST);
   setMatch_alnSet('m','u',def_MEqlT_alnDefs,alnSetST);

   setMatch_alnSet('m','w',def_MEqlW_alnDefs,alnSetST);
   setMatch_alnSet('m','s',def_MEqlS_alnDefs,alnSetST);
   setMatch_alnSet('m','m',def_MEqlM_alnDefs,alnSetST);
   setMatch_alnSet('m','k',def_MEqlK_alnDefs,alnSetST);
   setMatch_alnSet('m','r',def_MEqlR_alnDefs,alnSetST);
   setMatch_alnSet('m','y',def_MEqlY_alnDefs,alnSetST);
   setMatch_alnSet('m','b',def_MEqlB_alnDefs,alnSetST);
   setMatch_alnSet('m','d',def_MEqlD_alnDefs,alnSetST);
   setMatch_alnSet('m','h',def_MEqlH_alnDefs,alnSetST);
   setMatch_alnSet('m','v',def_MEqlV_alnDefs,alnSetST);
   setMatch_alnSet('m','n',def_MEqlN_alnDefs,alnSetST);
   setMatch_alnSet('m','x',def_MEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub09:
   *   - k (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('k','a',def_KEqlA_alnDefs,alnSetST);
   setMatch_alnSet('k','c',def_KEqlC_alnDefs,alnSetST);
   setMatch_alnSet('k','g',def_KEqlG_alnDefs,alnSetST);
   setMatch_alnSet('k','t',def_KEqlT_alnDefs,alnSetST);
   setMatch_alnSet('k','u',def_KEqlT_alnDefs,alnSetST);

   setMatch_alnSet('k','w',def_KEqlW_alnDefs,alnSetST);
   setMatch_alnSet('k','s',def_KEqlS_alnDefs,alnSetST);
   setMatch_alnSet('k','m',def_KEqlM_alnDefs,alnSetST);
   setMatch_alnSet('k','k',def_KEqlK_alnDefs,alnSetST);
   setMatch_alnSet('k','r',def_KEqlR_alnDefs,alnSetST);
   setMatch_alnSet('k','y',def_KEqlY_alnDefs,alnSetST);
   setMatch_alnSet('k','b',def_KEqlB_alnDefs,alnSetST);
   setMatch_alnSet('k','d',def_KEqlD_alnDefs,alnSetST);
   setMatch_alnSet('k','h',def_KEqlH_alnDefs,alnSetST);
   setMatch_alnSet('k','v',def_KEqlV_alnDefs,alnSetST);
   setMatch_alnSet('k','n',def_KEqlN_alnDefs,alnSetST);
   setMatch_alnSet('k','x',def_KEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub10:
   *   - r (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('r','a',def_REqlA_alnDefs,alnSetST);
   setMatch_alnSet('r','c',def_REqlC_alnDefs,alnSetST);
   setMatch_alnSet('r','g',def_REqlG_alnDefs,alnSetST);
   setMatch_alnSet('r','t',def_REqlT_alnDefs,alnSetST);
   setMatch_alnSet('r','u',def_REqlT_alnDefs,alnSetST);

   setMatch_alnSet('r','w',def_REqlW_alnDefs,alnSetST);
   setMatch_alnSet('r','s',def_REqlS_alnDefs,alnSetST);
   setMatch_alnSet('r','m',def_REqlM_alnDefs,alnSetST);
   setMatch_alnSet('r','k',def_REqlK_alnDefs,alnSetST);
   setMatch_alnSet('r','r',def_REqlR_alnDefs,alnSetST);
   setMatch_alnSet('r','y',def_REqlY_alnDefs,alnSetST);
   setMatch_alnSet('r','b',def_REqlB_alnDefs,alnSetST);
   setMatch_alnSet('r','d',def_REqlD_alnDefs,alnSetST);
   setMatch_alnSet('r','h',def_REqlH_alnDefs,alnSetST);
   setMatch_alnSet('r','v',def_REqlV_alnDefs,alnSetST);
   setMatch_alnSet('r','n',def_REqlN_alnDefs,alnSetST);
   setMatch_alnSet('r','x',def_REqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub11:
   *   - y (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('y','a',def_YEqlA_alnDefs,alnSetST);
   setMatch_alnSet('y','c',def_YEqlC_alnDefs,alnSetST);
   setMatch_alnSet('y','g',def_YEqlG_alnDefs,alnSetST);
   setMatch_alnSet('y','t',def_YEqlT_alnDefs,alnSetST);
   setMatch_alnSet('y','u',def_YEqlT_alnDefs,alnSetST);

   setMatch_alnSet('y','w',def_YEqlW_alnDefs,alnSetST);
   setMatch_alnSet('y','s',def_YEqlS_alnDefs,alnSetST);
   setMatch_alnSet('y','m',def_YEqlM_alnDefs,alnSetST);
   setMatch_alnSet('y','k',def_YEqlK_alnDefs,alnSetST);
   setMatch_alnSet('y','r',def_YEqlR_alnDefs,alnSetST);
   setMatch_alnSet('y','y',def_YEqlY_alnDefs,alnSetST);
   setMatch_alnSet('y','b',def_YEqlB_alnDefs,alnSetST);
   setMatch_alnSet('y','d',def_YEqlD_alnDefs,alnSetST);
   setMatch_alnSet('y','h',def_YEqlH_alnDefs,alnSetST);
   setMatch_alnSet('y','v',def_YEqlV_alnDefs,alnSetST);
   setMatch_alnSet('y','n',def_YEqlN_alnDefs,alnSetST);
   setMatch_alnSet('y','x',def_YEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub12:
   *   - b (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('b','a',def_BEqlA_alnDefs,alnSetST);
   setMatch_alnSet('b','c',def_BEqlC_alnDefs,alnSetST);
   setMatch_alnSet('b','g',def_BEqlG_alnDefs,alnSetST);
   setMatch_alnSet('b','t',def_BEqlT_alnDefs,alnSetST);
   setMatch_alnSet('b','u',def_BEqlT_alnDefs,alnSetST);

   setMatch_alnSet('b','w',def_BEqlW_alnDefs,alnSetST);
   setMatch_alnSet('b','s',def_BEqlS_alnDefs,alnSetST);
   setMatch_alnSet('b','m',def_BEqlM_alnDefs,alnSetST);
   setMatch_alnSet('b','k',def_BEqlK_alnDefs,alnSetST);
   setMatch_alnSet('b','r',def_BEqlR_alnDefs,alnSetST);
   setMatch_alnSet('b','y',def_BEqlY_alnDefs,alnSetST);
   setMatch_alnSet('b','b',def_BEqlB_alnDefs,alnSetST);
   setMatch_alnSet('b','d',def_BEqlD_alnDefs,alnSetST);
   setMatch_alnSet('b','h',def_BEqlH_alnDefs,alnSetST);
   setMatch_alnSet('b','v',def_BEqlV_alnDefs,alnSetST);
   setMatch_alnSet('b','n',def_BEqlN_alnDefs,alnSetST);
   setMatch_alnSet('b','x',def_BEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub13:
   *   - d (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('d','a',def_DEqlA_alnDefs,alnSetST);
   setMatch_alnSet('d','c',def_DEqlC_alnDefs,alnSetST);
   setMatch_alnSet('d','g',def_DEqlG_alnDefs,alnSetST);
   setMatch_alnSet('d','t',def_DEqlT_alnDefs,alnSetST);
   setMatch_alnSet('d','u',def_DEqlT_alnDefs,alnSetST);

   setMatch_alnSet('d','w',def_DEqlW_alnDefs,alnSetST);
   setMatch_alnSet('d','s',def_DEqlS_alnDefs,alnSetST);
   setMatch_alnSet('d','m',def_DEqlM_alnDefs,alnSetST);
   setMatch_alnSet('d','k',def_DEqlK_alnDefs,alnSetST);
   setMatch_alnSet('d','r',def_DEqlR_alnDefs,alnSetST);
   setMatch_alnSet('d','y',def_DEqlY_alnDefs,alnSetST);
   setMatch_alnSet('d','b',def_DEqlB_alnDefs,alnSetST);
   setMatch_alnSet('d','d',def_DEqlD_alnDefs,alnSetST);
   setMatch_alnSet('d','h',def_DEqlH_alnDefs,alnSetST);
   setMatch_alnSet('d','v',def_DEqlV_alnDefs,alnSetST);
   setMatch_alnSet('d','n',def_DEqlN_alnDefs,alnSetST);
   setMatch_alnSet('d','x',def_DEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub14:
   *   - h (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('h','a',def_HEqlA_alnDefs,alnSetST);
   setMatch_alnSet('h','c',def_HEqlC_alnDefs,alnSetST);
   setMatch_alnSet('h','g',def_HEqlG_alnDefs,alnSetST);
   setMatch_alnSet('h','t',def_HEqlT_alnDefs,alnSetST);
   setMatch_alnSet('h','u',def_HEqlT_alnDefs,alnSetST);

   setMatch_alnSet('h','w',def_HEqlW_alnDefs,alnSetST);
   setMatch_alnSet('h','s',def_HEqlS_alnDefs,alnSetST);
   setMatch_alnSet('h','m',def_HEqlM_alnDefs,alnSetST);
   setMatch_alnSet('h','k',def_HEqlK_alnDefs,alnSetST);
   setMatch_alnSet('h','r',def_HEqlR_alnDefs,alnSetST);
   setMatch_alnSet('h','y',def_HEqlY_alnDefs,alnSetST);
   setMatch_alnSet('h','b',def_HEqlB_alnDefs,alnSetST);
   setMatch_alnSet('h','d',def_HEqlD_alnDefs,alnSetST);
   setMatch_alnSet('h','h',def_HEqlH_alnDefs,alnSetST);
   setMatch_alnSet('h','v',def_HEqlV_alnDefs,alnSetST);
   setMatch_alnSet('h','n',def_HEqlN_alnDefs,alnSetST);
   setMatch_alnSet('h','x',def_HEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub15:
   *   - v (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('v','a',def_VEqlA_alnDefs,alnSetST);
   setMatch_alnSet('v','c',def_VEqlC_alnDefs,alnSetST);
   setMatch_alnSet('v','g',def_VEqlG_alnDefs,alnSetST);
   setMatch_alnSet('v','t',def_VEqlT_alnDefs,alnSetST);
   setMatch_alnSet('v','u',def_VEqlT_alnDefs,alnSetST);

   setMatch_alnSet('v','w',def_VEqlW_alnDefs,alnSetST);
   setMatch_alnSet('v','s',def_VEqlS_alnDefs,alnSetST);
   setMatch_alnSet('v','m',def_VEqlM_alnDefs,alnSetST);
   setMatch_alnSet('v','k',def_VEqlK_alnDefs,alnSetST);
   setMatch_alnSet('v','r',def_VEqlR_alnDefs,alnSetST);
   setMatch_alnSet('v','y',def_VEqlY_alnDefs,alnSetST);
   setMatch_alnSet('v','b',def_VEqlB_alnDefs,alnSetST);
   setMatch_alnSet('v','d',def_VEqlD_alnDefs,alnSetST);
   setMatch_alnSet('v','h',def_VEqlH_alnDefs,alnSetST);
   setMatch_alnSet('v','v',def_VEqlV_alnDefs,alnSetST);
   setMatch_alnSet('v','n',def_VEqlN_alnDefs,alnSetST);
   setMatch_alnSet('v','x',def_VEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub16:
   *   - n (anonymous) as first base
   \*****************************************************/

   setMatch_alnSet('n','a',def_NEqlA_alnDefs,alnSetST);
   setMatch_alnSet('n','c',def_NEqlC_alnDefs,alnSetST);
   setMatch_alnSet('n','g',def_NEqlG_alnDefs,alnSetST);
   setMatch_alnSet('n','t',def_NEqlT_alnDefs,alnSetST);
   setMatch_alnSet('n','u',def_NEqlT_alnDefs,alnSetST);

   setMatch_alnSet('n','w',def_NEqlW_alnDefs,alnSetST);
   setMatch_alnSet('n','s',def_NEqlS_alnDefs,alnSetST);
   setMatch_alnSet('n','m',def_NEqlM_alnDefs,alnSetST);
   setMatch_alnSet('n','k',def_NEqlK_alnDefs,alnSetST);
   setMatch_alnSet('n','r',def_NEqlR_alnDefs,alnSetST);
   setMatch_alnSet('n','y',def_NEqlY_alnDefs,alnSetST);
   setMatch_alnSet('n','b',def_NEqlB_alnDefs,alnSetST);
   setMatch_alnSet('n','d',def_NEqlD_alnDefs,alnSetST);
   setMatch_alnSet('n','h',def_NEqlH_alnDefs,alnSetST);
   setMatch_alnSet('n','v',def_NEqlV_alnDefs,alnSetST);
   setMatch_alnSet('n','n',def_NEqlN_alnDefs,alnSetST);
   setMatch_alnSet('n','x',def_NEqlX_alnDefs,alnSetST);

   /*****************************************************\
   * Fun13 Sec05 Sub17:
   *   - x (anonymous) as first base (technically aa)
   \*****************************************************/

   setMatch_alnSet('x','a',def_XEqlA_alnDefs,alnSetST);
   setMatch_alnSet('x','c',def_XEqlC_alnDefs,alnSetST);
   setMatch_alnSet('x','g',def_XEqlG_alnDefs,alnSetST);
   setMatch_alnSet('x','t',def_XEqlT_alnDefs,alnSetST);
   setMatch_alnSet('x','u',def_XEqlT_alnDefs,alnSetST);

   setMatch_alnSet('x','w',def_XEqlW_alnDefs,alnSetST);
   setMatch_alnSet('x','s',def_XEqlS_alnDefs,alnSetST);
   setMatch_alnSet('x','m',def_XEqlM_alnDefs,alnSetST);
   setMatch_alnSet('x','k',def_XEqlK_alnDefs,alnSetST);
   setMatch_alnSet('x','r',def_XEqlR_alnDefs,alnSetST);
   setMatch_alnSet('x','y',def_XEqlY_alnDefs,alnSetST);
   setMatch_alnSet('x','b',def_XEqlB_alnDefs,alnSetST);
   setMatch_alnSet('x','d',def_XEqlD_alnDefs,alnSetST);
   setMatch_alnSet('x','h',def_XEqlH_alnDefs,alnSetST);
   setMatch_alnSet('x','v',def_XEqlV_alnDefs,alnSetST);
   setMatch_alnSet('x','n',def_XEqlN_alnDefs,alnSetST);
   setMatch_alnSet('x','x',def_XEqlX_alnDefs,alnSetST);

   return;
} /*init_alnSet*/

/*-------------------------------------------------------\
| Fun14: pDefMatchMatrix_alnSet
|   - print out the default match matrix
| Input:
|   - outFILE:
|     o file to print match matrix to
| Output:
|   - Prints:
|     o default match matrix in alnDefs.h to outFILE
\-------------------------------------------------------*/
void
pDefMatchMatrix_alnSet(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun14 TOC:
   '   - print out the default match matrix
   '   o fun14 sec01:
   '     - print out header
   '   o fun14 sec02:
   '     - print out match matrix
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec01:
   ^  - print out header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   fprintf(
      (FILE *) outFILE,
      "// default match matrix\n"
   );

   fprintf(
      (FILE *) outFILE,
      "// legend:\n"
   );

   fprintf(
      (FILE *) outFILE,
      "//   - %i is match\n",
      def_ntEql_alnDefs
   );

   fprintf(
      (FILE *) outFILE,
      "//   - %i is mismatch\n",
      def_ntNotEql_alnDefs
   );

   fprintf(
      (FILE *) outFILE,
      "//   - %i is anonymous match\n",
      def_anonMatch_alnDefs
   );

   fprintf(
      (FILE *) outFILE,
      "//   - %i flag for anonymous (alone is mismatch)\n",
      def_anonymous_alnDefs
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec02:
   ^  - print out match matrix
   ^  o fun14 sec02 sub01:
   ^    - a as first base
   ^  o fun14 sec02 sub02:
   ^    - t as first base
   ^  o fun14 sec02 sub03:
   ^    - u (t) as first base
   ^  o fun14 sec02 sub04:
   ^    - g as first base
   ^  o fun14 sec02 sub05:
   ^    - c as first base
   ^  o fun14 sec02 sub06:
   ^    - w (anonymous) as first base
   ^  o fun14 sec02 sub07:
   ^    - s (anonymous) as first base
   ^  o fun14 sec02 sub08:
   ^    - m (anonymous) as first base
   ^  o fun14 sec02 sub09:
   ^    - k (anonymous) as first base
   ^  o fun14 sec02 sub10:
   ^    - r (anonymous) as first base
   ^  o fun14 sec02 sub11:
   ^    - y (anonymous) as first base
   ^  o fun14 sec02 sub12:
   ^    - b (anonymous) as first base
   ^  o fun14 sec02 sub13:
   ^    - d (anonymous) as first base
   ^  o fun14 sec02 sub14:
   ^    - h (anonymous) as first base
   ^  o fun14 sec02 sub15:
   ^    - v (anonymous) as first base
   ^  o fun14 sec02 sub16:
   ^    - n (anonymous) as first base
   ^  o fun14 sec02 sub17:
   ^    - x (anonymous) as first base (technically aa)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun14 Sec02 Sub01:
   *   - a as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"a a %i\n",def_AEqlA_alnDefs);
   fprintf((FILE *) outFILE,"a t %i\n",def_AEqlT_alnDefs);
   fprintf((FILE *) outFILE,"a u %i\n",def_AEqlU_alnDefs);
   fprintf((FILE *) outFILE,"a g %i\n",def_AEqlG_alnDefs);
   fprintf((FILE *) outFILE,"a c %i\n",def_AEqlC_alnDefs);

   /*anonymous matches*/
   fprintf((FILE *) outFILE,"a w %i\n",def_AEqlW_alnDefs);
   fprintf((FILE *) outFILE,"a s %i\n",def_AEqlS_alnDefs);
   fprintf((FILE *) outFILE,"a m %i\n",def_AEqlM_alnDefs);
   fprintf((FILE *) outFILE,"a k %i\n",def_AEqlK_alnDefs);
   fprintf((FILE *) outFILE,"a r %i\n",def_AEqlR_alnDefs);
   fprintf((FILE *) outFILE,"a y %i\n",def_AEqlY_alnDefs);
   fprintf((FILE *) outFILE,"a b %i\n",def_AEqlB_alnDefs);
   fprintf((FILE *) outFILE,"a d %i\n",def_AEqlD_alnDefs);
   fprintf((FILE *) outFILE,"a h %i\n",def_AEqlH_alnDefs);
   fprintf((FILE *) outFILE,"a v %i\n",def_AEqlV_alnDefs);
   fprintf((FILE *) outFILE,"a n %i\n",def_AEqlN_alnDefs);
   fprintf((FILE *) outFILE,"a x %i\n",def_AEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub02:
   *   - t as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"t a %i\n",def_TEqlA_alnDefs);
   fprintf((FILE *) outFILE,"t t %i\n",def_TEqlT_alnDefs);
   fprintf((FILE *) outFILE,"t u %i\n",def_TEqlU_alnDefs);
   fprintf((FILE *) outFILE,"t g %i\n",def_TEqlG_alnDefs);
   fprintf((FILE *) outFILE,"t c %i\n",def_TEqlC_alnDefs);

   /*anonymous matches*/
   fprintf((FILE *) outFILE,"t w %i\n",def_TEqlW_alnDefs);
   fprintf((FILE *) outFILE,"t s %i\n",def_TEqlS_alnDefs);
   fprintf((FILE *) outFILE,"t m %i\n",def_TEqlM_alnDefs);
   fprintf((FILE *) outFILE,"t k %i\n",def_TEqlK_alnDefs);
   fprintf((FILE *) outFILE,"t r %i\n",def_TEqlR_alnDefs);
   fprintf((FILE *) outFILE,"t y %i\n",def_TEqlY_alnDefs);
   fprintf((FILE *) outFILE,"t b %i\n",def_TEqlB_alnDefs);
   fprintf((FILE *) outFILE,"t d %i\n",def_TEqlD_alnDefs);
   fprintf((FILE *) outFILE,"t h %i\n",def_TEqlH_alnDefs);
   fprintf((FILE *) outFILE,"t v %i\n",def_TEqlV_alnDefs);
   fprintf((FILE *) outFILE,"t n %i\n",def_TEqlN_alnDefs);
   fprintf((FILE *) outFILE,"t x %i\n",def_TEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub03:
   *   - u (t) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"u a %i\n",def_UEqlA_alnDefs);
   fprintf((FILE *) outFILE,"u g %i\n",def_UEqlG_alnDefs);
   fprintf((FILE *) outFILE,"u c %i\n",def_UEqlC_alnDefs);
   fprintf((FILE *) outFILE,"u t %i\n",def_UEqlT_alnDefs);
   fprintf((FILE *) outFILE,"u u %i\n",def_UEqlU_alnDefs);

   /*Set u & t to same scores (U is RNA version of T)*/
   fprintf((FILE *) outFILE,"u w %i\n",def_UEqlW_alnDefs);
   fprintf((FILE *) outFILE,"u s %i\n",def_UEqlS_alnDefs);
   fprintf((FILE *) outFILE,"u m %i\n",def_UEqlM_alnDefs);
   fprintf((FILE *) outFILE,"u k %i\n",def_UEqlK_alnDefs);
   fprintf((FILE *) outFILE,"u r %i\n",def_UEqlR_alnDefs);
   fprintf((FILE *) outFILE,"u y %i\n",def_UEqlY_alnDefs);
   fprintf((FILE *) outFILE,"u b %i\n",def_UEqlB_alnDefs);
   fprintf((FILE *) outFILE,"u d %i\n",def_UEqlD_alnDefs);
   fprintf((FILE *) outFILE,"u h %i\n",def_UEqlH_alnDefs);
   fprintf((FILE *) outFILE,"u v %i\n",def_UEqlV_alnDefs);
   fprintf((FILE *) outFILE,"u n %i\n",def_UEqlN_alnDefs);
   fprintf((FILE *) outFILE,"u x %i\n",def_UEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub04:
   *   - g as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"g a %i\n",def_GEqlA_alnDefs);
   fprintf((FILE *) outFILE,"g t %i\n",def_GEqlT_alnDefs);
   fprintf((FILE *) outFILE,"g u %i\n",def_GEqlU_alnDefs);
   fprintf((FILE *) outFILE,"g g %i\n",def_GEqlG_alnDefs);
   fprintf((FILE *) outFILE,"g c %i\n",def_GEqlC_alnDefs);

   /*anonymous matches*/
   fprintf((FILE *) outFILE,"g w %i\n",def_GEqlW_alnDefs);
   fprintf((FILE *) outFILE,"g s %i\n",def_GEqlS_alnDefs);
   fprintf((FILE *) outFILE,"g m %i\n",def_GEqlM_alnDefs);
   fprintf((FILE *) outFILE,"g k %i\n",def_GEqlK_alnDefs);
   fprintf((FILE *) outFILE,"g r %i\n",def_GEqlR_alnDefs);
   fprintf((FILE *) outFILE,"g y %i\n",def_GEqlY_alnDefs);
   fprintf((FILE *) outFILE,"g b %i\n",def_GEqlB_alnDefs);
   fprintf((FILE *) outFILE,"g d %i\n",def_GEqlD_alnDefs);
   fprintf((FILE *) outFILE,"g h %i\n",def_GEqlH_alnDefs);
   fprintf((FILE *) outFILE,"g v %i\n",def_GEqlV_alnDefs);
   fprintf((FILE *) outFILE,"g n %i\n",def_GEqlN_alnDefs);
   fprintf((FILE *) outFILE,"g x %i\n",def_GEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub05:
   *   - c as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"c a %i\n",def_CEqlA_alnDefs);
   fprintf((FILE *) outFILE,"c t %i\n",def_CEqlT_alnDefs);
   fprintf((FILE *) outFILE,"c u %i\n",def_CEqlT_alnDefs);
   fprintf((FILE *) outFILE,"c g %i\n",def_CEqlG_alnDefs);
   fprintf((FILE *) outFILE,"c c %i\n",def_CEqlC_alnDefs);

   /*anonymous matches*/
   fprintf((FILE *) outFILE,"c w %i\n",def_CEqlW_alnDefs);
   fprintf((FILE *) outFILE,"c s %i\n",def_CEqlS_alnDefs);
   fprintf((FILE *) outFILE,"c m %i\n",def_CEqlM_alnDefs);
   fprintf((FILE *) outFILE,"c k %i\n",def_CEqlK_alnDefs);
   fprintf((FILE *) outFILE,"c r %i\n",def_CEqlR_alnDefs);
   fprintf((FILE *) outFILE,"c y %i\n",def_CEqlY_alnDefs);
   fprintf((FILE *) outFILE,"c b %i\n",def_CEqlB_alnDefs);
   fprintf((FILE *) outFILE,"c d %i\n",def_CEqlD_alnDefs);
   fprintf((FILE *) outFILE,"c h %i\n",def_CEqlH_alnDefs);
   fprintf((FILE *) outFILE,"c v %i\n",def_CEqlV_alnDefs);
   fprintf((FILE *) outFILE,"c n %i\n",def_CEqlN_alnDefs);
   fprintf((FILE *) outFILE,"c x %i\n",def_CEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub06:
   *   - w (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"w a %i\n",def_WEqlA_alnDefs);
   fprintf((FILE *) outFILE,"w c %i\n",def_WEqlC_alnDefs);
   fprintf((FILE *) outFILE,"w g %i\n",def_WEqlG_alnDefs);
   fprintf((FILE *) outFILE,"w t %i\n",def_WEqlT_alnDefs);
   fprintf((FILE *) outFILE,"w u %i\n",def_WEqlT_alnDefs);

   fprintf((FILE *) outFILE,"w w %i\n",def_WEqlW_alnDefs);
   fprintf((FILE *) outFILE,"w s %i\n",def_WEqlS_alnDefs);
   fprintf((FILE *) outFILE,"w m %i\n",def_WEqlM_alnDefs);
   fprintf((FILE *) outFILE,"w k %i\n",def_WEqlK_alnDefs);
   fprintf((FILE *) outFILE,"w r %i\n",def_WEqlR_alnDefs);
   fprintf((FILE *) outFILE,"w y %i\n",def_WEqlY_alnDefs);
   fprintf((FILE *) outFILE,"w b %i\n",def_WEqlB_alnDefs);
   fprintf((FILE *) outFILE,"w d %i\n",def_WEqlD_alnDefs);
   fprintf((FILE *) outFILE,"w h %i\n",def_WEqlH_alnDefs);
   fprintf((FILE *) outFILE,"w v %i\n",def_WEqlV_alnDefs);
   fprintf((FILE *) outFILE,"w n %i\n",def_WEqlN_alnDefs);
   fprintf((FILE *) outFILE,"w x %i\n",def_WEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub07:
   *   - s (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"s a %i\n",def_SEqlA_alnDefs);
   fprintf((FILE *) outFILE,"s c %i\n",def_SEqlC_alnDefs);
   fprintf((FILE *) outFILE,"s g %i\n",def_SEqlG_alnDefs);
   fprintf((FILE *) outFILE,"s t %i\n",def_SEqlT_alnDefs);
   fprintf((FILE *) outFILE,"s u %i\n",def_SEqlT_alnDefs);

   fprintf((FILE *) outFILE,"s w %i\n",def_SEqlW_alnDefs);
   fprintf((FILE *) outFILE,"s s %i\n",def_SEqlS_alnDefs);
   fprintf((FILE *) outFILE,"s m %i\n",def_SEqlM_alnDefs);
   fprintf((FILE *) outFILE,"s k %i\n",def_SEqlK_alnDefs);
   fprintf((FILE *) outFILE,"s r %i\n",def_SEqlR_alnDefs);
   fprintf((FILE *) outFILE,"s y %i\n",def_SEqlY_alnDefs);
   fprintf((FILE *) outFILE,"s b %i\n",def_SEqlB_alnDefs);
   fprintf((FILE *) outFILE,"s d %i\n",def_SEqlD_alnDefs);
   fprintf((FILE *) outFILE,"s h %i\n",def_SEqlH_alnDefs);
   fprintf((FILE *) outFILE,"s v %i\n",def_SEqlV_alnDefs);
   fprintf((FILE *) outFILE,"s n %i\n",def_SEqlN_alnDefs);
   fprintf((FILE *) outFILE,"s x %i\n",def_SEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub08:
   *   - m (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"m a %i\n",def_MEqlA_alnDefs);
   fprintf((FILE *) outFILE,"m c %i\n",def_MEqlC_alnDefs);
   fprintf((FILE *) outFILE,"m g %i\n",def_MEqlG_alnDefs);
   fprintf((FILE *) outFILE,"m t %i\n",def_MEqlT_alnDefs);
   fprintf((FILE *) outFILE,"m u %i\n",def_MEqlT_alnDefs);

   fprintf((FILE *) outFILE,"m w %i\n",def_MEqlW_alnDefs);
   fprintf((FILE *) outFILE,"m s %i\n",def_MEqlS_alnDefs);
   fprintf((FILE *) outFILE,"m m %i\n",def_MEqlM_alnDefs);
   fprintf((FILE *) outFILE,"m k %i\n",def_MEqlK_alnDefs);
   fprintf((FILE *) outFILE,"m r %i\n",def_MEqlR_alnDefs);
   fprintf((FILE *) outFILE,"m y %i\n",def_MEqlY_alnDefs);
   fprintf((FILE *) outFILE,"m b %i\n",def_MEqlB_alnDefs);
   fprintf((FILE *) outFILE,"m d %i\n",def_MEqlD_alnDefs);
   fprintf((FILE *) outFILE,"m h %i\n",def_MEqlH_alnDefs);
   fprintf((FILE *) outFILE,"m v %i\n",def_MEqlV_alnDefs);
   fprintf((FILE *) outFILE,"m n %i\n",def_MEqlN_alnDefs);
   fprintf((FILE *) outFILE,"m x %i\n",def_MEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub09:
   *   - k (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"k a %i\n",def_KEqlA_alnDefs);
   fprintf((FILE *) outFILE,"k c %i\n",def_KEqlC_alnDefs);
   fprintf((FILE *) outFILE,"k g %i\n",def_KEqlG_alnDefs);
   fprintf((FILE *) outFILE,"k t %i\n",def_KEqlT_alnDefs);
   fprintf((FILE *) outFILE,"k u %i\n",def_KEqlT_alnDefs);

   fprintf((FILE *) outFILE,"k w %i\n",def_KEqlW_alnDefs);
   fprintf((FILE *) outFILE,"k s %i\n",def_KEqlS_alnDefs);
   fprintf((FILE *) outFILE,"k m %i\n",def_KEqlM_alnDefs);
   fprintf((FILE *) outFILE,"k k %i\n",def_KEqlK_alnDefs);
   fprintf((FILE *) outFILE,"k r %i\n",def_KEqlR_alnDefs);
   fprintf((FILE *) outFILE,"k y %i\n",def_KEqlY_alnDefs);
   fprintf((FILE *) outFILE,"k b %i\n",def_KEqlB_alnDefs);
   fprintf((FILE *) outFILE,"k d %i\n",def_KEqlD_alnDefs);
   fprintf((FILE *) outFILE,"k h %i\n",def_KEqlH_alnDefs);
   fprintf((FILE *) outFILE,"k v %i\n",def_KEqlV_alnDefs);
   fprintf((FILE *) outFILE,"k n %i\n",def_KEqlN_alnDefs);
   fprintf((FILE *) outFILE,"k x %i\n",def_KEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub10:
   *   - r (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"r a %i\n",def_REqlA_alnDefs);
   fprintf((FILE *) outFILE,"r c %i\n",def_REqlC_alnDefs);
   fprintf((FILE *) outFILE,"r g %i\n",def_REqlG_alnDefs);
   fprintf((FILE *) outFILE,"r t %i\n",def_REqlT_alnDefs);
   fprintf((FILE *) outFILE,"r u %i\n",def_REqlT_alnDefs);

   fprintf((FILE *) outFILE,"r w %i\n",def_REqlW_alnDefs);
   fprintf((FILE *) outFILE,"r s %i\n",def_REqlS_alnDefs);
   fprintf((FILE *) outFILE,"r m %i\n",def_REqlM_alnDefs);
   fprintf((FILE *) outFILE,"r k %i\n",def_REqlK_alnDefs);
   fprintf((FILE *) outFILE,"r r %i\n",def_REqlR_alnDefs);
   fprintf((FILE *) outFILE,"r y %i\n",def_REqlY_alnDefs);
   fprintf((FILE *) outFILE,"r b %i\n",def_REqlB_alnDefs);
   fprintf((FILE *) outFILE,"r d %i\n",def_REqlD_alnDefs);
   fprintf((FILE *) outFILE,"r h %i\n",def_REqlH_alnDefs);
   fprintf((FILE *) outFILE,"r v %i\n",def_REqlV_alnDefs);
   fprintf((FILE *) outFILE,"r n %i\n",def_REqlN_alnDefs);
   fprintf((FILE *) outFILE,"r x %i\n",def_REqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub11:
   *   - y (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"y a %i\n",def_YEqlA_alnDefs);
   fprintf((FILE *) outFILE,"y c %i\n",def_YEqlC_alnDefs);
   fprintf((FILE *) outFILE,"y g %i\n",def_YEqlG_alnDefs);
   fprintf((FILE *) outFILE,"y t %i\n",def_YEqlT_alnDefs);
   fprintf((FILE *) outFILE,"y u %i\n",def_YEqlT_alnDefs);

   fprintf((FILE *) outFILE,"y w %i\n",def_YEqlW_alnDefs);
   fprintf((FILE *) outFILE,"y s %i\n",def_YEqlS_alnDefs);
   fprintf((FILE *) outFILE,"y m %i\n",def_YEqlM_alnDefs);
   fprintf((FILE *) outFILE,"y k %i\n",def_YEqlK_alnDefs);
   fprintf((FILE *) outFILE,"y r %i\n",def_YEqlR_alnDefs);
   fprintf((FILE *) outFILE,"y y %i\n",def_YEqlY_alnDefs);
   fprintf((FILE *) outFILE,"y b %i\n",def_YEqlB_alnDefs);
   fprintf((FILE *) outFILE,"y d %i\n",def_YEqlD_alnDefs);
   fprintf((FILE *) outFILE,"y h %i\n",def_YEqlH_alnDefs);
   fprintf((FILE *) outFILE,"y v %i\n",def_YEqlV_alnDefs);
   fprintf((FILE *) outFILE,"y n %i\n",def_YEqlN_alnDefs);
   fprintf((FILE *) outFILE,"y x %i\n",def_YEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub12:
   *   - b (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"b a %i\n",def_BEqlA_alnDefs);
   fprintf((FILE *) outFILE,"b c %i\n",def_BEqlC_alnDefs);
   fprintf((FILE *) outFILE,"b g %i\n",def_BEqlG_alnDefs);
   fprintf((FILE *) outFILE,"b t %i\n",def_BEqlT_alnDefs);
   fprintf((FILE *) outFILE,"b u %i\n",def_BEqlT_alnDefs);

   fprintf((FILE *) outFILE,"b w %i\n",def_BEqlW_alnDefs);
   fprintf((FILE *) outFILE,"b s %i\n",def_BEqlS_alnDefs);
   fprintf((FILE *) outFILE,"b m %i\n",def_BEqlM_alnDefs);
   fprintf((FILE *) outFILE,"b k %i\n",def_BEqlK_alnDefs);
   fprintf((FILE *) outFILE,"b r %i\n",def_BEqlR_alnDefs);
   fprintf((FILE *) outFILE,"b y %i\n",def_BEqlY_alnDefs);
   fprintf((FILE *) outFILE,"b b %i\n",def_BEqlB_alnDefs);
   fprintf((FILE *) outFILE,"b d %i\n",def_BEqlD_alnDefs);
   fprintf((FILE *) outFILE,"b h %i\n",def_BEqlH_alnDefs);
   fprintf((FILE *) outFILE,"b v %i\n",def_BEqlV_alnDefs);
   fprintf((FILE *) outFILE,"b n %i\n",def_BEqlN_alnDefs);
   fprintf((FILE *) outFILE,"b x %i\n",def_BEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub13:
   *   - d (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"d a %i\n",def_DEqlA_alnDefs);
   fprintf((FILE *) outFILE,"d c %i\n",def_DEqlC_alnDefs);
   fprintf((FILE *) outFILE,"d g %i\n",def_DEqlG_alnDefs);
   fprintf((FILE *) outFILE,"d t %i\n",def_DEqlT_alnDefs);
   fprintf((FILE *) outFILE,"d u %i\n",def_DEqlT_alnDefs);

   fprintf((FILE *) outFILE,"d w %i\n",def_DEqlW_alnDefs);
   fprintf((FILE *) outFILE,"d s %i\n",def_DEqlS_alnDefs);
   fprintf((FILE *) outFILE,"d m %i\n",def_DEqlM_alnDefs);
   fprintf((FILE *) outFILE,"d k %i\n",def_DEqlK_alnDefs);
   fprintf((FILE *) outFILE,"d r %i\n",def_DEqlR_alnDefs);
   fprintf((FILE *) outFILE,"d y %i\n",def_DEqlY_alnDefs);
   fprintf((FILE *) outFILE,"d b %i\n",def_DEqlB_alnDefs);
   fprintf((FILE *) outFILE,"d d %i\n",def_DEqlD_alnDefs);
   fprintf((FILE *) outFILE,"d h %i\n",def_DEqlH_alnDefs);
   fprintf((FILE *) outFILE,"d v %i\n",def_DEqlV_alnDefs);
   fprintf((FILE *) outFILE,"d n %i\n",def_DEqlN_alnDefs);
   fprintf((FILE *) outFILE,"d x %i\n",def_DEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub14:
   *   - h (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"h a %i\n",def_HEqlA_alnDefs);
   fprintf((FILE *) outFILE,"h c %i\n",def_HEqlC_alnDefs);
   fprintf((FILE *) outFILE,"h g %i\n",def_HEqlG_alnDefs);
   fprintf((FILE *) outFILE,"h t %i\n",def_HEqlT_alnDefs);
   fprintf((FILE *) outFILE,"h u %i\n",def_HEqlT_alnDefs);

   fprintf((FILE *) outFILE,"h w %i\n",def_HEqlW_alnDefs);
   fprintf((FILE *) outFILE,"h s %i\n",def_HEqlS_alnDefs);
   fprintf((FILE *) outFILE,"h m %i\n",def_HEqlM_alnDefs);
   fprintf((FILE *) outFILE,"h k %i\n",def_HEqlK_alnDefs);
   fprintf((FILE *) outFILE,"h r %i\n",def_HEqlR_alnDefs);
   fprintf((FILE *) outFILE,"h y %i\n",def_HEqlY_alnDefs);
   fprintf((FILE *) outFILE,"h b %i\n",def_HEqlB_alnDefs);
   fprintf((FILE *) outFILE,"h d %i\n",def_HEqlD_alnDefs);
   fprintf((FILE *) outFILE,"h h %i\n",def_HEqlH_alnDefs);
   fprintf((FILE *) outFILE,"h v %i\n",def_HEqlV_alnDefs);
   fprintf((FILE *) outFILE,"h n %i\n",def_HEqlN_alnDefs);
   fprintf((FILE *) outFILE,"h x %i\n",def_HEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub15:
   *   - v (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"v a %i\n",def_VEqlA_alnDefs);
   fprintf((FILE *) outFILE,"v c %i\n",def_VEqlC_alnDefs);
   fprintf((FILE *) outFILE,"v g %i\n",def_VEqlG_alnDefs);
   fprintf((FILE *) outFILE,"v t %i\n",def_VEqlT_alnDefs);
   fprintf((FILE *) outFILE,"v u %i\n",def_VEqlT_alnDefs);

   fprintf((FILE *) outFILE,"v w %i\n",def_VEqlW_alnDefs);
   fprintf((FILE *) outFILE,"v s %i\n",def_VEqlS_alnDefs);
   fprintf((FILE *) outFILE,"v m %i\n",def_VEqlM_alnDefs);
   fprintf((FILE *) outFILE,"v k %i\n",def_VEqlK_alnDefs);
   fprintf((FILE *) outFILE,"v r %i\n",def_VEqlR_alnDefs);
   fprintf((FILE *) outFILE,"v y %i\n",def_VEqlY_alnDefs);
   fprintf((FILE *) outFILE,"v b %i\n",def_VEqlB_alnDefs);
   fprintf((FILE *) outFILE,"v d %i\n",def_VEqlD_alnDefs);
   fprintf((FILE *) outFILE,"v h %i\n",def_VEqlH_alnDefs);
   fprintf((FILE *) outFILE,"v v %i\n",def_VEqlV_alnDefs);
   fprintf((FILE *) outFILE,"v n %i\n",def_VEqlN_alnDefs);
   fprintf((FILE *) outFILE,"v x %i\n",def_VEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub16:
   *   - n (anonymous) as first base
   \*****************************************************/

   fprintf((FILE *) outFILE,"n a %i\n",def_NEqlA_alnDefs);
   fprintf((FILE *) outFILE,"n c %i\n",def_NEqlC_alnDefs);
   fprintf((FILE *) outFILE,"n g %i\n",def_NEqlG_alnDefs);
   fprintf((FILE *) outFILE,"n t %i\n",def_NEqlT_alnDefs);
   fprintf((FILE *) outFILE,"n u %i\n",def_NEqlT_alnDefs);

   fprintf((FILE *) outFILE,"n w %i\n",def_NEqlW_alnDefs);
   fprintf((FILE *) outFILE,"n s %i\n",def_NEqlS_alnDefs);
   fprintf((FILE *) outFILE,"n m %i\n",def_NEqlM_alnDefs);
   fprintf((FILE *) outFILE,"n k %i\n",def_NEqlK_alnDefs);
   fprintf((FILE *) outFILE,"n r %i\n",def_NEqlR_alnDefs);
   fprintf((FILE *) outFILE,"n y %i\n",def_NEqlY_alnDefs);
   fprintf((FILE *) outFILE,"n b %i\n",def_NEqlB_alnDefs);
   fprintf((FILE *) outFILE,"n d %i\n",def_NEqlD_alnDefs);
   fprintf((FILE *) outFILE,"n h %i\n",def_NEqlH_alnDefs);
   fprintf((FILE *) outFILE,"n v %i\n",def_NEqlV_alnDefs);
   fprintf((FILE *) outFILE,"n n %i\n",def_NEqlN_alnDefs);
   fprintf((FILE *) outFILE,"n x %i\n",def_NEqlX_alnDefs);

   /*****************************************************\
   * Fun14 Sec02 Sub17:
   *   - x (anonymous) as first base (technically aa)
   \*****************************************************/

   fprintf((FILE *) outFILE,"x a %i\n",def_XEqlA_alnDefs);
   fprintf((FILE *) outFILE,"x c %i\n",def_XEqlC_alnDefs);
   fprintf((FILE *) outFILE,"x g %i\n",def_XEqlG_alnDefs);
   fprintf((FILE *) outFILE,"x t %i\n",def_XEqlT_alnDefs);
   fprintf((FILE *) outFILE,"x u %i\n",def_XEqlT_alnDefs);

   fprintf((FILE *) outFILE,"x w %i\n",def_XEqlW_alnDefs);
   fprintf((FILE *) outFILE,"x s %i\n",def_XEqlS_alnDefs);
   fprintf((FILE *) outFILE,"x m %i\n",def_XEqlM_alnDefs);
   fprintf((FILE *) outFILE,"x k %i\n",def_XEqlK_alnDefs);
   fprintf((FILE *) outFILE,"x r %i\n",def_XEqlR_alnDefs);
   fprintf((FILE *) outFILE,"x y %i\n",def_XEqlY_alnDefs);
   fprintf((FILE *) outFILE,"x b %i\n",def_XEqlB_alnDefs);
   fprintf((FILE *) outFILE,"x d %i\n",def_XEqlD_alnDefs);
   fprintf((FILE *) outFILE,"x h %i\n",def_XEqlH_alnDefs);
   fprintf((FILE *) outFILE,"x v %i\n",def_XEqlV_alnDefs);
   fprintf((FILE *) outFILE,"x n %i\n",def_XEqlN_alnDefs);
   fprintf((FILE *) outFILE,"x x %i\n",def_XEqlX_alnDefs);
} /*pDefMatchMatrix_alnSet*/

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
