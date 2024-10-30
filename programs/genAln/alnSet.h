/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' alnSet SOF: Start Of File
'  - holds settings structures and supporting functions
'    for my pairwise aligners
'  o header:
'    - Header guards and defined variables
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
'  o fun13 init_alnSet:
'    - Set all values in altSet (alingment settings)
'      structure to defaults
'  o fun12: maxScore_alnSet
'    - finds maximum score possible for a sequence
'  o fun14: pDefMatchMatrix_alnSet
'    - print out the default match matrix
'  o license:
'    - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - header guards and defined variables
\-------------------------------------------------------*/

#ifndef ALIGNMENT_SETTINGS_STRUCT_H
#define ALIGNMENT_SETTINGS_STRUCT_H

#define def_lenCol_alnSet 27 /*one column/row of matrix*/

/*For seq conversion to index functions*/
#define def_ntToCode_alnSet (1 | 2 | 4 | 8 | 16)
  /*clear 64 bit and case*/

/*-------------------------------------------------------\
| ST01: alnSet
| Use: Holds settings for my alignment program
\-------------------------------------------------------*/
typedef struct alnSet
{ /*alnSet*/
   /*Line wrap for printing out an alignment*/
   unsigned int lineWrapUI;

   signed char pBasePosBl; /*1 Print out base numbers*/
   signed char pFullAlnBl;
     /*1: Print out the full alignmnet
     ` 0: Print out the aligned region
     */

   signed char formatFlag;
     /*defExpandCig: is default format (S D I = X)
     ` defEMBOSS: is EMBOSS format (| space)
     ` defClustal: is clustal format (* space)
     */

   signed long minScoreSL; /*Min score to keep alignment*/

   /*scoring variables*/
   signed char noGapBl;   /*1: use gapExtendC; 0 do not*/
   signed short gapSS;    /*Penalty for starting indel*/
   signed short extendSS; /* Penalty for extending indel*/

   signed short delArySS[4]; /*gap penalties*/
   signed short insArySS[4]; /*gap penalties*/

   signed short scoreMatrixSS[
         def_lenCol_alnSet
      ][
         def_lenCol_alnSet
      ]; /*scoring matrix*/

   signed char matchMatrixSC[
         def_lenCol_alnSet
      ][
         def_lenCol_alnSet
      ]; /*match matrix*/

     /* Size of both non-WORDS matrixies (27 or 26) is due
     ` to wanting a look up table that can handle
     ` anonymous bases.  Most cells will be set to 0. 
     */
}alnSet;

/*-------------------------------------------------------\
| Fun01: setScore_alnSet
|  - Sets the score for a base pair (reference/query)
| Input:
|  - qryBase:
|    o query base to set score fore
|  - refBase:
|    o reference base to set score fore
|  - score:
|    o Score to have for this query/reference combination
|      (max is a short)
|  - alnSetSTPtr:
|    o pointer to alnSet (settings) structer with matrix
|      to add the new query reference comparison score to.
| Output:
|  o Modifies:
|    - one score in an snp/match scoring matrix
|  - default
|    o This assumes that the sequences are converted to
|      indexes
\-------------------------------------------------------*/
#define setScore_alnSet(qryBase, refBase, score, alnSetSTPtr) ((alnSetSTPtr)->scoreMatrixSS[(unsigned char) ((qryBase) & def_ntToCode_alnSet)][(unsigned char) ((refBase) & def_ntToCode_alnSet)] = (score))
/*setScore_alnSet, sequences will be lookup indexes*/

/*-------------------------------------------------------\
| Fun02: setMatch_alnSet
|  - Sets if two bases are a match or mismtach
| Input:
|  - qryBase:
|    o query base to set score fore
|  - refBase:
|    o reference base to set score fore
|  - match:
|    o 1: the query and reference bases are matches
|    o 0: the query and reference bases are mismatches
|  - alnSetSTPtr:
|    o pointer to alnSet (settings) structer with matrix
|      to add the new query reference comparison score to.
| Output:
|  - Modifies:
|    o one match in the matchMatrix
| Variations:
|  - default
|    o This assumes that the sequences are converted to
|      indexes
\-------------------------------------------------------*/
/*The user will not be providing converted bases*/
#define setMatch_alnSet(qryBase, refBase, match, alnSetSTPtr) ((alnSetSTPtr)->matchMatrixSC[(unsigned char) ((qryBase) & def_ntToCode_alnSet)][(unsigned char) ((refBase) & def_ntToCode_alnSet)] = (match))
/*setBasePairMatch, sequences will be lookup indexes*/


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
);

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
);

/*-------------------------------------------------------\
| Fun05: getScore_alnSet:
|  - Get the score for a pair of bases from an alignment
| Input:
|  - qryBase:
|    o character with query base to get score for
|  - refBase:
|    o character with reference base to get score for
|  - alnSetPtr:
|    o alnSet (settings) structer pionter with scoring
|      matrix
| Output:
|  - Returns:
|    o score of the input base pair
|  - default
|    o This assumes that the sequences are converted to
|      indexes
\-------------------------------------------------------*/
#define getScore_alnSet(qryBase, refBase, alnSetSTPtr) ((alnSetSTPtr)->scoreMatrixSS[(unsigned char) (qryBase)][(unsigned char) (refBase)])
/*getBasePairScore, with sequence converted to index*/

/*-------------------------------------------------------\
| Fun06: getMatch_alnSet
|  - Check if two bases were a match or mismatch
| Input:
|  - qryBase:
|    o character with query base to compare
|  - refBase:
|    o character with reference base to compare to
|  - alnSetPtr:
|    o alnSet (settings) structer pionter with scoring
|      matrix
| Output:
|  - Returns:
|    o def_ntEql_alnDefs if the bases matched (same)
|    o def_ntNotEql_alnDefs if bases were not a match
|    o def_anonymous_alnDefs for anonymous bases
|    o def_anonymous_alnDefs | def_ntEql_alnDefs for
|      anonymous matches (anoymous base could be same)
|  - default
|    o This assumes that the sequences are converted to
|      indexes
\-------------------------------------------------------*/
#define getMatch_alnSet(qryBase, refBase, alnSetSTPtr) ((alnSetSTPtr)->matchMatrixSC[(unsigned char) (qryBase)][(unsigned char) (refBase)])
/*getBasePairMatch, with sequence converted to index*/

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
    void *scoreFILE  /*File scoring matrix scores*/
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

#endif

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
