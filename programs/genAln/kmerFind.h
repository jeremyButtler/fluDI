/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' kmerFind SOF: Start Of File
'   - holds functions for kmer detection of primers
'   o header:
'     - guards, forward declerations, & defined variables
'   o tbl01: alnNtToBit_kmerFind
'     - converts an nucleotide alignment code form
'       alnSetStruct.h to an two bit value, with an 3rd
'       bit being used for anonymous bases and the 4th bit
'       for errors
'   o .h st01: tblST_kmerFind
'     - holds the kmer tables for detecting spoligytpes
'   o .h st02: refST_kmerFind
'     - holds the kmer pattern for the reference
'   o fun01: blank_tblST_kmerFind
'     - blanks all stored values in an tblST_kmerFind
'   o fun02: init_tblST_kmerFind
'     - initializes an tblST_kmerFind structure
'   o fun03: setup_tblST_kmerFind
'     - allocates memory for variables in a tblST_kmerFind
'       structure
'   o fun04: freeStack_tblST_kmerFind
'     - frees all variables in an tblST_kmerFind structure
'   o fun05: freeHeap_tblST_kmerFind
'     - frees and tblST_kmerFind structure
'   o fun06: blank_refST_kmerFind
'     - sets the counting and kmer values in an
'       refST_kmerFind to 0 or def_noKmer_kmerFind
'   o fun07: init_refST_kmerFind
'     - initializes an refST_kmerFind structure
'   o fun08: setup_refST_kmerFind
'     - allocates memory for a refST_kmerFind structure
'   o fun09: freeStack_refST_kmerFind
'     - frees the variables in an refST_kmerFind structure
'   o fun10: freeHeap_refST_kmerFind
'     - frees an refST_kmerFind structure
'   o fun11: freeHeapAry_refST_kmerFind
'     - frees an array of refST_kmerFind structure
'   o fun12: uiTwinShort_kmerFind
'     - sorts a unsigned int array from least to greatest
'       and keeps it linked to a second unsigned int array
'   o fun13: addSeqToRefST_kmerFInd
'     - adds a sequence to a refST_kmerFind structure
'   o fun14: prep_tblST_kmerFind
'     - sets up an tblST_kmerFind structure for primer
'       searching
'   o fun15: tsvToAry_refST_kmerFind
'     - makes an array of refST_kmerFind structures from a
'       tsv file
'   o fun16: faToAry_refST_kmerFind
'     - makes an array of refST_kmerFind structures
'   o fun17: nextSeqChunk_tblST_kmerFind
'     - adds a new set of kmers from an sequence to an
'       tblST_kmerFind structure
'   o fun18: forCntMatchs_kmerFind
'     - finds the number of kmers that are in both the
'       kmer table (query) and the pattern (reference)
'   o fun19: revCntMatchs_kmerFind
'     - finds the number of kmers that are shared in the
'       kmer table (query) and the reverse pattern
'       (reference)
'   o fun20: matchCheck_kmerFind
'     - tells if the  match meets the min requirements to
'       do an alignment or not
'   o fun21: findRefInChunk_kmerFind
'     - does an kmer check and alings an single sequence
'       in an refST_kmerFind structure to see if there is
'       an match
'   o fun22: waterFindPrims_kmerFind
'     - finds primers in an sequence (from fastx file)
'       using a slower, but more percise waterman
'   o fun23: fxFindPrims_kmerFind
'     - finds spoligotype spacers in an sequence (from
'       fastx file) using an faster kmer search followed
'       by an slower waterman to finalize alignments
'   o fun24: phit_kmerFind
'     - prints out the primer hits for a sequence
'   o fun25: pHeaderHit_kmerFind
'      - prints header for phit_kmerFind (fun24)
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards, forward declerations, and defined variables
\-------------------------------------------------------*/

#ifndef KMER_FIND_H
#define KMER_FIND_H

typedef struct seqST seqST;
typedef struct alnSet alnSet;

#define def_noKmer_kmerFind -1
#define def_endKmers_kmerFind -2
#define def_bitsPerKmer_kmerFind 2 /*do not change*/

#define def_noPrim_kmerFind 1
#define def_fileErr_kmerFind 2
#define def_memErr_kmerFind 4

#define def_minPercScore_kmerFind 0.85f
   /*this setting here seems to be the biggest deal*/

/*this is about twice as fast, but also misses half the
`   results
*/
#define def_minKmerPerc_kmerFind 0.4f
#define def_percShift_kmerFind 0.5f
#define def_extraNtInWin_kmerFind 1.0f
#define def_lenKmer_kmerFind 5
   /*3mers are slower than the waterman alone*/

/*-------------------------------------------------------\
| ST01: tblST_kmerFind
|   - holds the kmer tables for detecting spoligytpes
\-------------------------------------------------------*/
typedef struct tblST_kmerFind
{
   unsigned char lenKmerUC; /*number nucleotides in kmer*/

   signed int *tblSI; /*table (array) of kmer counts*/
   unsigned int lenTblUI; /*length of tblSI*/

   signed int *kmerArySI; /*converted kmer sequence*/
                          /*make sure there is always one
                          `   extra element for the ending
                          `  -2
                          */
   unsigned int numKmerUI; /*number of kmers kmerArySI
                           `  should store
                           */

   unsigned int ntInWinUI; /*number kmers in window*/
   unsigned int rmNtUI; /*number bases to remove per add*/
   unsigned long kmerMaskUL; /*masks extra bits in long*/

   unsigned long lastKmerUL;    /*last kmer in array*/
   unsigned long lenLastKmerUL;
      /*number bases since last anonymous base*/

   unsigned long seqPosUL; /*position at in sequence*/
   struct seqST *seqSTPtr; /*sequence working on*/
}tblST_kmerFind;

/*-------------------------------------------------------\
| ST02: refST_kmerFind
|   - holds the kmer pattern for the reference
\-------------------------------------------------------*/
typedef struct refST_kmerFind
{
   unsigned char lenKmerUC;  /*bases in one kmer*/
   unsigned int minKmersUI;  /*min matches to count*/

   signed int *forKmerArySI;/*array of forward kmers*/
   unsigned int *forRepAryUI;/*number times kmer repeats*/

   signed int *revKmerArySI;/*array of reverse kmers*/
   unsigned int *revRepAryUI;/*number times kmer repeats*/

   float maxForScoreF;       /*maximum score possible*/
   float maxRevScoreF;       /*maximum score possible*/
   unsigned int lenAryUI;    /*length of both arrays*/

   struct seqST *forSeqST;
   struct seqST *revSeqST;   /*reverse sequence*/
      /*the logic here is that I am only counting the
      `  number of kmers and that I am doing short
      `  sequences. so, it is ok to just count the number
      `  of times each kmer repeats
      */

   signed int mateSI; /*matching primer*/
}refST_kmerFind;

/*-------------------------------------------------------\
| Fun01: blank_tblST_kmerFind
|   - blanks all stored values in an tblST_kmerFind
|     structure
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to blank
|   - blankSeqBl:
|     o 1: blank the seqST (sequence) in tblSTPtr
|     o 0: do not blank the seqST in tblSTPtr
| Output:
|   - Modifies:
|     o tblSI in tlbSTPtr to be full of zeros
|     o kmerArySI to be full of -1's and to end in -2
|     o lenLastKmer to be 0
|     o seqPosUL to be 0
|     o if requested seqST with blaknSeqST from
|       ../memwater/seqST.h
\-------------------------------------------------------*/
void
blank_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   signed char blankSeqBl
);

/*-------------------------------------------------------\
| Fun02: init_tblST_kmerFind
|   - initializes a tblST_kmerFind structure
| Input:
|   - tblSTPtr:
|     o pointer to tblST_kmerFind structure to initialize
| Output:
|   - Modifies:
|     o all varaibles in an tblST_kmerFind find structure
|       to be 0
\-------------------------------------------------------*/
void
init_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr
);

/*-------------------------------------------------------\
| Fun03: setup_tblST_kmerFind
|   - allocates memory for variables in a tblST_kmerFind
|     structure
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to blank
|   - lenKmerUC:
|     o length of one kmer
| Output:
|   - Modifies:
|     o sets lenKmerUC and allocates memory for arrays
|       and seqST structures
|   - Returns:
|     o 0 for no errors
|     o def_memErr_kmerFind for an memory error
\-------------------------------------------------------*/
signed char
setup_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   unsigned char lenKmerUC
);

/*-------------------------------------------------------\
| Fun04: freeStack_tblST_kmerFind
|   - frees all variables in an tblST_kmerFind structure
|     and sets other values to defaults (calls blank)
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to free
|       variables in
| Output:
|   - Frees:
|     o all allocated varialbes in tlbStPtr
\-------------------------------------------------------*/
void
freeStack_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr
);

/*-------------------------------------------------------\
| Fun05: freeHeap_tblST_kmerFind
|   - frees and tblST_kmerFind structure
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to free
| Output:
|   - Frees:
|     o tblSTPtr (does not set to null)
\-------------------------------------------------------*/
void
freeHeap_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr
);

/*-------------------------------------------------------\
| Fun06: blank_refST_kmerFind
|   - sets the counting and kmer values in an
|     refST_kmerFind to 0 or def_noKmer_kmerFind
| Input:
|   - refSTPtr:
|     o pointer to an refST_kmerFind structure to blank
|       variables in
| Output:
|   - Modifies:
|     o forKmerArySI to be full of -1's
|     o forRefAryUS to be full of 0's
|     o revKmerArySI to be full of -1's
|     o revRefAryUS to be full of 0's
|     o forSeqST to be blanked
|     o revSeqST to be blanked
\-------------------------------------------------------*/
void
blank_refST_kmerFind(
   struct refST_kmerFind *refSTPtr
);

/*-------------------------------------------------------\
| Fun07: init_refST_kmerFind
|   - initializes an refST_kmerFind structure
| Input:
|   - refSTPtr:
|     o pointer to refST_kmerFind structure to initialize
| Output:
|   - Modifies:
|     o refSTPtr to be initialized
\-------------------------------------------------------*/
void
init_refST_kmerFind(
   struct refST_kmerFind *refSTPtr
);

/*-------------------------------------------------------\
| Fun08: setup_refST_kmerFind
|   - allocates memory for a refST_kmerFind structure
| Input:
|   o refSTPtr,
|     - pointer to a refST_kmerFind structure to setup
|   o lenKmerUC:
|     - length of one kmer
| Output:
|   - Modifies:
|     o refSTPtr to have memory allocated to variables and
|       kmer size
|   - Returns:
|     o 0 if no errors
|     o def_memErr_kmerFind for memory errors
\-------------------------------------------------------*/
signed char
setup_refST_kmerFind(
   struct refST_kmerFind *refSTPtr,
   unsigned char lenKmerUC
);

/*-------------------------------------------------------\
| Fun09: freeStack_refST_kmerFind
|   - initializes an refST_kmerFind structure
| Input:
|   - refSTPtr:
|     o pointer to refST_kmerFind structure to free
|       varialbes in
| Output:
|   - Frees:
|     o all pointers in refSTPtr
|   - Sets
|     o pointers in refSTPtr to 0
\-------------------------------------------------------*/
void
freeStack_refST_kmerFind(
   struct refST_kmerFind *refSTPtr
);

/*-------------------------------------------------------\
| Fun10: freeHeap_refST_kmerFind
|   - frees an refST_kmerFind structure
| Input:
|   - refSTPtr:
|     o pointer to refST_kmerFind structure to free
|       varialbes in
| Output:
|   - Frees:
|     o refSTPtr
\-------------------------------------------------------*/
void
freeHeap_refST_kmerFind(
   struct refST_kmerFind *refSTPtr
);

/*-------------------------------------------------------\
| Fun11: freeHeapAry_refST_kmerFind
|   - frees an array of refST_kmerFind structure
| Input:
|   - refSTAry:
|     o pointer to refST_kmerFind structure array to free
|   - lenArySI:
|     o number of refST_kmerFind structures in refSTPtr
| Output:
|   - Frees:
|     o all refST_kmerFind structures in refSTAry
|     o refSTAry (does not set to 0)
\-------------------------------------------------------*/
void
freeHeapAry_refST_kmerFind(
   struct refST_kmerFind *refSTAry,
   signed int lenArySI
);

/*-------------------------------------------------------\
| Fun12: uiTwinShort_kmerFind
|   - sorts a unsigned int array from least to greatest
|     and keeps it linked to a second unsigned int array
| Input:
|   - uiAry:
|     o array to sort
|   - uiSecAry:
|     o second array to keep in order with uiAry
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o uiAry to be sorted form least to greatest
|     o uiSecAry to be sorted by uiAry
\-------------------------------------------------------*/
void
uiTwinSort_kmerFind(
   unsigned int *uiAry,
   unsigned int *uiSecAry,
   unsigned long startUL,
   unsigned long endUL
);

/*-------------------------------------------------------\
| Fun13: addSeqToRefST_kmerFInd
|   - adds a sequence to a refST_kmerFind structure
| Input:
|   - tblSTPtr:
|     o pointer to a tblST_kmerFind structure with
|       settings, such as the kmer length, mask, and
|       maximum number of kmers
|   - refSTPtr:
|     o pionter to the refST_kmerFind structure to add the
|       sequence to
|   - seqSTPtr:
|     o pointer to the seqST to get the sequence from
|   - minPercKmersF:
|     o float with minimum percentage of kmers to start
|       considering an window supports an spacer
|   - longestSeqUI:
|     o length of the longest sequence in a refSTPtr
|       structure
|   - alnSetPtr:
|     o pointer to alnSet structure with score matrix
|       (for finding max possible score)
| Output:
|   - Returns:
|     o 0 for memory error
|     o length of longest sequence in an refST_kmerFind
|       structure
\-------------------------------------------------------*/
unsigned int
addSeqToRefST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   struct refST_kmerFind *refSTPtr,
   struct seqST *seqSTPtr,
   float minPercKmersF,
   unsigned int longestSeqUI,
   struct alnSet *alnSetPtr
);

/*-------------------------------------------------------\
| Fun14: prep_tblST_kmerFind
|   - sets up an tblST_kmerFind structure for primer
|     searching
| Input:
|   - tblSTPtr:
|     o pointer to a tblST_kmerFind structure to set up
|   - percExtraNtInWinF:
|     o float with percentage of extra nucleotides to
|       store in one window (beyond reference length)
|   - percWinShiftF:
|     o float with percentage of bases to shift for each
|       new window in tblSTPtr
|   - longestSeqUI:
|     o longest sequence to map against. Used to find the
|       maximum window size
| Output:
|   - Returns:
|     o 0 for no errors
|     o def_memErr_kmerFind for memory errors
\-------------------------------------------------------*/
signed char
prep_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   float percExtraNtInWinF,
   float percWinShiftF,
   unsigned long longestSeqUI
);

/*-------------------------------------------------------\
| Fun15: tsvToAry_refST_kmerFind
|   - makes an array of refST_kmerFind structures from a
|     tsv file
| Input:
|   - tsvFileStr:
|     o c-string with path to tsv file with sequences
|       to build an refST_kmerFind array from
|       - column 1: primer id
|       - column 2: paired (True, Yes, 1) or unpaired
|       - column 3: forward primer
|       - column 4: reverse primer
|   - lenKmerUC:
|     o length of one kmer
|   - lenArySIPtr:
|     o will hold the number of refST_kmerFind structures
|       made
|   - minPercKmersF:
|     o float with minimum percentage of kmers to start
|       considering an window supports an spacer
|   - maxNtUIPtr:
|     o pointer to unsigned int to hold the largest
|       reference sequence
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to link
|       with this reference set
|   - percExtraNtInWinF:
|     o float with percentage of extra nucleotides to
|       store in one window (beyond reference length)
|   - percWinShiftF:
|     o float with percentage of bases to shift for each
|       new window in tblSTPtr
|   - alnSetSTPtr:
|     o pointer to an alnSet structure with the alignment
|       settings (for finding max possible score)
|   - errSC:
|     o holds error type
| Output:
|   - Modifies:
|     o lenArySIPtr to hold the number of sequences in
|       refST_kmerFind array
|     o errSC to hold the error type
|       - 0 for no errors
|       - def_memErr_kmerFind for memory errors
|       - def_fileErr_kmerFind for file errors
|     o ntInWinUI in tblSTPtr to have the number of bases
|       in one window (2 * longest primer length)
|     o rmNtUI in tblSTPtr to have the number of bases
|       to remove with each window shift
|       (percWinShiftF * tblSTPtr->ntINWinUI)
|   - Returns:
|     o an array of refST_kmerFind structures with
|       sequences to check
|     o 0 for an error
\-------------------------------------------------------*/
refST_kmerFind *
tsvToAry_refST_kmerFind(
   signed char *tsvFileStr,
   unsigned char lenKmerUC,
   signed int *lenArySIPtr,
   float minPercKmersF,
   struct tblST_kmerFind *tblSTPtr,
   float percExtraNtInWinF,
   float percWinShiftF,
   struct alnSet *alnSetSTPtr,
   signed char *errSC
);

/*-------------------------------------------------------\
| Fun16: faToAry_refST_kmerFind
|   - makes an array of refST_kmerFind structures
| Input:
|   - faFileStr:
|     o c-string with path to fasta file with sequences
|       to build an refST_kmerFind array from
|   - lenKmerUC:
|     o length of one kmer
|   - lenArySIPtr:
|     o will hold the number of refST_kmerFind structures
|       made
|   - minPercKmersF:
|     o float with minimum percentage of kmers to start
|       considering an window supports an spacer
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to link
|       with this reference set
|   - percExtraNtInWinF:
|     o float with percentage of extra nucleotides to
|       store in one window (beyond reference length)
|   - percWinShiftF:
|     o float with percentage of bases to shift for each
|       new window in tblSTPtr
|   - alnSetSTPtr:
|     o pointer to an alnSet structure with the alignment
|       settings (for finding max possible score)
|   - errSC:
|     o holds error type
| Output:
|   - Modifies:
|     o lenArySIPtr to hold the number of sequences in
|       th refST_kmerFind array
|     o errSC to hold the error type
|       - 0 for no errors
|       - def_memErr_kmerFind for memory errors
|       - def_fileErr_kmerFind for file errors
|     o ntInWinUI in tblSTPtr to have the number of bases
|       in one window (2 * longest primer length)
|     o rmNtUI in tblSTPtr to have the number of bases
|       to remove with each window shift
|       (percWinShiftF * tblSTPtr->ntINWinUI)
|   - Returns:
|     o an array of refST_kmerFind structures with
|       sequences to check
|     o 0 for an error
\-------------------------------------------------------*/
refST_kmerFind *
faToAry_refST_kmerFind(
   signed char *faFileStr,
   unsigned char lenKmerUC,
   signed int *lenArySIPtr,
   float minPercKmersF,
   struct tblST_kmerFind *tblSTPtr,
   float percExtraNtInWinF,
   float percWinShiftF,
   struct alnSet *alnSetSTPtr,
   signed char *errSC
);

/*-------------------------------------------------------\
| Fun17: nextSeqChunk_tblST_kmerFind
|   - adds a new set of kmers from an sequence to an
|     tblST_kmerFind structure
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure to add
|       kmers to
|   - firstTimeBl:
|     o 1: first time adding sequence (blank kmer array)
|     o 0: updating the kmer window
| Output:
|   - Modifies:
|     o tblSI and seqAryUS in tblSTPtr to have the old
|       kmers (number specified by rmNtUI in tblSI)
|       remove and the new kmers added in
|       - for end of sequence it sets an index to
|         def_endKmers_kmerFind
|    o firstTimeBl:
|      o to be 0 if it is 1
|   - Returns:
|     o 0 for not end of sequence
|     o 1 for end of sequence
\-------------------------------------------------------*/
signed char
nextSeqChunk_tblST_kmerFind(
   struct tblST_kmerFind *tblSTPtr,/*table to add seq to*/
   signed char *firstTimeBl /*1: first set of kmers*/
);

/*-------------------------------------------------------\
| Fun18: forCntMatchs_kmerFind
|   - finds the number of kmers that are in both the
|     kmer table (query) and the pattern (reference)
| Input:
|   - tblST_kmerFindPtr:
|     o pointer to an tblST_kmerFind structure with an
|       kmer table to get the number of matches from
|   - refST_kmerFindPtr:
|     o pointer to an refST_kmerFind structure with an
|       array of kmers to get the number of matches for
| Output:
|   - Returns:
|     o number of matching kmers
\-------------------------------------------------------*/
signed long
forCntMatchs_kmerFind(
   struct tblST_kmerFind *tblST_kmerFindPtr,
   struct refST_kmerFind *refST_kmerFindPtr
);

/*-------------------------------------------------------\
| Fun19: revCntMatchs_kmerFind
|   - finds the number of kmers that are shared in the
|     kmer table (query) and the reverse pattern
|     (reference)
| Input:
|   - tblST_kmerFindPtr:
|     o pointer to an tblST_kmerFind structure with an
|       kmer table to get the number of matches from
|   - refST_kmerFindPtr:
|     o pointer to an refST_kmerFind structure with an
|       array of kmers to get the number of matches rev
| Output:
|   - Returns:
|     o number of matching kmers
\-------------------------------------------------------*/
signed long
revCntMatchs_kmerFind(
   struct tblST_kmerFind *tblST_kmerFindPtr,
   struct refST_kmerFind *refST_kmerFindPtr
);

/*-------------------------------------------------------\
| Fun20: matchCheck_kmerFind
|   - tells if the  match meets the min requirements to
|     do an alignment or not
| Input:
|   - tblST_kmerFindPtr:
|     o pointer to an tblST_kmerFind structure with the
|       chunk of query (kmer table) to check
|   - refST_kmerFindPtr:
|     o pointer to an refST_kmerFind structure with the
|       reference (primers) kmers to check
| Output:
|   - Returns:
|     o 0 (foward) or 1 (best) if there is no match
|     o 1 | 0 if match and the forward direction is best
|     o 1 | 2 if match and the reverse direction is best
\-------------------------------------------------------*/
signed char
matchCheck_kmerFind(  
   struct tblST_kmerFind *tblST_kmerFindPtr, 
   struct refST_kmerFind *refST_kmerFindPtr  
); 

/*-------------------------------------------------------\
| Fun21: findRefInChunk_kmerFind
|   - does an kmer check and alings an single sequence
|     in an refST_kmerFind structure to see if there is
|     an match
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure with the
|       chunk of query (kmer table) to check
|   - refSTPtr:
|     o pointer to an refST_kmerFind structure with the
|       reference (primers) kmers to check
|   - alnSetSTPtr:
|     o pointer to an alnSet structure with the alignment
|       settings
|   - minPerScoreF:
|     o float with minimum percent score to keep an
|       alingment
|   - scoreSL:
|     o pointer to an signed long to hold the alingment
|       score
|   - qryStartUL:
|     o pointer to an unsigned long to hold the first base
|       in the query alingment
|   - qryEndUL:
|     o pointer to an unsigned long to hold the last base
|       in the query alingment
|   - refStartUL:
|     o pointer to an unsigned long to hold the first base
|       in the reference alingment
|   - refEndUL:
|     o pointer to an unsigned long to hold the last base
|       in the reference alingment
| Output:
|   - Modifies:
|     o scoreSL
|       - 0 if no alignment done
|       - score if an alignment was done
|     o qryStartUL
|       - 0 if no alignment done
|       - first aligned query base if alignment done
|     o qryEndtUL
|       - 0 if no alignment done
|       - last aligned query base if alignment done
|     o refStartUL
|       - 0 if no alignment done
|       - first aligned reference base if alignment done
|     o refEndtUL
|       - 0 if no alignment done
|       - last aligned reference base if alignment done
|   - Returns:
|     o 1 if the reference sequence was found in the
|       kmer table (query) sequence
|     o 0 if reference sequence not found
\-------------------------------------------------------*/
signed char
findRefInChunk_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   struct refST_kmerFind *refSTPtr,
   struct alnSet *alnSetSTPtr,
   float minPercScoreF,
   signed long *scoreSL,
   unsigned long *qryStartUL,
   unsigned long *qryEndUL,
   unsigned long *refStartUL,
   unsigned long *refEndUL
);

/*-------------------------------------------------------\
| Fun22: waterFindPrims_kmerFind
|   - finds primers in an sequence (from fastx file) using
|     a slower, but more percise waterman
| Input:
|   - refAryST
|     o array of refST_kmerFind structures with reference
|       (primer) sequences to search for
|   - lenRefAryUI:
|     o number of refST_kmerFind structures in refSTAry
|   - seqSTPtr:
|     o pointer to an seqST structer with the query
|       sequence
|   - minPerScoreF:
|     o float with minimum percent score to keep an
|       alingment
|   - codeAryUI:
|     o unsigned integer array to hold the number of times
|       each spacer was detected
|   - dirArySC:
|     o pointer to a signed char array to hold mapped
|       primer directions
|   - scoreArySL:
|     o array of signed longs with the best score for each
|       matched primer
|   - seqStartAryUL:
|     o array of unsigned longs with the starting position
|       one the sequence for each score in scoreArySL
|   - endAyUL:
|     o array of unsigned longs with the ending position
|       on the sequence for each score in scoreArySL
|   - primStartAryUL:
|     o array of unsigned longs with the starting position
|       on the primer for each score in scoreArySL
|   - primEndAyUL:
|     o array of unsigned longs with the ending position
|       on the primer for score in scoreArySL
|   - alnSetSTPtr:
|     o pointer to an alnSet structure with the alignment
|       settings
| Output:
|   - Modifies:
|     o each position codeStr
|       - 1 if detected at one primer
|       - 2 if detected at both primers
|       - 0 if no primers were found
|     o dirAryStr to have direction of best alignment
|       - F for foward
|       - R for reverse
|     o scoreArySL score of best alignment for each primer
|     o seqStartAryUL starting sequence position of best
|       alignment for each primer
|     o seqEndAryUL ending sequence position of best
|       alignment for each primer
|     o primgStartAryUL starting primer position of best
|       alignment for each primer
|     o primgEndAryUL ending primer position of best
|       alignment for each primer
|   - Returns:
|     o 0 if the reference sequence was found in the kmer
|       table (query) sequence
|     o def_noPrim_kmerFind if no primers were found
\-------------------------------------------------------*/
signed char
waterFindPrims_kmerFind(
   struct refST_kmerFind *refSTAry,
   unsigned int lenRefAryUI,
   struct seqST *seqSTPtr,
   float minPercScoreF,
   unsigned int codeAryUI[],
   signed char dirArySC[],
   signed long scoreArySL[],
   unsigned long seqStartAryUL[],
   unsigned long seqEndAryUL[],
   unsigned long primStartAryUL[],
   unsigned long primEndAryUL[],
   struct alnSet *alnSetSTPtr
);

/*-------------------------------------------------------\
| Fun23: fxFindPrims_kmerFind
|   - finds primers in an sequence (from fastx file) using
|     an faster kmer search followed by an slower waterman
|     to finalize alignments
| Input:
|   - tblSTPtr:
|     o pointer to an tblST_kmerFind structure with
|       settings and to hold the query sequence to check
|   - refAryST
|     o array of refST_kmerFind structures with reference
|       (primer) sequences to search for
|   - lenRefAryUI:
|     o number of refST_kmerFind structures in refSTAry
|   - seqSTPtr:
|     o pointer to an seqST structure with the
|       sequence to check for primers in
|   - minPerScoreF:
|     o float with minimum percent score to keep an
|       alingment
|   - codeAryUI:
|     o unsigned integer array to hold the number of times
|       each spacer was detected
|   - dirArySC:
|     o pointer to a signed char array to hold mapped
|       primer directions
|   - scoreArySL:
|     o array of signed longs with the best score for each
|       matched primer
|   - seqStartAryUL:
|     o array of unsigned longs with the starting position
|       one the sequence for each score in scoreArySL
|   - endAyUL:
|     o array of unsigned longs with the ending position
|       on the sequence for each score in scoreArySL
|   - primStartAryUL:
|     o array of unsigned longs with the starting position
|       on the primer for each score in scoreArySL
|   - primEndAyUL:
|     o array of unsigned longs with the ending position
|       on the primer for score in scoreArySL
|   - alnSetSTPtr:
|     o pointer to an alnSet structure with the alignment
|       settings
| Output:
|   - Modifies:
|     o each position codeStr to have number of times I
|       detected each primer
|     o dirAryStr to have direction of best alignment
|       - F for foward
|       - R for reverse
|     o scoreArySL score of best alignment for each primer
|     o seqStartAryUL starting sequence position of best
|       alignment for each primer
|     o seqEndAryUL ending sequence position of best
|       alignment for each primer
|     o primgStartAryUL starting primer position of best
|       alignment for each primer
|     o primgEndAryUL ending primer position of best
|       alignment for each primer
|   - Returns:
|     o 0 if the reference sequence was found in the kmer
|       table (query) sequence
|     o def_noPrim_kmerFind if no primers were found
\-------------------------------------------------------*/
signed char
fxFindPrims_kmerFind(
   struct tblST_kmerFind *tblSTPtr,
   struct refST_kmerFind *refSTAry,
   unsigned int lenRefAryUI,
   struct seqST *seqSTPtr,
   float minPercScoreF,
   unsigned int codeAryUI[],
   signed char dirArySC[],
   signed long scoreArySL[],
   unsigned long seqStartAryUL[],
   unsigned long seqEndAryUL[],
   unsigned long primStartAryUL[],
   unsigned long primEndAryUL[],
   struct alnSet *alnSetSTPtr
);

/*-------------------------------------------------------\
| Fun24: phit_kmerFind
|   - prints out the primer hits for a sequence
| Input:
|   - refAryST:
|     o Array of refST_kmerFind structures with the primer
|       sequences
|   - numRefsSI:
|     o number of refST_kmerFind (primers) structures in
|       refAryST
|   - seqSTPtr:
|     o pointer to an seqST structer with the read id
|   - codeAryUI:
|     o pointer to a unsigned int array with the number of
|       times the sequence had each primer
|   - dirArySC:
|     o pointer to a signed char array with the mapped
|       direction for each primer
|   - scoreArySL:
|     o pointer to a signed long array with the score bes
|       score for each primer
|   - seqStartArySL:
|     o pointer to a unsigned long array with the primers
|       starting position on the sequence for the best
|       score
|   - seqEndArySL:
|     o pointer to a unsigned long array with the primers
|       ending position on the sequence for the best
|       score
|   - primStartArySL:
|     o pointer to a unsigned long array with the
|       sequences starting position on the primer for the
|       best score
|   - primEndArySL:
|     o pointer to a unsigned long array with the
|       sequences ending position on the primer for the
|       best score
|   - outFILE:
|     o file to print the stats to
\-------------------------------------------------------*/
void
phit_kmerFind(
   struct refST_kmerFind *refAryST,
   signed int numRefsSI,
   struct seqST *seqSTPtr,
   unsigned int *codeAryUI,
   signed char *dirArySC,
   signed long *scoreArySL,
   unsigned long *seqStartAryUL,
   unsigned long *seqEndAryUL,
   unsigned long *primStartAryUL,
   unsigned long *primEndAryUL,
   void *outFILE
);

/*-------------------------------------------------------\
| Fun25: pHeaderHit_kmerFind
|    - prints header for phit_kmerFind (fun24)
| Input:
|   - outFILE:
|     o file to print header to
| Output:
|   - Prints:
|     o header to outFILE
\-------------------------------------------------------*/
void
pHeaderHit_kmerFind(
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
