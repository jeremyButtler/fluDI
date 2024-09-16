/*########################################################
# Name: diScan
#   - has functions to scan a sequence to see if it has DI
#     events (large deletions)
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - included libraries
'   o fun01: findSeg_diScan
'     - uses kmer profiling to find the best matching
'       segment
'   o fun02: waterScan_diScan
'     - scan for DI fragments using kmer profiling and
'       Watermen alignment
'   o fun03: phead_diScan
'     - print out header for fragment tsv
'   o fun04: pfrag_diScan
'     - print out the fragment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - forward declerations, defined variables, and guards
\-------------------------------------------------------*/

#ifndef DEFECTIVE_INFLUENZA_SCAN_H
#define DEFECTIVE_INFLUENZA_SCAN_H

typedef struct seqST seqST;
typedef struct samEntry samEntry;
typedef struct dirMatrix dirMatrix;
typedef struct alnSet alnSet;
typedef struct kmerCnt kmerCnt;

/*these should always be less than 0*/
#define def_noMatch_diScan -1
#define def_memErr_diScan -2

#define def_matchScore_diScan 5
   /*This appears to disagree with alnDefs, but after
   `   division by 10, this is true
   */

/*-------------------------------------------------------\
| Fun01: findSeg_diScan
|   - uses kmer profiling to find the best matching
|     segment
| Input:
|   - seqSTPtr:
|     o pointer to seqST with sequence to find segment
|       for
|   - refAryST:
|     o pointer to kmerCnt structure array with segments
|       to compare againts
|   - numSegUI:
|     o number segments (references) in refAryST
|   - kmerArySI:
|     o signed int array to hold the kmers in seqSTPtr
|   - cntArySI:
|     o signed int array to hold kmer counts for seqSTPtr
|   - lenKmerUC:
|     o length of one kmer
|   - minKmerPercF:
|     o min percent of shared kmers to keep a read
|   - maxKmersSI:
|     o pointer to signed long to hold number of kmers
|       that mapped (- for reverse, + for forward)
| Output:
|   - Modifies:
|     o maxKmersSI to have the number of matches
|       - is positive if foward reference was best
|       - is negative if reverse reference was best
|   - Returns:
|     o -1 if no segment mapped
|     o > 0 (segment index) if a segment mapped
\-------------------------------------------------------*/
signed int
findSeg_diScan(
   struct seqST *seqSTPtr, /*sequence to check*/
   struct kmerCnt *refAryST,  /*array of references*/
   unsigned int numSegUI,      /*length of ref array*/
   signed int *kmerArySI,      /*holds sequence kmers*/
   signed int *cntArySI,       /*holds kmer counts*/
   unsigned char lenKmerUC,    /*length of one kmer*/
   float minKmerPercF,         /*min perc kmers to keep*/
   signed int *maxKmersSI      /*will hold kmer count*/
);

/*-------------------------------------------------------\
| Fun02: waterScan_diScan
|   - scan for DI fragments using kmer profiling and
|     Watermen alignment
| Input:
|    - seqSTPtr:
|      o pointer to seqST with sequence to scan
|   - refAryST:
|     o pointer to kmerCnt structure array with segments
|       to compare againts
|   - lenRefUI:
|     o number of kmerCnt structures to scan in refAryST
|   - kmerArySI:
|     o signed int array to hold the kmers in seqSTPtr
|   - cntArySI:
|     o signed int array to hold kmer counts for seqSTPtr
|   - minPercScoreF:
|     o minimum percent score from waterman alingment to
|       check for DI (count as mapped)
|   - minKmerPercF:
|     o min percent of shared kmers to keep a read
|   - lenKmerUC:
|     o length of one kmer
|   - samSTPtr:
|     o pointer to samEntry structure to hold aligned
|       sequence
|   - segSIPtr:
|     o pointer to signed int ot hold the segment mapped
|       to or a negative number
|   - minDIDelUI:
|     o minimum deletion size to flag as a DI sequence
|   - minEndNtUI:
|     o how many bases in a DI event must be to be a DI
|       event
|   - numKmersSIPtr:
|     o pointer to signed int to hold the number of kmers
|       shared with the reference
|   - alnSetSTPtr:
|     o pointer to alnSet struct with alignment settings
|   - matrixSTPtr:
|     o pointer to dirMatrix struct to use in alignment
| Output:
|   - Modifies:
|     o samSTPtr to hold the new alignment
|     o matrixSTPtr to have the direction matrix, score,
|       and index from alignment
|     o kmerArySI to have the sequences kmers list
|     o cntArySI to have the sequences kmer counts
|     o segSIPtr to have the segment mapped to or -1
|   - Returns:
|     o 0 if there were no DI events
|     o > 0 if found DI events (number of events returned)
|     o def_noMatch_diScan (-1) if not a DI sequence
|     o def_memErr_diScan (-2) if memory error
\-------------------------------------------------------*/
signed int
waterScan_diScan(
   struct seqST *seqSTPtr,/*sequence to check*/
   struct kmerCnt *refAryST,  /*array of references*/
   unsigned int lenRefUI,     /*# of references*/
   signed int *kmerArySI,     /*holds sequence kmers*/
   signed int *cntArySI,      /*holds kmer counts*/
   float minPercScoreF,       /*min % score to check DIs*/
   float minKmerPercF,         /*min perc kmers to keep*/
   unsigned char lenKmerUC,   /*length of one kmer*/
   unsigned int minDIDelUI,   /*min del size in DI*/
   unsigned int minPadNtUI,   /*min start/end length*/
   struct samEntry *samSTPtr, /*holds alignment*/
   signed int *segSIPtr,      /*segment mapped to*/
   signed int *numKmersSIPtr, /*number kmers shared*/
   struct alnSet *alnSetSTPtr,/*alignment settings*/
   struct dirMatrix *matrixSTPtr /*matrix for alignment*/
);

/*-------------------------------------------------------\
| Fun03: phead_diScan
|   - print out header for fragment tsv
| Input:
|   - outFILE:
|     o file to print header to
|   - lenKmerUC:
|     o length of one kmer
| Output:
|   - Prints:
|     o header to outFILE
\-------------------------------------------------------*/
void
phead_diScan(
   void * outFILE,
   unsigned char lenKmerUC
);

/*-------------------------------------------------------\
| Fun04: pfrag_diScan
|   - print out the fragment
| Input:
|   - samSTPtr:
|     o pointer to samEntry structure with sequence to
|       print
|   - numDISI:
|     o number of DI events detected
|   - segLenSI:
|     o the length of the mapped segment
|   - scoreSL:
|     o score for waterman alignment
|   - numKmersSI:
|     o number of kmers shared between ref and read
|   - outFILE:
|     o file to print read to
| Output:
|   - Prints:
|     o read id, mapped segment, classifaction, and other
|       information to a tsv file
\-------------------------------------------------------*/
void
pfrag_diScan(
   struct samEntry *samSTPtr,
   signed int numDISI,
   signed int segLenSI,
   signed long scoreSL,
   signed int numKmersSI,
   void * outFILE
);

#endif
