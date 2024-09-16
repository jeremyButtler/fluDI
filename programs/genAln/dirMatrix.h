/*########################################################
# Name: dirMatrix
#  - holds functions for dealing with the dirMatrix
#    returned by water and needle
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - forward declerations and guards
'   o .h st01: alnMatrixStruct
'     - Holds the direction matrix and best score(s) for a
'       single aligment
'   o fun01: blank_dirMatrix
'     - blanks all score, index, and length variables in a
'       dirMatrix structure
'   o fun02: init_dirMatrix
'     - initializes a dirMatrix structure
'   o fun03: freeStack_dirMatrix
'     - frees heap allocated variables in a dirMatrix
'       structure
'   o fun04: freeHeap_dirMatrix
'     - frees a dirMatrix structure
'   o fun05: getAln_dirMatrix
'     - gets a sam file entry (alignment) from a direction
'       matrix (inside the dirMatrix structure)
'   o fun06: getCig_dirMatrix
'     - gets a cigar for an alignment
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - forward declerations and guards
\-------------------------------------------------------*/

#ifndef DIRECTION_MATRIX_H
#define DIRECTION_MATRIX_H

typedef struct samEntry samEntry;
typedef struct seqST seqST;
typedef struct alnSet alnSet;

/*-------------------------------------------------------\
| St01: dirMatrix
|  - Holds the direction matrix and best score(s) for a
|    single aligment
\-------------------------------------------------------*/
typedef struct dirMatrix
{ /*alnStruct*/
  signed char *dirMatrixSC; /*directional matrix*/
  unsigned long lenMatrixUL;/*size of directional matrix*/

  signed long scoreSL;       /*best score in alignment*/
  unsigned long indexUL;     /*index for best score*/

  signed long *scoreArySL;   /*score array used to score*/
  unsigned long lenScoreUL;  /*length of score array*/

  unsigned long lenRefUL;    /*reference length*/
  unsigned long refOffsetUL; /*first ref base to align*/
  unsigned long refEndUL;    /*last ref base to align*/

  unsigned long lenQryUL;    /*length of query*/
  unsigned long qryOffsetUL; /*first query base to align*/
  unsigned long qryEndUL;    /*last query base to align*/

  signed char errSC;         /*holds any errors*/
}dirMatrix;

/*-------------------------------------------------------\
| Fun01: blank_dirMatrix
|   - blanks all score, index, and length variables in a
|     dirMatrix structure
| Input:
|   - dirMatrixSTPtr:
|     o a pointer to a dirMatrix to blank
| Output:
|   - Sets:
|     o all score, length, and index variables are set to
|       0 (direction matrix is not touched)
\-------------------------------------------------------*/
void
blank_dirMatrix(
  struct dirMatrix *dirMatrixSCSTPtr
);

/*-------------------------------------------------------\
| Fun02: init_dirMatrix
|   - initializes a dirMatrix structure
| Input:
|   - dirMatrixSTPtr:
|     o a pointer to a dirMatrix structure initialize
| Output:
|   - Sets:
|     o all variables in matrixSTPtr to 0
\-------------------------------------------------------*/
void
init_dirMatrix(
  struct dirMatrix *dirMatrixSTPtr
);

/*-------------------------------------------------------\
| Fun03: freeStack_dirMatrix
|   - frees heap allocated variables in a dirMatrix
|     structure
| Input:
|   - matrixSTPtr
|     o pointer to dirMatrix structure with variables to
|       free
| Output:
|   - Frees:
|     o dirMatrix->dirMatrixSC
|   - Sets:
|     o all non-freeded variables to 0
\-------------------------------------------------------*/
void
freeStack_dirMatrix(
   struct dirMatrix *matrixSTPtr
);

/*-------------------------------------------------------\
| Fun04: freeHeap_dirMatrix
|   - Frees a dirMatrix structure
| Input:
|   - matrixSTPtr
|     o pointer to a dirMatrix structure to free
| Output:
|   - Frees:
|     o matrixSTPtr
\-------------------------------------------------------*/
void
freeHeap_dirMatrix(
   struct dirMatrix *matrixSTPtr
);

/*-------------------------------------------------------\
| Fun05: getAln_dirMatrix
|   - gets a sam file entry (alignment) from a direction
|     matrix (inside the dirMatrix structure)
| Input:
|   - matrixSTPtr
|     o pointer to a dirMatrix structure to get alignment
|       from
|   - indexUL:
|     o index of last base in the alignment
|     o 0 to use index from matirxSTPtr
|   - revBl:
|     o 1: reverse alignment (sam flag is 16)
|     o 0: foward alignment (sam flag is 0)
|   - qrySTPtr:
|     o pointer to a seqST with the query sequence
|   - refSTPtr:
|     o pointer to a seqST with the reference sequence
|   - samSTPtr:
|     o pointer to a samEntry struct to hold the alignment
|   - numAnonUI:
|     o pointer to unsigned in to hold the number of
|       anonymous bases (matches only)
|   - alnSetSTPtr:
|     o pointer to alnSet structure with the match matrix
| Output:
|   - Modifies:
|     o samSTPtr to have the aligned sequence
|   - Returns:
|     o 0 for no errors
|     o 1 for memory error (only error possible)
\-------------------------------------------------------*/
signed char
getAln_dirMatrix(
   struct dirMatrix *matrixSTPtr,
   unsigned long indexUL,
   signed char revBl,
   struct seqST *qrySTPtr,
   struct seqST *refSTPtr,
   struct samEntry *samSTPtr,
   unsigned int *numAnonUI,
   struct alnSet *alnSetSTPtr
);

/*-------------------------------------------------------\
| Fun06: getCig_dirMatrix
|   - gets a cigar for an alignment
| Input:
|   - matrixSTPtr
|     o pointer to a dirMatrix structure to get alignment
|       from
|   - indexUL:
|     o index of last base in the alignment
|     o 0 to use index from matirxSTPtr
|   - revBl:
|     o 1: reverse alignment (sam flag is 16)
|       - this means I had to reverse complement the
|         reference sequence
|     o 0: foward alignment (sam flag is 0)
|   - qrySTPtr:
|     o pointer to a seqST with the query sequence
|   - refSTPtr:
|     o pointer to a seqST with the reference sequence
|   - cigTypeStr:
|     o pointer to c-string with cigar entry type array
|   - cigArySI:
|     o pointer to signed int array with the length of
|       each cibar entry
|   - cigPosUI:
|     o position at in cigar
|   - lenCigUI:
|     o pointer to unsigned int with length of the cigar
|       buffer
|   - refStartUI:
|     o unsigned int pointer to point to first reference
|       base in cigar
|   - numAnonUI:
|     o pointer to unsigned in to hold the number of
|       anonymous bases (matches only)
|   - numMaskUI:
|     o pointer to unsigned long to hold number of
|       maksed bases
|   - delAtEndBl:
|     o 1: add deltions if reference is short at end
|     o 0: ignore
|   - alnSetSTPtr:
|     o pointer to alnSet structure with the match matrix
| Output:
|   - Modifies:
|     o cigTypeStr to have the cigar entry types
|     o cigArySI to have the length of each cigar entry
|     o lenCigUI if cigTypeStr and cigArySI needed to be
|       resized
|     o refStartUI to have first reference base in cigar
|     o numAnonUI to have number of matching anonymous
|       bases
|     o numMaskUI to have number of masked bases
|   - Returns:
|     o new position in cigar
|     o -1 for memory error (only error possible)
\-------------------------------------------------------*/
signed long
getCig_dirMatrix(
   struct dirMatrix *matrixSTPtr,
   unsigned long indexUL,
   signed char revBl,
   struct seqST *qrySTPtr,
   struct seqST *refSTPtr,
   signed char **cigTypeStr,
   signed int **cigArySI,
   unsigned int cigPosUI,
   unsigned int *lenCigUI,
   unsigned int *refStartUI,
   unsigned int *numAnonUI,
   unsigned int *numMaskUI,
   signed char delAtEndBl,
   struct alnSet *alnSetSTPtr
);

#endif
