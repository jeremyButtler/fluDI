/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' seqST SOF: Start Of File
'   - holds functions for reading in or manipulating
'     sequences
'   o header:
'     - guards and defined variables
'   o .h st01: seqST
'     - holds an single sequence (fasta/fastq)
'   o .c fun01 addLine_seqST:
'     - add characters from file to buffer, if needed 
'       resize. This will only read in till the end of the
'       line
'   o fun02 getFqSeq_seqST:
'     - reads a fastq sequence from a fastq file
'   o fun03 getFaSeq_seqST:
'     - grabs the next read in the fasta file
'   o fun04 revComp_seqST:
'     - reverse complement a sequence
'   o fun05 blank_seqST:
'     - sets values in seqST to zero
'   o fun06 init_seqST:
'     - sets values in seqST to blank values
'   o fun07 freeStack_seqST:
'     - frees variables in an seqST (calls blank_seqST)
'   o fun08 freeHeap_seqST:
'     - frees an seqST strucuter (calls fredSeqSTStack)
'   o fun09: freeHeapAry_seqST
'     - frees an array of seqST's
'   o fun10 cpIdEndPad_seqST:
'      - copies read id to a buffer and adds in endIdC to
'        the end. If needed, this function will add right
'        padding of spaces to the end.
'   o fun11: cp_seqST
'     - copies an seqST structure
'   o fun12: swap_seqST
'     - swaps values in two seqST structs
'   o fun13: sort_seqST
'     - sorts an array of seqST structs
'   o fun14: search_seqST
'     - searchs for sequence in seqST struct array
'   o fun15: mkAry_seqST
'     - makes a seqST struct array
'   o fun16: realloc_seqST
'     - add more memory to a seqST struct array
'   o fun17: readFaFile_seqST
'     - get all sequences from a fasta file
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - guards and defined variables
\-------------------------------------------------------*/

#ifndef SEQUENCE_STRUCTURE_H
#define SEQUENCE_STRUCTURE_H

#define def_EOF_seqST 1
#define def_fileErr_seqST 2
#define def_badLine_seqST 16 /*invald fastx entry*/
#define def_memErr_seqST 4

/*--------------------------------------------------------\
| ST01: seqST
|  - holds an single sequence (fasta/fastq)
\--------------------------------------------------------*/
typedef struct seqST
{
  signed char *idStr;         /*sequence id*/
  unsigned long  lenIdUL;     /*length of sequence id*/
  unsigned long  lenIdBuffUL; /*length of Id buffer*/

  signed char *seqStr;        /*sequence*/
  unsigned long  lenSeqUL;    /*length of the sequence*/
  unsigned long  lenSeqBuffUL;/*length sequence buffer*/

  signed char *qStr;          /*q-score entry*/
  unsigned long  lenQUL;      /*length of the q-score*/
  unsigned long  lenQBuffUL;  /*length of q-score buffer*/

  unsigned long  offsetUL;    /*offset for an alignment*/
  unsigned long  endAlnUL;    /*marks end of alignment*/
}seqST;

/*-------------------------------------------------------\
| Fun02: readRefFqSeq
|  - grabs the next read in the fastq file
| Input:
|  - fqFILE:
|    o pointer to FILE handle to a fastq file to get the
|      next read from
|  - seqST:
|    o pointer to seqST structure to hold next read
| Output:
|  - Modifies:
|    o refStruct to hold the read in fastq entry & sets
|      its pointers
|  - Returns:
|     o 0: for no errors
|     o def_EOF_seqST: if EOF
|     o def_fileErr_seqST If file was not a fastq file
|     o def_badLine_seqST | def_fileErr_seqST for  a
|       invalid fastq entry
|     o def_memErr_seqST If malloc failed to find memory
\-------------------------------------------------------*/
signed char
getFqSeq_seqST(
  void *fqFILE,           /*fastq file with sequence*/
  struct seqST *seqSTPtr /*will hold one fastq entry*/
);

/*-------------------------------------------------------\
| Fun03 TOC: getFaSeq_seqST
|  -  grabs the next read in the fasta file
| Input:
|  - faFILE:
|    o Pointer to FILE handle to a fasta file to get the
|      next read from
|  - seqST:
|    o pointer to seqST structure to hold next read
| Output:
|  - Modifies:
|    o seqST to hold one fasta entry
|  - Returns:
|     o 0: for no errors
|     o def_EOF_seqST: if EOF
|     o def_fileErr_seqST If file was not a fastq file
|     o def_badLine_seqST | def_fileErr_seqST for  a
|       invalid fasta entry
|     o def_memErr_seqST If malloc failed to find memory
| Note:
|   - this will remove new lines from the sequence.
|   - this will only remove spaces or tabs at the end of
|     each sequence entry, so "atg atc \n" will got to
|     "atc atc".
\-------------------------------------------------------*/
signed char
getFaSeq_seqST(
  void *faFILE,           /*fasta file with sequence*/
  struct seqST *seqSTPtr  /*will hold one fastq entry*/
);

/*-------------------------------------------------------\
| Fun04: revComp_seqST
|  - reverse complement a sequence
| Input:
|  - seqST:
|    o pointer to seqST structure with a sequence to
|      reverse complement
| Output:
|  - Modfies
|    o reverse complements the sequence entry and if
|      reverses the Q-score entry
\-------------------------------------------------------*/
void
revComp_seqST(
  struct seqST *seqSTPtr  /*sequence to reverse comp*/
);


/*-------------------------------------------------------\
| Fun05: blank_seqST
|  - sets values in seqST to blank values
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer to blank
| Output:
|  - Modifies
|    o sets id, sequence, and Q-score entreis to start
|      with '\0' and the id, sequence, and Q-score lengths
|      to 0. This does not change the buffer lengths.
\-------------------------------------------------------*/
void
blank_seqST(
   struct seqST *seqSTPtr
);

/*-------------------------------------------------------\
| Fun06: init_seqST
|  - sets vlues in seqST to zero
| Input:
|  - seqSTPtr:
|    o Pointer to an seqST structuer to initialize
| Output:
|  - Modifies
|    o all values in seqST to be 0
\-------------------------------------------------------*/
void
init_seqST(
   struct seqST *seqSTPtr
);

/*-------------------------------------------------------\
| Fun07: freeStack_seqST
|  - frees the variables in an seqST strucuter
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer that has variables
|      to free
| Output:
|  - Frees
|    o seqSTPtr-idStr and sets to 0
|    o seqSTPtr-seqStr and sets to 0
|    o seqSTPtr-qStr and sets to 0
\-------------------------------------------------------*/
void
freeStack_seqST(
   struct seqST *seqSTPtr
);

/*-------------------------------------------------------\
| Fun08: freeHeap_seqST
|  - frees the seqST strucuter
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer to free
| Output:
|  - Frees
|    o seqST from memory and sets to 0
\-------------------------------------------------------*/
void
freeHeap_seqST(
   struct seqST *seqSTPtr
);

/*-------------------------------------------------------\
| Fun09: freeHeapAry_seqST
|  - frees an array of seqST's
| Input:
|  - seqAryST:
|    o pointer to an arrya of seqST structuer to free
|  - numSeqSI:
|    o number of seqSTs in the seqST array
|     (seqSTAryPtr)
| Output:
|  - Frees
|    o seqST from memory
\-------------------------------------------------------*/
void
freeHeapAry_seqST(
   struct seqST *seqAryST,
   signed long  numSeqSL
);

/*-------------------------------------------------------\
| Fun10: cpIdEndPad_seqST
|  - copies read id to a buffer and adds in endIdC to
|    the end. If needed, this function will add right
|    padding of spaces to the end.
| Input:
|  - seqST:
|    o pointer to seqST with read id to copy
|  - buffStr:
|    o c-string to copy read id to
|  - endIdC:
|    o character to add in to mark the end of the read id
|  - padRI:
|    o amount of padding to add to the right of the read
|      id
| Output:
|  - Modifies:
|    o buffStr to hold the sequence id + endIdC. If the
|      id + endIdC is shorter than padRI, the copied id is
|      padded with spaces on the right till it is the same
|      size as padRI.
|  - Returns:
|    o pointer to end of copied id or padding if padding
|      was applied.
\-------------------------------------------------------*/
signed char *
cpIdEndPad_seqST(
   struct seqST *seqSTPtr,  /*has read id to copy*/
   signed char *buffStr, /*buffer to add read id to*/
   signed char endIdC,   /*padding; id end (0 to skip)*/
   signed int padRI      /*padding; id start*/
);

/*-------------------------------------------------------\
| Fun11: cp_seqST
|   - Copies an seqST structure
| Input:
|   - dupSeqST:
|     o Pointer to an seqST to copy to (duplicate)
|   - cp_seqST:
|     o Pointer to an seqST to copy
| Output:
|   - Modifies:
|     o dupSeqST to be an copy of cp_seqST
|   - returns:
|     o 0 for no errors
|     o def_memErr_seqST for an memory error
\-------------------------------------------------------*/
signed char
cp_seqST(
   struct seqST *dupSeqST, /*copy to (duplicate)*/
   struct seqST *cpSeqST   /*seqSTPtr to copy*/
);

/*-------------------------------------------------------\
| Fun12: swap_seqST
|   - swaps values in two seqST structs
| Input:
|   - firstSTPtr:
|     o pointer first seqST struct to swap
|   - secSTPtr:
|     o pointer second seqST struct to swap
| Output:
|   - Modifies:
|     o firstSTPtr to have values from secSTPtr
|     o secSTPtr to have values from firstSTPtr
\-------------------------------------------------------*/
void
swap_seqST(
   struct seqST *firstSTPtr,
   struct seqST *secSTPtr
);

/*-------------------------------------------------------\
| Fun13: sort_seqST
|   - sorts an array of seqST structs
| Input:
|   - arySTPtr:
|     o pointer to array of seqST struct to sort
|   - lenSL:
|     o length to sort in array
| Output:
|   - Modifies:
|     o arySTPtr to be sorted by read/ref id
\-------------------------------------------------------*/
void
sort_seqST(
   struct seqST *arySTPtr, /*array to sort*/
   signed long lenSL       /*length of array*/
);

/*-------------------------------------------------------\
| Fun14: search_seqST
|  - searchs for sequence in seqST struct array
| Input:
|  - seqArySTPtr:
|    o pionter to seqST stuct array to search
|  - qryStr:
|    o query to find
|  - lenSL:
|    o length of seqArySTPtr (index 1)
| Output:
|  - Returns:
|    o index of qryStr in seqArySTPtr
|    o -1 if qryUL is not in ulAry
\-------------------------------------------------------*/
signed long
search_seqST(
   struct seqST *seqArySTPtr, /*sequence array to search*/
   signed char *qryStr,       /*query to hunt for*/
   signed long lenSL          /*length to search in*/
);

/*-------------------------------------------------------\
| Fun15: mkAry_seqST
|  - makes a seqST struct array
| Input:
|  - sizeSL:
|    o size of seqST struct array
| Output:
|  - Returns:
|    o pointer to seqST struct array
|    o 0 for memory errors
\-------------------------------------------------------*/
struct seqST *
mkAry_seqST(
   signed long sizeSL         /*new length (size)*/
);

/*-------------------------------------------------------\
| Fun16: realloc_seqST
|  - add more memory to a seqST struct array
| Input:
|  - seqArySTPtr:
|    o pionter to seqST stuct array to expand
|  - lenSL:
|    o index of last item
|    o this marks the start of initialization
|  - sizeSL:
|    o new size of seqST struct array
| Output:
|  - Returns:
|    o 0 for no errors
|    o def_memErr_seqST for memory errors
\-------------------------------------------------------*/
signed char
realloc_seqST(
   struct seqST **seqArySTPtr,/*sequence array to expand*/
   signed long lenSL,         /*index of last item*/
   signed long sizeSL         /*new length (size)*/
);

/*-------------------------------------------------------\
| Fun17: readFaFile_seqST
|   - get all sequences from a fasta file
| Input:
|   - faFILE:
|     o FILE pointer to file to get sequences from
|   - lenSLPtr:
|     o pointer to signed long to hold number sequences
|       in returned array (also in faFILE)
|   - sizeSLPtr:
|     o pointer to signed long to hold returned array size
|   - errSCPtr:
|     o pointer to signed char to hold errors
| Output:
|   - Modifies:
|     o faFILE to be at end of file
|     o lenSLPtr to have number sequences in faFILE/array
|     o sizeSLPtr to have size of returned array
|     o errSCPtr to have errors
|       * 0 for no errors
|       * def_fileErr_seqST if had file error
|       * def_badLine_seqST | def_fileErr_seqST for
|         invalid fasta entry
|       * def_memErr_seqST for memory errors
|   - Returns:
|     o seqST struct array with sequences
|     o 0 for errors
\-------------------------------------------------------*/
struct seqST *
readFaFile_seqST(
   void *faFILE,           /*file to get reads from*/
   signed long *lenSLPtr,  /*number items in array*/
   signed long *sizeSLPtr, /*length of array*/
   signed char *errSCPtr   /*holds error values*/
);

#endif
