/*########################################################
# Name: kmerCnt
#   - holds functions to count number of matching kmers in
#     sequences
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - forward declerations, defined variables, guards
'   o fun01: blank_kmerCnt
'     - blanks all variables in a kmerCnt structure
'   o fun02: init_kmerCnt
'     - initialize a kmerCnt structure
'   o fun04: freeStack_kmerCnt
'     - frees variables inside a kmerCnt structure
'   o fun05: freeHeap_kmerCnt
'     - frees a kmerCnt structure (you set to 0)
'   o fun06: freeStackAry_kmerCnt
'     - frees internal memory form all kmerCnt structure
'       in an array
'   o fun07: freeHeapAry_kmerCnt
'     - frees an array of kmerCnt structure (you set to 0)
'   o .c fun08: mkKmerMask_kmerCnt
'     - makes a kmer mask for removing extra bases from
'       kmer
'   o fun09: addSeq_kmerCnt
'     - adds a sequence to a kmerCnt structure
'   o .c fun10: dualSort_kmerCnt
'     - sorts first unsigned int array and keeps second
'       signed int array in order with first
'   o fun11: ntToKmerAry_kmerCnt
'     - converts a nucleotide sequence to a array of kmer
'       counts
'   o fun12: get_kmerCnt
'     - gets number of matching kmers between a kmerCnt
'       structure and the kmer arrays
'   o fun13: faToKmerCnt_kmerCnt
'     - converts a the sequences in a fasta file to a
'       kmerCnt array
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - forward declerations, defined variables, guards
\-------------------------------------------------------*/

#ifndef KMER_COUNT_H
#define KMER_COUNT_H

typedef struct seqST seqST;

#define def_noKmer_kmerCnt -1
#define def_endKmers_kmerCnt -2
#define def_bitsPerKmer_kmerCnt 2 /*do not change*/

#define def_noMatch_kmerCnt 1
#define def_fileErr_kmerCnt 2
#define def_memErr_kmerCnt 4

/*-------------------------------------------------------\
| ST01: kmerCnt
|   - holds the number of kmers in each sequence
\-------------------------------------------------------*/
typedef struct kmerCnt
{
   unsigned char lenKmerUC;  /*bases in one kmer*/
   unsigned int maxKmersUI;  /*maximum number of kmers*/

   signed int *forKmerArySI; /*array of forward kmers*/
   unsigned int forKmersUI;  /*kmers in sequence*/

   signed int *revKmerArySI; /*array of reverse kmers*/
   unsigned int revKmersUI;  /*kmers in sequence*/

   /*costs a bit more to have both, but saves a little
   `  time later
   */
   struct seqST *forSeqST; /*forward sequence*/
   struct seqST *revSeqST; /*reverse seqeunce*/
}kmerCnt;

/*-------------------------------------------------------\
| Fun01: blank_kmerCnt
|   - blanks all variables in a kmerCnt structure
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCnt structure to blank
| Output:
|   - Modifies:
|     o all variables in kmerCntSTPtr to be blank settings
\-------------------------------------------------------*/
void
blank_kmerCnt(
   struct kmerCnt *kmerCntSTPtr
);

/*-------------------------------------------------------\
| Fun02: init_kmerCnt
|   - initialize a kmerCnt structure
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCnt structure to initialize
| Output:
|   - Modifies:
|     o all values in kmerCntSTPtr to be 0 or null
\-------------------------------------------------------*/
void
init_kmerCnt(
   struct kmerCnt *kmerCntSTPtr
);

/*-------------------------------------------------------\
| Fun03: setup_kmerCnt
|   - sets up kmerCnt structure for use
|     (memory allocations)
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCnt structure to initialize
|   - lenKmerUC:
|     o length of one kmer
| Output:
|   - Modifies:
|     o lenKmerUC to be one kmer length
|     o maxKmersUI to be the maximum number of kmers
|     o totalKmersUI to be 0 (by blank_kmerCnt)
|   - Allocates:
|     o maxKmersUI uints for forKmerArySI and revKMerArySI
|     o seqST for forSeqST and revSeqST
|   - Returns:
|     o 0 for no problems
|     o def_memErr_kmerCnt for memory errors
\-------------------------------------------------------*/
signed char
setup_kmerCnt(
   struct kmerCnt *kmerCntSTPtr,
   unsigned char lenKmerUC
);

/*-------------------------------------------------------\
| Fun04: freeStack_kmerCnt
|   - frees variables inside a kmerCnt structure
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCnt structure with variables to
|       free
| Output:
|   - Frees:
|     o forKmerArySI, revKmerArySI, forSeqST, and revSeqST
|   - Sets:
|     o blanks structure and sets pointers to 0
|       - you will needed to call init_kmerCnt again
|       - only lenKmerUC and maxKMerUC are not set to 0
\-------------------------------------------------------*/
void
freeStack_kmerCnt(
   struct kmerCnt *kmerCntSTPtr
);

/*-------------------------------------------------------\
| Fun05: freeHeap_kmerCnt
|   - frees a kmerCnt structure (you set to 0)
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCnt structure to free
| Output:
|   - Frees:
|     o kmerCntSTPtr
\-------------------------------------------------------*/
void
freeHeap_kmerCnt(
   struct kmerCnt *kmerCntSTPtr
);

/*-------------------------------------------------------\
| Fun06: freeStackAry_kmerCnt
|   - frees internal memory form all kmerCnt structure
|     in an array
| Input:
|   - kmerCntAryST:
|     o pointer to an array of kmerCnt structures to
|       all variables for (goes to uninitialzed state)
|   - numSTUI:
|     o number of initialized kmerCnt structures in
|       kmerCntAryST (actually have internal memory)
| Output:
|   - Frees:
|     o forKmerArySI, revKmerArySI, forSeqST, and revSeqST
|       for numSTUI kmerCnt structures in the array
|   - Sets:
|     o blanks structure and sets pointers to 0
|       - you will needed to call init_kmerCnt again
\-------------------------------------------------------*/
void
freeStackAry_kmerCnt(
   struct kmerCnt *kmerCntAryST,
   unsigned int numSTUI
);

/*-------------------------------------------------------\
| Fun07: freeHeapAry_kmerCnt
|   - frees an array of kmerCnt structure (you set to 0)
| Input:
|   - kmerCntAryST:
|     o pointer to an array of kmerCnt structures to free
|   - numSTUI:
|     o number of initialized kmerCnt structures in
|       kmerCntAryST (actually have internal memory)
| Output:
|   - Frees:
|     o kmerCntAryST
\-------------------------------------------------------*/
void
freeHeapAry_kmerCnt(
   struct kmerCnt *kmerCntAryST,
   unsigned int numSTUI
);

/*-------------------------------------------------------\
| Fun09: addSeq_kmerCnt
|   - adds a sequence to a kmerCnt structure
| Input:
|   - kmerCntSTPtr:
|     o pointer to kmerCnt structure to add sequence to
|   - seqSTPtr:
|     o pointer to a seqST with sequence to add
| Output:
|   - Modifies:
|     o kmerCntSTPtr to have the new sequence
|   - Returns:
|     o 0 for no errors
|     o def_memErr_kmerCnt for memory errors
\-------------------------------------------------------*/
signed char
addSeq_kmerCnt(
   struct kmerCnt *kmerCntSTPtr,/*will hold new sequence*/
   struct seqST *seqSTPtr   /*sequence to copy*/
);

/*-------------------------------------------------------\
| Fun11: ntToKmerAry_kmerCnt
|   - converts a nucleotide sequence to a array of kmer
|     counts
| Input:
|   - seqSTPtr:
|     o pointer to a seqST with sequence to get kmer
|       counts for
|   - lenKmerUC:
|     o length of one kmer
|   - kmerArySI:
|     o pointer to a signed int array to add kmer so
|     o needs to be size of ((lenKmerUC ^ 4) + 1)
|   - cntArySI:
|     o pointer to a signed int array to add kmer counts
|       to
|     o needs to be size of ((lenKmerUC ^ 4) + 1)
| Output:
|   - Modifies:
|     o kmerArySI to hold kmers in seqSTPtr
|     o cntArySI to hold number times each kmer happened
|     o sorts kmerArySI and cntArySI by kmer
|       - this converts the hash table to a list of kmers
|       - end will be marked with a -2
|   - Returns:
|     o number of kmers in sequence
\-------------------------------------------------------*/
signed int
ntToKmerAry_kmerCnt(
   struct seqST *seqSTPtr, /*sequence to convert*/
   unsigned char lenKmerUC,    /*length of one kmer*/
   signed int *kmerArySI,      /*will hold uniqe kmers*/
   signed int *cntArySI        /*gets # kmer duplicates*/
);

/*-------------------------------------------------------\
| Fun12: get_kmerCnt
|   - gets number of matching kmers between a kmerCnt
|     structure and the kmer arrays
| Input:
|   - kmerCntSTPtr:
|     o pointer to a kmerCntSTPtr structure with the kmer
|       table to compare to
|   - kmerArySI:
|     o pointer to a singed in array with the unique kmers
|       in the sequence (use fun11; ntToKmerAry_kmerCnt)
|   - cntArySI:
|     o pointer to a singed in array with the number of
|       duplicates per unique kmers in the sequence (use
|       fun11; ntToKmerAry_kmerCnt)
| Output:
|   - Returns:
|     o highes number of kmers found
|       - + signed long if forward sequence more kmers
|       - - signed long if reverse sequence more kmers
\-------------------------------------------------------*/
signed int
get_kmerCnt(
   struct kmerCnt *kmerCntSTPtr, /*table to get counts*/
   signed int *kmerArySI,        /*sequence uniqe kmers*/
   signed int *cntArySI          /*sequence kmer counts*/
);

/*-------------------------------------------------------\
| Fun13: faToKmerCnt_kmerCnt
|   - converts a the sequences in a fasta file to a
|     kmerCnt array
| Input:
|   - faFileStr:
|     o path to and name of fasta file with sequnces
|   - lenKmerUC:
|     o length of one kmer
|   - numSeqUI:
|     o pointer to unsigned int to hold the number of
|       sequences in the kmerCnt array
|   - errSCPtr:
|     o pointer to signed char to hold the errors
| Output:
|   - Modfies:
|     o numSeqUI to have the number of sequences in the
|       returned kmerCnt array
|     o errSCPtr to hold errors
|       - 0 for no erors
|       - def_fileErr_kmerCnt if could not open faFileStr
|       - def_memErr_kmerCnt for memory errors
\-------------------------------------------------------*/
struct kmerCnt *
faToKmerCnt_kmerCnt(
   signed char *faFileStr,/*path to fasta file*/
   unsigned char lenKmerUC, /*length of one kmer*/
   unsigned int *numSeqUI,/*will hold number of sequnces*/
   signed char *errSCPtr  /*holds error message*/
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
