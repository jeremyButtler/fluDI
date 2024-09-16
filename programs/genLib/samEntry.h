/*########################################################
# Name: samEntryStruct
#  - Holds structer to hold a sam file entry. This also
#    includes the functions needed to support this
#    structure.
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - header guards and definitions
'   o .h st01 samEntry:
'     - Holds a single samfile entry
'   o fun01 blank_samEntry:
'     - Sets all non-alloacted variables in samEntryST to 0
'   o fun02 init_samEntry:
'     - Initalize a samEntry struct to 0's
'   o fun03: setup_samEntry
'     - allocates memory for a samEntry structure (call
'       after init_samEntry (fun02))
'   o fun04 freeStack_samEntry:
'     - Frees heap allocations in a stack allocated
'       samEntry struct
'   o fun05 freeHeap_samEntry:
'     - Frees a samEntry structer (and sets to null)
'   o fun06: makeSamEntry
'     - Makes an heap allocated samEntry structure
'   o fun07: qhistToMed_samEntry
'     - Gets the median q-score for an histogram of
'        q-scores in a samStruct
'   o fun08: findQScores_samEntry
'      - Gets the median and mean q-scores from a samEntry
'        Structure.
'   o fun09: cpQEntry_samEntry
'     - Copies q-scores from a string into a samEntry
'       structure
'   o fun10: get_samEntry
'     - Reads in a single line from a sam file
'   o fun11: findRefPos_samEntry
'     - Find an reference coordinate in an sequence in
'       an sam entry structure
'   o fun12: p_samEntry
'     - Prints the sam file entry to a file. This does not
'       print any extra stats that were found.
'   o fun13: p_samEntryAsFastq
'     - Prints the sam entry as a fastq entry to a fastq
'       file
'   o fun14: pfa_samEntry
'     - Prints the sam entry as a fasta entry to a fasta
'       file
'   o fun15: pstats_samEntry
'     - Prints out the stats in a samEntry struct to a file
'   o fun16: revCmp_samEntry
'      - reverse complements a sam file sequence entry
'   o .h note01:
'      - Notes about the sam file format from the sam file
'        pdf
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Header gaurds and definitions
\-------------------------------------------------------*/

#ifndef SAMENTRYSTRUCT_H
#define SAMENTRYSTRUCT_H

#define def_adjQ_samEntry 33 /*offest to get q-score of 0*/
#define def_maxQ_samEntry 94 /*highest possible q-score*/

#define def_EOF_samEntry 1
#define def_memErr_samEntry 2
#define def_fileErr_samEntry 4

/*-------------------------------------------------------\
| ST01: samEntry
|  - Holds a single samfile entry
\-------------------------------------------------------*/
typedef struct samEntry
{ /*samEntry*/
    /*Buffers for storing strings*/
    signed char qryIdStr[128]; /*Holds query id/name*/
    unsigned char lenQryIdUC;  /*Length of query id/name*/

    signed char refIdStr[128]; /*Holds reference id/name*/
    unsigned char lenRefIdUC;  /*Length of ref id/name*/

    signed char *cigTypeStr;   /*Holds cigar type entry*/
    signed int *cigArySI;      /*Holds the cigar number*/
    unsigned int lenCigUI;     /*Length of cigar entry*/
    unsigned int lenCigBuffUI;/*# bytes malloc to cigStr*/

    signed char rNextStr[128]; /*Holds rNext entry*/
    unsigned char lenRNextUC;  /*Length of rNextStr*/

    /*The sequence and q-scores entries have lengths
    ` stored in readLenUI
    */
    signed char *seqStr;       /*Holds sequence entry*/
    unsigned int lenSeqBuffUI;/*bytes malloc seqStr/qStr*/

    signed char *qStr;         /*Holds the q-score entry*/
    unsigned int lenQBuffUI;   /*number q-score entries*/

    signed char *extraStr;     /*Extra (columns 12-end)*/
    unsigned int lenExtraUI;   /*Length of extraStr*/
    unsigned int lenExtraBuffUI;/*# bytes malloc to qStr*/

    /*Flags/single numeric values in the sam entry*/
    unsigned char mapqUC;  /*Holds mapping quality*/
    unsigned short flagUS; /*Holds the flag*/
    signed int pNextSI;  /*PNEXT mate/next read position*/
    signed int tLenSI;   /*TLEN observed template length*/

    /*Stats*/

    /*The q-scores; includes the softmasked regions*/
    float medianQF;    /*holds median read q-score*/
    float meanQF;      /*holds mean read q-score*/
 
    unsigned int refStartUI;  /*First reference base*/
    unsigned int refEndUI;    /*Last reference base*/
    unsigned int readLenUI;   /*Holds read length*/
    unsigned int alnReadLenUI;/*Number ref bases aligned*/

    unsigned int numMatchUI;/*Holds number of matches*/
    unsigned int numSnpUI;  /*Holds number of mismatches*/
    unsigned int numInsUI;  /*Holds number of insertions*/
    unsigned int numDelUI;  /*number of deletions*/
    unsigned int numMaskUI; /*number soft masked bases*/

    /*These variables are used in finding the q-scores*/
    unsigned int qHistUI[def_maxQ_samEntry + 8];
    unsigned long sumQUL;             /*Total for mean Q*/
}samEntry;

/*-------------------------------------------------------\
| Fun01: blank_samEntry
|   - Sets all values to 0, or for c-strings to '\0'
| Input:
|   - samSTPtr:
|     o Pointer to samEntry structure to blank
| Output:
|   - Modifies:
|     o Sets every variable to 0
|     o the first character in every c-string to be '\0'
\-------------------------------------------------------*/
void
blank_samEntry(
   struct samEntry *samSTPtr
);

/*-------------------------------------------------------\
| Fun02: init_samEntry
|   - Initializes a samEntry structure for use. This 
|     function should only ever be called once per
|     structure or after freeStack_samEntry has been used.
|     Use blank_samEntry to reset the values in a samEntry
|     struct after initialization.
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to initialize
| Output:
|   - Modifies:
|     o all variables in samEntry to be 0 or null
|   - Returns:
|     o 0 for success
|     o def_memErr_samEntry for memory error
\-------------------------------------------------------*/
void
init_samEntry(
   struct samEntry *samSTPtr
);

/*-------------------------------------------------------\
| Fun03: setup_samEntry
|   - allocates memory for a samEntry structure (call
|     after init_samEntry (fun02))
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to allocate memory to
| Output:
|   - Allocates:
|     o memory for seqStr, qStr, cigTypStr, cigValAryI,
|       and extraStr
|   - Returns:
|     o 0 for success
|     o def_memErr_samEntry for memory error
\-------------------------------------------------------*/
signed char
setup_samEntry(
   struct samEntry *samSTPtr
);

/*-------------------------------------------------------\
| Fun04: freeStack_samEntry
| Use:
|   - Frees all variables in samEntry, but not samEntry
| Input:
|   - samSTPtr
|     o Pointer to samEntry struct to free the interal
|       memory of
| Output:
|   - Frees:
|     o allocated memory in samSTPtr
\-------------------------------------------------------*/
void
freeStack_samEntry(
   struct samEntry *samSTPtr
);

/*-------------------------------------------------------\
| Fun05: freeHeap_samEntry
|   - Frees a samEntry struct
| Input:
|   - samSTPtr
|     o Pointer to Sam entry struct to free
| Output:
|   - Frees:
|     o samSTPtr and its variables with allocated memory
\-------------------------------------------------------*/
void
freeHeap_samEntry(
   struct samEntry **samSTPtr
);

/*-------------------------------------------------------\
| Fun06: makeSamEntry
|   - Makes an heap allocated samEntry structure
| Input:
| Output:
|   - Returns:
|     o Pointer to a samEntry structure
|     o 0 if had an memory error
\-------------------------------------------------------*/
struct samEntry *
makeSamEntry();

/*-------------------------------------------------------\
| Fun07: qhistToMed_samEntry
|   - Gets the median q-score for an histogram of q-scores
|     in a samStruct
| Input:
|   - samSTPTr:
|     o Pointer to samEntry struct with a q-score
|       histogram to find the median q-score for
| Output:
|   - Modifies:
|     o samSTPtr->medianQF to have the median q-score
\-------------------------------------------------------*/
void
qhistToMed_samEntry(
   struct samEntry *samSTPtr 
);

/*-------------------------------------------------------\
| Fun08: findQScores_samEntry
|   - Gets the median and mean q-scores from a samEntry
|     Structure.
| Input:
|   - samSTPTr:
|     o Pointer to samEntry struct to find the median and
|       mean q-scores for
| Output:
|   - Modifies:
|     o samSTPtr->medianQF to have the median q-score
|     o samSTPtr->meanQF to have the mean q-score
|     o samSTPtr->qHistUI to have a histogram of q-scores
|     o samSTPtr->sumQUL to have the sum of all q-scores
\-------------------------------------------------------*/
void
findQScores_samEntry(
   struct samEntry *samSTPtr
);

/*-------------------------------------------------------\
| Fun09: cpQEntry_samEntry
|   - Copies q-scores from a string into a samEntry
|     structure
| Input:
|   - samSTPtr:
|     o Pionter to sam entry struct to copys q-scores to
|   _ cpQStr:
|     o C-string with q-scores to copy to samSTPtr
|   - blankQHistBl:
|     o 1: Blank q-score vars (histogram/sum/mean/median)
|     o 0: do not blank the q-score variables
| Output:
|   - Mofidies:
|     o qStr in samSTPtry to have the q-scores
|     o medianQF in samSTPtr to have the median q-score
|     o meanQF in samSTPtr to have the mean q-score
|     o qHistUI in samSTPtr to have histogram of q-scores
|     o samQUL in samSTPtr to have sum off all q-scores
|   - Returns
|     o The value in samSTPtr->readLenUI
\-------------------------------------------------------*/
unsigned int
cpQEntry_samEntry(
   struct samEntry *samSTPtr, /*Copy q-scores to*/
   signed char *cpQStr,       /*q-scores to copy*/
   signed char blankQHistBl   /*1: to blank q-score hist*/
);

/*-------------------------------------------------------\
| Fun10: get_samEntry
|  - Reads in a single line from a sam file
| Input:
|  - samSTPtr:
|    o Pionter to samEntry structure to store the sam file
|      line in. This structure should be initialized
|  - buffStr:
|    o Buffer to read the sam file line temporarly into.
|      This is resized if needed. You can input NULL to
|      create a new buffer.
|  - lenBuffUL:
|    o Length of buffStr (updated if buffStr is resized)
|  - samVoidFILE:
|    o Sam file to read a line from. The void is so that
|      I can use samFILE in the function.
| Output:
|  - Modifies:
|    o samSTPtr to have the next line
|    o samFILE to be on the next line
|    o buffStr to hold a sam file line (resized if needed)
|    o lenBuffUL to hold the resized length of buffStr
|  - Returns:
|    o 0 for success
|    o def_EOF_samEntry for EOF (End Of File)
|    o def_memErr_samEntry for memory errors
\-------------------------------------------------------*/
signed char
get_samEntry(
   struct samEntry *samSTPtr,
   signed char **buffStr,
   unsigned long *lenBuffUL,
   void *samVoidFILE
);

/*-------------------------------------------------------\
| Fun11: findRefPos_samEntry
|   - Find an reference coordinate in an sequence in
|     an sam entry structure
| Input:
|   - samSTPtr:
|     o Pointer to samEntry structure with sequence
|       (query) to find target reference position in
|   - siCig:
|     o pionter to signed int with index of the cigar
|       entry I am on
|   - cigNtSI:
|     o pionter to a singed int with the number of bases
|       left in the cigar entry
|   - targPosSI:
|     o Reference position I wish to find
|   - refPosSI:
|     o pointer to signed int with the current position
|       in reference sequence
|   - seqPosSI:
|     o pointer to signed int with the current position
|       in query sequence
| Output:
|   - Modifies:
|     o siCig to point to the next open cigar entry
|       - will be > samSTPtr->lenCigUI when the sequence
|         does not end at at targPosSI
|     o cigNtSI to have the number of bases remianing
|       in the current siCig entry
|     o refPosSI to be the index to the position in targSI
|     o seqPosSI to be the index of the sequence base at
|       refPosSI (check cigar to see if deletion)
|   - Returns:
|     o Index of last cigar entry I incurmented
\-------------------------------------------------------*/
signed int
findRefPos_samEntry(
   struct samEntry *samSTPtr,/*has cigar and sequence*/
   signed int *siCig,        /*current cigar entry on*/
   signed int *cigNtSI,     /*# of bases left in cigar*/
   signed int targPosSI,     /*position looking for*/
   signed int *refPosSI,     /*current ref position*/
   signed int *seqPosSI      /*current sequence position*/
);

/*-------------------------------------------------------\
| Fun12: p_samEntry
|   - Prints the sam file entry to a file. This does not
|     print any extra stats that were found.
| Input:
|   - samST
|     o Pointer to samEntry struct with sam entry to print
|   - buffStr:
|     o Pointer to c-string buffer to temporarly hold the
|       cigar entry (speeds things up)
|   - lenBuffUL:
|     o Current length of buffer, adjusted if buffStr is
|       expanded
|   - pNoNewLineBl:
|     o 1: do not print a new line after; you will do this
|     o 0: end of sam entry, print a new line
|   - outFILE:
|     o File to print the sam entry to
| Output:
|   - Prints:
|     o Sam file entry in samST to outFILE.
|   - Modifies:
|     o inceases the size of buffStr, if is not 8x the
|       cigar length
|     o Sets lenBuffUL to the new buffStr size when
|       buffStr is resized
|   - Returns:
|     o 0 for no problems
|     o def_memErr_samEntry for memory errors
\-------------------------------------------------------*/
signed char
p_samEntry(
   struct samEntry *samSTPtr,
   signed char **buffStr,
   unsigned long *lenBuffUL,
   signed char pNoNewLineBl,
   void *outFILE
);

/*-------------------------------------------------------\
| Fun13: pfq_samEntry
|   - Prints the sam entry as fastq entry to a fastq file
| Input:
|   - samST:
|     o Pointer to samEntry structure with fastq entry to
|       print out
|   - outFILE:
|     o Fastq file to print the new fastq entry to
| Output:
|   - Prints:
|     o fastq entry from samST to outFILE
\-------------------------------------------------------*/
void
pfq_samEntry(
   struct samEntry *samSTPtr,
   void *outFILE
);

/*-------------------------------------------------------\
| Fun14: pfa_samEntry
| Use:
|   - Prints the sam entry as fasta entry to a fasta file
| Input:
|   - samST:
|     o Pointer to samEntry structure with fastq entry to
|       print out
|   - outFILE:
|     o Fasta file to print the new fasta entry to
| Output:
|   - Prints:
|     o fasta entry from samST to outFILE
\-------------------------------------------------------*/
void
pfa_samEntry(
   struct samEntry *samSTPtr,
   void *outFILE
);

/*-------------------------------------------------------\
| Fun15: pstats_samEntry
| Use:
|  - Prints out the stats in a samEntry struct to a file
| Input:
|  - samEntryST:
|    o Pointer to samEntry struct to print stats for
|  - pHeadBl:
|    o 1: Print the header for the stats tsv file
|    o 0: Do not print the header
|  - pNsBl:
|    o 1: find and print out the anonymous base counts
|    o 0: do not print out anonymous base counts
|  - outFILE:
|    o TSV (tab deliminated) file to print stats to
| Output:
|  - Prints:
|    o line with stats from samEntryStruct to file (tsv)
|  - Modifies:
|    o printHeaderChar to be 0 if set to 1
\-------------------------------------------------------*/
void
pstats_samEntry(
   struct samEntry *samSTPtr,
   signed char *pHeadBl,
   signed char pNsBl,
   void *outFILE
);

/*-------------------------------------------------------\
| Fun16: revCmp_samEntry
|   - reverse complements a sam file sequence entry
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct to reverse complement
| Output:
|   - Modifies:
|     o samSTPtr (includes flag) to be reverse complment
\-------------------------------------------------------*/
void
revCmp_samEntry(
   struct samEntry *samSTPtr
);

#endif

/*-------------------------------------------------------\
| Note01:
|   - Notes about the sam file format from the sam file
|     pdf
\-------------------------------------------------------*/

/*
Sam file table for first 11 columns (all sam files have)
| Col | Field |  Type  |        Brief description              |
|:---:|:-----:|:------:|:-------------------------------------:|
|  1  | QNAME | String |       Query template NAME             |
|  2  | FLAG  |  Int   |          bitwise FLAG                 |
|  3  | RNAME | String |     Reference sequence NAME           |
|  4  |  POS  |  Int   |  1-based leftmost mapping POSition    |
|  5  | MAPQ  |  Int   |          MAPping Quality              |
|  6  | CIGAR | String |            CIGAR string               |
|  7  | RNEXT | String | Reference name of the mate/next read  |
|  8  | PNEXT |  Int   |   Position of the mate/next read      |
|  9  | TLEN  |  Int   |      observed Template LENgth         |
| 10  |  SEQ  | String |          segment SEQuence             |
| 11  | QUAL  | String | ASCII of Phred-scaled base Quality+33 |
*/

/*eqx cigar entry
`   matches are '='
`   mimsatches are 'X'
'   non-eqx merges matches and mismatches into 'M'
`   deletions are 'D'
`   insertions are 'I'
`   soft masks are 'S'
`   Numbers come before entry (=, X, D, I, or S) & show
`     number of times the entry is repeated
`   Everything else is a hard mask (was removed) or
'      removed the bases from the sequence
`   EX: 10S5=1X701 (10 soft masked bases, 5 matches,
`       1 SNP, 701 matches)
*/

/* Sam file flag values tables
| Bit  | FLAG  |                        Description                                 |
|:----:|:-----:|:------------------------------------------------------------------:|
| 1    | 0x1   | template having multiple segments in sequencing                    |
| 2    | 0x2   | each segment properly aligned according to the aligner             |
| 4    | 0x4   | segment unmapped                                                   |
| 8    | 0x8   | next segment in the template unmapped                              |
| 16   | 0x10  | SEQ being reverse complemented                                     |
| 32   | 0x20  | SEQ of the next segment in the template being reverse complemented |
| 64   | 0x40  | the first segment in the template                                  |
| 128  | 0x80  | the last segment in the template                                   |
| 256  | 0x100 | secondary alignment                                                |
| 512  | 0x200 | not passing filters, such as platform/vendor quality controls      |
| 1024 | 0x400 | PCR or optical duplicate                                           |
| 2048 | 0x800 | supplementary alignment                                            |
*/

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
