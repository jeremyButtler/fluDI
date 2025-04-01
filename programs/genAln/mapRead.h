/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' mapRead SOF: Start Of File
'   - maps a read to a reference
'   o header:
'     - defined variables and guards
'   o .h st01: chains_mapRead
'     - holds kmer chains build from kmerChain_mapRead
'   o fun01: blank_chains_mapRead
'     - blanks a chains_mapRead struct
'   o fun02: init_chains_mapRead
'     - initialize a chains_mapRead struct
'   o fun03: freeStack_chains_mapRead
'     - frees all variables in a chains_mapRead struct
'   o fun04: freeHeap_chains_mapRead
'     - frees a chains_mapRead struct
'   o fun05: setup_chains_mapRead
'     - sets up memory for a chains_mapRead struct
'   o fun06: mk_chains_mapRead
'     - makes a heap allocated chains_mapRead struct
'   o fun07: realloc_chains_mapRead
'     - reallocates memroy in a chains_mapRead struct
'   o .c fun08: uiFind_mapRead
'     - finds first signed int in a range
'   o .c fun09: findChain_chains_mapRead
'     - finds kmer in range from kmerChain_mapRead return
'   o fun10: kmerChain_mapRead
'     - merge kmers into chains (longer kmers)
'   o .c fun11: numQryChains_mapRead
'     - gets the number of chains assigned to a single
'       query kmer
'   o fun12: mergeChains_mapRead
'     - merges kmer chains from fun10 (kmerChain) into
'       the chain the maximizes base count
'   o fun13: mergeToSam_mapRead
'     - convert output from mergeChains_map read to sam
'       file entry (mapq is set to 0)
'   o fun14: flipRef_mapRead
'     - flips samEntry reference to query
'   o fun15: scoreSubAln_mapRead
'     - scores the aligment in samEntry and finds best
'       sub-aligment
'   o fun16: blank_set_mapRead
'     - blank a set_mapRead struct
'   o fun17: init_set_mapRead
'     - initializes a set_mapRead struct
'   o fun18: freeStack_set_mapRead
'     - frees variables in a set_mapRead struct
'   o fun19: freeHeap_set_mapRead
'     - frees a set_mapRead struct
'   o fun20: setup_set_mapRead
'     - allocates memory for structs in a set_mapRead
'   o fun21: addQryKmers_mapRead
'     - get forward and reverse kmer sequences
'   o fun22: aln_mapRead
'     - maps a read to a reference
'   o fun23: qckAln_mapRead
'     - maps a read to a reference; returns best alignment
'     - qckAln_mapRead differs from aln_mapRead in that
'       it uses kmerScan to find best window in entire
'       reference, and then after scan, makes alignment
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header
|   - defined variables and guards
\-------------------------------------------------------*/

#ifndef MAP_READ_H
#define MAP_READ_H

struct seqST;
struct samEntry;
struct alnSet;
struct refST_kmerFind;
struct tblST_kmerFind;

/*keep these negative, merge Chains returns postive/0 on
`  success
*/
#define def_memErr_mapRead -1
#define def_noChains_mapRead -2
#define def_noScore_mapRead -4
#define def_noAln_mapRead -8

/*-------------------------------------------------------\
| ST01: chains_mapRead
|   - holds kmer chains build from kmerChain_mapRead
\-------------------------------------------------------*/
typedef struct chains_mapRead
{
   /*chain arrays (start coords, length, and score)*/
   signed int *refArySI;   /*ref kmer starting coords*/
   signed int *qryArySI;   /*query kmer start coords*/
   signed int *chainArySI; /*length of each chain*/
   signed int *scoreArySI; /*maximal bases in a chain*/
      /*has best scores foun in previous round for
      `  a chain when merged with other chain; ant method?
      */
   signed int *nextArySI;
      /*breadcrum trail used to get best score*/

   /*20bytes per chain*/

   /*limiting to ints, because if I really have
   `    over 2 billion items I am in trouble.
   `    At that point I need at least 40Gb of ram and I
   `    excedes the largest solved traveling salesman
   `    problem by many folds (< 5 million cities)
   */
   signed int lenChainSI;  /*length of chain arrays*/
   signed int sizeChainSI; /*size of chain arrays*/

   /*these allow quick jumping/indexing*/
   signed int *indexArySI;/*refAryUI ref kmers index's*/
   signed int lenIndexSI; /*length match/index arrays*/
   signed int sizeIndexSI;/*size of match/index arrays*/
   signed int lastIndexSI;/*last query nt to start chain*/
}chains_mapRead;

/*-------------------------------------------------------\
| ST02: set_mapRead
|   - settings for a aln_mapRead
\-------------------------------------------------------*/
typedef struct set_mapRead
{
   /*for kmer scaning (sort by lenArySI)*/
   signed int lenArySI[16];  /*length to swith kmer size*/
   unsigned char kmerAryUC[16]; /*kmer size per length*/
   float minKmersAryF[16];   /*% kmers needed per length*/
   signed int chainLenArySI[16]; /*minimum chain length
                                 `  to keep a chain
                                 */

   unsigned char lenKmersUC;    /*number kmer lengths*/

   float percExtraNtF;      /*% extra nt in window*/
   float percShiftF;        /*% to shift window by*/

   signed char numWinSC;    /*number windows to scan
                             `  before doing alignment
                             */
   /*for gap filling*/
   struct alnSet *alnSetST;

   /*for scoreSub_mapRead*/
   signed char subBl;     /*1: do sub-aligment score*/
   float minScoreF;       /*min % score to keep*/
   float minPercLenF;     /*min % score to keep*/

   /*this part is here to reduce memory allocation steps*/
   struct tblST_kmerFind *tblSTPtr; /*for scan/chaining*/
   unsigned char maxKmerLenUC;      /*max tbl kmer len*/
      /*keeps track of previous alignments*/
   struct refST_kmerFind *qrySTPtr; /*for scan/chaining*/
   struct chains_mapRead *chainsSTPtr; /*for chaining*/

   /*holds kmer lists needed for chaining*/
   signed int *forQryKmersArySI;/*holds query kmer list*/
   signed int *revQryKmersArySI;/*holds query kmer list*/
   signed int lenQrySI;        /*length of qryKmersArySI*/
   signed int sizeQrySI;       /*size of qryKmersArySI*/

   signed int *refKmersArySI;/*has reference kmer list*/
   signed int *refSortArySI; /*sorted ref kmer array*/
                            /*for resorting refKmerArySI*/
   signed int lenRefSI;       /*length of refKmersArySI*/
   signed int sizeRefSI;      /*size of refKmersArySI*/
}set_mapRead;

/*-------------------------------------------------------\
| Fun01: blank_chains_mapRead
|   - blanks a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct to blank
| Output:
|   - Modifies:
|     o refKmerAryUI and qryKmerAryUI in chainsSTPtr to
|       have all values set to -1
|     o chainLenArySI and scoreAryF in chainsSTPtr to have 
|       all values set to 0
|     o indexArySL in chainsSTPtr to be filled with -1's
|     o lenChainSL and lenIndexSL in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
blank_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to blank*/
);

/*-------------------------------------------------------\
| Fun02: init_chains_mapRead
|   - initialize a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct to initialize
| Output:
|   - Modifies:
|     o refKmerAryUI, qryKmerAryUI, chainLenArySI,
|       scoreAryF, and indexArySL to be 0 (null)
|     o lenChainSL, sizeChainSL, lenIndexSL, and
|       sizeIndexSL in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
init_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to initialize*/
);

/*-------------------------------------------------------\
| Fun03: freeStack_chains_mapRead
|   - frees all variables in a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead with variables to free
| Output:
|   - Modifies:
|     o frees refKmerAryUI, qryKmerAryUI, chainLenArySI,
|       scoreAryF, and indexArySL and sets to 0 (null)
|     o lenChainSL, sizeChainSL, lenIndexSL, and
|       sizeIndexSL in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
freeStack_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to free vars*/
);

/*-------------------------------------------------------\
| Fun04: freeHeap_chains_mapRead
|   - frees a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to free
| Output:
|   - Modifies:
|     o frees chainsSTPtr (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to free*/
);

/*-------------------------------------------------------\
| Fun05: setup_chains_mapRead
|   - sets up memory for a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to setup memory
|   - lenChainsSL:
|     o length of chains arrays
|   - lenRefUI:
|     o length of reference
| Output:
|   - Modifies:
|     o allocates lenChainsSL items of memory to
|       refCkmerArySI, qryKmerArySI, chainLenArySI,
|       and scoreAryF in chainsSTPtr
|     o allocates lenRefUI memory to indexArySL in
|       chainsSTPtr
|     o sizeChainSL in chainsSTPtr to be lenChainSL 
|     o sizeIndexSL in chainsSTPtr to be lenIndexSL 
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*to setup*/
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
);

/*-------------------------------------------------------\
| Fun06: mk_chains_mapRead
|   - makes a heap allocated chains_mapRead struct
| Input:
|   - lenChainsSL:
|     o length of chains arrays
|   - lenRefUI:
|     o length of reference
| Output:
|   - Returns:
|     o pointer to initialized and setup chains_mapRead
|       struct
|     o 0 for memory errors
\-------------------------------------------------------*/
struct chains_mapRead *
mk_chains_mapRead(
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
);

/*-------------------------------------------------------\
| Fun07: realloc_chains_mapRead
|   - reallocates memroy in a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to rellocate memory
|   - lenChainsSL:
|     o new length of chains arrays (0 for do not resize)
|   - lenRefUI:
|     o new length of reference (0 for do not resize)
| Output:
|   - Modifies:
|     o refKmerArySI, qryKmerArySI, chainLenArySI,
|       scoreAryF in chainsSTPtr to be reallocated to
|       lenChainSL if lenChainSL >chainsSTPtr->sizeChainSL
|     o indexArySL chainsSTPtr to be reallocated to
|       lenRefUI if lenRefUI  > chainsSTPtr->sizeIndexSL
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
realloc_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*reallocate*/
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
);

/*-------------------------------------------------------\
| Fun10: kmerChain_mapRead
|   - merge kmers into chains (longer kmers)
| Input:
|   - chainsSTPtr:
|     o pointer to chainsSTPtr struct to add kmer chains
|       to
|   - minChainLenSI:
|     o minimum number of kmers to keep a chain
|   - qryKmerArySI:
|     o signed int array with query kmers (by position)
|   - lenQrySI:
|     o number kmers in qryKmerArySI
|   - refKmerArySI:
|     o signed int array with reference kmers to chain
|   - refPosArySI:
|     o signed int array with index of reference kmers
|       (starts with 0, ends with lenRefSI - 1)
|     o for resorting refKmerArySI
|   - refSortArySI:
|     o signed int array that is same length as
|       refKmerArySI
|     o modified to have sorted kmers in refKmerArySI
|   - lenRefSI:
|     o number kmers in refKmerArySI
|   - refStartSI:
|     o first refrence base in window
| Output:
|   - Modifies:
|     o refKmerArySI, qryKmerArySI, and chainLenArySI in
|       chainsSTPtr to have chain starting coordinates and
|       length
|     o scoreArySI and nextArySI in chainsSTPtr to have
|       all detected kmer positions filled with 0's
|     o indexArySI in chainsSTPtr to have index or -1 for
|       each query kmer start position in qryKmerArySI
|     o lenChainSI and lenIndexSI to have length of
|       indexArySI (num ref kmers) and number of chains
|     o sizeChainSI and sizeIndexSI in chainsSTPtr if
|       needed to rellocate memory for chain/index arrays
|     o refSortArySI to have kmers in refKmerArySI
|       - kmers sorted, with negatives at end (unsigned)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|   - Note:
|     o sorting order of chains is query coordinate and
|       then reference coordinate
\-------------------------------------------------------*/
signed char
kmerChain_mapRead(
   struct chains_mapRead *chainsSTPtr,/*gets kmer chains*/
   signed int minChainLenSI,          /*min chain length*/
   signed int *qryKmerArySI,          /*query kmers*/
   signed int lenQrySI,               /*query length*/
   signed int *refKmerArySI,          /*reference kmers*/
   signed int *refSortArySI,          /*to sort kmers*/
   signed int lenRefSI,               /*number ref kmers*/
   signed int refStartSI              /*first ref base*/
);

/*-------------------------------------------------------\
| Fun12: mergeChains_mapRead
|   - merges kmer chains from fun10 (kmerChain) into
|     the chain the maximizes base count
|   - not optimal solution, but also reduce complexity
|     by not counting for gaps
| Input:
|   - chainAryUI:
|     o unsigned int array with kmer chains to merge into
|       a "best" chain
|   - numChainsSL:
|     o number of kmer chains in chainAryUI
|   - penaltyF:
|     o penalty for gap or mismatch (should be positive)
| Output:
|   - Modifies:
|     o nextArySI in chains_mapRead to have best scoring
|       trail for each chain
|     o scoreArySI to have most bases possible in a chain
|   - Returns:
|     o index of best merged chain
|     o def_noChains_mapRead if no chains in array
\-------------------------------------------------------*/ 
signed int
mergeChains_mapRead(
   struct chains_mapRead *chainsSTPtr  /*chains to merge*/
);

/*-------------------------------------------------------\
| Fun13: mergeToSam_mapRead
|   - convert output from mergeChains_map read to sam
|     file entry (mapq is set to 0)
| Input:
|   - chainsSTPtr:
|     o chains_mapRead struct pointer with chains to use
|       in making the alignment
|   - indexSI:
|     o starting index of merged chain to align
|   - refSTPtr:
|     o tbl_kmerFind struct pointer with reference
|       sequence
|   - qrySTPtr:
|     o ref_kmerFind struct used with kmerChain_mapRead
|       (is query)
|   - revBl:
|     o query direction (same as kmerChain_mapRead)
|       * 1 is reverse complement
|       * 9 is foward
|   - chainLenSI:
|     o minimum length of a chain
|     o used to align ends before chains, alignment size
|       is 1.5 chains (for chainLenSI = 20, is 30 bases)
|       * 1.5 is off the top of my head
|   - settings:
|     o settings for alignment
|   - samSTPtr:
|     o samEntry struct pointer to add alignment to
| Output:
|   - Modifies:
|     o samSTPtr to have alignment; everything in samSTPtr
|       is modified
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
mergeToSam_mapRead(
   struct chains_mapRead *chainsSTPtr,/*chains to use*/
   signed int indexSI,            /*chain index to start*/
   struct tblST_kmerFind *refSTPtr,/*reference sequence*/
   struct refST_kmerFind *qrySTPtr,/*has query sequence*/
   signed char revBl,            /*1 is reverse*/
   signed int chainLenSI,      /*min length of one chain*/
   struct alnSet *settings,      /*aligment settings*/
   struct samEntry *samSTPtr     /*gets final alignment*/
);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
| Fun14: flipRef_mapRead
|   - flips samEntry reference to query
| Input:
|   - qrySTPtr:
|     o seqST struct pointer with query sequence to add
|   - samSTPtr:
|     o samEntry struct pointer with sam entry to flip
|       * qrySTPtr must be reference mapped to in samSTPtr
| Output:
|   - Modifies:
|     o samSTPtr to have querySTPtr as sequence
|       * most variables are modified
|       * reverse complements if samSTPtr->flagUS & 16
|         - put in foward query sequence, this function
|           reverse complements the sequence copied to
|           samSTPtr
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
signed char
flipRef_mapRead(
   struct seqST *qrySTPtr,   /*query to add in*/
   struct samEntry *samSTPtr /*to flip*/
);

/*-------------------------------------------------------\
| Fun15: scoreSubAln_mapRead
|   - scores the aligment in samEntry and finds best
|     sub-aligment
| Input:
|   - samSTPtr:
|     o samEntry struct pointer with alignment to score
|     o do not trim
|   - subBl:
|     o modifiy alignment to best sub-alignment
|       * 1: look for best sub-alignment
|       * 0: do not look for sub-alignments
|   - refSTPtr:
|     o seqST struct pointer with reference sequence used
|       in alignment
|     o seqStr must be as lookup index (seqToIndex_alnSet)
|   - qrySTPtr:
|     o seqST struct pointer with query sequence used in
|       alignment
|     o seqStr must be as lookup index (seqToIndex_alnSet)
|   - alnSetSTPTr:
|     o alnSet struct pointer with scoring settings
|   - maxScoreSLPtr:
|     o signed long pointer to hold maximum score possible
|       for the aligned part of the query
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o extraStr in samSTPtr to have score "AS:i:score"
|       entry
|     o cigTypeStr and cigArySI in samSTPtr to be setup
|       for best sub-alignment (if subBl is 1)
|     o lenCigUI in samSTPtr to be new cigar length if
|       looking for best sub-alignment
|     o maxScoreSLPtr to have maximum possible score for
|       returned alignment (exlcudes unaligned bases)
|     o errSCPtr:
|       * 0 for no errors
|       * def_noScore_mapRead if no sub-score greater
|         than 0 (unlikley)
|       * def_memErr_mapRead for memory errors
|   - Returns:
|     o score for best sub-alignment or if subSl is 0, for
|       the entire samEntry
\-------------------------------------------------------*/
signed long
scoreSub_mapRead(
   struct samEntry *samSTPtr, /*alignment to score*/
   signed char subBl,         /*1:look for sub-alingment*/
   struct seqST *refSTPtr,    /*reference sequence*/
   struct seqST *qrySTPtr,    /*query sequence*/
   struct alnSet *alnSetSTPtr,/*scoring settings*/
   signed long *maxScoreSLPtr,/*gets max possible score*/
   signed char *errSCPtr      /*errors*/
);

/*-------------------------------------------------------\
| Fun16: blank_set_mapRead
|   - blank a set_mapRead struct
| Input:
|   - setSTPtr:
|     o pointer to set_mapRead struct to blank
|   - quickBl:
|     o tells if doing quick blank (used kmers only)
|       * 1: quick blank; only blank kmers in tblSTPtr
|       * 0: slow blank; blank all of tblSTPtr
|   - tblSeqBlankBl:
|     o tells if blanking sequence in table
|       * 1: blank sequence in tblSTPtr
|       * 0: do not blank sequece in tblSTPtr
| Output:
|   - Modifies:
|     o blanks structs in setSTPtr
\-------------------------------------------------------*/
void
blank_set_mapRead(
   struct set_mapRead *setSTPtr,
   signed char quickBl,/*1: do quick blank of tblSTPtr*/
   signed char tblSeqBlankBl /*blank sequences*/
);

/*-------------------------------------------------------\
| Fun17: init_set_mapRead
|   - initializes a set_mapRead struct
| Input:
|   - setSTPtr:
|     o pointer to set_mapRead struct to initializes
| Output:
|   - Modifies:
|     o sets variables to defaults and pointers to null
\-------------------------------------------------------*/
void
init_set_mapRead(
   struct set_mapRead *setSTPtr
);

/*-------------------------------------------------------\
| Fun18: freeStack_set_mapRead
|   - frees variables in a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer with variables to free
| Output:
|   - Modifies:
|     o frees all pointers and intitializes; defaults/null
\-------------------------------------------------------*/
void
freeStack_set_mapRead(
   struct set_mapRead *setSTPtr
);

/*-------------------------------------------------------\
| Fun19: freeHeap_set_mapRead
|   - frees a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o frees setSTPtr (you must set to null/0)
\-------------------------------------------------------*/
void
freeHeap_set_mapRead(
   struct set_mapRead *setSTPtr
);

/*-------------------------------------------------------\
| Fun20: setup_set_mapRead
|   - allocates memory for structs in a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer to free
|   - kmerUC:
|     o size of one kmer
|   - maxRefLenSL:
|     o size of expected longest reference
| Output:
|   - Modifies:
|     o all arrays in setSTPtr to have memory
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_set_mapRead(
   struct set_mapRead *setSTPtr, /*to setup*/
   unsigned char kmerUC,         /*size of one kmer*/
   signed long maxRefLenSL       /*size of longest ref*/
);

/*-------------------------------------------------------\
| Fun21: addQryKmers_mapRead
|   - get forward and reverse kmer sequences
| Input:
|   - forSeqStr:
|     o c-string with forward query sequence
|   - refSeqStr:
|     o c-string with reverse query sequence
|   - lenSeqSI:
|     o sequence length
|   - forKmerArySI:
|     o signed int array pionter to hold forward kmers
|   - revKmerArySI:
|     o signed int array pionter to hold reverse kmers
|   - sizeSIPtr:
|     o signed int pointer with current for/rev array size
|   - kmerLenUC:
|     o length of one kmer (index 1)
| Output:
|   - Modifies:
|     o forKmerArySI to have forward kmers (resized if
|       needed)
|     o revKmerArySI to have reverse kmers (resized if
|       needed)
|     o sizeSIPtr to have resized array length
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
addQryKmers_mapRead(
   signed char *forSeqStr,     /*forward sequence*/
   signed char *revSeqStr,     /*reverse sequence*/
   signed int lenSeqSI,        /*sequence length*/
   signed int **forKmerArySI,  /*gets forward kmers*/
   signed int **revKmerArySI,  /*gets reverse kmers*/
   signed int *sizeSIPtr,      /*current array sizes*/
   unsigned char kmerLenUC     /*size of one kmer*/
);

/*-------------------------------------------------------\
| Fun22: aln_mapRead
|   - maps a read to a reference
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to hold alignment
|   - qrySTPtr:
|     o seqST struct pointer with query sequence
|   - refSTPtr:
|     o seqST struct pointer with reference sequence
|   - setSTPtr:
|     o pointer to set_mapRead struct with settings
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o samSTPtr to have alignment (everything modified)
|     o errSCPtr to have error message
|       * 0 for no errors
|       * def_memErr_mapRead for no errors
|       * def_noChains_mapRead if chould not get chains
|       * def_noAln_mapRead if could not find alignment
|     o these are input to avoid making each time (reduce
|       memory allcation/deallocation), so no need to
|       check
|       * everything tblSTPtr
|       * everything in refKmerFindSTPtr
|       * everything in chainsSTPtr
|   - Returns:
|     o score for alignment
\-------------------------------------------------------*/
signed char
aln_mapRead(
   struct samEntry *samSTPtr,     /*gets aligment*/
   struct seqST *qrySTPtr,        /*query sequence*/
   struct seqST *refSTPtr,        /*reference sequence*/
   struct set_mapRead *setSTPtr,  /*settings*/
   signed char *errSCPtr          /*gets errors*/
);

/*-------------------------------------------------------\
| Fun23: qckAln_mapRead
|   - maps a read to a reference; returns best alignment
|   - qckAln_mapRead differs from aln_mapRead in that
|     it uses kmerScan to find best window in entire
|     reference, and then after scan, makes alignment
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to hold alignment
|   - qrySTPtr:
|     o seqST struct pointer with query sequence
|   - refSTPtr:
|     o seqST struct pointer with reference sequence
|   - setSTPtr:
|     o pointer to set_mapRead struct with settings
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o samSTPtr to have alignment (everything modified)
|     o errSCPtr to have error message
|       * 0 for no errors
|       * def_memErr_mapRead for no errors
|       * def_noChains_mapRead if chould not get chains
|       * def_noAln_mapRead if could not find alignment
|     o these are input to avoid making each time (reduce
|       memory allcation/deallocation), so no need to
|       check
|       * everything tblSTPtr
|       * everything in refKmerFindSTPtr
|       * everything in chainsSTPtr
|   - Returns:
|     o score for alignment
\-------------------------------------------------------*/
signed long
qckAln_mapRead(
   struct samEntry *samSTPtr,     /*gets aligment*/
   struct seqST *qrySTPtr,        /*query sequence*/
   struct seqST *refSTPtr,        /*reference sequence*/
   struct set_mapRead *setSTPtr,  /*settings*/
   signed char *errSCPtr          /*gets errors*/
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
: Copyright (c) 2025 jeremyButtler
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
