# Goal:

Explain how to include this code in your own C programs.

Not a very good guide, but it will tell you what most
  of the stuff does.

**TODO: update this**

# Structures:

## fluST:

### fluST structure

The fluST structure has the flu segment names, lengths,
  and the patterns to identify flu segments by primer
  sites. It also has the settings for the maximum mvRNA
  size and minimum number of deletions (how much smaller
  the read needs to be) to be a diRNA.

- fluST general variables:
  - fluST->segIDAryStr: has segment ids
  - fluST->lenSegArySI: has expected segment lengths
  - fluST->minDIDelSI: has minimum deletions to be diRNA
  - fluST->maxMVLenSI: has maximum length to be mvRNA

- fluST segment identification variables
  - fluST->forSearchST: linked list with patterns for the
    forward primers (5'->3' direction)
  - fluST->forCompSearchST: linked list with patterns for
    the forward primers when they are reverse complement
    (3'->5' direction)
  - fluST->revSearchST: linked list with patterns for the
    reverse primers (5'->3' direction)
  - fluST->revCompSearchST: linked list with patterns for
    the reverse primers when they are reverse complement
    (3'->5' direction)

### initializing fluST

You should always initialize a fluST with init_fluST
  (fun12 fluST.c/h). This will set all variables to the
  defaults in fluSeg.h.

### adding patterns for segment ids

You can add segment id sequences (patterns) or change the
  expected segment lengths with addSegTo_fluST
  (fun08 fluST.c/h).

- Input:
  - pointer to fluST structure to add pattern to
  - sequence (pattern) to add in
  - new expected length of the segment
  - 0 or 1
    - 0 is for a forward primer 
      - forSearchST and forCompSearchST are updated
    - 1 is for a reverse primer 
      - revSearchST and revCompSearchST are updated
  - segment number the pattern/sequence belongs to

### removing patterns for segment ids

You can remove segment id sequences (patterns) with
  rmSegTo_fluST (fun09 fluST.c/h). This does a lazy
  delete, so the pattern remains in the list, but its
  id is set to -1 (no match).

- Input:
  - pointer to fluST structure to remove pattern from
  - sequence (pattern) to remove
  - 0 or 1
    - 0 is for a forward primer 
      - forSearchST and forCompSearchST are updated
    - 1 is for a reverse primer 
      - revSearchST and revCompSearchST are updated
  - segment number the pattern/sequence belongs to

### searching for segment ids with fluST

When searching for ids it is best to use findSeg_fluST
  (fun10 fluST.c/h). The output is the segment number or
  -1 (def_noSeg_fluST from fluST.h). The segment number
  can be used to get the segment id and segment length
  from fluST. This is done in the detectDI_fluST step.

- Input:
  - fluST structure with patterns to search
  - last/first base in primer
    - last base when id sequence comes before primer
    - first base when id sequence comes after primer
  - 0 or 1
    - 0 if sequence is not reverse complement
      - searches forSearchST or revSearchST
    - 1 if sequence is reverse complement
      - searches forCompSearchST or revCompSearchST
  - 0 or 1
    - 0 if sequence is a forward primer
      - searches forSearchST or forCompSearchST
    - 1 if sequence is a reverse primer
      - searches revSearchST or revCompSearchST

### detecting DI segments with fluST

You can detect diRNA, mvRNA, or vRNA segments using
  detectDI_fluST (fun15). This returns a set of flags
  have two parts. The first part is how many primers
  supported the segment id. The second flag tells if the
  read was classified as diRNA, mvRNA, or vRNA.

- Return flags; segment id support:
  - 0 (no flags): could not id segment
  - def_segFound_fluST: both primers supported the same
    segment
  - def_partSeg_fluST: segment supported by one primer
    (the other primer supported no segments)
  - def_revSegSup_fluST: the segment was supported by the
    reverse primer
    - def_partSeg_fluST means id was from forward primer
    - def_revSegSup_fluST | def_partSeg_fluST means id was
      from reverser primer only
    - def_segFound_fluST always has def_revSegSup_fluST 
  - def_diffSeg_fluST both primers supported different
    segments.
    - only forward and reverse segment ids are returned
- Return flags; classification:
  - def_diFound_fluST: read is diRNA
  - def_mvFound_fluST: read is mvRNA
  - def_fullFound_fluST: read is vRNA (full length)

Some of the input for detectDI_fluST is from the primer
  scanning step, which can be done with
  waterFindPrims_kmerFind (fun20; kmerFind.c/h) or
  fxFindPrims_kmerFind (fun21; kmerFind.c/h).

- Input:
  - fluST structure to search for segment ids with and
    having thresholds for diRNA and mvRNA
  - last/first base in primer
    - last base when id sequence comes before primer
    - first base when id sequence comes after primer
  - hits; array of length 2 with 0 for no hits
    - is checked to see if both primers were found
    - return from the primer scanning step 
  - direction for each primer (array of length 2)
    - 'F' is forward
    - 'R' is reverse
    - returned from the primer scanning step 
  - first primer base in the sequence for each primer
    (array of length 2)
    - returned from the primer scanning step 
  - last primer base in the sequence for each primer
    (array of length 2)
    - returned from the primer scanning step 
  - a variable to hold the segment number (&variable)
    - set to def_noSeg_fluST if a read was assigned to
      two segments
  - a variable to hold the forwards primer segment number
    (&variable)
  - a variable to hold the reverse primer segment number
    (&variable)
  - a variable to hold the mapped length (&variable)
    - length from first primer base to last primer base
      

### printing results

To print the header use pidHeader_fluST (fun16 fluST.c/h).

- Input:
  - file to output to

To print the results of detectDI_fluST use pid_fluST
  (fun17 fluST.c/h). For disagreements this will also find
  if the longest supported segment is diRNA.

- Input:
  - fluST structure with segment lengths, ids, and
    settings for DI checking
  - read id to print out (first character is ignored)
  - segment number for read (both primers in agreement)
  - segment number for forward primer of read
  - segment number for reverse primer of read
  - returned flag from detectDI_fluST (classification)
  - direction primer mapped onto sequence (array length 2)
    - returned from the primer scanning step 
  - score for both primers (array of length 2)
    - returned from the primer scanning step 
  - maximum score possible for forward primer
  - maximum score possible for reverse primer
  - primers starting position on sequence
    (array of length 2)
    - returned from the primer scanning step 
  - primers ending position on sequence
    (array of length 2)
    - returned from the primer scanning step 
  - sequences starting position on primer
    (array of length 2)
    - returned from the primer scanning step 
  - sequences ending position on primer
    (array of length 2)
    - returned from the primer scanning step 
  - length of sequence
  - mapped length from detectDI_fluST
  - file to print id to

Returns nothing.

### Freeing fluST structure

You can free the fluST structure using the freeHeap_fluST
  (fun14 fluST.c/h) or freeStack_fluST (fun13 fluST.c/h)
  functions. Use freeHeap_fluST when you allocated the
  fluST structure to the heap (malloc/calloc). Make sure
  to set the fluST pointer to null afterwards. Use
  freeStack_fluST when fluST is on the stack
  (struct fluST variable;).

- Input (same for both):
  - fluST structure to free (as pointer)

## seqStruct

The seqStruct structure is set up to hold a fastq or fasta
  sequence.

- Variables:
  - seqStruct->idStr has the read id
    - this includes the starting '>' or '@' and all
       white space
  - seqStruct->lenIdUL length of read id
  - seqStruct->seqStr has the sequence with white space
    removed (as much as possible)
  - seqStruct->lenSeqUL length of sequence
  - seqStruct->qStr has the q-score entry or is 0 or if
    buffer set to '\0'
  - seqStruct->lenQUL length of q-score entry
  - seqStruct->offsetUL is for memwater
    - tells were to start the alignment
  - seqStruct->endAlnUL is for memwater
    - tells were to end the alignment

### initializing seqStruct

You should always call init_seqST (fun07 seqST.h) before
  using a seqStruct for the first time. This sets all
  pointers to null.

- Input:
  - seqStruct to initialize

### Reading in fastq/fasta sequences

You can read in fastq sequences with getFqSeq_seqST (fun02
  src/seqST.c/h). For fasta sequences use getFaSeq_sedqST
  (fun03 src/seqST.c/h).

- Input:
  - fastq (getFq) or fasta (getFa) file to get sequence
    from
  - seqStruct to hold the sequence

### Reverse complement

You can reverse complement a sequence using revComp_seqST
  (fun05 src/seqST.c/h).

- Input:
  - seqStruct to reverse complement

### converting sequence to index for alignment

For the waterman alignment (or any aligner from alnSeq)
  you will need to convert the sequence to look up
  index's. You can do this with seqToIndex_alnSetST
  (fun09 memwater/alnSetST.h).

To convert the index back use indexToSeq_alnSetST (fun10
  memwater/alnSetST.h).

- Input (both are same):
  - sequence to convert (seqStruct->seqStr).

### Freeing seqStruct

For seqStructs on the heap (malloc/calloc) use
  freeHeap_seqST (fun09 src/seqST.c/h). For arrays of
  seqStructs allocated on the heap use freeAry_seqST
  (fun10 src/seqST.c/h). For seqStructs on the stack use
  freeStack_seqST (fun08 src/seqST.c/h).

- Input:
  - seqStruct to free
  - for freeAry_seqST the number of seqStructs in the
    array

## alnSet structure

The alnSet structure has the settings for the memwater
  alignment. Most of the settings in it were designed
  for alnSeq and do not affect kmer find.

You can change the score for a match with
  setScore_alnSetST (fun01 alnSetST.h).

You can change the gap extension penalty with
  alnSet->gapExtendC. The gap opening penalty can be
  changed with alnSet->gapOpenC.

All scores and gap penalties are signed chars, meaning
  they can be between -127 to 128.

The defaults for the alnSet structure are in
  "memwater/alnDefaults.h". Though both primFind and
  getDIIds modify the gapopen penalty to be -4 instead of
  -10.

### initialize and freeing a alnSet structure

You should initialize the alnSet structure using
  init_alnSetST (fun11 memwater/alnSetST.c/h).

There are no heap variables in a alnSet structure, but if
  you want to make sure it is free you can call
  freeHeap_alnSetST (fun04 memwater/alnSetST.c/h) for
  heap allocated structures or freeStack_alnST (fun03
  memwater/alnSetST.c/h) for stack structures.

For heap frees make sure to set the alnSet structure
  pointer to null after freeHeap_alnSetST.

## refST_kmerFind structure

The refST_kmerFind structure holds the reference (primer)
  sequences to use in primer searching.

- Variables:
  - refST_kmerFind->mateSI stores the index of the
    reverse/forward primer to primer is paired to
    - Only modify it if you want the read to have both
      primers
    - requires a refST_kmerFind array
    - Set it to -1 for primers not paired
  - refST_kmerFind->maxForScoreF maximum possible score
    for a Waterman alignment of stored sequence (forward)
  - refST_kmerFind->maxRevScoreF maximum possible score
    for a Waterman alignment of stored sequence (reverse
    complement)
- refST_kmerFind->forSeqST is a seqStruct with the forward
  sequence
- refST_kmerFind->revSeqST is a seqStruct with the reverse
  complement sequence

### add sequences to refST_kmerFind

First off you need to initialize the refST_kmerFind
  structure with init_refST_kmerFind (fun06 kmerFind.c/h).
  The refST_kmerFind is not set up to change sequences in,
  so you will have to do a stack free (will set variables
  to defaults).

You can add sequences to the refST_kmerFind structure by
  hand using addSeqToRefST_kmerFind (fun11 kemrFind.c/h).
  However, if you do this you will have to set the
  refST_kmerFind->mateSI flag by hand (if you want it).

- Output:
  - Returns the length of the longest sequence added to
    the window
    - You will need this later to set up tblST_kmerFind
    - is 0 if there was an error

- Input:
  - tblST_kmerFind structure (next structure) to use with
    the refST_kmerFind structure/group of structures
  - refST_kmerFind structure to add sequence to
  - sequence in seqStruct to add to the refST_kmerFind
    structure
    - this is copied, so you can free this afterwards
  - minimum percentage of kmers needed to consider a
    window a match
  - the length of the longest sequence added to a previous
    refST_kmerFind structure
  - settings for the Waterman alignment (alnSet structure)
    - use to find the maximum possible score

### getting sequences from file

You can get sequences from a tsv file using
  tsvToAry_refST_kmerFind (fun13 kmerFind.c/h) or you 
  can use a fasta file with faToAry_refST_kmerFind
  (fun14 kmerFind.c/h). The advantage about the tsv is
  that you can specify if the primers are paired or not.
  The fasta file method assumes nothing is paired.
 
- tsv format:
  - column 1 is the primer id
  - column 2 has "True", "Yes", or "1" for paired primers
    - It really is "T", "Y", "1"
  - column 3 has the forward primer sequence
    - I think you can use "NA" for none
  - column 4 has the reverse primer sequence
    - I think you can use "NA" for none

- Output (both are same):
  - an array of refST_kmerFind structures with sequences
    - 0 for an error
  - sets up the tblST_kmerFind structure

- Input:
  - file to read in (fasta/tsv depending on function)
  - kmer length for kmer search
  - variable to hold the refST_kmerFind array length
    - used with freeHeapAry_kmerFind
  - minimum percentage of kmers needed in one read window
    to consider a potential match
  - a tblST_kmerFind structure to set up (next structure)
  - the percentage of extra nucleotides in window
    - (longest reference length) * percentage
  - the percentage of nucleotides to shift window by 
    - (window size) * percentage
  - settings for the Waterman alignment (alnSet structure)
    - use to find the maximum possible score
  - a variable to hold an errors; check to find error type

### Freeing a refST_kmerFind structure

You can free a refST_kmerFind structure with
  freeHeap_refST_kmerFind (heap allocated; fun08
  kmerFind.c/n), freeHeapAry_refST_kmerFind (heap
  allocated array (tsvToAry or faToAry);
  fun09 kmerFind.c/h), or freeStack_refST_kmerFind (not
  on heap; fun07 kmerFind.c/h).

- Input:
  - refST_kmerFind to free
  - for array free; number of structures in array

## tblST_kmerFind

The tblST_kmerFind structure is used to do the kmer search
  on the sequence.

- Variables:
  - tblST_kmerFind->seqST seqStruct that holds the read
    being checked
  - tblST_kmerFind->seqPosUl is the position at in the
    read

### initialization and set up tblST_kmerFind

You can initialize a new tblST_kmerFind structure using
  init_tlbST_kmerFind (fun02 kmerFind.c/h).

You can set up a tblST_kmerFind structure using
  setUpTblST_kmerFind (fun12 kmerFind.c/h). However, the
  tsvToAry_refST_kmerFind and faToAry_refST_kmerFind will
  do this for you.

- Input (setUpTblST_kmerFind)
  - tblST_kmerFind structure to initialize
  - the percentage of extra nucleotides in window
    - (longest reference length) * percentage
  - the percentage of nucleotides to shift window by 
    - (window size) * percentage
  - the longest sequence in your refST_kmerFind structure
    array

### freeing tblST_kmerFind

For freeing a tblST_kmerFind on the heap use
  freeHeap_tblST_kmerFind (fun04 kmerFind.c/h). For
  freeing from the stack use freeStack_tblST_kmerFind
  (fun03 kmerFind.c/h). For

For heap frees, remember to set you tblST_kmerFind
  structure pointer to null.

# Functions:

## Primer searching

For primmer searching you have two functions,
  waterFindPrims_kmerFind (fun20 kmerFind.c/h) and
  fxFindPrims_kmerFind (fun21 kmerFind.c/h). The
  waterFindPrims_kmerFind function uses a memwater to find
  the primer score and coordinates. While
  fxFindPrims_kmerFind uses a kmer scan to find potential
  windows (I like double longest reference size) and then 
  does a waterman alignment on potential windows.

- Returns:
  - 0 if no primers found
  - 1 at least one primer or primer pair was found

- Input:
  - tblST_kmerFind structure to scan with
    - fxFindPrims_kmerFind only
  - refST_kmerFind array of primers to search for
  - number of primers in the refST_kmerFind array
  - sequence to scan (in seqStruct)
  - minimum percent score to count a primer as mapped
    - found by score / (maximum score)
  - array to hold the number of times each primer mapped
    to a read
    - same length as refST_kmerFind array
  - array to hold the direction of each mapped primer
    - same length as refST_kmerFind array
  - array to hold the score of each mapped primer
    - same length as refST_kmerFind array
  - array to hold the first sequence base that a primer
    mapped to for each mapped primer
    - same length as refST_kmerFind array
  - array to hold the last sequence base that a primer
    mapped to for each mapped primer
    - same length as refST_kmerFind array
  - array to hold the first primer base that a sequence
    mapped to for each mapped primer
    - same length as refST_kmerFind array
  - array to hold the last primer base that a sequence
    mapped to for each mapped primer
    - same length as refST_kmerFind array
  - settings for the alignment (alnSet structure)

## printing primer results

Use pHEaderHit_kmerFind (fun23 kmerFind.c/h) to print the
  header.

- Input
  - file to print header to

Use phit_kmerFind (fun22 kmerFind.c/h) to print the
  results of a single read

- Input
  - refST_kmerFind array with primer ids (in seqStruct)
  - number of primers in the refST_kmerFind array
  - sequence to print out in a seqStruct
  - array having the number of times each primer mapped
    to a read
    - same length as refST_kmerFind array
  - array having the direction of each mapped primer
    - same length as refST_kmerFind array
  - array having the score of each mapped primer
    - same length as refST_kmerFind array
  - array having the first sequence base that a primer
    mapped to for each mapped primer
    - same length as refST_kmerFind array
  - array having the last sequence base that a primer
    mapped to for each mapped primer
    - same length as refST_kmerFind array
  - array having the first primer base that a sequence
    mapped to for each mapped primer
    - same length as refST_kmerFind array
  - array having the last primer base that a sequence
    mapped to for each mapped primer
    - same length as refST_kmerFind array
  - file to print the hit to

## Waterman alignment

primFind and getDIIds use the memwater alignment. It only
  returns the score and alignment coordinates (start and
  end). This takes a little longer than a traditional
  Waterman alignment, so it is slower than the striped
  Waterman.  You can call memwater with memWater (fun09
  memwater/memwater.c/h).

The mem in memwater is for memory efficient. It comes from
  the memwater alignment using about twice as much memory
  as a Hirschberg.

- Input:
  - query sequence in a seqStruct
    - offsetUL set to first base to align in query
    - endAlnUL set to last base to align in query
  - reference sequence in a seqStruct
    - offsetUL set to first base to align in reference
    - endAlnUL set to last base to align in reference
  - variable to hold first mapped base in reference
  - variable to hold last mapped base in reference
  - variable to hold first mapped base in query
  - variable to hold last mapped base in query

- Output:
  - score for the alignment

# Side libraries (generalLib):

These libraries are here for general support. In this
  case they are only .h files.

- shellShort: has a binary search, swap function, and
  shell sort functions (min to max)
- base10str: converts string with base 10 numbers to
  integers
- ulCp: has long bases copy functions that allow you to
  change the deliminator. Not as fast as strcpy, but not
  limited to '\0'.
  - also has a find end of line function ('\n' or '\0')
  - returns number of bytes copied
- charCp: byte version of ulCp, but also includes a
  function for comparing strings
  - returns number of bytes copied
- dataTypeShortHand:
  - my typedef short hand for variable names
  - unsigned type goes to utype
  - signed type goes to stype
- genMath: does branch less max, mins, and absolute values
  - from
    [https://graphics.stanford.edu/~seander/bithacks.html](
    https://graphics.stanford.edu/~seander/bithacks.html)
- ntToTwoBit: look up table to convert nucleotide to a two
  bit value. A third bit is used to report error or
  anonymous bases.
