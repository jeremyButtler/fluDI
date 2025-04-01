/*-------------------------------------------------------\
' fluST SOF: Start Of File
'   - holds fluST structure and the supporting structures
'     for finding flu segment once primer location is
'     found
'   o header:
'     - guards and defined variables
'   o .h st01: segIdSearch_fluST
'     - is a single link in a flu segment search
'   o .h st02: fluST
'     - has the segement id list to search and segments
'   o fun01: blank_segIdSearch_fluST
'     - sets a segIdSearch_fluST structure to defaults
'   o fun02: init_segIdSearch_fluST
'     - initializes variables in a segIdSearch_fluST
'       structure
'   o fun03: freeStack_segIdSearch_fluST
'     - free variables in a segIdSearch_fluST structure
'     - This will free all nodes, except the previous node
'       in the linked list
'   o fun04: freeHeap_segIdSearch_fluST
'     - free variables in a segIdSearch_fluST structure
'     - This will free all nodes, except previous nodes
'       in the linked list
'   o .c fun05: addSeqTo_segIdSeach_fluST
'     - adds a new pattern (sequence) to a
'       segIdSearch_fluST linked list
'   o .c fun06: findSeq_segIdSeach_fluST
'     - finds a sequence in segIdSearch_fluST linked list
'   o .c fun07: cpCompSeq_fluST
'     - copies the complement sequence to a new string.
'       anonymous bases and '-'s are convereted to 'N's
'   o fun08: addSegTo_fluST
'     - adds a segment to a fluST structure
'   o fun09: rmSegSeqFrom_fluST
'     - removes a segment sequence from a fluST structure
'   o fun10: findSeg_fluST
'     - find the segement a sequence belongs to
'   o fun11: blank_fluST
'     - sets a fluST stucture to defaults
'   o fun13: setup_fluST
'     - allocates and sets up a fluST struture to default
'       values (call init_fluST first)
'   o fun14: freeStack_fluST
'     - frees all variables inside a fluST structure
'   o fun15: freeHeap_fluST
'     - frees a fluST structure
'   o fun16: detectDI_fluST
'     - finds segment number and then detects if read is
'       a diRNA, mvRNA, or a vRNA (full)
'   o fun17: pidHeader_fluST
'     - prints out the header for pid_fluST
'   o fun18: pid_fluST
'     - prints out read id and stats for a flu segment
'   o license:
'     - licensing for this code (public domain / mit)
\-------------------------------------------------------*/

/*-------------------------------------------------------\
| Header:
|   - guards and defined variables
\-------------------------------------------------------*/

#ifndef FLU_STRUCTURE_H
#define FLU_STRUCTURE_H

#define def_segFound_fluST 1   /*prims support same seg*/
#define def_partSeg_fluST 2    /*1 prim has segment type*/
#define def_revSegSup_fluST 8  /*rev primer supports*/
#define def_diFound_fluST 16   /*diRNA found*/
#define def_fullFound_fluST 32 /*found a full segment*/
#define def_multiSeg_fluST 64  /*to many segments*/

#define def_diMinDel_fluST 50 /*50 nucleotides off full*/
#define def_maxMvRNALen_fluST 127

#define def_lenId_fluST 4
#define def_numSeg_fluST 8

/*numbers for each segment*/
#define def_noSeg_fluST -1 /*DO NOT CHANGE FROM -1*/

/*-------------------------------------------------------\
| ST01: segIdSearch_fluST
|   - is a single link in a flu segment search
\-------------------------------------------------------*/
typedef struct segIdSearch_fluST
{
   signed char segArySC[def_numSeg_fluST];
       /*matching segment numbers*/
   unsigned char numSegsUC; /*# segments with pattern*/

   /*the four next possible bases*/
   struct segIdSearch_fluST *aNtST;
   struct segIdSearch_fluST *cNtST;
   struct segIdSearch_fluST *gNtST;
   struct segIdSearch_fluST *tNtST;

   /*previous structure for backtracking*/
   struct segIdSearch_fluST *lastNtST; 
}segIdSearch_fluST;

/*-------------------------------------------------------\
| ST02: fluST
|   - has the segement id list to search and the segments
\-------------------------------------------------------*/
typedef struct fluST
{
   /*list to id flu segments with*/
   struct segIdSearch_fluST *forSearchST;
   struct segIdSearch_fluST *forCompSearchST;

   struct segIdSearch_fluST *revSearchST;
   struct segIdSearch_fluST *revCompSearchST;

   /*c-string array of segment names; use segArySC from
   `  segIdSearch to get id from segIdAryStr
   */
   signed char
      segIdAryStr[def_numSeg_fluST][def_lenId_fluST];

   /*array of segment lengths; use segArySC from
   `   segIdSearch to get length of segment
   */
   signed int lenSegArySI[def_numSeg_fluST];

   /*variables for detecting di/mcRNAs*/
   signed int minDiDelSI; /*min deletions for di*/
   signed int maxMvLenSI; /*max length for for mvRNA*/
}fluST;

/*-------------------------------------------------------\
| Fun01: blank_segIdSearch_fluST
|   - sets a segIdSearch_fluST structure to defaults
| Input:
|   - segIdSearch_fluST_ptr:
|     o pointer to a segIdSearch_fluST to set to defaults
| Output:
|   - Modifies:
|     o values in segIdSearch_fluST_ptr to be defaults
\-------------------------------------------------------*/
void
blank_segIdSearch_fluST(
   struct segIdSearch_fluST *segIdSearch_fluST_ptr
);

/*-------------------------------------------------------\
| Fun02: init_segIdSearch_fluST
|   - initializes variables in a segIdSearch_fluST
|     structure
| Input:
|   - segIdSearch_fluST_ptr:
|     o pointer to a segIdSearch_fluST to initialize
| Output:
|   - Modifies:
|     o values in segIdSearch_fluST_ptr to be defaults or
|       null
\-------------------------------------------------------*/
void
init_segIdSearch_fluST(
   struct segIdSearch_fluST *segIdSearch_fluST_ptr
);

/*-------------------------------------------------------\
| Fun03: freeStack_segIdSearch_fluST
|   - free variables in a segIdSearch_fluST structure
|   - This will free all nodes, except the previous node
|     in the linked list
| Input:
|   - segSTPtr:
|     o pointer to a segIdSearch_fluST to free
| Output:
|   - Frees:
|     o the A, C, G, and T linked lists in segSTPtr
|   - Modifies:
|     o sets segArySC to 0
|     o sets the A, C, G, and T pointers to 0
\-------------------------------------------------------*/
void
freeStack_segIdSearch_fluST( 
   struct segIdSearch_fluST *segSTPtr
);

/*-------------------------------------------------------\
| Fun04: freeHeap_segIdSearch_fluST
|   - free variables in a segIdSearch_fluST structure
|   - This will free all nodes, except previous nodes
|     in the linked list
| Input:
|   - segSTPtr:
|     o pointer to a segIdSearch_fluST to free
| Output:
|   - Frees:
|     o segSTPtr and all values in it
|       - you must set segSTPtr to null
|   - Modifies:
|     o segSTPtr->lastNtST to have the link to segSTPtr
|       set to null
\-------------------------------------------------------*/
void
freeHeap_segIdSearch_fluST( 
   struct segIdSearch_fluST *segSTPtr
);

/*-------------------------------------------------------\
| Fun08: addSegTo_fluST
|   - adds a segment to a fluST structure
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to add segment to
|   - seqStr:
|     o c-string with the sequence to add
|   - lenSegSI:
|     o expected length of the segment
|   - revPrimBl:
|     o direction of primer (1 = reverse; 0 = foward)
|   - segSC:
|     o segment id to add; use def_XXXNum_fluSeg variables
| Output:
|   - Modifies:
|     o all variables in fluST to support the new pattern
|   - Returns:
|     o 0 for no errors
|     o 1 for memory error
|     o 2 for invalid base in a seqStr
\-------------------------------------------------------*/
signed char
addSegTo_fluST(
   struct fluST *fluSTPtr,
   signed char *seqStr,
   signed int lenSegSI,
   signed char revPrimBl,
   signed char segSC
);

/*-------------------------------------------------------\
| Fun09: rmSegSeqFrom_fluST
|   - removes a segment sequence from a fluST structure
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to remove segment
|       sequence
|   - seqStr:
|     o c-string with the sequence to remove
|   - revPrimBl:
|     o direction of primer (1 = reverse; 0 = foward)
|   - segArySC:
|     o segment id to add; use def_XXXNum_fluSeg variables
| Output:
|   - Modifies:
|     o removes the sequence from the segment linked list
|   - Returns:
|     o 0 for no errors
|     o 1 for sequence not in list
\-------------------------------------------------------*/
signed char
rmSegSeqFrom_fluST(
   struct fluST *fluSTPtr,
   signed char *seqStr,
   signed char revPrimBl,
   signed char rmSegNumSC
);

/*-------------------------------------------------------\
| Fun10: findSeg_fluST
|   - find the segement a sequence belongs to
| Input:
|   - fluSTPtr:
|     o pointer to fluST structure to search for segment
|       in
|   - seqStr:
|     o c-sting pointing to first base that is part of
|       segment pattern in the sequence
|   - cmpBl:
|     o 1: sequence is a complement sequence
|     o 0: sequence is normal
|   - revBl:
|     o 1: is a reverse primer
|     o 0: forward primer
|   - numSegsUCPtr:
|     o pointer to unsigned char to hold number of
|       segments with the pattern
| Output:
|   - Modifies:
|     o numSegUCPtr to have the number of segments in the
|       returned array
|   - Returns:
|     o pointer to array with the segment number
|     o 0 if no segment found
\-------------------------------------------------------*/
signed char *
findSeg_fluST(
   struct fluST *fluSTPtr,
   signed char *seqStr,
   signed char cmpBl,
   signed char revBl,
   unsigned char *numSegsUCPtr
);

/*-------------------------------------------------------\
| Fun11: blank_fluST
|   - sets a fluST stucture to defaults, but really does
|     nothing, since no variables in fluST have defaults
|     that should be changed
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to blank
| Output:
|   - Nothing:
|     o this is here for future changes
\-------------------------------------------------------*/
void
blank_fluST(
   struct fluST *fluSTPtr
);

/*-------------------------------------------------------\
| Fun12: init_fluST
|   - initializes all values in a fluST struct to 0
| Input:
|   - fluSTPtr:
|     o pointer to fluST structure to initialize
| Output:
|   - Modifies:
|     o all values in fluSTPtr to be defaults or 0
\-------------------------------------------------------*/
void
init_fluST(
  struct fluST *fluSTPtr
);

/*-------------------------------------------------------\
| Fun13: setup_fluST
|   - allocates and sets up a fluST struture to default
|     values (call init_fluST first)
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to initialize
| Output:
|   - Modifies:
|     o all values in fluSTPtr to be defaults.
|   - Returns:
|     o 0 for no errors
|     o 1 for memory error
|     o 2 for invalid base in a seqStr
\-------------------------------------------------------*/
signed char
setup_fluST(
  struct fluST *fluSTPtr
);

/*-------------------------------------------------------\
| Fun14: freeStack_fluST
|   - frees all variables inside a fluST structure
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to free variables in
| Output:
|   - Frees:
|     o all variables in fluSTPtr and sets pointers to 0
\-------------------------------------------------------*/
void
freeStack_fluST(
   struct fluST *fluSTPtr
);

/*-------------------------------------------------------\
| Fun15: freeHeap_fluST
|   - frees a fluST structure
| Input:
|   - fluSTPtr:
|     o pointer to a fluST structure to free
| Output:
|   - Frees:
|     o fluSTPtr (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_fluST(
   struct fluST *fluSTPtr
);

/*-------------------------------------------------------\
| Fun16: detectDI_fluST
|   - finds segment number and then detects if read is
|     a diRNA or a vRNA (full)
| Input:
|   - fluSTPtr:
|     o pointer to fluST structure with thresholds and
|       segment number
|   - seqStr:
|     o c-string with sequence to see if is diRNA
|   - hitAryUI:
|     o unsigned in array with number of hits for both
|       primers
|   - dirArySC:
|     o signed char array with the direction "F" for
|       foward of the best mapped primer
|   - seqStartAryUL:
|     o unsigned long array with the primers starting
|       coordinate on the sequence
|   - seqEndAryUL:
|     o unsigned long array with the primers ending
|       coordinate on the sequence
|   - segSC:
|     o pointer to signed char to hold the segment number
|       decided on
|   - mappedLenUL:
|     o pionter to unsigned long to hold the length of
|       region from start to end of primers
| Output:
|   - Modifies:
|     o segArySC to have the segment number
|   - Returns:
|     o 0 if both primers were not detected, if primers
|       are in same direction, or if no segments were
|       found or if have multiple segments for one primer
|     o segement stats | genome type
|       - segment stats:
|         o 0 for no segment found
|         o def_segFound_fluST | def_revSegSup_fluST if
|           both primers support the same segment
|         o def_partSeg_fluST if only forward primer
|           supported a segment
|         o def_partSeg_fluST | def_revSegSup_fluST if
|           only reverse primer supported a segment
|         o def_noSeg_fluST if neither primer supports a
|           segment
|         o def_multiSeg_fluST if more than one segment
|           was supported
|       - genome type:
|         o def_diFound_fluST if meets diRNA length
|         o def_mvFound_fluST if is mvRNA length
|         o def_fullFound_fluST if is larger then diRNA
|           length (full genome)
\-------------------------------------------------------*/
signed char
detectDI_fluST(
   struct fluST *fluSTPtr,
   signed char *seqStr,
   unsigned int hitAryUI[],       /*primer hits*/
   signed char dirArySC[],        /*primer direction*/
   unsigned long seqStartAryUL[], /*seq map coordiantes*/
   unsigned long seqEndAryUL[],   /*seq map coordiantes*/
   signed char *segSC,
   unsigned long *mappedLenUL
);

/*-------------------------------------------------------\
| Fun17: pidHeader_fluST
|   - prints out the header for pid_fluST
| Input:
|   - outFILE:
|     o file to print header to
| Output:
|   - Prints:
|     o header to outFILE
\-------------------------------------------------------*/
void
pidHeader_fluST(
   void *outFILE
);

/*-------------------------------------------------------\
| Fun18: pid_fluST
|   - prints out read id and stats for a flu segment
| Input:
|   - fluSTPtr:
|     o pointer to a fluSTPtr structure with lengths and
|       segment names
|   - idStr:
|     o c-string with read id
|   - segArySC:
|     o segment number found
|   - forSegSC:
|     o segment number of the foward primer
|   - revSegSC:
|     o segment number of the reverse primer
|   - diFlagSC:
|     o return from detect_DI_fluST (fun16)
|   - dirArySC:
|     o array of primer directions
|   - scoreArySL:
|     o array of alignment scores
|   - maxForScoreF:
|     o maximum score for foward primer
|   - maxRevScoreF:
|     o maximum score for reverse primer
|   - seqStartAryUL:
|     o array with primer starting positions on sequence
|   - seqEndAryUL:
|     o array with primer ending positions on sequence
|   - primStartAryUL:
|     o array with sequence starting positions on primer
|   - primEndAryUL:
|     o array with sequence ending positions on primer
|   - seqLenUL:
|     o length of sequence
|   - mappedLenUL:
|     o length of region between primers (start to end)
|   - outFILE:
|     o file to print id and stats to
| Output:
|   - Prints:
|     o id and other stats to outFILE
\-------------------------------------------------------*/
void
pid_fluST(
   struct fluST *fluSTPtr,
   signed char *idStr,
   signed char segArySC,
   signed char diFlagSC,
   signed char dirArySC[],
   signed long scoreArySL[],
   float maxForScoreF,
   float maxRevScoreF,
   unsigned long seqStartAryUL[],
   unsigned long seqEndAryUL[],
   unsigned long primStartAryUL[],
   unsigned long primEndAryUL[],
   unsigned long seqLenUL,
   unsigned long mappedLenUL,
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
