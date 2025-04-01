/*#######################################################\
# Name: diCoords
#   - has functions to get DI coordinates from a sam
#     file
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header
'     - guards
'   o fun01: get_diCoords
'     - gets start and ending coordinates for DI events
'   o fun02: scan_diCoords
'     - scans to see if DI events and returns the number
'       found (get_diCoords without memory allocation)
'   o fun03: pDIHead_diCoords:
'     - print out diCoords tsv header
'   o fun04: pDI_diCoords:
'     - print out a di entry to a file
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards
\-------------------------------------------------------*/

#ifndef GET_DEFECTIVE_INTERFEARING_COORDINATES_H
#define GET_DEFECTIVE_INTERFEARING_COORDINATES_H

/*-------------------------------------------------------\
| Fun01: get_diCoords
|   - gets start and ending coordinates for DI events
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct with sequencer to scan
|   - minDILenUI:
|     o minimum deletion length to classify cigar entry as
|       a DI entry
|   - minEndNtUI:
|     o how many bases in a DI event must be to be a DI
|       event
|       event
|   - indelFlagSC:
|     o holds what type of indels targeting
|       - 1: dels
|       - 2: ins
|       - 3: del/ins
|   - diStartAryUI:
|     o pointer to unsigned long to hold the starting
|       coordinates of each DI event
|   - diEndAryUI:
|     o pointer to unsigned long to hold the ending
|       coordinates of each DI event
|   - delSizeAryUI:
|     o pointer to unsigned int array to hold the deletion
|       size of each DI event
|   - indelArySC:
|     o holds indels types found ("D" = del; "I" = ins)
|   - lenArysUIPtr:
|     o pointer to unsigned long to hold diStartAryUI and
|       diEndAryUI lengths
| Output:
|   - Modifies:
|     o diStartAryUI to have DI starting coordinates
|     o diEndAryUI to have DI ending coordinates
|     o lenArysUIPtr to have new array lengths if
|       diStartAryUI and diEndAryUI are changed
|   - Returns:
|     o number of DI events detected
|     o < 0 for memory errors
\-------------------------------------------------------*/
signed int
get_diCoords(
   struct samEntry *samSTPtr,  /*sam entry to scan*/
   unsigned int minDILenUI,    /*min del length for DI*/
   unsigned int minEndNtUI,    /*max NT at ends*/
   signed char indelFlagSC,    /*holds indel type*/
   unsigned int **diStartAryUI,/*gets start coordinates*/
   unsigned int **diEndAryUI,  /*gets DI end coordinates*/
   unsigned int **delSizeAryUI,/*gets DI deletion sizes*/
   signed char **indelArySC,   /*indel types found*/
   unsigned int *lenArysUIPtr  /*current arrays lengths*/
);

/*-------------------------------------------------------\
| Fun02: scan_diCoords
|   - scans to see if DI events and returns the number
|     found (get_diCoords without memory allocation)
| Input:
|   - samSTPtr:
|     o pointer to samEntry struct with sequencer to scan
|   - minDILenUI:
|     o minimum deletion length to classify cigar entry as
|       a DI entry
|   - minEndNtUI:
|     o how many bases in a DI event must be to be a DI
|       event
| Output:
|   - Returns:
|     o number of DI events detected
|     o < 0 for memory errors
\-------------------------------------------------------*/
signed int
scan_diCoords(
   struct samEntry *samSTPtr,  /*sam entry to scan*/
   unsigned int minDILenUI,    /*min del length for DI*/
   unsigned int minEndNtUI     /*max NT at ends*/
);

/*-------------------------------------------------------\
| Fun03: pDIHead_diCoords:
|   - print out diCoords tsv header
| Input:
|   - outFILE:
|     o file to print to
| Output:
|   - Prints:
|     o header for pDI_diCoords to outFILE
\-------------------------------------------------------*/
void
pDIHead_diCoords(
   void *outFILE
);

/*-------------------------------------------------------\
| Fun04: pDI_diCoords:
|   - print out a di entry to a file
| Input:
|   - qryIdStr:
|     o c-string with query id to print out
|   - refIdStr:
|     o c-string with reference id to print out
|   - diStartAryUI:
|     o unsigned int array with starting coordinates
|       for each DI event (get_diCoords fun01)
|   - diEndAryUI:
|     o unsigned int array with ending coordinates
|       for each DI event (get_diCoords fun01)
|   - delSizeAryUI:
|     o unsigned int array with the deletion size for each
|       DI event
|   - indelArySC:
|     o holds indels types found ("D" = del; "I" = ins)
|   - numDIsSI:
|     o number of DI events in diStartAryUI/diEndAryUI
|   - refStartUI:
|     o frist reference base in alignment
|   - refEndUI:
|     o last reference base in alignment
|   - outFILE:
|     o file to print to
| Output:
|   - Prints:
|     o read id, number of DI events, and the coordinates
|       for each event as a tsv to outFILE
\-------------------------------------------------------*/
void
pDI_diCoords(
   signed char *qryIdStr,      /*read id*/
   signed char *refIdStr,      /*reference id*/
   unsigned int *diStartAryUI, /*starting coordiantes*/
   unsigned int *diEndAryUI,   /*ending coordiantes*/
   unsigned int *delSizeAryUI, /*size of each deletion*/
   signed char *indelArySC,    /*indel types found*/
   signed int numDIsSI,        /*number of DI events*/
   unsigned int refStartUI,    /*first mapped ref base*/
   unsigned int refEndUI,      /*last mapped ref base*/
   void *outFILE               /*filt to print to*/
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
