/*#######################################################\
# Name: diCoords
#   - has functions to get DI coordinates from a sam
#     file
\#######################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header
'     - included libraries
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
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>
#include "../genBio/samEntry.h"

/*.h files only*/
#include "../genLib/dataTypeShortHand.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   o .h  #include "../genLib/base10str.h"
!   o .h  #include "../genLib/numToStr.h"
!   o .h  #include "../genLib/ulCp.h"
!   o .h  #include "../genLib/strAry.h"
!   o .h  #include "../genBio/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - gets start and ending coordinates for DI events
   '   o fun01 sec01:
   '     - variable declerations
   '   o fun01 sec02:
   '     - allocate memory
   '   o fun01 sec03:
   '     - scan for DI coordinates
   '   o fun01 sec04:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int numDISI = 0;
   unsigned int coordUI = 0;
   unsigned int uiCig = 0;
   unsigned int curNonDIUI = 0;
   unsigned int startUI = 0;    /*start of DI event*/

   unsigned int stopScanUI =
      samSTPtr->readLenUI - minEndNtUI;
      /*point in query to stop scaning for DI's*/
  
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
         ! *diStartAryUI
      || ! *diEndAryUI
      || ! *indelArySC
      || *lenArysUIPtr < samSTPtr->lenCigUI
   ){ /*If: need to resize the arrays*/
      if(*diStartAryUI)
         free(*diStartAryUI);

      *diStartAryUI = 0;

      *lenArysUIPtr = samSTPtr->lenCigUI;

      *diStartAryUI= malloc(*lenArysUIPtr * sizeof(uint));

      if(! *diStartAryUI)
         goto memErr_sec04;


      if(*diEndAryUI)
         free(*diEndAryUI);

      *diEndAryUI = 0;

      *diEndAryUI = malloc(*lenArysUIPtr * sizeof(uint));

      if(! *diEndAryUI)
         goto memErr_sec04;


      if(*delSizeAryUI)
         free(*delSizeAryUI);

      *delSizeAryUI = 0;

      *delSizeAryUI =
         malloc(*lenArysUIPtr * sizeof(uint));

      if(! *delSizeAryUI)
         goto memErr_sec04;


     *indelArySC = 0;
     *indelArySC = malloc(*lenArysUIPtr * sizeof(schar));

     if(! *indelArySC)
        goto memErr_sec04;
   } /*If: need to resize the arrays*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - scan for DI coordinates
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(samSTPtr->cigTypeStr[samSTPtr->lenCigUI - 1] == 'S')
      stopScanUI -=
         samSTPtr->cigArySI[samSTPtr->lenCigUI - 1];
         /*account for softmasking at end*/

   if(samSTPtr->cigTypeStr[0] == 'S')
      stopScanUI -= samSTPtr->cigArySI[0];
      /*account for softmasking at start*/
      
   coordUI = samSTPtr->refStartUI;


   for(
      uiCig = 0;
      uiCig < samSTPtr->lenCigUI;
      ++uiCig
   ){ /*Loop: find di events*/ 

      if(samSTPtr->cigTypeStr[uiCig] == 'S')
         continue;

      else if(samSTPtr->cigTypeStr[uiCig] == 'H')
         continue;

      else if( samSTPtr->cigTypeStr[uiCig] == 'D' )
      { /*If: have deletion*/
         startUI = coordUI; /*start of deletion*/
         coordUI += samSTPtr->cigArySI[uiCig];

         if(! (indelFlagSC & 1) )
            continue; /*ignoring deletions*/
      } /*If: have deletion*/
  
      else if(samSTPtr->cigTypeStr[uiCig] == 'I')
      { /*If: have insertion*/
         startUI = coordUI;

         if(curNonDIUI < minEndNtUI)
         { /*If: to early for recording DI*/
            curNonDIUI += samSTPtr->cigArySI[uiCig];
            continue; /*ignoring insertions*/
         } /*If: to early for recording DI*/

         curNonDIUI += samSTPtr->cigArySI[uiCig];

         if(! (indelFlagSC & 2) )
            continue; /*ignoring insertions*/
      } /*If: have insertion*/

      else
      { /*If: not indel*/
         curNonDIUI += samSTPtr->cigArySI[uiCig];
         coordUI += samSTPtr->cigArySI[uiCig];
         continue;
      } /*If: not indel*/


      if(curNonDIUI >= stopScanUI)
         break; /*finished*/

      else if(curNonDIUI < minEndNtUI)
         continue; /*to early for DI event*/

      else if(
         samSTPtr->cigArySI[uiCig] < (sint) minDILenUI
      ) continue; /*DI event to small*/


      (*delSizeAryUI)[numDISI]= samSTPtr->cigArySI[uiCig];
      (*indelArySC)[numDISI]= samSTPtr->cigTypeStr[uiCig];
      (*diStartAryUI)[numDISI] = startUI;/*first DI base*/

      if(samSTPtr->cigTypeStr[uiCig] == 'D')
         (*diEndAryUI)[numDISI] = coordUI - 1;/*end del*/
      else
         (*diEndAryUI)[numDISI] = coordUI + 1;
         /*base after insertion*/

      ++numDISI; /*number DI events*/

   } /*Loop: find di events*/ 

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto cleanUp_sec04;

   memErr_sec04:;
      numDISI = -1;
      goto cleanUp_sec04;

   cleanUp_sec04:;
      return(numDISI);
} /*get_diCoords*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - gets start and ending coordinates for DI events
   '   o fun02 sec02:
   '     - variable declerations
   '   o fun02 sec02:
   '     - scan for DI coordinates
   '   o fun02 sec03:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   sint numDISI = 0;
   uint coordUI = 0;
   uint uiCig = 0;

   unsigned int stopScanUI =
      samSTPtr->readLenUI - minEndNtUI;
      /*point in query to stop scaning for DI's*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - scan for DI coordinates
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(samSTPtr->cigTypeStr[samSTPtr->lenCigUI - 1] == 'S')
      stopScanUI -=
         samSTPtr->cigArySI[samSTPtr->lenCigUI - 1];
         /*account for softmasking at end*/

   if(samSTPtr->cigTypeStr[0] == 'S')
      stopScanUI -= samSTPtr->cigArySI[0];
      /*account for softmasking at start*/

   coordUI = samSTPtr->refStartUI;

   for(
      uiCig = 0;
      uiCig < samSTPtr->lenCigUI;
      ++uiCig
   ){ /*Loop: find di events*/ 
      if(samSTPtr->cigTypeStr[uiCig] == 'S')
         continue; /*soft masked, do not care*/

      else if(samSTPtr->cigTypeStr[uiCig] != 'D')
      { /*Else If: not a deletion*/
         coordUI += samSTPtr->cigArySI[uiCig];
         continue;
      } /*Else If: not a deletion*/

      else if(
         samSTPtr->cigArySI[uiCig] < (sint) minDILenUI
      ){ /*Else If: deletion is to small*/
         coordUI += samSTPtr->cigArySI[uiCig];
         continue; /*if not a deletion*/
      }  /*Else If: deletion is to small*/

      else if(coordUI >= stopScanUI)
         break; /*finished (end of query scan area)*/

      ++numDISI; /*number DI events*/
   } /*Loop: find di events*/ 

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return(numDISI);
} /*scan_diCoords*/

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
){
   fprintf(
      (FILE *) outFILE,
      "read\tref\tstart_di\tend_di\tnum_del\ttype"
   );

   fprintf(
      (FILE *) outFILE,
      "\tref_start\tref_end\tnum_di\n"
   );
} /*pDIHead_diCoords*/

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
){
   sint siDI = 0;

   for(
      siDI = 0;
      siDI < numDIsSI;
      ++siDI
   ){ /*Loop: print out coordinates*/
      fprintf(
         (FILE *) outFILE,
         "%s\t%s\t%u\t%u\t%u\t%c\t%u\t%u\t%i\n",
         qryIdStr,
         refIdStr,
         diStartAryUI[siDI] + 1,
         diEndAryUI[siDI] + 1,
         delSizeAryUI[siDI],
         indelArySC[siDI],
         refStartUI + 1,
         refEndUI + 1,
         numDIsSI
      );
   } /*Loop: print out coordinates*/

   if(! siDI)
   { /*If: no DI events*/
      fprintf(
         (FILE *) outFILE,
         "%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\t0\n",
         qryIdStr,
         refIdStr
      );
   } /*If: no DI events*/
} /*pDI_diCoords*/

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
