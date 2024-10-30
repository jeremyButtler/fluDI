/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' adjCoords SOF: Start Of File
'   - holds functions for adjusting sam file coordinates
'     from an gene mapping to an reference mapping
'   o header:
'     - included libraries
'   o fun01: adjCoords
'     - adjust coordinates in sam file to reflect the 
'       actual reference position
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

#include "adjCoords.h"

#include <stdio.h>

#include "geneCoord.h"
#include "samEntry.h"
#include "../genLib/ulCp.h"

/*These have no .c files*/
#include "../genLib/dataTypeShortHand.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden includes:
!   - .c  #include "../genLib/numToStr.h"
!   - .c  #include "../genLib/base10StrToNum.h"
!   - .c  #include "../genLib/charCp.h"
!   - .c  #include "../genLib/strAry.h"
!   - .h  #include "../genLib/genMath.h" .h macro only
!   - .h  #include "ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: adjCoords
|   - addjust coordinates in sam file to reflect the 
|     actual reference position
| Input:
|   - samSTPtr:
|     o pointer to samEntry structure with the read to
|       adjust coordinates for
|   - coordsSTPtr:
|     o pointer to an geneCoord structure with the
|       coordinates for each gene
|   - numGenesSI:
|     o number of genes in coordsSIPtr
| Output:
|   - Returns:
|     o 0 for no problems
|     o 1 if gene is not in the list
\-------------------------------------------------------*/
signed char
adjCoords(
   struct samEntry *samSTPtr,
   struct geneCoord *coordsSTPtr,
   signed int numGenesSI
){
   schar errSC = 0;
   sint siIndx = 0;

   siIndx =
      findName_geneCoord(
         coordsSTPtr,
         (schar *) samSTPtr->refIdStr,
         numGenesSI
      );

   if(siIndx == -1)
      goto noGene_fun01;

   /*sam files will have gene if reference orientaiton,
   '  while the coordinates will be in reading frame
   */
   if(coordsSTPtr->dirAryUC[siIndx])
   { /*If: this was an reverse complement gene*/
      revCmp_samEntry(samSTPtr);
         
      samSTPtr->refStartUI = 
           coordsSTPtr->endAryUI[siIndx]
         - samSTPtr->refStartUI
         - samSTPtr->alnReadLenUI;

      /*Account for subtracting an index 1 value*/
      ++samSTPtr->refStartUI;

      samSTPtr->refEndUI = samSTPtr->refStartUI;
      samSTPtr->refEndUI += samSTPtr->alnReadLenUI;

      /*Account for adding an index one value*/
      --samSTPtr->refEndUI;
   } /*If: this was an reverse complement gene*/

   else
   { /*Else: this is a foward gene*/
      samSTPtr->refStartUI +=
         coordsSTPtr->startAryUI[siIndx];

      samSTPtr->refEndUI +=
         coordsSTPtr->startAryUI[siIndx];
   } /*Else: this is a foward gene*/

   samSTPtr->lenRefIdUC =
      (uchar)
      cpDelim_ulCp(
         samSTPtr->refIdStr,
         coordsSTPtr->refAryStr[siIndx],
         0,
         '\0'
      );

   errSC = 0;
   goto ret_fun01;

   noGene_fun01:;
   errSC = 1;
   goto ret_fun01;

   ret_fun01:;
   return errSC;
} /*adjCoords*/

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
