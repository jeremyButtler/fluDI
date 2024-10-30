/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' shellSort SOF: Start Of File
'   - holds functions for doing shell sorts on arrays
'   o note:
'     - most of this code is from
'        o Adam Drozdek. 2013. Data Structures and
'          Algorithims in c++. Cengage Leraning. fourth
'          edition. pages 505-508
'   o note: this is long, but is really a set of three
'     functions per data type
'     o xxSearch_shellSort
'       - search for query in a sorted array
'     o xxRange_shellSort
'       - search for start of range in a sorted array
'     o xx_shellSort
'       - sort an array by lowest value first
'   o header:
'     - guards
'   o fun01: ulSearch_shellSort
'     - search for query in unsigned long array
'   o fun02: ulRange_shellSort
'     - searches for range in unsinged long array
'   o fun03: ul_shellSort
'     - sorts a unsigned long array from least to greatest
'   o fun04: uiSearch_shellSort
'     - search for query in unsigned int array
'   o fun05: uiRange_shellSort
'     - searches for range in unsinged int array
'   o fun06: ui_shellSort
'     - sorts a unsigned int array from least to greatest
'   o fun07: us_Search_shellSort
'     - search for query in unsigned short array
'   o fun08: us_Range_shellSort
'     - searches for range in unsinged short array
'   o fun09: us_shellSort
'     - sorts unsigned short array from least to greatest
'   o fun10: uc_Search_shellSort
'     - search for query in unsigned char array
'   o fun11: uc_Range_shellSort
'     - searches for range in unsinged char array
'   o fun12: sl_shellSort
'     - sorts signed long array from least to greatest
'   o fun13: sl_Search_shellSort
'     - search for query in signed long array
'   o fun14: sl_Range_shellSort
'     - searches for range in singed long array
'   o fun15: sl_shellSort
'     - sorts signed long array from least to greatest
'   o fun16: si_Search_shellSort
'     - search for query in signed int array
'   o fun17: si_Range_shellSort
'     - searches for range in singed int array
'   o fun18: si_shellSort
'     - sorts signed int array from least to greatest
'   o fun19: ss_Search_shellSort
'     - search for query in signed short array
'   o fun20: ss_Range_shellSort
'     - searches for range in singed short array
'   o fun21: ss_shellSort
'     - sorts signed short array from least to greatest
'   o fun22: sc_Search_shellSort
'     - search for query in signed char array
'   o fun23: sc_Range_shellSort
'     - searches for range in singed char array
'   o fun24: sc_shellSort
'     - sorts signed char array from least to greatest
'   o fun25: dbl_shellSort
'     - sorts a double array from least to greatest
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards
\-------------------------------------------------------*/

#ifndef GENERAL_SHELL_SHORT_H
#define GENERAL_SHELL_SHORT_H

/*-------------------------------------------------------\
| Fun01: ulSearch_shellSort
|  - search for query in unsigned long array
| Input:
|  - ulAry:
|    o to unsigned int array
|  - qryUL:
|    o number to find
|  - lenArySL:
|    o length of ulAry (index 1)
| Output:
|  - Returns:
|    o index of qryUL in ulAry
|    o -1 if qryUL is not in ulAry
\-------------------------------------------------------*/
signed long
ulSearch_shellSort(
   unsigned long *ulAry,
   unsigned long qryUL,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun02: ulRange_shellSort
|   - searches for range in unsinged long array
| Input:
|   - ulAry:
|     o unsigned long array
|   - startUL:
|     o start of the range to look for
|   - endUL:
|     o end of the range to look for
|   - lenArySL:
|     o length of qrySIPtr (index 1)
| Output:
|   - Returns:
|     o first index in the range
|     o -1 if no value is in the range
\-------------------------------------------------------*/
signed long
ulRange_shellSort(
   unsigned long *ulAry,
   unsigned long startUL,
   unsigned long endUL,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun03: ul_shellSort
|   - sorts a unsigned long array from least to greatest
| Input:
|   - ulAry:
|     o array to sort
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o ulAry to be sorted form least to greatest
\-------------------------------------------------------*/
void
ul_shellSort(
   unsigned long *ulAry,
   unsigned long startUL,
   unsigned long endUL
);

/*-------------------------------------------------------\
| Fun04: uiSearch_shellSort
|  - search for query in unsigned int array
| Input:
|  - uiAry:
|    o unsigned int array
|  - qryUI:
|    o number to find
|  - lenArySL:
|    o length of uiAry (index 1)
| Output:
|  - Returns:
|    o index of qryUI in uiAry
|    o -1 if qryUI is not in uiAry
\-------------------------------------------------------*/
signed long
uiSearch_shellSort(
   unsigned int *uiAry,
   unsigned int qryUI,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun05: uiRange_shellSort
|   - searches for range in unsinged int array
| Input:
|   - uiAry:
|     o unsigned int array
|   - startUI:
|     o start of the range to look for
|   - endUI:
|     o end of the range to look for
|   - lenArySL:
|     o length of qrySIPtr (index 1)
| Output:
|   - Returns:
|     o first index in the range
|     o -1 if no value is in the range
\-------------------------------------------------------*/
signed long
uiRange_shellSort(
   unsigned int *uiAry,
   unsigned int startUI,
   unsigned int endUI,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun06: ui_shellSort
|   - sorts a unsigned int array from least to greatest
| Input:
|   - uiAry:
|     o array to sort
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o uiAry to be sorted form least to greatest
\-------------------------------------------------------*/
void
ui_shellSort(
   unsigned int *uiAry,
   unsigned long startUL,
   unsigned long endUL
);

/*-------------------------------------------------------\
| Fun07: usSearch_shellSort
|  - search for query in unsigned short array
| Input:
|  - usAry:
|    o unsigned short array
|  - qryUS:
|    o number to find
|  - lenArySL:
|    o length of usAry (index 1)
| Output:
|  - Returns:
|    o index of qryUI in usAry
|    o -1 if qryUI is not in usAry
\-------------------------------------------------------*/
signed long
usSearch_shellSort(
   unsigned short *usAry,
   unsigned short qryUS,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun08: usRange_shellSort
|   - searches for range in unsinged short array
| Input:
|   - usAry:
|     o unsigned short array
|   - startUS:
|     o start of the range to look for
|   - endUS:
|     o end of the range to look for
|   - lenArySL:
|     o length of qrySIPtr (index 1)
| Output:
|   - Returns:
|     o first index in the range
|     o -1 if no value is in the range
\-------------------------------------------------------*/
signed long
usRange_shellSort(
   unsigned short *usAry,
   unsigned short startUS,
   unsigned short endUS,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun09: us_shellSort
|   - sorts a unsigned short array from least to greatest
| Input:
|   - usAry:
|     o array to sort
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o usAry to be sorted form least to greatest
\-------------------------------------------------------*/
void
us_shellSort(
   unsigned short *usAry,
   unsigned long startUL,
   unsigned long endUL
);

/*-------------------------------------------------------\
| Fun10: ucSearch_shellSort
|  - search for query in unsigned char array
| Input:
|  - ucAry:
|    o unsigned char array
|  - qryUC:
|    o number to find
|  - lenArySL:
|    o length of ucAry (index 1)
| Output:
|  - Returns:
|    o index of qryUI in ucAry
|    o -1 if qryUI is not in ucAry
\-------------------------------------------------------*/
signed long
ucSearch_shellSort(
   unsigned char *ucAry,
   unsigned char qryUC,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun11: ucRange_shellSort
|   - searches for range in unsinged char array
| Input:
|   - ucAry:
|     o unsigned char array
|   - startUC:
|     o start of the range to look for
|   - endUC:
|     o end of the range to look for
|   - lenArySL:
|     o length of qrySIPtr (index 1)
| Output:
|   - Returns:
|     o first index in the range
|     o -1 if no value is in the range
\-------------------------------------------------------*/
signed long
ucRange_shellSort(
   unsigned char *ucAry,
   unsigned char startUC,
   unsigned char endUC,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun12: uc_shellSort
|   - sorts a unsigned char array from least to greatest
| Input:
|   - ucAry:
|     o array to sort
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o ucAry to be sorted form least to greatest
\-------------------------------------------------------*/
void
uc_shellSort(
   unsigned char *ucAry,
   unsigned long startUL,
   unsigned long endUL
);

/*-------------------------------------------------------\
| Fun13: slSearch_shellSort
|  - search for query in signed long array
| Input:
|  - slAry:
|    o signed long array
|  - qrySL:
|    o number to find
|  - lenArySL:
|    o length of slAry (index 1)
| Output:
|  - Returns:
|    o index of qryUI in slAry
|    o -1 if qryUI is not in slAry
\-------------------------------------------------------*/
signed long
slSearch_shellSort(
   signed long *slAry,
   signed long qrySL,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun14: slRange_shellSort
|   - searches for range in singed long array
| Input:
|   - slAry:
|     o signed long array
|   - startSL:
|     o start of the range to look for
|   - endSL:
|     o end of the range to look for
|   - lenArySL:
|     o length of qrySIPtr (index 1)
| Output:
|   - Returns:
|     o first index in the range
|     o -1 if no value is in the range
\-------------------------------------------------------*/
signed long
slRange_shellSort(
   signed long *slAry,
   signed long startSL,
   signed long endSL,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun15: sl_shellSort
|   - sorts a signed long array from least to greatest
| Input:
|   - slAry:
|     o array to sort
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o slAry to be sorted form least to greatest
\-------------------------------------------------------*/
void
sl_shellSort(
   signed long *slAry,
   unsigned long startUL,
   unsigned long endUL
);

/*-------------------------------------------------------\
| Fun16: siSearch_shellSort
|  - search for query in signed int array
| Input:
|  - siAry:
|    o signed int array
|  - qrySI:
|    o number to find
|  - lenArySL:
|    o length of siAry (index 1)
| Output:
|  - Returns:
|    o index of qrySI in siAry
|    o -1 if qrySI is not in siAry
\-------------------------------------------------------*/
signed long
siSearch_shellSort(
   signed int *siAry,
   signed int qrySI,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun17: siRange_shellSort
|   - searches for range in singed int array
| Input:
|   - siAry:
|     o signed int array
|   - startSI:
|     o start of the range to look for
|   - endSI:
|     o end of the range to look for
|   - lenArySL:
|     o length of qrySIPtr (index 1)
| Output:
|   - Returns:
|     o first index in the range
|     o -1 if no value is in the range
\-------------------------------------------------------*/
signed long
siRange_shellSort(
   signed int *siAry,
   signed int startSI,
   signed int endSI,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun18: si_shellSort
|   - sorts a signed int array from least to greatest
| Input:
|   - siAry:
|     o array to sort
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o siAry to be sorted form least to greatest
\-------------------------------------------------------*/
void
si_shellSort(
   signed int *siAry,
   unsigned long startUL,
   unsigned long endUL
);

/*-------------------------------------------------------\
| Fun19: ssSearch_shellSort
|  - search for query in signed short array
| Input:
|  - ssAry:
|    o signed short array
|  - qrySS:
|    o number to find
|  - lenArySL:
|    o length of ssAry (index 1)
| Output:
|  - Returns:
|    o index of qrySS in ssAry
|    o -1 if qrySS is not in ssAry
\-------------------------------------------------------*/
signed long
ssSearch_shellSort(
   signed short *ssAry,
   signed short qrySS,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun20: ssRange_shellSort
|   - searches for range in unsinged char array
| Input:
|   - ssAry:
|     o signed short array
|   - startSS:
|     o start of the range to look for
|   - endSS:
|     o end of the range to look for
|   - lenArySL:
|     o length of qrySIPtr (index 1)
| Output:
|   - Returns:
|     o first index in the range
|     o -1 if no value is in the range
\-------------------------------------------------------*/
signed long
ssRange_shellSort(
   signed short *ssAry,
   signed short startSS,
   signed short endSS,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun21: ss_shellSort
|   - sorts a signed short array from least to greatest
| Input:
|   - ssAry:
|     o array to sort
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o ssAry to be sorted form least to greatest
\-------------------------------------------------------*/
void
ss_shellSort(
   signed short *ssAry,
   unsigned long startUL,
   unsigned long endUL
);

/*-------------------------------------------------------\
| Fun22: scSearch_shellSort
|  - search for query in signed char array
| Input:
|  - scAry:
|    o signed char array
|  - qrySC:
|    o number to find
|  - lenArySL:
|    o length of scAry (index 1)
| Output:
|  - Returns:
|    o index of qrySC in scAry
|    o -1 if qrySC is not in scAry
\-------------------------------------------------------*/
signed long
scSearch_shellSort(
   signed char *scAry,
   signed char qrySC,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun23: scRange_shellSort
|   - searches for range in unsinged char array
| Input:
|   - scAry:
|     o signed char array
|   - startSC:
|     o start of the range to look for
|   - endSC:
|     o end of the range to look for
|   - lenArySL:
|     o length of qrySIPtr (index 1)
| Output:
|   - Returns:
|     o first index in the range
|     o -1 if no value is in the range
\-------------------------------------------------------*/
signed long
scRange_shellSort(
   signed char *scAry,
   signed char startSC,
   signed char endSC,
   signed long lenArySL
);

/*-------------------------------------------------------\
| Fun24: sc_shellSort
|   - sorts a signed char array from least to greatest
| Input:
|   - scAry:
|     o array to sort
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o scAry to be sorted form least to greatest
\-------------------------------------------------------*/
void
sc_shellSort(
   signed char *scAry,
   unsigned long startUL,
   unsigned long endUL
);

/*-------------------------------------------------------\
| Fun25: dbl_shellSort
|   - sorts a double array from least to greatest
| Input:
|   - dblAry:
|     o array to sort
|   - startUL:
|     o first element to start sorting at
|   - endUL:
|     o last element to sort (index 0)
| Output:
|   - Modifies:
|     o scAry to be sorted form least to greatest
\-------------------------------------------------------*/
void
dbl_shellSort(
   double *dblAry,
   unsigned long startUL,
   unsigned long endUL
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
