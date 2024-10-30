/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' strAry SOF:
'   - functions to maintain and organize an array of
'     strings (max length is 63 characters)
'   o header:
'     - defined variables and guards
'   o fun01: mk_strAry
'     - make a string array
'   o fun02: realloc_strAry
'     - rellocates memory for a string array
'   o fun04: add_strAry
'     - adds a string to a string array
'   o .h fun03: get_strAry
'     - finds pointer of a string in a string array
'   o fun05: swap_strAry
'     - swaps two strings in a string array
'   o fun06: cmp_strAry
'     - compares two strings in a string array
'   o fun07: sort_strAry
'     - sorts a string array from least to greatest; is
'       case sensitive
'   o fun08: sortSync_strAry
'     - sorts a string array from least to greatest, but
'       keeps the unsigned int array in sync with strings
'   o fun09: find_strAry
'     - search for query in string array (must be sorted)
'   o fun10: findNoSort_strAry
'     - search for query in string array (dumb search)
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - defined variables and guards
\-------------------------------------------------------*/

#ifndef STRING_ARRAY_H
#define STRING_ARRAY_H

#define def_lenStr_strAry 64 /*max length of a string*/

/*-------------------------------------------------------\
| Fun01: mk_strAry
|   - make a string array
| Input:
|   - sizeUL:
|     o number of strings to store
| Output:
|   - Returns:
|     o pointer to string array
|     o 0 for memory error
\-------------------------------------------------------*/
signed char *
mk_strAry(
   unsigned long sizeUL
);

/*-------------------------------------------------------\
| Fun02: realloc_strAry
|   - rellocates memory for a string array
| Input:
|   - strAry:
|     o string array to reallocate memory for
|   - sizeUL:
|     o new number of strings to store
| Output:
|   - Returns:
|     o pointer to string array (if succeded)
|       - realloc frees only on success
|     o 0 for memory error
\-------------------------------------------------------*/
signed char *
realloc_strAry(
   signed char *strAry,
   unsigned long sizeUL
);

/*-------------------------------------------------------\
| Fun03: get_strAry
|   - finds pointer of a string in a string array
| Input:
|   - strAry:
|     o string array to get string index from
|   - indexUL:
|     o index of target string
| Output:
|   - Returns:
|     o pointer to string at indexUL
\-------------------------------------------------------*/
#define get_strAry(strAry, indexUL) (strAry + (indexUL * def_lenStr_strAry))

/*-------------------------------------------------------\
| Fun04: add_strAry
|   - adds a string to a string array
| Input:
|   - newStr:
|     o string to add to array; must be 63 char or shorter
|   - strAry:
|     o string array to add string to
|   - indexUL:
|     o index to add string at
| Output:
|   - Modifies:
|     o strAry to have newStr at indexUL
\-------------------------------------------------------*/
void
add_strAry(
   signed char *newStr,
   signed char *strAry,
   unsigned long indexUL
);


/*-------------------------------------------------------\
| Fun05: swap_strAry
|   - swaps two strings in a string array
| Input:
|   - strAry:
|     o string array to sort
|   - firstUL:
|     o first index to swap
|   - secUL:
|     o second index to swap
| Output:
|   - Modifies:
|     o strAry to have strings at firstUL and secUL swaped
\-------------------------------------------------------*/
void
swap_strAry(
   signed char *strAry,
   unsigned long firstUL,
   unsigned long secUL
);

/*-------------------------------------------------------\
| Fun06: cmp_strAry
|   - compares two strings in a string array
| Input:
|   - strAry:
|     o string array with strings to compare
|   - qryUL:
|     o index of query to compare
|   - refUL:
|     o index of reference to compare
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query is greater
|     o < 0 if reference is greater
\-------------------------------------------------------*/
signed long
cmp_strAry(
   signed char *strAry,
   unsigned long qryUL,
   unsigned long refUL
);

/*-------------------------------------------------------\
| Fun07: sort_strAry
|   - sorts a string array from least to greatest; is case
|     sensitive
| Input:
|   - strAry:
|     o string array to sort
|   - lenUL:
|     o length of strAry (index 1)
| Output:
|   - Modifies:
|     o strAry to be sorted
\-------------------------------------------------------*/
void
sort_strAry(
   signed char *strAry,
   unsigned long lenUL
);

/*-------------------------------------------------------\
| Fun08: sortSync_strAry
|   - sorts a string array from least to greatest, but
|     keeps the unsigned int array in sync with strings
| Input:
|   - strAry:
|     o string array to sort
|   - uiAry:
|     o unsigned int array to keep in sync with strAry
|   - lenUL:
|     o length of strAry (index 1)
| Output:
|   - Modifies:
|     o strAry to be sorted
\-------------------------------------------------------*/
void
sortSync_strAry(
   signed char *strAry,
   unsigned int *uiAry,
   unsigned long lenUL
);

/*-------------------------------------------------------\
| Fun09: find_strAry
|  - search for query in string array (must be sorted)
| Input:
|  - strAry:
|    o string array
|  - qryStr:
|    o string to find
|  - lenUL:
|    o length of strAry (index 1)
| Output:
|  - Returns:
|    o index of qryStr in strAry
|    o -1 if qryStr is not in strAry
\-------------------------------------------------------*/
signed long
find_strAry(
   signed char *strAry,
   signed char *qryStr,
   signed long lenSL
);

/*-------------------------------------------------------\
| Fun10: findNoSort_strAry
|  - search for query in string array (dumb search)
| Input:
|  - strAry:
|    o string array
|  - qryStr:
|    o string to find
|  - lenUL:
|    o length of strAry (index 1)
| Output:
|  - Returns:
|    o index of qryStr in strAry
|    o -1 if qryStr is not in strAry
\-------------------------------------------------------*/
signed long
findNoSort_strAry(
   signed char *strAry,
   signed char *qryStr,
   signed long lenSL
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
