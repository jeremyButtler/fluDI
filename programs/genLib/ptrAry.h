/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' ptrAry SOF:
'   - functions to maintain and organize an array of
'     string pointers
'   o header:
'     - defined variables and guards
'   o .h st01: strPtrAry
'     - struct holding string pointer array & string sizes
'   o st01: str_ptrAry
'     - struct holding string pointer array and string sizes
'   o fun01: blank_str_ptrAry
'     - blanks a str_ptrAry struct
'   o fun02: init_str_ptrAry
'     - initialize a str_ptrAry struct
'   o fun03: freeStack_str_ptrAry
'     - frees c-string array in a str_ptrAry struct
'   o fun04: freeHeap_str_ptrAry
'     - frees a str_ptrAry struct
'   o fun05: resize_str_ptrAry
'     - resizes arrays in a str_ptrAry struct
'   o fun06: mk_str_ptrAry
'     - make initialized heap allocated str_ptrAry struct
'   o fun07: add_str_ptrAry
'     - adds a string to a str_ptrAry struct
'   o fun09: swap_str_ptrAry
'     - swaps two strings in a str_ptrAry struct
'   o fun10: addSort_str_ptrAry
'     - adds a string to a str_ptrAry in sorted positon
'     - assumes you do not manipulate struct in other ways
'   o fun11: cmp_str_ptrAry
'     - compares two strings in a str_ptrAry struct
'   o fun12: sort_str_ptrAry
'     - sorts string pointer array in str_ptrAry struct
'   o fun13: sort_str_ptrAry
'     - sorts string pointer array in str_ptrAry struct
'       while keeping an unsigned int array in sync
'   o fun14: find_str_ptrAry
'     - search for query in str_ptrAry (must be sorted)
'   o fun15: findNoSort_str_ptrAry
'     - search for query in str_ptrAry (dumb search)
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - defined variables and guards
\-------------------------------------------------------*/

#ifndef STRING_POINTER_ARRAY_H
#define STRING_POINTER_ARRAY_H

#define def_memErr_ptrAry 1

/*-------------------------------------------------------\
| ST01: str_ptrAry
|   - struct holding string pointer array and string sizes
\-------------------------------------------------------*/
typedef struct str_ptrAry
{
   signed char **strAry;     /*string pointer array*/
   unsigned int *lenAryUI;  /*length of each string*/
   unsigned int *sizeAryUI; /*size of each string*/

   signed long sizeSL;     /*maximum strings in ary*/
   signed long lenSL;      /*number filled strings*/
}str_ptrAry;

/*-------------------------------------------------------\
| Fun01: blank_str_ptrAry
|   - blanks a str_ptrAry struct
| Input:
|   - strPtrAry:
|     o pointer to strPtrAry structure to blank
| Output:
|   - Modifies:
|     o all strings to be set to null, all lengths to be
|       set to 0, and lenArySL to 0
\-------------------------------------------------------*/
void
blank_str_ptrAry(
   struct str_ptrAry *strPtrAry
);

/*-------------------------------------------------------\
| Fun02: init_str_ptrAry
|   - initialize a str_ptrAry struct
| Input:
|   - strPtrAry:
|     o pointer to str_ptrAry structure to initializes
| Output:
|   - Modifies:
|     o all values set to null
\-------------------------------------------------------*/
void
init_str_ptrAry(
  struct str_ptrAry *strPtrAry
);

/*-------------------------------------------------------\
| Fun03: freeStack_str_ptrAry
|   - frees c-string array in a str_ptrAry struct
| Input:
|   - strPtrAry:
|     o pointer to str_ptrAry struct with vars to free
| Output:
|   - Frees:
|     o all arrays and then initializes
\-------------------------------------------------------*/
void
freeStack_str_ptrAry(
  struct str_ptrAry *strPtrAry
);

/*-------------------------------------------------------\
| Fun04: freeHeap_str_ptrAry
|   - frees a str_ptrAry struct
| Input:
|   - strPtrAry:
|     o pointer to str_ptrAry struct to free
| Output:
|   - Frees:
|     o strPtrAry (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_str_ptrAry(
  struct str_ptrAry *strPtrAry
);

/*-------------------------------------------------------\
| Fun05: resize_str_ptrAry
|   - resizes arrays in a str_ptrAry struct
| Input:
|   - strAryPtr:
|     o pointer to a str_ptrAry struct to add string to
|   - numElmSL:
|     o number elements to have after resizing
| Output:
|   - Modifies:
|     o resizes all arrays in strAryPtr
|     o updates sizeSL in strAryPtr to numElmSL
|   - Returns:
|     o 0 for no errors
|     o def_memErr_ptrAry for memory errors
\-------------------------------------------------------*/
signed char
resize_str_ptrAry(
   struct str_ptrAry *strAryPtr,
   signed long numElmSL
);

/*-------------------------------------------------------\
| Fun06: mk_str_ptrAry
|   - make an initialized heap allocated str_ptrAry struct
| Input:
|   - numElmSL:
|     o number elements to have
| Output:
|   - Returns:
|     o pointer to str_ptrAry struct
|     o 0 for memory error
\-------------------------------------------------------*/
struct str_ptrAry *
mk_str_ptrAry(
   signed long numElmSL
);

/*-------------------------------------------------------\
| Fun07: add_str_ptrAry
|   - adds a string to a str_ptrAry struct
| Input:
|   - newStr:
|     o string to add to array
|   - strAryPtr:
|     o pointer to a str_ptrAry struct to add string to
|   - indexUL:
|     o index to add string at
| Output:
|   - Modifies:
|     o copys newStr to indexUL in strAryPtr
|   - Returns:
|     o 0 for no errors
|     o def_memErr_ptrAry for memory errors
\-------------------------------------------------------*/
signed char
add_str_ptrAry(
   signed char *newStr,
   struct str_ptrAry *strAryPtr,
   signed long indexSL
);

/*-------------------------------------------------------\
| Fun08: findClose_str_ptraRy
|  - finds location to insert query in string ptr array
|    (must be sorted)
| Input:
|  - strAryPtr:
|    o string array
|  - qryStr:
|    o string to find insert location
| Output:
|  - Returns:
|    o index to insert qryStr at
\-------------------------------------------------------*/
signed long
findClose_str_ptrAry(
   struct str_ptrAry *strAryPtr,
   signed char *qryStr
);

/*-------------------------------------------------------\
| Fun09: swap_str_ptrAry
|   - swaps two strings in a str_ptrAry struct
| Input:
|   - strAryPtr:
|     o pointer to str_ptrAry struct with strings to swap
|   - firstSL:
|     o first index to swap
|   - secSL:
|     o second index to swap
| Output:
|   - Modifies:
|     o arrays in str_ptrAry to have swapped strings
\-------------------------------------------------------*/
void
swap_str_ptrAry(
   struct str_ptrAry *strAryPtr,
   signed long firstSL,
   signed long secSL
);

/*-------------------------------------------------------\
| Fun10: addSort_str_ptrAry
|   - adds a string to a str_ptrAry in sorted positon
|   - assumes you do not manipulate struct in other ways
| Input:
|   - newStr:
|     o string to add to array
|   - strAryPtr:
|     o pointer to str_ptrAry struct to add string to
| Output:
|   - Modifies:
|     o strAryPtr to have newStr at closest index
|   - Returns:
|     o index of added string
|     o -1 for memory errors
\-------------------------------------------------------*/
unsigned long
addSort_str_ptrAry(
   signed char *newStr,          /*string to add*/
   struct str_ptrAry *strAryPtr /*add newStr to*/
);

/*-------------------------------------------------------\
| Fun11: cmp_str_ptrAry
|   - compares two strings in a str_ptrAry struct
| Input:
|   - strAryPtr:
|     o str_ptrAry stuct with strings to compare
|   - qrySL:
|     o index of query to compare
|   - refSL:
|     o index of reference to compare
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query is greater
|     o < 0 if reference is greater
\-------------------------------------------------------*/
signed long
cmp_str_ptrAry(
   str_ptrAry *strAryPtr,
   signed long qrySL,
   signed long refSL
);

/*-------------------------------------------------------\
| Fun12: sort_str_ptrAry
|   - sorts string pointer array in str_ptrAry struct
| Input:
|   - strPtrAry:
|     o pointer to str_ptrAry struct to sort
| Output:
|   - Modifies:
|     o arrays in strPtrAry to be sorted by values in
|       strAry array in strPtrAry
\-------------------------------------------------------*/
void
sort_str_ptrAry(
   struct str_ptrAry *strPtrAry
);

/*-------------------------------------------------------\
| Fun13: sort_str_ptrAry
|   - sorts string pointer array in str_ptrAry struct
|     while keeping an unsigned int array in sync
| Input:
|   - strPtrAry:
|     o pointer to str_ptrAry struct to sort
|   - uiAry:
|     o unsigned int array to keep in sync
| Output:
|   - Modifies:
|     o arrays in strPtrAry to be sorted by values in
|       strAry array in strPtrAry
|     o uiAry to be sorted with strPtrAry
\-------------------------------------------------------*/
void
sortSync_str_ptrAry(
   struct str_ptrAry *strPtrAry,
   unsigned int *uiAry
);

/*-------------------------------------------------------\
| Fun14: find_str_ptrAry
|  - search for query in str_ptrAry (must be sorted)
| Input:
|  - strPtrAry:
|    o pointer to str_ptrAry struct to search for query
|  - qryStr:
|    o string to find
| Output:
|  - Returns:
|    o index of qryStr in strAry
|    o -1 if qryStr is not in strAry
\-------------------------------------------------------*/
signed long
find_str_ptrAry(
   struct str_ptrAry *strPtrAry,
   signed char *qryStr
);

/*-------------------------------------------------------\
| Fun15: findNoSort_str_ptrAry
|  - search for query in str_ptrAry (dumb search)
| Input:
|  - strPtrAry:
|    o pointer to str_ptrAry struct to search for query
|  - qryStr:
|    o string to find
| Output:
|  - Returns:
|    o index of qryStr in strAry
|    o -1 if qryStr is not in strAry
\-------------------------------------------------------*/
signed long
findNoSort_str_ptrAry(
   struct str_ptrAry *strPtrAry,
   signed char *qryStr
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
