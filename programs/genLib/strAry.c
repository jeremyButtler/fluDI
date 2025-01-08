/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' strAry SOF:
'   - functions to maintain and organize an array of
'     strings (max length is 63 characters)
'   o header:
'     - included libraries
'   o fun01: mk_strAry
'     - make a string array
'   o fun02: realloc_strAry
'     - rellocates memory for a string array
'   o .h fun03: get_strAry
'     - finds pointer of a string in a string array
'   o fun04: add_strAry
'     - adds a string to a string array
'   o fun05: swap_strAry
'     - swaps two strings in a string array
'   o fun06: cmp_strAry
'     - compares a string to a string array
'   o fun07: cmpIndex_strAry
'     - compares two strings in a string array
'   o fun08: findInsert_strAry
'     - finds location to insert query in string array
'       (must be sorted)
'   o fun09: addSort_strAry
'     - adds a string to string array in a sorted positon
'   o fun10: sort_strAry
'     - sorts a string array from least to greatest; is
'       case sensitive
'   o fun11: sortSync_strAry
'     - sorts a string array from least to greatest, but
'       keeps the unsigned int array in sync with strings
'   o fun12: find_strAry
'     - search for query in string array (must be sorted)
'   o fun13: findNoSort_strAry
'     - search for query in string array (dumb search)
'   o license:
'     - licensing for this code (public domain / mit)
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

#include "strAry.h"
#include "ulCp.h"

/*-------------------------------------------------------\
| Fun01: mk_strAry
|   - make a string array
| Input:
<|   - sizeUL:
|     o number of strings to store
| Output:
|   - Returns:
|     o pointer to string array
|     o 0 for memory error
\-------------------------------------------------------*/
signed char *
mk_strAry(
   unsigned long sizeUL
){
   return
     calloc(
         sizeUL * def_lenStr_strAry,
         sizeof(signed char)
      );
} /*mk_strAry*/

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
){
   return
      realloc(
         strAry,
         sizeUL * def_lenStr_strAry * sizeof(signed char)
      );
} /*realloc_strAry*/

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
){
   signed char *tmpStr = get_strAry(strAry, indexUL);

   tmpStr +=
      cpStr_ulCp(
         tmpStr,
         newStr
      );
} /*add_strAry*/

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
){
   swapDelim_ulCp(
      get_strAry(strAry, firstUL),
      get_strAry(strAry, secUL),
      0,
      '\0'
   ); /*swap strings*/
} /*swap_strAry*/

/*-------------------------------------------------------\
| Fun06: cmp_strAry
|   - compares a string to a string array
| Input:
|   - qryStr:
|     o c-string with query to compare
|   - strAry:
|     o c-string array with strings to compare
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
   signed char *qryStr,
   signed char *strAry,
   unsigned long refUL
){
   return
      eqlNull_ulCp(
         qryStr,
         get_strAry(strAry, refUL)
      ); /*compare strings*/
} /*cmp_strAry*/


/*-------------------------------------------------------\
| Fun07: cmpIndex_strAry
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
cmpIndex_strAry(
   signed char *strAry,
   unsigned long qryUL,
   unsigned long refUL
){
   return
      eqlNull_ulCp(
         get_strAry(strAry, qryUL),
         get_strAry(strAry, refUL)
      ); /*compare strings*/
} /*cmpIndex_strAry*/

/*-------------------------------------------------------\
| Fun08: findInsert_strAry
|  - finds location to insert query in string array
|    (must be sorted)
| Input:
|  - strAry:
|    o string array
|  - qryStr:
|    o string to find insert location
|  - lenUL:
|    o length of strAry (index 1)
| Output:
|  - Returns:
|    o index to insert qryStr at
\-------------------------------------------------------*/
signed long
findClose_strAry(
   signed char *strAry,
   signed char *qryStr,
   signed long lenSL
){
   signed long midSL = 0;
   signed long rightSL = lenSL - 1;
   signed long leftSL = 0;
   signed long cmpSL = 0;

   while(leftSL <= rightSL)
   { /*Loop: Search for the querys index*/
      midSL = (leftSL + rightSL) >> 1;

      cmpSL =
         eql_ulCp(
            qryStr,
            get_strAry(strAry, midSL),
            0,
            '\0'
         ); /*compare query to array value*/

      if(cmpSL > 0)
          leftSL = midSL + 1;

      else if(cmpSL < 0)
          rightSL = midSL - 1;

      else
      { /*Else: found matching string*/
         if(midSL + 1 >= rightSL)
            return midSL;

         cmpSL =
            eql_ulCp(
               qryStr,
               get_strAry(strAry, midSL + 1),
               0,
               '\0'
            ); /*compare query to array value*/

         if(! cmpSL)
            leftSL = midSL + 1; /*find end of duplicates*/
         else
            return midSL;
      } /*Else: found matching string*/
   } /*Loop: Search for the querys index*/

   if(midSL == lenSL)
      return midSL;     /*at end of array, insert here*/
   if(cmpSL > 0)
      return midSL + 1; /*query still greater*/
   else
      return midSL;     /*query less than, insert here*/
} /*findClose_strAry*/

/*-------------------------------------------------------\
| Fun09: addSort_strAry
|   - adds a string to a string array in a sorted positon
|   - assumes you do not manipulate array in other ways
| Input:
|   - newStr:
|     o string to add to array; must be 63 char or shorter
|   - strAry:
|     o string array to add string to
|   - lenUL:
|     o has number of elements in array
| Output:
|   - Modifies:
|     o strAry to have newStr at closest index
|   - Returns:
|     o index of added string
\-------------------------------------------------------*/
unsigned long
addSort_strAry(
   signed char *newStr, /*string to add to array*/
   signed char *strAry, /*string array to add newStr to*/
   unsigned long lenUL  /*number of elements in array*/
){
   unsigned long indexUL = 0;
   signed long diffSL = 0;


   if(! lenUL)
   { /*If: first string in array*/
      cpStr_ulCp(
         strAry,
         newStr
      );

      indexUL = 0;
      goto ret_fun09;
   } /*If: first string in array*/


   indexUL =
      (unsigned long)
      findClose_strAry(
         strAry,
         newStr,
         (signed long) lenUL
      ); /*find insert location*/

   diffSL = (signed long) lenUL;

   if(diffSL > (signed long) indexUL)
   { /*If: can copy last value*/
      cpStr_ulCp(
         get_strAry(strAry, diffSL),
         get_strAry(strAry, diffSL - 1)
      ); /*prevents copying unitialized last value*/

      --diffSL;
   } /*If: can copy last value*/

   while(diffSL > (signed long) indexUL)
   { /*Loop: move strings up (clear position)*/
      swap_strAry(
         strAry,
         diffSL,
         diffSL - 1
      );

      --diffSL;
   } /*Loop: move strings up (clear position)*/
   

   add_strAry(
      newStr,
      strAry,
      indexUL
   ); /*add new string in*/

   ret_fun09:;
      return indexUL;
} /*addSort_strAry*/

/*-------------------------------------------------------\
| Fun10: sort_strAry
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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun10 TOC:
   '   - sorts a string array from least to greatest; is
   '     case sensitive
   '   o fun10 sec01:
   '     - variable declerations
   '   o fun10 sec02:
   '     - find the number of rounds to sort for
   '   o fun10 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of sorting rounds*/
   unsigned long subUL = 0;
   unsigned long nextUL = 0;
   unsigned long lastUL = 0;
   unsigned long onUL = 0;

   /*Variables to incurment loops*/
   unsigned long ulIndex = 0;
   unsigned long ulElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Recursion formuia: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subUL = 1; /*Initialzie first array*/

   while(subUL < lenUL - 1)
      subUL = (3 * subUL) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec03:
   ^   - sort arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(subUL > 0)
   { /*Loop: all rounds*/
      for(
         ulIndex = 0;
         ulIndex <= subUL;
         ++ulIndex
      ){ /*Loop: though sub array*/
         ulElm = ulIndex;

         for(
            ulElm = ulIndex;
            ulElm + subUL < lenUL;
            ulElm += subUL
         ){ /*Loop: swap elements in subarray*/
            nextUL = ulElm + subUL;

            if(
               cmpIndex_strAry(
                  strAry,
                  ulElm,
                  nextUL
               ) > 0
            ){ /*If I need to swap an element*/
               swap_strAry(
                  strAry,
                  ulElm,
                  nextUL
               );

               lastUL = ulElm;
               onUL = ulElm;

               while(lastUL >= subUL)
               { /*loop: move swapped element back*/
                  lastUL -= subUL;

                  if(
                     cmpIndex_strAry(
                        strAry,
                        onUL,
                        lastUL
                     ) > 0
                  ) break; /*moved onUL to its position*/

                  swap_strAry(
                     strAry,
                     onUL,
                     lastUL
                  );

                  onUL = lastUL;
               } /*Loop: move swapped element back*/
            } /*If I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subUL = (subUL - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*sort_strAry*/

/*-------------------------------------------------------\
| Fun11: sortSync_strAry
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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC:
   '   - sorts a string array from least to greatest, but
   '     keeps the unsigned int array in sync with strings
   '   o fun11 sec01:
   '     - variable declerations
   '   o fun11 sec02:
   '     - find the number of rounds to sort for
   '   o fun11 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of sorting rounds*/
   unsigned long subUL = 0;
   unsigned long nextUL = 0;
   unsigned long lastUL = 0;
   unsigned long onUL = 0;

   /*Variables to incurment loops*/
   unsigned long ulIndex = 0;
   unsigned long ulElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(lenUL <= 1)
      return; /*nothing to do*/

   /*Recursion formuia: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subUL = 1; /*Initialzie first array*/

   while(subUL < lenUL - 1)
      subUL = (3 * subUL) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec03:
   ^   - sort arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(subUL > 0)
   { /*Loop: all rounds*/
      for(
         ulIndex = 0;
         ulIndex <= subUL;
         ++ulIndex
      ){ /*Loop: though sub array*/
         ulElm = ulIndex;

         for(
            ulElm = ulIndex;
            ulElm + subUL < lenUL;
            ulElm += subUL
         ){ /*Loop: swap elements in subarray*/
            nextUL = ulElm + subUL;

            if(
               cmpIndex_strAry(
                  strAry,
                  ulElm,
                  nextUL
               ) > 0
            ){ /*If I need to swap an element*/
               swap_strAry(
                  strAry,
                  ulElm,
                  nextUL
               );

               uiAry[ulElm] ^= uiAry[nextUL];
               uiAry[nextUL] ^= uiAry[ulElm];
               uiAry[ulElm] ^= uiAry[nextUL];

               lastUL = ulElm;
               onUL = ulElm;

               while(lastUL >= subUL)
               { /*loop: move swapped element back*/
                  lastUL -= subUL;

                  if(
                     cmpIndex_strAry(
                        strAry,
                        onUL,
                        lastUL
                     ) > 0
                  ) break; /*moved onUL to its position*/

                  swap_strAry(
                     strAry,
                     onUL,
                     lastUL
                  );

                  uiAry[onUL] ^= uiAry[lastUL];
                  uiAry[lastUL] ^= uiAry[onUL];
                  uiAry[onUL] ^= uiAry[lastUL];

                  onUL = lastUL;
               } /*Loop: move swapped element back*/
            } /*If I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subUL = (subUL - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*sortSync_strAry*/

/*-------------------------------------------------------\
| Fun12: find_strAry
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
){
   signed long midSL = 0;
   signed long rightSL = lenSL - 1;
   signed long leftSL = 0;
   signed long cmpSL = 0;

   while(leftSL <= rightSL)
   { /*Loop: Search for the querys index*/
      midSL = (leftSL + rightSL) >> 1;

      cmpSL =
         eql_ulCp(
            qryStr,
            get_strAry(strAry, midSL),
            0,
            '\0'
         ); /*compare query to array value*/

      if(cmpSL > 0)
          leftSL = midSL + 1;

      else if(cmpSL < 0)
          rightSL = midSL - 1;

      else
         return midSL;
   } /*Loop: Search for the querys index*/

   return -1; /*query not found*/
} /*find_strAry*/

/*-------------------------------------------------------\
| Fun13: findNoSort_strAry
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
){
   signed long indexSL = 0;

   while(indexSL < lenSL)
   { /*Loop: Search for the querys index*/
      if(
        ! eql_ulCp(
            qryStr,
            get_strAry(strAry, indexSL),
            0,
            '\0'
         ) /*compare query to array value*/
     ) return indexSL;
 
     ++indexSL;
   } /*Loop: Search for the querys index*/

   return -1; /*query not found*/
} /*findNoSort_strAry*/

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
