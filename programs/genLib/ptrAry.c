/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' ptrAry SOF:
'   - functions to maintain and organize an array of
'     string pointers
'   o header:
'     - included libraries
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
'   o fun08: findClose_str_ptrAry
'    - finds location to insert query in string ptr array
'      (must be sorted)
'   o fun09: swap_str_ptrAry
'     - swaps two strings in a str_ptrAry struct
'   o fun10: addSort_str_ptrAry
'     - adds a string to a str_ptrAry in sorted positon
'     - assumes you do not manipulate struct in other ways
'   o fun11: cmp_str_ptrAry
'     - compares two strings in a str_ptrAry struct
'   o fun12: sort_str_ptrAry
'     - sorts string pointer array in str_ptrAry struct
'   o fun13: sortSync_str_ptrAry
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
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "ptrAry.h"
#include "ulCp.h"

/*-------------------------------------------------------\
| Fun01: blank_str_ptrAry
|   - blanks a str_ptrAry struct
| Input:
|   - strPtrAry:
|     o pointer to strPtrAry structure to blank
| Output:
|   - Modifies:
|     o all strings to be set to null, all lengths to be
|       set to 0, and lenAryUI to 0
\-------------------------------------------------------*/
void
blank_str_ptrAry(
   struct str_ptrAry *strPtrAry
){
   signed long slPos = 0;

   for(
      slPos = 0;
      slPos < strPtrAry->lenSL;
      ++slPos
   ){ /*Loop: blank the string array*/
      strPtrAry->strAry[slPos][0] = '\0';
      strPtrAry->lenAryUI[slPos] = 0;
      strPtrAry->sizeAryUI[slPos] = 0;
   } /*Loop: blank the string array*/

   strPtrAry->lenSL = 0;
} /*blank_str_ptrAry*/

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
){
   strPtrAry->strAry = 0;
   strPtrAry->lenAryUI = 0;
   strPtrAry->sizeAryUI = 0;

   strPtrAry->sizeSL = 0;
   strPtrAry->lenSL = 0;

   blank_str_ptrAry(strPtrAry);
} /*init_str_ptrAry*/

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
){
   signed long slPos = 0;

   for(
      slPos = 0;
      slPos < strPtrAry->sizeSL;
      ++slPos
   ){ /*Loop: free all strings*/
      if(strPtrAry->strAry[slPos])
         free(strPtrAry->strAry[slPos]);
      strPtrAry->strAry[slPos] = 0;
   } /*Loop: free all strings*/

   free(strPtrAry->strAry);
   free(strPtrAry->lenAryUI);
   free(strPtrAry->sizeAryUI);

   init_str_ptrAry(strPtrAry);
} /*freeStack_str_ptrAry*/

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
){
   if(! strPtrAry)
      return ; /*nothing to free*/

   freeStack_str_ptrAry(strPtrAry);
   free(strPtrAry);
} /*freeHeap_str_ptrAry*/

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
){
   unsigned int *tmpUIPtr = 0;
   signed char **tmpAryStr = 0;
   signed long oldSizeSL = strAryPtr->sizeSL;

   if(strAryPtr->strAry)
      tmpAryStr =
         realloc(
            strAryPtr->strAry,
            numElmSL * sizeof(signed char *)
         );

   else
      tmpAryStr= malloc(numElmSL * sizeof(signed char *));

   if(! tmpAryStr)
      return def_memErr_ptrAry;

   strAryPtr->strAry = tmpAryStr;
   strAryPtr->sizeSL = numElmSL;


   if(strAryPtr->lenAryUI)
      tmpUIPtr =
         realloc(
            strAryPtr->lenAryUI,
            numElmSL * sizeof(unsigned int)
         );
   else
      tmpUIPtr = malloc(numElmSL * sizeof(unsigned int));

   if(! tmpUIPtr)
      return def_memErr_ptrAry;

   strAryPtr->lenAryUI = tmpUIPtr;


   if(strAryPtr->sizeAryUI)
      tmpUIPtr =
         realloc(
            strAryPtr->sizeAryUI,
            numElmSL * sizeof(unsigned int)
         );
   else
      tmpUIPtr = malloc(numElmSL * sizeof(unsigned int));

   if(! tmpUIPtr)
      return def_memErr_ptrAry;

   strAryPtr->sizeAryUI = tmpUIPtr;

   for(
      numElmSL = oldSizeSL;
      numElmSL < strAryPtr->sizeSL;
      ++numElmSL
   ){ /*Loop: blank sizes*/
      strAryPtr->sizeAryUI[numElmSL] = 0;
      strAryPtr->lenAryUI[numElmSL] = 0;
      strAryPtr->strAry[numElmSL] = 0;
   } /*Loop: blank sizes*/

   return 0;
} /*resize_str_ptrAry*/

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
){
   struct str_ptrAry *retHeapST = 0;
   retHeapST = malloc(sizeof(struct str_ptrAry));

   if(! retHeapST)
      goto memErr_fun05;

   init_str_ptrAry(retHeapST);

   if(
      resize_str_ptrAry(
         retHeapST,
         numElmSL
      )
   ) goto memErr_fun05;

   goto ret_fun05;

   memErr_fun05:;
      if(retHeapST)
         freeHeap_str_ptrAry(retHeapST);
      retHeapST = 0;
      goto ret_fun05;

   ret_fun05:;
      return retHeapST;
} /*mk_str_ptrAry*/


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
){
   unsigned int lenUI = 0;

   if(indexSL >= strAryPtr->sizeSL)
   { /*If: need to expand the array*/
      lenUI =
         (unsigned int)
         resize_str_ptrAry(
            strAryPtr,
            indexSL + 16
         );

      if(lenUI)
         return def_memErr_ptrAry;
   } /*If: need to expand the array*/


   lenUI = endStr_ulCp(newStr); /*find length of string*/


   if(! strAryPtr->sizeAryUI[indexSL])
   { /*If: need to allocate memory*/
      strAryPtr->strAry[indexSL] =
         malloc((lenUI + 8) * sizeof(signed char));

      if(! strAryPtr->strAry[indexSL])
         return def_memErr_ptrAry;

      strAryPtr->sizeAryUI[indexSL] = lenUI;
   } /*If: need to allocate memory*/

   if(strAryPtr->sizeAryUI[indexSL] < lenUI)
   { /*If: need to resize array*/
       free(strAryPtr->strAry[indexSL]);
       strAryPtr->strAry[indexSL] = 0;

       strAryPtr->strAry[indexSL] =
         malloc((lenUI + 8) * sizeof(signed char));

      if(! strAryPtr->strAry[indexSL])
         return def_memErr_ptrAry;

      strAryPtr->sizeAryUI[indexSL] = lenUI;
   } /*If: need to resize array*/

   cpLen_ulCp(
      strAryPtr->strAry[indexSL],
      newStr,
      lenUI
   ); /*copy string*/

   strAryPtr->lenAryUI[indexSL] = lenUI;

   if(indexSL >= strAryPtr->lenSL)
      strAryPtr->lenSL = indexSL + 1;

   return 0;
} /*add_str_ptrAry*/

/*-------------------------------------------------------\
| Fun08: findClose_str_ptraRy
|  - finds location to insert query in string ptr array
|    (must be sorted)
| Input:
|  - strAryPtr:
|    o pointer to str_ptrAry to find query in
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
){
   signed long midSL = 0;
   signed long rightSL = strAryPtr->lenSL - 1;
   signed long leftSL = 0;
   signed long cmpSL = 0;

   while(leftSL <= rightSL)
   { /*Loop: Search for the querys index*/
      midSL = (leftSL + rightSL) >> 1;

      cmpSL =
         eqlNull_ulCp(
            qryStr,
            strAryPtr->strAry[midSL]
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
            eqlNull_ulCp(
               qryStr,
               strAryPtr->strAry[midSL + 1]
            ); /*compare query to array value*/

         if(! cmpSL)
            leftSL = midSL + 1; /*find end of duplicates*/
         else
            return midSL;
      } /*Else: found matching string*/
   } /*Loop: Search for the querys index*/

   if(midSL == strAryPtr->lenSL)
      return midSL;     /*at end of array, insert here*/
   if(cmpSL > 0)
      return midSL + 1; /*query still greater*/
   else
      return midSL;     /*query less than, insert here*/
} /*findClose_strAry*/

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
){
   signed char *tmpStr = 0;

   tmpStr = strAryPtr->strAry[firstSL];
   strAryPtr->strAry[firstSL] = strAryPtr->strAry[secSL];
   strAryPtr->strAry[secSL] = tmpStr;

   strAryPtr->lenAryUI[firstSL] ^=
      strAryPtr->lenAryUI[secSL];
   strAryPtr->lenAryUI[secSL] ^=
      strAryPtr->lenAryUI[firstSL];
   strAryPtr->lenAryUI[firstSL] ^=
      strAryPtr->lenAryUI[secSL];

   strAryPtr->sizeAryUI[firstSL] ^=
      strAryPtr->sizeAryUI[secSL];
   strAryPtr->sizeAryUI[secSL] ^=
      strAryPtr->sizeAryUI[firstSL];
   strAryPtr->sizeAryUI[firstSL] ^=
      strAryPtr->sizeAryUI[secSL];
} /*swap_str_ptrAry*/

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
){
   signed long indexSL = 0;
   signed long diffSL = 0;

   if(! strAryPtr->lenSL)
   { /*If: first string in array*/
      if(
         add_str_ptrAry(
            newStr,
            strAryPtr,
            0
         ) /*add new string in*/
      ) return -1; /*memory error*/

      indexSL = 0;
      ++strAryPtr->lenSL;
      goto ret_fun10;
   } /*If: first string in array*/


   indexSL =
      findClose_str_ptrAry(
         strAryPtr,
         newStr
      ); /*find insert location*/

   diffSL = strAryPtr->lenSL;

   if(diffSL > indexSL)
   { /*If: can copy last value*/
      strAryPtr->strAry[diffSL] = 0;
      strAryPtr->lenAryUI[diffSL] = 0;
      strAryPtr->sizeAryUI[diffSL] = 0;

      /*prevents copying unitialized last value*/
   } /*If: can copy last value*/

   while(diffSL > indexSL)
   { /*Loop: move strings up (clear position)*/
      swap_str_ptrAry(
         strAryPtr,
         diffSL,
         diffSL - 1
      );

      --diffSL;
   } /*Loop: move strings up (clear position)*/
   

   if(
      add_str_ptrAry(
         newStr,
         strAryPtr,
         (unsigned long) indexSL
      ) /*add new string in*/
   ) return -1; /*memory error*/

   ret_fun10:;
      return (unsigned long) indexSL;
} /*addSort_strAry*/

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
){
   return
      eqlNull_ulCp(
         strAryPtr->strAry[qrySL],
         strAryPtr->strAry[refSL]
      ); /*compare strings*/
} /*cmp_str_ptrAry*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun12 TOC:
   '   - sorts string pointer array in str_ptrAry struct
   '   o fun12 sec01:
   '     - variable declerations
   '   o fun12 sec02:
   '     - find the number of rounds to sort for
   '   o fun12 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of sorting rounds*/
   signed long subSL = 0;
   signed long nextSL = 0;
   signed long lastSL = 0;
   signed long onSL = 0;
   signed long lenSL = strPtrAry->lenSL;

   /*Variables to incurment loops*/
   signed long slIndex = 0;
   signed long slElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(lenSL <= 1)
      return; /*nothing to do*/

   /*Recursion formuia: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subSL = 1; /*Initialzie first array*/

   while(subSL < lenSL - 1)
      subSL = (3 * subSL) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec03:
   ^   - sort arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(subSL > 0)
   { /*Loop: all rounds*/
      for(
         slIndex = 0;
         slIndex <= subSL;
         ++slIndex
      ){ /*Loop: though sub array*/
         slElm = slIndex;

         for(
            slElm = slIndex;
            slElm + subSL < lenSL;
            slElm += subSL
         ){ /*Loop: swap elements in subarray*/
            nextSL = slElm + subSL;

            if(
               cmp_str_ptrAry(
                  strPtrAry,
                  slElm,
                  nextSL
               ) > 0
            ){ /*If I need to swap an element*/
               swap_str_ptrAry(
                  strPtrAry,
                  slElm,
                  nextSL
               );

               lastSL = slElm;
               onSL = slElm;

               while(lastSL >= subSL)
               { /*loop: move swapped element back*/
                  lastSL -= subSL;

                  if(
                     cmp_str_ptrAry(
                        strPtrAry,
                        onSL,
                        lastSL
                     ) > 0
                  ) break; /*moved onSL to its position*/

                  swap_str_ptrAry(
                     strPtrAry,
                     onSL,
                     lastSL
                  );

                  onSL = lastSL;
               } /*Loop: move swapped element back*/
            } /*If I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subSL = (subSL - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*sort_str_ptrAry*/

/*-------------------------------------------------------\
| Fun13: sortSync_str_ptrAry
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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun13 TOC:
   '   - sorts string pointer array in str_ptrAry struct
   '     while keeping an unsigned int array in sync
   '   o fun13 sec01:
   '     - variable declerations
   '   o fun13 sec02:
   '     - find the number of rounds to sort for
   '   o fun13 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of sorting rounds*/
   signed long subSL = 0;
   signed long nextSL = 0;
   signed long lastSL = 0;
   signed long onSL = 0;

   signed long lenSL = strPtrAry->lenSL;

   /*Variables to incurment loops*/
   signed long slIndex = 0;
   signed long slElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(lenSL <= 1)
      return; /*nothing to do*/

   /*Recursion formuia: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subSL = 1; /*Initialzie first array*/

   while(subSL < lenSL - 1)
      subSL = (3 * subSL) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun13 Sec03:
   ^   - sort arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(subSL > 0)
   { /*Loop: all rounds*/
      for(
         slIndex = 0;
         slIndex <= subSL;
         ++slIndex
      ){ /*Loop: though sub array*/
         slElm = slIndex;

         for(
            slElm = slIndex;
            slElm + subSL < lenSL;
            slElm += subSL
         ){ /*Loop: swap elements in subarray*/
            nextSL = slElm + subSL;

            if(
               cmp_str_ptrAry(
                  strPtrAry,
                  slElm,
                  nextSL
               ) > 0
            ){ /*If I need to swap an element*/
               swap_str_ptrAry(
                  strPtrAry,
                  slElm,
                  nextSL
               );

               uiAry[slElm] ^= uiAry[nextSL];
               uiAry[nextSL] ^= uiAry[slElm];
               uiAry[slElm] ^= uiAry[nextSL];

               lastSL = slElm;
               onSL = slElm;

               while(lastSL >= subSL)
               { /*loop: move swapped element back*/
                  lastSL -= subSL;

                  if(
                     cmp_str_ptrAry(
                        strPtrAry,
                        onSL,
                        lastSL
                     ) > 0
                  ) break; /*moved onSL to its position*/

                  swap_str_ptrAry(
                     strPtrAry,
                     onSL,
                     lastSL
                  );

                  uiAry[onSL] ^= uiAry[lastSL];
                  uiAry[lastSL] ^= uiAry[onSL];
                  uiAry[onSL] ^= uiAry[lastSL];

                  onSL = lastSL;
               } /*Loop: move swapped element back*/
            } /*If I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subSL = (subSL - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*sortSync_str_ptrAry*/

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
){
   signed long midSL = 0;
   signed long rightSL = strPtrAry->lenSL - 1;
   signed long leftSL = 0;
   signed long cmpSL = 0;

   while(leftSL <= rightSL)
   { /*Loop: Search for the querys index*/
      midSL = (leftSL + rightSL) >> 1;

      cmpSL =
         eqlNull_ulCp(
            qryStr,
            strPtrAry->strAry[midSL]
         ); /*compare query to array value*/

      if(cmpSL > 0)
          leftSL = midSL + 1;

      else if(cmpSL < 0)
          rightSL = midSL - 1;

      else
         return midSL;
   } /*Loop: Search for the querys index*/

   return -1; /*query not found*/
} /*find_str_ptrAry*/

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
){
   signed long indexSL = 0;

   while(indexSL < strPtrAry->lenSL)
   { /*Loop: Search for the querys index*/
      if(
        ! eqlNull_ulCp(
            qryStr,
            strPtrAry->strAry[indexSL]
         ) /*compare query to array value*/
     ) return indexSL;
 
     ++indexSL;
   } /*Loop: Search for the querys index*/

   return -1; /*query not found*/
} /*findNoSort_str_ptrAry*/

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
