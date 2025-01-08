/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' numToStr SOF: Start Of File
'   - function(s) to convert a number to a string
'   o header:
'     - included libraries
'   o fun01: numToStr
'     - converts a number to a c-string
'   o fun02: backNumToStr
'     - converts a number to a c-string backwards
'   o fun03: double_numToStr
'     - converts a double to a c-string
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
#endif

#include "numToStr.h"

/*-------------------------------------------------------\
| Fun01: numToStr
|   - converts a number to a c-string
| Input:
|   - cstr:
|     o c-string to hold the converted number
|   - numUL:
|     o number to convert
| Output:
|   - Modifies:
|     o cstr to have the number (adds a '\0' at the end)
|   - Returns:
|     o number of characters copied to cstr
\-------------------------------------------------------*/
unsigned int
numToStr(
   signed char *cstr,
   unsigned long numUL
){
   unsigned int retUI = 0;
   signed char *startStr = 0;
   signed char *endStr = 0;

   do{
      cstr[retUI] = (numUL % 10) + 48;
      ++retUI;
      numUL /= 10;
   } while(numUL);

   cstr[retUI] = '\0';

   endStr = &cstr[retUI - 1];
   startStr = cstr;

   while(startStr < endStr)
   { /*Loop: Reverse array*/
      *endStr ^= *startStr;
      *startStr ^= *endStr;
      *endStr-- ^= *startStr++;
   } /*Loop: Reverse array*/

   return retUI;
} /*numToStr*/

/*-------------------------------------------------------\
| Fun02: backNumToStr
|   - converts a number to a c-string backwards
| Input:
|   - cstr:
|     o c-string to hold the converted number
|   - numUL:
|     o number to convert
| Output:
|   - Modifies:
|     o cstr to have the backwards number (adds a '\0' at
|       the end)
|   - Returns:
|     o number of characters copied to cstr
\-------------------------------------------------------*/
unsigned int
backwards_numToStr(
   signed char *cstr,
   unsigned long numUL
){
   unsigned int retUI = 0;

   do{
      cstr[retUI] = (numUL % 10) + 48;
      ++retUI;
      numUL /= 10;
   } while(numUL);

   cstr[retUI] = '\0';
   return retUI;
} /*backNumToStr*/

/*-------------------------------------------------------\
| Fun03: double_numToStr
|   - converts a double to a c-string
| Input:
|   - cstr:
|     o c-string to hold the converted number
|   - numDbl:
|     o number to convert
|   - decUC:
|     o number decimal digits to keep
|       (max_dblDecimal_numToStr)
| Output:
|   - Modifies:
|     o cstr to have the number (adds a '\0' at the end)
|   - Returns:
|     o number of characters copied to cstr
\-------------------------------------------------------*/
unsigned int
double_numToStr(
   signed char *cstr,
   double numDbl,
   unsigned char decUC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - converts a double to a c-string
   '   o fun03 sec01:
   '     - variable declarations
   '   o fun03 sec02:
   '     - convert non-decimal part to c-string
   '   o fun03 sec03:
   '     - convert decimal part to number
   '   o fun03 sec04:
   '     - check if need to round
   '   o fun03 sec05:
   '     - add null and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned int retUI = 0;
   signed char *startStr = 0;
   signed char *endStr = 0;
   signed long numSL = 0;   /*has non-decimal number*/
   double decDbl = 0;       /*has decimal fraction*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - convert non-decimal part to c-string
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numDbl < 0)
   { /*If: negative number*/
      numDbl *= -1;
      cstr[retUI++] = '-'; /*negative number*/
   } /*If: negative number*/

   numSL = (signed long) numDbl;
   decDbl = numDbl - numSL;

   do{
      cstr[retUI++] = (numSL % 10) + 48;
      numSL /= 10;
   } while(numSL);

   endStr = &cstr[retUI - 1];
   startStr = cstr;

   while(startStr < endStr)
   { /*Loop: Reverse array (is bacwards)*/
      *endStr ^= *startStr;
      *startStr ^= *endStr;
      *endStr-- ^= *startStr++;
   } /*Loop: Reverse array (is backwards)*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - convert decimal part to number
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(decUC > max_dblDecimal_numToStr)
      decUC = max_dblDecimal_numToStr;
      /*make sure not out of percision range*/

   if(decDbl != 0)
   { /*If: have decimal part*/
      cstr[retUI++] = '.';

      while(decDbl != 0)
      { /*Loop: add in decimals*/
         if(! decUC)
            break;

         decDbl *= 10; /*convert decimal to whole number*/
         cstr[retUI++] = ((signed char) decDbl) + 48;

         decDbl -= ((signed char) decDbl);
            /*remove decimal*/

         --decUC;
      } /*Loop: add in decimals*/
   } /*If: have decimal part*/

   cstr[retUI] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - check if need to round
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(decDbl * 10 >= 5)
   { /*If: need to round up*/
      startStr = &cstr[retUI - 1];

      if(*startStr != '9')
         ++(*startStr);

      else
      { /*Else: need to round multiple values*/
         while(*startStr == '9')
         { /*Loop: round up*/
            if(
                  startStr == cstr
               || *startStr == '-'
            ){ /*If: need to round up on frist digit*/
               endStr = startStr + 1;

               while(*endStr != '\0')        
               { /*Loop: copy values*/
                  *endStr ^= *startStr;
                  *startStr ^= *endStr;
                  *endStr++ ^= *startStr;
               } /*Loop: copy values*/

               *endStr = *startStr; /*copy last value*/
               *startStr = '1';
               ++retUI; /*account for one more digit*/

               goto rounded_fun03_sec03;
            } /*If: need to round up on frist digit*/

            if(*startStr != '.')
               *startStr = '0';

            --startStr;
         } /*Loop: round up*/

         ++(*startStr); /*round up last digit*/

         rounded_fun03_sec03:;
      } /*Else: need to round multiple values*/
   } /*If: need to round up*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec05:
   ^   - add null and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   return retUI;
} /*double_numToStr*/

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
