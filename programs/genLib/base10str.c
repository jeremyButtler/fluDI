/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' base10str SOF: Start Of Header
'   - functions to converts c-strings to numbers
'   o header:
'     - included libraries
'   o fun01: strToUL_base10str
'     - converts a base 10 c-string to an unsigned long
'   o fun02: strToUI_base10str
'     - converts a base 10 c-string to an unsigned int
'   o fun03: strToUS_base10str
'     - converts a base 10 c-string to an unsigned short
'   o fun04: strToUC_base10str
'     - converts a base 10 c-string to an unsigned char
'   o fun05: strToSL_base10str
'     - converts base 10 c-string to a signed long
'   o fun06: strToSI_base10str
'     - converts base 10 c-string to a signed int
'   o fun07: strToSS_base10str
'     - converts base 10 c-string to a signed short
'   o fun08: strToSC_base10str
'     - converts base 10 c-string to a signed char
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "base10str.h"

/*-------------------------------------------------------\
| Fun01: strToUL_base10str
|   - converts a base 10 c-string to an unsigned long
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retUL:
|     o pointer to unisnged long to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retUL to hold the unsigned long
| Note:
|   - only converts till retUL goes out off bounds
\-------------------------------------------------------*/
unsigned int
strToUL_base10str(
   signed char *inStr,
   unsigned long *retUL
){
   unsigned int uiDig = 0;
   unsigned long tmpUL = 0;

   *retUL = 0;

   for(
      uiDig = 0;
      uiDig < def_ulDig_base10str;
      ++uiDig
   ){ /*Loop: Convert digits with no overflow concern*/
      if(
            inStr[uiDig] > 47
         && inStr[uiDig] < 58
      ) *retUL = *retUL * 10 + inStr[uiDig] - 48;

      else
         break;
   } /*Loop: Convert digits with no overflow concern*/

   if(
         inStr[uiDig] > 47
      && inStr[uiDig] < 58
   ){ /*If have a last digit (could overflow)*/
      if(*retUL <= def_maxUL_base10str / 10)
      { /*If: number is < the max unsigned int size/10*/
         tmpUL = 10 * *retUL + inStr[uiDig] - 48;

         if(tmpUL > 9)
         { /*If: I have room for 10th digit*/
            *retUL = tmpUL;
            ++uiDig;
         } /*If: I have room for 10th digit*/
      } /*If: number is < the max unsigned int size/10*/
   } /*If have a last digit (could overflow)*/

   return uiDig;
} /*strToUL_base10str*/

/*-------------------------------------------------------\
| Fun02: strToUI_base10str
|   - converts a base 10 c-string to an unsigned int
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retUI:
|     o pointer to unisnged int to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retUI to hold the unsigned int
| Note:
|   - only converts till retUI goes out off bounds
\-------------------------------------------------------*/
unsigned int
strToUI_base10str(
   signed char *inStr,
   unsigned int *retUI
){
   unsigned int uiDig = 0;
   unsigned int tmpUI = 0;

   *retUI = 0;

   for(
      uiDig = 0;
      uiDig < def_uiDig_base10str;
      ++uiDig
   ){ /*Loop: Convert digits with no overflow concern*/
      if(
            inStr[uiDig] > 47
         && inStr[uiDig] < 58
      ) *retUI = *retUI * 10 + inStr[uiDig] - 48;

      else
         break;
   } /*Loop: Convert digits with no overflow concern*/

   if(
         inStr[uiDig] > 47
      && inStr[uiDig] < 58
   ){ /*If have a last digit (couid overflow)*/
      if(*retUI <= def_maxUI_base10str / 10)
      { /*If: number is < the max unsigned int size/10*/
         tmpUI = 10 * *retUI + inStr[uiDig] - 48;

         if(tmpUI > 9)
         { /*If: I have room for 10th digit*/
            *retUI = tmpUI;
            ++uiDig;
         } /*If: I have room for 10th digit*/
      } /*If: number is < the max unsigned int size/10*/
   } /*If have a last digit (couid overflow)*/

   return uiDig;
} /*strToUI_base10str*/

/*-------------------------------------------------------\
| Fun03: strToUS_base10str
|   - converts a base 10 c-string to an unsigned short
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retUS:
|     o pointer to unisnged short to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retUS to hold the unsigned short
| Note:
|   - only converts till retUS goes out off bounds
\-------------------------------------------------------*/
unsigned int
strToUS_base10str(
   signed char *inStr,
   unsigned short *retUS
){
   unsigned int uiDig = 0;
   unsigned short tmpUS = 0;

   *retUS = 0;

   for(
      uiDig = 0;
      uiDig < def_usDig_base10str;
      ++uiDig
   ){ /*Loop: Convert digits with no overflow concern*/
      if(
            inStr[uiDig] > 47
         && inStr[uiDig] < 58
      ) *retUS = *retUS * 10 + inStr[uiDig] - 48;

      else
         break;
   } /*Loop: Convert digits with no overflow concern*/

   if(
         inStr[uiDig] > 47
      && inStr[uiDig] < 58
   ){ /*If have a last digit (cousd overflow)*/
      if(*retUS <= def_maxUS_base10str / 10)
      { /*If: number is < the max unsigned int size/10*/
         tmpUS = 10 * *retUS + inStr[uiDig] - 48;

         if(tmpUS > 9)
         { /*If: I have room for 10th digit*/
            *retUS = tmpUS;
            ++uiDig;
         } /*If: I have room for 10th digit*/
      } /*If: number is < the max unsigned int size/10*/
   } /*If have a last digit (cousd overflow)*/

   return uiDig;
} /*strToUS_base10str*/

/*-------------------------------------------------------\
| Fun04: strToUC_base10str
|   - converts a base 10 c-string to an unsigned char
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retUC:
|     o pointer to unisnged char to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retUC to hold the unsigned char
| Note:
|   - only converts till retUC goes out off bounds
\-------------------------------------------------------*/
unsigned int
strToUC_base10str(
   signed char *inStr,
   unsigned char *retUC
){
   unsigned int uiDig = 0;
   unsigned char tmpUC = 0;

   *retUC = 0;

   for(
      uiDig = 0;
      uiDig < def_ucDig_base10str;
      ++uiDig
   ){ /*Loop: Convert digits with no overflow concern*/
      if(
            inStr[uiDig] > 47
         && inStr[uiDig] < 58
      ) *retUC = *retUC * 10 + inStr[uiDig] - 48;

      else
         break;
   } /*Loop: Convert digits with no overflow concern*/

   if(
         inStr[uiDig] > 47
      && inStr[uiDig] < 58
   ){ /*If have a last digit (coucd overflow)*/
      if(*retUC <= def_maxUC_base10str / 10)
      { /*If: number is < the max unsigned int size/10*/
         tmpUC = 10 * *retUC + inStr[uiDig] - 48;

         if(tmpUC > 9)
         { /*If: I have room for 10th digit*/
            *retUC = tmpUC;
            ++uiDig;
         } /*If: I have room for 10th digit*/
      } /*If: number is < the max unsigned int size/10*/
   } /*If have a last digit (coucd overflow)*/

   return uiDig;
} /*strToUC_base10str*/

/*-------------------------------------------------------\
| Fun05: strToSL_base10str
|   - converts base 10 c-string to a signed long
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retSL:
|     o pointer to signed long to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retSL to hold the signed long
| Note:
|   - only converts until retSL overflows
\-------------------------------------------------------*/
unsigned int
strToSL_base10str(
   signed char *inStr,
   signed long *retSL
){
   unsigned int uiDig = 0;
   signed long tmpSL = 0;

   *retSL = 0;

   for(
      uiDig = (inStr[0] == '-');
      uiDig < def_slDig_base10str;
      ++uiDig
   ){ /*Loop: Convert digits with no overflow concern*/
      if(
            inStr[uiDig] > 47
         && inStr[uiDig] < 58
      ) *retSL = *retSL * 10 + inStr[uiDig] - 48;
      else
         break;
   } /*Loop: Convert digits with no overflow concern*/

   if(
         inStr[uiDig] > 47
      && inStr[uiDig] < 58
   ){ /*If have a last digit (could overflow)*/
      if(*retSL <= def_maxSL_base10str / 10)
      { /*If: number is < the max unsigned int size/10*/
         tmpSL = 10 * *retSL + inStr[uiDig] - 48;

         if(tmpSL > 9)
         { /*If: I have room for 10th digit*/
            *retSL = tmpSL;
            ++uiDig;
         } /*If: I have room for 10th digit*/
      } /*If: number is < the max unsigned int size/10*/
   } /*If have a last digit (could overflow)*/

   if(inStr[0] == '-')
      *retSL = -(*retSL);

   return uiDig;
} /*strToSL_base10str*/

/*-------------------------------------------------------\
| Fun06: strToSI_base10str
|   - converts base 10 c-string to a signed int
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retSI:
|     o pointer to signed int to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retSI to hold the signed int
| Note:
|   - only converts until retSI overflows
\-------------------------------------------------------*/
unsigned int
strToSI_base10str(
   signed char *inStr,
   signed int *retSI
){
   unsigned int uiDig = 0;
   signed int tmpSI = 0;

   *retSI = 0;

   for(
      uiDig = (inStr[0] == '-');
      uiDig < def_siDig_base10str;
      ++uiDig
   ){ /*Loop: Convert digits with no overflow concern*/
      if(
            inStr[uiDig] > 47
         && inStr[uiDig] < 58
      ) *retSI = *retSI * 10 + inStr[uiDig] - 48;
      else
         break;
   } /*Loop: Convert digits with no overflow concern*/

   if(
         inStr[uiDig] > 47
      && inStr[uiDig] < 58
   ){ /*If have a last digit (could overflow)*/
      if(*retSI <= def_maxSI_base10str / 10)
      { /*If: number is < the max unsigned int size/10*/
         tmpSI = 10 * *retSI + inStr[uiDig] - 48;

         if(tmpSI > 9)
         { /*If: I have room for 10th digit*/
            *retSI = tmpSI;
            ++uiDig;
         } /*If: I have room for 10th digit*/
      } /*If: number is < the max unsigned int size/10*/
   } /*If have a last digit (could overflow)*/

   if(inStr[0] == '-')
      *retSI = -(*retSI);

   return uiDig;
} /*strToSI_base10str*/

/*-------------------------------------------------------\
| Fun07: strToSS_base10str
|   - converts base 10 c-string to a signed short
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retSS:
|     o pointer to signed short to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retSS to hold the signed short
| Note:
|   - only converts until retSS overflows
\-------------------------------------------------------*/
unsigned short
strToSS_base10str(
   signed char *inStr,
   signed short *retSS
){
   unsigned short uiDig = 0;
   signed short tmpSS = 0;

   *retSS = 0;

   for(
      uiDig = (inStr[0] == '-');
      uiDig < def_ssDig_base10str;
      ++uiDig
   ){ /*Loop: Convert digits with no overflow concern*/
      if(
            inStr[uiDig] > 47
         && inStr[uiDig] < 58
      ) *retSS = *retSS * 10 + inStr[uiDig] - 48;
      else
         break;
   } /*Loop: Convert digits with no overflow concern*/

   if(
         inStr[uiDig] > 47
      && inStr[uiDig] < 58
   ){ /*If have a last digit (could overflow)*/
      if(*retSS <= def_maxSS_base10str / 10)
      { /*If: number is < the max unsigned short ssze/10*/
         tmpSS = 10 * *retSS + inStr[uiDig] - 48;

         if(tmpSS > 9)
         { /*If: I have room for 10th digit*/
            *retSS = tmpSS;
            ++uiDig;
         } /*If: I have room for 10th digit*/
      } /*If: number is < the max unsigned short ssze/10*/
   } /*If have a last digit (could overflow)*/

   if(inStr[0] == '-')
      *retSS = -(*retSS);

   return uiDig;
} /*strToSS_base10str*/

/*-------------------------------------------------------\
| Fun08: strToSC_base10str
|   - converts base 10 c-string to a signed char
| Input:
|   - inStr:
|     o c-string with base 10 number to convert
|   - retSC:
|     o pointer to signed char to hold converted number
| Output:
|   - Returns:
|     o number of characters converted
|   - Modifies:
|     o retSC to hold the signed char
| Note:
|   - only converts until retSC overflows
\-------------------------------------------------------*/
unsigned char
strToSC_base10str(
   signed char *inStr,
   signed char *retSC
){
   unsigned char uiDig = 0;
   signed char tmpSC = 0;

   *retSC = 0;

   for(
      uiDig = (inStr[0] == '-');
      uiDig < def_scDig_base10str;
      ++uiDig
   ){ /*Loop: Convert digits with no overflow concern*/
      if(
            inStr[uiDig] > 47
         && inStr[uiDig] < 58
      ) *retSC = *retSC * 10 + inStr[uiDig] - 48;
      else
         break;
   } /*Loop: Convert digits with no overflow concern*/

   if(
         inStr[uiDig] > 47
      && inStr[uiDig] < 58
   ){ /*If have a last digit (could overflow)*/
      if(*retSC <= def_maxSC_base10str / 10)
      { /*If: number is < the max unsigned char scze/10*/
         tmpSC = 10 * *retSC + inStr[uiDig] - 48;

         if(tmpSC > 9)
         { /*If: I have room for 10th digit*/
            *retSC = tmpSC;
            ++uiDig;
         } /*If: I have room for 10th digit*/
      } /*If: number is < the max unsigned char scze/10*/
   } /*If have a last digit (could overflow)*/

   if(inStr[0] == '-')
      *retSC = -(*retSC);

   return uiDig;
} /*strToSC_base10str*/

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
