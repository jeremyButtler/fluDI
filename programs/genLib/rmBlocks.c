/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' rmBlocks SOF: Start Of File
'   - remove blocks (start ---, end ---) from file and
'     merge block contents onto a line
'   - also removes all comments
'     - formats: #.*\n; #---.*---#; or #---.*---#
'   o header:
'     - included libraries
'   o fun01: str_rmBlocks
'     - converts blocks in string (---) to single lines
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

/*-------------------------------------------------------\
| Fun01: str_rmBlocks
|   - converts blocks in string (---) to single lines
| Input:
|   - textStr:
|     o c-string with blocks to convert to lines
|   - lenUL:
|     o has length of textStr
|     o 0 if modifying original string
|   - retLenULPtr:
|     o pointer to unsigned long with new string length
|     o can set equal to lenUL
|     o if length not waned, use null
|   - keepCommentBl:
|     o 1: keep all non-block comments
|     o 0: remove all comments
| Output:
|   - Modifies:
|     o retLenULPtr to have length of modified string
|     0 textStr to have not blocks if lenUL is 0
|   - Returns:
|     o modified c-string with textStr with all blocks
|       (start and end with ---) symbols removed and
|       block entries converted to a single line
|       - is new pointer if lenUL is not 0
|       - pointer to textStr if lenUL is 0
|     o 0 for memory errors
\-------------------------------------------------------*/
signed char *
str_rmBlocks(
   signed char *textStr,
   unsigned long lenUL,
   unsigned long *retLenULPtr,
   signed char keepCommentBl
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - converts blocks (#--- ---) to single lines
   '   o fun01 sec01:
   '     - variable declarations
   '   o fun01 sec02:
   '     - allocated return string memory
   '   o fun01 sec03:
   '     - remove blocks
   '   o fun01 sec04:
   '     - return string without blocks
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *tmpStr = textStr;
   signed char *scanStr = textStr;

   signed char *retStr = 0;
   signed char *dupStr = 0;
   signed char blockBl = 0;
   signed char stringBl = 0;

   unsigned long tmpUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - allocated return string memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! retLenULPtr)
      retLenULPtr = &tmpUL;

   *retLenULPtr = 0;

   if(lenUL)
   { /*If: can not modify input string*/
      retStr = malloc((lenUL + 9) * sizeof(signed char));

      if(! retStr)
         return 0;

      dupStr = retStr;
   } /*If: can not modify input string*/

   else
      dupStr = textStr; /*modifing input string*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - remove blocks
   ^   o fun01 sec03 sub01:
   ^     - get past shebang (if present)
   ^   o fun01 sec03 sub02:
   ^     - deal with escaped charaters + start loop
   ^   o fun01 sec03 sub03:
   ^     - deal with strings + start loop
   ^   o fun01 sec03 sub04:
   ^     - deal with comments
   ^   o fun01 sec03 sub05:
   ^     - deal with extra white space
   ^   o fun01 sec03 sub06:
   ^     - detect if working on a block
   ^   o fun01 sec03 sub07:
   ^     - merge lines in the bock + debuging white space
   ^   o fun01 sec03 sub08:
   ^     - move to next characer in block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec03 Sub01:
   *   - get past shebang (if present)
   \*****************************************************/

   if(
         tmpStr[0] == '#' 
      && tmpStr[1] == '!' 
   ){ /*If: have shebang*/
      while(*tmpStr != '\n')
      { /*Loop: copy the shebang*/
          if(*tmpStr == '\0')
             break;

          *dupStr++ = *tmpStr++;
         ++(*retLenULPtr);
      } /*Loop: copy the shebang*/

      if(*tmpStr == '\n')
      { /*If: not at end of string*/
         *dupStr++ = *tmpStr++;
         ++(*retLenULPtr);
      } /*If: not at end of string*/
   } /*If: have shebang*/

   /*****************************************************\
   * Fun01 Sec03 Sub02:
   *   - deal with escaped charaters + start loop
   \*****************************************************/

   while(*tmpStr != '\0')
   { /*Loop: scan for blocks*/
      if(tmpStr[0] == '\\')
      { /*If: have escaped value*/
         ++tmpStr;

         if(*tmpStr == '"')
         { /*If: double quote string*/
            *dupStr++ = '"';
            ++tmpStr;
            ++(*retLenULPtr);

            continue;
         } /*If: double qoute string*/

         else if(stringBl)
         { /*Else: in string, no special symbols*/
            --tmpStr; /*keep \\*/
            goto nextChar_fun01_sec03_sub04;
         } /*Else: in string, no special symbols*/

         else if(*tmpStr == '\\')
         { /*Else If: slash escape*/
            *dupStr++ = '\\';
            ++tmpStr;
            ++(*retLenULPtr);

            continue;
         } /*Else If: slash escape*/

         else if(*tmpStr == '#')
         { /*Else If: comment escape (no comment)*/
            *dupStr++ = '#';
            ++tmpStr;
            ++(*retLenULPtr);

            if(
                  tmpStr[0] == '-'
               && tmpStr[1] == '-'
               && tmpStr[2] == '-'
               && tmpStr[3] != '-'
            ) goto blockEscape_fun01_sec03_sub01;

            continue;
         } /*Else If: comment escape (no comment)*/

         else if(
               tmpStr[0] == '-'
            && tmpStr[1] == '-'
            && tmpStr[2] == '-'
            && tmpStr[3] != '-'
         ){ /*Else If: block escape*/
            blockEscape_fun01_sec03_sub01:;

            *dupStr++ = '-';
            ++tmpStr;
            ++(*retLenULPtr);

            *dupStr++ = '-';
            ++tmpStr;
            ++(*retLenULPtr);

            *dupStr++ = '-';
            ++tmpStr;
            ++(*retLenULPtr);

            continue;
         } /*Else If: block escape*/
      } /*Else If: have escaped value*/

      /**************************************************\
      * Fun01 Sec03 Sub03:
      *   - deal with strings
      \**************************************************/

      if(tmpStr[0] == '"')
      { /*If: have string*/
         stringBl = ! stringBl;
         goto nextChar_fun01_sec03_sub04;
      } /*If: have string*/

      if(stringBl)
         goto nextChar_fun01_sec03_sub04;
         /*middle of string, nothing to check here*/

      /**************************************************\
      * Fun01 Sec03 Sub04:
      *   - deal with comments
      \**************************************************/

      if(tmpStr[0] == '#')
      { /*If: have comment*/

         ++tmpStr;

         if(
               tmpStr[0] == '-'
            && tmpStr[1] == '-'
            && tmpStr[2] == '-'
            && tmpStr[3] != '-'
         ){ /*If: have block comment*/
            tmpStr += 3; /*got off block start*/

            while(*tmpStr != '\0')
            { /*Loop: find end of block*/

               if( *tmpStr == '\\' )
               { /*If: may have an escaped block*/
                  ++tmpStr;

                  while(*tmpStr == '-')
                     ++tmpStr;

                  continue;
               } /*If: may have an escaped block*/

               if(
                     tmpStr[0] == '-'
                  && tmpStr[1] == '-'
                  && tmpStr[2] == '-'
               ){ /*If: have comment block*/
                  tmpStr += 3;

                  if(tmpStr[0] != '-')
                  { /*If: found end block*/
                     if(*tmpStr == '#')
                        ++tmpStr; /*if using ---# end*/

                     break;
                  } /*If: found end block*/

                  while(*tmpStr == '-')
                    ++tmpStr; /*move to end of dashes*/
               } /*If: have comment block*/

               else
                  ++tmpStr;

            } /*Loop: find end of block*/

            goto blockEnd_rmSpace_fun01_sec03_sub05;
         } /*If: have block comment*/

         else
         { /*Else: line comment (single line)*/
            if(
                  keepCommentBl
               && ! blockBl
            ) --tmpStr; /*get back onto comment*/

            while(*tmpStr != '\n')
            { /*Loop: find end of comment*/
               if(*tmpStr == '\0') 
                  break;

               if(
                     keepCommentBl
                  && ! blockBl
               ){ /*If: keeping non-block comments*/
                  *dupStr++ = *tmpStr;
                  ++(*retLenULPtr);

                  if(*retLenULPtr == lenUL)
                     break;
               } /*If: keeping non-block comments*/

               ++tmpStr;
            } /*Loop: find end of comment*/
         } /*Else: line comment (single line)*/

         continue;
      } /*If: have comment*/

      /**************************************************\
      * Fun01 Sec03 Sub05:
      *   - deal with extra white space
      \**************************************************/

      #ifndef BUG
         if(tmpStr[0] == '\n')
         { /*If: new line*/ 
            if(! blockBl)
            { /*If: not on block, keep new line*/
               *dupStr++ = '\n';
               ++(*retLenULPtr);
            } /*If: not on block, keep new line*/

            else
            { /*Else: block, keep one space*/
               *dupStr++ = ' ';
               ++(*retLenULPtr);
            } /*Else: block, keep one space*/

            while(*tmpStr < 33)
            { /*Loop: remove white space*/
               if(*tmpStr == '\0')
                  break;
               ++tmpStr;
            } /*Loop: remove white space*/

            continue;
         } /*If: new line*/ 

         else if(
               tmpStr[0] <= 32
            && tmpStr[1] <= 32
         ){ /*If: extra white space*/

            if(tmpStr[1] == '\n')
               ++tmpStr;
            else
            { /*Else: not end line; remove extra space*/
               while(tmpStr[1] <= 32)
               { /*Loop: remove extra white space*/
                  if(tmpStr[1] == '\0')
                     break;

                  if(tmpStr[1] == '\n')
                     break;

                  ++tmpStr;
               } /*Loop: remove extra white space*/
            } /*Else: not end line; remove extra space*/

            continue;
         } /*If: extra white space*/
      #endif

      /**************************************************\
      * Fun01 Sec03 Sub06:
      *   - detect if working on a block
      \**************************************************/

      if(
            *(tmpStr) != '\\'
         && tmpStr[0] == '-'
         && tmpStr[1] == '-'
         && tmpStr[2] == '-'
         && tmpStr[3] != '-'
      ){ /*If: found a block marker*/
         blockBl = ! blockBl;
         tmpStr += 3; /*get off block*/

         if(! blockBl)
         { /*If: at end of block*/
            blockEnd_rmSpace_fun01_sec03_sub05:;

            scanStr = tmpStr;

            while(*tmpStr < 33)
            { /*Loop: remove extra white space*/
               if(*tmpStr == '\0')
                  break;

               if(*tmpStr == '\n')
                  break;

               ++tmpStr;
            } /*Loop: remove extra white space*/

            if(scanStr == tmpStr)
               continue; /*no white space at end*/

            scanStr = tmpStr;
            --tmpStr;

            while(*scanStr < 33)
            { /*Loop: check if have empty lines*/
               if(*scanStr == '\0')
               { /*If: end of string*/
                  tmpStr = scanStr;
                  break;
               } /*If: end of string*/

               if(*scanStr == '\n')
                  tmpStr = scanStr;

               ++scanStr;
            } /*Loop: check if have empty lines*/
         } /*If: at end of block*/

         continue;
      } /*If: found a block marker*/

      else if(! blockBl)
         goto nextChar_fun01_sec03_sub04;

      /**************************************************\
      * Fun01 Sec03 Sub07:
      *   - merge lines in the bock + debuging white space
      \**************************************************/

      if(
            tmpStr[0] == '\n'
         || tmpStr[0] == '\r'
      ){ /*If: is new line in a block*/
         ++tmpStr;
         continue;
      } /*If: is new line in a block*/

      #ifdef BUG
         /*preserves indent white space in non-blocks*/

         if(
               tmpStr[0] <= 32
            && tmpStr[1] <= 32
            && tmpStr[1] != '\n'
         ){ /*If: extra white space*/
            while(tmpStr[1] <= 32)
            { /*Loop: remove extra white space*/
               if(tmpStr[0] == '\0')
                  break;

               if(tmpStr[0] == '\n')
                  break;

               ++tmpStr;
            } /*Loop: remove extra white space*/

            continue;
         } /*If: extra white space*/
      #endif

      /**************************************************\
      * Fun01 Sec03 Sub08:
      *   - move to next characer in block
      \**************************************************/

      nextChar_fun01_sec03_sub04:;
         *dupStr++ = *tmpStr++;
         ++(*retLenULPtr);

         if(*retLenULPtr == lenUL)
            break;
   } /*Loop: scan for blocks*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^   - return string without blocks
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *dupStr = '\0';
   return retStr;
} /*str_rmBlocks*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconveint / not possible, this code is under the
:   MIT license
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
