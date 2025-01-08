/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' ulCp SOF: Start Of File
'   - uses longs to copy contents of a string
'   - these functions are slower than strcpy and strlen,
'     but they allow deliminators to be used.
'   - DOES NOT WORK ON **BIG ENDIN** SYSTEM
'   o header:
'     - included libraries
'   deliminator detection/making:
'     o .h fun01: rshiftByte_ulCp
'       - shift byte right by x bytes; over shifts go to 0
'     o .h fun02: mkDelim_ulCp
'       - make ulong_ulCp delimintor from character
'     o .h fun03: ifDelim_ulCp
'       - check if have deliminator in long
'     o .h fun04: ifNullDelim_ulCp
'       - checks for null or deliminator in ulong_ulCp
'     o .h fun05: ifNull_ulCp
'       - checks to see if null in ulong_ulCp
'     o .h fun06: ifLineUnix_ulCp
'       - checks if unix line end ('\n' or '\0')
'     o .h fun07: ifEndLine_ulCp
'       - checks if input long is end of line
'     o .h fun08: ifWhite_ulCp
'       - check if input long has white space
'         (' ', '\r', '\n', '\t', '\0')
'   String copying:
'     o fun09: cpLen_ulCp
'       - copies cpStr into dupStr using ulong_ulCps
'     o fun10: cpDelim_ulCp
'       - copies string until deliminator is found
'     o fun11: cpStr_ulCp
'       - copies string until \0
'     o fun12: cpLineUnix_ulCp
'       - copies string until end of line (\0, \n)
'     o fun13: cpLine_ulCp
'       - copies string until end of line (\0, \r, \n)
'       - ingores unused ascii characters (> 32, not '\t')
'     o fun14: cpWhite_ulCp
'       - copies string until white space
'   String length:
'     o fun15: lenStr_ulCp
'       - finds the length of a string using ulong_ulCps
'     o fun16: lenStrNull_ulCp
'       - finds the length of a string using ulong_ulCps
'         but also stops at null '\0'
'     o fun17: endStr_ulCp
'       - finds the end of a c-string ('\0')
'     o fun18: endLineUnix_ulCp
'      - finds the end of a c-string. This assumes that
'        the line ends in '\0' or '\n'
'     o fun19: endLine_ulCp
'       - finds the end of a c-string (all OS's; looks for
'         '\0', '\n', and '\r')
'     o fun20: endWhite_ulCp
'       - finds 1st white space ('\t',' ','\n','\0','\r')
'   String comparision:
'     o fun21: eql_ulCp
'       - compares two strings until deliminator is found
'     o fun22: eqlNull_ulCp
'       - compares two strings until null is found
'     o fun23: eqlWhite_ulCp
'       - compares two strings until white space is found
'   String cleanup:
'     o fun24: rmWhite_ulCp
'       - removes white space from c-string
'   String swap:
'      o fun25: swapDelim_ulCp
'        - swaps two strings until deliminator is found
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

#include "ulCp.h"

/*-------------------------------------------------------\
| Fun09: cpLen_ulCp
|   - copies cpStr into dupStr using ulong_ulCps
| Input:
|   - dupStr:
|     o pointer to string to copy cpStr into
|   - cpStr:
|     o pointer to string to copy
|   - lenUI:
|     o length of cpStr
| Output:
|   - Modifies:
|     o  dupStr to hold lenUI characters from cpStr
\-------------------------------------------------------*/
void
cpLen_ulCp(
   signed char *dupStr,
   signed char *cpStr,
   unsigned int lenUI
){
   ulong_ulCp *cpUL = (ulong_ulCp *) cpStr;
   ulong_ulCp *dupUL = (ulong_ulCp *) dupStr;
   unsigned int uiChar = 0;

   for(
      uiChar= 0;
      uiChar < (lenUI >> def_shiftULBy_ulCp);
      ++uiChar
   ) dupUL[uiChar] = cpUL[uiChar];

   for(
      uiChar = ( lenUI - (lenUI & def_modUL_ulCp) );
      uiChar < lenUI;
     ++uiChar
   ) dupStr[uiChar] = cpStr[uiChar];

   dupStr[lenUI] = '\0';
} /*cpLen_ulCp*/

/*-------------------------------------------------------\
| Fun10: cpDelim_ulCp
|   - copies string until deliminator is found
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
|   - delimUL:
|     o delminator to end at (as long). Use makeULDelim
|       to build this deliminator
|   - delimSC:
|     o delminator (as char) to stop copying at
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
unsigned int
cpDelim_ulCp(
   signed char *dupStr,
   signed char *cpStr,
   ulong_ulCp delimUL,
   signed char delimSC
){
   ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
   ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
   signed char *dupTmpStr = 0;

   while( ! ifDelim_ulCp(*cpUL, delimUL) )
      *dupUL++ = *cpUL++;

   cpStr = (signed char *) cpUL;
   dupTmpStr = (signed char *) dupUL;

   while( *cpStr != delimSC )
      *dupTmpStr++ = *cpStr++;

   *dupTmpStr = '\0';
   return dupTmpStr - dupStr; /*number of chars copied*/
} /*cpDelim_ulCp*/

/*-------------------------------------------------------\
| Fun11: cpStr_ulCp
|   - copies string until \0
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
unsigned int
cpStr_ulCp(
   signed char *dupStr,
   signed char *cpStr
){
   ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
   ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
   signed char *dupTmpStr = 0;

   while( ! ifNull_ulCp(*cpUL) )
      *dupUL++ = *cpUL++;

   cpStr = (signed char *) cpUL;
   dupTmpStr = (signed char *) dupUL;

   while( *cpStr != '\0' )
      *dupTmpStr++ = *cpStr++;

   *dupTmpStr = '\0';
   return dupTmpStr - dupStr; /*number of char copied*/
} /*cpStr_ulCp*/

/*-------------------------------------------------------\
| Fun12: cpLineUnix_ulCp
|   - copies string until end of line (\0, \n)
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
unsigned int
cpLineUnix_ulCp(
   signed char *dupStr,
   signed char *cpStr
){
   ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
   ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
   signed char *dupTmpStr = 0;

   while( ! ifLineUnix_ulCp(*cpUL) )
      *dupUL++ = *cpUL++;

   cpStr = (signed char *) cpUL;
   dupTmpStr = (signed char *) dupUL;

   while( *cpStr & (~ '\n') )
      *dupTmpStr++ = *cpStr++;

   *dupTmpStr = '\0';
   return dupTmpStr - dupStr; /*number of char copied*/
} /*cpLineUnix_ulCp*/

/*-------------------------------------------------------\
| Fun13: cpLine_ulCp
|   - copies string until end of line (\0, \r, \n)
|   - ingores all unused ascii characters (> 32, not '\t')
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
unsigned int
cpLine_ulCp(
   signed char *dupStr,
   signed char *cpStr
){
   ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
   ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
   signed char *dupTmpStr = 0;

   while( ! ifEndLine_ulCp(*cpUL) )
      *dupUL++ = *cpUL++;

   cpStr = (signed char *) cpUL;
   dupTmpStr = (signed char *) dupUL;

   while(
         *cpStr > '\r'  /*\r > \t > \n > \0*/
      || *cpStr == '\t' /*so catch tab case*/
   ) *dupTmpStr++ = *cpStr++;

   *dupTmpStr = '\0';
   return dupTmpStr - dupStr; /*number of char copied*/
} /*cpLine_ulCp*/

/*-------------------------------------------------------\
| Fun14: cpWhite_ulCp
|   - copies string until white space
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
unsigned int
cpWhite_ulCp(
   signed char *dupStr,
   signed char *cpStr
){
   ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
   ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
   signed char *dupTmpStr = 0;

   while( ! ifWhite_ulCp(*cpUL) )
      *dupUL++ = *cpUL++;

   cpStr = (signed char *) cpUL;
   dupTmpStr = (signed char *) dupUL;

   while( *cpStr > 33 )
      *dupTmpStr++ = *cpStr++;

   *dupTmpStr = '\0';
   return dupTmpStr - dupStr; /*number of char copied*/
} /*cpWhite_ulCp*/


/*-------------------------------------------------------\
| Fun15: lenStr_ulCp
|   - finds the length of a string using ulong_ulCps
| Input:
|   - inStr:
|     o c-string or string with deliminator to find length
|       of
|   - delimUL:
|     o deliminator (as ulong_ulCp (fun02)) at end of
|       string
|   - delimSC:
|     o deliminator (as char) at end of string
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
lenStr_ulCp(
   signed char *inStr,
   ulong_ulCp delimUL,
   signed char delimSC
){
   ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
   signed char *tmpStr = 0;

   while( ! ifDelim_ulCp(*ptrUL, delimUL) )
      ++ptrUL;

   tmpStr = (signed char *) ptrUL;

   while(*tmpStr != delimSC)
      ++tmpStr;

   return tmpStr - inStr;
} /*lenStr_ulCp*/

/*-------------------------------------------------------\
| Fun16: lenStrNull_ulCp
|   - finds the length of a string using ulong_ulCps
|     but also stops at null '\0'
| Input:
|   - inStr:
|     o c-string or string with deliminator to find length
|       of
|   - delimUL:
|     o deliminator (as ulong_ulCp (fun02)) at end of
|       string
|   - delimSC:
|     o deliminator (as char) at end of string
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
lenStrNull_ulCp(
   signed char *inStr,
   ulong_ulCp delimUL,
   signed char delimSC
){
   ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
   signed char *tmpStr = 0;

   while( ! ifNullDelim_ulCp(*ptrUL, delimUL) )
      ++ptrUL;

   tmpStr = (signed char *) ptrUL;

   while( (*tmpStr ^ delimSC) )
      ++tmpStr;

   return tmpStr - inStr;
} /*lenStrNull_ulCp*/

/*-------------------------------------------------------\
| Fun17: endStr_ulCp
|   - finds the end of a c-string ('\0')
| Input:
|   - inStr:
|     o c-string or string to look for end in
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
endStr_ulCp(
   signed char *inStr
){
   ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
   signed char *tmpStr = 0;

   while( ! ifNull_ulCp(*ptrUL) )
      ++ptrUL;

   tmpStr = (signed char *) ptrUL;

   while( *tmpStr != '\0' )
      ++tmpStr;

   return tmpStr - inStr;
} /*endStr_ulCp*/

/*-------------------------------------------------------\
| Fun18: endLineUnix_ulCp
|   - finds the end of a c-string. This assumes that the
|     line ends in an '\0' or an '\n'
| Input:
|   - inStr:
|     o c-string or string to look for end in
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
endLineUnix_ulCp(
   signed char *inStr
){
   ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
   signed char *tmpStr = 0;

   while( ! ifLineUnix_ulCp(*ptrUL) )
      ++ptrUL;

   tmpStr = (signed char *) ptrUL;

   while( *tmpStr & (~ '\n') )
      ++tmpStr;

   return tmpStr - inStr;
} /*endLineUnix_ulCp*/

/*-------------------------------------------------------\
| Fun19: endLine_ulCp
|   - finds the end of a c-string (all OS's; looks for
|     '\0', '\n', and '\r')
|   - ingores all unused ascii characters (> 32, not '\t')
| Input:
|   - inStr:
|     o c-string or string to look for end in
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
endLine_ulCp(
   signed char *inStr
){
   ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
   signed char *tmpStr = 0;

   while( ! ifEndLine_ulCp(*ptrUL) )
      ++ptrUL;

   tmpStr = (signed char *) ptrUL;

   while(
         *tmpStr > '\r'  /*\r > \t > \n > \0*/
      || *tmpStr == '\t' /*so catch tab case*/
   ) ++tmpStr;

   return tmpStr - inStr;
} /*endLine_ulCp*/

/*-------------------------------------------------------\
| Fun20: endWhite_ulCp
|   - finds 1st white space ('\t', ' ', '\n', '\0', '\r')
| Input:
|   - inStr:
|     o c-string or string to look for whitespace in
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
unsigned int
endWhite_ulCp(
   signed char *inStr
){
   ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
   signed char *tmpStr = 0;

   while( ! ifWhite_ulCp(*ptrUL) )
      ++ptrUL;

   tmpStr = (signed char *) ptrUL;

   while( *tmpStr > 33)
      ++tmpStr;

   return tmpStr - inStr;
} /*endWhite_ulCp*/

/*-------------------------------------------------------\
| Fun21: eql_ulCp
|   - compares two strings until deliminator is found
| Input:
|   - qryStr:
|     o Pointer to query string
|   - refStr:
|     o Pointer to reference strin
|   - delimUL:
|     o delminator to end at (as long). Use makeULDelim
|       to build this deliminator
|   - delimSC:
|     o delminator (as char) to stop copying at
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query > reference
|     o < 0 if query < reference
| Note:
|   - This will likely not be very good at comparing
|     short strings.
\-------------------------------------------------------*/
signed long
eql_ulCp(
   signed char *qryStr,
   signed char *refStr,
   ulong_ulCp delimUL,
   signed char delimSC
){
   ulong_ulCp *qryUL = (ulong_ulCp *) qryStr;
   ulong_ulCp *refUL = (ulong_ulCp *) refStr;

   while( ! ifDelim_ulCp(*qryUL, delimUL) )
   { /*Loop: Copy cpStr to dupStr*/
      if(*qryUL != *refUL)
         break;

      ++qryUL;
      ++refUL;
   } /*Loop: Copy cpStr to dupStr*/

   qryStr = (signed char *) qryUL;
   refStr = (signed char *) refUL;

   while( *qryStr != delimSC )
   { /*Loop: find difference*/
      if(*qryStr != *refStr)
         break;

      ++qryStr;
      ++refStr;
   } /*Loop: find difference*/

   return *qryStr - *refStr;
} /*eql_ulCp*/

/*-------------------------------------------------------\
| Fun22: eqlNull_ulCp
|   - compares two strings until null is found
| Input:
|   - qryStr:
|     o Pointer to query string
|   - refStr:
|     o Pointer to reference strin
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query > reference
|     o < 0 if query < reference
| Note:
|   - this will likely not be very good at comparing
|     short strings.
\-------------------------------------------------------*/
signed long
eqlNull_ulCp(
   signed char *qryStr,
   signed char *refStr
){
   ulong_ulCp *qryUL = (ulong_ulCp *) qryStr;
   ulong_ulCp *refUL = (ulong_ulCp *) refStr;

   while( ! ifNull_ulCp(*qryUL) )
   { /*Loop: Copy cpStr to dupStr*/
      if(*qryUL != *refUL)
         break;

      ++qryUL;
      ++refUL;
   } /*Loop: Copy cpStr to dupStr*/

   qryStr = (signed char *) qryUL;
   refStr = (signed char *) refUL;

   while( *qryStr != (signed char) '\0' )
   { /*Loop: find difference*/
      if(*qryStr != *refStr)
         break;

      ++qryStr;
      ++refStr;
   } /*Loop: find difference*/

   return *qryStr - *refStr;
} /*eqlNull_ulCp*/

/*-------------------------------------------------------\
| Fun23: eqlWhite_ulCp
|   - compares two strings until white space is found
| Input:
|   - qryStr:
|     o Pointer to query string
|   - refStr:
|     o Pointer to reference strin
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query > reference
|     o < 0 if query < reference
| Note:
|   - this will likely not be very good at comparing
|     short strings.
\-------------------------------------------------------*/
signed long
eqlWhite_ulCp(
   signed char *qryStr,
   signed char *refStr
){
   ulong_ulCp *qryUL = (ulong_ulCp *) qryStr;
   ulong_ulCp *refUL = (ulong_ulCp *) refStr;

   while( ! ifWhite_ulCp(*qryUL) )
   { /*Loop: Copy cpStr to dupStr*/
      if(*qryUL != *refUL)
         break;

      ++qryUL;
      ++refUL;
   } /*Loop: Copy cpStr to dupStr*/

   qryStr = (signed char *) qryUL;
   refStr = (signed char *) refUL;

   while( *qryStr > 32 )
   { /*Loop: find difference*/
      if(*qryStr != *refStr)
         break;

      ++qryStr;
      ++refStr;
   } /*Loop: find difference*/

   if(*refStr > 32)
      return *qryStr - *refStr;

   return 0;
} /*eqlWhite_ulCp*/

/*-------------------------------------------------------\
| Fun24: rmWhite_ulCp
|   - removes white space from c-string
| Input:
|   - inStr:
|     o c-string to remove whites space
| Output:
|   - Modifies:
|     o  inStr to have white space removed
|   - Returns:
|     o number of characters in inStr
\-------------------------------------------------------*/
unsigned int
rmWhite_ulCp(
   signed char *inStr
){
   ulong_ulCp *cpUL = (ulong_ulCp *) inStr;
   ulong_ulCp *dupUL = (ulong_ulCp *) inStr;

   signed char *dupStr = inStr;
   signed char *cpStr = inStr;

   unsigned int uiChar = 0;

   while(*cpStr != '\0')
   { /*Loop: remove white space*/

      if( ! ifWhite_ulCp(*cpUL) )
         *dupUL++ = *cpUL++;

      else
      { /*Else: have white space, manually copy*/
         cpStr = (signed char *) cpUL;
         dupStr = (signed char *) dupUL;
         
         for(
            uiChar = 0;
            uiChar < sizeof(ulong_ulCp);
            ++uiChar
         ){ /*Loop: remove white space*/
            if(*cpStr == '\0')
               goto done_fun24;
            else if(*cpStr < 33)
               ++cpStr;
            else
               *cpStr++ = *dupStr++;
         } /*Loop: remove white space*/

         cpUL = (ulong_ulCp *) cpStr;
         dupUL = (ulong_ulCp *) dupStr;
      } /*Else: have white space, manually copy*/
   } /*Loop: remove white space*/

   done_fun24:;
      *dupStr = '\0';
      return dupStr - inStr;
} /*rmWhite_ulCp*/

/*-------------------------------------------------------\
| Fun25: swapDelim_ulCp
|   - swaps two strings until deliminator is found
| Input:
|   - firstStr:
|     o Pointer to string to first string to swap
|   - secStr:
|     o Pointer to second string to swap
|   - delimUL:
|     o delminator to end at (as long). Use makeULDelim
|       to build this deliminator
|   - delimSC:
|     o delminator (as char) to stop copying at
| Output:
|   - Modifies:
|     o firstStr to have secStr string
|     o secStr to have firstStr string
| Note:
|   - This will likely not be very good at swapping short
|     strings.
\-------------------------------------------------------*/
void
swapDelim_ulCp(
   signed char *firstStr,
   signed char *secStr,
   ulong_ulCp delimUL,
   signed char delimSC
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun15: swapDelim_ulCp
   '   - swaps two strings until deliminator is found
   '   o fun15 sec01:
   '     - variable declarations
   '   o fun15 sec02:
   '     - swap until first deliminator
   '   o fun15 sec03:
   '     - if 1st string ends early, finsh swapping second
   '   o fun15 sec04:
   '     - else 2nd string ends early, finsh swapping 1st
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ulong_ulCp *firstUL = (ulong_ulCp *) firstStr;
   ulong_ulCp *secUL = (ulong_ulCp *) secStr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec02:
   ^   - swap until first deliminator
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(
      ! (
            ifDelim_ulCp(*firstUL, delimUL)
          & ifDelim_ulCp(*secUL, delimUL) 
        )
   ){ /*Loop: Copy cpStr to dupStr*/
      *firstUL ^= *secUL;
      *secUL ^= *firstUL;
      *firstUL++ ^= *secUL++;
   } /*Loop: Copy cpStr to dupStr*/

   firstStr = (signed char *) firstUL;
   secStr = (signed char *) secUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec03:
   ^   - if first string ends early, finsh swapping second
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(ifDelim_ulCp(*firstUL, delimUL))
   { /*If: first string ended*/

      while(*firstStr != delimSC)
      { /*Loop: copy first string*/
         if(*secStr == '\0')
            break;

         *firstStr ^= *secStr;
         *secStr ^= *firstStr;
         *firstStr++ ^= *secStr++;
      } /*Loop: copy first string*/

      if(*secStr != '\0')
      { /*If: second string has more values*/
         *firstStr++ = *secStr;
         *secStr++ = '\0';
      } /*If: second string has more values*/

      firstUL = (ulong_ulCp *) firstStr;
      secUL = (ulong_ulCp *) secStr;

      while( ! ifDelim_ulCp(*secUL, delimUL) )
         *firstUL++ = *secUL++;

      firstStr = (signed char *) firstUL;
      secStr = (signed char *) secUL;

      while(*secStr != delimSC)
         *firstStr++ = *secStr++;

      *firstStr = '\0';
   } /*If: first string ended*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun15 Sec04:
   ^   - else 2nd string ends early, finsh swapping 1st
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   else
   { /*Else: second string ended*/

      while(*secStr != delimSC)
      { /*Loop: copy first string*/
         if(*firstStr == '\0')
            break;

         *firstStr ^= *secStr;
         *secStr ^= *firstStr;
         *firstStr++ ^= *secStr++;
      } /*Loop: copy first string*/

      if(*firstStr != '\0')
      { /*If: first string has more values*/
         *secStr++ = *firstStr;
         *firstStr++ = '\0';
      } /*If: first string has more values*/

      firstUL = (ulong_ulCp *) firstStr;
      secUL = (ulong_ulCp *) secStr;

      while( ! ifDelim_ulCp(*firstUL, delimUL) )
         *secUL++ = *firstUL++;

      firstStr = (signed char *) firstUL;
      secStr = (signed char *) secUL;

      while(*firstStr != delimSC)
         *secStr++ = *firstStr++;

      *secStr = '\0';
   } /*Else: second string ended*/
} /*swapDelim_ulCp*/

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
