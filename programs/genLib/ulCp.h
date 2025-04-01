/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' ulCp SOF: Start Of File
'   - uses longs to copy contents of a string
'   - these functions are slower than strcpy and strlen,
'     but they allow deliminators to be used.
'   - DOES NOT WORK ON **BIG ENDIN** SYSTEM
'   o header:
'     - Included libraries, definitions and preprocessor
'       checks
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
'      o fun26: swapNull_ulCp
'        - swaps two strings until null
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - Included libraries, definitions and preprocessor
|     checks
\-------------------------------------------------------*/

#ifndef UNSINGED_LONG_COPY_H
#define UNSINGED_LONG_COPY_H

/*plan9 defines long as 32 bit, so need an extra check*/
#ifdef PLAN9_64
   typedef unsigned long long ulong_ulCP;   
#else
   typedef unsigned long ulong_ulCp;   
#endif

/*ulong_ulCp shortcuts*/
/*the compiler should truncate these as needed*/
#define def_null_ulCp 0
#define def_newline_ulCp (ulong_ulCp)0x0a0a0a0a0a0a0a0a
#define def_carriage_ulCp (ulong_ulCp) 0x0d0d0d0d0d0d0d0d
#define def_tab_ulCp (ulong_ulCp) 0x0909090909090909
#define def_space_ulCp (ulong_ulCp) 0x2020202020202020
#define def_comma_ulCp (ulong_ulCp) 0x2c2c2c2c2c2c2c2c

#define def_one_ulCp (ulong_ulCp) 0x0101010101010101
#define def_highBit_ulCp (ulong_ulCp) 0x8080808080808080

#define def_four_ulCp (ulong_ulCp) 0x0404040404040404
#define def_eight_ulCp (ulong_ulCp) 0x0808080808080808

#define def_hi3Bits_ulCp (ulong_ulCp) 0xe0e0e0e0e0e0e0e0
   /*gets the high three bits in a char*/

/*for white space checks*/
#define def_31_ulCp (ulong_ulCp) 0x0f0f0f0f0f0f0f0f
#define def_gt63_ulCp (ulong_ulCp) 0xc0c0c0c0c0c0c0c0
   /*0xc0 is 12*/
#define def_32Or64_ulCp (ulong_ulCp) 0x6060606060606060
   /*32 | 64 = (0x20 | 0x40) = (0010 0000 | 0100 0000)*/


#define def_bitsPerChar_ulCp 8
#define def_charInUL_ulCp sizeof(ulong_ulCp)
#define def_modUL_ulCp (sizeof(ulong_ulCp ) -1)
   /*The mod values are 7 (longs), 3 (ints), and 1
   `  (shorts), so, 1 - number of bytes will give me the
   `  mod value
   */

/*This should be compliled in.
`  This will only work for up to 64 bit numbers, after
`  that it will fail. However, it is very unlikely a
`  long will by 128 bits.
*/
#define def_shiftULBy_ulCp ((sizeof(ulong_ulCp) >> 1) - (sizeof(ulong_ulCp) >> 3))

/*Logic:
`   - sizeof(ulong_ulCp) >> 1:
`     o This gets the number of bytes to shift to get
`       an 32 bit or 16 bit number to one bytes. Basically
`       it is divide by 2
`   - (sizeof(ulong_ulCp) >> 3)
`     o This adds in 1 if I have eight bytes. It is needed
`       as the 64 bit correction step
`   - (sizeof(long) >> 1) - (sizeof(long) >> 3)
`     o Gives me the number of bits I need to shift to
`       get the number of characters in an long
`   - (sizeof(long) >> 3) << (sizeof(long) >> 4)
`     o This is not used, but would correct for a 128 bit
`       number. Each time I increase the bit size I need
`       to add a correction in
*/

/*-------------------------------------------------------\
| Fun01: rshiftByte_ulCp
|   - shift a byte right by x bytes (over shifts go to 0)
| Input:
|   - val:
|     o value to shift right
|   - byte:
|     o number of bytes to shift by
| Output:
|   - Returns:
|     o shifted value
|     o 0 if over shifted
\-------------------------------------------------------*/
#define rshiftByte_ulCp(val, byte) ( ( ((ulong_ulCp) (val)) << ( ( ((ulong_ulCp) (byte)) & (-(( ((signed int) sizeof(ulong_ulCp)) - ((signed int) (byte)) ) > 0 ))) << 3 ) ) & (-(( ((signed int) sizeof(ulong_ulCp)) - ((signed int) (byte)) ) > 0 )))
/* Logic:
`    - right shift by bytes
`    - boundsCheck: (sizeof(ulong_ulCp)) - byte) ) > 0:
`      o check if shift is out of bounds
`      o 0 if out of bounds (ulong_ulCp - byte) <= 0
`      o 1 if in bounds of ulong_ulCp
`    - boundsMask: - boundsCheck: 
`      o 0 if out of bounds
`      o -1 if in of bounds
`    - getBytes: byte & boundsMask:
`      o check if keep bytes
`      o 0 if out of bounds
`      o bytes if in of bounds
`    - shiftVal: getBytes << 3:
`      o get the number of bits to shift
`      o 0 if out of bounds
`      0 number bytes * 8
`    - shift: val << shiftVal:
`      o shift value to correct position
`      0 value if out of bounds
`      0 shifted value
`    - finalVal: shift & boundsMask:
`      o sets out of bounds values to 0
`      o 0 if out of bounds
`      o shifted value if in bounds
*/

/*-------------------------------------------------------\
| Fun02: mkDelim_ulCp
|   - make ulong_ulCp delimintor from character
| Input:
|   - delim:
|     o deliminator (char); make to ulong_ulCp deliminator
| Output:
|   - Returns:
|     o ulong_ulCp as new deliminator
\-------------------------------------------------------*/
#define mkDelim_ulCp(delim) ( rshiftByte_ulCp((delim), 7) | rshiftByte_ulCp((delim), 6) | rshiftByte_ulCp((delim), 5) | rshiftByte_ulCp((delim), 4) | rshiftByte_ulCp((delim), 3) | rshiftByte_ulCp((delim), 2) | rshiftByte_ulCp((delim), 1) | rshiftByte_ulCp((delim), 0) )
/*Logic:
`  - make a ulong_ulCp delminator
`  - addDelim: rshiftByte_ulCp(delim, 7):
`    o add char deliminator to ulong_ulCp returned
`  - nextDelim: addDelim | rshiftByte_ulCp(delim, byte):
`    o add deliminator to other byte positions in
`      ulong_ulCp deliminator
`    o byte ranges from 6 to 0
*/

/*-------------------------------------------------------\
| Fun03: ifDelim_ulCp
|   - check if have deliminator in long
| Input:
|   - ulStr:
|     o ulong_ulCp (ulong or ulong long) to check for the
|       deliminator
|   - delim:
|     o delimintor in a ulong_ulCp to search for
| Output:
|   - Returns:
|     o 0 if deliminator not found
|     o not 0 if deliminator found
\-------------------------------------------------------*/
#define ifDelim_ulCp(ulStr, delim) (  ((((ulStr) ^ (delim)) - def_one_ulCp) & def_highBit_ulCp) )
/*Logic:
`  - see if deliminator is in long
`  - delimToNull: ulCp ^ delimUL:
`    o delim goes to 0
`    o non-delim values (including null) are > 0
`  - delimToNeg: delimToNull - def_one_ulCp
`    o convets delimiator to -1
`  - getNeg: delmToneg & def_highBit_ulCp:
`     o get negative flags (only in delimiator)
`     o non-delim values go to 0
*/

/*-------------------------------------------------------\
| Fun04: ifNullDelim_ulCp
|   - checks to see if null or deliminator in ulong_ulCp
| Input:
|   - ulStr:
|     o ulong_ulCp to check for null or deliminator
| Output:
|   - Returns:
|     o 0 if null or deliminator not found
|     o not 0 if null or deliminator found
\-------------------------------------------------------*/
#define ifNullDelim_ulCp(ulStr, delim) ( ( ( (ulStr) & (~(delim)) ) - def_one_ulCp ) & def_highBit_ulCp )
/*Logic
`  - checks if null or deliminator in ulong_ulCp
`  - rmDelim: ulStr & ~delim:
`    o clears all deliminator bits from ulStr
`    o null goes to 0
`    o delim values goes 0
`    o everything else is at least 1
`  - negDelim: rmDelim - def_one_ulCp
`    o 0 (deliminator / null) goes to -1
`    o everything else is >= 0
`  - getNeg: negDelim - def_highBit_ulCp
`      o clear everything except negative bit
`      o null and deliminator go to 1000 000
`      o everything else goes to 0
*/

/*-------------------------------------------------------\
| Fun05: ifNull_ulCp
|   - checks to see if null in ulong_ulCp
| Input:
|   - ulStr:
|     o ulong_ulCp (ulong or ulong long) to check for null
| Output:
|   - Returns:
|     o 0 if null not found
|     o not 0 if null found
\-------------------------------------------------------*/
#define ifNull_ulCp(ulStr) (((ulStr) - def_one_ulCp) & def_highBit_ulCp)
/* Logic
`    - finds if null is in a ulong_ulCp
`    - nullToNeg: ulStr - def_one_ulCp:
`      o null ('\0') goes to -1
`      o everything else goes to 0 or 1
`    - detectNull: nullToNeg & def_highBit_ulCp:
`      o clear everything except negative bit
`      o null goes to 1000 000
`      o everything else goes to 0
*/

/*-------------------------------------------------------\
| Fun06: ifLineUnix_ulCp
|   - checks if unix line end ('\n' or '\0')
| Input:
|   - ulStr:
|     o ulong_ulCp to check for unix end of line
| Output:
|   - Returns:
|     o 0 if unix end of line not found
|     o not 0 if unix end of line found
\-------------------------------------------------------*/
#define ifLineUnix_ulCp(ulStr) ( ( ( (ulStr) & (~def_newline_ulCp) ) - def_one_ulCp ) & def_highBit_ulCp )

/* Logic
`   - find unix line end ('\n' or '\0')
`   - blankLine: ulStr & ~def_newline_ulCp:
`     o '\n' goes to 0
`     o '\0' stays at 0
`     o everything else is > 0
`   - negLine: blankLine - def_one_upCp:
`     o converts newlines and nulls to -1
`     o all other values are >= 0
`   - findLine: negLine & def_highBit_ulCp
`     o clears all non-negative bit values
`     o '\n' and '\0' go to -127 (1000 000)
`     o evertying else goes to 0
*/

/*-------------------------------------------------------\
| Fun07: ifEndLine_ulCp
|   - checks if input long is end of line
| Input:
|   - strUL:
|     o ulong_ulCp to check for end of line
| Output:
|   - Returns:
|     o > 0 if have end of line
|     o 0 if no '\0', '\n', '\r'
| Note:
|   - asccii table values 1 to 8, 11, 12, 14, and 15 will
|     also be set to > 0, however, these asccii charaters
|     are not used in text files, so should not appear
\-------------------------------------------------------*/
#define ifEndLine_ulCp(strUL) ( (( ((strUL) & ~(def_eight_ulCp | def_four_ulCp)) + def_one_ulCp) - def_four_ulCp) & def_highBit_ulCp )

/*```````````````````````````````````````````````````````\
` Fun06 Logic:
`   - rmShared: strUL & ~(0x08... | 0x04...)
`     o removes variables making new line and \r unique
`       - only care about tab, \n, \r, \0 and values >= 32
`         (assuming all other values are never input)
`     o \0 (0) stays 0
`     o \n (8 + 2) goes to 2
`     o \r (8 + 4 + 1) goes to 1
`     o \t (8 + 2 + 1) goes to 2 + 1
`   - addOne: rmShared + 0x01...:
`     o converts tab to 4, but new line to 3
`     o \0 (0) goes to 1
`     o \n (2) goes to 2 + 1
`     o \r (1) goes to 2
`     o \t (2 + 1) goes to 4
`   - subFour: addOne - 0x04...
`     o converts \n, \r, and \0 to negative number
`     o \0 (1) goes to -3
`     o \n (2 + 1) goes to -1
`     o \r (2) goes to -2
`     o \t (4) goes to 1
`   - subFour & 0x80...:
`     o converts clears everthing except negative bit
`     o \0 (-3) goes to 0x80
`     o \n (-1) goes to 0x80
`     o \n (-2) goes to 0x80
`     o \t goes to 0
`     o 16 or higher goest to 0
\```````````````````````````````````````````````````````*/

/*-------------------------------------------------------\
| Fun08: ifWhite_ulCp
|   - check if input long has white space
|     (' ', '\r', '\n', '\t', '\0')
| Input:
|   - ulStr:
|     o ulong_ulCp (ulong or ulong long) to check for
|       white space
| Output:
|   - Returns:
|     o 0 if white space not found
|     o not 0 if white space found
\-------------------------------------------------------*/
#define ifWhite_ulCp(ulStr) ( ( (((ulStr) - def_one_ulCp) & def_hi3Bits_ulCp) - def_one_ulCp ) & def_highBit_ulCp )
/*Logic:
`   - checks if have ' ', '\t', '\n', '\r', or '\0' (< 33)
`   - cvtSpace: ulStr - def_one_ulCp:
`     o convert space to 31 [no 32 flag] (null goes to -1)
`     o values > 32 (non-white) still have 32 flag
`   - rmLowBit: cvtSpace & def_hi3Bits_ulCp:
`     o removes bits < 32/64/128 (1, 2, 4, 8, 16 flags)
`   - whiteToNeg: rmLowBit - def_one_ulCp:
`     o converts white space to negative -1
`     o null goes to -2
`     o everything else goes to 0 or higher
`   - getNeg: rmLowBit & def_highBit_ulCp:
`     o get negative flags (only in white space)
`     o non-white space (>32) goes to 0
*/


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
);

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
);


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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

/*-------------------------------------------------------\
| Fun26: swapNull_ulCp
|   - swaps two strings until null
| Input:
|   - firstStr:
|     o Pointer to string to first string to swap
|   - secStr:
|     o Pointer to second string to swap
| Output:
|   - Modifies:
|     o firstStr to have secStr string
|     o secStr to have firstStr string
| Note:
|   - This will likely not be very good at swapping short
|     strings.
\-------------------------------------------------------*/
void
swapNull_ulCp(
   signed char *firstStr,
   signed char *secStr
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
