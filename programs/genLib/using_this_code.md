# Use:

Here to give you an orientation of the C code files and
  functions.

# Libraries:

- Table Of Contents (TOC)
  - [base10str](#base10str)
    - c-string to base10 number conversoin
  - [charCp](#charcp)
    - byte copy functions
  - [dataTypeShortHand](#datatypeshorthand)
    - shorthand values for datatypes
  - [genMath](#genmath)
    - math (branchless) functions
  - [numToStr](#numtostr)
    - number to c-string
  - [ptrAry](#ptrary)
    - c-string array (using pointers)
  - [shellSort](#shellsort)
    - shell short and their functions
  - [strAry](#strary)
    - c-string array (no pointers, fixed, 60 byte width)
  - [ulCp](#ulcp)
    - long byte copy and string manipulation functions

I gues for me libraries include .c and .h files and are
  not limited to .so/.dylib files. These are sometimes
  redunent, but often are more nuianced then the base
  C code.

To find a list of functions see the SOF (Start Of File)
  block (really a table of contens, but already use TOC)
  at the start of each file. They should be mostly up to
  date.

## base10str

base10str is here to convert c-strings into nubmers. Has
  both a .c and .h file.

Functions are in form of strTo<DataType>_base10str. With
  <DataType> being the data type to convert to. Input is
  a c-string and a pointer to the data type to convert.

| short hand | data Type      | function          |
|:----------:|:---------------|:------------------|
|     UL     | unsigned long  | strToUL_base10str |
|     UI     | unsigned int   | strToUI_base10str |
|     US     | unsigned short | strToUS_base10str |
|     UC     | unsigned char  | strToUC_base10str |
|     SL     | signed long    | strToSL_base10str |
|     SI     | signed int     | strToSI_base10str |
|     SS     | signed short   | strToSS_base10str |
|     SC     | signed char    | strToSC_base10str |

Table: of base10str functions.


## charCp

Has functions for manipulating c-strings. Not very
  usefull, since they work on one byte at a time.

Has both .c and .h files.

## dataTypeShortHand:

Has my shorthand notations for data types. Mainly because
  I do not like to type `unsigned` or `signed` all the
  time. These mirror the plan9 short hands.

.h file only.

## genMath

These are my general math functions. Most are branchless
  functions from the stanfford bithacking guide. Has
  min, max, abloxlute values (SL for long, regular int).
  Also includes and ifMax/ifMin function to find the
  maximum/minimum and copy its paired input
  (is branchless).

The min/max ifMin/ifMax exist are .h file macros. The
  other functions are in the .c file.

## numToStr:

Converts unsigned numbers (numToStr) or doubles
  (double_numToStr) to c-strings. For doubles, negatives
  values are allowed.

Has both a .c and a .h file.

## ptrAry

Depends on ulCp.h/ulCp.c

Has my functions to manage a string pointer array.
  Includes setup (Fun05 mk_str_ptrAry), adding strings
  (Fun07 add_str_ptrAry), sorting strings
  (Fun10 sort_str_ptrAry), finding strings from a
  sorted array (Fun12 find_str_ptrAry), and comparing
  strings in the array (Fun08 cmp_str_ptrAry).

Inorder to use this you first must create a str_ptrAry
  structure with mk_str_ptrAry (provide starting array
  size). Then you add strings with add_str_ptrAry (memory
  is reallocated as needed). Finally, when finished, you
  can free the array with freeHeap_str_ptrAry (fun04).

## shellSort

Has shell sorting functions. There are a lot (one for
  each data type). The pattern is <dataType>_shellSort to
  sort, <dataType>Search_shellSort to search for an item
  in a sorted array, and <dataType>Range_shellSort to
  find the first value in a range in a sorted array. Each
  type of function takes the same input, but is designed
  for a paticular datatype.

For example, the unsigned long sort is ul_shellSort(),
  while the signed long sort is sl_shellSort(). To find
  and item in a unsigned long array I would do
  ulSearch_shellSort();.

Shell sort has both a .c file and .h file. I had an older
  verison wich was just macros, but it turns out those
  relied on GNU extensions, and so were not c89
  compatable. That was a reall pain and I miss the
  convience of macros.

## strAry

Depends on ulCp.h/ulCp.c

This is a set of functions to manage a string array in
  a single buffer (no pointers) [think matrix]. Each
  string is 64 bytes longs.

You can make a c-string buffer with mk_strAry(size)
  (Fun01). You can then add strings with add_strAry
  (Fun04) (WARNING DOES NOT RELLOC FOR YOU). To add more
  memory use realloc_strAry(). You can get pointers to
  strings with get_strAry (Fun03; .h macro). You can sort
  strings with sort_strAry (Fun07) and find strings in the
  sorted array with find_strAry (Fun09). To compare
  strings use cmp_strAry (Fun06).

The string is justs a c-string, so use free() to free you
  strAry.

## ulCp

I have no idea if this will work on big endin systems. Not
  a big deal since most user OSs are little endin. I have
  tested one or two functions using qemu to emulate a
  mk68040 processor on a debian image I found online
  ([https://people.debian.org/~gio/dqib/](
   https://people.debian.org/~gio/dqib/)
  ). It seems to work, but then I also had mixed results
  with a browser emulator.

ulCp is short for unsigned long copy. It is like charCp,
  except that it uses unsigned long. As a general rule it
  is slower thang strcpy, which uses SIMD. However, it
  is often faster than charCp or other byte by byte
  methods (at least on 64 bit). It also allows for you to
  use your own deliminators for some functions. Use
  mkDelim_ulCp (fun02) to make a deliminator (store it in
  an ulong_ulCp dataType).

I had to define my own data type for ulCp because plan9
  treats longs as 32 bits. So, their is a special flag
  (-DPLAN9_64) which will use `long long` as ulong_ulCp to
  get 64 bit dataTypes (only use on 64 bit systems).

I have defined several deliminators at the top of the
  ulCp.h file. To see a list of functions vew the table
  of contents (Start Of File/SOF) at the start of ulCp.c
  or ulCp.h. I define 64 bit values and then use casting
  to resize it to 32 bit or whatever size is needed.

| deliminator     | defined macro        |
|:----------------|:---------------------|
| null            | def_null_ulCp (or 0) |
| new line        | def_newline_ulCp     |
| carriage return | def_carriage_ulCp    |
| tab             | def_tab_ulCp         |
| space           | def_space_ulCp       |
| comma           | def_comma_ulCp       |

Table: of deliminators in ulCp.h

There are also several long check macros to see if a value
  is in a long (fun03 to fun08 in ulCp.h).

| deliminator     | function          | function id |
|:----------------|:-----------------:|:------------|
| user defined    | ifDelim_ulCp      | Fun03       |
| null + user     | ifDelimNull_ulCp  | Fun04       |
| null ('\\0')    | ifNull_ulCp       | Fun05       |
| newline + null  | ifLineUnix_ulCp   | Fun06       |
| \\n + \\r + \\0 | ifEndLine_ulCp    | Fun07       |
| white space     | ifWhite_ulCp      | Fun08       |

Table: of macros that check a long for certain charaters.
  for user or user defined deliminators, use a pre-defined
  deliminator or `ulong_ulCp delim = mkDelim_ulCp(value)`
  to make a deliminator. By white space I mean any
  asccii value less then 33 (this includes space, tab,
  carriage return, newline, and null).
