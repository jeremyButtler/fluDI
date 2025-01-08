/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' rmBlocks SOF: Start Of File
'   - has functions to remove blocks (start ---, end ---)
'     from file and merge block contents onto a line
'   - also removes all comments
'     - formats: #.*\n; #---.*---#; or #---.*---#
'   o header:
'     - guards
'   o fun01: str_rmBlocks
'     - converts blocks in string (---) to single lines
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - gaurds
\-------------------------------------------------------*/

#ifndef REMOVE_BLOCKS_AND_COMMENTS_H
#define REMOVE_BLOCKS_AND_COMMENTS_H

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
);

#endif

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
