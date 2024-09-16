/*#########################################################
# Name: codonTbl
#    - has codon table and functions to support the
#      codon table
#########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - guards
'   o tbl01 codonLkTbl:
'     - Table to convert three bases to codons
'   o fun01: codonToAA_codonTbl
'     - converts codon into an amino acid
'   o fun02: revCodonToAA_codonTbl
'     - reverses complements input codon & gets amino acid
'   o fun03: bacStartCode_codonTbl
'     - checks to see if input bases is a bacterial start
'       codon (convert nt with ntTo2Bit.h)
'   o fun04: bacStart_codonTbl
'     - checks if codon (no conversion) is a start codon
'   o fun05: bactRevStart_codonTbl
'     - checks if codon (no conversion) is a start codon
'       for a reverse complement sequence
'   o fun06: aaTripToChar_codonTbl
'     - converts a three letter amino acid idenity to its
'       single letter amino acid identity
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards
\-------------------------------------------------------*/

#ifndef CODON_TABLE_H
#define CODON_TABLE_H

/*--------------------------------------------------------\
| Tbl01 codonTbl:
|  - table to convert three bases to codons
\--------------------------------------------------------*/
static signed char codonTbl[5][5][5] =
   {  /*codonTbl*/
      { /*First element is an T*/
         {'f', 'f', 'l', 'l', 'x'}, /*2nd element is T*/
         {'s', 's', 's', 's', 's'}, /*2nd element is C*/
         {'y', 'y', '*', '*', 'x'}, /*2nd element is A*/
         {'c', 'c', '*', 'w', 'x'}, /*2nd element is G*/
         {'x', 'x', 'x', 'x', 'x'}  /*2nd element is N*/
      }, /*First element is an T*/

      { /*First element is an C*/
         {'l', 'l', 'l', 'l', 'l'}, /*2nd element is T*/
         {'p', 'p', 'p', 'p', 'p'}, /*2nd element is C*/
         {'h', 'h', 'q', 'q', 'x'}, /*2nd element is A*/
         {'r', 'r', 'r', 'r', 'r'}, /*2nd element is G*/
         {'x', 'x', 'x', 'x', 'x'}  /*2nd element is N*/
      }, /*First element is an C*/

      { /*First element is an A*/
         {'i', 'i', 'i', 'm', 'x'}, /*2nd element is T*/
         {'t', 't', 't', 't', 't'}, /*2nd element is C*/
         {'n', 'n', 'k', 'k', 'x'}, /*2nd element is A*/
         {'s', 's', 'r', 'r', 'x'}, /*2nd element is G*/
         {'x', 'x', 'x', 'x', 'x'}  /*2nd element is N*/
      }, /*First element is an A*/

      { /*First element is an G*/
         {'v', 'v', 'v', 'v', 'v'}, /*2nd element is T*/
         {'a', 'a', 'a', 'a', 'a'}, /*2nd element is C*/
         {'d', 'd', 'e', 'e', 'x'}, /*2nd element is A*/
         {'g', 'g', 'g', 'g', 'g'}, /*2nd element is G*/
         {'x', 'x', 'x', 'x', 'x'}  /*2nd element is N*/
      }, /*First element is an G*/

      { /*First element is an N*/
         {'x', 'x', 'x', 'x', 'x'}, /*2nd element is T*/
         {'x', 'x', 'x', 'x', 'x'}, /*2nd element is C*/
         {'x', 'x', 'x', 'x', 'x'}, /*2nd element is A*/
         {'x', 'x', 'x', 'x', 'x'}, /*2nd element is G*/
         {'x', 'x', 'x', 'x', 'x'}  /*2nd element is N*/
      }, /*First element is an N*/
   }; /*codonTbl*/

/*-------------------------------------------------------\
| Fun01: codonToAA_codonTbl
|   - converts codon into an amino acid
| Input:
|   - firstNtSC:
|     o first base in codon
|   - secNtSC:
|     o second base in codon
|   - thirdNtSC:
|     o last base in codon
| Output:
|   - Returns:
|     o amino acid of the input codon
\-------------------------------------------------------*/
signed char
codonToAA_codonTbl(
   signed char firstNtSC,
   signed char secNtSC,
   signed char thirdNtSC
);

/*-------------------------------------------------------\
| Fun02: revCodonToAA_codonTbl
|   - reverses complements input codon & finds amino acid
| Input:
|   - firstNtSC:
|     o first base in the foward codon, but third (last)
|       base in the reverse complement codon
|   - secNtSC:
|     o second base in forward codon, but second (middle)
|       base in the reverse complement codon
|   - thirdNtSC:
|     o last base in forward codon, but first base in the
|       reverse complement codon
| Output:
|   - Returns:
|     o reverse complement amino acid of the input
|       codon
\-------------------------------------------------------*/
signed char
revCodonToAA_codonTbl(
   signed char firstNtSC,
   signed char secNtSC,
   signed char thirdNtSC
);

/*-------------------------------------------------------\
| Fun03: bacStartCode_codonTbl
|   - checks to see if input bases is a bacterial start
|     codon (convert nt with ntTo2Bit.h)
| Input:
|   - firstNtSC:
|     o first base to check; it needs to be the output of
|       an base to code look up table
|   - secBaseC:
|     o second base to check; it needs to be the output
|       of an base to code look up table
|   - thirdBaseC:
|     o third base to check; it needs to be the output
|       of an base to code look up table
| Output:
|   - Returns:
|     o 1 if this is an bacterial start codon
|     o 0 if this is not an bacterial start codon
\-------------------------------------------------------*/
signed char
bactStartCode_codonTbl(
   unsigned char firstNtSC,
   unsigned char secNtSC,
   unsigned char thirdNtSC
);

/*-------------------------------------------------------\
| Fun04: bacStart_codonTbl
|   - checks if codon (no conversion) is a start codon
| Input:
|   - firstNtSC:
|     o first base to check
|   - secBaseC:
|     o second base to check
|   - thirdBaseC:
|     o third base to check
| Output:
|   - Returns:
|     o 1 if this is a bacterial start codon
|     o 0 if this is not a bacterial start codon
\-------------------------------------------------------*/
signed char
bactStart_codonTbl(
   signed char firstNtSC,
   signed char secNtSC,
   signed char thirdNtSC
);

/*-------------------------------------------------------\
| Fun05: bactRevStart_codonTbl
|   - checks if codon (no conversion) is a start codon
|     for a reverse complement sequence
| Input:
|   - firstNtSC:
|     o first base to check (last in codon)
|   - secBaseC:
|     o second base to check (middle in codon)
|   - thirdBaseC:
|     o third base to check (first in codon)
| Output:
|   - Returns:
|     o 1 if this is a reverse bacterial start codon
|     o 0 if this is not a reverse bacterial start codon
\-------------------------------------------------------*/
signed char
bactRevStart_codonTbl(
   signed char firstNtSC,
   signed char secNtSC,
   signed char thirdNtSC
);

/*-------------------------------------------------------\
| Fun06: aaTripToChar_codonTbl
|   - converts a three letter amino acid idenity to its
|     single letter amino acid identity
| Input:
|   - codonStr:
|     o c-string with the three letter amino acid to
|       convert
| Output:
|   - Returns:
|     o one letter code for the amino acid
|     o 0 for an invalid code
\-------------------------------------------------------*/
signed char
aaTripToChar_codonTbl(
   signed char *codonStr
);

#endif

/*
   Ala A Alanine
   Arg R Arginine
   Asn N Asparagine
   Asp D Aspartic acid
   Cys C Cysteine
   Gln Q Glutamine
   Glu E Glutamic acid
   Gly G Glycine
   His H Histidine
   Ile I Isoleucine
   Lys K Lysine
   Leu L Luecine
   Met M Methionine
   Phe F Phenylalanine
   Pro P Proline
   Ser S Serine
   Thr T Threonine
   Trp W Tryptophan
   Try Y Tyrosine
   Val V Valine
*/

/*Table
     T         C        A        G
  +----------+----------+----------+----------+
  | TTT  F . | TCT  S . | TAT  Y . | TGT  C . | T
T | TTC  F . | TCC  S . | TAC  Y . | TGC  C . | C
  | TTA  L . | TCA  S . | TAA  * . | TGA  * . | A
  | TTG  L 2 | TCG  S . | TAG  * . | TGG  W . | G
  +----------+----------+----------+----------+
  | CTT  L . | CCT  P . | CAT  H . | CGT  R . | T
C | CTC  L . | CCC  P . | CAC  H . | CGC  R . | C
  | CTA  L . | CCA  P . | CAA  Q . | CGA  R . | A
  | CTG  L . | CCG  P . | CAG  Q . | CGG  R . | G
  +----------+----------+----------+----------+
  | ATT  I . | ACT  T . | AAT  N . | AGT  S . | T
A | ATC  I . | ACC  T . | AAC  N . | AGC  S . | C
  | ATA  I . | ACA  T . | AAA  K . | AGA  R . | A
  | ATG  M 1 | ACG  T . | AAG  K . | AGG  R . | G
  +----------+----------+----------+----------+
  | GTT  V . | GCT  A . | GAT  D . | GGT  G . | T
G | GTC  V . | GCC  A . | GAC  D . | GGC  G . | C
  | GTA  V . | GCA  A . | GAA  E . | GGA  G . | A
  | GTG  V 2 | GCG  A . | GAG  E . | GGG  G . | G
  +----------+----------+----------+----------+

  1 = Canonical start codon
  2 = Bacterial non-canconical start codon
  . = Nothing
*/

/*ASCII table
   Dec Hex  Binary   Value  Description
   0    00  00000000 NUL    Null character
   1    01  00000001 SOH    Start of Heading
   2    02  00000010 STX    Start of Text
   3    03  00000011 ETX    End of Text
   4    04  00000100 EOT    End of Transmission
   5    05  00000101 ENQ    Enquiry
   6    06  00000110 ACK    Acknowledge
   7    07  00000111 BEL    Bell, Alert
   8    08  00001000 BS     Backspace
   9    09  00001001 HT     Horizontal Tab
   10   0A  00001010 LF     Line Feed
   11   0B  00001011 VT     Vertical Tabulation
   12   0C  00001100 FF     Form Feed
   13   0D  00001101 CR     Carriage Return
   14   0E  00001110 SO     Shift Out
   15   0F  00001111 SI     Shift In
   16   10  00010000 DLE    Data Link Escape
   17   11  00010001 DC1    Device Control One (XON)
   18   12  00010010 DC2    Device Control Two
   19   13  00010011 DC3    Device Control Three (XOFF)
   20   14  00010100 DC4    Device Control Four
   21   15  00010101 NAK    Negative Acknowledge
   22   16  00010110 SYN    Synchronous Idle
   23   17  00010111 ETB    End of Transmission Block
   24   18  00011000 CAN    Cancel
   25   19  00011001 EM     End of medium
   26   1A  00011010 SUB    Substitute
   27   1B  00011011 ESC    Escape
   28   1C  00011100 FS     File Separator
   29   1D  00011101 GS     Group Separator
   30   1E  00011110 RS     Record Separator
   31   1F  00011111 US     Unit Separator

   32   20  00100000 SP     Space
   33   21  00100001 !      Exclamation mark
   34   22  00100010 "      Double quotes
   35   23  00100011 #      Number sign
   36   24  00100100 $      Dollar
   37   25  00100101 %      Per cent sign
   38   26  00100110 &      Ampersand
   39   27  00100111 '      Single quote
   40   28  00101000 (      Open parenthesis
   41   29  00101001 )      Close parenthesis
   42   2A  00101010 *      Asterisk
   43   2B  00101011 +      Plus
   44   2C  00101100 ,      Comma
   45   2D  00101101 -      Hyphen-minus
   46   2E  00101110 .      Period, dot or full stop
   47   2F  00101111 /      Slash or divide
   48   30  00110000 0      Zero
   49   31  00110001 1      One
   50   32  00110010 2      Two
   51   33  00110011 3      Three
   52   34  00110100 4      Four
   53   35  00110101 5      Five
   54   36  00110110 6      Six
   55   37  00110111 7      Seven
   56   38  00111000 8      Eight
   57   39  00111001 9      Nine
   58   3A  00111010 :      Colon
   59   3B  00111011 ;      Semicolon
   60   3C  00111100 <      Less than
   61   3D  00111101 =      Equals
   62   3E  00111110 >      Greater than
   63   3F  00111111 ?      Question mark
   64   40  01000000 @      At sign
   65   41  01000001 A      Uppercase A
   66   42  01000010 B      Uppercase B
   67   43  01000011 C      Uppercase C
   68   44  01000100 D      Uppercase D
   69   45  01000101 E      Uppercase E
   70   46  01000110 F      Uppercase F
   71   47  01000111 G      Uppercase G
   72   48  01001000 H      Uppercase H
   73   49  01001001 I      Uppercase I
   74   4A  01001010 J      Uppercase J
   75   4B  01001011 K      Uppercase K
   76   4C  01001100 L      Uppercase L
   77   4D  01001101 M      Uppercase M
   78   4E  01001110 N      Uppercase N
   79   4F  01001111 O      Uppercase O
   80   50  01010000 P      Uppercase P
   81   51  01010001 Q      Uppercase Q
   82   52  01010010 R      Uppercase R
   83   53  01010011 S      Uppercase S
   84   54  01010100 T      Uppercase T
   85   55  01010101 U      Uppercase U
   86   56  01010110 V      Uppercase V
   87   57  01010111 W      Uppercase W
   88   58  01011000 X      Uppercase X
   89   59  01011001 Y      Uppercase Y
   90   5A  01011010 Z      Uppercase Z
   91   5B  01011011 [      Opening bracket
   92   5C  01011100 \      Backslash
   93   5D  01011101 ]      Closing bracket
   94   5E  01011110 ^      Caret - circumflex
   95   5F  01011111 _      Underscore
   96   60  01100000 `      Grave accent
   97   61  01100001 a      Lowercase a
   98   62  01100010 b      Lowercase b
   99   63  01100011 c      Lowercase c
   100  64  01100100 d      Lowercase d
   101  65  01100101 e      Lowercase e
   102  66  01100110 f      Lowercase f
   103  67  01100111 g      Lowercase g
   104  68  01101000 h      Lowercase h
   105  69  01101001 i      Lowercase i
   106  6A  01101010 j      Lowercase j
   107  6B  01101011 k      Lowercase k
   108  6C  01101100 l      Lowercase l
   109  6D  01101101 m      Lowercase m
   110  6E  01101110 n      Lowercase n
   111  6F  01101111 o      Lowercase o
   112  70  01110000 p      Lowercase p
   113  71  01110001 q      Lowercase q
   114  72  01110010 r      Lowercase r
   115  73  01110011 s      Lowercase s
   116  74  01110100 t      Lowercase t
   117  75  01110101 u      Lowercase u
   118  76  01110110 v      Lowercase v
   119  77  01110111 w      Lowercase w
   120  78  01111000 x      Lowercase x
   121  79  01111001 y      Lowercase y
   122  7A  01111010 z      Lowercase z
   123  7B  01111011 {      Opening brace
   124  7C  01111100 |      Vertical bar
   125  7D  01111101 }      Closing brace
   126  7E  01111110 ~      Equivalency sign - tilde
   127  7F  01111111 DEL    Delete

There is an extended ascii table, but I do not use it
  here since fasta/fastq use regular ascii
*/

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
