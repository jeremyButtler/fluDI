/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' revNtTo2Bit SOF: Start Of File
'   - had tables to convert nucleodtides to two bits
'   - the third bit holds a error bit that marks anonymous
'     base or non-nucleotide characters. it can be ingored
'     after checks
'   o header:
'     - guards and defined variables
'   o tbl01: revNtTo2Bit
'     - table to convert bases to four bit (a/t/g/c only)
'       values to two bits with the 3rd bit saved for an
'       error/anonymous base
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards and defined variables
\-------------------------------------------------------*/

#ifndef REVERSE_NUCLEOTIDE_TO_TWO_BIT_H
#define REVERSE_NUCLEOTIDE_TO_TWO_BIT_H

#define def_t_revNtTo2Bit 0
#define def_c_revNtTo2Bit 1
#define def_a_revNtTo2Bit 2
#define def_g_revNtTo2Bit 3
#define def_err3rdBit_revNtTo2Bit 4

/*-------------------------------------------------------\
| Tbl01: revNtTo2Bit
|  - table to convert bases to four bit (a/t/g/c only)
|    values to two bits with the 3rd bit saved for an
|    error/anonymous base
|    o 1st bit is A (a_fourBit)
|    o 2nd bit is C (c_fourBit)
|    o 3rd bit is G (g_fourBit)
|    o 4th bit is T/U (t_fourBit)
\-------------------------------------------------------*/
static
unsigned char revNtTo2Bit[] =
{  /*revNtTo2Bit*/
   /*White space/invisible charactes block*/
   def_err3rdBit_revNtTo2Bit, /*0   = Null character*/

   def_err3rdBit_revNtTo2Bit, /*1   = Start of Heading*/
   def_err3rdBit_revNtTo2Bit, /*2   = Start of Text*/
   def_err3rdBit_revNtTo2Bit, /*3   = End of Text*/
   def_err3rdBit_revNtTo2Bit, /*4   = End Transmission*/
   def_err3rdBit_revNtTo2Bit, /*5   = Enquiry*/
   def_err3rdBit_revNtTo2Bit, /*6   = Acknowledge*/
   def_err3rdBit_revNtTo2Bit, /*7   = Bell*/
   def_err3rdBit_revNtTo2Bit, /*8   = Backspace*/

   def_err3rdBit_revNtTo2Bit, /*9   =  tab (horizontal)*/
   def_err3rdBit_revNtTo2Bit, /*10  = New line*/

   def_err3rdBit_revNtTo2Bit, /*11  = Vertical Tab*/
   def_err3rdBit_revNtTo2Bit, /*12  = Form Feed*/

   def_err3rdBit_revNtTo2Bit, /*13  = Carriage Return*/

   def_err3rdBit_revNtTo2Bit, /*14  = Shift Out*/
   def_err3rdBit_revNtTo2Bit, /*15  = Shift In*/
   def_err3rdBit_revNtTo2Bit, /*16  = Data Link Escape*/
   def_err3rdBit_revNtTo2Bit, /*17  = Device Control One*/
   def_err3rdBit_revNtTo2Bit, /*18  = Device Control Two*/
   def_err3rdBit_revNtTo2Bit, /*19  = Device Contol 3*/
   def_err3rdBit_revNtTo2Bit, /*20  = Device Control 4*/
   def_err3rdBit_revNtTo2Bit, /*21  = Neg Acknowledge*/
   def_err3rdBit_revNtTo2Bit, /*22  = Synchronous Idle*/
   def_err3rdBit_revNtTo2Bit, /*23  = End Transmission*/
   def_err3rdBit_revNtTo2Bit, /*24  = Cancle*/
   def_err3rdBit_revNtTo2Bit, /*25  = End of medium*/
   def_err3rdBit_revNtTo2Bit, /*26  = Substitute*/
   def_err3rdBit_revNtTo2Bit, /*27  = escape*/
   def_err3rdBit_revNtTo2Bit, /*28  = File Separator*/
   def_err3rdBit_revNtTo2Bit, /*29  = Group Separator*/
   def_err3rdBit_revNtTo2Bit, /*30  = Record Separator*/
   def_err3rdBit_revNtTo2Bit, /*31  = Unit Separator*/

   /*symbol/number block*/
   def_err3rdBit_revNtTo2Bit, /*32  = space*/
   def_err3rdBit_revNtTo2Bit, /*33  = !*/
   def_err3rdBit_revNtTo2Bit, /*34  = "*/
   def_err3rdBit_revNtTo2Bit, /*35  = #*/
   def_err3rdBit_revNtTo2Bit, /*36  = $*/
   def_err3rdBit_revNtTo2Bit, /*37  = %*/
   def_err3rdBit_revNtTo2Bit, /*38  = &*/
   def_err3rdBit_revNtTo2Bit, /*39  = '*/
   def_err3rdBit_revNtTo2Bit, /*40  = (*/
   def_err3rdBit_revNtTo2Bit, /*41  = )*/
   def_err3rdBit_revNtTo2Bit, /*42  = **/
   def_err3rdBit_revNtTo2Bit, /*43  = +*/
   def_err3rdBit_revNtTo2Bit, /*44  = ,*/
   def_err3rdBit_revNtTo2Bit, /*45  = -*/
   def_err3rdBit_revNtTo2Bit, /*46  = .*/
   def_err3rdBit_revNtTo2Bit, /*47  = /*/
   def_err3rdBit_revNtTo2Bit, /*48  = 0*/
   def_err3rdBit_revNtTo2Bit, /*49  = 1*/
   def_err3rdBit_revNtTo2Bit, /*50  = 2*/
   def_err3rdBit_revNtTo2Bit, /*51  = 3*/
   def_err3rdBit_revNtTo2Bit, /*52  = 4*/
   def_err3rdBit_revNtTo2Bit, /*53  = 5*/
   def_err3rdBit_revNtTo2Bit, /*54  = 6*/
   def_err3rdBit_revNtTo2Bit, /*55  = 7*/
   def_err3rdBit_revNtTo2Bit, /*56  = 8*/
   def_err3rdBit_revNtTo2Bit, /*57  = 9*/
   def_err3rdBit_revNtTo2Bit, /*58  = :*/
   def_err3rdBit_revNtTo2Bit, /*59  = ;*/
   def_err3rdBit_revNtTo2Bit, /*60  = <*/
   def_err3rdBit_revNtTo2Bit, /*61  = =*/
   def_err3rdBit_revNtTo2Bit, /*62  = >*/
   def_err3rdBit_revNtTo2Bit, /*63  = ?*/
   def_err3rdBit_revNtTo2Bit, /*64  = @*/

   /*Uppercase letters*/
   def_t_revNtTo2Bit,   /*65  = A*/

   def_err3rdBit_revNtTo2Bit, /*66=B; C/G/T*/

   def_g_revNtTo2Bit,   /*67  = C*/


   def_err3rdBit_revNtTo2Bit, /*68=D; A/G/T*/
   def_err3rdBit_revNtTo2Bit, /*69  = E not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*70  = F not nucleotide*/

   def_c_revNtTo2Bit,     /*71  = G*/

   def_err3rdBit_revNtTo2Bit, /*72=H; A/G/T*/
   def_err3rdBit_revNtTo2Bit, /*73  = I not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*74  = J not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*75=K; G/T*/
   def_err3rdBit_revNtTo2Bit, /*76  = L not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*77=M; A/C*/

   def_err3rdBit_revNtTo2Bit, /*78=N; any base*/
   def_err3rdBit_revNtTo2Bit, /*79  = O not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*80  = P not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*81  = Q not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*82=R A/G*/
   def_err3rdBit_revNtTo2Bit, /*83=S C/G*/

   def_a_revNtTo2Bit,   /*84  = T*/
   def_a_revNtTo2Bit,   /*85  = U*/

   def_err3rdBit_revNtTo2Bit, /*86  = V (ACG),treat as N*/
   def_err3rdBit_revNtTo2Bit, /*87=W A\T*/
   def_err3rdBit_revNtTo2Bit, /*88=X; any aa*/
   def_err3rdBit_revNtTo2Bit,   /*89=Y C/T*/
   def_err3rdBit_revNtTo2Bit, /*90  = Z not nucleotide*/

   /*Special characters after upercase letters*/
   def_err3rdBit_revNtTo2Bit, /*91  = [*/
   def_err3rdBit_revNtTo2Bit, /*92  = \*/
   def_err3rdBit_revNtTo2Bit, /*93  = ]*/
   def_err3rdBit_revNtTo2Bit, /*94  = ^*/
   def_err3rdBit_revNtTo2Bit, /*95  = _*/
   def_err3rdBit_revNtTo2Bit, /*96  = `*/

   /*lower case letters*/
   def_t_revNtTo2Bit,   /*97=a*/

   def_err3rdBit_revNtTo2Bit, /*98=b; C/G/T*/

   def_g_revNtTo2Bit,   /*99  = c*/

   def_err3rdBit_revNtTo2Bit, /*100=d; A/G/T*/
   def_err3rdBit_revNtTo2Bit, /*101 = e not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*102 = f not nucleotide*/
   def_c_revNtTo2Bit,   /*103 = g*/
   def_err3rdBit_revNtTo2Bit, /*104=h; A/G/T*/
   def_err3rdBit_revNtTo2Bit, /*105 = i not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*106 = j not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*107=k; G/T*/
   def_err3rdBit_revNtTo2Bit, /*108  = l not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*109=m; A/C*/

   def_err3rdBit_revNtTo2Bit, /*110=n; any base*/

   def_err3rdBit_revNtTo2Bit, /*111 = o not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*112 = p not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*113 = q not nucleotide*/
   def_err3rdBit_revNtTo2Bit, /*114=r A/G*/
   def_err3rdBit_revNtTo2Bit, /*115=s C/G*/
   def_a_revNtTo2Bit,     /*116 = t*/
   def_a_revNtTo2Bit,     /*117 = u*/

   
   def_err3rdBit_revNtTo2Bit, /*118 = v (ACG),treat as N*/
   def_err3rdBit_revNtTo2Bit, /*119=w A\T*/
   def_err3rdBit_revNtTo2Bit, /*120 = x*/
   def_err3rdBit_revNtTo2Bit,   /*121=y C/T*/
   def_err3rdBit_revNtTo2Bit, /*122 = z not nucleotide*/

   /*Special characters after lowercase letters*/
   def_err3rdBit_revNtTo2Bit, /*123 = {*/
   def_err3rdBit_revNtTo2Bit, /*124 = |*/
   def_err3rdBit_revNtTo2Bit, /*125 = }*/
   def_err3rdBit_revNtTo2Bit, /*126 = ~*/
   def_err3rdBit_revNtTo2Bit, /*127 = Del*/
}; /*revNtTo2Bit*/

#endif

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the publikc domain
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
