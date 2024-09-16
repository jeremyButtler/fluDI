/*########################################################
# Name: ntTo5Bit
#   - has look up table to convert nucleodtides to a five
#     bit value
#   - the six bit is an error bit for non-nucleotide
#     characters. it can be ignored after checking.
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   o header:
'     - guards and defined variables
'   o tbl01: ntTo5Bit
'     - table to convert bases to five bit values, with an
'       extra sixth bit acting as an error
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards and defined variables
\-------------------------------------------------------*/

#ifndef NUCLEOTIDE_TO_FIVE_BIT_H
#define NUCLEOTIDE_TO_FIVE_BIT_H

#define def_a_fourBit_ntTo5Bit 1
#define def_c_fourBit_ntTo5Bit 1 << 1
#define def_g_fourBit_ntTo5Bit 1 << 2
#define def_t_fourBit_ntTo5Bit 1 << 3
#define def_n_fithBit_ntTo5Bit 1 << 4 /*any anonymous nt*/
#define def_err_sixBit_ntTo5Bit 1 << 5

/*-------------------------------------------------------\
| Tbl-01: ntTo5Bit
|  - table to convert bases to five bit values (a/c/g/t/
|    is anonymous), with an extra sixth bit acting as an
|    error flag
|    o 1st bit is A (def_a_fourBit_ntTo5Bit)
|    o 2nd bit is C (def_c_fourBit_ntTo5Bit)
|    o 3rd bit is G (def_g_fourBit_ntTo5Bit)
|    o 4th bit is T/U (def_t_fourBit_ntTo5Bit)
|    o 5th bit marks anonymous base (def_n_fithBit_ntTo5Bit)
|    o 6th bit is error (def_err_sixBit_ntTo5Bit)
\-------------------------------------------------------*/
static
unsigned char ntTo5Bit[] =
{ /*ntTo5Bit*/
   /*White space/invisible charactes block*/
   def_err_sixBit_ntTo5Bit, /*0   = Null character*/

   def_err_sixBit_ntTo5Bit, /*1   = Start of Heading*/
   def_err_sixBit_ntTo5Bit, /*2   = Start of Text*/
   def_err_sixBit_ntTo5Bit, /*3   = End of Text*/
   def_err_sixBit_ntTo5Bit, /*4   = End of Transmission*/
   def_err_sixBit_ntTo5Bit, /*5   = Enquiry*/
   def_err_sixBit_ntTo5Bit, /*6   = Acknowledge*/
   def_err_sixBit_ntTo5Bit, /*7   = Bell*/
   def_err_sixBit_ntTo5Bit, /*8   = Backspace*/

   def_err_sixBit_ntTo5Bit, /*9   =  tab (horizontal)*/
   def_err_sixBit_ntTo5Bit, /*10  = New line*/

   def_err_sixBit_ntTo5Bit, /*11  =Vertical Tab (not key)*/
   def_err_sixBit_ntTo5Bit, /*12  = Form Feed*/

   def_err_sixBit_ntTo5Bit, /*13  = Carriage Return*/

   def_err_sixBit_ntTo5Bit, /*14  = Shift Out*/
   def_err_sixBit_ntTo5Bit, /*15  = Shift In*/
   def_err_sixBit_ntTo5Bit, /*16  = Data Link Escape*/
   def_err_sixBit_ntTo5Bit, /*17  = Device Control One*/
   def_err_sixBit_ntTo5Bit, /*18  = Device Control Two*/
   def_err_sixBit_ntTo5Bit, /*19  = Device Contol Three*/
   def_err_sixBit_ntTo5Bit, /*20  = Device Control Four*/
   def_err_sixBit_ntTo5Bit, /*21  = Negative Acknowledge*/
   def_err_sixBit_ntTo5Bit, /*22  = Synchronous Idle*/
   def_err_sixBit_ntTo5Bit, /*23  =End Transmission Block*/
   def_err_sixBit_ntTo5Bit, /*24  = Cancle*/
   def_err_sixBit_ntTo5Bit, /*25  = End of medium*/
   def_err_sixBit_ntTo5Bit, /*26  = Substitute*/
   def_err_sixBit_ntTo5Bit, /*27  = escape*/
   def_err_sixBit_ntTo5Bit, /*28  = File Separator*/
   def_err_sixBit_ntTo5Bit, /*29  = Group Separator*/
   def_err_sixBit_ntTo5Bit, /*30  = Record Separator*/
   def_err_sixBit_ntTo5Bit, /*31  = Unit Separator*/

   /*symbol/number block*/
   def_err_sixBit_ntTo5Bit, /*32  = space*/
   def_err_sixBit_ntTo5Bit, /*33  = !*/
   def_err_sixBit_ntTo5Bit, /*34  = "*/
   def_err_sixBit_ntTo5Bit, /*35  = #*/
   def_err_sixBit_ntTo5Bit, /*36  = $*/
   def_err_sixBit_ntTo5Bit, /*37  = %*/
   def_err_sixBit_ntTo5Bit, /*38  = &*/
   def_err_sixBit_ntTo5Bit, /*39  = '*/
   def_err_sixBit_ntTo5Bit, /*40  = (*/
   def_err_sixBit_ntTo5Bit, /*41  = )*/
   def_err_sixBit_ntTo5Bit, /*42  = **/
   def_err_sixBit_ntTo5Bit, /*43  = +*/
   def_err_sixBit_ntTo5Bit, /*44  = ,*/
   def_err_sixBit_ntTo5Bit, /*45  = -*/
   def_err_sixBit_ntTo5Bit, /*46  = .*/
   def_err_sixBit_ntTo5Bit, /*47  = /*/
   def_err_sixBit_ntTo5Bit, /*48  = 0*/
   def_err_sixBit_ntTo5Bit, /*49  = 1*/
   def_err_sixBit_ntTo5Bit, /*50  = 2*/
   def_err_sixBit_ntTo5Bit, /*51  = 3*/
   def_err_sixBit_ntTo5Bit, /*52  = 4*/
   def_err_sixBit_ntTo5Bit, /*53  = 5*/
   def_err_sixBit_ntTo5Bit, /*54  = 6*/
   def_err_sixBit_ntTo5Bit, /*55  = 7*/
   def_err_sixBit_ntTo5Bit, /*56  = 8*/
   def_err_sixBit_ntTo5Bit, /*57  = 9*/
   def_err_sixBit_ntTo5Bit, /*58  = :*/
   def_err_sixBit_ntTo5Bit, /*59  = ;*/
   def_err_sixBit_ntTo5Bit, /*60  = <*/
   def_err_sixBit_ntTo5Bit, /*61  = =*/
   def_err_sixBit_ntTo5Bit, /*62  = >*/
   def_err_sixBit_ntTo5Bit, /*63  = ?*/
   def_err_sixBit_ntTo5Bit, /*64  = @*/

   /*Uppercase letters*/
   def_a_fourBit_ntTo5Bit,   /*65  = A*/

   (
       def_n_fithBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*66=B; C/G/T*/

   def_c_fourBit_ntTo5Bit,   /*67  = C*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*68=D; A/G/T*/

   def_err_sixBit_ntTo5Bit, /*69  = E not nucleotide*/
   def_err_sixBit_ntTo5Bit, /*70  = F not nucleotide*/
   def_g_fourBit_ntTo5Bit,   /*71  = G*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*72=H; A/G/T*/

   def_err_sixBit_ntTo5Bit, /*73  = I not nucleotide*/
   def_err_sixBit_ntTo5Bit, /*74  = J not nucleotide*/

   (
       def_n_fithBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*75=K; G/T*/

   def_err_sixBit_ntTo5Bit, /*76  = L not nucleotide*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
   ), /*77=M; A/C*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*78=N; any base*/

   def_err_sixBit_ntTo5Bit, /*79  = O not nucleotide*/
   def_err_sixBit_ntTo5Bit, /*80  = P not nucleotide*/
   def_err_sixBit_ntTo5Bit, /*81  = Q not nucleotide*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
   ), /*82=R A/G*/

   (
       def_n_fithBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
   ), /*83=S C/G*/

   def_t_fourBit_ntTo5Bit,   /*84  = T*/
   def_t_fourBit_ntTo5Bit,   /*85  = U*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
   ), /*86  = V (ACG), treat as N*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*87=W A\T*/

   ( 
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*88=X; any aa (I will hold same for nucelotides)*/

   (
       def_n_fithBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ),   /*89=Y C/T*/

   def_err_sixBit_ntTo5Bit, /*90  = Z not nucleotide*/

   /*Special characters after upercase letters*/
   def_err_sixBit_ntTo5Bit, /*91  = [*/
   def_err_sixBit_ntTo5Bit, /*92  = \*/
   def_err_sixBit_ntTo5Bit, /*93  = ]*/
   def_err_sixBit_ntTo5Bit, /*94  = ^*/
   def_err_sixBit_ntTo5Bit, /*95  = _*/
   def_err_sixBit_ntTo5Bit, /*96  = `*/

   /*lower case letters*/
   def_a_fourBit_ntTo5Bit,   /*97=a*/

   (
       def_n_fithBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*98=b; C/G/T*/

   def_c_fourBit_ntTo5Bit,   /*99  = c*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*100=d; A/G/T*/

   def_err_sixBit_ntTo5Bit, /*101 = e not nucleotide*/
   def_err_sixBit_ntTo5Bit, /*102 = f not nucleotide*/
   def_g_fourBit_ntTo5Bit,   /*103 = g*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*104=h; A/G/T*/

   def_err_sixBit_ntTo5Bit, /*105 = i not nucleotide*/
   def_err_sixBit_ntTo5Bit, /*106 = j not nucleotide*/

   (
       def_n_fithBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*107=k; G/T*/

   def_err_sixBit_ntTo5Bit, /*108  = l not nucleotide*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
   ), /*109=m; A/C*/

   /*110=n; any base*/
   (  
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*Any base (n)*/

   def_err_sixBit_ntTo5Bit, /*111 = o not nucleotide*/
   def_err_sixBit_ntTo5Bit, /*112 = p not nucleotide*/
   def_err_sixBit_ntTo5Bit, /*113 = q not nucleotide*/

   (
        def_n_fithBit_ntTo5Bit
      | def_a_fourBit_ntTo5Bit
      | def_g_fourBit_ntTo5Bit
   ), /*114=r A/G*/

   (
        def_n_fithBit_ntTo5Bit
      | def_c_fourBit_ntTo5Bit
      | def_g_fourBit_ntTo5Bit
   ),/*115=s C/G*/

   def_t_fourBit_ntTo5Bit,   /*116 = t*/
   def_t_fourBit_ntTo5Bit,   /*117 = u*/

   (
       def_n_fithBit_ntTo5Bit 
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
   ), /*118 = v (ACG), treat as N*/

   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*119=w A\T*/

   /*120=x; any aa (I will hold same for nucelotides)*/
   (
       def_n_fithBit_ntTo5Bit
     | def_a_fourBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_g_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*x*/

   (
       def_n_fithBit_ntTo5Bit
     | def_c_fourBit_ntTo5Bit
     | def_t_fourBit_ntTo5Bit
   ), /*121=y C/T*/

   def_err_sixBit_ntTo5Bit, /*122 = z not nucleotide*/

   /*Special characters after lowercase letters*/
   def_err_sixBit_ntTo5Bit, /*123 = {*/
   def_err_sixBit_ntTo5Bit, /*124 = |*/
   def_err_sixBit_ntTo5Bit, /*125 = }*/
   def_err_sixBit_ntTo5Bit, /*126 = ~*/
   def_err_sixBit_ntTo5Bit, /*127 = Del*/
}; /*ntTo5Bit*/

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
