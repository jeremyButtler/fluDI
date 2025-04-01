/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' kmerBit SOF: Start Of File
'   - converts nucleotide index (from seqToIndex_alnSet)
'     to a kmer bit
'   o .h tbl01: alnNtTo_kmerBit
'     - converts an nucleotide alignment code form
'       alnSetStruct.h to an two bit value, with an 3rd
'       bit being used for anonymous bases and the 4th bit
'       for errors
'   o fun01: mkMask_kmerBit
'     - makes a mask for clearing extra bits from a kmer
'   o fun02: ntIndexToKmer_kmerBit
'     - adds a nucleotide (in index format
'      [seqToIndex_alnSet]) to a kmer
'   o fun03: ntBitToKmer_kmerBit
'     - adds a nucleotide bit (by alNtTo_kmerBit) to kmer
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef KMER_BIT_H
#define KMER_BIT_H

#define def_bitsPerKmer_kmerBit 2 /*do not change*/

#define def_a_kmerBit 0
#define def_c_kmerBit 1
#define def_g_kmerBit 2
#define def_t_kmerBit 3
#define def_anonNt_kmerBit 4 /*anonmyous nucleotide*/
#define def_noNt_kmerBit 8

/*-------------------------------------------------------\
| Tbl01: alnNtTo_kmerBit
|   - converts an nucleotide alignment code form
|     alnSetStruct.h to an two bit value, with an 3rd
|     bit being used for anonymous bases and the 4th bit
|     for errors
\-------------------------------------------------------*/
static unsigned char alnNtTo_kmerBit[] =
{
   def_noNt_kmerBit,    /*00 = 64 = ,*/  
   def_a_kmerBit,       /*01 = 65 = A/a*/  
   def_anonNt_kmerBit,  /*02 = 66 = B/b*/  
   def_c_kmerBit,       /*03 = 67 = C/c*/  
   def_anonNt_kmerBit,  /*04 = 68 = D/d*/  
   def_noNt_kmerBit,    /*05 = 69 = E/e*/  
   def_noNt_kmerBit,    /*06 = 70 = F/f*/  
   def_g_kmerBit,       /*07 = 71 = G/g*/  
   def_anonNt_kmerBit,  /*08 = 72 = H/h*/  
   def_noNt_kmerBit,    /*09 = 73 = I/i*/  
   def_noNt_kmerBit,    /*10 = 74 = J/j*/  
   def_anonNt_kmerBit,  /*11 = 75 = K/k*/  
   def_noNt_kmerBit,    /*12 = 76 = L/l*/  
   def_anonNt_kmerBit,  /*13 = 77 = M/m*/  
   def_anonNt_kmerBit,  /*14 = 78 = N/n*/  
   def_noNt_kmerBit,    /*15 = 79 = O/o*/  
   def_noNt_kmerBit,    /*16 = 80 = P/p*/  
   def_noNt_kmerBit,    /*17 = 81 = Q/q*/  
   def_anonNt_kmerBit,  /*18 = 82 = R/r*/  
   def_anonNt_kmerBit,  /*19 = 83 = S/s*/  
   def_t_kmerBit,       /*20 = 84 = T/t*/  
   def_t_kmerBit,       /*21 = 85 = U/u*/  
   def_anonNt_kmerBit,  /*22 = 86 = V/v*/  
   def_anonNt_kmerBit,  /*23 = 87 = W/w*/  
   def_anonNt_kmerBit,  /*24 = 88 = X/x (amino acids)*/  
   def_anonNt_kmerBit,  /*25 = 89 = Y/y*/  
   def_noNt_kmerBit     /*26 = 90 = Z/z*/  
}; /*alnNtTo_kmerBit*/


/*-------------------------------------------------------\
| Fun01: mkMask_kmerBit
|   - makes a mask for clearing extra bits from a kmer
| Input:
|   - lenKmerMac:
|     o number nucleotides in one kmer
| Output:
|   - Returns:
|     o mask to clear kmer (0's for unused bits and 1's
|       for used bits)
\-------------------------------------------------------*/
#define mkMask_kmerBit(lenKmerMac) ( ((unsigned long) -1) >> (( sizeof(unsigned long) << 3 ) - ( (lenKmerMac) * def_bitsPerKmer_kmerBit )) )
/*Logic:
`  - kmerBits: lenKmerMac * def_bitsPerKmer_kmerBit:
`    o gets number of bits needed to fill one kmer
`  - sizeLong: sizeof(unsigned long) << 3
`    o number of bits in unsigned long
`  - emptyBits: sizeLong - kmerBits:
`    o number of extra bits in unsigned long (not used)
`  - ( (unsinged long) -1 ) >> emptyBits
`    o shifts -1 till all unsed bits are set to 0, and all
`      used bits are set to 1
`     
*/

/*-------------------------------------------------------\
| Fun02: ntIndexToKmer_kmerBit
|   - adds a nucleotide (in index format
|    [seqToIndex_alnSet]) to a kmer
| Input:
|   - ntMac:
|     o nucleotide to add
|   - kmerMac:
|     o kmer with previous nucleotides
|   - maskMac:
|     o mask to clear extr bits with (use mkMask_kmerBit)
| Output:
|   - Returns:
|     o kmer with added bits (does not modify kmerMac)
\-------------------------------------------------------*/
#define ntIndexToKmer_kmerBit(ntMac, kmerMac, maskMac) ( ( ((kmerMac) << def_bitsPerKmer_kmerBit) | (alnNtTo_kmerBit[ (unsigned char) (ntMac) ]) ) & (maskMac) )
/*Logic:
`   - mkRoom: kmerMac << def_bitsPerKmer_kmerBit:
`     o add room for new bit
`   - getBit: alnNtTo_kmerBit[(unsigned char) ntMac]:
`     o convert nucleotide index to bits for kmer
`   - kmer: mkRoom | getBit:
`     o add converted nucleotide index to kmer
`   - kmer & mascMac:
`     o clear any extra bits (not in kmer)
*/

/*-------------------------------------------------------\
| Fun03: ntBitToKmer_kmerBit
|   - adds a nucleotide bit (by alNtTo_kmerBit) to a kmer
| Input:
|   - ntBitMac:
|     o nucleotide bit to add; covert by alnNtTo_kmerBit
|   - kmerMac:
|     o kmer with previous nucleotides
|   - maskMac:
|     o mask to clear extr bits with (use mkMask_kmerBit)
| Output:
|   - Returns:
|     o kmer with added bits (does not modify kmerMac)
\-------------------------------------------------------*/
#define ntBitToKmer_kmerBit(ntBitMac, kmerMac, maskMac) ( ( ((kmerMac) << def_bitsPerKmer_kmerBit) | ((unsigned char) (ntBitMac)) ) & (maskMac) )
/*Logic:
`   - mkRoom: kmerMac << def_bitsPerKmer_kmerBit:
`     o add room for new bit
`   - kmer: mkRoom | ntBitMac:
`     o add converted nucleotide index to kmer
`   - kmer & mascMac:
`     o clear any extra bits (not in kmer)
*/

#endif
