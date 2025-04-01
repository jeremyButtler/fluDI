/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' diScore SOF: Start Of File:
'   - sets scoring matrix up for diFrag
'   o header:
'     - included libraries and defaults
'   o fun01: set_diScore
'     - sets up scoring matrix for diFrag
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries and defaults
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif 

#include "diScore.h"
#include "../genAln/alnSet.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries:
!   - .c  #include "../genLib/base10str.h"
!   - .c  #include "../genLib/ulCp.h"
!   - .h  #include "alnDefs.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: set_diScore
|   - sets up scoring matrix for diFrag
| Input:
|   - alnSetSTPtr:
|     o pointer to alnSet struct to setup
| Output:
|   - Modifies:
|     o scoreMatrixSS in alnSetSTPtr to have all
|       non-anoymous matches and snps to def_match_diScore
|       and def_snp_diScore
\-------------------------------------------------------*/
void
set_diScore(
   struct alnSet *alnSetSTPtr
){
   setScore_alnSet('A','A',def_match_diScore,alnSetSTPtr);
   setScore_alnSet('C','C',def_match_diScore,alnSetSTPtr);
   setScore_alnSet('G','G',def_match_diScore,alnSetSTPtr);
   setScore_alnSet('T','T',def_match_diScore,alnSetSTPtr);
   setScore_alnSet('U','U',def_match_diScore,alnSetSTPtr);

   setScore_alnSet('A','C',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('A','G',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('A','T',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('A','U',def_snp_diScore,alnSetSTPtr);

   setScore_alnSet('C','A',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('C','G',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('C','T',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('C','U',def_snp_diScore,alnSetSTPtr);

   setScore_alnSet('G','A',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('G','C',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('G','T',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('G','U',def_snp_diScore,alnSetSTPtr);

   setScore_alnSet('T','A',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('T','C',def_snp_diScore,alnSetSTPtr);
   setScore_alnSet('T','G',def_snp_diScore,alnSetSTPtr);

   changeGap_alnSet(
      alnSetSTPtr,
      def_open_diScore,
      def_extend_diScore
   );
} /*set_diScore*/
