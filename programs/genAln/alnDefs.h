/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' alnDefs SOF: Start Of File
'   - holds default settings and global definitions for
'     sequence alignment
'   o header:
'     - guards
'   o sec02:
'     - matrix variables (DO NOT CHANGE)
'   o sec03:
'     - score variables (gap, ajdustments, and general)
'   o fun01: maxScore_alnDefs
'     - finds maximum possible score for a sequence
'       (assumes all matches are a score of 5)
'   o sec04:
'     - score matrix defaults (nucleotide only)
'   o sec05:
'     - match matrix defaults (nucleotide only)
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards
\-------------------------------------------------------*/

#ifndef ALIGNMENT_SEQUENCES_DEFAULTS_H
#define ALIGNMENT_SEQUENCES_DEFAULTS_H

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec01:
^   - errors and return variables
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define def_memErr_alnDefs 1

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec02:
^   - matrix variables (DO NOT CHANGE)
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*Alignment matrix movements*/
/*Do not change these values*/
#define def_mvStop_alnDefs 0  /*Stop*/
#define def_mvDel_alnDefs 1   /*Move left (deletion)*/
#define def_mvIns_alnDefs 3   /*Move up (insertion)*/
#define def_mvSnp_alnDefs 2   /*Move diagnol (snp/match)*/
  /*so a gap can be found by & 1*/

/*for matching matrix (each flag is a different bit)*/
#define def_ntEql_alnDefs 1
#define def_ntNotEql_alnDefs 2

#define def_anonymous_alnDefs 4

#define def_anonMatch_alnDefs (def_ntEql_alnDefs | def_anonymous_alnDefs)

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec03:
^   - score variables (gap, ajdustments, and general)
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*Scoring variables*/
#define def_gapOpen_alnDefs -1000
   /*Penalty for starting indel*/

#define def_gapExtend_alnDefs -100
   /*Penalty for extading an indel*/

#define def_scoreAdj_alnDefs 100
   /*adjust score for eDNA matrix. 100 allows -1 to be
   `   0.01
   */

/*-------------------------------------------------------\
| Fun01: maxScore_alnDefs
|   - finds maximum possible score for a sequence
|     (assumes all matches are a score of 5)
| Input:
|   - lenReadUI:
|     o length of sequence
| Output:
|   - Returns:
|     o maximum possible score
\-------------------------------------------------------*/
/*quick way to get maximum possible score for a read*/
#define maxScore_alnDefs(lenReadUI) ((lenReadUI) * 5 * def_scoreAdj_alnDefs)

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec04:
^   - score matrix defaults (nucleotide only)
^   o sec04 sub01:
^     - a as first base
^   o sec04 sub02:
^     - t as first base
^   o sec04 sub03:
^     - u (t) as first base
^   o sec04 sub04:
^     - g as first base
^   o sec04 sub05:
^     - c as first base
^   o sec04 sub06:
^     - w (anonymous) as first base
^   o sec04 sub07:
^     - s (anonymous) as first base
^   o sec04 sub08:
^     - m (anonymous) as first base
^   o sec04 sub09:
^     - k (anonymous) as first base
^   o sec04 sub10:
^     - r (anonymous) as first base
^   o sec04 sub11:
^     - y (anonymous) as first base
^   o sec04 sub12:
^     - b (anonymous) as first base
^   o sec04 sub13:
^     - d (anonymous) as first base
^   o sec04 sub14:
^     - h (anonymous) as first base
^   o sec04 sub15:
^     - v (anonymous) as first base
^   o sec04 sub16:
^     - n (anonymous) as first base
^   o sec04 sub17:
^     - x (anonymous) as first base (technically aa)
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/* Scoring matrix is EDNAFULL or something close to it. I
`   think I have the anonymous bases with the correct
`   scores.
*/

/*Alignment scoring matrix (non-anonymous); this is
'  in multiples of 100 to allow for a 0.01 gap extend
'  penalty
*/

/********************************************************\
* Sec04 Sub01:
*   - a as first base
\********************************************************/

#define def_AToA_alnDefs 500    /*score: A and A*/
#define def_AToT_alnDefs -400   /*score: A and T*/
#define def_AToU_alnDefs -400   /*score: A and T*/
#define def_AToG_alnDefs -400   /*score: A and G*/
#define def_AToC_alnDefs -400   /*score: A and C*/

#define def_AToW_alnDefs 100    /*score: A and W (AT)*/
#define def_AToS_alnDefs -400   /*score: A and S (CG)*/
#define def_AToM_alnDefs 100    /*score: A and M (AC)*/
#define def_AToK_alnDefs -400   /*score: A and K (GT)*/
#define def_AToR_alnDefs 100    /*score: A and R (AG)*/
#define def_AToY_alnDefs -400   /*score: A and Y (CT)*/
#define def_AToB_alnDefs -400   /*score: A and B (CGT)*/
#define def_AToD_alnDefs -100   /*score: A and D (AGT)*/
#define def_AToH_alnDefs -100   /*score: A and H (ACT)*/
#define def_AToV_alnDefs -100   /*score: A and V (ACG)*/
#define def_AToN_alnDefs -200   /*score: A and N (ACGT)*/
#define def_AToX_alnDefs -200   /*score: A and X (ACGT)*/

/********************************************************\
* Sec04 Sub02:
*   - t as first base
\********************************************************/

#define def_TToA_alnDefs -400   /*score: T and A*/
#define def_TToT_alnDefs 500    /*score: T and T*/
#define def_TToU_alnDefs 500    /*score: T and T*/
#define def_TToG_alnDefs -400   /*score: T and G*/
#define def_TToC_alnDefs -400   /*score: T and C*/

#define def_TToW_alnDefs 100    /*score: T and W (AT)*/
#define def_TToS_alnDefs -400   /*score: T and S (CG)*/
#define def_TToM_alnDefs -400   /*score: T and M (AC)*/
#define def_TToK_alnDefs 100    /*score: T and K (GT)*/
#define def_TToR_alnDefs -400   /*score: T and R (AG)*/
#define def_TToY_alnDefs 100    /*score: T and Y (CT)*/
#define def_TToB_alnDefs -100   /*score: T and B (CGT)*/
#define def_TToD_alnDefs -100   /*score: T and D (AGT)*/
#define def_TToH_alnDefs -100   /*score: T and H (ACT)*/
#define def_TToV_alnDefs -400   /*score: T and V (ACG)*/
#define def_TToN_alnDefs -200   /*score: T and N (ACGT)*/
#define def_TToX_alnDefs -200   /*score: T and X (ACGT)*/

/********************************************************\
* Sec04 Sub03:
*   - u as first base
\********************************************************/

#define def_UToA_alnDefs -400   /*score: U and A*/
#define def_UToT_alnDefs 500    /*score: U and T*/
#define def_UToU_alnDefs 500    /*score: U and U*/
#define def_UToG_alnDefs -400   /*score: U and G*/
#define def_UToC_alnDefs -400   /*score: U and C*/

#define def_UToW_alnDefs 100    /*score: T and W (AT)*/
#define def_UToS_alnDefs -400   /*score: T and S (CG)*/
#define def_UToM_alnDefs -400   /*score: T and M (AC)*/
#define def_UToK_alnDefs 100    /*score: T and K (GT)*/
#define def_UToR_alnDefs -400   /*score: T and R (AG)*/
#define def_UToY_alnDefs 100    /*score: T and Y (CT)*/
#define def_UToB_alnDefs -100   /*score: T and B (CGT)*/
#define def_UToD_alnDefs -100   /*score: T and D (AGT)*/
#define def_UToH_alnDefs -100   /*score: T and H (ACT)*/
#define def_UToV_alnDefs -400   /*score: T and V (ACG)*/
#define def_UToN_alnDefs -200   /*score: T and N (ACGT)*/
#define def_UToX_alnDefs -200   /*score: T and X (ACGT)*/

/********************************************************\
* Sec04 Sub04:
*   - g as first base
\********************************************************/

#define def_GToA_alnDefs -400   /*score: G and A*/
#define def_GToT_alnDefs -400   /*score: G and T*/
#define def_GToU_alnDefs -400   /*score: G and T*/
#define def_GToG_alnDefs 500    /*score: G and G*/
#define def_GToC_alnDefs -400   /*score: G and C*/

#define def_GToW_alnDefs -400   /*score: G and W (AT)*/
#define def_GToS_alnDefs 100    /*score: G and S (CG)*/
#define def_GToM_alnDefs -400   /*score: G and M (AC)*/
#define def_GToK_alnDefs 100    /*score: G and K (GT)*/
#define def_GToR_alnDefs 100    /*score: G and R (AG)*/
#define def_GToY_alnDefs -400   /*score: G and Y (CT)*/
#define def_GToB_alnDefs -100   /*score: G and B (CGT)*/
#define def_GToD_alnDefs -100   /*score: G and D (AGT)*/
#define def_GToH_alnDefs -400   /*score: G and H (ACT)*/
#define def_GToV_alnDefs -100   /*score: G and V (ACG)*/
#define def_GToN_alnDefs -200   /*score: G and N (ACGT)*/
#define def_GToX_alnDefs -200   /*score: G and X (ACGT)*/

/********************************************************\
* Sec04 Sub05:
*   - c as first base
\********************************************************/

#define def_CToA_alnDefs -400   /*score: C and A*/
#define def_CToT_alnDefs -400   /*score: C and T*/
#define def_CToU_alnDefs -400   /*score: C and T*/
#define def_CToG_alnDefs -400   /*score: C and G*/
#define def_CToC_alnDefs 500    /*score: C and C*/

#define def_CToW_alnDefs -400   /*score: C and W (AT)*/
#define def_CToS_alnDefs 100    /*score: C and S (CG)*/
#define def_CToM_alnDefs 100    /*score: C and M (AC)*/
#define def_CToK_alnDefs -400   /*score: C and K (GT)*/
#define def_CToR_alnDefs -400   /*score: C and R (AG)*/
#define def_CToY_alnDefs 100    /*score: C and Y (CT)*/
#define def_CToB_alnDefs -100   /*score: C and B (CGT)*/
#define def_CToD_alnDefs -400   /*score: C and D (AGT)*/
#define def_CToH_alnDefs -100   /*score: C and H (ACT)*/
#define def_CToV_alnDefs -100   /*score: C and V (ACG)*/
#define def_CToN_alnDefs -200   /*score: C and N (ACGT)*/
#define def_CToX_alnDefs -200   /*score: C and X (ACGT)*/

/********************************************************\
* Sec04 Sub06:
*   - w as first base
\********************************************************/

#define def_WToA_alnDefs 100    /*score: W and A (AT)*/
#define def_WToT_alnDefs 100    /*score: T and W (AT)*/
#define def_WToU_alnDefs 100    /*score: T and W (AT)*/
#define def_WToG_alnDefs -400   /*score: G and W (AT)*/
#define def_WToC_alnDefs -400   /*score: W and C (AT)*/

#define def_WToW_alnDefs -100   /*score: W and W (AT)*/
#define def_WToS_alnDefs -400   /*score: W and S (CG)*/
#define def_WToM_alnDefs -100   /*score: W and M (AC)*/
#define def_WToK_alnDefs -100   /*score: W and K (GT)*/
#define def_WToR_alnDefs -100   /*score: W and R (AG)*/
#define def_WToY_alnDefs -100   /*score: W and Y (CT)*/
#define def_WToB_alnDefs -100   /*score: W and B (CGT)*/
#define def_WToD_alnDefs -100   /*score: W and D (AGT)*/
#define def_WToH_alnDefs -100   /*score: W and H (ACT)*/
#define def_WToV_alnDefs -100   /*score: W and V (ACG)*/
#define def_WToN_alnDefs -100   /*score: W and N (ACGT)*/
#define def_WToX_alnDefs -100   /*score: W and X (ACGT)*/

/********************************************************\
* Sec04 Sub07:
*   - s as first base
\********************************************************/

#define def_SToA_alnDefs -400   /*score: S and A (CG)*/
#define def_SToT_alnDefs -400   /*score: T and S (CG)*/
#define def_SToU_alnDefs -400   /*score: T and S (CG)*/
#define def_SToG_alnDefs 100    /*score: G and S (CG)*/
#define def_SToC_alnDefs 100    /*score: S and C (CG)*/

#define def_SToW_alnDefs -400   /*score: S and W (AT)*/
#define def_SToS_alnDefs -100   /*score: S and S (CG)*/
#define def_SToM_alnDefs -200   /*score: S and M (AC)*/
#define def_SToK_alnDefs -200   /*score: S and K (GT)*/
#define def_SToR_alnDefs -200   /*score: S and R (AG)*/
#define def_SToY_alnDefs -200   /*score: S and Y (CT)*/
#define def_SToB_alnDefs -100   /*score: S and B (CGT)*/
#define def_SToD_alnDefs -300   /*score: S and D (AGT)*/
#define def_SToH_alnDefs -300   /*score: S and H (ACT)*/
#define def_SToV_alnDefs -100   /*score: S and V (ACG)*/
#define def_SToN_alnDefs -100   /*score: S and N (ACGT)*/
#define def_SToX_alnDefs -100   /*score: S and X (ACGT)*/

/********************************************************\
* Sec04 Sub08:
*   - m as first base
\********************************************************/

#define def_MToA_alnDefs 100    /*score: M and A (AC)*/
#define def_MToT_alnDefs -400   /*score: T and M (AC)*/
#define def_MToU_alnDefs -400   /*score: T and M (AC)*/
#define def_MToG_alnDefs -400   /*score: G and M (AC)*/
#define def_MToC_alnDefs 100    /*score: M and C (AC)*/

#define def_MToW_alnDefs -200   /*score: M and W (AT)*/
#define def_MToS_alnDefs -200   /*score: M and S (CG)*/
#define def_MToM_alnDefs -100   /*score: M and M (AC)*/
#define def_MToK_alnDefs -400   /*score: M and K (GT)*/
#define def_MToR_alnDefs -200   /*score: M and R (AG)*/
#define def_MToY_alnDefs -200   /*score: M and Y (CT)*/
#define def_MToB_alnDefs -300   /*score: M and B (CGT)*/
#define def_MToD_alnDefs -300   /*score: M and D (AGT)*/
#define def_MToH_alnDefs -100   /*score: M and H (ACT)*/
#define def_MToV_alnDefs -100   /*score: M and V (ACG)*/
#define def_MToN_alnDefs -100   /*score: M and N (ACGT)*/
#define def_MToX_alnDefs -100   /*score: M and X (ACGT)*/

/********************************************************\
* Sec04 Sub09:
*   - k as first base
\********************************************************/

#define def_KToA_alnDefs -400   /*score: K and A (GT)*/
#define def_KToT_alnDefs 100    /*score: T and K (GT)*/
#define def_KToU_alnDefs 100    /*score: T and K (GT)*/
#define def_KToG_alnDefs 100    /*score: G and K (GT)*/
#define def_KToC_alnDefs -400   /*score: K and C (GT)*/

#define def_KToW_alnDefs -200   /*score: K and W (AT)*/
#define def_KToS_alnDefs -200   /*score: K and S (CG)*/
#define def_KToM_alnDefs -400   /*score: K and M (AC)*/
#define def_KToK_alnDefs -100   /*score: K and K (GT)*/
#define def_KToR_alnDefs -200   /*score: K and R (AG)*/
#define def_KToY_alnDefs -200   /*score: K and Y (CT)*/
#define def_KToB_alnDefs -100   /*score: K and B (CGT)*/
#define def_KToD_alnDefs -100   /*score: K and D (AGT)*/
#define def_KToH_alnDefs -300   /*score: K and H (ACT)*/
#define def_KToV_alnDefs -300   /*score: K and V (ACG)*/
#define def_KToN_alnDefs -100   /*score: K and N (ACGT)*/
#define def_KToX_alnDefs -100   /*score: K and X (ACGT)*/

/********************************************************\
* Sec04 Sub10:
*   - r as first base
\********************************************************/

#define def_RToA_alnDefs 100    /*score: R and A (AG)*/
#define def_RToT_alnDefs -400   /*score: T and R (AG)*/
#define def_RToU_alnDefs -400   /*score: T and R (AG)*/
#define def_RToG_alnDefs 100    /*score: G and R (AG)*/
#define def_RToC_alnDefs -400   /*score: R and C (AG)*/

#define def_RToW_alnDefs -100   /*score: R and W (AT)*/
#define def_RToS_alnDefs -100   /*score: R and S (CG)*/
#define def_RToM_alnDefs -100   /*score: R and M (AC)*/
#define def_RToK_alnDefs -200   /*score: R and K (GT)*/
#define def_RToR_alnDefs -100   /*score: R and R (AG)*/
#define def_RToY_alnDefs -400   /*score: R and Y (CT)*/
#define def_RToB_alnDefs -300   /*score: R and B (CGT)*/
#define def_RToD_alnDefs -100   /*score: R and D (AGT)*/
#define def_RToH_alnDefs -100   /*score: R and H (ACT)*/
#define def_RToV_alnDefs -300   /*score: R and V (ACG)*/
#define def_RToN_alnDefs -100   /*score: R and N (ACGT)*/
#define def_RToX_alnDefs -100   /*score: R and X (ACGT)*/

/********************************************************\
* Sec04 Sub11:
*   - y as first base
\********************************************************/

#define def_YToA_alnDefs -400   /*score: Y and A (CT)*/
#define def_YToT_alnDefs -400   /*score: T and Y (CT)*/
#define def_YToU_alnDefs -400   /*score: T and Y (CT)*/
#define def_YToG_alnDefs -400   /*score: G and Y (CT)*/
#define def_YToC_alnDefs 100    /*score: Y and C (CT)*/

#define def_YToW_alnDefs -200   /*score: Y and W (AT)*/
#define def_YToS_alnDefs -200   /*score: Y and S (CG)*/
#define def_YToM_alnDefs -200   /*score: Y and M (AC)*/
#define def_YToK_alnDefs -200   /*score: Y and K (GT)*/
#define def_YToR_alnDefs -400   /*score: Y and R (AG)*/
#define def_YToY_alnDefs -100   /*score: Y and Y (CT)*/
#define def_YToB_alnDefs -200   /*score: Y and B (CGT)*/
#define def_YToD_alnDefs -300   /*score: Y and D (AGT)*/
#define def_YToH_alnDefs -100   /*score: Y and H (ACT)*/
#define def_YToV_alnDefs -300   /*score: Y and V (ACG)*/
#define def_YToN_alnDefs -100   /*score: Y and N (ACGT)*/
#define def_YToX_alnDefs -100   /*score: Y and X (ACGT)*/

/********************************************************\
* Sec04 Sub12:
*   - b as first base
\********************************************************/

#define def_BToA_alnDefs -400   /*score: B and A (CGT)*/
#define def_BToT_alnDefs -100   /*score: T and B (CGT)*/
#define def_BToU_alnDefs -100   /*score: T and B (CGT)*/
#define def_BToG_alnDefs -100   /*score: G and B (CGT)*/
#define def_BToC_alnDefs -100   /*score: B and C (CGT)*/

#define def_BToW_alnDefs -300   /*score: B and W (AT)*/
#define def_BToS_alnDefs -100   /*score: B and S (CG)*/
#define def_BToM_alnDefs -300   /*score: B and M (AC)*/
#define def_BToK_alnDefs -100   /*score: B and K (GT)*/
#define def_BToR_alnDefs -300   /*score: B and R (AG)*/
#define def_BToY_alnDefs -100   /*score: B and Y (CT)*/
#define def_BToB_alnDefs -100   /*score: B and B (CGT)*/
#define def_BToD_alnDefs -200   /*score: B and D (AGT)*/
#define def_BToH_alnDefs -200   /*score: B and H (ACT)*/
#define def_BToV_alnDefs -200   /*score: B and V (ACG)*/
#define def_BToN_alnDefs -100   /*score: B and N (ACGT)*/
#define def_BToX_alnDefs -100   /*score: B and X (ACGT)*/

/********************************************************\
* Sec04 Sub13:
*   - d as first base
\********************************************************/

#define def_DToA_alnDefs -100   /*score: D and A (AGT)*/
#define def_DToT_alnDefs -100   /*score: T and D (AGT)*/
#define def_DToU_alnDefs -100   /*score: T and D (AGT)*/
#define def_DToG_alnDefs -100   /*score: G and D (AGT)*/
#define def_DToC_alnDefs -400   /*score: D and C (AGT)*/

#define def_DToW_alnDefs -100   /*score: D and W (AT)*/
#define def_DToS_alnDefs -300   /*score: D and S (CG)*/
#define def_DToM_alnDefs -300   /*score: D and M (AC)*/
#define def_DToK_alnDefs -100   /*score: D and K (GT)*/
#define def_DToR_alnDefs -100   /*score: D and R (AG)*/
#define def_DToY_alnDefs -300   /*score: D and Y (CT)*/
#define def_DToB_alnDefs -200   /*score: D and B (CGT)*/
#define def_DToD_alnDefs -100   /*score: D and D (AGT)*/
#define def_DToH_alnDefs -200   /*score: D and H (ACT)*/
#define def_DToV_alnDefs -200   /*score: D and V (ACG)*/
#define def_DToN_alnDefs -100   /*score: D and N (ACGT)*/
#define def_DToX_alnDefs -100   /*score: D and X (ACGT)*/


/********************************************************\
* Sec04 Sub14:
*   - h as first base
\********************************************************/

#define def_HToA_alnDefs -100   /*score: H and A (ACT)*/
#define def_HToT_alnDefs -100   /*score: T and H (ACT)*/
#define def_HToU_alnDefs -100   /*score: T and H (ACT)*/
#define def_HToG_alnDefs -400   /*score: G and H (ACT)*/
#define def_HToC_alnDefs -100   /*score: H and C (ACT)*/

#define def_HToW_alnDefs -100   /*score: H and W (AT)*/
#define def_HToS_alnDefs -300   /*score: H and S (CG)*/
#define def_HToM_alnDefs -100   /*score: H and M (AC)*/
#define def_HToK_alnDefs -300   /*score: H and K (GT)*/
#define def_HToR_alnDefs -300   /*score: H and R (AG)*/
#define def_HToY_alnDefs -100   /*score: H and Y (CT)*/
#define def_HToB_alnDefs -200   /*score: H and B (CGT)*/
#define def_HToD_alnDefs -200   /*score: H and D (AGT)*/
#define def_HToH_alnDefs -100   /*score: H and H (ACT)*/
#define def_HToV_alnDefs -200   /*score: H and V (ACG)*/
#define def_HToN_alnDefs -100   /*score: H and N (ACGT)*/
#define def_HToX_alnDefs -100   /*score: H and X (ACGT)*/

/********************************************************\
* Sec04 Sub15:
*   - v as first base
\********************************************************/

#define def_VToA_alnDefs -100   /*score: V and A (ACG)*/
#define def_VToT_alnDefs -400   /*score: T and V (ACG)*/
#define def_VToU_alnDefs -400   /*score: T and V (ACG)*/
#define def_VToG_alnDefs -100   /*score: G and V (ACG)*/
#define def_VToC_alnDefs -100   /*score: V and C (ACG)*/

#define def_VToW_alnDefs -300   /*score: V and W (AT)*/
#define def_VToS_alnDefs -100   /*score: V and S (CG)*/
#define def_VToM_alnDefs -100   /*score: V and M (AC)*/
#define def_VToK_alnDefs -300   /*score: V and K (GT)*/
#define def_VToR_alnDefs -100   /*score: V and R (AG)*/
#define def_VToY_alnDefs -300   /*score: V and Y (CT)*/
#define def_VToB_alnDefs -200   /*score: V and B (CGT)*/
#define def_VToD_alnDefs -200   /*score: V and D (AGT)*/
#define def_VToH_alnDefs -200   /*score: V and H (ACT)*/
#define def_VToV_alnDefs -100   /*score: V and V (ACG)*/
#define def_VToN_alnDefs -100   /*score: V and N (ACGT)*/
#define def_VToX_alnDefs -100   /*score: V and X (ACGT)*/

/********************************************************\
* Sec04 Sub16:
*   - n as first base
\********************************************************/

#define def_NToA_alnDefs -200   /*score: N and A (ACGT)*/
#define def_NToT_alnDefs -200   /*score: T and N (ACGT)*/
#define def_NToU_alnDefs -200   /*score: T and N (ACGT)*/
#define def_NToG_alnDefs -200   /*score: G and N (ACGT)*/
#define def_NToC_alnDefs -200   /*score: N and C (ACGT)*/

#define def_NToW_alnDefs -200   /*score: N and W (AT)*/
#define def_NToS_alnDefs -200   /*score: N and S (CG)*/
#define def_NToM_alnDefs -200   /*score: N and M (AC)*/
#define def_NToK_alnDefs -200   /*score: N and K (GT)*/
#define def_NToR_alnDefs -200   /*score: N and R (AG)*/
#define def_NToY_alnDefs -200   /*score: N and Y (CT)*/
#define def_NToB_alnDefs -200   /*score: N and B (CGT)*/
#define def_NToD_alnDefs -200   /*score: N and D (AGT)*/
#define def_NToH_alnDefs -200   /*score: N and H (ACT)*/
#define def_NToV_alnDefs -200   /*score: N and V (ACG)*/
#define def_NToN_alnDefs -100   /*score: N and N (ACGT)*/
#define def_NToX_alnDefs -100   /*score: N and X (ACGT)*/

/********************************************************\
* Sec04 Sub17:
*   - x as first base
\********************************************************/

#define def_XToA_alnDefs -200   /*score: X and A (ACGT)*/
#define def_XToT_alnDefs -200   /*score: T and X (ACGT)*/
#define def_XToU_alnDefs -200   /*score: T and X (ACGT)*/
#define def_XToG_alnDefs -200   /*score: G and X (ACGT)*/
#define def_XToC_alnDefs -200   /*score: X and C (ACGT)*/

#define def_XToW_alnDefs -100   /*score: X and W (AT)*/
#define def_XToS_alnDefs -100   /*score: X and S (CG)*/
#define def_XToM_alnDefs -100   /*score: X and M (AC)*/
#define def_XToK_alnDefs -100   /*score: X and K (GT)*/
#define def_XToR_alnDefs -100   /*score: X and R (AG)*/
#define def_XToY_alnDefs -100   /*score: X and Y (CT)*/
#define def_XToB_alnDefs -100   /*score: X and B (CGT)*/
#define def_XToD_alnDefs -100   /*score: X and D (AGT)*/
#define def_XToH_alnDefs -100   /*score: X and H (ACT)*/
#define def_XToV_alnDefs -100   /*score: X and V (ACG)*/
#define def_XToN_alnDefs -100   /*score: X and N (ACGT)*/
#define def_XToX_alnDefs -100   /*score: X and X (ACGT)*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec05:
^   - match matrix defaults (nucleotide only)
^   o sec05 sub01:
^     - a as first base
^   o sec05 sub02:
^     - t as first base
^   o sec05 sub03:
^     - u (t) as first base
^   o sec05 sub04:
^     - g as first base
^   o sec05 sub05:
^     - c as first base
^   o sec05 sub06:
^     - w (anonymous) as first base
^   o sec05 sub07:
^     - s (anonymous) as first base
^   o sec05 sub08:
^     - m (anonymous) as first base
^   o sec05 sub09:
^     - k (anonymous) as first base
^   o sec05 sub10:
^     - r (anonymous) as first base
^   o sec05 sub11:
^     - y (anonymous) as first base
^   o sec05 sub12:
^     - b (anonymous) as first base
^   o sec05 sub13:
^     - d (anonymous) as first base
^   o sec05 sub14:
^     - h (anonymous) as first base
^   o sec05 sub15:
^     - v (anonymous) as first base
^   o sec05 sub16:
^     - n (anonymous) as first base
^   o sec05 sub17:
^     - x (anonymous) as first base (technically aa)
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*using 1 for match and 2 for anonymous. So an anonymous
`  match becomes 3
*/

/********************************************************\
* Sec05 Sub01:
*   - a as first base
\********************************************************/

#define def_AEqlA_alnDefs def_ntEql_alnDefs    /*A and A*/
#define def_AEqlT_alnDefs 0                    /*A and T*/
#define def_AEqlU_alnDefs 0                    /*A and T*/
#define def_AEqlG_alnDefs 0                    /*A and G*/
#define def_AEqlC_alnDefs 0                    /*A and C*/

#define def_AEqlW_alnDefs def_anonMatch_alnDefs
#define def_AEqlS_alnDefs def_anonymous_alnDefs
#define def_AEqlM_alnDefs def_anonMatch_alnDefs
#define def_AEqlK_alnDefs def_anonymous_alnDefs
#define def_AEqlR_alnDefs def_anonMatch_alnDefs
#define def_AEqlY_alnDefs def_anonymous_alnDefs
#define def_AEqlB_alnDefs def_anonymous_alnDefs
#define def_AEqlD_alnDefs def_anonMatch_alnDefs
#define def_AEqlH_alnDefs def_anonMatch_alnDefs
#define def_AEqlV_alnDefs def_anonMatch_alnDefs
#define def_AEqlN_alnDefs def_anonMatch_alnDefs
#define def_AEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub02:
*   - t as first base
\********************************************************/

#define def_TEqlA_alnDefs 0                    /*T and A*/
#define def_TEqlT_alnDefs def_ntEql_alnDefs    /*T and T*/
#define def_TEqlU_alnDefs def_ntEql_alnDefs    /*T and U*/
#define def_TEqlG_alnDefs 0                    /*T and G*/
#define def_TEqlC_alnDefs 0                    /*T and C*/

#define def_TEqlW_alnDefs def_anonMatch_alnDefs
#define def_TEqlS_alnDefs def_anonymous_alnDefs
#define def_TEqlM_alnDefs def_anonymous_alnDefs
#define def_TEqlK_alnDefs def_anonMatch_alnDefs
#define def_TEqlR_alnDefs def_anonymous_alnDefs
#define def_TEqlY_alnDefs def_anonMatch_alnDefs
#define def_TEqlB_alnDefs def_anonMatch_alnDefs
#define def_TEqlD_alnDefs def_anonMatch_alnDefs
#define def_TEqlH_alnDefs def_anonMatch_alnDefs
#define def_TEqlV_alnDefs def_anonymous_alnDefs
#define def_TEqlN_alnDefs def_anonMatch_alnDefs
#define def_TEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub03:
*   - u as first base
\********************************************************/

#define def_UEqlA_alnDefs 0                    /*U and A*/
#define def_UEqlT_alnDefs def_ntEql_alnDefs    /*U and T*/
#define def_UEqlU_alnDefs def_ntEql_alnDefs    /*U and U*/
#define def_UEqlG_alnDefs 0                    /*U and G*/
#define def_UEqlC_alnDefs 0                    /*U and C*/

#define def_UEqlW_alnDefs def_anonMatch_alnDefs
#define def_UEqlS_alnDefs def_anonymous_alnDefs
#define def_UEqlM_alnDefs def_anonymous_alnDefs
#define def_UEqlK_alnDefs def_anonMatch_alnDefs
#define def_UEqlR_alnDefs def_anonymous_alnDefs
#define def_UEqlY_alnDefs def_anonMatch_alnDefs
#define def_UEqlB_alnDefs def_anonMatch_alnDefs
#define def_UEqlD_alnDefs def_anonMatch_alnDefs
#define def_UEqlH_alnDefs def_anonMatch_alnDefs
#define def_UEqlV_alnDefs def_anonymous_alnDefs
#define def_UEqlN_alnDefs def_anonMatch_alnDefs
#define def_UEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub04:
*   - g as first base
\********************************************************/

#define def_GEqlA_alnDefs 0                    /*G and A*/
#define def_GEqlT_alnDefs 0                    /*G and T*/
#define def_GEqlU_alnDefs 0                    /*G and T*/
#define def_GEqlG_alnDefs def_ntEql_alnDefs    /*G and G*/
#define def_GEqlC_alnDefs 0                    /*G and C*/

#define def_GEqlW_alnDefs def_anonymous_alnDefs
#define def_GEqlS_alnDefs def_anonMatch_alnDefs
#define def_GEqlM_alnDefs def_anonymous_alnDefs
#define def_GEqlK_alnDefs def_anonMatch_alnDefs
#define def_GEqlR_alnDefs def_anonMatch_alnDefs
#define def_GEqlY_alnDefs def_anonymous_alnDefs
#define def_GEqlB_alnDefs def_anonMatch_alnDefs
#define def_GEqlD_alnDefs def_anonMatch_alnDefs
#define def_GEqlH_alnDefs def_anonymous_alnDefs
#define def_GEqlV_alnDefs def_anonMatch_alnDefs
#define def_GEqlN_alnDefs def_anonMatch_alnDefs
#define def_GEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub05:
*   - c as first base
\********************************************************/

#define def_CEqlA_alnDefs 0                    /*C and A*/
#define def_CEqlT_alnDefs 0                    /*C and T*/
#define def_CEqlU_alnDefs 0                    /*C and T*/
#define def_CEqlG_alnDefs 0                    /*C and G*/
#define def_CEqlC_alnDefs def_ntEql_alnDefs    /*C and C*/

#define def_CEqlW_alnDefs def_anonymous_alnDefs
#define def_CEqlS_alnDefs def_anonMatch_alnDefs
#define def_CEqlM_alnDefs def_anonMatch_alnDefs
#define def_CEqlK_alnDefs def_anonymous_alnDefs
#define def_CEqlR_alnDefs def_anonymous_alnDefs
#define def_CEqlY_alnDefs def_anonMatch_alnDefs
#define def_CEqlB_alnDefs def_anonMatch_alnDefs
#define def_CEqlD_alnDefs def_anonymous_alnDefs
#define def_CEqlH_alnDefs def_anonMatch_alnDefs
#define def_CEqlV_alnDefs def_anonMatch_alnDefs
#define def_CEqlN_alnDefs def_anonMatch_alnDefs
#define def_CEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub06:
*   - w as first base
\********************************************************/

#define def_WEqlA_alnDefs def_anonMatch_alnDefs
#define def_WEqlT_alnDefs def_anonMatch_alnDefs
#define def_WEqlU_alnDefs def_anonMatch_alnDefs
#define def_WEqlG_alnDefs def_anonymous_alnDefs
#define def_WEqlC_alnDefs def_anonymous_alnDefs

#define def_WEqlW_alnDefs def_anonMatch_alnDefs
#define def_WEqlS_alnDefs def_anonymous_alnDefs
#define def_WEqlM_alnDefs def_anonMatch_alnDefs
#define def_WEqlK_alnDefs def_anonMatch_alnDefs
#define def_WEqlR_alnDefs def_anonMatch_alnDefs
#define def_WEqlY_alnDefs def_anonMatch_alnDefs
#define def_WEqlB_alnDefs def_anonMatch_alnDefs
#define def_WEqlD_alnDefs def_anonMatch_alnDefs
#define def_WEqlH_alnDefs def_anonMatch_alnDefs
#define def_WEqlV_alnDefs def_anonMatch_alnDefs
#define def_WEqlN_alnDefs def_anonMatch_alnDefs
#define def_WEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub07:
*   - s as first base
\********************************************************/

#define def_SEqlA_alnDefs def_anonymous_alnDefs
#define def_SEqlT_alnDefs def_anonymous_alnDefs
#define def_SEqlU_alnDefs def_anonymous_alnDefs
#define def_SEqlG_alnDefs def_anonMatch_alnDefs
#define def_SEqlC_alnDefs def_anonMatch_alnDefs

#define def_SEqlW_alnDefs def_anonymous_alnDefs
#define def_SEqlS_alnDefs def_anonMatch_alnDefs
#define def_SEqlM_alnDefs def_anonMatch_alnDefs
#define def_SEqlK_alnDefs def_anonMatch_alnDefs
#define def_SEqlR_alnDefs def_anonMatch_alnDefs
#define def_SEqlY_alnDefs def_anonMatch_alnDefs
#define def_SEqlB_alnDefs def_anonMatch_alnDefs
#define def_SEqlD_alnDefs def_anonMatch_alnDefs
#define def_SEqlH_alnDefs def_anonMatch_alnDefs
#define def_SEqlV_alnDefs def_anonMatch_alnDefs
#define def_SEqlN_alnDefs def_anonMatch_alnDefs
#define def_SEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub08:
*   - m as first base
\********************************************************/

#define def_MEqlA_alnDefs def_anonMatch_alnDefs
#define def_MEqlT_alnDefs def_anonymous_alnDefs
#define def_MEqlU_alnDefs def_anonymous_alnDefs
#define def_MEqlG_alnDefs def_anonymous_alnDefs
#define def_MEqlC_alnDefs def_anonMatch_alnDefs

#define def_MEqlW_alnDefs def_anonMatch_alnDefs
#define def_MEqlS_alnDefs def_anonMatch_alnDefs
#define def_MEqlM_alnDefs def_anonMatch_alnDefs
#define def_MEqlK_alnDefs def_anonymous_alnDefs
#define def_MEqlR_alnDefs def_anonMatch_alnDefs
#define def_MEqlY_alnDefs def_anonMatch_alnDefs
#define def_MEqlB_alnDefs def_anonMatch_alnDefs
#define def_MEqlD_alnDefs def_anonMatch_alnDefs
#define def_MEqlH_alnDefs def_anonMatch_alnDefs
#define def_MEqlV_alnDefs def_anonMatch_alnDefs
#define def_MEqlN_alnDefs def_anonMatch_alnDefs
#define def_MEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub09:
*   - k as first base
\********************************************************/

#define def_KEqlA_alnDefs def_anonymous_alnDefs
#define def_KEqlT_alnDefs def_anonMatch_alnDefs
#define def_KEqlU_alnDefs def_anonMatch_alnDefs
#define def_KEqlG_alnDefs def_anonMatch_alnDefs
#define def_KEqlC_alnDefs def_anonymous_alnDefs

#define def_KEqlW_alnDefs def_anonMatch_alnDefs
#define def_KEqlS_alnDefs def_anonMatch_alnDefs
#define def_KEqlM_alnDefs def_anonymous_alnDefs
#define def_KEqlK_alnDefs def_anonMatch_alnDefs
#define def_KEqlR_alnDefs def_anonMatch_alnDefs
#define def_KEqlY_alnDefs def_anonMatch_alnDefs
#define def_KEqlB_alnDefs def_anonMatch_alnDefs
#define def_KEqlD_alnDefs def_anonMatch_alnDefs
#define def_KEqlH_alnDefs def_anonMatch_alnDefs
#define def_KEqlV_alnDefs def_anonMatch_alnDefs
#define def_KEqlN_alnDefs def_anonMatch_alnDefs
#define def_KEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub10:
*   - r as first base
\********************************************************/

#define def_REqlA_alnDefs def_anonMatch_alnDefs
#define def_REqlT_alnDefs def_anonymous_alnDefs
#define def_REqlU_alnDefs def_anonymous_alnDefs
#define def_REqlG_alnDefs def_anonMatch_alnDefs
#define def_REqlC_alnDefs def_anonymous_alnDefs

#define def_REqlW_alnDefs def_anonMatch_alnDefs
#define def_REqlS_alnDefs def_anonMatch_alnDefs
#define def_REqlM_alnDefs def_anonMatch_alnDefs
#define def_REqlK_alnDefs def_anonMatch_alnDefs
#define def_REqlR_alnDefs def_anonMatch_alnDefs
#define def_REqlY_alnDefs def_anonymous_alnDefs
#define def_REqlB_alnDefs def_anonMatch_alnDefs
#define def_REqlD_alnDefs def_anonMatch_alnDefs
#define def_REqlH_alnDefs def_anonMatch_alnDefs
#define def_REqlV_alnDefs def_anonMatch_alnDefs
#define def_REqlN_alnDefs def_anonMatch_alnDefs
#define def_REqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub11:
*   - y as first base
\********************************************************/

#define def_YEqlA_alnDefs def_anonymous_alnDefs
#define def_YEqlT_alnDefs def_anonMatch_alnDefs
#define def_YEqlU_alnDefs def_anonMatch_alnDefs
#define def_YEqlG_alnDefs def_anonymous_alnDefs
#define def_YEqlC_alnDefs def_anonMatch_alnDefs

#define def_YEqlW_alnDefs def_anonMatch_alnDefs
#define def_YEqlS_alnDefs def_anonMatch_alnDefs
#define def_YEqlM_alnDefs def_anonMatch_alnDefs
#define def_YEqlK_alnDefs def_anonMatch_alnDefs
#define def_YEqlR_alnDefs def_anonymous_alnDefs
#define def_YEqlY_alnDefs def_anonMatch_alnDefs
#define def_YEqlB_alnDefs def_anonMatch_alnDefs
#define def_YEqlD_alnDefs def_anonMatch_alnDefs
#define def_YEqlH_alnDefs def_anonMatch_alnDefs
#define def_YEqlV_alnDefs def_anonMatch_alnDefs
#define def_YEqlN_alnDefs def_anonMatch_alnDefs
#define def_YEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub12:
*   - b as first base
\********************************************************/

#define def_BEqlA_alnDefs def_anonymous_alnDefs
#define def_BEqlT_alnDefs def_anonMatch_alnDefs
#define def_BEqlU_alnDefs def_anonMatch_alnDefs
#define def_BEqlG_alnDefs def_anonMatch_alnDefs
#define def_BEqlC_alnDefs def_anonMatch_alnDefs

#define def_BEqlW_alnDefs def_anonMatch_alnDefs
#define def_BEqlS_alnDefs def_anonMatch_alnDefs
#define def_BEqlM_alnDefs def_anonMatch_alnDefs
#define def_BEqlK_alnDefs def_anonMatch_alnDefs
#define def_BEqlR_alnDefs def_anonMatch_alnDefs
#define def_BEqlY_alnDefs def_anonMatch_alnDefs
#define def_BEqlB_alnDefs def_anonMatch_alnDefs
#define def_BEqlD_alnDefs def_anonMatch_alnDefs
#define def_BEqlH_alnDefs def_anonMatch_alnDefs
#define def_BEqlV_alnDefs def_anonMatch_alnDefs
#define def_BEqlN_alnDefs def_anonMatch_alnDefs
#define def_BEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub13:
*   - d as first base
\********************************************************/

#define def_DEqlA_alnDefs def_anonMatch_alnDefs
#define def_DEqlT_alnDefs def_anonMatch_alnDefs
#define def_DEqlU_alnDefs def_anonMatch_alnDefs
#define def_DEqlG_alnDefs def_anonMatch_alnDefs
#define def_DEqlC_alnDefs def_anonymous_alnDefs

#define def_DEqlW_alnDefs def_anonMatch_alnDefs
#define def_DEqlS_alnDefs def_anonMatch_alnDefs
#define def_DEqlM_alnDefs def_anonMatch_alnDefs
#define def_DEqlK_alnDefs def_anonMatch_alnDefs
#define def_DEqlR_alnDefs def_anonMatch_alnDefs
#define def_DEqlY_alnDefs def_anonMatch_alnDefs
#define def_DEqlB_alnDefs def_anonMatch_alnDefs
#define def_DEqlD_alnDefs def_anonMatch_alnDefs
#define def_DEqlH_alnDefs def_anonMatch_alnDefs
#define def_DEqlV_alnDefs def_anonMatch_alnDefs
#define def_DEqlN_alnDefs def_anonMatch_alnDefs
#define def_DEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub14:
*   - h as first base
\********************************************************/

#define def_HEqlA_alnDefs def_anonMatch_alnDefs
#define def_HEqlT_alnDefs def_anonMatch_alnDefs
#define def_HEqlU_alnDefs def_anonMatch_alnDefs
#define def_HEqlG_alnDefs def_anonymous_alnDefs
#define def_HEqlC_alnDefs def_anonMatch_alnDefs

#define def_HEqlW_alnDefs def_anonMatch_alnDefs
#define def_HEqlS_alnDefs def_anonMatch_alnDefs
#define def_HEqlM_alnDefs def_anonMatch_alnDefs
#define def_HEqlK_alnDefs def_anonMatch_alnDefs
#define def_HEqlR_alnDefs def_anonMatch_alnDefs
#define def_HEqlY_alnDefs def_anonMatch_alnDefs
#define def_HEqlB_alnDefs def_anonMatch_alnDefs
#define def_HEqlD_alnDefs def_anonMatch_alnDefs
#define def_HEqlH_alnDefs def_anonMatch_alnDefs
#define def_HEqlV_alnDefs def_anonMatch_alnDefs
#define def_HEqlN_alnDefs def_anonMatch_alnDefs
#define def_HEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub15:
*   - v as first base
\********************************************************/

#define def_VEqlA_alnDefs def_anonMatch_alnDefs
#define def_VEqlT_alnDefs def_anonymous_alnDefs
#define def_VEqlU_alnDefs def_anonymous_alnDefs
#define def_VEqlG_alnDefs def_anonMatch_alnDefs
#define def_VEqlC_alnDefs def_anonMatch_alnDefs

#define def_VEqlW_alnDefs def_anonMatch_alnDefs
#define def_VEqlS_alnDefs def_anonMatch_alnDefs
#define def_VEqlM_alnDefs def_anonMatch_alnDefs
#define def_VEqlK_alnDefs def_anonMatch_alnDefs
#define def_VEqlR_alnDefs def_anonMatch_alnDefs
#define def_VEqlY_alnDefs def_anonMatch_alnDefs
#define def_VEqlB_alnDefs def_anonMatch_alnDefs
#define def_VEqlD_alnDefs def_anonMatch_alnDefs
#define def_VEqlH_alnDefs def_anonMatch_alnDefs
#define def_VEqlV_alnDefs def_anonMatch_alnDefs
#define def_VEqlN_alnDefs def_anonMatch_alnDefs
#define def_VEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub16:
*   - n as first base
\********************************************************/

#define def_NEqlA_alnDefs def_anonMatch_alnDefs
#define def_NEqlT_alnDefs def_anonMatch_alnDefs
#define def_NEqlU_alnDefs def_anonMatch_alnDefs
#define def_NEqlG_alnDefs def_anonMatch_alnDefs
#define def_NEqlC_alnDefs def_anonMatch_alnDefs

#define def_NEqlW_alnDefs def_anonMatch_alnDefs
#define def_NEqlS_alnDefs def_anonMatch_alnDefs
#define def_NEqlM_alnDefs def_anonMatch_alnDefs
#define def_NEqlK_alnDefs def_anonMatch_alnDefs
#define def_NEqlR_alnDefs def_anonMatch_alnDefs
#define def_NEqlY_alnDefs def_anonMatch_alnDefs
#define def_NEqlB_alnDefs def_anonMatch_alnDefs
#define def_NEqlD_alnDefs def_anonMatch_alnDefs
#define def_NEqlH_alnDefs def_anonMatch_alnDefs
#define def_NEqlV_alnDefs def_anonMatch_alnDefs
#define def_NEqlN_alnDefs def_anonMatch_alnDefs
#define def_NEqlX_alnDefs def_anonMatch_alnDefs

/********************************************************\
* Sec05 Sub17:
*   - x as first base
\********************************************************/

#define def_XEqlA_alnDefs def_anonMatch_alnDefs
#define def_XEqlT_alnDefs def_anonMatch_alnDefs
#define def_XEqlU_alnDefs def_anonMatch_alnDefs
#define def_XEqlG_alnDefs def_anonMatch_alnDefs
#define def_XEqlC_alnDefs def_anonMatch_alnDefs

#define def_XEqlW_alnDefs def_anonMatch_alnDefs
#define def_XEqlS_alnDefs def_anonMatch_alnDefs
#define def_XEqlM_alnDefs def_anonMatch_alnDefs
#define def_XEqlK_alnDefs def_anonMatch_alnDefs
#define def_XEqlR_alnDefs def_anonMatch_alnDefs
#define def_XEqlY_alnDefs def_anonMatch_alnDefs
#define def_XEqlB_alnDefs def_anonMatch_alnDefs
#define def_XEqlD_alnDefs def_anonMatch_alnDefs
#define def_XEqlH_alnDefs def_anonMatch_alnDefs
#define def_XEqlV_alnDefs def_anonMatch_alnDefs
#define def_XEqlN_alnDefs def_anonMatch_alnDefs
#define def_XEqlX_alnDefs def_anonMatch_alnDefs

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
: Copyright (c) 3034 jeremyButtler
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
