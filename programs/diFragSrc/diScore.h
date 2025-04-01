/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' diScore SOF: Start Of File:
'   - sets scoring matrix up for diFrag
'   o header:
'     - guards and defined variables
'   o fun01: set_diScore
'     - sets up scoring matrix for diFrag
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards and defined variables
\-------------------------------------------------------*/

#ifndef DEFICTIVE_INTERFERING_SCORING_H
#define DEFICTIVE_INTERFERING_SCORING_H

struct alnSet;

#define def_match_diScore 500  /*really 5*/
/*#define def_snp_diScore -400*/   /*really -4*/
#define def_snp_diScore -40   /*really -0.04*/
#define def_open_diScore -1000 /*really -10*/
#define def_extend_diScore -10 /*really 0.1*/
#define def_adjust_diScore 100 /*adjust to reall score*/

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
);

#endif
