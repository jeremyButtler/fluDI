
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' varFilt SOF:
'   - has functions for filtering reads with variants
'   o header:
'     - defined variables, forward declarations, guards
'   o fun02: getProf_varFilt
'     - builds a profile for filtering
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef VARIANT_FITER_H
#define VARIANT_FITER_H

typedef struct samEntry samEntry;
typedef struct set_tbCon set_tbCon;
typedef struct conNt_tbCon conNt_tbCon;
typedef struct index_clustST index_clustST;

#define def_memErr_varFilt 1
#define def_fileErr_varFilt 2

/*for filtering variants*/
#define def_rmA_varFilt 1
#define def_rmC_varFilt 2
#define def_rmG_varFilt 4
#define def_rmT_varFilt 8

/*-------------------------------------------------------\
| Fun02: getProf_varFilt
|   - builds a profile for filtering
| Input:
|   - minSupF:
|     o minimum percent support to keep variant
|   - varUIPtr:
|     o pointer to unsigned int to hold number of variants
|       kept
|   - conAryST:
|     o conNt_tbCon array with the uncollapsed consensus
|       to use in variant filtering
|   - lenConUI:
|     o length (number of bases) in conAryST
| Output:
|   - Modifies:
|     o varUIPtr to have number of found variants
|   - Returns:
|     o pointer to signed char array with variants to
|       filter out
\-------------------------------------------------------*/
signed char *
getProf_varFilt(
   float minSupF,                /*min % sup to filter*/
   unsigned int *varUIPtr,       /*number variants found*/
   struct conNt_tbCon *conAryST, /*uncollapsed con*/
   unsigned int lenConUI         /*consensus length*/
);

#endif
