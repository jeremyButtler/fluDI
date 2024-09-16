/*########################################################
# Name indexToCoord
#  o Holds functions to convert alignment matrix index's
#    to coordinates
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of Functions
'   o fun01: qryCoord_water
'     - gets the query coordinates of the query sequence
'       in an matrix.
'   o fun02: refCoord_water
'     - gets the coordinates of the reference sequence in
'       in an matrix.
'   o fun03: indexToCoord_water
'     - gets the coordinates of the reference and query
'       sequence in an matrix.
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "indexToCoord.h"

/*-------------------------------------------------------\
| Fun01: qryCoord_indexToCoord
|   - Gets the query coordinates of the query sequence in
|     an matrix.
| Input
|   - refLenUL:
|     o length of the reference sequence
|   - indexUL:
|     o indexUL to convert to cooridnates
| Output:
|  - Returns:
|    o the query coordinate that was in the indexUL
\-------------------------------------------------------*/
unsigned long
qryCoord_indexToCoord(
   unsigned long refLenUL,
   unsigned long indexUL
){
     indexUL = ( indexUL / (refLenUL + 1) );
     indexUL -= (!! indexUL);
     return indexUL;
} /*refCoord_indexToCoord*/
   /*Logic:
   `  - indexUL = indexUL / (refLenUL + 1):
   `    o gives the row i ended on
   `    o The + 1 accounts for the gap row
   `  - indexUL - (!! indexUL)
   `    o indexUL - 0 if on gap row (indexUL = 0)
   `    o indexUL - 1 if not gap row (indexUL > 0)
   `    o converts non-gap row values to index 0
   */

/*-------------------------------------------------------\
| Fun02: refCoord_indexToCoord
|   - gets the coordinates of the reference sequence in
|     an matrix.
| Input
|   - refLen:
|     o length of the reference sequence
|   - index:
|     o index to convert to cooridnates
|   - refCoord:
|     o will hold the coordinate of the reference sequence
| Output:
|  - Returns
|    o reference coordinate that was in the index
\-------------------------------------------------------*/
unsigned long
refCoord_indexToCoord( 
   unsigned long refLenUL,
   unsigned long indexUL
){
     indexUL = ( indexUL % (refLenUL + 1) );
     indexUL -= (!! indexUL);
     return indexUL;
} /*refCoord_indexToCoord*/
   /*Logic:
   `  - indexUL = indexUL % (refLenUL + 1):
   `    o gives the column i ended on
   `    o The + 1 accounts for the gap column
   `  - indexUL - (!! indexUL)
   `    o indexUL - 0 if on gap colunm (indexUL = 0)
   `    o indexUL - 1 if not gap colunm (indexUL > 0)
   `    o converts non-gap column values to index 0
   */

/*-------------------------------------------------------\
| Fun03: indexToCoord
|   - gets the coordinates of the reference and query
|     sequence in an matrix.
| Input
|   - refLenUL:
|     o length of the reference sequence
|   - indexUL:
|     o index to convert to cooridnates
|   - refCoordULPtr:
|     o pointer to unsigned long to hold reference
|       coordinate
|   - qryCoordULPtr:
|     o pointer to unsigned long to hold query coordinate
| Output:
|  - Modifies:
|    o refCoordULPtr to reference coordinate
|    o qryCoordULPtr to query coordinate
\-------------------------------------------------------*/
void
indexToCoord(
   unsigned long refLenUL,
   unsigned long indexUL,
   unsigned long *refCoordULPtr,
   unsigned long *qryCoordULPtr
){
   *refCoordULPtr =
      refCoord_indexToCoord(
         refLenUL,
         indexUL
      ); /*get reference coordinate*/

   *qryCoordULPtr =
      qryCoord_indexToCoord(
         refLenUL,
         indexUL
      ); /*get reference coordinate*/
} /*indexToCoord*/
