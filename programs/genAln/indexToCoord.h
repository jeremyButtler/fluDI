/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' indexToCoord SOF: Start Of Functions
'  - holds functions to convert alignment matrix index's
'    to coordinates
'   o header:
'     - guards
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
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards
\-------------------------------------------------------*/

#ifndef INDEX_TO_COORDINATES_H
#define INDEX_TO_COORDINATES_H

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
|    o query coordinate that was in the indexUL
\-------------------------------------------------------*/
unsigned long
qryCoord_indexToCoord(
   unsigned long refLenUL,
   unsigned long indexUL
);

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
);

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
);

#endif
