/*########################################################
# Name: fluSeg
#   - holds flu segment information (numbers, length,
#     ids/names, and id sequences)
########################################################*/

#ifndef FLU_SEGMENT_H
#define FLU_SEGMENT_H

#define def_PB2Num_fluSeg 0
#define def_PB1Num_fluSeg 1
#define def_PANum_fluSeg 2
#define def_HANum_fluSeg 3
#define def_NPNum_fluSeg 4
#define def_NANum_fluSeg 5
#define def_MNum_fluSeg 6
#define def_NSNum_fluSeg 7 /*always last segment*/

static signed char
   segIdAryStr_fluSeg[def_NSNum_fluSeg + 1][4] =
   {
      "PB2",
      "PB1",
      "PA", 
      "HA", 
      "NP", 
      "NA", 
      "M",  
      "NS"  
   }; /*forIdAryStr_fluSeg*/

static signed short
   segLenArySS_fluSeg[def_NSNum_fluSeg + 1] =
   {
      2341, /*PB2*/
      2341, /*PB1*/
      2233, /*PA*/
      1778, /*HA*/
      1565, /*NP*/
      1413, /*NA*/
      1027, /*M*/
       890  /*NS*/
   }; /*segLenArySS_fluSeg*/

static signed char
   forSeqAryStr_fluSeg[def_NSNum_fluSeg + 1][10] =
   {
      "TC",   /*PB2*/
      "CA",   /*PB1*/
      "TAC",  /*PA*/
      "GG",   /*HA*/
      "GTA",  /*NP*/
      "AGT",  /*NA*/
      "TAG",  /*M*/
      "GTG"   /*NS*/
   }; /*forSeqAryStr_fluSeg*/

static signed char
   revSeqAryStr_fluSeg[def_NSNum_fluSeg + 1][10] =
   {
      "TCGT",  /*PB2*/
      "CATT",  /*PB1*/
      "TACT",  /*PA*/
      "GTGT",  /*HA*/
      "GTAT",  /*NP*/
      "AGTT",  /*NA*/
      "TAGT",  /*M*/
      "GTGT"   /*NS*/
   }; /*revSeqsAryStr_fluSeg*/

#endif
