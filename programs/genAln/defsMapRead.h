/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' defsMapRead SOF: Start Of File
'   - default values for mapRead
'   o header:
'     - gaurds
'   o sec01:
'     - kmer scan variables
'   o sec02:
'     - chaining and scoring variables
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - gaurds
\-------------------------------------------------------*/

#ifndef DEFAULTS_MAP_READ_H
#define DEFAULTS_MAP_READ_H

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec01:
^   - kmer scan variables
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*kmer scan defaults*/
#define def_minLen_defsMapRead 0
#define def_secLen_defsMapRead 100
#define def_thirdLen_defsMapRead 800
#define def_fourthLen_defsMapRead 3000
#define def_fifthLen_defsMapRead 12000
#define def_sixthLen_defsMapRead 50000
#define def_finalLen_defsMapRead 200000

#define def_minKmer_defsMapRead 5 
#define def_secKmer_defsMapRead 6 
#define def_thirdKmer_defsMapRead 7
#define def_fourthKmer_defsMapRead 8
#define def_fifthKmer_defsMapRead 9
#define def_sixthKmer_defsMapRead 10
#define def_finalKmer_defsMapRead 11

#define def_minPerc_defsMapRead 0.7f
#define def_secPerc_defsMapRead 0.7f
#define def_thirdPerc_defsMapRead 0.8f
#define def_fourthPerc_defsMapRead 0.8f
#define def_fifthPerc_defsMapRead 0.8f
#define def_sixthPerc_defsMapRead 0.9f
#define def_finalPerc_defsMapRead 0.9f

#define def_minChain_defsMapRead 2
#define def_secChain_defsMapRead 15
#define def_thirdChain_defsMapRead 15
#define def_fourthChain_defsMapRead 20
#define def_fifthChain_defsMapRead 20
#define def_sixthChain_defsMapRead 25
#define def_finalChain_defsMapRead 25


/* What these defaults mean:
` - as a rule of thumb I am trying to keep the kmer
`   space double the maximum widnow size, that way
`   there is plenty of empty room
` - for the minimum length I want more sparcity, since
`   there are even less kmers
` - after 11mers, the kmer space gets massive and I need
`   16mb (unless used hash method), so best to limit,
`   - > 1mb reads are rare, so reasnable limit
` - the percent kmer values are just a random guess
` - min: < 100
`   - 5mers; kmer space is 1024
`   - max window 250 - 3
`   - need at least 70% matching kmers to do aligment
`   - chains are at least 2 kmers (6 bases) long
` - two: 100 to 799:
`   - 6mers, which have kmer space of 4098
`   - maximum window length 2000 - 3
`   - need at least 70% matching kmers to do aligment
`   - chains are at least 15 kmers (20 bases) long
` - third: 880 to 2999:
`   - 7mers; kmer space is ~ 16000
`   - max window 7500 - 3
`   - need at least 80% matching kmers to do aligment
`   - chains are at least 15 kmers (21 bases) long
` - fourth: 3000 to 11999:
`   - 8mers; kmer space is ~ 64000
`   - max window lenth 30000 - 3
`   - need at least 80% matching kmers to do aligment
`   - chains are at least 20 kmers (27 bases) long
` - fifth: 12000 to 49999:
`   - 9mers; kmer space is ~ 256000
`   - max window length 175000 - 3
`   - need at least 80% matching kmers to do aligment
`   - chains are at least 20 kmers (28 bases) long
` - sixth: 50000 to 199999
`   - 10mers kmer space is ~ 1000000
`   - max window length 500000 - 3
`   - need at least 90% matching kmers to do aligment
`   - chains are at least 25 kmers (34 bases) long
` - last: > 199999
`   - 11mers kmer space is ~ 4000000
`   - no max window length
`   - need at least 90% matching kmers to do aligment
`   - chains are at least 25 kmers (35 bases) long
*/

#define def_percExtraNt_defsMapRead 1.5f
   /*this is (query length) + 150% * (query length)=250%*/
#define def_percShift_defsMapRead 1.25f
   /*shift window by 1.25 query length
   `   works well because have 2.5 query lengths for my
   `   window
   */

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec02:
^   - chaining and scoring variables
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define def_minScore_defsMapRead 0.9f
    /*min % score to keep alignment*/
#define def_minPercLen_defsMapRead 0.75f
    /*min % of query length mapped*/

#define def_subAln_defsMapRead 1 /*keep best sub-aln*/

/*the Needleman/Waterman matrix's can be more complex,
`  however, in this case keep same as core.
`  Also, their scores are only used in gap filling,
`  not in the final score
*/
#endif
