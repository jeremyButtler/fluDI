/*########################################################
# Name: seqST
#   - Holds functions for reading in or manipulating
#     sequences
########################################################*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'  o header:
'    - included libraries
'  o .h st01: seqST
'    - holds an single sequence (fasta/fastq)
'  o .c fun01 addLine_seqST:
'    - add characters from file to buffer, if needed 
'      resize. This will only read in till the end of the
'      line
'  o fun02 getFqSeq_seqST:
'    - reads a fastq sequence from a fastq file
'  o fun03 getFaSeq_seqST:
'    - grabs the next read in the fasta file
'  o fun04 revComp_seqST:
'    - reverse complement a sequence
'  o fun05 blank_seqST:
'    - sets values in seqST to zero
'  O fun06 init_seqST:
'    - sets values in seqST to blank values
'  o fun07 freeStack_seqST:
'    - frees variables in an seqST (calls blank_seqST)
'  o fun08 freeHeap_seqST:
'    - frees an seqST strucuter (calls fredSeqSTStack)
'  o fun09: freeHeapAry_seqST
'    - frees an array of seqST's
'  o fun10 cpIdEndPad_seqST:
'     - copies read id to a buffer and adds in endIdC to
'       the end. If needed, this function will add right
'       padding of spaces to the end.
'  o fun11: cp_seqST
'    - copies an seqST structure
'  o license:
'    - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "seqST.h"

#include <stdio.h>

#include "ulCp.h"

/*no .c files*/
#include "dataTypeShortHand.h"

/*-------------------------------------------------------\
| Fun01: addLine_seqST
|  - read line of characters into the buffer.If needed
|    this will resize the buffer.
| Input:
|  - buffStr:
|    o pointer to c-string to add buffer into.
|  - lenBuffUL:
|    o pointer to the size of buffStr
|  - curBuffUL:
|    o pointer to the number of elements in buffStr  
|  - resBuffUL:
|    o how much to expand buffer by when the buffer is
|      full
|  - inFILE:
|    o pointer to FILE handle to get the next lein from
| Output:
|   - Modifies:
|     o buffStr to hold the next line.
|       - buffStr is resizied if it is to small to hold
|         the next line.
|       - buffStr + lenBuffUL - 2 will be '\0' or '\n'
|       - buffStr set to 0 for memory allocation errors
|     o curBuffUL: Has the number of chars in the buffer
|     o lenBuffUL: Has the buffer size
|     o inFILE: Points to next line in file
|   - Returns:
|     o 0 if no errors
|     o def_fileErr_seqST if was end of file (EOF)
|     o def_memErr_seqST if had a memory allocation error
\-------------------------------------------------------*/
unsigned char
addLine_seqST(
    signed char **buffStr,  /*Buffer to add data to*/
    unsigned long *lenBuffUL,/*Size of the buffer*/
    unsigned long *curBuffUL,/*Number of chars in buffer*/
    unsigned long resBuffUL, /*How much to resize buff by*/
    unsigned long *filePosUL, /*position at in file*/
    void *inFILE             /*File to grab data from*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC: addLine_seqST
   '  - read line of characters into the buffer.If needed
   '    this will resize the buffer.
   '   o fun01 sec01:
   '     - variable declerations
   '   o fun01 sec02:
   '     - check if need to resize the buffer
   '   o fun01 sec03:
   '     - read in the next line in the buffer
   '   o fun01 sec04:
   '     - if at end of file, update read in lengths
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec01:
    ^  - variable declerations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    schar *tmpStr = 0;
    unsigned long spareBuffUL = 0;
    unsigned long tmpUL = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec02:
    ^  - check if need to resize the buffer
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(*curBuffUL == 0 && *lenBuffUL > 0)
        spareBuffUL = *lenBuffUL;

    else if(*lenBuffUL == 0 || *curBuffUL-1 >= *lenBuffUL)
    { /*If need to resize the buffer (-1 for 0 index)*/
        *lenBuffUL += resBuffUL;

        if(*lenBuffUL == 0)
          tmpStr =malloc(sizeof(schar) * (*lenBuffUL+9));
        else
        { /*Else I need to resize the buffer*/
            tmpStr =
              realloc(
                *buffStr,
                sizeof(schar) * (*lenBuffUL + 9)
            ); /*Resize hte buffer*/
        } /*Else I need to resize the buffer*/
       
        if(tmpStr == 0)
            return def_memErr_seqST; /*Memory error*/

        *buffStr = tmpStr;

        /*Amount of extra space in the buffer*/
        spareBuffUL = resBuffUL;
    } /*If need to resize the buffer*/

    else
        spareBuffUL = *lenBuffUL - *curBuffUL;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec03:
    ^  - read in the next line in the buffer
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    /*Set up marker to catch when read in the entire line*/
    tmpStr = *buffStr + *curBuffUL;
 
    while(
       fgets((char *) tmpStr,spareBuffUL,(FILE *) inFILE)
    ){ /*While I have lines to read*/
        tmpUL = endLine_ulCp(tmpStr);
        *curBuffUL += tmpUL;
        tmpStr = *buffStr + *curBuffUL;
        *filePosUL += tmpUL;
        spareBuffUL -= tmpUL;

        if(spareBuffUL < 16)
        { /*If: need to resize buffer*/
            *lenBuffUL += resBuffUL;

            tmpStr =
                realloc(
                    *buffStr,
                    sizeof(char) * (*lenBuffUL + 9)
                ); /*Resize the buffer*/

            if (tmpStr == 0)
                return def_memErr_seqST; /*Memory error*/

            *buffStr = tmpStr;

            spareBuffUL = resBuffUL;
            tmpStr = *buffStr + *curBuffUL;
        } /*If: need to resize buffer*/

        if (*tmpStr == '\n')
        { /*If: found end of line*/
            ++(*filePosUL);
            return 0; /*finshed with line*/
        } /*If: found end of line*/
    } /*While I have lines to read*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun01 Sec04:
    ^  - if at end of file, update read in lengths
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    
    return def_EOF_seqST; /*End of file*/
} /*addLine_seqST*/

/*-------------------------------------------------------\
| Fun02: getFqSeq_seqST
|  - grabs the next read in the fastq file
| Input:
|  - fqFILE:
|    o pointer to FILE handle to a fastq file to get the
|      next read from
|  - seqSTPtr:
|    o pointer to seqSTPtr structure to hold next read
| Output:
|  - Modifies:
|    o refStruct to hold the read in fastq entry & sets
|      its pointers
|  - Returns:
|     o 0: for no errors
|     o def_EOF_seqST: if EOF
|     o def_fileErr_seqST If file was not a fastq file
|     o def_badLine_seqST | def_fileErr_seqST for  a
|       invalid fastq entry
|     o def_memErr_seqST If malloc failed to find memory
\-------------------------------------------------------*/
signed char
getFqSeq_seqST(
  void *fqFILE,           /*fastq file with sequence*/
  struct seqST *seqSTPtr  /*will hold one fastq entry*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC: readRefFqSeq
   '  -  grabs the next read in the fastq file
   '  o fun02 sec01:
   '    - variable declarations
   '  o fun02 sec02:
   '    - check if need to allocate memory for buffer
   '  o fun02 sec03:
   '    - read in the first data
   '  o fun02 sec04:
   '    - if not at file end, see if have the full entry
   '  o fun02 sec05:
   '    - read till end of file/check if valid fastq entry
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun02 Sec01:
    ^  - variable declarations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    schar tmpSC = 'C';         /*initilize sequence loop*/

    unsigned short extraBuffUS = 1024;
    unsigned long tmpBuffUL = 0;
    unsigned long filePosUL = 0;

    /*Holds number of lines in sequence entry*/
    unsigned char numLinesUC = 0;

    unsigned char errUC = 0;
    schar *oldIterStr = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun02 Sec02:
    ^  - read in the header
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if((FILE *) fqFILE == 0)
        return def_fileErr_seqST;    /*no file*/

    seqSTPtr->lenIdUL = 0;

    errUC =
        addLine_seqST(
            &seqSTPtr->idStr,   /*C-string; hold header*/
            &seqSTPtr->lenIdBuffUL,/*length header buff*/
            &seqSTPtr->lenIdUL, /*number bytes in buffer*/
            extraBuffUS,        /*amount to increase buf*/
            &filePosUL,
            (FILE *) fqFILE     /*fq file with header*/
    ); /*get the header (will resize as needed)*/


    if(errUC)
       return errUC; /*EOF or memory allocation error*/

    /*Account for the new line*/
    --seqSTPtr->lenIdUL;
    seqSTPtr->idStr[seqSTPtr->lenIdUL] = '\0';
    
    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun02 Sec03:
    ^  - read in the sequence & spacer
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    seqSTPtr->lenSeqUL = 0;
    /*need to set this up so the loop does not error out*/
    oldIterStr = &tmpSC;

    while(*oldIterStr != '+')
    { /*While I have not reached the spacer entry*/
        /*So can get to this position later*/
        tmpBuffUL = seqSTPtr->lenSeqUL;

        errUC =
            addLine_seqST(
                &seqSTPtr->seqStr,    /*buff to copy seq*/
                &seqSTPtr->lenSeqBuffUL,/*size; buffer*/
                &seqSTPtr->lenSeqUL,/*length of sequence*/
                extraBuffUS,     /*resize buff by*/
                &filePosUL,
                (FILE *) fqFILE  /*fq file with sequence*/
        ); /*Get the header*/

        /*set up new lines to be removed on next read*/
        oldIterStr =
           seqSTPtr->seqStr + seqSTPtr->lenSeqUL - 1;

        while(*oldIterStr < 33)
        { /*While removing new lines*/
            --seqSTPtr->lenSeqUL;
            --oldIterStr;
        } /*While removing new lines*/

        /*Check for memory errors (64) & invalid entries*/
        if(errUC & def_memErr_seqST)
           return errUC; 

        if(errUC & def_EOF_seqST)
           return
              def_badLine_seqST | def_fileErr_seqST;

        /*Get on first character in the new buffer*/
        oldIterStr = seqSTPtr->seqStr + tmpBuffUL;
        ++numLinesUC; /*number of '\n' in sequence entry*/
    } /*While I have not reached the spacer entry*/

    --numLinesUC; /*Account for the overcounting*/

    oldIterStr = seqSTPtr->seqStr + seqSTPtr->lenSeqUL;

    while(*oldIterStr != '+')
    { /*While not at start of spacer entry*/
        --oldIterStr;
        --seqSTPtr->lenSeqUL; /*Acount for white space*/
    } /*While not at start of spacer entry*/

    *oldIterStr = '\0';  /*Set spacer to null*/
    --oldIterStr;
    --seqSTPtr->lenSeqUL;    /*Acount for white space*/

    while(*oldIterStr < 33)
    { /*While have white space at end to remove*/
        *oldIterStr = '\0';
        --oldIterStr;
        --seqSTPtr->lenSeqUL; /*acount for white space*/
    } /*While have white space at end to remove*/
    
    ++seqSTPtr->lenSeqUL; /*account being one base off*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun02 Sec04:
    ^  - read in the q-score entry
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    seqSTPtr->lenQUL = 0;

    while(numLinesUC > 0)
    { /*Loop: I need to read in the Q-score entry*/
        errUC =
            addLine_seqST(
              &seqSTPtr->qStr,  /*buffer; q-score entry*/
              &seqSTPtr->lenQBuffUL,/*size of buffer*/
              &seqSTPtr->lenQUL, /*length; Q-score entry*/
              extraBuffUS,      /*for reszing buffer*/
              &filePosUL,
              (FILE *) fqFILE   /*fq file with q-score*/
        ); /*get the header*/

        /*set up white space to removed on next read*/
        oldIterStr =
           seqSTPtr->qStr + seqSTPtr->lenQUL - 1;

        while(*oldIterStr < 33) /*removes white space*/
        { /*While removing new lines*/
            --seqSTPtr->lenQUL;
            --oldIterStr;
        } /*While removing new lines*/

        /*Check for errors; memory=64/invalid file=130*/
        if(errUC & def_memErr_seqST)
           return errUC; 

        if(errUC & def_EOF_seqST)
           return
              def_badLine_seqST | def_fileErr_seqST;

        --numLinesUC; /*Number of lines to read in*/
    } /*Loop: I need to read in the Q-score entry*/

    /*Remove any white space at the end*/
    oldIterStr = seqSTPtr->qStr + seqSTPtr->lenQUL;

    while(*oldIterStr < 33)
    { /*While have white space at end to remove*/
        *oldIterStr = '\0';
        --oldIterStr;
        --seqSTPtr->lenQUL; /*Acount for white space*/
    } /*While have white space at end to remove*/

    seqSTPtr->endAlnUL = seqSTPtr->lenSeqUL - 1;
    seqSTPtr->offsetUL = 0;

    return errUC; /*0 = more reads or def_EOF_seqST*/
} /*getFqSeq_seqST*/

/*-------------------------------------------------------\
| Fun03: getFaSeq_seqST
|  -  grabs the next read in the fasta file
| Input:
|  - faFILE:
|    o pointer to FILE handle to a fasta file to get the
|      next read from
|  - seqSTPtr:
|    o pointer to seqSTPtr structure to hold next read
| Output:
|  - Modifies:
|    o seqSTPtr to hold one fasta entry
|  - Returns:
|     o 0: for no errors
|     o def_EOF_seqST: if EOF
|     o def_fileErr_seqST If file was not a fastq file
|     o def_badLine_seqST | def_fileErr_seqST for  a
|       invalid fasta entry
|     o def_memErr_seqST If malloc failed to find memory
| Note:
|   - This will remove new lines from the sequence.
|   - This will only remove spaces or tabs at the end of
|     each sequence entry, so "atg atc \n" will got to
|     "atc atc".
\-------------------------------------------------------*/
signed char
getFaSeq_seqST(
  void *faFILE,           /*fasta file with sequence*/
  struct seqST *seqSTPtr  /*will hold one fastq entry*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC: getFaSeq_seqST
   '  - grabs the next read in the fasta file
   '  o fun03 sec01:
   '    - variable declarations
   '  o fun03 sec02:
   '    - check if need to allocate memory for buffer
   '  o fun03 sec03:
   '    - read in the sequence
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun03 Sec01:
    ^  - variable declarations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /*Used To initilize the sequence loop*/
    schar tmpSC = 'C';
    unsigned long tmpUL = 0;
    unsigned long filePosUL = 0; /*reset file lengths*/

    unsigned short extraBuffUS = 1024;
    signed char errSC = 0;
    schar *tmpStr= 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun03 Sec02:
    ^  - read in header
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    seqSTPtr->lenIdUL = 0;

    errSC =
        (signed char)
        addLine_seqST(
            &seqSTPtr->idStr,   /*buffer to hold header*/
            &seqSTPtr->lenIdBuffUL,/*length of buffer*/
            &seqSTPtr->lenIdUL, /*length of header*/
            extraBuffUS,        /*resize buff by*/
            &filePosUL,
            (FILE *) faFILE     /*Fasta file with header*/
    ); /*Get the header (will resize as needed)*/

    if(errSC)
    { /*If: had error*/
       if(errSC == def_EOF_seqST)
          goto eof_fun03_sec04_sub02;
       else
          goto memErr_fun03_sec04_sub04;
    } /*If: had error*/

    if(seqSTPtr->idStr[0] == '>')
    { /*If: have header*/
       cpLen_ulCp(
          seqSTPtr->idStr,
          &seqSTPtr->idStr[1],
          seqSTPtr->lenIdUL - 1
       );

       --seqSTPtr->lenIdUL;
    } /*If: have header*/

    while(seqSTPtr->idStr[seqSTPtr->lenIdUL - 1] < 33)
    { /*Loop: remove white space at end*/
       if(seqSTPtr->idStr[seqSTPtr->lenIdUL - 1] == '\0')
          break;

       --seqSTPtr->lenIdUL;
    } /*Loop: remove white space at end*/

    seqSTPtr->idStr[seqSTPtr->lenIdUL] = '\0';

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun03 Sec03:
    ^   - read in the sequence
    ^   o fun03 sec03 sub01:
    ^     - get input + start loop
    ^   o fun03 sec03 sub02:
    ^     - remove white space & check if end of entry
    ^   o fun03 sec03 sub03:
    ^     - resize (if needed) and add next character
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /****************************************************\
    * Fun03 Sec03 Sub01:
    *   - get input + start loop
    \****************************************************/

    seqSTPtr->lenSeqUL = 0;

    if(seqSTPtr->seqStr != 0)
       *seqSTPtr->seqStr = '\0';
    else
    { /*Else: allocate memory*/
       seqSTPtr->seqStr =
          malloc((extraBuffUS + 9) * sizeof(schar));

       seqSTPtr->lenSeqBuffUL = extraBuffUS;
    } /*Else: allocate memory*/

    while(1 == 1)
    { /*Loop: read in sequence*/
        tmpStr =
           (schar *)
           fgets(
              (char *)
                 seqSTPtr->seqStr + seqSTPtr->lenSeqUL,
              seqSTPtr->lenSeqBuffUL - seqSTPtr->lenSeqUL,
              (FILE *) faFILE
           );

        if(! tmpStr)
           break; /*will catch EOF on header read in*/

        /************************************************\
        * Fun03 Sec03 Sub02:
        *   - remove white space & check if end of entry
        \************************************************/

        seqSTPtr->lenSeqUL +=
           rmWhite_ulCp(
              seqSTPtr->seqStr + seqSTPtr->lenSeqUL
           );

        /*peek ahead in the file*/
        tmpSC = 0;

        while(tmpSC < 32)
        { /*Loop: find next non-white space character*/
           tmpUL = fread(&tmpSC, 1, 1, (FILE *) faFILE);

           if(tmpUL == 0)
              break;
        } /*Loop: find next non-white space character*/

        if(tmpUL == 0)
           break; /*will catch EOF on header read in*/

        if(tmpSC == '>')
           break; /*finished*/

        /************************************************\
        * Fun03 Sec03 Sub03:
        *   - resize (if needed) and add next character
        \************************************************/

        if(
              seqSTPtr->lenSeqUL
           >= seqSTPtr->lenSeqBuffUL - 10
        ){ /*If: need to allocate more memory*/
           tmpStr =
              realloc(
                 seqSTPtr->seqStr,
                   (
                        seqSTPtr->lenSeqBuffUL
                      + 9
                      + extraBuffUS
                   )
                 * sizeof(schar)
              );

           if(! tmpStr)
              goto memErr_fun03_sec04_sub04;

           seqSTPtr->seqStr = tmpStr;
           seqSTPtr->lenSeqBuffUL +=  extraBuffUS;
        } /*If: need to allocate more memory*/

        /*add peek character to buffer*/
        if(tmpSC > 32)
        { /*If: not white space*/
           seqSTPtr->seqStr[seqSTPtr->lenSeqUL] = tmpSC;
           ++seqSTPtr->lenSeqUL;
        } /*If: not white space*/

        seqSTPtr->seqStr[seqSTPtr->lenSeqUL] = '\0';
    } /*Loop: read in sequence*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun03 Sec04:
    ^   - check for errors and return
    ^   o fun03 sec04 sub01:
    ^     - not EOF or error
    ^   o fun03 sec04 sub02:
    ^     - EOF
    ^   o fun03 sec04 sub03:
    ^     - check if valid fasta sequence
    ^   o fun03 sec04 sub04:
    ^     - memory error
    ^   o fun03 sec04 sub05:
    ^     - return result/error type
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    /****************************************************\
    * Fun03 Sec04 Sub01:
    *   - not EOF or error
    \****************************************************/

    errSC = 0;
    goto seqChecks_fun03_sec04_sub03;

    /****************************************************\
    * Fun03 Sec04 Sub02:
    *   - EOF
    \****************************************************/

    eof_fun03_sec04_sub02:;
    errSC = def_EOF_seqST;
    goto seqChecks_fun03_sec04_sub03;

    /****************************************************\
    * Fun03 Sec04 Sub03:
    *   - check if valid fasta sequence
    \****************************************************/

    seqChecks_fun03_sec04_sub03:;

    if(
          seqSTPtr->seqStr == 0
       || *seqSTPtr->seqStr == '\0'
    ) return def_badLine_seqST | def_fileErr_seqST;

    seqSTPtr->endAlnUL = seqSTPtr->lenSeqUL - 1;
    seqSTPtr->offsetUL = 0;
    goto ret_fun03_sec04_sub05;

    /****************************************************\
    * Fun03 Sec04 Sub04:
    *   - memory error
    \****************************************************/

    memErr_fun03_sec04_sub04:;
    errSC = def_memErr_seqST;
    goto ret_fun03_sec04_sub05;

    /****************************************************\
    * Fun03 Sec04 Sub05:
    *   - return result/error type
    \****************************************************/

    ret_fun03_sec04_sub05:;
    return errSC;
} /*getFaSeq_seqST*/

/*-------------------------------------------------------\
| Fun04: revComp_seqST
|  - reverse complement a sequence
| Input:
|  - seqSTPtr:
|    o pointer to seqSTPtr structure with a sequence to
|      reverse complement
| Output:
|  - Modfies
|    o reverse complements the sequence entry and if
|      reverses the Q-score entry
\-------------------------------------------------------*/
void
revComp_seqST(
  struct seqST *seqSTPtr  /*sequence to reverse comp*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' SOF: Start Of File
    '   o fun04 sec01:
    '     - variable declarations
    '   o fun04 sec02:
    '     - check if have q-score entry
    '   o fun04 sec03:
    '     - reverse complement sequence
    '   o fun04 sec04:
    '     - reverse complement final base and return
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun04 Sec01:
    ^   - variable declarations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    schar *seqStr = seqSTPtr->seqStr;

    schar *endStr =
       seqSTPtr->seqStr + seqSTPtr->lenSeqUL - 1;

    schar *qStr = 0;
    schar *qEndStr = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun04 Sec02:
    ^   - check if have q-score entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(seqSTPtr->qStr != 0 && *seqSTPtr->qStr != '\0')
    { /*If have a Q-score entry*/
       qStr = seqSTPtr->qStr;
       qEndStr = seqSTPtr->qStr + seqSTPtr->lenQUL - 1;
    } /*If have a Q-score entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun04 Sec03:
    ^   - reverse complement sequence
    ^   o fun04 sec03 sub01:
    ^     - reverse complement start base + start loop
    ^   o fun04 sec03 sub02:
    ^     - reverse complement end base
    ^   o fun04 sec03 sub03:
    ^     - swap start and end bases
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /****************************************************\
    * Fun04 Sec03 Sub01:
    *   - reverse complement start base + start loop
    \****************************************************/

    while(endStr > seqStr)
    { /*While I have bases to reverse complement*/
       switch(*seqStr & (~32)) /*Make sure upper case*/
       { /*switch, reverse complement*/
           case 'A': *seqStr = 'T'; break;
           case 'C': *seqStr = 'G'; break;
           case 'G': *seqStr = 'C'; break;
           case 'T': *seqStr = 'A'; break;
           case 'U': *seqStr = 'A'; break;
           case 'W': *seqStr = 'W'; break;
           case 'S': *seqStr = 'S'; break;
           case 'M': *seqStr = 'K'; break;  /*A/C->T/G*/
           case 'K': *seqStr = 'M'; break;  /*T/G->A/C*/
           case 'R': *seqStr = 'Y'; break;  /*A/G->T/C*/
           case 'Y': *seqStr = 'R'; break;  /*T/C->A/G*/
           case 'B': *seqStr = 'V'; break;/*C/G/T->G/C/A*/
           case 'D': *seqStr = 'H'; break;/*G/T/A->C/A/T*/
           case 'H': *seqStr = 'D'; break;/*C/A/T->G/T/A*/
           case 'V': *seqStr = 'B'; break;/*A/C/G->T/G/C*/
           case 'N': *seqStr = 'N'; break;/*A/C/G/T*/
           default: *seqStr = 0;
       } /*switch, reverse complement*/

       /*************************************************\
       * Fun04 Sec03 Sub02:
       *   - reverse complement end base
       \*************************************************/

       switch(*endStr & (~32)) /*Make sure upper case*/
       { /*switch, reverse complement*/
           case 'A': *endStr = 'T'; break;
           case 'C': *endStr = 'G'; break;
           case 'G': *endStr = 'C'; break;
           case 'T': *endStr = 'A'; break;
           case 'U': *endStr = 'A'; break;
           case 'W': *endStr = 'W'; break;
           case 'S': *endStr = 'S'; break;
           case 'M': *endStr = 'K'; break;  /*A/C->T/G*/
           case 'K': *endStr = 'M'; break;  /*T/G->A/C*/
           case 'R': *endStr = 'Y'; break;  /*A/G->T/C*/
           case 'Y': *endStr = 'R'; break;  /*T/C->A/G*/
           case 'B': *endStr = 'V'; break;/*C/G/T->G/C/A*/
           case 'D': *endStr = 'H'; break;/*G/T/A->C/A/T*/
           case 'H': *endStr = 'D'; break;/*C/A/T->G/T/A*/
           case 'V': *endStr = 'B'; break;/*A/C/G->T/G/C*/
           case 'N': *endStr = 'N'; break;/*A/C/G/T*/
           default: *endStr = 0;
       } /*switch, reverse complement*/

       /*************************************************\
       * Fun04 Sec03 Sub03:
       *   - swap start and end bases
       \*************************************************/

       *seqStr ^= *endStr;
       *endStr ^= *seqStr;
       *seqStr ^= *endStr;
       
       if(qStr)
       { /*If I also need to swap Q-scores*/
         *qStr ^= *qEndStr;
         *qEndStr ^= *qStr;
         *qStr ^= *qEndStr;

         ++qStr;
         --qEndStr;
       } /*If I also need to swap Q-scores*/
       
       ++seqStr;
       --endStr;
    } /*While I have bases to reverse complement*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun04 Sec04:
    ^   - reverse complement final base and return
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if ended on same base; if so complement base*/
    if(endStr == seqStr)
    { /*If: ending on same base (not swaped)*/
       switch(*seqStr & (~32)) /*Make sure upper case*/
       { /*switch, reverse complement*/
           case 'A': *seqStr = 'T'; break;
           case 'C': *seqStr = 'G'; break;
           case 'G': *seqStr = 'C'; break;
           case 'T': *seqStr = 'A'; break;
           case 'U': *seqStr = 'A'; break;
           case 'W': *seqStr = 'W'; break;
           case 'S': *seqStr = 'S'; break;
           case 'M': *seqStr = 'K'; break;  /*A/C->T/G*/
           case 'K': *seqStr = 'M'; break;  /*T/G->A/C*/
           case 'R': *seqStr = 'Y'; break;  /*A/G->T/C*/
           case 'Y': *seqStr = 'R'; break;  /*T/C->A/G*/
           case 'B': *seqStr = 'V'; break;/*C/G/T->G/C/A*/
           case 'D': *seqStr = 'H'; break;/*G/T/A->C/A/T*/
           case 'H': *seqStr = 'D'; break;/*C/A/T->G/T/A*/
           case 'V': *seqStr = 'B'; break;/*A/C/G->T/G/C*/
           case 'N': *seqStr = 'N'; break;/*A/C/G/T*/
           default: *seqStr = 0;
       } /*switch, reverse complement*/
    } /*If: ending on same base (not swaped)*/

    return;
} /*revComp_seqST*/

/*-------------------------------------------------------\
| Fun05: blank_seqST
|  - sets values in seqST to blank values
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer to blank
| Output:
|  - Modifies
|    o sets id, sequence, and q-score entries to start
|      with '\0' and the id, sequence, and q-score lengths
|      to 0. This does not change the buffer lengths.
\-------------------------------------------------------*/
void
blank_seqST(
   struct seqST *seqSTPtr
){
   if(seqSTPtr->idStr)
      *seqSTPtr->idStr = '\0';

   seqSTPtr->lenIdUL = 0;

   if(seqSTPtr->seqStr)
      *seqSTPtr->seqStr='\0';

   seqSTPtr->lenSeqUL = 0;

   if(seqSTPtr->qStr)
      *seqSTPtr->qStr = '\0';

   seqSTPtr->lenQUL = 0;

   seqSTPtr->offsetUL = 0;
   seqSTPtr->endAlnUL = 0;
} /*blank_seqST*/

/*-------------------------------------------------------\
| Fun06: init_seqST
|  - sets vlues in seqST to zero
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer to initialize
| Output:
|  - Modifies
|    o all values in seqST to be 0
\-------------------------------------------------------*/
void
init_seqST(
   struct seqST *seqSTPtr
){
   seqSTPtr->idStr = 0;
   seqSTPtr->lenIdBuffUL = 0;

   seqSTPtr->seqStr = 0;
   seqSTPtr->lenSeqBuffUL = 0;

   seqSTPtr->qStr = 0;
   seqSTPtr->lenQBuffUL = 0;

   blank_seqST(seqSTPtr);
} /*init_seqST*/

/*-------------------------------------------------------\
| Fun07: freeStack_seqST
|  - frees the variables in an seqST strucuter
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer that has variables
|      to free
| Output:
|  - Frees
|    o seqSTPtr-idStr and sets to 0
|    o seqSTPtr-seqStr and sets to 0
|    o seqSTPtr-qStr and sets to 0
\-------------------------------------------------------*/
void
freeStack_seqST(
   struct seqST *seqSTPtr
){
   if(seqSTPtr)
   { /*If: I have something to free*/
      free(seqSTPtr->idStr);
      seqSTPtr->idStr = 0;

      free(seqSTPtr->seqStr);
      seqSTPtr->seqStr = 0;

      free(seqSTPtr->qStr);
      seqSTPtr->qStr = 0;

      blank_seqST(seqSTPtr);
   } /*If: I have something to free*/\
} /*freeStack_seqST*/

/*-------------------------------------------------------\
| Fun08: freeHeap_seqST
|  - frees the seqST strucuter
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer to free
| Output:
|  - frees
|    o seqST from memory
\-------------------------------------------------------*/
void
freeHeap_seqST(
   struct seqST *seqSTPtr
){
   freeStack_seqST(seqSTPtr); /*checks if have null*/
   free(seqSTPtr);
} /*freeHeap_seqST*/

/*-------------------------------------------------------\
| Fun09: freeHeapAry_seqST
|  - frees an array of seqST's
| Input:
|  - seqAryST:
|    o pointer to an arrya of seqST structuer to free
|  - numSeqSI:
|    o number of seqSTs in the seqST array
|     (seqSTAryPtr)
| Output:
|  - Frees
|    o seqST from memory
\-------------------------------------------------------*/
void
freeHeapAry_seqST(
   struct seqST *seqAryST,
   signed int numSeqSI
){
   sint siFree = 0;

   for(
      siFree = 0;
      siFree < numSeqSI;
      ++siFree
   ){ /*Loop: Free the individual structures*/
      freeStack_seqST(&seqAryST[siFree]);
   } /*Loop: Free the individual structures*/

   free(seqAryST);
} /*freeHeapAry_seqST*/

/*-------------------------------------------------------\
| Fun10: cpIdEndPad_seqST
|  - copies read id to a buffer and adds in endIdC to
|    the end. If needed, this function will add right
|    padding of spaces to the end.
| Input:
|  - seqSTPtr:
|    o pointer to seqST with read id to copy
|  - buffStr:
|    o c-string to copy read id to
|  - endIdC:
|    o character to add in to mark the end of the read id
|  - padRI:
|    o amount of padding to add to the right of the read
|      id
| Output:
|  - Modifies:
|    o buffStr to hold the sequence id + endIdC. If the
|      id + endIdC is shorter than padRI, the copied id is
|      padded with spaces on the right till it is the same
|      size as padRI.
|  - Returns:
|    o pointer to end of copied id or padding if padding
|      was applied.
\-------------------------------------------------------*/
signed char *
cpIdEndPad_seqST(
   struct seqST *seqSTPtr,/*has read id to copy*/
   signed char *buffStr, /*buffer to add read id to*/
   signed char endIdC,   /*padding; id end (0 to skip)*/
   signed int padRI      /*padding; id start*/
){
   sint idCntI = 0;
   schar *tmpStr = seqSTPtr->idStr;

   if(*tmpStr == '>') ++tmpStr; /*Get of header*/

   while(*tmpStr != '\0')
   { /*While I have a read id to copy over*/
      *buffStr = *tmpStr;
      ++buffStr;
      ++tmpStr;
      ++idCntI;
   } /*While I have a read id to copy over*/

   if(endIdC)
   { /*If adding in a character after the read id*/
      *buffStr = endIdC;
      ++buffStr; /*If not a null*/
      ++idCntI;
   } /*If adding in a character after the read id*/

   while(idCntI < padRI)
   { /*While I have padding to add to the end*/
      *buffStr = ' ';
      ++buffStr;
      ++idCntI;
   } /*While I have padding to add to the end*/

   return buffStr;
} /*cpIdEndPad_seqST*/

/*-------------------------------------------------------\
| Fun11: cp_seqST
|   - copies an seqST structure
| Input:
|   - dupSeqST:
|     o pointer to an seqST to copy to (duplicate)
|   - cpSeqST:
|     o pointer to an seqST to copy
| Output:
|   - Modifies:
|     o dupSeqST to be an copy of cp_seqST
|   - returns:
|     o 0 for no errors
|     o def_memErr_seqST for an memory error
\-------------------------------------------------------*/
signed char
cp_seqST(
   struct seqST *dupSeqST, /*copy to (duplicate)*/
   struct seqST *cpSeqST   /*seqST to copy*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC:
   '   - copies an seqST structure
   '   o fun11 sec01:
   '     - blank and copy aligment settings
   '   o fun11 sec02:
   '     - copy q-score entry
   '   o fun11 sec03:
   '     - copy sequence entry
   '   o fun11 sec04:
   '     - copy read id
   '   o fun11 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^   - blank and copy aligment settings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_seqST(dupSeqST);
   dupSeqST->offsetUL = cpSeqST->offsetUL;
   dupSeqST->endAlnUL = cpSeqST->endAlnUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec02:
   ^   - copy q-score entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   if(cpSeqST->qStr)
   { /*If: have q-score entry*/
      if(
         ! dupSeqST->qStr
         || dupSeqST->lenQBuffUL < cpSeqST->lenQUL
      ){ /*If: have to resize buffer*/

         if(dupSeqST->qStr)
         { /*If: need to resize buffer*/
            free(dupSeqST->qStr);
            dupSeqST->qStr = 0;
         } /*If: need to resize buffer*/

         dupSeqST->lenQBuffUL = cpSeqST->lenQUL;

         dupSeqST->qStr =
            malloc(
               (dupSeqST->lenQBuffUL + 1) * sizeof(schar)
            );

         if(! dupSeqST)
            goto memErr_fun11;
      } /*If: have to resize buffer*/
      
      dupSeqST->lenQUL = cpSeqST->lenQUL;

      cpLen_ulCp(
         dupSeqST->qStr,
         cpSeqST->qStr,
         dupSeqST->lenSeqUL
      );
   } /*If: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec03:
   ^   - copy sequence entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
      ! dupSeqST->seqStr
      || dupSeqST->lenSeqBuffUL < cpSeqST->lenSeqUL
   ){ /*If: have to resize buffer*/

      if(dupSeqST->seqStr)
      { /*If: need to resize buffer*/
         free(dupSeqST->seqStr);
         dupSeqST->seqStr = 0;
      } /*If: need to resize buffer*/

      dupSeqST->lenSeqBuffUL = cpSeqST->lenSeqUL;

      dupSeqST->seqStr =
         malloc(
            (dupSeqST->lenSeqBuffUL + 1) * sizeof(schar)
         );

      if(! dupSeqST)
         goto memErr_fun11;
   } /*If: have to resize buffer*/

   dupSeqST->lenSeqUL = cpSeqST->lenSeqUL;

   cpLen_ulCp(
      dupSeqST->seqStr,
      cpSeqST->seqStr,
      dupSeqST->lenSeqUL
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec04:
   ^   - copy read id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
      ! dupSeqST->idStr
      || dupSeqST->lenIdBuffUL < cpSeqST->lenIdUL
   ){ /*If: have to resize buffer*/

      if(dupSeqST->idStr)
      { /*If: need to resize buffer*/
         free(dupSeqST->idStr);
         dupSeqST->idStr = 0;
      } /*If: need to resize buffer*/

      dupSeqST->lenIdBuffUL = cpSeqST->lenIdUL;

      dupSeqST->idStr =
         malloc(
            (dupSeqST->lenIdBuffUL + 1) * sizeof(schar)
         );

      if(! dupSeqST)
         goto memErr_fun11;
   } /*If: have to resize buffer*/
   
   dupSeqST->lenIdUL = cpSeqST->lenIdUL;

   cpLen_ulCp(
      dupSeqST->idStr,
      cpSeqST->idStr,
      dupSeqST->lenIdUL
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec05:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun11:;
   return def_memErr_seqST;
} /*cp_seqST*/
  
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
