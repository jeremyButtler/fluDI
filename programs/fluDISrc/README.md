# USE:

Scan a fastq file for vRNA and diRNA reads using diIDs and
  diFrag.

# License:

Dual licensed under Unlicense and MIT. Unlicense is
  the primary lincense. However, this defaults to MIT
  if the Unlicense is not recognized or not wanted (for
  any reason).

# Install:

## Unix:

```
if [[ ! -d "/usr/local/bin" ]]; then
   sudo mkdir -p "/usr/local/bin";
fi # Mac does not have directory, but is still in path

make -f mkfile.unix
sudo make -f mkfile.unix install
```

## Windows:

You can do this, but it is not recomended (is command
  line program).

1. Install the windows developer console (comes with
   visual studio (20Gb install; so if can just get
   developer console is better)).
2. Download fluDI from github and unzip (do Downloads)
3. Open a developer console
   - in visual studio: view->developer console
   - start menu->visusal studio->developer console
4. enter `cd $HOME$\Downloads\fluDI\programs\fluDISrc`
5. enter `nmake \F mkfile.win`
6. Close developer (`exit) console and copy the fluDI exe
   file in `fluDI\programs\fluDISrc` to your desired
   location.

# Run:

You can get the help message with `fluDI -h`.

To run fluDI you will need a fastq file with reads
  (preferably long reads) and a set of flu references for
  each segment (see ../diFragSrc/fluGenbankRef.fa) in a
  fasta file. The header of segment needs to start with
  the segment id/number (ex `>PB2_ref` or `1_ref`). You
  can then run fluDI with
  `fluDI -fq reads.fastq -ref refs.fasta -prefix goodName`

- Segments:
  - PB2 is 1
  - PB1 is 2
  - PA is 3
  - HA is 4
  - NP is 5
  - NA is 6
  - M is 7
  - NS is 8

- The output of fluDI will be sam files and tsv files.
  - delete-report shows vRNA and diRNA reads found by
    fluDI, the identified segment, and the classification
    (diRNA or vRNA) called by diFrag and diIds.
  - prefix-frag.tsv is the diFrag report
    - A `*` means the value could not be found. It happens
     when only one program can classify a read
  - prefix-di-IDs.tsv is the diIds report
  - prefix-diRna.sam is a sam file with reads classifed as
    a DI read
    - the call:s: flag tells if the classification was
      from diIds (diID), diFrag (frag), or both programs
  - prefix-vRna.sam is a sam file with reads classified
    as vRna
    - the call:s: flag tells if the classification was
      from diFrag (frag), or both diFrag and diIds (diIds
      is not allowed to classify vRna on its own)

# How works:

This is not perfect and you will likey have to scan
  through the results for mis-classified reads. However,
  it is a start.

## diFrag:

Read 7mer counts are compared to the list of references
  and the reference with the closest 7mer count is
  selected. The read is then aligned to the selected
  reference with a Waterman alignment. Reads with less
  than 90% of the maximum score are discareded, but the
  alingment is kept for diFag and diIds merging. Each
  kept read is then scaned for large deletions (>20). If
  large deletions are present then the read is classified
  as diRNA, else it is vRNA.

## diIds:

The forward omni primer sites are found in reads with a
  kmer scan using 5mers. Reads with only one primer site
  are discarded. The flu segment the read is mapped to is
  then found using both primers sites. The aligned length
  is found by finding the length between the primer sites.
  If the aligned length is less than the expected segment 
  length by 20 bases, then the reads is classified as DI.
  Otherwise, it is classifed as vRNA.

## Combing reports:

The results diFrag and diIds are combined together to
  make a final classification. Reads with different
  flu segments reported by diFrag and diIds are discarded.
  While reads with that idFrag did not classify, but did
  align, are kept, so long as the found segment agrees
  with diIds. Reads that were not aligned by diFrag are
  discarded. Though the diIds classification is still
  reported.

Reads that could only be classifed by diFrag are printed
  to their respective sam files (prefix-diRna.sam or
  prefix-vRna.sam). However, for cases when only diIds
  could classify reads, but diFrag could align the read,
  only the diRNA classifications are printed. For cases
  were diFrag and diIds could classify a read, but
  disagreeded about the vRNA or diRNA classifaction, then
  the read is considered diRNA.

Finally, edClust is run on the diRNA classified reads and
  vRNA classified reads to build any consensuses.

# Updates:
