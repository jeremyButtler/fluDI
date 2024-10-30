# Use:

Detects DI (defective influenza) reads in sequenced
  samples.

diIds (segment length reads only) and diFrag both
  use different methods to identify DI reads. So, they
  will likely not agree.

# License:

This coded is dual licensed under public domain and MIT
  license. Public domain is the default license, unless
  you, your institution/company, your government, or 
  any one else has a problem with it. Then it defaults to
  MIT.

# Install:

# General unix

## Local install

The disavatange of a local install is that it is only
  available to you. However, it does not need any extra
  permissions.

First make sure your local enviroment is set up.

```
# make sure have local install and Downloads directory

if [[ ! -d "$HOME/Downloads" ]]; then
   mkdir "$HOME/Downloads";
fi

if [[ ! -d "$HOME/local/bin" ]]; then
   mkdir -p "$HOME/local/bin";
fi

# check if have path setup to local (add to path if not)

localBl="$(grep "$HOME/local/bin" <<< "$PATH")"

if [[ "$localBl" == "" ]]; then
   printf "export \"PATH=$HOME/local/bin:$PATH\"\n" >> $HOME/.bashrc;
fi

export "PATH=$HOME/local/bin:$PATH"

# you can remove $HOME/local/bin from your path by
# deleting the "export PATH=/home/<user name>/local/bin:"
# from your .bashrc file
```

```
cd $HOME/Downloads
git clone https://github.com/jeremybuttler/fluDI
cd fluDI/programs/diCoordsSrc
make -f mkfile.unix
make -f mkfile.unix install PREFIX="$HOME/local/bin"

cd ../diFragSrc
make -f mkfile.unix
make -f mkfile.unix install PREFIX="$HOME/local/bin"

cd ../diIdsSrc
make -f mkfile.unix
make -f mkfile.unix install PREFIX="$HOME/local/bin"
```

## global install

Everyone has access, but you need to install as root.

Initial setup

```
if [[ ! -d "$HOME/Downloads" ]]; then
   mkdir "$HOME/Downloads";
fi
```

build and install programs

```
cd $HOME/Downloads
git clone https://github.com/jeremybuttler/fluDI
cd ../diCoordsSrc
make -f mkfile.unix
sudo make -f mkfile.unix install

cd ../diFragSrc
make -f mkfile.unix
sudo make -f mkfile.unix install

cd ../diIdsSrc
make -f mkfile.unix
sudo make -f mkfile.unix install
```

# Core programs:

## diIds

### diIds Overview and usage

diIds uses 12 to 13 conserved nucleotides at 5' and 3'
  ends to identify segments. DI reads are called if the
  region between the primers has 50 or fewer nucleotides
  than the expected segment length.

Problem could come from id sights identify wrong segments.
  So, chimeric reads will be a nightmare. Also the method
  of calling DI reads is rough.

Run by `diIds -fq reads.fastq > out.tsv`. You can
  get a help message with `diIds -h`.

This outputs a tsv (tab deliminated file) with the read
  ids that were kept.

### diIds Output

- Columns:
  - Column 1 has read id
  - Column 2 has segment read is from
  - Column 3 has classification (vRNA/diRNA/mvRNA)
  - column 4 is the mapped length
    - distance between primer ends
  - column 5 is the read length
  - column 6 is the segments expected length
  - column 7 is
    - both for both primers supporting
    - for for forward primer support only
    - rev for forward primer support only
  - remaining columns are primers stats

### How diIds works

1. Find universal 5' and 3' primer sites in reads
   - reads missing one of these sites are discarded
   - uses the modified (updated/bug fixed) memwater
     from my alnSeq repository
   - or uses kmerFind code (uses waterman code) from
     tbSpoligo, which is a module from freezeTB
2. The segment the read is from is found by comparing the
   bases at the end of the primer sites
   - Reads that have no segment id are discarded
   - Reads that have differing segment ids are discarded
     (unless `-diff` is used).
   - Reads that segments that could only be identified by
     one site are kept (unless `-no-part` is used)
   - Reads were both segments agree are kept
3. The length of the mapping region between primers
   (end of reverse primer - start of forward primer) is
   found
   - reads having a mapping length less then 85% of the
     read are discarded
   - reads with a mapping length that is 10% larger (110%)
     than the expected segment length are discarded
4. Kept read ids are printed out

## diFrag

### diFrag overview

diFrag aligns reads to a set of user provided
  references using a waterman and then searches
  alignments to find number of DIs.

Problems include that this may not always pick the true
  reference. That the low gap extension penalty may
  sometimes create false deletions or may add the last
  few bases to the end. And that it is really slow (this
  one I can probably make better with some time).

Run by `diFrag -fq reads.fastq -ref refs.fasta -out-sam out.sam > out.tsv`.
  You can get a help message with `findIDFrag -h`.

The default output is a tsv file, but you can also output
  a sam file (includes tsv file) with -out-sam. The sam
  file can be passed to stdout by using `-out file.tsv`
  and `-sam-out -`.

The sam file does not have mapping qualities, but does
  include the alignment score as the `AS` (AS:i:) entry
  (number 12).

### diFrag output

- tsv columns:
  - Column 1 has read id
  - Column 2 has reference name
  - Column 3 has classification (diRNA/vRNA)
  - Column 4 has number of detected DI events
  - Column 5 has the direction (forward/reverse)
  - Column 6 has the read length
  - Column 7 has the aligned length (number reference
    bases covered)
  - Columns 8 and 9 have the mapping coordinates
  - Column 10 has the reference length
  - Column 11 has the alignment score
  - Column 12 has the maximum possible alignment score
    for the read
  - Columns 11 to 15 has the number of matches, SNPs,
    insertions, deletions, and masked bases
  - Column 16 has the number of shared kmers with the
    reference (kmer size is in column name)
  - Column 17 and 18 have the mean and median Q-scores
    - this is often over inflated (not sure if I am off or
      if it really is overinflated)

### How diFrag works

1. Counts number of kmers shared between references and
   read
   - The reference (and direction) sharing the most kmers
     with the read is kept
2. Checks to see if (shared kmers) / (reads total kmers)
   is over 40% (otherwise assumes no alignment)
3. Kept reads are mapped to best reference using a
   Waterman Smith alignment with a gap extension penalty
   of -0.05 (match = 5, snp = -4, gap open = -10)
4. Reads with alignment scores under 50% of maximum score
   for a read are discarded
5. The alignment is extracted from the directional matrix
   for the waterman as sam file entry
6. The sam entry cigar is scanned for large deletions (DI)
   (>= 20) that are at least 12 nucleotides away from
   either end
7. The number of detected DI events and other alignment
   information is printed out (+ sam file if requested)

## diCoords

### diCoords overview

getDiCoords is the DI detection part of diFrag, only
    the output prints the coordinates 

One problem here is that this can not detect missing
  reads. So, make sure your read mapper can find DI reads.

Run by `diCoords -sam reads.sam > out.tsv`. You can get
  a help message with `diCoords -h`.

diCoords outputs a tsv with the coordinates of detected
  DI events. Each row represents one detected DI event,
  so there can be multiple rows per read.

### diCoords output

- Columns:
  - Column 1 has the read id
  - Column 2 has the reference id
  - Column 3 has the number of DI events (rows) for this
    read
  - Column 4 is the start of one DI event
  - Column 5 is the end of one DI event

### How diCoords works:

Basically step 6 of diFrag, with coordinate recording.

1. Read entry from sam file
2. Entry cigar is scanned for large deletions (DI)
   (>= 20) that are at least 12 nucleotides away from
   either end
3. The start and end of all DI events are recorded
4. Print read id, number DI events, and coordinates
   for each DI event

# Updates:

- 2024-10-30
  - updateing to keep in sync with my bioTools
    organization
- 2024-09-15:
  - modified code base so c89 complaint
  - shorted names of several programs
    - getDIIds is now diIds
    - findDIFrag is now diFrag
    - getDICoords is now diCoords
  - some minor bug fixes in alignment coordinates
  - side programs moved to bioTools repository, so this
    repository is focused on DI detection programs only
    - [https://github.com/jeremybuttler/bioTools](
       https://github.com/jeremybuttler/bioTools)
- 2024-07-16:
  - Added in a binning program to bin reads in a sam file
    by reference. This will be my first step in consensus
    building.
- 2024-07-15:
  - various small changes to tsv outputs (program
    versions not changed)
  - diFrag (only update):
    - changed references from being hardcoded to being
      user provided
      - providing multiple references fixed the problem
        with HA segments (reference was to distant)
    - improved speed (still slow) by filtering out reads
      with few kmers shared with references (less
      alignments)
    - added in options to control min score (as percent),
      min shared kmers (as percent), and kmer length

# Credits

- The idea of making a DI detection program was from
  Eric Bortz.
- The TB crew (Tara Ness and Bryce Inman) for being part
  of the TB program. This is the freezeTB project that
  contributed the source code for the primer search step.
- Family for being there.
