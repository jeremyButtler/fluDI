# Use:

Finds DI fragments in a fastq file.

# License

Primary is unlicense (public domain), but if publice
  domain does not work or is not desired it defaults to
  MIT (alternate).

# Build

Each make file (mkfile) is prefixed by the OS type
  (mkfile.OS) it builds.

- Make files
  - mkfile.unix: general unix (including Mac)
  - mkfile.static: static builds for Linux/BSD
  - mkfile.win: for windows

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

Then install diFrag

```
# install idFrag locally

cd $HOME/Downloads
git clone https://github.com/jeremybuttler/fluDI
cd fluDI/programs/diFragSrc
make -f mkfile.unix
make -f mkfile.unix install PREFIX="$HOME/local/bin"
```

## global install

Everyone has access, but you need to install as root.

```
if [[ ! -d "$HOME/Downloads" ]]; then
   mkdir "$HOME/Downloads";
fi

# install diFrag globally

cd $HOME/Downloads
git clone https://github.com/jeremybuttler/fluDI
cd fluDI/programs/diFragSrc
make -f mkfile.unix
sudo make -f mkfile.unix install
```

# Run:

To get the help message for diFrag do `diFrag -h`.

To run diFrag you will need a set of references and a
  fastq file with reads to map.

```
diFrag -ref references.fasta -fq reads.fastq
```

The output will be a tsv with the list of read ids,
  reference, if vRNA or diRNA, and mapping inforamtion.
  If you want to extract the reads use seqById from
  bioTools
  [https://github.com/jeremybuttler/bioTools](
   https://github.com/jeremybuttler/bioTools) or seqkit
  (use `cut -f 1 out.tsv | tail -n+2` to prepare file for
  seqkit).

You removing everything but diRNA
  using `grep "diRNA" stats.tsv > diReads.tsv`.

You removing everything but vRNA (full length)
  using `grep "vRNA" stats.tsv > diReads.tsv`.

# overview

findDIFrag aligns reads to a set of user provided
  references using a waterman and then searches
  alignments to find number of DIs.

Problems include that this may not always pick the true
  reference. That the low gap extension penalty may
  sometimes create false deletions or may add the last
  few bases to the end. And that it is really slow (this
  one I can probably make better with some time).

The default output is a tsv file, but you can also output
  a sam file (includes tsv file) with -out-sam. The sam
  file can be passed to stdout by using `-out file.tsv`
  and `-sam-out -`.

The sam file does not have mapping qualities, but does
  include the alignment score as the `AS` (AS:i:) entry
  (number 12).

## output

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

## how works

1. Counts number of kmers shared between references and
   read
   - The reference (and direction) sharing the most kmers
     with the read is kept
2. Checks to see if (shared kmers) / (reads total kmers)
   is over 40% (otherwise assumes no alignment)
3. Kept reads are mapped to best reference using a
   Waterman Smith alignment with a gap extension penalty
   of -0.1 (match = 5, snp = -4, gap open = -10)
4. Reads with alignment scores under 50% of maximum score
   for a read are discarded
5. The alignment is extracted from the directional matrix
   for the waterman as sam file entry
6. The sam entry cigar is scanned for large deletions (DI)
   (>= 20) that are at least 12 nucleotides away from
   either end
7. The number of detected DI events and other alignment
   information is printed out (+ sam file if requested)
