# Use:

Detects DI (defective influenza) reads in sequenced
  samples by primer locations and segment ids.

# License:

This coded is dual licensed under public domain and MIT
  license. Public domain is the default license, unless
  you, your institution/company, your government, or 
  any one else has a problem with it. Then it defaults to
  MIT.

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

Then install diIds

```
# install diIds locally

cd $HOME/Downloads
git clone https://github.com/jeremybuttler/fluDI
cd fluDI/programs/diIdsSrc
make -f mkfile.unix
make -f mkfile.unix install PREFIX="$HOME/local/bin" ```

## global install

Everyone has access, but you need to install as root.

```
if [[ ! -d "$HOME/Downloads" ]]; then
   mkdir "$HOME/Downloads";
fi

# install diIds globally

cd $HOME/Downloads
git clone https://github.com/jeremybuttler/fluDI
cd fluDI/programs/diIdsSrc
make -f mkfile.unix
sudo make -f mkfile.unix install
```

# Run:

You can print the help message with `diIds -h`.

To find the read ids of DI reads
  do `diIds -fq reads.fastq > out.tsv`.

If you want the sequences use seqByID from
  bioTools (fastq files only) or seqkit (fastq/fasta).

# Overview

Uses 12 to 13 conserved nucleotides at 5' and 3'
  ends to identify segments. DI reads are called if the
  region between the primers has 50 or fewer nucleotides
  than the expected segment length.

Problem could come from id sights identify wrong segments.
  So, chimeric reads will be a nightmare. Also the method
  of calling DI reads is rough.

This outputs a tsv (tab deliminated file) with the read
  ids that were kept.

## Output

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

## works

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
