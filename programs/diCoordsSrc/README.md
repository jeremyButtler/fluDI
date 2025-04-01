# Use:

Detects lare deletions events in each entry in a sam file.

When readding the output tsv file, start_di is the position
  of the first deletion (on reference). end_di is the
  position of the last deletion. For insertions, start_di
  and end_di are the two bases around the insertion.

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

Then install diCoords

```
# install diCoords locally

cd $HOME/Downloads
git clone https://github.com/jeremybuttler/fluDI
cd fluDI/programs/diCoordsSrc
make -f mkfile.unix
make -f mkfile.unix install PREFIX="$HOME/local/bin"
```

## global install

Everyone has access, but you need to install as root.

```
if [[ ! -d "$HOME/Downloads" ]]; then
   mkdir "$HOME/Downloads";
fi

# install diCoords globally

cd $HOME/Downloads
git clone https://github.com/jeremybuttler/fluDI
cd fluDI/programs/diCoordsSrc
make -f mkfile.unix
sudo make -f mkfile.unix install
```

# Run:

You can get the help message with `diCoords -h`. To check
  a sam file for deletions use `diCoords -sam reads.sam`.

# overview

This is DI detection part of diFrag, only the output
  prints the coordinates 

One problem here is that this can not detect missing
  reads. So, make sure your read mapper can find DI reads.

Run by `diCoords -sam reads.sam > out.tsv`. You can get
  a help message with `diCoords -h`.

The output is a tsv with the coordinates of detected
  DI events. Each row represents one detected DI event,
  so there can be multiple rows per read.

## output

- Columns:
  - Column 1 has the read id
  - Column 2 has the reference id
  - Column 3 has the number of DI events (rows) for this
    read
  - Column 4 is the start of one DI event
  - Column 5 is the end of one DI event

### How works:

Basically step 6 of diFrag, with coordinate recording.

1. Read entry from sam file
2. Entry cigar is scanned for large deletions (DI)
   (>= 20) that are at least 12 nucleotides away from
   either end
3. The start and end of all DI events are recorded
4. Print read id, number DI events, and coordinates
   for each DI event

# Updates:

- 2024-08-29:
  - updated to C89 standard, should be windows freindly
- 2024-07-16:
  - Added in a binning program to bin reads in a sam file
    by reference. This will be my first step in consensus
    building.
- 2024-07-15:
  - various small changes to tsv outputs (program
    versions not changed)
  - findDIFrag (only update):
    - changed references from being hardcoded to being
      user provided
      - providing multiple references fixed the problem
        with HA segments (reference was to distant)
    - improved speed (still slow) by filtering out reads
      with few kmers shared with references (less
      alignments)
    - added in options to control min score (as percent),
      min shared kmers (as percent), and kmer length

