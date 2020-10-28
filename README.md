# SAMToFasta
Convert a SAM and reference files into a paired FASTA alignment file. This is most useful when you've aligned an amplicon or some small fragments to a known, small contig. This simple program pulls out all aligned reads, having an aligned contig name and position in their SAM record, and uses their CIGAR string to write a equally spaced reference and read pair for each alignment in FASTA format.

## Usage

```
Convert a SAM/BAM to fasta file, given a reference

USAGE:
    samtofasta [FLAGS] --input <FILE> --output <FILE> --ref <FILE>

FLAGS:
    -f, --full_reference    should we try to output the full alignment, i.e. the gaps in the read from the beginning of
                            the contig to the end
    -h, --help              Prints help information
    -V, --version           Prints version information

OPTIONS:
    -i, --input <FILE>     An input SAM/BAM file containing reads already aligned to the reference
    -o, --output <FILE>    the FASTA alignment output file
    -r, --ref <FILE>       The reference we've align our reads to, can have multiple contigs
```

## Build

samtofasta is written in rust, and you can build (after install rust and cloning the repo) by typing:

```
cargo build --release
```

From the main directory. A binary will be put into ```./target/release/samtofasta```
