#!/bin/bash

threads=32
genome_size=1.32g

reads=lavagull_filt_primary.fastq

flye.sh --nano-hq $reads --out-dir . --scaffold -g $genome_size




###################################################################################################
: "
  --nano-raw path [path ...]
                        ONT regular reads, pre-Guppy5 (<20% error)
  --nano-corr path [path ...]
                        ONT reads that were corrected with other methods (<3% error)
  --nano-hq path [path ...]
                        ONT high-quality reads: Guppy5+ SUP or Q20 (<5% error)
  -g size, --genome-size size
                        estimated genome size (for example, 5m or 2.6g)
  -o path, --out-dir path
                        Output directory
  -t int, --threads int
                        number of parallel threads [1]
  -i int, --iterations int
                        number of polishing iterations [1]
  --scaffold            enable scaffolding using graph [disabled by default]
  --polish-target path  run polisher on the target sequence
  --resume              resume from the last completed stage
  --debug               enable debug output
  -v, --version         show program's version number and exit

Input reads can be in FASTA or FASTQ format, uncompressed or compressed with gz.
Currently, PacBio (CLR, HiFi, corrected) and ONT reads (regular, HQ, corrected) are supported.
Expected error rates are <15% for PB CLR/regular ONT; <5% for ONT HQ, <3% for corrected, and <1% for HiFi.
Note that Flye was primarily developed to run on uncorrected reads.
You may specify multiple files with reads (separated by spaces).

Mixing different read types is not yet supported.
"
###################################################################################################
