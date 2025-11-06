#!/bin/bash

: BRAKER3 run with BRAKER2 mode -- no rnaseq -- with compleasm BUSCO completion of Lava Gull repeats softmasked assembly

threads=96

asm=input/lavagull_v0.C.softmasked.fa
proteins=input/Vertebrata_OrthoDB_11.fa

species=LeucoFulig
lineage=aves

outdir="."  # output in current dir

: BRAKER_run.sh puts Augustus and other dirs into the PATH and sets things up to call braker.pl and calls it with the given args
: cores option was changed to threads with BRAKER3

BRAKER_run.sh --genome=$asm        \
              --species=$species   \
              --prot_seq=$proteins \
              --workingdir=$outdir \
              --softmasking        \
              --gff3               \
              --busco_lineage=$lineage         \
              --threads=$threads               |&
tee braker_run.log
