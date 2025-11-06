#!/bin/bash

: use new RepeatMasker version from Nov 2024 4.1.7 and aves repeats from Dfam38 -- many more than Dfam37
: use much smaller curated aves set to compare results

threads=96
asm=input/lavagull_v0.C.fasta
repeats_fasta=input/denovo_and_aves_curated_repeats_Dfam38.fa

repeatmasker_run.sh $asm $repeats_fasta $threads
