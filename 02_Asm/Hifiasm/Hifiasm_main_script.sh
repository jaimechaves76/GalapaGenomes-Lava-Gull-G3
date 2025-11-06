#!/bin/bash

# change as desired to more threads if they can be spared
threads=64

# lex lower of TTAGGG vertebrate telomere motif
telo=CCCTAA

# this version of hifiasm will do its own ONT read error-correction with --ont arg set
adapter_trimmed_ONT_reads=LavaGull_primary.fastq

hifiasm.sh -t $threads --telo-m $telo --ont $adapter_trimmed_ONT_reads
