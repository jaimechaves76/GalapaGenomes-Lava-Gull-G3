#!/bin/bash

threads=32
genome_size=1.32G

input=lavagull_pass_min_250.fastq.gz
output=lavagull_chopped.fastq
readid_len_qscore=$(replace_ext lavagull_chopped.fastq readid_len_qscore.txt)

log=lavagull_porechop.log

# we will have porechop write to stdout so we can add len and qscore to the comment

porechop.sh --version | tee $log
porechop.sh --threads $threads --verbosity 1 -i $input |
add_len_nanoplot_qscore_to_fastq.sh > $output

if [ -s $output ]; then
   fastq_readid_qscore_len_fields.sh $output > $readid_len_qscore
   seq_summary_pass_qscores_by_thous.sh $readid_len_qscore $genome_size > qual_breakdown.txt
fi
