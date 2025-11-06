#!/bin/bash

: dorado correct requires bgzip not gzip file

threads=96
genome_size=1.32G
min_len=250
concated=lavagull_pass_min_250.fastq.gz

function format_seqsumm_qscore {
   # $2 is read id, $15 length $16 quality score in sequencing_summary.txt
   awk '{print ">"$2 " " $16}' sequencing_summary.txt
}

function add_len_qscore {  # uses sequencing_summary.txt and reads fastq recs from stdin
   bawk -v min_len=$min_len '
      FNR==NR { ar[$name] = $comment; next }  # get read qscore mean

      { slen = length($seq) }

      slen >= min_len {
         print "@" $name " " $comment " " slen " " ar[$name]
         print $seq
         print "+"
         print $qual
      }
   ' <(format_seqsumm_qscore) -
}


############################################################
#                      decompress all                      #
#                add len, qscore to comment                #
#                bgzip as single fastq file                #
############################################################

bgzip -dc basecall_pass/*.fastq.gz |
add_len_qscore                     |
bgzip -@ $threads --output $concated

seq_summary_pass_qscores_by_thous.sh $genome_size $min_len > qual_breakdown.txt
