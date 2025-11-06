#!/bin/bash

: Lmic1_RagTag    1       39018472        1       W       contig_216      1       39018472        -

awk '
   $5 != "W" { next }
   cur != $1 {
      if (cur != "") printf("\t%s in %d\n", len_last, parts)
      printf("%s", $1)
      cur = $1
      parts = 0
   }
   {
      strand = ($9 == "-") ? "-" : ""
      printf("\t%s%s", $6, strand)
      len_last = $3
      parts++
   }
   END { print "" }
' $1 | sort -k1,1V
