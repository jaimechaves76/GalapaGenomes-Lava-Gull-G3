#!/bin/bash

bawk '
   FILENUM == 1 {
      for (f = 5; f <= NF-1; f++) {
         map[$name] = map[$name] " " $f
         ctg = $f; l = length(ctg)
         rc = substr(ctg, l, 1)=="-" ? 1 : 0
         if (rc) { ctg = substr(ctg, 1, l-1) }
         used_in_scaffolding[ctg]++
      }
      next
   }

   FILENUM == 2 {
      recs[$name] = $seq
      if ( ! ($name in used_in_scaffolding) ) {
         print ">" $name " " length($seq)
         print $seq
      }
   }

   END {
      for (n = 1; n <= 100; n++) Ns = Ns "N"  # for gap between scaffold components

      for (scaff in map) {
         n = split(map[scaff], ar, " ")
         print ">" scaff " " map[scaff]
         for (c = 1; c <= n; c++) {
            ctg = ar[c]; l = length(ctg)
            rc = substr(ctg, l, 1)=="-" ? 1 : 0
            if (rc) { ctg = substr(ctg, 1, l-1) }
            ctg_seq = toupper(recs[ctg])
            if (rc) revcomp(ctg_seq)

            print ctg_seq
            if (c < n)
               print Ns

            # print scaff, ar[c], ctg, rc, length(ctg_seq)
         }
      }
   }
' <(prefix_gt new_scaffold_map.tsv) assembly.fasta |
bawk '{print length($seq) " " $0}' |
sort -k2,2V                        |
cawk '{ print ">" $2 " " $1 " " fldcat(4, NF, " "); print $3 }
'
