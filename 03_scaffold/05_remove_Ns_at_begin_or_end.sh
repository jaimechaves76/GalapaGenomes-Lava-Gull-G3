function remove_Ns_at_begin_or_end {
   bawk '
      # looks like purged_asm left Ns at seq begin and/or seq end for scaffolds it modified, remove them
      { sub("N+$", "", $seq); sub("^N+", "", $seq) }

      { print ">" $name; print $seq }
   ' $asm
}

asm=$1
remove_Ns_at_begin_or_end $asm | seqfold 120
