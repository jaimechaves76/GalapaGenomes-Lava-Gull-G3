#!/bin/bash

ragtag.sh scaffold Larus_mich_1.1_genomic.fna assembly.fasta -t 96 -o ./assembly_Lmich_ragtag_output
ragtag.sh scaffold Larus_mich_1.1_genomic.fna purged_asm.fa -t 96 -o ./purged-asm_Lmich_ragtag_output

ragtag.sh scaffold Larus_fuscus_1.1_genomic.fna assembly.fasta -t 96 -o ./assembly_Lfuscus_ragtag_output
ragtag.sh scaffold Larus_fuscus_1.1_genomic.fna purged_asm.fa -t 96 -o ./purged-asm_Lfuscus_ragtag_output
