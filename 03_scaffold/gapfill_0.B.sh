#!/bin/bash

asm=assembly_v0.B.fasta

reads=lavagull_chopped.fasta
threads=32

quartet_gapfiller.py -d $asm -g $reads -t $threads
