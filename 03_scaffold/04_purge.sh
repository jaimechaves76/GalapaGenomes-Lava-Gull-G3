#!/bin/bash

ref=assembly_v0.9.fasta
reads=lavagull_chopped.fasta

purge_dups.sh $ref $reads
