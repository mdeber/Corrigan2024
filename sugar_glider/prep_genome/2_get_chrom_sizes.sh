#!/usr/bin/env bash
# run with command line arg of fasta file path
in=$1
[ -e ${in}.fai ] || samtools faidx $in
bname=$(basename "${in%.*}")
cut -f1,2 ${in}.fai > ${bname}_chrom.sizes
