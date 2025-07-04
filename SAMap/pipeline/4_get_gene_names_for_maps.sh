#!/usr/bin/env bash
#
# For each transcriptome, produce a csv file: first column with the identifier found in the blast mapping, 
# second column with the gene name (symbol) as they would be usually found in the RNA h5ad files. 
# It's not necessary to subset them to match the h5ad files, SAMap will do that when it runs.
# 
DIR_SAMAP=$(dirname $0)/..
cd $DIR_SAMAP
trhg=${DIR_SAMAP}/transcriptomes/human/GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna
trcj=${DIR_SAMAP}/transcriptomes/marm_hmba/GCF_011100555.1_mCalJa1.2.pat.X_cds_from_genomic.fna
trmm=${DIR_SAMAP}/transcriptomes/mouse/GCF_000001635.26_GRCm38.p6_cds_from_genomic.fna

[ -e maps ] || exit 1
[ -e genes_for_maps ] || mkdir genes_for_maps

grep "^>" $trhg | sed -E 's/>(.*) \[gene=(.*)\] (.*)/\1,\2/' | sed 's/\].*//' > genes_for_maps/hg.csv
grep "^>" $trcj | sed -E 's/>(.*) \[gene=(.*)\] (.*)/\1,\2/' | sed 's/\].*//' > genes_for_maps/cj.csv
grep "^>" $trmm | sed -E 's/>(.*) \[gene=(.*)\] (.*)/\1,\2/' | sed 's/\].*//' > genes_for_maps/mm.csv
