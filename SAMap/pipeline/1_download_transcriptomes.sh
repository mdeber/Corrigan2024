#!/usr/bin/env bash
DIR_SAMAP=$(dirname $0)/..
cd $DIR_SAMAP
[ -e transcriptomes ] || mkdir transcriptomes
cd transcriptomes
[ -e human ]       || mkdir human
[ -e mouse ]       || mkdir mouse
[ -e marm_census ] || mkdir marm_census
[ -e marm_hmba ]   || mkdir marm_hmba

cd human
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna.gz"
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.pc_transcripts.fa.gz"
gunzip *gz

cd ../mouse 
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_cds_from_genomic.fna.gz"
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.pc_transcripts.fa.gz"
gunzip *gz

cd ../marm_census
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/663/435/GCF_009663435.1_Callithrix_jacchus_cj1700_1.1/GCF_009663435.1_Callithrix_jacchus_cj1700_1.1_cds_from_genomic.fna.gz"
gunzip *gz

cd ../marm_hmba
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/100/555/GCF_011100555.1_mCalJa1.2.pat.X/GCF_011100555.1_mCalJa1.2.pat.X_cds_from_genomic.fna.gz"
gunzip *gz
