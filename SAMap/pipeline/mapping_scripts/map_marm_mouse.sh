#!/usr/bin/env bash
#SBATCH -o bl_cjmm.out
#SBATCH -p all
#SBATCH --time=96:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=16       # 
#SBATCH --mem-per-cpu=2G         # 

n_threads=16
SCRIPT_DIR="$(dirname "$(realpath "$0")")"
DIR_SAMAP=${SCRIPT_DIR}/../..
tr1=${DIR_SAMAP}/transcriptomes/marm_hmba/GCF_011100555.1_mCalJa1.2.pat.X_cds_from_genomic.fna
n1='cj'
tr2=${DIR_SAMAP}/transcriptomes/mouse/GCF_000001635.26_GRCm38.p6_cds_from_genomic.fna
n2='mm'

cd $DIR_SAMAP
tblastx -query $tr1 -db $tr2 -outfmt 6 -out "maps/${n1}_to_${n2}.txt" -num_threads $n_threads -max_hsps 1 -evalue 1e-6
exit 0
