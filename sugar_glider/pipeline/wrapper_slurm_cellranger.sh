#!/usr/bin/env bash
#SBATCH -o slurm-cellranger_%j.out
#SBATCH -p all
#SBATCH --time=72:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=8        # 
#SBATCH --mem-per-cpu=4G         # 

# Make sure to run this script within the directory made for each sample, for instance: 
#   .../alignment/data/THE_SAMPLE_NAME;
# 
# and run this script with positional arguments like so:
#   sbatch this_script.sh THE_SAMPLE_NAME FASTQ_DIRECTORY CELLR_SAMPLE_NAME

cellranger count \
--id=$1 \
--fastqs=$2 \
--sample=$3 \
--transcriptome=/PATH/TO/Pbrev_HiC_3kbExt \
--include-introns=true \
--create-bam=true \
--localcores=8 \
--localmem=30;
