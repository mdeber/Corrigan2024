#!/usr/bin/env bash
DIR_SAMAP=$(dirname $0)/..
cd $DIR_SAMAP

# for generating the original blastdb files, don't want to launch multiple jobs that act on the same
# files (learned that the hard way), and this only takes a second:
trhg=${DIR_SAMAP}/transcriptomes/human/GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna
trcj=${DIR_SAMAP}/transcriptomes/marm_hmba/GCF_011100555.1_mCalJa1.2.pat.X_cds_from_genomic.fna
trmm=${DIR_SAMAP}/transcriptomes/mouse/GCF_000001635.26_GRCm38.p6_cds_from_genomic.fna
makeblastdb -in $trhg -dbtype 'nucl'
makeblastdb -in $trcg -dbtype 'nucl'
makeblastdb -in $trmm -dbtype 'nucl'

# but for actually generating the tblastx maps, we want to launch all these jobs in parallel
[ -e maps ] || mkdir maps
cd scripts
ls mapping_scripts | while read SCRIPT; do
    # sbatch mapping_scripts/${SCRIPT} ## what i did (and why it's laid out like this)
    bash mapping_scripts/${SCRIPT}
done
