#!/usr/bin/env bash
cd $(dirname $0)
PIPELINE_DIR=`pwd`
cd $PIPELINE_DIR/..
PROJECT_DIR=`pwd`
ALIGN_DIR=${PROJECT_DIR}/cellranger
FQ_DIR=${PROJECT_DIR}/fastq

# generate samples.txt
ls -1 $FQ_DIR | sed 's/_S[0-9]*_L[0-9]*_.*//' | uniq > samples.txt

[ -e $ALIGN_DIR ] || mkdir $ALIGN_DIR
cd $ALIGN_DIR

# for each sample, slurm launch cellranger
cat ../samples.txt | while read SAMPLE; do
    [ -d "$SAMPLE" ] || mkdir "$SAMPLE"
    cd "$SAMPLE"
    CELLR_SNAME=$(ls -1 "$FQ_DIR" | grep "$SAMPLE" | sed 's/_S[0-9]*_L[0-9]*_.*//' | uniq)
    if [ ! -e "${SAMPLE}/outs" ]; then
        [ -e "${SAMPLE}/_lock" ] && rm "${SAMPLE}/_lock"
        sbatch -J "cellranger_${SAMPLE}" "${PIPELINE_DIR}/wrapper_slurm_cellranger.sh" "$SAMPLE" $FQ_DIR $CELLR_SNAME
    fi
    cd ..
done
