#!/usr/bin/env bash
cd $(dirname $0)
PIPELINE_DIR=`pwd`
cd $PIPELINE_DIR/..
PROJECT_DIR=`pwd`
ALIGN_DIR=${PROJECT_DIR}/cellranger
[ -e cellbender ] || mkdir cellbender
cd cellbender
cat ../samples.txt | while read SAMPLE; do

    # skip samples that already failed
    [ -e ${ALIGN_DIR}/${SAMPLE}/QC_FAIL ] && continue
    [ -e $SAMPLE ] || mkdir $SAMPLE
    cd $SAMPLE

    PATH_H5="cellbender_output_filtered.h5"
    PATH_SEUR_H5="cellbender_output_filtered_seurat.h5"

    if [ -e $PATH_SEUR_H5 ]; then
        cd ..
        continue
    fi

    cellbender remove-background \
    --input ${ALIGN_DIR}/${SAMPLE}/${SAMPLE}/outs/raw_feature_bc_matrix.h5 \
    --output cellbender_output.h5 \
    --exclude-feature-types Peaks \
    --cuda;
    
    [ -e "$PATH_H5" ] || exit 1
    ptrepack --complevel 5 "${PATH_H5}:/matrix" "${PATH_SEUR_H5}:/matrix"
    cd ..
done
exit 0
