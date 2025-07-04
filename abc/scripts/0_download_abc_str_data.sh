#!/usr/bin/env bash
ABC_SCRIPT_DIR=$(dirname $0)
ABC_DIR=${ABC_SCRIPT_DIR}/..

[ -e ${ABC_DIR}/data ] || mkdir ${ABC_DIR}/data

cd ${ABC_DIR}/data
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/WMB-10X/20230830/views/cell_metadata_with_cluster_annotation.csv

# most analysis only uses cells from the striatal dataset:
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-STR-raw.h5ad

# however, for later comparative integration analyses, also getting all cells from class8 (regardless of annotated dissection region);
# these datasets were determined in script: 6_get_all_abc_class8.py
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-PAL-raw.h5ad
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-HY-raw.h5ad
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-1-raw.h5ad
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-CTXsp-raw.h5ad
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-OLF-raw.h5ad
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-Isocortex-2-raw.h5ad
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-HPF-raw.h5ad
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-TH-raw.h5ad
