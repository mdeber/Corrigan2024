#!/usr/bin/env bash
# 
# Args:
# 1) Path to input Seurat RDS file
# 2) Path for the output h5ad file
#
# This script will make a temporary folder for the intermediate files and 
# delete them when the script finishes. It only converts the RNA assay
# (see MatrixMarket_to_AnnData.py for details).

PATH_SR_TO_MM=$(dirname $0)/Seurat_to_MatrixMarket.R;
PATH_MM_TO_AD=$(dirname $0)/MatrixMarket_to_AnnData.py;

PATH_TMP=$(dirname $2)/tmp;
[ -e $PATH_TMP ] || mkdir $PATH_TMP;

Rscript $PATH_SR_TO_MM $1 $PATH_TMP;
python3 $PATH_MM_TO_AD $PATH_TMP ${PATH_TMP}/h5ad;

mv ${PATH_TMP}/h5ad/RNA.h5ad $2;
[ -e "$2" ] || exit 1;
rm -r $PATH_TMP
exit 0;
