#!/usr/bin/env bash
# 
# Args:
# 1) Path to input h5ad file
# 2) Path for the output Seurat RDS file
#
# This script will make a temporary folder for the intermediate files and 
# delete them when the script finishes.

PATH_AD_TO_MM=$(dirname $0)/AnnData_to_MatrixMarket.py;
PATH_MM_TO_SR=$(dirname $0)/MatrixMarket_to_Seurat.R;

PATH_TMP=$(dirname $2)/tmp;
[ -e $PATH_TMP ] || mkdir $PATH_TMP;

python3 $PATH_AD_TO_MM $1 $PATH_TMP;
Rscript $PATH_MM_TO_SR $PATH_TMP $2;

[ -e "$2" ] || exit 1;
rm -r $PATH_TMP
exit 0;
