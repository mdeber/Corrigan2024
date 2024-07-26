#!/usr/bin/env bash
DEVPIG_SCRIPT_DIR=$(dirname $0);
DEVPIG_DIR=${DEVPIG_SCRIPT_DIR}/..;
TOP_DIR=${DEVPIG_DIR}/..;

WRAPPER_PATH=${TOP_DIR}/scripts/conversion_scripts/anndata_to_seurat/wrapper_AnnData_to_Seurat.sh;
bash $WRAPPER_PATH ${DEVPIG_DIR}/data/pig_processed.h5ad ${DEVPIG_DIR}/data/pig_processed.rds;
