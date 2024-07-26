#!/usr/bin/env bash
DEVFERRET_SCRIPT_DIR=$(dirname $0);
DEVFERRET_DIR=${DEVFERRET_SCRIPT_DIR}/..;
TOP_DIR=${DEVFERRET_DIR}/..;

WRAPPER_PATH=${TOP_DIR}/scripts/conversion_scripts/anndata_to_seurat/wrapper_AnnData_to_Seurat.sh;
bash $WRAPPER_PATH ${DEVFERRET_DIR}/data/ferret_processed.h5ad ${DEVFERRET_DIR}/data/ferret_processed.rds;
