#!/usr/bin/env bash
DEVMOUSE_SCRIPT_DIR=$(dirname $0);
DEVMOUSE_DIR=${DEVMOUSE_SCRIPT_DIR}/..;
TOP_DIR=${DEVMOUSE_DIR}/..;

WRAPPER_PATH=${TOP_DIR}/scripts/conversion_scripts/anndata_to_seurat/wrapper_AnnData_to_Seurat.sh;
bash $WRAPPER_PATH ${DEVMOUSE_DIR}/data/mouse_processed.h5ad ${DEVMOUSE_DIR}/data/mouse_processed.rds;
