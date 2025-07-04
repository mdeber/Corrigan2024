#!/usr/bin/env bash
DEVMAC_SCRIPT_DIR=$(dirname $0);
DEVMAC_DIR=${DEVMAC_SCRIPT_DIR}/..;
TOP_DIR=${DEVMAC_DIR}/..;

WRAPPER_PATH=${TOP_DIR}/scripts/conversion_scripts/anndata_to_seurat/wrapper_AnnData_to_Seurat.sh;
# run with counts_only=True
bash $WRAPPER_PATH ${DEVMAC_DIR}/data/MacaqueDevInhibitoryNeurons.h5ad ${DEVMAC_DIR}/data/MacaqueDevInhibitoryNeurons.rds 'True';
