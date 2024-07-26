#!/usr/bin/env bash
DEVMAC_SCRIPT_DIR=$(dirname $0);
DEVMAC_DIR=${DEVMAC_SCRIPT_DIR}/..;
TOP_DIR=${DEVMAC_DIR}/..;

DATA_DIR=${DEVMAC_DIR}/data;
PATH_RUN_HOTSPOT=${TOP_DIR}/scripts/run_hotspot.py;
python3 $PATH_RUN_HOTSPOT ${DATA_DIR}/MacaqueDevInhibitoryNeurons_subsample_50perc.h5ad ${DATA_DIR}/hotspot;
