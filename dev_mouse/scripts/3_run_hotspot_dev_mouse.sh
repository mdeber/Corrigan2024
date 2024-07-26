#!/usr/bin/env bash
DEVMOUSE_SCRIPT_DIR=$(dirname $0);
DEVMOUSE_DIR=${DEVMOUSE_SCRIPT_DIR}/..;
TOP_DIR=${DEVMOUSE_DIR}/..;

DATA_DIR=${DEVMOUSE_DIR}/data;
PATH_RUN_HOTSPOT=${TOP_DIR}/scripts/run_hotspot.py;
python3 $PATH_RUN_HOTSPOT ${DATA_DIR}/mouse_processed_subsample_50perc.h5ad ${DATA_DIR}/hotspot;
