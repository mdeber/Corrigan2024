#!/usr/bin/env bash
DEVPIG_SCRIPT_DIR=$(dirname $0);
DEVPIG_DIR=${DEVPIG_SCRIPT_DIR}/..;
TOP_DIR=${DEVPIG_DIR}/..;

DATA_DIR=${DEVPIG_DIR}/data;
PATH_RUN_HOTSPOT=${TOP_DIR}/scripts/run_hotspot.py;
python3 $PATH_RUN_HOTSPOT ${DATA_DIR}/pig_processed.h5ad ${DATA_DIR}/hotspot;
