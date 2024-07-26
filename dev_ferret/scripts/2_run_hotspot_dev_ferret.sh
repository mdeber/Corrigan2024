#!/usr/bin/env bash
DEVFERRET_SCRIPT_DIR=$(dirname $0);
DEVFERRET_DIR=${DEVFERRET_SCRIPT_DIR}/..;
TOP_DIR=${DEVFERRET_DIR}/..;

DATA_DIR=${DEVFERRET_DIR}/data;
PATH_RUN_HOTSPOT=${TOP_DIR}/scripts/run_hotspot.py;
python3 $PATH_RUN_HOTSPOT ${DATA_DIR}/ferret_processed.h5ad ${DATA_DIR}/hotspot;
