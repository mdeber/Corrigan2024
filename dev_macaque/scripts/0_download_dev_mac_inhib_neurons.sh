#!/usr/bin/env bash
DEVMAC_SCRIPT_DIR=$(dirname $0);
DEVMAC_DIR=${DEVMAC_SCRIPT_DIR}/..;

cd ${DEVMAC_DIR}/data;
wget https://cells-test.gi.ucsc.edu/dev-inhibitory-neurons/macaque/MacaqueDevInhibitoryNeurons.h5ad;