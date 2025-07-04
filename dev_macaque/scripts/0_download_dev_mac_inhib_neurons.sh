#!/usr/bin/env bash
DEVMAC_SCRIPT_DIR=$(dirname $0);
DEVMAC_DIR=${DEVMAC_SCRIPT_DIR}/..;

[ -e cd ${DEVMAC_DIR}/data ] || mkdir cd ${DEVMAC_DIR}/data;

cd ${DEVMAC_DIR}/data;
wget https://cells-test.gi.ucsc.edu/dev-inhibitory-neurons/macaque/MacaqueDevInhibitoryNeurons.h5ad;