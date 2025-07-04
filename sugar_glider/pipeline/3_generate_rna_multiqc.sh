#!/usr/bin/env bash
cd $(dirname $0)/..
ALIGN_DIR=`pwd`/cellranger
[ -e qc ] || mkdir qc
cd qc
multiqc $ALIGN_DIR
exit 0
