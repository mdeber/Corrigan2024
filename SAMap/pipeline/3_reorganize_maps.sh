#!/usr/bin/env bash
# actually want to keep these organized like SAMap expects
DIR_SAMAP=$(dirname $0)/..
cd $DIR_SAMAP/maps
[ -e cjmm ] && exit 1
mkdir cjmm
mkdir hgcj
mkdir hgmm
mv cj_to_mm.txt cjmm
mv mm_to_cj.txt cjmm
mv hg_to_cj.txt hgcj
mv cj_to_hg.txt hgcj
mv hg_to_mm.txt hgmm
mv mm_to_hg.txt hgmm
