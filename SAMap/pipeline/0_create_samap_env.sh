#!/usr/bin/env bash
echo "Don't run this script without editing"
exit 0
CDIR=/PATH/TO/envs/SAMap
conda create -p $CDIR -c conda-forge python=3.9 numpy=1.23.5 pip pybind11 h5py=3.8.0 leidenalg python-igraph texttable
conda activate $CDIR
conda install -c bioconda blast=2.9.0
mkdir ${CDIR}/samap
git clone https://github.com/atarashansky/SAMap.git ${CDIR}/samap
python -m pip install ${CDIR}/samap
conda install sympy
