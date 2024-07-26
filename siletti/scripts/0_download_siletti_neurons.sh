SILETTI_SCRIPT_DIR=$(dirname $0);
SILETTI_DIR=${SILETTI_SCRIPT_DIR}/..;

cd ${SILETTI_DIR}/data;
wget https://storage.googleapis.com/linnarsson-lab-human/Neurons.h5ad;
wget https://github.com/linnarsson-lab/adult-human-brain/raw/main/tables/cluster_annotation.xlsx
wget https://github.com/linnarsson-lab/adult-human-brain/raw/main/tables/subcluster_annotation.xlsx
