ABC_SCRIPT_DIR=$(dirname $0);
ABC_DIR=${ABC_SCRIPT_DIR}/..;

[ -e ${ABC_DIR}/data ] || mkdir ${ABC_DIR}/data;

cd ${ABC_DIR}/data;
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/WMB-10Xv3/20230630/WMB-10Xv3-STR-raw.h5ad;
wget https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/WMB-10X/20230830/views/cell_metadata_with_cluster_annotation.csv;
