ABC_SCRIPT_DIR=$(dirname $0)
ABC_DIR=${ABC_SCRIPT_DIR}/..;
TOP_DIR=${ABC_DIR}/..;

WRAPPER_PATH=${TOP_DIR}/scripts/conversion_scripts/anndata_to_seurat/wrapper_AnnData_to_Seurat.sh;
bash $WRAPPER_PATH ${ABC_DIR}/data/WMB-10Xv3-STR-raw.h5ad ${ABC_DIR}/data/WMB-10Xv3-STR-raw.rds;
