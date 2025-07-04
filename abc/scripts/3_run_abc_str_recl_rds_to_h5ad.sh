ABC_SCRIPT_DIR=$(dirname $0)
ABC_DIR=${ABC_SCRIPT_DIR}/..;
TOP_DIR=${ABC_DIR}/..;

CLUST_DIR=${ABC_DIR}/data/abc_seurat_cl08_str_dropSubc057_0p8;
WRAPPER_PATH=${TOP_DIR}/scripts/conversion_scripts/seurat_to_anndata/wrapper_Seurat_to_AnnData.sh;
bash $WRAPPER_PATH ${CLUST_DIR}/abc_seurat_cl08_str_dropSubc057_0p8.rds ${CLUST_DIR}/abc_seurat_cl08_str_dropSubc057_0p8.h5ad;
