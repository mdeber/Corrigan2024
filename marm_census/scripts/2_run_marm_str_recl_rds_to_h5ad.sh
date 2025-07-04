MARM_SCRIPT_DIR=$(dirname $0)
MARM_DIR=${MARM_SCRIPT_DIR}/..;
TOP_DIR=${MARM_DIR}/..;

CLUST_DIR=${MARM_DIR}/data/marm_str_gad_drop2_dropGlut_0p8;
WRAPPER_PATH=${TOP_DIR}/scripts/conversion_scripts/seurat_to_anndata/wrapper_Seurat_to_AnnData.sh;
bash $WRAPPER_PATH ${CLUST_DIR}/marm_str_gad_drop2_dropGlut_0p8.rds ${CLUST_DIR}/marm_str_gad_drop2_dropGlut_0p8.h5ad;
