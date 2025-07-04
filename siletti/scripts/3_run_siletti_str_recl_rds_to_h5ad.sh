SILETTI_SCRIPT_DIR=$(dirname $0);
SILETTI_DIR=${SILETTI_SCRIPT_DIR}/..;
TOP_DIR=${SILETTI_DIR}/..;

CLUST_DIR=${SILETTI_DIR}/data/siletti_neur_roi_sclust_dropAllMeis2_0p8;
WRAPPER_PATH=${TOP_DIR}/scripts/conversion_scripts/seurat_to_anndata/wrapper_Seurat_to_AnnData.sh;
bash $WRAPPER_PATH ${CLUST_DIR}/siletti_neur_roi_sclust_dropAllMeis2_0p8.rds ${CLUST_DIR}/siletti_neur_roi_sclust_dropAllMeis2_0p8.h5ad;
