SILETTI_SCRIPT_DIR=$(dirname $0);
SILETTI_DIR=${SILETTI_SCRIPT_DIR}/..;
TOP_DIR=${SILETTI_DIR}/..;

CLUST_DIR=${SILETTI_DIR}/data/siletti_neur_roi_sclust_dropAllMeis2_0p8;
PATH_RUN_HOTSPOT=${TOP_DIR}/scripts/run_hotspot.py;
python3 $PATH_RUN_HOTSPOT ${CLUST_DIR}/siletti_neur_roi_sclust_dropAllMeis2_0p8.h5ad ${CLUST_DIR}/hotspot;
