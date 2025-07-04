MARM_SCRIPT_DIR=$(dirname $0)
MARM_DIR=${MARM_SCRIPT_DIR}/..;
TOP_DIR=${MARM_DIR}/..;

CLUST_DIR=${MARM_DIR}/data/marm_str_gad_drop2_dropGlut_0p8;
PATH_RUN_HOTSPOT=${TOP_DIR}/scripts/run_hotspot.py;
python3 $PATH_RUN_HOTSPOT ${CLUST_DIR}/marm_str_gad_drop2_dropGlut_0p8.h5ad ${CLUST_DIR}/hotspot;
