ABC_SCRIPT_DIR=$(dirname $0)
ABC_DIR=${ABC_SCRIPT_DIR}/..;
TOP_DIR=${ABC_DIR}/..;

CLUST_DIR=${ABC_DIR}/data/abc_seurat_cl08_str_dropSubc057_0p8;
PATH_RUN_HOTSPOT=${TOP_DIR}/scripts/run_hotspot.py;
python3 $PATH_RUN_HOTSPOT ${CLUST_DIR}/abc_seurat_cl08_str_dropSubc057_0p8.h5ad ${CLUST_DIR}/hotspot;
