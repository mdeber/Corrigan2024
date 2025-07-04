# running hotspot on the developing mouse data (n=71023 cells) takes days;
# in this script, will cluster the mouse data at high resolution and then
# randomly sample 50% of cells from each cluster (output n=35507 cells)
#
# (the same is also done for the developing macaque data in another script)

import anndata as ad
import scanpy as sc
import pandas as pd
from pathlib import Path

devmouse_script_dir = Path(__file__).parent.resolve()
devmouse_data_dir = Path(devmouse_script_dir.parent, "data")
path_h5ad = Path(devmouse_data_dir, "mouse_processed.h5ad")
path_h5ad_out = Path(devmouse_data_dir, "mouse_processed_subsample_50perc.h5ad")

adata = ad.read_h5ad(path_h5ad, backed = 'r')

### data & existing clusters
# >>> adata.obs.columns
# Index(['n_genes', 'percent_mito', 'percent_ribo', 'n_counts', 'doublet_score',
#        'predicted_doublet', 'doublet_info', 'batch_name', 'batch', 'timepoint',
#        'leiden', 'leiden_subcluster_CRABPY', 'leiden_subcluster_CRABPY2',
#        'class'],
#       dtype='object')
# >>> adata.obs.leiden.nunique()
# 14

sc.tl.leiden(adata, resolution=10, key_added='leiden_hires', n_iterations=-1)
# >>> adata.obs.leiden_hires.nunique()
# 157

cell_sample = adata.obs.groupby(by='leiden_hires').sample(frac=0.5).index
adata_sample = adata[adata.obs.index.isin(cell_sample), ]
adata_sample.write_h5ad(path_h5ad_out)
