# on this macaque data:
# running hotspot on n=109111 cells (without functioning parallelization...) 
# would take ~8 days, and it's doubtful that's worthwhile.
# 
# so instead, will subset to use only half that number of cells. to ensure
# good representation of all cell types, will subsample the cells within 
# high-resolution clusters.
#
# (the same is done for the developing mouse data in another script)
import anndata as ad
import pandas as pd
from pathlib import Path

devmac_script_dir = Path(__file__).parent.resolve()
devmac_data_dir = Path(devmac_script_dir.parent, "data")
path_h5ad = Path(devmac_data_dir, "MacaqueDevInhibitoryNeurons.h5ad")
path_h5ad_out = Path(devmac_data_dir, "MacaqueDevInhibitoryNeurons_subsample_50perc.h5ad")

adata = ad.read_h5ad(path_h5ad, backed = 'r')

### Want to use the embedding/clustering to guide sampling:
# 
# >>> adata.obs.columns
# Index(['batch', 'batch_name', 'file_name', 'timepoint', 'region', 'class',
#        'hires_leiden', 'leiden', 'n_genes', 'latent_cell_probability', 'phase',
#        'latent_time', 'n_counts'],
#       dtype='object')
#
# >>> adata.obs['hires_leiden'].nunique()
# 124
# >>> adata.obs['hires_leiden'].min()
# 0
# >>> adata.obs['hires_leiden'].max()
# 138
#
###
#
# Finest cluster resolution has 124 clusters (with values 0-138);
# For each cluster, randomly sample half the cells
cell_sample = adata.obs.groupby(by='hires_leiden').sample(frac=0.5).index
adata_sample = adata[adata.obs.index.isin(cell_sample), ]
adata_sample.write_h5ad(path_h5ad_out)
