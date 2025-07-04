#!/usr/bin/env python
#
# Adding annotations and doing some subsetting before scVI integration. 
# However, output file will not yet be homolog mapped (will do that in the integration script).
# 
import pandas as pd
import anndata as ad
from pathlib import Path

siletti_script_dir = Path(__file__).parent.resolve()
siletti_data_dir   = Path(siletti_script_dir.parent, "data")
path_h5ad          = Path(siletti_data_dir, "siletti_neur_roi_sclust_dropAllMeis2_0p8/siletti_neur_roi_sclust_dropAllMeis2_0p8.h5ad")
path_clust_anno    = Path(siletti_data_dir, "siletti_neur_roi_sclust_dropAllMeis2_0p8/clust_anno.csv")
path_out           = Path(siletti_data_dir, "siletti_for_scvi.h5ad")


# --------------------------------------------------
# Imports
# --------------------------------------------------

adata      = ad.read_h5ad(path_h5ad)
clust_anno = pd.read_csv(path_clust_anno)


# --------------------------------------------------
# Merge 'clust_anno'
# --------------------------------------------------

adata.obs = adata.obs.join(
    clust_anno.rename({"seurat_clusters_anno" : "clust_anno"}, axis = 1).set_index("seurat_clusters"),
    on = "seurat_clusters"
)

# --------------------------------------------------
# Some subsetting
# --------------------------------------------------

#
# The siletti "Hypothalamus" cells don't belong, and can integrate with other hypothalamus cells.
# However, the siletti "Mixed" cells integrate with marmoset HMBA 'STR SST RSPO2 GABA', so were kept in another version.
# But this time we're also dropping Siletti 'Mixed'. Even though we wanted those HMBA cells, they do just integrate with
# that Siletti cluster and nothing else in the other species, and they created a problem where they dragged over a few 
# cells from 'STR SST ADARB2 GABA' and created disjoint leiden clusters.
#
# We don't have any space or bandwidth to be discussing all these minor variations that don't affect any of our conclusions,
# we will try to present the most honest and representative version of integration analysis.
#
# >>> adata.shape[0]
# 6617
# 
# dropping 339 cells
adata = adata[~adata.obs.clust_anno.isin(['Hypothalamus', 'Mixed']), :].copy()

# >>> adata.shape[0]
# 6042

# --------------------------------------------------
# Export
# --------------------------------------------------

adata.write_h5ad(path_out)
