#!/usr/bin/env python
#
# Adding annotations and doing some subsetting before scVI integration. 
# However, output file will not yet be homolog mapped (will do that in the integration script).
# 
import pandas as pd
import anndata as ad
from pathlib import Path

marm_script_dir = Path(__file__).parent.resolve()
marm_data_dir   = Path(marm_script_dir.parent, "data")
path_h5ad       = Path(marm_data_dir, "marm_str_gad_drop2_dropGlut_0p8/marm_str_gad_drop2_dropGlut_0p8.h5ad")
path_clust_anno = Path(marm_data_dir, "marm_str_gad_drop2_dropGlut_0p8/clust_anno.csv")
path_out        = Path(marm_data_dir, "marm_for_scvi.h5ad")


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

# The marmoset "Mixed" cells don't integrate with anything, they form a tight singular cluster.
# >>> adata.shape[0]
# 3379
# 
# dropping 98 cells
adata = adata[adata.obs.clust_anno != 'Mixed', :].copy()

# >>> adata.shape[0]
# 3281

### add annotations that split the LHX8+ from other TAC3
adata.obs['clust_anno_lhx8'] = adata.obs['clust_anno'].astype(str)
adata.obs.loc[adata.obs.seurat_clusters == 8, "clust_anno_lhx8"] = "TAC3/LHX8"


# --------------------------------------------------
# Export
# --------------------------------------------------

adata.write_h5ad(path_out)
