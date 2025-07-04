#!/usr/bin/env python
#
# This script needs to be run after scVI integration.
# The scVI integration was used to identify relevant cells from the previously
# unlabelled HMBA dataset (which contains cells from diverse subcortical structures). 
# 
# The scVI integration script exported cell barcodes + our new annotations for the HMBA cells,
# which are used here to create the h5ad file. This file is in turn used for correlation and
# SAMap analysis.
#
# During initial analysis, no public annotations were available for the HMBA data,
# but they became available shortly before resubmission. In the scVI integration,
# we found strong agreement between our results and the labels in the HMBA taxonomy.
# 
import pandas as pd
import anndata as ad
from pathlib import Path

hmba_script_dir = Path(__file__).parent.resolve()
hmba_data_dir   = Path(hmba_script_dir.parent, "data")
path_h5ad       = Path(hmba_data_dir, "rna_merged_qc_pass_raw.h5ad")   # [!]
path_cells_anno = Path(hmba_data_dir, "cells_annotated.csv")
path_out        = Path(hmba_data_dir, "marm_hmba.h5ad")

###
### [!] HMBA is prerelease from: 240905_reprocess_and_recluster/data/rna_merged_qc_pass_raw.h5ad
### 
### The barcodes present in the HMBA data folder should make these cells findable in any future public release
### using this script.
###

# --------------------------------------------------
# Imports
# --------------------------------------------------

adata      = ad.read_h5ad(path_h5ad)
cells_anno = pd.read_csv(path_cells_anno)


# --------------------------------------------------
# Subset data and join annotations
# --------------------------------------------------

## Select only cells from the scVI integration, and also join in the cluster annotations created there
cells_anno.set_index(cells_anno.columns[0], inplace=True)
cells_anno.drop('leiden', axis=1, inplace=True)

# >>> adata.shape
# (282806, 35787)
assert all(cells_anno.index.isin(adata.obs.index))
adata = adata[cells_anno.index, :].copy()

# >>> adata.shape
# (12182, 35787)

adata.obs = adata.obs.join(cells_anno)


# --------------------------------------------------
# Export
# --------------------------------------------------

adata.write_h5ad(path_out)
