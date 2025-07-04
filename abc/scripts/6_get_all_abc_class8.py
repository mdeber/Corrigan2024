#!/usr/bin/env python
#
# Find all cells that are in ABC class 8, regardless of dissection region.
# This will be used as part of an expanded dataset for comparative integration analysis.
# 
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path

# Paths
abc_script_dir = Path(__file__).parent.resolve()
abc_data_dir   = Path(abc_script_dir.parent, "data")
path_abc_anno  = Path(abc_data_dir, "cell_metadata_with_cluster_annotation.csv")
path_out       = Path(abc_data_dir, "abc_cl08_merge.h5ad")

# Imports
abc_anno = pd.read_csv(path_abc_anno)
abc_anno.set_index('cell_label', inplace=True)


# --------------------------------------------------
# Quick check class8 cells
# --------------------------------------------------

# >>> abc_anno.loc[abc_anno['class'] == '08 CNU-MGE GABA', 'feature_matrix_label'].value_counts()
# feature_matrix_label
# WMB-10Xv3-PAL            7974
# WMB-10Xv3-STR            4785
# WMB-10Xv3-HY             1850
# WMB-10Xv2-Isocortex-4    1276
# WMB-10Xv2-Isocortex-1     814
# WMB-10Xv2-HY              496
# WMB-10Xv3-Isocortex-1     432
# WMB-10Xv3-CTXsp           425
# WMB-10Xv2-CTXsp           338
# WMB-10Xv2-OLF             133
# WMB-10Xv2-HPF             118
# WMB-10Xv3-OLF              88
# WMB-10Xv3-Isocortex-2      49
# WMB-10Xv3-HPF              43
# WMB-10Xv3-TH               16
# WMB-10Xv2-TH               12
# Name: count, dtype: int64

### SUBC055, I know the only other datasets we end up with are PAL and HY (v3)

# >>> abc_anno.loc[abc_anno['subclass'] == '055 STR Lhx8 Gaba', 'feature_matrix_label'].value_counts()
# feature_matrix_label
# WMB-10Xv3-PAL      684
# WMB-10Xv3-STR      646
# WMB-10Xv3-HY       161
# WMB-10Xv2-HY        54
# WMB-10Xv3-CTXsp     11
# WMB-10Xv3-TH         9
# WMB-10Xv2-TH         4
# WMB-10Xv2-CTXsp      3
# WMB-10Xv2-OLF        1
# Name: count, dtype: int64

# 
# Just note with cell counts, several places they end up varying in analysis:
# 
# In scVI integration, some differences can come from dropping donors with low cell counts (as donor is modeled). See bottom of script.
# Some cells in the initial totals loacked metadata, implying they were filtered in the atlas analysis set.
#

# --------------------------------------------------
# Export a merged class8 dataset
# --------------------------------------------------

# >>> abc_anno.loc[(abc_anno['class'] == '08 CNU-MGE GABA') & (abc_anno['library_method'] == '10Xv3'), 'feature_matrix_label'].value_counts()
# feature_matrix_label
# WMB-10Xv3-PAL            7974
# WMB-10Xv3-STR            4785
# WMB-10Xv3-HY             1850
# WMB-10Xv3-Isocortex-1     432
# WMB-10Xv3-CTXsp           425
# WMB-10Xv3-OLF              88
# WMB-10Xv3-Isocortex-2      49
# WMB-10Xv3-HPF              43
# WMB-10Xv3-TH               16
# Name: count, dtype: int64

cell_idx = abc_anno.index[(abc_anno['class'] == '08 CNU-MGE GABA') & (abc_anno['library_method'] == '10Xv3')]
ds_list  = list(abc_anno.loc[cell_idx, 'feature_matrix_label'].drop_duplicates().sort_values())

# >>> ds_list
# ['WMB-10Xv3-CTXsp', 'WMB-10Xv3-HPF', 'WMB-10Xv3-HY', 'WMB-10Xv3-Isocortex-1', 'WMB-10Xv3-Isocortex-2', 'WMB-10Xv3-OLF', 'WMB-10Xv3-PAL', 'WMB-10Xv3-STR', 'WMB-10Xv3-TH']

ad_dict = {}
for ds in ds_list:
    path_ad = Path(abc_data_dir, f"{ds}-raw.h5ad")
    assert path_ad.is_file()
    ad_dict[ds] = ad.read_h5ad(path_ad)
    # AnnData object with n_obs × n_vars = __ × __
    # obs: 'cell_barcode', 'library_label', 'anatomical_division_label'
    # var: 'gene_symbol'
    # uns: 'normalization', 'parent', 'parent_layer', 'parent_rows'
    ad_dict[ds] = ad_dict[ds][ad_dict[ds].obs.index.isin(cell_idx), :].copy()

adata     = ad.concat(ad_dict)
adata.var = ad_dict[ds].var

# >>> adata
# AnnData object with n_obs × n_vars = 15662 × 32285
#     obs: 'cell_barcode', 'library_label', 'anatomical_division_label'
#     var: 'gene_symbol'

adata.write_h5ad(path_out)


# --------------------------------------------------
# Cells dropped later for low donor counts
# --------------------------------------------------

# cell_sub = abc_anno.loc[cell_idx, :]
# donor_counts = cell_sub['donor_label'].value_counts()
# cell_sub_sub = cell_sub.loc[cell_sub['donor_label'].isin(donor_counts.index[donor_counts >= 30]), :]

# >>> cell_sub_sub.loc[cell_sub_sub.subclass == '055 STR Lhx8 Gaba', 'feature_matrix_label'].value_counts()
# feature_matrix_label
# WMB-10Xv3-PAL      684
# WMB-10Xv3-STR      638
# WMB-10Xv3-HY       149
# WMB-10Xv3-CTXsp     11

