#!/usr/bin/env python3
#
# Create the expanded ABC dataset used for scVI integration. This has other useful classes for comparative
# analysis that aren't in the other targeted subset of striatal GABAergic interneurons.
#
# However, output file will not yet be homolog mapped (will do that in the integration script).
# 
import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path


# --------------------------------------------------
# Paths
# --------------------------------------------------

abc_script_dir = Path(__file__).parent.resolve()
abc_data_dir   = Path(abc_script_dir.parent, "data")

# file with all ABC striatal cells (of any class)
path_abc_str   = Path(abc_data_dir, "WMB-10Xv3-STR-raw.h5ad")

# file with all ABC cells in class8 (including striatal, which are also present in abc_str file);
# will be including the non-striatal cells from subclass 055
path_abc_cl8   = Path(abc_data_dir, "abc_cl08_merge.h5ad")

# file with the ABC cells in our initial analysis set (re-processed and re-annotated)
path_abc_cl8_str_anno = Path(abc_data_dir, "abc_seurat_cl08_str_dropSubc057_0p8/abc_seurat_cl08_str_dropSubc057_0p8.h5ad")

# original complete annotations from ABC atlas
path_abc_anno = Path(abc_data_dir, "cell_metadata_with_cluster_annotation.csv")

# our re-annotated clusters (only striatal GABAergic interneurons)
path_clust_anno = Path(abc_data_dir, "abc_seurat_cl08_str_dropSubc057_0p8/clust_anno.csv")


# --------------------------------------------------
# Imports
# --------------------------------------------------

print("Importing data...")
abc_str          = ad.read_h5ad(path_abc_str)
abc_cl8          = ad.read_h5ad(path_abc_cl8)
abc_cl8_str_anno = ad.read_h5ad(path_abc_cl8_str_anno)
clust_anno       = pd.read_csv(path_clust_anno)


# --------------------------------------------------
# ABC metadata
# --------------------------------------------------

print("Importing, formatting metadata...")
# Import, join ABC cell metadata
abc_anno = pd.read_csv(path_abc_anno)
abc_anno.set_index('cell_label', inplace=True)
anno_col_keep = [
    "feature_matrix_label",
    "region_of_interest_acronym",
    "library_method",
    "donor_label",
    "donor_genotype",
    "donor_sex",
    "cluster_alias",
    "neurotransmitter",
    "class",
    "subclass", 
    "supertype",
    "cluster"
]
abc_anno = abc_anno.loc[:, abc_anno.columns.isin(anno_col_keep)]


# --------------------------------------------------
# Before joining ABC metadata, combine non-striatal cells from ABC SUBC055 into the ABC striatal data
# --------------------------------------------------

print("Before joining ABC cell metadata, adding in the subc055 cells from outside the striatum...")
cbc_subc055 = abc_anno.index[(abc_anno['subclass'] == '055 STR Lhx8 Gaba') & (abc_anno.index.isin(abc_cl8.obs.index))]
cbc_subc055_add = cbc_subc055[~cbc_subc055.isin(abc_str.obs.index)]

# >>> len(cbc_subc055_add)
# 865

# >>> abc_anno.loc[cbc_subc055_add, 'feature_matrix_label'].value_counts()
# feature_matrix_label
# WMB-10Xv3-PAL      684
# WMB-10Xv3-HY       161
# WMB-10Xv3-CTXsp     11
# WMB-10Xv3-TH         9
# Name: count, dtype: int64

### check before merge

# >>> abc_cl8
# AnnData object with n_obs × n_vars = 15662 × 32285
#     obs: 'cell_barcode', 'library_label', 'anatomical_division_label'
#     var: 'gene_symbol'
# >>> abc_str
# AnnData object with n_obs × n_vars = 285167 × 32285
#     obs: 'cell_barcode', 'library_label', 'anatomical_division_label'
#     var: 'gene_symbol'
#     uns: 'normalization', 'parent', 'parent_layer', 'parent_rows'

### merge
assert all(abc_str.var == abc_cl8.var)
abc_str = ad.concat([abc_str, abc_cl8[cbc_subc055_add, :]])
abc_str.var = abc_cl8.var

# >>> abc_str
# AnnData object with n_obs × n_vars = 286032 × 32285
#     obs: 'cell_barcode', 'library_label', 'anatomical_division_label'
#     var: 'gene_symbol'


### now join in the cell metadata
print("Joining ABC cell metadata...")
abc_str.obs = abc_str.obs.join(abc_anno)


# --------------------------------------------------
# Fix ABC ENSEMBL
# --------------------------------------------------

# doing it the same way it was done before

idx_keep = np.where(~abc_str.var['gene_symbol'].duplicated())[0]
# >>> abc_str.shape[1] - len(idx_keep)
# 40
# this is true:
if abc_str.raw is None:
    abc_str = ad.AnnData(
        X = abc_str.X[:, idx_keep],
        obs = abc_str.obs,
        var = abc_str.var.iloc[idx_keep, :],
        uns = abc_str.uns,
        obsm = abc_str.obsm,
        varm = abc_str.varm
    )
else:
    abc_str = ad.AnnData(
        X = abc_str.X[:, idx_keep],
        obs = abc_str.obs,
        var = abc_str.var.iloc[idx_keep, :],
        uns = abc_str.uns,
        obsm = abc_str.obsm,
        varm = abc_str.varm,
        raw = abc_str.raw[:, idx_keep]
    )

abc_str.var['ENSEMBL'] = abc_str.var.index
abc_str.var.set_index('gene_symbol', inplace=True)

# --------------------------------------------------
# Some ABC subsetting
# --------------------------------------------------

### Must at least get rid of the MSNs; and other Meis2+ cells

# >>> abc_str.obs['class'].value_counts()
# class
# 09 CNU-LGE GABA      78485
# 11 CNU-HYa GABA      40170
# 30 Astro-Epen        39192
# 31 OPC-Oligo         35165
# 10 LSX GABA          29175
# 13 CNU-HYa Glut      20296
# 33 Vascular          11979
# 34 Immune            11385
# 08 CNU-MGE GABA       5650
# 01 IT-ET Glut         5106
# 05 OB-IMN GABA        4493
# 06 CTX-CGE GABA       1462
# 03 OB-CR Glut          488
# 07 CTX-MGE GABA        485
# 02 NP-CT-L6b Glut      361
# 04 DG-IMN Glut         287
# 14 HY Glut             240
# 12 HY GABA             116
# 18 TH Glut              63
# 16 HY MM Glut           15
# 19 MB Glut              10
# 32 OEC                  10
# 29 CB Glut               5
# 20 MB GABA               5
# 17 MH-LH Glut            2
# 21 MB Dopa               1
# 15 HY Gnrh1 Glut         1
# Name: count, dtype: int64

# >>> abc_str.obs.neurotransmitter.value_counts()
# neurotransmitter
# GABA         156419
# Glut          26763
# Chol            326
# Glut-GABA       189
# Dopa             16
# Name: count, dtype: int64

# >>> abc_str.obs.region_of_interest_acronym.value_counts()
# region_of_interest_acronym
# sAMY     120764
# STRd      55626
# LSX       53819
# STRv      53573
# PAL         684
# HY          161
# CTXsp        11
# TH            9
# Name: count, dtype: int64

### Will start off with doing the full subsetting (for most cell types), and then just add back in all cells from subclass055

# will just drop a couple classes here: class9=LGE/MSNs, class5=Meis2+;
# note the Chat+ cells are have neurotransmitter given as Chol
idx_gaba     = np.where(abc_str.obs['neurotransmitter'].isin(['GABA', 'Chol']))[0]
idx_roi_keep = np.where(abc_str.obs['region_of_interest_acronym'].isin(['STRd', 'STRv']))[0]
idx_cls_drop = np.where(abc_str.obs['class'].isin(['09 CNU-LGE GABA', '05 OB-IMN GABA']))[0]
idx_subc55   = np.where(abc_str.obs['subclass'] == '055 STR Lhx8 Gaba')[0]

# >>> len(idx_gaba)
# 156745
# >>> len(idx_roi_keep)
# 109199
# >>> len(idx_cls_drop)
# 82978
# >>> len(idx_subc55)
# 1511

idx_keep = np.intersect1d(idx_gaba, idx_roi_keep)
# >>> len(idx_keep)
# 74597
idx_keep = np.setdiff1d(idx_keep, idx_cls_drop)
# >>> len(idx_keep)
# 2763
### not totally unexpected; vast majority of striatal interneurons are MSNs
idx_keep = np.union1d(idx_keep, idx_subc55)
# >>> len(idx_keep)
# 3770

abc_str  = abc_str[idx_keep, :].copy()


# --------------------------------------------------
# NOT USED: Subset ABC to striatum PLUS all other subclass 055 (supertypes 0236, 0237) MINUS any subclass 057
# --------------------------------------------------

# idx_sc55 = np.where(abc_str.obs['subclass'] == "055 STR Lhx8 Gaba")[0]
# idx_sc57 = np.where(abc_str.obs['subclass'] == "057 NDB-SI-MA-STRv Lhx8 Gaba")[0]
# idx_str  = np.where(abc_str.obs['anatomical_division_label'] == 'STR')[0]
# idx_keep = np.union1d(idx_sc55, idx_str)
# idx_keep = np.setdiff1d(idx_keep, idx_sc57)

# # >>> len(idx_sc55)
# # 1573
# # >>> len(idx_sc57)
# # 8429
# # >>> len(idx_str)
# # 4785
# # >>> len(idx_keep) # before remove subclass57
# # 5712
# # >>> len(idx_keep)
# # 4513

# abc_str = abc_str[idx_keep, :].copy()


# --------------------------------------------------
# Merge consistent 'clust_anno' for adult datasets (marm, abc, siletti)
# --------------------------------------------------

# Export cell barcodes + annotations (clust_anno) for the cells we annotated
abc_cl8_str_anno.obs = abc_cl8_str_anno.obs.join(
    clust_anno['abc'].rename({"seurat_clusters_anno" : "clust_anno"}, axis = 1).set_index("seurat_clusters"),
    on = "seurat_clusters"
)

# fill NaN for cells outside of this annotation set
abc_str.obs = abc_str.obs.join(abc_cl8_str_anno.obs.loc[:, ["clust_anno"]])

# >>> abc_cl8_str_anno.shape
# (3586, 32245)
# >>> abc_str.shape
# (3770, 32245)
# >>> sum(abc_cl8_str_anno.obs.index.isin(abc_str.obs.index))
# 2097


# --------------------------------------------------
# Add Tac2+/Tac2- for ABC Th class
# --------------------------------------------------

idx_gene_tac3    = np.where(abc_str.var.index == "Tac2")[0][0]
idx_cell_tac3    = abc_str.X.getcol(idx_gene_tac3).nonzero()[0] # index nz rows
idx_cell_th      = np.where(abc_str.obs['clust_anno'] == 'Th')[0]
idx_cell_th_tac3 = np.intersect1d(idx_cell_tac3, idx_cell_th)

abc_str.obs["clust_anno_tac3"] = abc_str.obs.loc[:, "clust_anno"]
idx_col_clustanno = abc_str.obs.columns.get_loc('clust_anno_tac3')
abc_str.obs.iloc[idx_cell_th,      idx_col_clustanno] = 'Th: Tac2-'
abc_str.obs.iloc[idx_cell_th_tac3, idx_col_clustanno] = 'Th: Tac2+'

# >>> abc_str.shape
# (3770, 32245)

# >>> abc_str.obs.clust_anno.value_counts()
# clust_anno
# Pthlh/Pvalb    757
# Th             634
# Sst/Npy        588
# Chat           107
# Lhx6/Prox1      11
# Name: count, dtype: int64
# 
# >>> abc_str.obs.clust_anno_tac3.value_counts()
# clust_anno_tac3
# Pthlh/Pvalb    757
# Sst/Npy        588
# Th: Tac2-      583
# Chat           107
# Th: Tac2+       51
# Lhx6/Prox1      11
# Name: count, dtype: int64


#### --- this is what this was in old analyses --- ####
# >>> abc_str.obs.clust_anno.value_counts()
# clust_anno
# Sst/Npy        1222
# Pthlh/Pvalb     969
# Th: Tac2-       625
# Lhx6/Prox1      434
# Chat            285
# Th: Tac2+        51
# Name: count, dtype: int64
#### --- --- ####


# --------------------------------------------------
# Also have a column where all types are clust_anno, except subclass055 cells are labelled by the ABC cluster
# --------------------------------------------------

abc_str.obs['clust_anno_splitsubc55'] = abc_str.obs['clust_anno'].astype(str)
abc_str.obs.loc[abc_str.obs.subclass == '055 STR Lhx8 Gaba', 'clust_anno_splitsubc55'] = abc_str.obs.loc[abc_str.obs.subclass == '055 STR Lhx8 Gaba', 'cluster']

# >>> abc_str.obs['clust_anno_splitsubc55'].value_counts()
# clust_anno_splitsubc55
# nan                     808
# Pthlh/Pvalb             748
# Sst/Npy                 565
# 0841 STR Lhx8 Gaba_1    527
# 0844 STR Lhx8 Gaba_2    479
# 0843 STR Lhx8 Gaba_2    226
# 0840 STR Lhx8 Gaba_1    118
# 0842 STR Lhx8 Gaba_1    113
# Chat                    107
# 0839 STR Lhx8 Gaba_1     48
# Th                       20
# Lhx6/Prox1               11
# Name: count, dtype: int64

### Also make a column where we just use the subclass (except split 55)
abc_str.obs['subclass_split55'] = abc_str.obs['subclass'].astype(str)
abc_str.obs.loc[abc_str.obs.subclass == '055 STR Lhx8 Gaba', 'subclass_split55'] = abc_str.obs.loc[abc_str.obs.subclass == '055 STR Lhx8 Gaba', 'cluster']


# --------------------------------------------------
# Write out
# --------------------------------------------------

abc_str.write_h5ad(Path(abc_data_dir, "abc_for_scvi.h5ad"))

abc_subset = abc_str[~abc_str.obs.clust_anno.isna(), ].copy()
abc_subset.write_h5ad(Path(abc_data_dir, "abc_for_scvi_anno_subset.h5ad"))
