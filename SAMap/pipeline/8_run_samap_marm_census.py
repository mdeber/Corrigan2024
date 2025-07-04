#!/usr/bin/env python

import pandas as pd
import numpy as np
import anndata as ad
from pathlib import Path
import os
import matplotlib.pyplot as plt
import seaborn as sns
import re

from samap.mapping import SAMAP
from samap.analysis import get_mapping_scores
from samalg import SAM
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

# --------------------------------------------------
# Paths
# --------------------------------------------------

script_dir = Path(__file__).parent.resolve()
dir_samap  = script_dir.parent.resolve()
pdir       = dir_samap.parent.resolve()
#
# [!] For ABC, want to have them split by cluster, like in the correlation analysis. 
#     And want to avoid NaN in annotations, so that's why using the subset_anno file
#
paths_h5ad = {
    "marm"      : Path(pdir, "marm_census/data/marm_for_scvi.h5ad"),
    "abc"       : Path(pdir, "abc/data/abc_for_scvi_anno_subset.h5ad"),
    "siletti"   : Path(pdir, "siletti/data/siletti_for_scvi.h5ad"),
    "marm_hmba" : Path(pdir, "marm_hmba/data/marm_hmba.h5ad") 
}
pdir_maps = Path(dir_samap, "maps_by_gene")
outdir    = Path(dir_samap, "results_marm_census")

if not outdir.is_dir():
    os.mkdir(outdir)

#
# Apparently I didn't have to remake the maps; could have loaded the full 
# original maps and passed a 'names' arg to SAMAP, see:
# https://github.com/atarashansky/SAMap/blob/main/samap/mapping.py
#


# --------------------------------------------------
# Run SAMAP
# --------------------------------------------------

# run from raw
sams = {
    'hg' : str(paths_h5ad['siletti']),
    'cj' : str(paths_h5ad['marm']),
    'mm' : str(paths_h5ad['abc'])
}
anno_keys = {
    'hg' : 'clust_anno',
    'cj' : 'clust_anno_lhx8', # where tac3 is split by lhx8+ vs lhx8-
    'mm' : 'clust_anno_splitsubc55'
}

if False:
    sm_blind = SAMAP(sams, str(pdir_maps) + "/")
    sm_blind.run(pairwise=True)

sm = SAMAP(sams, str(pdir_maps) + "/", keys = anno_keys)
# 19176 `hg` gene symbols match between the datasets and the BLAST graph.
# 14459 `cj` gene symbols match between the datasets and the BLAST graph. # vs 21900 for HMBA data
# 20812 `mm` gene symbols match between the datasets and the BLAST graph.

sm.run(pairwise=True)


# ===

# >>> sm.samap.adata
# AnnData object with n_obs × n_vars = 11420 × 110259
#     obs: 'hg_orig.ident', 'hg_ROIGroup', 'hg_ROIGroupCoarse', 'hg_ROIGroupFine', 'hg_roi', 'hg_organism_ontology_term_id', 'hg_disease_ontology_term_id', 'hg_self_reported_ethnicity_ontology_term_id', 'hg_assay_ontology_term_id', 'hg_sex_ontology_term_id', 'hg_development_stage_ontology_term_id', 'hg_donor_id', 'hg_dissection', 'hg_sample_id', 'hg_supercluster_term', 'hg_cell_type_ontology_term_id', 'hg_tissue_ontology_term_id', 'hg_clust_anno', 'cj_orig.ident', 'cj_class_label', 'cj_class_name', 'cj_subclass_label', 'cj_subclass_name', 'cj_supertype_label', 'cj_supertype_name', 'cj_cluster_label', 'cj_cluster_name', 'cj_CLUSTER.SUBCLUSTER', 'cj_organism_ontology_term_id', 'cj_self_reported_ethnicity_ontology_term_id', 'cj_disease_ontology_term_id', 'cj_cell_type_ontology_term_id', 'cj_cell_type_ontology_term_name', 'cj_assay_ontology_term_id', 'cj_suspension_type', 'cj_seq_pool', 'cj_donor_id', 'cj_structure', 'cj_region', 'cj_anatomical_name', 'cj_sex', 'cj_age', 'cj_tissue_ontology_term_id', 'cj_tissue_ontology_term_name', 'cj_development_stage_ontology_term_id', 'cj_sex_ontology_term_id', 'cj_clust_anno', 'mm_cell_barcode', 'mm_library_label', 'mm_anatomical_division_label', 'mm_feature_matrix_label', 'mm_library_method', 'mm_region_of_interest_acronym', 'mm_donor_label', 'mm_donor_genotype', 'mm_donor_sex', 'mm_neurotransmitter', 'mm_class', 'mm_subclass', 'mm_supertype', 'mm_cluster', 'mm_clust_anno', 'mm_clust_anno_tac3', 'mm_clust_anno_splitsubc55', 'mm_subclass_split55', 'batch', 'species'
#     uns: 'neighbors', 'gnnm_corr', 'mapping_K', 'umap', 'homology_gene_names_dict'
#     obsm: 'X_umap'
#     layers: 'X_disp'
#     obsp: 'connectivities'
#     varp: 'homology_graph_reweighted', 'homology_graph'



# --------------------------------------------------
# Calculate cell type mapping scores
# --------------------------------------------------

# n_top uses all cells from clusters; setting n_top can use just the top n aligning cells
_, mt = get_mapping_scores(sm, anno_keys, n_top=0)

# >>> mt
#                            hg_CCK  hg_CCK/VIP   hg_CHAT  hg_PTHLH/PVALB  hg_SST/NPY   hg_TAC3    cj_CCK  ...  mm_0843 STR Lhx8 Gaba_2  mm_0844 STR Lhx8 Gaba_2   mm_Chat  mm_Lhx6/Prox1  mm_Pthlh/Pvalb  mm_Sst/Npy     mm_Th
# hg_CCK                   0.000000    0.000000  0.000000        0.000000    0.000000  0.000000  0.988750  ...                 0.000000                 0.000125  0.000000       0.000000        0.003209    0.003635  0.000000
# hg_CCK/VIP               0.000000    0.000000  0.000000        0.000000    0.000000  0.000000  0.022279  ...                 0.000000                 0.002267  0.000000       0.000000        0.000367    0.021128  0.000000
# hg_CHAT                  0.000000    0.000000  0.000000        0.000000    0.000000  0.000000  0.000000  ...                 0.000000                 0.000000  0.994023       0.000000        0.000000    0.000000  0.000000
# hg_PTHLH/PVALB           0.000000    0.000000  0.000000        0.000000    0.000000  0.000000  0.002675  ...                 0.000000                 0.529931  0.000000       0.083401        0.974391    0.000637  0.306467
# hg_SST/NPY               0.000000    0.000000  0.000000        0.000000    0.000000  0.000000  0.000000  ...                 0.000000                 0.008333  0.000000       0.000000        0.000000    0.983615  0.000000
# hg_TAC3                  0.000000    0.000000  0.000000        0.000000    0.000000  0.000000  0.000483  ...                 0.004501                 0.210072  0.000000       0.000000        0.011105    0.001045  0.078980
# cj_CCK                   0.988750    0.022279  0.000000        0.002675    0.000000  0.000483  0.000000  ...                 0.000755                 0.025418  0.000000       0.024559        0.118170    0.000012  0.000000
# cj_CCK/VIP               0.007379    0.968854  0.000000        0.004445    0.001651  0.020492  0.000000  ...                 0.496354                 0.332662  0.000000       0.000008        0.012708    0.005868  0.010247
# cj_CHAT                  0.000000    0.000000  0.997585        0.000000    0.001550  0.002022  0.000000  ...                 0.000000                 0.000000  0.979505       0.000000        0.000000    0.001505  0.000000
# cj_LAMP5/LHX6            0.008267    0.000000  0.000000        0.342712    0.000000  0.003539  0.000000  ...                 0.000000                 0.210926  0.000000       0.067018        0.167965    0.000000  0.091446
# cj_PTHLH/PVALB           0.000325    0.000000  0.000000        0.938887    0.001948  0.021801  0.000000  ...                 0.001119                 0.192712  0.000000       0.273759        0.920364    0.004463  0.221733
# cj_SST/NPY               0.000000    0.000000  0.000660        0.011429    0.997072  0.036004  0.000000  ...                 0.001236                 0.044517  0.000186       0.000000        0.001274    0.978130  0.000718
# cj_TAC3                  0.000000    0.008601  0.000000        0.013406    0.000496  0.954184  0.000000  ...                 0.000000                 0.030924  0.000000       0.000000        0.018500    0.000386  0.092754
# cj_TAC3/LHX8             0.000000    0.002911  0.000000        0.031981    0.000000  0.820682  0.000000  ...                 0.012363                 0.030093  0.000000       0.000000        0.012273    0.004237  0.001761
# mm_0839 STR Lhx8 Gaba_1  0.000000    0.000000  0.000000        0.168138    0.000000  0.272708  0.019182  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_0840 STR Lhx8 Gaba_1  0.000000    0.000000  0.000000        0.000000    0.000000  0.865738  0.000000  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_0841 STR Lhx8 Gaba_1  0.000000    0.000000  0.000000        0.000000    0.000000  0.875593  0.000330  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_0842 STR Lhx8 Gaba_1  0.000000    0.000000  0.000000        0.000000    0.000000  0.960839  0.000024  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_0843 STR Lhx8 Gaba_2  0.000000    0.000000  0.000000        0.000000    0.000000  0.004501  0.000755  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_0844 STR Lhx8 Gaba_2  0.000125    0.002267  0.000000        0.529931    0.008333  0.210072  0.025418  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_Chat                  0.000000    0.000000  0.994023        0.000000    0.000000  0.000000  0.000000  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_Lhx6/Prox1            0.000000    0.000000  0.000000        0.083401    0.000000  0.000000  0.024559  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_Pthlh/Pvalb           0.003209    0.000367  0.000000        0.974391    0.000000  0.011105  0.118170  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_Sst/Npy               0.003635    0.021128  0.000000        0.000637    0.983615  0.001045  0.000012  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000
# mm_Th                    0.000000    0.000000  0.000000        0.306467    0.000000  0.078980  0.000000  ...                 0.000000                 0.000000  0.000000       0.000000        0.000000    0.000000  0.000000


mt.to_csv(Path(outdir, "mapping_table_clust_anno_splitsubc55.csv"))


### Heatmap [no row-ordering, diagonals set to 0]

plt.figure(figsize=(8, 6))
ax = sns.heatmap(mt, cmap='viridis', vmin=0, vmax=1)
ax.set_frame_on(True)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.8)
    spine.set_color('black')

cbar = ax.collections[0].colorbar
cbar.outline.set_linewidth(0.8)
plt.title('SAMap')
# plt.xlabel('')
# plt.ylabel('')
plt.tight_layout()
plt.savefig(Path(outdir, "mapping_table_clust_anno_splitsubc55_heatmap.pdf"))
plt.clf()
del ax, cbar

# --------------------------------------------------
# hclust mapping scores
# --------------------------------------------------



# convert similarity matrix to condensed format distance matrix
# -> it looks like the diagonal is zeroed out by default (which is weird)...
distance_matrix = 1 - mt
for i in range(distance_matrix.shape[0]):
    distance_matrix.iloc[i,i] = 0
condensed_dist = ssd.squareform(distance_matrix.values)

# hclust (make linkage matrix) and get row order
linkage = sch.linkage(condensed_dist, method='average')
dendro = sch.dendrogram(linkage, no_plot=True)
ordered_indices = dendro['leaves']
ordered_labels = mt.index[ordered_indices]

# make reordered df
mt_reorder = mt.loc[ordered_labels, ordered_labels]

### Heatmap

plt.figure(figsize=(8, 6))
ax = sns.heatmap(mt_reorder, cmap='viridis', vmin=0, vmax=1)
ax.set_frame_on(True)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.8)
    spine.set_color('black')

cbar = ax.collections[0].colorbar
cbar.outline.set_linewidth(0.8)
plt.title('SAMap')
# plt.xlabel('')
# plt.ylabel('')
plt.tight_layout()
plt.savefig(Path(outdir, "mapping_table_clust_anno_splitsubc55_heatmap_hclust.pdf"))
plt.clf()
del ax, cbar


### Redo this without zeroing out the diagonal...
mt_diag = mt
for i in range(mt_diag.shape[0]):
    mt_diag.iloc[i,i] = 1
mt_diag_reord = mt_diag.loc[ordered_labels, ordered_labels]


plt.figure(figsize=(8, 6))
ax = sns.heatmap(mt_diag_reord, cmap='viridis', vmin=0, vmax=1)
ax.set_frame_on(True)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.8)
    spine.set_color('black')

cbar = ax.collections[0].colorbar
cbar.outline.set_linewidth(0.8)
plt.title('SAMap')
# plt.xlabel('')
# plt.ylabel('')
plt.tight_layout()
plt.savefig(Path(outdir, "mapping_table_clust_anno_splitsubc55_heatmap_hclust_fixdiag.pdf"))
plt.clf()
del ax, cbar


