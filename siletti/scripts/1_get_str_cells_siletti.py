"""
Siletti atlas datasets for striatal interneuron reclustering:
    Basal Nuclei (BN) - Body of the Caudate - CaB
    Basal Nuclei (BN) - Putamen - Pu
    Basal Nuclei (BN) - Nucleus Accumbens - NAC

In these datasets, select cells by supercluster_term:
    CGE interneuron
    MGE interneuron
    Splatter

Note:
Due to dissection artifacts that seem to be present in these datasets, 
I expect most of the striatal interneurons are in the "Splatter" type and that 
many of those classified as "CGE interneuron" and "MGE interneuron" type 
actually come from neighboring brain structures and not striatum

After selecting the cells of interest, this script writes out MatrixMarket
files for the data and csv files for the metadata (for subsequent import
into R for reprocessing).
"""

import anndata as ad
import pandas as pd
import numpy as np
from scipy.io import mmwrite
from pathlib import Path

# Various object summaries at bottom of script
siletti_script_dir = Path(__file__).parent.resolve()
siletti_data_dir = Path(siletti_script_dir.parent, "data")
path_h5ad = Path(siletti_data_dir, "Neurons.h5ad")
adata = ad.read_h5ad(path_h5ad, backed = 'r')

idx_roi = np.where(adata.obs.roi.values.isin(['Human CaB', 'Human Pu', 'Human NAC']))[0]
# n=88044 [not the same number as the loom (n=97226), but this is just neurons]
# (-> implies there are 9182 non-neurons in those tissues)

## [not used for anything]
# idx_supercl = np.where(adata.obs.supercluster_term.values.isin(['CGE interneuron', 'MGE interneuron', 'Splatter']))[0]
## n=741938

# >>> len(np.intersect1d(idx_roi, idx_supercl))
# 7665 # -> not many cells in common

# >>> len(np.union1d(idx_roi, idx_supercl))
# 822317

# Write out genes
adata.var.to_csv(Path(siletti_data_dir, 'Neurons_var.csv'))

# Write out files for roi-selected cells
# (note: adata[idx_roi, ].X uses less memory than adata.X[idx_roi, ])
adata.obs.iloc[idx_roi, :].to_csv(Path(siletti_data_dir, 'Neurons_select_roi_obs.csv'))
mmwrite(Path(siletti_data_dir, 'Neurons_select_roi_X.mtx'), adata[idx_roi, ].X)

## [not used] Write out files for supercluster-term-selected cells
# adata.obs.iloc[idx_supercl, :].to_csv(Path(siletti_data_dir, 'Neurons_select_supercluster_obs.csv'))
# mmwrite(Path(siletti_data_dir, 'Neurons_select_supercluster_X.mtx'), adata[idx_supercl, ].X)

# ----------

# >>> adata
# AnnData object with n_obs x n_vars = 2480956 x 59480 backed at '/Users/md6347/Library/CloudStorage/Dropbox/krienen_lab/240314_ReclusterSilettiStriatum/Neurons.h5ad'
#     obs: 'ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi', 'organism_ontology_term_id', etc. etc. [see next line]
#     var: 'Biotype', 'Chromosome', 'End', 'Gene', 'Start'
#     uns: 'batch_condition', 'schema_version', 'title'
#     obsm: 'X_UMAP', 'X_tSNE'

# >>> adata.obs.columns
# Index(['ROIGroup', 'ROIGroupCoarse', 'ROIGroupFine', 'roi',
#        'organism_ontology_term_id', 'disease_ontology_term_id',
#        'self_reported_ethnicity_ontology_term_id', 'assay_ontology_term_id',
#        'sex_ontology_term_id', 'development_stage_ontology_term_id',
#        'donor_id', 'dissection', 'cell_cycle_score', 'sample_id', 'cluster_id',
#        'subcluster_id', 'supercluster_term', 'cell_type_ontology_term_id',
#        'tissue_ontology_term_id'],
#       dtype='object')

# >>> adata.obs.roi.values
# ['Human MoAN', 'Human MoSR', 'Human MoAN', 'Human PnAN', 'Human MoAN', ..., 'Human A1C', 'Human SEP', 'Human GPe', 'Human PnRF', 'Human ITG']
# Length: 2480956
# Categories (108, object): ['Human A1C', 'Human A5-A7', 'Human A13', 'Human A14', ..., 'Human V2',
#                            'Human VA', 'Human VLN', 'Human VPL']
#
# 108 categories, when in loom file there were 105 'Roi' terms and 107 'Tissue' terms...

# >>> list(adata.obs.roi.values.categories)
# ['Human A1C', 'Human A5-A7', 'Human A13', 'Human A14', 'Human A19', 'Human A23', 
#  'Human A24', 'Human A25', 'Human A29-A30', 'Human A32', 'Human A35-36', 
#  'Human A35-A36', 'Human A35r', 'Human A38', 'Human A40', 'Human A43', 
#  'Human A44-A45', 'Human A46', 'Human ANC', 'Human AON', 'Human BL', 'Human BM', 
#  'Human BNST', 'Human CA1-3', 'Human CA1C-CA3C', 'Human CA1R-CA2R', 
#  'Human CA1R-CA2R-CA3R', 'Human CA1U', 'Human CA1U-CA2U-CA3U', 'Human CA2U-CA3U',
#  'Human CA3R', 'Human CA4C-DGC', 'Human CBL', 'Human CBV', 'Human CEN', 'Human CM',
#  'Human CM-Pf', 'Human CMN', 'Human CaB', 'Human CbDN', 'Human Cla', 'Human CoA', 
#  'Human DGR-CA4Rpy', 'Human DGU-CA4Upy', 'Human DTg', 'Human ETH', 'Human FI', 
#  'Human GPe', 'Human GPi', 'Human Gpe', 'Human HTHma', 'Human HTHma-HTHtub', 
#  'Human HTHpo', 'Human HTHpo-HTHso', 'Human HTHso', 'Human HTHso-HTHtub', 
#  'Human HTHtub', 'Human IC', 'Human IO', 'Human ITG', 'Human Idg', 'Human Ig', 
#  'Human LEC', 'Human LG', 'Human LP', 'Human LP-VPL', 'Human La', 'Human M1C', 
#  'Human MD', 'Human MD-Re', 'Human MEC', 'Human MG', 'Human MN', 'Human MTG', 
#  'Human MoAN', 'Human MoRF-MoEN', 'Human MoSR', 'Human NAC', 'Human PAG', 
#  'Human PAG-DR', 'Human PB', 'Human PN', 'Human PTR', 'Human Pir', 'Human PnAN', 
#  'Human PnEN', 'Human PnRF', 'Human Pro', 'Human Pu', 'Human Pul', 'Human RN', 
#  'Human S1C', 'Human SC', 'Human SEP', 'Human SI', 'Human SN', 'Human SN-RN', 
#  'Human STG', 'Human STH', 'Human SpC', 'Human Sub', 'Human TF', 'Human TH-TL', 
#  'Human V1C', 'Human V2', 'Human VA', 'Human VLN', 'Human VPL']
#
# Includes our roi terms of interest:
# 'Human CaB'
# 'Human Pu'
# 'Human NAC'

# >>> adata.obs.supercluster_term.values.categories
# Index(['Amygdala excitatory', 'CGE interneuron', 'Cerebellar inhibitory',
#        'Deep-layer corticothalamic and 6b', 'Deep-layer intratelencephalic',
#        'Deep-layer near-projecting', 'Eccentric medium spiny neuron',
#        'Hippocampal CA1-3', 'Hippocampal CA4', 'Hippocampal dentate gyrus',
#        'LAMP5-LHX6 and Chandelier', 'Lower rhombic lip', 'MGE interneuron',
#        'Mammillary body', 'Medium spiny neuron', 'Midbrain-derived inhibitory',
#        'Miscellaneous', 'Splatter', 'Thalamic excitatory', 'Upper rhombic lip',
#        'Upper-layer intratelencephalic'],
#       dtype='object')
#
# Includes our supercluster terms of interest:
# 'CGE interneuron'
# 'MGE interneuron'
# 'Splatter'
