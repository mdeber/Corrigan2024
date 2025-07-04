#!/usr/bin/env python

# 
# Integrating the following datasets. In all cases, see their data directories for the
# actual data selection and subsetting, but summarized here:
# 
# - marmoset HMBA
#   - beginning with all unlabelled cells and subsetting based on cross-species integration
# - ABC (mouse) atlas
#   - beginning with all cells from striatal dissections PLUS all cells from subclass 055,
#     but then subsetting that by cross-species integration
# - Siletti (human) atlas
#   - using our existing analysis set with some additional filtering e.g., no "mixed" or hypothalamic cells
# - marmoset census
#   - using our existing analysis set with some additional filtering, e.g., no "mixed" cells
#
# One complication here is that we now have two poorly controlled inputs (from ABC and marmoset HMBA),
# which may complicate curation of cells that are truly striatal interneurons, so careful attention
# will have to be given to any clusters that only contain cells from those datasets (and which to keep).
#

import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import scvi
from pathlib import Path
from scipy.sparse import csr_matrix
import scipy.stats
import os
import matplotlib.pyplot as plt
import seaborn as sns
import re

# --------------------------------------------------
# Paths
# --------------------------------------------------

pdir       = Path(__file__).parent.resolve()
path_ortho = Path(pdir, "ortholog_tables/ens_112_ortho_human_ferret_pig_mouse_marmoset.csv")
outdir     = Path(pdir, "scvi_integrations")

if not Path.is_dir(outdir):
    os.mkdir(outdir)

paths_h5ad = {
    "marm"      : Path(pdir, "marm_census/data/marm_for_scvi.h5ad"),
    "abc"       : Path(pdir, "abc/data/abc_for_scvi.h5ad"),
    "siletti"   : Path(pdir, "siletti/data/siletti_for_scvi.h5ad"),
    "marm_hmba" : Path(pdir, "marm_hmba/data/rna_merged_qc_pass_raw.h5ad") 
}

### HMBA is prerelease from: 240905_reprocess_and_recluster/data/rna_merged_qc_pass_raw.h5ad
### --> see that data directory for cell barcodes to recreate the starting (and final) sets

# --------------------------------------------------
# Imports
# --------------------------------------------------

print("Importing data...")
all_h5ad = {
    k: ad.read_h5ad(v) for k, v in paths_h5ad.items()
}

df_ortho = pd.read_csv(path_ortho)

# human is the reference:
#
# >>> df_ortho.columns
# Index(['Gene.stable.ID', 'Gene.stable.ID.version', 'Gene.name',
#        'Ferret.gene.stable.ID', 'Ferret.gene.name',
#        'X.id..target.Ferret.gene.identical.to.query.gene',
#        'X.id..query.gene.identical.to.target.Ferret.gene',
#        'Ferret.orthology.confidence..0.low..1.high.', 'Pig.gene.stable.ID',
#        'Pig.gene.name', 'X.id..target.Pig.gene.identical.to.query.gene',
#        'X.id..query.gene.identical.to.target.Pig.gene',
#        'Pig.orthology.confidence..0.low..1.high.', 'Mouse.gene.stable.ID',
#        'Mouse.gene.name', 'X.id..target.Mouse.gene.identical.to.query.gene',
#        'X.id..query.gene.identical.to.target.Mouse.gene',
#        'Mouse.orthology.confidence..0.low..1.high.',
#        'White.tufted.ear.marmoset.gene.stable.ID',
#        'White.tufted.ear.marmoset.gene.name',
#        'X.id..target.White.tufted.ear.marmoset.gene.identical.to.query.gene',
#        'X.id..query.gene.identical.to.target.White.tufted.ear.marmoset.gene',
#        'White.tufted.ear.marmoset.orthology.confidence..0.low..1.high.'],
#       dtype='object')


# --------------------------------------------------
# Get single ortholog table
# --------------------------------------------------

# mapping between dataset names and their gene name columns in the ortholog table
dict_dat_colname = {
    'marm'      : 'White.tufted.ear.marmoset.gene.name',
    'abc'       : 'Mouse.gene.name',
    'siletti'   : 'Gene.name',
    'marm_hmba' : 'White.tufted.ear.marmoset.gene.name'
}

print("Number of dataset genes found in ortholog table:")
for dname, cname in dict_dat_colname.items():
    ng_dat = all_h5ad[dname].shape[1]
    ng_ortho = sum(all_h5ad[dname].var.index.isin(df_ortho[cname]))
    print(f"{dname}: {ng_ortho} / {ng_dat}")
    del ng_dat, ng_ortho

# Number of dataset genes found in ortholog table:
# marm: 11857 / 19762
# abc: 17383 / 32245
# siletti: 37297 / 58252
# marm_hmba: 14736 / 35787

# (!) Fix marmoset TAC3. Missing in this version of ensembl (although even mouse correctly has it as "Tac2")
df_ortho.loc[df_ortho['Gene.name'] == "TAC3", 'White.tufted.ear.marmoset.gene.name'] = "TAC3"

# Get single ortholog table
df_ortho_3way = df_ortho
for dname, cname in dict_dat_colname.items():
    df_ortho_3way = df_ortho_3way[df_ortho_3way[cname].isin(all_h5ad[dname].var.index)]
    df_ortho_3way = df_ortho_3way.drop_duplicates(subset = cname)

# columns to keep, but without duplicates and keeping in this order (so not using 'set')
cols_keep = list(dict.fromkeys(
    ['Gene.name', 'Gene.stable.ID', 'Gene.stable.ID.version'] + list(dict_dat_colname.values())
))
# >>> cols_keep
# ['Gene.name', 'Gene.stable.ID', 'Gene.stable.ID.version', 'White.tufted.ear.marmoset.gene.name', 'Mouse.gene.name']
df_ortho_3way = df_ortho_3way.loc[:, cols_keep]

print(f"Genes in a single 3-way ortholog table (well, 3 species but 4 datasets): {df_ortho_3way.shape[0]}")
# Genes in a single 3-way ortholog table: 11280


# --------------------------------------------------
# Subset & rename genes in anndata objects
# --------------------------------------------------

# AnnData deliberately prevents subsetting genes in raw, and because
# scVI only uses raw, I will just copy raw into X and delete raw
#
# in hindsight not the best idea b/c i end up duplicating the anndata to normalize
for dname, cname in dict_dat_colname.items():
    # re-assign raw, then delete, and subset genes
    if not isinstance(all_h5ad[dname].X, csr_matrix):
        all_h5ad[dname].X = csr_matrix(all_h5ad[dname].X)
    # reassigning raw as X (and using indexing in case X has a subset of genes in raw)
    if not all_h5ad[dname].raw is None:
        all_h5ad[dname].X = all_h5ad[dname].raw[:, all_h5ad[dname].var.index].X
        del all_h5ad[dname].raw
    all_h5ad[dname] = all_h5ad[dname][:, df_ortho_3way[cname]].copy()
    # join in gene ortholog table (to keep gene name mappings, human ENSEMBL ID, etc.)
    all_h5ad[dname].var = all_h5ad[dname].var.join(df_ortho_3way.set_index(cname), how = "left")
    # for any species not currently using the human gene names:
    if "Gene.name" in all_h5ad[dname].var.columns:
        all_h5ad[dname].var[cname] = all_h5ad[dname].var.index
        all_h5ad[dname].var = all_h5ad[dname].var.set_index("Gene.name")

print("Finished subsetting.")


# --------------------------------------------------
# Add somewhat consistent 'donor_coding' and 'species' columns for integration
# --------------------------------------------------

all_h5ad['marm_hmba'].obs['donor'] = [re.sub("^.*_HMBA_", "", x) for x in all_h5ad['marm_hmba'].obs.sample_name]
all_h5ad['marm_hmba'].obs['donor'] = [re.sub("_.*", "", x)       for x in all_h5ad['marm_hmba'].obs['donor']]

dict_dat_donorcol = {
    'marm'      : 'donor_id',
    'abc'       : 'donor_label',
    'siletti'   : 'donor_id',
    'marm_hmba' : 'donor'
}

for dname, cname in dict_dat_donorcol.items():
    if cname == '':
        all_h5ad[dname].obs['donor_coding'] = 'single_donor'
    else:
        all_h5ad[dname].obs['donor_coding'] = all_h5ad[dname].obs[cname]

# not used for integration (instead, integrating across dataset/experiment)
dict_dat_speciescol = {
    'marm'      : 'marmoset',
    'abc'       : 'mouse',
    'siletti'   : 'human',
    'marm_hmba' : 'marmoset'
}

for dname, sname in dict_dat_speciescol.items():
    all_h5ad[dname].obs['species'] = sname


# --------------------------------------------------
# Merge for integration
# --------------------------------------------------

print("Merging data...")

# I think ad.concat generates layer 'counts', which is all nan
adata = ad.concat(all_h5ad, join = "outer", label = "experiment")
del all_h5ad


# --------------------------------------------------
# Remove problematic columns
# --------------------------------------------------

# can't write out the h5ad with this NaN-containing column 'is_primary_data':
adata.obs.drop(['is_primary_data'], axis = 1, inplace = True)


# --------------------------------------------------
# Drop donors with <30 cells
# --------------------------------------------------

# Note that there are some donors with very few counts...; we may want to drop donors with maybe <10 or <50 cells...
donor_counts = adata.obs['donor_coding'].value_counts()

#
# >>> len(donor_counts.index)
# 59
# >>> sum(donor_counts < 50)
# 33
# >>> sum(donor_counts < 30)
# 31

###
### will just inspect the 2 donors that are added if we set threshold at 50 cells per donor
###

# chk_donor = list(donor_counts.index[(donor_counts < 50) & (donor_counts > 30)])
# >>> adata.obs.loc[adata.obs.donor_coding.isin(chk_donor), 'region_of_interest_acronym'].value_counts()
# region_of_interest_acronym
# HY       58
# sAMY     15

### well that's a puzzle, neither fit my search criteria, unless...
# 
# >>> adata.obs.loc[adata.obs.donor_coding == chk_donor[0], 'region_of_interest_acronym'].value_counts()
# region_of_interest_acronym
# HY       40
# >>> adata.obs.loc[adata.obs.donor_coding == chk_donor[1], 'region_of_interest_acronym'].value_counts()
# region_of_interest_acronym
# HY       18
# sAMY     15
# 
### ...yea the same donor had multiple dissections apparently. will keep

adata = adata[adata.obs['donor_coding'].isin(donor_counts.index[donor_counts >= 30]), :].copy()


# --------------------------------------------------
# scVI model
# --------------------------------------------------

# Finding HVGs and modeling integration by experiment.
# Saving counts for all 9428 homologs in raw (and only HVGs in X for integration) 
adata.raw = adata
sc.pp.highly_variable_genes(
    adata, 
    flavor      = "seurat_v3", 
    n_top_genes = 4000, 
    batch_key   = "experiment", 
    subset      = True
)

print("Saving merged h5ad with HVGs selected...")
adata.write_h5ad(Path(outdir, "adata.h5ad"))

## re-import
# adata = ad.read_h5ad(Path(outdir, "adata.h5ad"))

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key = "experiment",
    categorical_covariate_keys = ['donor_coding']
)

# default gene_likelihood='zinb'
scvi_model_initial_int = scvi.model.SCVI(
    adata, 
    n_hidden   = 256,
    n_latent   = 64,
    n_layers   = 3,
    dispersion = "gene-batch"
)

print("Train scVI (initial_int)...")
scvi_model_initial_int.train(accelerator = "gpu", max_epochs = 250)

print("Save model (initial_int)...")
scvi_model_initial_int.save(Path(outdir, "scvi_model_initial_int"))

## re-import
# scvi_model_initial_int = scvi.model.SCVI.load(Path(outdir, "scvi_model_initial_int"), adata) 


# --------------------------------------------------
# Get, save scVI embeddings
# --------------------------------------------------

adata.obsm["X_scVI"] = scvi_model_initial_int.get_latent_representation()
sc.pp.neighbors(adata, use_rep = "X_scVI")

# setting min_dist=0.25 and keeping spread=1; will see, we generally will 
# prefer more clumping because our interest is primarily label transfer
adata.obsm["X_umap_integrated"] = sc.tl.umap(adata, min_dist = 0.25, copy = True).obsm["X_umap"]

np.savetxt(
    Path(outdir, "scvi_model_initial_int_X_scVI.csv"), 
    adata.obsm["X_scVI"], 
    delimiter = ","
)

np.savetxt(
    Path(outdir, "scvi_model_initial_int_X_umap_integrated.csv"), 
    adata.obsm["X_umap_integrated"], 
    delimiter = ","
)

## re-import
# adata.obsm["X_scVI"] = np.genfromtxt(Path(outdir, "scvi_model_initial_int_X_scVI.csv"), delimiter = ',')
# adata.obsm["X_umap_integrated"] = np.genfromtxt(Path(outdir, "scvi_model_initial_int_X_umap_integrated.csv"), delimiter = ',')


# --------------------------------------------------
# Add some metadata columns
# --------------------------------------------------

# experiment + clust_anno
adata.obs['experiment_clust_anno'] = adata.obs['experiment'].astype(str) + '_' + adata.obs['clust_anno'].astype(str)

# ABC atlas clusters
adata.obs['experiment_cluster'] = adata.obs['experiment'].astype(str) + '_' + adata.obs['cluster'].astype(str)

### add column for experiment + dissection_region...
# 'structure' [that's marm, all str], 
# 'region' [marm], 
# 'anatomical_name' [marm], 
# 'anatomical_division_label' [abc, all str], 
# 'ROIGroup' [siletti, all cerebral nuclei], 
# 'ROIGroupCoarse' [ditto], 
# 'ROIGroupFine' [siletti, all basal nuclei], 
# 'roi' [siletti, CaB, Pu, and NAC], 
# 'dissection' [same but not contracted],
# 
# so it's 'roi' for siletti, and 'region_of_interest_acronym' for abc
# 
adata.obs['experiment_roi'] = adata.obs['experiment'].astype(str) + '_' + adata.obs['roi'].astype(str)
adata.obs.loc[adata.obs.experiment == 'abc', 'experiment_roi'] = 'abc_' + adata.obs.loc[adata.obs.experiment == 'abc', 'region_of_interest_acronym'].astype(str)


# --------------------------------------------------
# Export h5ad
# --------------------------------------------------

# adata.write_h5ad(Path(outdir, "adata_model_initial_int.h5ad"))

## re-import
# adata = ad.read_h5ad(Path(outdir, "adata_model_initial_int.h5ad"))


# --------------------------------------------------
# Cluster entropy: calculation
# --------------------------------------------------

### Changing things up a little here, will have to normalize adata at this point

adata_norm = ad.AnnData(
    X = adata.raw.X.copy(),
    obs = adata.obs,
    var = adata.raw.var,
    uns = adata.uns,
    obsm = adata.obsm,
    # raw = adata.raw
)
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)

# get an integrated clustering and cluster entropy ->
# this is following what Nelson Johansen at the Allen Institute has been doing for the ongoing Human and Mammalian Brain Atlas
sc.pp.neighbors(adata_norm, n_neighbors=10, use_rep="X_scVI")
sc.tl.leiden(adata_norm, flavor="igraph", n_iterations=2)

for anno in ["experiment"]:
    adata_norm.obs[anno + "_entropy"] = -1 ## Initialize with a value outside range of entropy.
    for cluster in np.unique(adata_norm.obs.leiden):
        adata_norm.obs.loc[adata_norm.obs.leiden == cluster, anno + "_entropy"] = scipy.stats.entropy(
            adata_norm.obs.loc[adata_norm.obs.leiden == cluster, anno].value_counts() / sum(adata_norm.obs.leiden == cluster)
        )


# --------------------------------------------------
# Cluster entropy: Set and check cutoff (start at 0.1)
# --------------------------------------------------

## Want to flag clusters from essentially one study (basically only marm_hmba in this case) using threshold of study entropy
adata_norm.obs["toRemove"] = adata_norm.obs["experiment_entropy"] <= 0.1

# for plotting
adata_norm.obs["experiment_toRemove"] = adata_norm.obs["experiment"].astype(str) + "_" + adata_norm.obs["toRemove"].astype(str)


# --------------------------------------------------
# Cluster entropy: Checks when cutoff=0.1
# --------------------------------------------------

# another 2.5k cells to remove when including all subc055

# >>> adata_norm.obs.toRemove.value_counts()
# toRemove
# True     270508
# False     25174
# Name: count, dtype: int64

# >>> adata_norm.obs.loc[:, ['experiment', 'toRemove']].value_counts()
# experiment  toRemove
# marm_hmba   True        267945
#             False        14861
# siletti     False         3618
# abc         False         3454
# marm        False         3241
# siletti     True          2424
# abc         True            99
# marm        True            40
# Name: count, dtype: int64


### it's cells from siletti; w/o all of subc055, only 8 siletti cells failed to integrate;
### but almost the same marm_hmba cells were flagged

# >>> adata_norm.obs.loc[:, ['experiment_clust_anno', 'toRemove']].value_counts()
# experiment_clust_anno  toRemove
# marm_hmba_nan          True        267945
#                        False        14861
# abc_nan                False         3454
# siletti_PTHLH/PVALB    False         2560
# siletti_TAC3           True          2413
# marm_PTHLH/PVALB       False         1225
# marm_TAC3              False         1097
# marm_SST/NPY           False          515
# siletti_SST/NPY        False          396
# marm_CHAT              False          269
# siletti_CCK/VIP        False          228
# siletti_CHAT           False          194
# siletti_CCK            False          190
# marm_CCK               False          135
# abc_nan                True            99
# siletti_TAC3           False           50
# marm_TAC3              True            14
# marm_CHAT              True            10
# marm_PTHLH/PVALB       True             9
# marm_SST/NPY           True             7
# siletti_PTHLH/PVALB    True             4
# siletti_SST/NPY        True             3
# siletti_CCK/VIP        True             2
# siletti_CHAT           True             2
# Name: count, dtype: int64

### and it's all the siletti TAC3; very surprising

# --------------------------------------------------
# Cluster entropy: Check HMBA cells being kept/dropped for cutoff=0.1
#
# -> went back to check using the HMBA taxonomy
# --------------------------------------------------

# hmba_bc_kept = adata_norm.obs.loc[(~adata_norm.obs.toRemove) & (adata_norm.obs.experiment == "marm_hmba"), :].index

# # import the subsetted & annotated version of the data
# # path_cell_meta = "{AIBS_S3}/241107_nutmeg_bg_RNAseq/241209_nutmeg_bg_clean_metadata.csv"
# cell_meta = pd.read_csv(path_cell_meta)
# cell_meta.set_index('cell_idx', inplace=True)

# hmba_anno_cells_kept = pd.DataFrame({
#     'Group': np.unique(cell_meta.Group),
#     'Kept': [sum(cell_meta.index.isin(hmba_bc_kept) & (cell_meta.Group == grp)) for grp in np.unique(cell_meta.Group)],
#     'Total': [sum(cell_meta.Group == grp) for grp in np.unique(cell_meta.Group)]
# })
# hmba_anno_cells_kept["Fraction"] = hmba_anno_cells_kept["Kept"] / hmba_anno_cells_kept["Total"]

# >>> hmba_anno_cells_kept.sort_values("Fraction", ascending = False)
#                                Group  Kept  Total  Fraction
# 38              STR Cholinergic GABA   533    533  1.000000 [*]
# 40                STR LHX8 ST18 GABA    78     78  1.000000 [*]
# 3               BN LAMP5 CXCL14 GABA   891    894  0.996644
# 39            STR FS PTHLH ST18 GABA  1573   1584  0.993056 [*]
# 45                      STR VIP GABA  3061   3090  0.990615 [*]
# 4                 BN LAMP5 LHX6 GABA  1387   1401  0.990007
# 44               STR TAC3 PLPP4 GABA  1263   1277  0.989037 [*]
# 42                STR SST Chodl GABA   805    817  0.985312 [*]
# 6                 CTX PVALB ChC GABA   407    427  0.953162
# 24                         OT D1-ICj   160    352  0.454545
# 5                      BN MEIS2 GABA   416   1889  0.220222
# 59                       low_quality    25    136  0.183824
# 53                     STRv D1 Shell   617   5035  0.122542
# 21             NDB SI LHX6 LHX8 GBX1    72   1189  0.060555
# 35                   SN STH VTA Glut   343   7943  0.043183
# 41               STR SST ADARB2 GABA     3     70  0.042857 [!]
# 15                     Hipp MGE GABA    80   2036  0.039293
# 30                       SN STH GABA   259   6713  0.038582
# 57                        Unassigned    60   6388  0.009393
# 43                STR SST RSPO2 GABA    21   2350  0.008936
# 22                  NDB SI LHX8 GABA     3    362  0.008287
# 28                      SLC17A6 Glut    23   4006  0.005741
# 11                GP MEIS2 SOX6 GABA     1    195  0.005128
# 52                     STRv D1 NUDAP    12   2685  0.004469
# 54                     STRv D2 Shell    11   2945  0.003735
# 33                  SN STH PAX8 GABA     2    788  0.002538
# 48                  STRd D1D2 Hybrid     5   2248  0.002224
# 36                        SN TH Dopa     1    586  0.001706
# 31          SN STH GATA3 TCF7L2 GABA     1   1828  0.000547
# 47                 STRd D1 Striosome     1   2611  0.000383
# 16         Hypothalamus SLC17A6 GAD1     1   2889  0.000346
# 46                    STRd D1 Matrix     2   9240  0.000216
# 18                 Lower Rhombic Lip     0   2938  0.000000
# 7                     Choroid Plexus     0    387  0.000000
# 51                 STRd D2 Striosome     0   1699  0.000000
# 49                    STRd D2 Matrix     0  11331  0.000000
# 50           STRd D2 StrioMat Hybrid     0    382  0.000000
# 2                       BN Astrocyte     0   4556  0.000000
# 55                TH MEIS2 OTX2 GABA     0   7067  0.000000
# 56                     Thalamus Glut     0    238  0.000000
# 58                              VLMC     0    255  0.000000
# 8                               Endo     0     61  0.000000
# 12                GP SOX6 KCNA1 GABA     0    306  0.000000
# 9                          Ependymal     0    244  0.000000
# 29               SN PVALB GATA3 GABA     0     85  0.000000
# 20                 Myelinating Oligo     0    219  0.000000
# 17         Hypothalamus SLC17A6 Glut     0    419  0.000000
# 23                               OPC     0   7645  0.000000
# 25                      Oligo OPALIN     0  32847  0.000000
# 26                      Oligo RBFOX1     0  15910  0.000000
# 27              Premyelinating Oligo     0    230  0.000000
# 1                      BG SKOR1 Glut     0    234  0.000000
# 10                         GAD1 GABA     0    345  0.000000
# 32  SN STH Hypothalamus SLC17A6 Glut     0   2921  0.000000
# 34           SN STH PVALB PITX2 Glut     0    567  0.000000
# 14                         GPi Shell     0     80  0.000000
# 37                     STR Astrocyte     0   6679  0.000000
# 13                          GPi Core     0    185  0.000000
# 19                         Microglia     0   6189  0.000000
# 0                                BAM     0    220  0.000000

### basically the same as when not all subc055 included, to a very surprising extent


### I'll show the follow up on 'STR SST RSPO2 GABA', but it looks the same as previous integrations

# >>> adata_norm.obs.loc[cell_meta.index[cell_meta.Group == 'STR SST ADARB2 GABA'], ['leiden', 'experiment_entropy']].value_counts()
# leiden  experiment_entropy
# 3       0.062095              41
# 13      0.023205              23
# 18      1.331140               3
# 43      0.094610               3
# Name: count, dtype: int64

### very similar cluster entropies as previously

# >>> adata_norm.obs.loc[adata_norm.obs.leiden == '3', 'experiment_clust_anno'].value_counts()
# experiment_clust_anno
# marm_hmba_nan          5063
# abc_nan                  40
# marm_PTHLH/PVALB          6
# marm_SST/NPY              3
# siletti_SST/NPY           2

# >>> cell_meta.loc[adata_norm.obs.index[(adata_norm.obs.leiden == '3') & (adata_norm.obs.index.isin(cell_meta.index))], "Group"].value_counts()
# Group
# STR SST RSPO2 GABA        2328
# Hipp MGE GABA             1956
# Unassigned                  60
# STR SST ADARB2 GABA         41
# low_quality                 25
# CTX PVALB ChC GABA          19
# NDB SI LHX6 LHX8 GBX1       16
# BN LAMP5 LHX6 GABA          14
# STR FS PTHLH ST18 GABA      10
# SN STH GABA                  4
# STR VIP GABA                 2
# STR SST Chodl GABA           2
# Name: count, dtype: int64


### again almost the exact same result

###
###
### === Will just keep these notes from previous integrations:
###
### There are a lot more marm_hmba cells in that cluster...
### and in previous integration 'STR SST ADARB2 GABA' actually didn't end up integrating with anything
### 
### And then (below) a lot of the other cells in that cluster are 'STR SST RSPO2 GABA', which only 
### integrated with siletti_Mixed (not that many cells) last time
###
### So even though, on appearance, we want those cells in principle, they don't integrate now and they won't integrate later
### [and that's been across a bunch of variations of these integrations]
###
### ===


# --------------------------------------------------
# Check cluster entropy optimization: Check HMBA cells being kept/dropped for higher/lower cutoffs than 0.1
# --------------------------------------------------

# df_entropy_counts = adata_norm.obs.loc[:, ['leiden', 'experiment_entropy']].drop_duplicates().set_index('leiden').sort_values('leiden')
# df_entropy_counts['non_hmba_cells'] = adata_norm.obs.loc[adata_norm.obs.experiment != 'marm_hmba', 'leiden'].value_counts(sort=False)
# df_entropy_counts['hmba_cells'] = adata_norm.obs.loc[adata_norm.obs.experiment == 'marm_hmba', 'leiden'].value_counts(sort=False)

# >>> df_entropy_counts.sort_values('experiment_entropy', ascending=False)
#         experiment_entropy  non_hmba_cells  hmba_cells
# leiden
# 18                1.331140            1481         953
# 0                 1.171024             575         663
# 21                0.758220             326        1100
# 19                0.736158            4694           8
# 17                0.698127            1097        1506
# 15                0.685156            1203        1965
# 12                0.469051             575        2647
# 22                0.289799             291        3874
# 16                0.148232              71        2145
# 43                0.094610              31        1590 [*]
# 3                 0.062095              51        5063 [*] already checked
# 24                0.038404            2412           8 
# 1                 0.038157              17        2932
# 13                0.023205              13        3998
# 20                0.016557               3        1275
# 2                 0.013465               3        1622
# 32                0.013182               4        2413
# 36                0.011850               1         627
# 23                0.011034              12        8964
# 39                0.008629               5        4519
# 25                0.005236               4        6993
# 28                0.003853               3        6791
# 31                0.003378               1        2626
# 40                0.001439               1        6827
# 47                0.001098               1        9221
# 46                0.000917               1       11262
# 58                0.000000               0        4818
# 61                0.000000               0         830
# 60                0.000000               0          25
# 59                0.000000               0        1509
# 14                0.000000               0         266
# 57                0.000000               0        6557
# 63                0.000000               0         242
# 56                0.000000               0         653
# 55                0.000000               0        2424
# 54                0.000000               0       18532
# 62                0.000000               0         704
# 65                0.000000               0         245
# 64                0.000000               0        4776
# 52                0.000000               0       34498
# 66                0.000000               0         360
# 67                0.000000               0         266
# 68                0.000000               0         452
# 69                0.000000               0         217
# 70                0.000000               0         792
# 71                0.000000               0         215
# 72                0.000000               0         117
# 73                0.000000               0         288
# 74                0.000000               0        4136
# 75                0.000000               0         397
# 76                0.000000               0         988
# 53                0.000000               0        1348
# 49                0.000000               0        3919
# 51                0.000000               0        7717
# 35                0.000000               0        3490
# 10                0.000000               0        3031
# 9                 0.000000               0        2665
# 8                 0.000000               0        1428
# 26                0.000000               0         626
# 27                0.000000               0          55
# 29                0.000000               0        6060
# 30                0.000000               0         726
# 7                 0.000000               0        2901
# 33                0.000000               0        2145
# 34                0.000000               0         663
# 6                 0.000000               0         168
# 50                0.000000               0        7448
# 37                0.000000               0         728
# 38                0.000000               0         637
# 5                 0.000000               0       12894
# 41                0.000000               0        5418
# 42                0.000000               0       15096
# 4                 0.000000               0        5834
# 44                0.000000               0        2375
# 45                0.000000               0        1548
# 48                0.000000               0         995
# 11                0.000000               0       16957
# 77                0.000000               0          35


### lower values, not much to check
# 
# cl 43 = 0.094610
# 
# >>> cell_meta.loc[adata_norm.obs.index[(adata_norm.obs.leiden == '43') & (adata_norm.obs.index.isin(cell_meta.index))], "Group"].value_counts()
# Group
# NDB SI LHX6 LHX8 GBX1       1055
# GP SOX6 KCNA1 GABA           306
# NDB SI LHX8 GABA              17
# SN STH GABA                   12
# Unassigned                     6
# STR SST Chodl GABA             6
# STR SST ADARB2 GABA            3
# SN STH GATA3 TCF7L2 GABA       3
# STR VIP GABA                   3
# SN STH VTA Glut                2
# STRv D1 NUDAP                  1
# Name: count, dtype: int64
#
### basically the same 'next lowest' cluster as previous integration

### higher vallues
#
# cl 16 = 0.148232
# 
# >>> cell_meta.loc[adata_norm.obs.index[(adata_norm.obs.leiden == '16') & (adata_norm.obs.index.isin(cell_meta.index))], "Group"].value_counts()
# Group
# BN LAMP5 LHX6 GABA        1383
# CTX PVALB ChC GABA         406
# STR FS PTHLH ST18 GABA       9
# Hipp MGE GABA                5
# low_quality                  5
# Unassigned                   3
# STR SST RSPO2 GABA           1
# STR VIP GABA                 1
# BN LAMP5 CXCL14 GABA         1
# Name: count, dtype: int64
#
### again effectively the same 'next higher' cluster as previous integration;
### that time I left in the 'BN LAMP5 LHX6 GABA', which actually did pull over some marm census CCK cells...
###

######## -> so this integration (and previous, that had only abc cells from STRv/STRd) could be redone with 
########    a higher initial entropy cutoff for excluding marm_hmba cells. but whatever, I'll keep this pair
########    being similar, it doesn't seem problematic that they're there...



# --------------------------------------------------
# Cluster entropy: plots before subsetting
# --------------------------------------------------

plotdir = Path(outdir, "scvi_model_initial_int_plots")
if not Path.is_dir(plotdir):
    os.mkdir(plotdir)

sc.settings.figdir = str(plotdir)


### each species ---
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment", 
    save = "_experiment.pdf",
    show = False
)
# highlighting the non-marm-hmba:
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment",
    groups = ['marm', 'abc', 'siletti'],
    save = "_experiment_nonhmba.pdf",
    show = False
)

sc.pl.embedding(
    adata_norm, 
    "X_umap_integrated", 
    color = "leiden", 
    save = "_leiden.pdf",
    show = False
)

sc.pl.embedding(
    adata_norm, 
    "X_umap_integrated", 
    color = "toRemove", 
    save = "_toRemove.pdf",
    show = False
)



### each species clust_anno ---

marm_clust = set(adata.obs.loc[
    adata.obs.index[adata.obs['experiment'].isin(["marm"])],
    'experiment_clust_anno'
])
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment_clust_anno",
    groups = marm_clust,
    palette = plt.get_cmap("tab20").colors[0:7],
    save = "_marm_clust_anno.pdf",
    show = False
)


# # sc.pl.embedding(
# #     adata, 
# #     "X_umap_integrated", 
# #     color = "experiment_clust_anno",
# #     groups = set(adata.obs['experiment_clust_anno'][adata.obs['experiment_clust_anno'].str.contains("marm_hmba")]),
# #     palette = plt.get_cmap("tab20").colors[0:9],
# #     save = "_marm_hmba_clust_anno.pdf",
# #     show = False
# # )
# adata_norm.obs = adata_norm.obs.join(pd.DataFrame(cell_meta.loc[:, 'Group']).rename(columns = {'Group': 'hmba_group'}))
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     # palette = plt.get_cmap("tab20b").colors[0:24],
#     save = "_hmba_group.pdf",
#     show = False
# )
# # one color is reused (BN LAMP5 CXCL14 GABA and 'Unassigned' (36 cells))
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     groups = ['Unassigned'],
#     # palette = plt.get_cmap("tab20b").colors[0:22],
#     save = "_hmba_group_unassigned.pdf",
#     show = False
# )


sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment_clust_anno",
    groups = set(adata.obs['experiment_clust_anno'][adata.obs['experiment_clust_anno'].str.contains("siletti")]),
    palette = plt.get_cmap("tab20").colors[0:9],
    save = "_siletti_clust_anno.pdf",
    show = False
)



# abc subclass
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "subclass",
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_subclass.pdf",
    show = False
)

# abc supertype
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "supertype",
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_supertype.pdf",
    show = False
)

# abc clusters within subclass 055
sc_clusters = set(adata.obs.loc[
    adata.obs.index[adata.obs['subclass'].isin(["055 STR Lhx8 Gaba"])],
    'experiment_cluster'
])

sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "experiment_cluster",
    groups = sc_clusters,
    # palette = "tab20b",
    palette = plt.get_cmap("tab20").colors[2:9],
    save = "_abc_cluster_in_subclass_055.pdf",
    show = False
)

# abc roi
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "region_of_interest_acronym",
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_roi.pdf",
    show = False
)

# siletti roi
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "roi",
    palette = plt.get_cmap("tab20b").colors,
    save = "_siletti_roi.pdf",
    show = False
)



# --------------------------------------------------
# Save normalized h5ad (has info for subsetting)
# --------------------------------------------------

# This contains the leiden clusters, the entropy calculations, etc, because i don't remember if that's all deterministic or not...

# adata_norm.write_h5ad(Path(outdir, "adata_norm_model_initial_int.h5ad"))

# reimport
# adata_norm = ad.read_h5ad(Path(outdir, "adata_norm_model_initial_int.h5ad"))


# --------------------------------------------------
# Cluster entropy: Use threshold of 0.1 to drop marm_hmba cells clustering with themselves
# --------------------------------------------------

idx_keep = np.where(~((adata_norm.obs.toRemove) & (adata_norm.obs.experiment == "marm_hmba")))[0]

### basically the same as w/o all subc055
# >>> len(idx_keep)
# 27737
# >>> adata_norm.shape[0]
# 295682
# >>> sum(adata_norm.obs.experiment == "marm_hmba")
# 282806
# >>> adata_norm.shape[0] - len(idx_keep)
# 267945
# 
### dropping 267945/282806 = 94.7% of marm_hmba cells

# >>> adata.shape
# (295682, 4000)
# >>> adata.raw.shape
# (295682, 11280)

# adata.raw.X is float64; make sure didn't just make ref to adata.raw when making adata_norm
assert (adata.raw.X.astype(int) != adata.raw.X).nnz == 0

adata = ad.AnnData(
    adata.raw.X[idx_keep, :].astype(int), 
    obs = adata.obs.iloc[idx_keep, :], 
    var = adata.raw.var
)

# >>> adata.shape
# (27737, 11280)


# --------------------------------------------------
# scVI model
# --------------------------------------------------

adata.raw = adata
sc.pp.highly_variable_genes(
    adata, 
    flavor      = "seurat_v3", 
    n_top_genes = 4000, 
    batch_key   = "experiment", 
    subset      = True
)

print("Saving merged h5ad with HVGs selected...")
adata.write_h5ad(Path(outdir, "adata_firstprune.h5ad"))

## re-import
# adata = ad.read_h5ad(Path(outdir, "adata_firstprune.h5ad"))

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key = "experiment",
    categorical_covariate_keys = ['donor_coding']
)
scvi_model_firstprune = scvi.model.SCVI(
    adata, 
    n_hidden   = 256,
    n_latent   = 64,
    n_layers   = 3,
    dispersion = "gene-batch"
)

print("Training model (firstprune)...")
scvi_model_firstprune.train(accelerator = "gpu", max_epochs = 250)
scvi_model_firstprune.save(Path(outdir, "scvi_model_firstprune"))

## re-import
# scvi_model_firstprune = scvi.model.SCVI.load(Path(outdir, "scvi_model_firstprune"), adata) 


# --------------------------------------------------
# Get, save scVI embeddings
# --------------------------------------------------

adata.obsm["X_scVI"] = scvi_model_firstprune.get_latent_representation()
sc.pp.neighbors(adata, use_rep = "X_scVI")
adata.obsm["X_umap_integrated"] = sc.tl.umap(adata, min_dist = 0.25, copy = True).obsm["X_umap"]
np.savetxt(
    Path(outdir, "scvi_model_firstprune_X_scVI.csv"), 
    adata.obsm["X_scVI"], 
    delimiter = ","
)
np.savetxt(
    Path(outdir, "scvi_model_firstprune_X_umap_integrated.csv"), 
    adata.obsm["X_umap_integrated"], 
    delimiter = ","
)

## re-import
# adata.obsm["X_scVI"] = np.genfromtxt(Path(outdir, "scvi_model_firstprune_X_scVI.csv"), delimiter = ',')
# adata.obsm["X_umap_integrated"] = np.genfromtxt(Path(outdir, "scvi_model_firstprune_X_umap_integrated.csv"), delimiter = ',')


# --------------------------------------------------
# Export h5ad 
# --------------------------------------------------

# adata.write_h5ad(Path(outdir, "adata_model_firstprune.h5ad"))

## re-import
# adata = ad.read_h5ad(Path(outdir, "adata_model_firstprune.h5ad"))


# --------------------------------------------------
# Normalize and recluster
# --------------------------------------------------

adata_norm = ad.AnnData(
    X    = adata.raw.X.copy(),
    obs  = adata.obs,
    var  = adata.raw.var,
    uns  = adata.uns,
    obsm = adata.obsm,
    # raw = adata.raw
)
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)
sc.pp.neighbors(adata_norm, n_neighbors=10, use_rep="X_scVI")
sc.tl.leiden(adata_norm, flavor="igraph", n_iterations=2)


# --------------------------------------------------
# Plots after first (HMBA-only) pruning
# --------------------------------------------------

plotdir = Path(outdir, "scvi_model_firstprune_plots")
if not Path.is_dir(plotdir):
    os.mkdir(plotdir)

sc.settings.figdir = str(plotdir)


### each species ---
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment", 
    save = "_experiment.pdf",
    show = False
)
# highlighting the non-marm-hmba:
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment",
    groups = ['marm', 'abc', 'siletti'],
    save = "_experiment_nonhmba.pdf",
    show = False
)

sc.pl.embedding(
    adata_norm, 
    "X_umap_integrated", 
    color = "leiden", 
    save = "_leiden.pdf",
    show = False
)


### each species clust_anno ---

marm_clust = set(adata.obs.loc[
    adata.obs.index[adata.obs['experiment'].isin(["marm"])],
    'experiment_clust_anno'
])
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment_clust_anno",
    groups = marm_clust,
    palette = plt.get_cmap("tab20").colors[5:12],
    save = "_marm_clust_anno.pdf",
    show = False
)

sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment_clust_anno",
    groups = set(adata.obs['experiment_clust_anno'][adata.obs['experiment_clust_anno'].str.contains("siletti")]),
    # palette = plt.get_cmap("tab20").colors[0:9],
    save = "_siletti_clust_anno.pdf",
    show = False
)


# adata_norm.obs = adata_norm.obs.join(pd.DataFrame(cell_meta.loc[:, 'Group']).rename(columns = {'Group': 'hmba_group'}))
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     # palette = plt.get_cmap("tab20b").colors[0:24],
#     save = "_hmba_group.pdf",
#     show = False
# )
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     groups = ['Unassigned'],
#     # palette = plt.get_cmap("tab20b").colors[0:22],
#     save = "_hmba_group_unassigned.pdf",
#     show = False
# )
# hmba_groups = adata_norm.obs.hmba_group.cat.categories
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     groups = list(hmba_groups[hmba_groups.str.startswith("STR")]),
#     # palette = plt.get_cmap("tab20b").colors[0:15],
#     save = "_hmba_group_striatal.pdf",
#     show = False
# )
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     groups = list(hmba_groups[~hmba_groups.str.startswith("STR")]),
#     # palette = plt.get_cmap("tab20b").colors[0:15],
#     save = "_hmba_group_nonstriatal.pdf",
#     show = False
# )


# abc class
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "class",
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_class.pdf",
    show = False
)

# abc subclass
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "subclass",
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_subclass.pdf",
    show = False
)

# abc subclasses within each abc class: 
# ['06 CTX-CGE GABA', '07 CTX-MGE GABA', '08 CNU-MGE GABA', '10 LSX GABA', '11 CNU-HYa GABA']
for abc_cl in list(adata.obs['class'].drop_duplicates().dropna().sort_values()):
    subc = set(adata.obs.loc[adata.obs['class'] == abc_cl, 'subclass'])
    sc.pl.embedding(
        adata,
        "X_umap_integrated", 
        color = "subclass",
        groups = subc,
        palette = plt.get_cmap("tab20b").colors,
        save = f"_abc_class_{abc_cl}_subclasses.pdf",
        show = False
    )

# abc supertype
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "supertype",
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_supertype.pdf",
    show = False
)

# abc supertypes within class8
st = set(adata.obs.loc[adata.obs['class'] == '08 CNU-MGE GABA', 'supertype'])
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "supertype",
    groups = st,
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_supertype_in_class8.pdf",
    show = False
)


# abc clusters within subclass 055
sc_clusters = set(adata.obs.loc[adata.obs['subclass'] == '055 STR Lhx8 Gaba', 'experiment_cluster'])
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "experiment_cluster",
    groups = sc_clusters,
    # palette = "tab20b",
    palette = plt.get_cmap("tab20").colors[5:13],
    save = "_abc_cluster_in_subclass_055.pdf",
    show = False
)


# anatomical division label
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "anatomical_division_label",
    palette = plt.get_cmap("tab20b").colors[6:9],
    save = "_abc_anatomical_division_label.pdf",
    show = False
)
# abc roi
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "region_of_interest_acronym",
    palette = plt.get_cmap("tab20b").colors[6:9],
    save = "_abc_roi.pdf",
    show = False
)


### siletti dissection region ---
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "experiment_roi",
    groups = set(adata.obs['experiment_roi'][adata.obs['experiment_roi'].str.contains("siletti")]),
    palette = plt.get_cmap("tab20b").colors[2:6],
    save = "_siletti_roi.pdf",
    show = False
)


### plot ABC donors
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "donor_coding",
    groups = set(adata.obs.loc[adata.obs.experiment == 'abc', 'donor_coding']),
    save = "_abc_donor_coding.pdf",
    show = False
)


### plot marm donors
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "donor_coding",
    groups = set(adata.obs.loc[adata.obs.experiment == 'marm', 'donor_coding']),
    palette = plt.get_cmap("tab20b").colors[2:6],
    save = "_marm_donor_coding.pdf",
    show = False
)

### plot siletti donors
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "donor_coding",
    groups = set(adata.obs.loc[adata.obs.experiment == 'siletti', 'donor_coding']),
    save = "_siletti_donor_coding.pdf",
    show = False
)


### plot marm_hmba donors
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "donor_coding",
    groups = set(adata.obs.loc[adata.obs.experiment == 'marm_hmba', 'donor_coding']),
    palette = plt.get_cmap("tab20b").colors[2:6],
    save = "_marm_hmba_donor_coding.pdf",
    show = False
)




# --------------------------------------------------
# Feature plots (requires normalized adata object)
# --------------------------------------------------

genedir = Path(plotdir, "genes")
if not Path.is_dir(genedir):
    os.mkdir(genedir)

sc.settings.figdir = str(genedir)

gn_to_plot = [
    "TH",
    "TAC3",
    "LHX6",
    "LHX8",
    "PLPP4",
    "PROX1",
    "PVALB",
    "PTHLH",
    "SST",
    "CHAT",
    "VIP",
    "LAMP5",
    "MEIS2",
    "SOX6"
]

for gn in gn_to_plot:
    sc.pl.embedding(
        adata_norm,
        "X_umap_integrated", 
        use_raw = False,
        color = gn,
        save = f"_{gn}.pdf",
        show = False
    )
    for exp in ['marm', 'abc', 'siletti', 'marm_hmba']:
        sc.pl.embedding(
            adata_norm,
            "X_umap_integrated", 
            mask_obs = adata_norm.obs['experiment'] == exp,
            use_raw = False,
            color = gn,
            save = f"_{gn}_{exp}.pdf",
            show = False
        )

sc.settings.figdir = str(plotdir)


# --------------------------------------------------
# Cluster pruning
# --------------------------------------------------

# Clearly some the ABC and marm_hmba data is Meis2+, which are these:
cl_meis2 = ['9', '16', '20', '21', '31', '32']

# sc.pl.heatmap(adata_norm, 'MEIS2', groupby = 'leiden')
# 
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "leiden", 
#     palette = plt.get_cmap("tab20b").colors,
#     groups = cl_meis2
# )
#
# -> those are also the source of abc cells from classes '10 LSX GABA' and '11 CNU-HYa GABA'


idx_keep = np.where(~adata_norm.obs.leiden.isin(cl_meis2))[0]

# >>> adata_norm.shape[0]
# 27737
# >>> len(idx_keep)
# 24441
# >>> adata_norm.shape[0] - len(idx_keep)
# 3296

adata = ad.AnnData(
    adata.raw.X[idx_keep, :].astype(int), 
    obs = adata.obs.iloc[idx_keep, :], 
    var = adata.raw.var
)

# >>> adata.shape
# (24441, 11280)


# --------------------------------------------------
# scVI model
# --------------------------------------------------

adata.raw = adata
sc.pp.highly_variable_genes(
    adata, 
    flavor = "seurat_v3", 
    n_top_genes = 4000, 
    batch_key = "experiment", 
    subset = True
)

print("Saving merged h5ad (secondprune) with HVGs selected...")
adata.write_h5ad(Path(outdir, "adata_secondprune.h5ad"))

## re-import
# adata = ad.read_h5ad(Path(outdir, "adata_secondprune.h5ad"))

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key = "experiment",
    categorical_covariate_keys = ['donor_coding']
)
scvi_model_secondprune = scvi.model.SCVI(
    adata, 
    n_hidden = 256,
    n_latent = 64,
    n_layers = 3,
    dispersion = "gene-batch"
)

print("Training model (secondprune)...")
scvi_model_secondprune.train(accelerator = "gpu", max_epochs = 250)
scvi_model_secondprune.save(Path(outdir, "scvi_model_secondprune"))

## re-import
# scvi_model_secondprune = scvi.model.SCVI.load(Path(outdir, "scvi_model_secondprune"), adata) 


# --------------------------------------------------
# Get, save scVI embeddings
# --------------------------------------------------

adata.obsm["X_scVI"] = scvi_model_secondprune.get_latent_representation()
sc.pp.neighbors(adata, use_rep = "X_scVI")
adata.obsm["X_umap_integrated"] = sc.tl.umap(adata, min_dist = 0.25, copy = True).obsm["X_umap"]
np.savetxt(
    Path(outdir, "scvi_model_secondprune_X_scVI.csv"), 
    adata.obsm["X_scVI"], 
    delimiter = ","
)
np.savetxt(
    Path(outdir, "scvi_model_secondprune_X_umap_integrated.csv"), 
    adata.obsm["X_umap_integrated"], 
    delimiter = ","
)

## re-import
# adata.obsm["X_scVI"] = np.genfromtxt(Path(outdir, "scvi_model_secondprune_X_scVI.csv"), delimiter = ',')
# adata.obsm["X_umap_integrated"] = np.genfromtxt(Path(outdir, "scvi_model_secondprune_X_umap_integrated.csv"), delimiter = ',')


# --------------------------------------------------
# Export h5ad 
# --------------------------------------------------

# adata.write_h5ad(Path(outdir, "adata_model_secondprune.h5ad"))

## re-import
# adata = ad.read_h5ad(Path(outdir, "adata_model_secondprune.h5ad"))


# --------------------------------------------------
# Normalize and recluster
# --------------------------------------------------

# reclustering this time may not be necessary
adata_norm = ad.AnnData(
    X    = adata.raw.X.copy(),
    obs  = adata.obs,
    var  = adata.raw.var,
    uns  = adata.uns,
    obsm = adata.obsm,
    # raw  = adata.raw
)
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)
sc.pp.neighbors(adata_norm, n_neighbors=10, use_rep="X_scVI")
sc.tl.leiden(adata_norm, flavor="igraph", n_iterations=2)


# --------------------------------------------------
# Plots after second pruning
# --------------------------------------------------

plotdir = Path(outdir, "scvi_model_secondprune_plots")
if not Path.is_dir(plotdir):
    os.mkdir(plotdir)

sc.settings.figdir = str(plotdir)


### each species ---
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment", 
    save = "_experiment.pdf",
    show = False
)
# highlighting the non-marm-hmba:
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment",
    groups = ['marm', 'abc', 'siletti'],
    save = "_experiment_nonhmba.pdf",
    show = False
)

sc.pl.embedding(
    adata_norm, 
    "X_umap_integrated", 
    color = "leiden", 
    save = "_leiden.pdf",
    show = False
)


### each species clust_anno ---

marm_clust = set(adata.obs.loc[
    adata.obs.index[adata.obs['experiment'].isin(["marm"])],
    'experiment_clust_anno'
])
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment_clust_anno",
    groups = marm_clust,
    palette = plt.get_cmap("tab20").colors[5:12],
    save = "_marm_clust_anno.pdf",
    show = False
)

sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "experiment_clust_anno",
    groups = set(adata.obs['experiment_clust_anno'][adata.obs['experiment_clust_anno'].str.contains("siletti")]),
    # palette = plt.get_cmap("tab20").colors[0:9],
    save = "_siletti_clust_anno.pdf",
    show = False
)


# adata_norm.obs = adata_norm.obs.join(pd.DataFrame(cell_meta.loc[:, 'Group']).rename(columns = {'Group': 'hmba_group'}))
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     # palette = plt.get_cmap("tab20b").colors[0:24],
#     save = "_hmba_group.pdf",
#     show = False
# )
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     groups = ['Unassigned'],
#     # palette = plt.get_cmap("tab20b").colors[0:22],
#     save = "_hmba_group_unassigned.pdf",
#     show = False
# )
# hmba_groups = adata_norm.obs.hmba_group.cat.categories
# hmba_groups_str = list(hmba_groups[hmba_groups.str.startswith("STR")])
# hmba_groups_nonstr = list(hmba_groups[~hmba_groups.str.startswith("STR")])
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     groups = hmba_groups_str,
#     palette = plt.get_cmap("tab20b").colors[0:15],
#     save = "_hmba_group_striatal.pdf",
#     show = False
# )
# sc.pl.embedding(
#     adata_norm, 
#     "X_umap_integrated", 
#     color = "hmba_group", 
#     groups = hmba_groups_nonstr,
#     palette = plt.get_cmap("tab20b").colors[0:15],
#     save = "_hmba_group_nonstriatal.pdf",
#     show = False
# )


# abc class
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "class",
    palette = plt.get_cmap("tab20b").colors[2:7],
    save = "_abc_class.pdf",
    show = False
)

# abc subclass
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "subclass",
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_subclass.pdf",
    show = False
)

# abc subclasses within each abc class: 
# ['06 CTX-CGE GABA', '07 CTX-MGE GABA', '08 CNU-MGE GABA', '10 LSX GABA', '11 CNU-HYa GABA']
for abc_cl in list(adata.obs['class'].drop_duplicates().dropna().sort_values()):
    subc = set(adata.obs.loc[adata.obs['class'] == abc_cl, 'subclass'])
    sc.pl.embedding(
        adata,
        "X_umap_integrated", 
        color = "subclass",
        groups = subc,
        palette = plt.get_cmap("tab20b").colors,
        save = f"_abc_class_{abc_cl}_subclasses.pdf",
        show = False
    )

# abc supertype
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "supertype",
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_supertype.pdf",
    show = False
)

# abc supertypes within class8
st = set(adata.obs.loc[adata.obs['class'] == '08 CNU-MGE GABA', 'supertype'])
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "supertype",
    groups = st,
    palette = plt.get_cmap("tab20b").colors,
    save = "_abc_supertype_in_class8.pdf",
    show = False
)


# abc clusters within subclass 055
sc_clusters = set(adata.obs.loc[adata.obs['subclass'] == '055 STR Lhx8 Gaba', 'experiment_cluster'])
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "experiment_cluster",
    groups = sc_clusters,
    # palette = "tab20b",
    palette = plt.get_cmap("tab20").colors[5:13],
    save = "_abc_cluster_in_subclass_055.pdf",
    show = False
)


# anatomical division label
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "anatomical_division_label",
    palette = plt.get_cmap("tab20b").colors[2:5],
    save = "_abc_anatomical_division_label.pdf",
    show = False
)
# abc roi
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "region_of_interest_acronym",
    palette = plt.get_cmap("tab20b").colors[1:6],
    save = "_abc_roi.pdf",
    show = False
)

### siletti dissection region ---
sc.pl.embedding(
    adata,
    "X_umap_integrated", 
    color = "experiment_roi",
    groups = set(adata.obs['experiment_roi'][adata.obs['experiment_roi'].str.contains("siletti")]),
    palette = plt.get_cmap("tab20b").colors[2:6],
    save = "_siletti_roi.pdf",
    show = False
)


### plot ABC donors
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "donor_coding",
    groups = set(adata.obs.loc[adata.obs.experiment == 'abc', 'donor_coding']),
    save = "_abc_donor_coding.pdf",
    show = False
)


### plot marm donors
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "donor_coding",
    groups = set(adata.obs.loc[adata.obs.experiment == 'marm', 'donor_coding']),
    palette = plt.get_cmap("tab20b").colors[2:6],
    save = "_marm_donor_coding.pdf",
    show = False
)

### plot siletti donors
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "donor_coding",
    groups = set(adata.obs.loc[adata.obs.experiment == 'siletti', 'donor_coding']),
    save = "_siletti_donor_coding.pdf",
    show = False
)


### plot marm_hmba donors
sc.pl.embedding(
    adata, 
    "X_umap_integrated", 
    color = "donor_coding",
    groups = set(adata.obs.loc[adata.obs.experiment == 'marm_hmba', 'donor_coding']),
    palette = plt.get_cmap("tab20b").colors[2:6],
    save = "_marm_hmba_donor_coding.pdf",
    show = False
)




# --------------------------------------------------
# Feature plots (requires normalized adata object)
# --------------------------------------------------

genedir = Path(plotdir, "genes")
if not Path.is_dir(genedir):
    os.mkdir(genedir)

sc.settings.figdir = str(genedir)

gn_to_plot = [
    "DLX5",
    "DLX1",   
    "DLX2",   
    "DLX6",   
    "GAD1",   
    "GAD2",   
    "SOX6",   
    "LHX6",   
    "NKX2-1",
    "MEF2C",  
    "PAX6",   
    "MEIS2",  
    "NPY1R",  
    "ISL1",   
    "FOXP1",  
    "FOXP2",  
    "TSHZ1",  
    "PENK",   
    "SCGN",   
    "LAMP5",  
    "ZIC1",   
    "ZIC2",   
    "ZIC4",   
    "HAP1",   
    "CHAT",   
    "NR2F2",  
    "PROX1",  
    "CCK",    
    "MKI67",  
    "ADARB2", 
    "PDE3A",  
    "TNS1",   
    "STXBP6", 
    "KIT",    
    "ETV1",   
    "CER1",   
    "TH",     
    "NXPH2",  
    "PLPP4",  
    "NXPH1",  
    "DRD2",   
    "MAF",    
    "CHRNA3", 
    "COL19A1",
    "TAC3",   
    "RBP4",   
    "LHX8",   
    "NPY",    
    "SST",    
    "TRHDE",  
    "SIX3",   
    "TBR1",   
    "VIP",    
    "HES5",   
    "ERBB4",  
    "CHRNA7", 
    "PVALB",  
    "PTHLH",  
    "VIPR2" 
]

# >>> np.setdiff1d(gn_to_plot, adata_norm.var.index)
# array(['CCK', 'CER1', 'CHRNA7', 'HAP1', 'HES5', 'MKI67'], dtype='<U7')

gn_to_plot = list(np.intersect1d(gn_to_plot, adata_norm.var.index))

for gn in gn_to_plot:
    if Path(genedir, f"X_umap_integrated_{gn}.pdf").is_file():
        continue
    sc.pl.embedding(
        adata_norm,
        "X_umap_integrated", 
        use_raw = False,
        color = gn,
        save = f"_{gn}.pdf",
        show = False
    )
    for exp in ['marm', 'abc', 'siletti', 'marm_hmba']:
        sc.pl.embedding(
            adata_norm,
            "X_umap_integrated", 
            mask_obs = adata_norm.obs['experiment'] == exp,
            use_raw = False,
            color = gn,
            save = f"_{gn}_{exp}.pdf",
            show = False
        )

sc.settings.figdir = str(plotdir)


# --------------------------------------------------
# Annotate HMBA data & export annotations/metadata
# --------------------------------------------------

# Could run scANVI but I'm going to manually assign this based on the UMAP and whatever cells co-embed

### Again, went back and verified how my clusters relate to the HMBA taxonomy:

# >>> adata_norm.obs.loc[adata_norm.obs.experiment == "marm_hmba", ["hmba_group", "leiden"]].groupby("hmba_group").value_counts()
# hmba_group                leiden
# BN LAMP5 CXCL14 GABA      14         871
#                           15          15
#                           22           2
#                           25           1
#                           24           1
#                           17           1
# BN LAMP5 LHX6 GABA        17         733
#                           19         640
#                           21           7
#                           6            2
#                           7            1
#                           4            1
#                           14           1
#                           23           1
#                           16           1
# CTX PVALB ChC GABA        21         350
#                           19          54
#                           17           2
#                           5            1
# Hipp MGE GABA             5           68
#                           21           5
#                           4            4
#                           2            2
#                           18           1
# NDB SI LHX6 LHX8 GBX1     18          51
#                           20           6
#                           0            3
#                           9            2
#                           8            2
#                           4            1
#                           25           1
# SN STH GABA               15           5
#                           18           2
#                           8            1
# SN STH GATA3 TCF7L2 GABA  18           1
# SN STH PAX8 GABA          18           1
# SN TH Dopa                0            1
# STR Cholinergic GABA      0          348
#                           20         185
# STR FS PTHLH ST18 GABA    3          831
#                           4          297
#                           5          225
#                           6           82
#                           7           82
#                           8           18
#                           11          14
#                           17          13
#                           2            6
#                           18           4
#                           19           1
# STR LHX8 ST18 GABA        8           78
# STR SST ADARB2 GABA       18           2
#                           9            1
# STR SST Chodl GABA        9          559
#                           10         141
#                           18         103
#                           4            2
# STR SST RSPO2 GABA        10          11
#                           5            6
#                           9            2
#                           19           1
#                           23           1
# STR TAC3 PLPP4 GABA       11         964
#                           12         132
#                           8          102
#                           13          56
#                           16           4
#                           2            3
#                           7            2
# STR VIP GABA              23        1184
#                           22         806
#                           24         535
#                           25         268
#                           15         256
#                           14           9
#                           19           1
#                           2            1
#                           9            1
# STRv D1 NUDAP             15           9
#                           20           1
#                           2            1
# Unassigned                15          31
#                           2            6
#                           22           2
#                           19           2
#                           18           2
#                           25           1
#                           5            1
#                           17           1
# low_quality               2            5
#                           23           4
#                           21           3
#                           24           2
#                           19           2
#                           14           2
#                           4            1
#                           5            1
#                           22           1
#                           3            1
#                           16           1
#                           17           1


### but what I originally used (but mostly deferred to the UMAP plots):

# >>> adata_norm.obs.loc[:, ['experiment_clust_anno', 'leiden']].groupby('leiden').value_counts()
# leiden  experiment_clust_anno
# 0       marm_hmba_nan             438
#         marm_CHAT                 257
#         siletti_CHAT              195
#         abc_nan                    95
#         marm_TAC3                   6
# 1       abc_nan                    13
# 2       abc_nan                    77
#         marm_hmba_nan              50
#         marm_PTHLH/PVALB           18
#         marm_TAC3                  10
#         siletti_PTHLH/PVALB         8
#         siletti_TAC3                6
#         marm_SST/NPY                1
# 3       siletti_PTHLH/PVALB      1262
#         marm_hmba_nan             934
#         marm_PTHLH/PVALB          687
#         abc_nan                   491
#         siletti_TAC3                7
#         marm_TAC3                   2
#         marm_CHAT                   1
# 4       siletti_PTHLH/PVALB       401
#         marm_hmba_nan             366
#         abc_nan                   278
#         marm_PTHLH/PVALB          255
#         siletti_TAC3                3
#         marm_SST/NPY                1
# 5       marm_hmba_nan             403
#         marm_PTHLH/PVALB          121
#         siletti_PTHLH/PVALB        67
#         abc_nan                    34
#         marm_TAC3                   1
# 6       marm_hmba_nan              86
#         siletti_PTHLH/PVALB        72
#         marm_PTHLH/PVALB           60
# 7       siletti_PTHLH/PVALB       748
#         abc_nan                   175
#         marm_hmba_nan              91
#         marm_PTHLH/PVALB           63
#         siletti_TAC3               27
#         marm_TAC3                   7
#         marm_CCK                    1
# 8       marm_TAC3                 237
#         marm_hmba_nan             212
#         siletti_TAC3              196
#         abc_nan                    61
#         marm_PTHLH/PVALB           11
# 9       marm_hmba_nan             655
#         abc_nan                   496
#         marm_SST/NPY              475
#         siletti_SST/NPY           346
#         marm_CHAT                   1
# 10      marm_hmba_nan             175
#         siletti_SST/NPY            32
#         abc_nan                    20
#         marm_SST/NPY               11
#         marm_PTHLH/PVALB            1
# 11      marm_hmba_nan            1063
#         marm_TAC3                 613
#         siletti_TAC3              323
#         abc_nan                    24
#         marm_PTHLH/PVALB            1
# 12      siletti_TAC3             1803
#         marm_TAC3                 189
#         marm_hmba_nan             138
#         abc_nan                    17
#         marm_PTHLH/PVALB            1
# 13      siletti_TAC3               86
#         marm_hmba_nan              58
#         marm_TAC3                  34
# 14      marm_hmba_nan            1111
#         siletti_CCK               190
#         marm_CCK                   70
#         abc_nan                    64
#         siletti_CCK/VIP             1
# 15      marm_hmba_nan             420
#         siletti_CCK/VIP            32
# 16      abc_nan                   930
#         siletti_TAC3               10
#         marm_TAC3                  10
#         marm_hmba_nan               8
#         siletti_PTHLH/PVALB         1
# 17      marm_hmba_nan             963
#         marm_CCK                   64
#         siletti_PTHLH/PVALB         2
#         abc_nan                     2
#         marm_SST/NPY                1
#         marm_PTHLH/PVALB            1
# 18      marm_hmba_nan             198
#         abc_nan                    96
#         marm_SST/NPY               29
#         siletti_SST/NPY            21
#         marm_PTHLH/PVALB           13
#         siletti_PTHLH/PVALB         2
# 19      marm_hmba_nan             717
#         abc_nan                     4
#         marm_PTHLH/PVALB            1
# 20      marm_hmba_nan             207
#         marm_CHAT                  20
# 21      marm_hmba_nan             463
#         abc_nan                     3
#         siletti_PTHLH/PVALB         1
#         marm_PTHLH/PVALB            1
# 22      marm_hmba_nan             937
#         siletti_CCK/VIP            20
#         abc_nan                     5
# 23      marm_hmba_nan            1552
#         siletti_CCK/VIP            76
#         abc_nan                    14
# 24      marm_hmba_nan             656
#         siletti_CCK/VIP           101
#         abc_nan                    41
#         siletti_TAC3                1
# 25      marm_hmba_nan             281
#         abc_nan                     4



# 1 has only some ABC cells
leiden_dict = {
    '14' : 'CCK',
    '15' : 'CCK/VIP',
    '22' : 'CCK/VIP',
    '23' : 'CCK/VIP',
    '24' : 'CCK/VIP',
    '25' : 'CCK/VIP',
     '0' : 'CHAT',
    '20' : 'CHAT',
    '17' : 'LAMP5/LHX6',
    '19' : 'LAMP5/LHX6',
    '21' : 'LAMP5/LHX6',
     '9' : 'SST/NPY',
    '10' : 'SST/NPY',
    '18' : 'SST/NPY',
     '2' : 'TAC3',
     '7' : 'TAC3',
    '11' : 'TAC3',
    '12' : 'TAC3',
    '13' : 'TAC3',
    '16' : 'TAC3',
     '8' : 'TAC3/LHX8',
     '3' : 'PTHLH/PVALB',
     '4' : 'PTHLH/PVALB',
     '5' : 'PTHLH/PVALB',
     '6' : 'PTHLH/PVALB',
     '7' : 'PTHLH/PVALB',
}

hmba_cells_selected = adata_norm.obs.loc[adata_norm.obs.experiment == "marm_hmba", ["leiden"]]
hmba_cells_selected['clust_anno'] = [leiden_dict[x] for x in hmba_cells_selected['leiden']]
hmba_cells_selected.to_csv(Path(pdir, "marm_hmba/data/cells_annotated.csv"))
