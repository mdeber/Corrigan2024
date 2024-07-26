"""
Script to run hotspot on an h5ad file. Everything of substance follows 
the hotspot tutorial workflow: https://hotspot.readthedocs.io/en/latest/

Args:
1) Path to h5ad input file
2) Path to an output directory for several output files, including the gene 
   modules. Folder will be created if it doesn't already exist. 

Note: If SAVE_LOCAL_CORRELATIONS = True, script outputs a file, 
'local_correlations.csv.gz', which is typically large and not necessary for 
anything, but takes a very long time to compute.
"""
## Options:
SAVE_LOCAL_CORRELATIONS = False
# NCORES=12 # error when passing to jobs arg in hotspot functions -> skip parallelization

print("##### Begin script run_hotspot.py #####")

import anndata as ad
import hotspot
from scipy.sparse import csc_matrix

from pathlib import Path
import sys
import os

in_path = Path(sys.argv[1])
out_dir = Path(sys.argv[2])

if not in_path.exists():
    sys.exit("Input h5ad file (first argument) not found.")

if not out_dir.is_dir():
    print(f"Creating output directory {out_dir.name}...")
    os.mkdir(out_dir)
else:
    print(f"Found existing output directory {out_dir.name}.")

print("Reading in h5ad file...")
adata = ad.read_h5ad(in_path)

print("Creating hotspot object...")

# verify that genes in adata.X are same as adata.raw.X
# (structuring this 'if statement' like this in case adata.raw.var doesn't exist)
if adata.raw.X.shape != adata.X.shape:
    print("Subsetting genes in 'adata.raw.X' to those found in 'adata.X'...")
    valid_genes_norm = adata.raw.var.index.isin(adata.var.index)
else:
    valid_genes_norm = adata.var.index.isin(adata.var.index) # (just array of True)

# (Hotspot will throw an error if any genes have 0 counts across all cells)
valid_genes_nonzero = adata.raw.X.sum(0).A1 > 0

valid_genes = valid_genes_norm & valid_genes_nonzero
if sum(valid_genes_nonzero) < adata.raw.shape[1]:
    # are any *additional* genes being removed for having 0 counts?
    if sum(valid_genes_norm == valid_genes) != len(valid_genes):        
        print("Removing genes with all 0 counts...")

if not 'X_pca' in adata.obsm.keys():
    if 'pca' in adata.obsm.keys():
        adata.obsm['X_pca'] = adata.obsm['pca']
    elif 'PCA' in adata.obsm.keys():
        adata.obsm['X_pca'] = adata.obsm['PCA']
    else:
        print("Did not find 'X_pca', 'pca' or 'PCA' in obsm. Will try running PCA...")
        import scanpy as sc
        sc.tl.pca(adata)

if sum(valid_genes) < adata.raw.shape[1]:
    # Subset genes and save everything needed from adata
    cell_meta = adata.obs
    gene_meta = adata.raw.var.iloc[valid_genes, :]
    X_raw = adata.raw.X[:, valid_genes]
    X_pca = adata.obsm['X_pca']
    # Re-assign adata
    adata = ad.AnnData(obs = cell_meta, var = gene_meta)
    adata.raw = ad.AnnData(X_raw, obs = cell_meta, var = gene_meta)
    adata.obsm['X_pca'] = X_pca
    del X_raw, X_pca, cell_meta, gene_meta

adata.obs['total_umi_counts'] = adata.raw.X.sum(1).A1
adata.layers['counts'] = csc_matrix(adata.raw.X) # csc significantly faster
hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model="danb",
    latent_obsm_key="X_pca",
    umi_counts_obs_key="total_umi_counts"
)

print("Creating the KNN graph & finding informative genes (by autocorrelation)...")
hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
hs_results = hs.compute_autocorrelations()

print("Evaluating pair-wise gene associations (local correlations). This may take a while...")
hs_genes = hs_results.loc[hs_results.FDR < 0.05].index
local_correlations = hs.compute_local_correlations(hs_genes) # parallelization fails

# (note: local_correlations are also stored in hs.local_correlations_z)
if SAVE_LOCAL_CORRELATIONS:
    print("Saving local_correlations.csv.gz...")
    local_correlations.to_csv(
        Path(out_dir, "local_correlations.csv.gz"),
        compression={'method': 'gzip'}
    )

# Group genes into modules (convenience function doing agglomerative clustering);
# unassigned genes have module of -1
# (note: min_gene_threshold prevents merging of modules with >30 genes)
print("Creating and saving gene modules...")
modules = hs.create_modules(
    min_gene_threshold=30, core_only=True, fdr_threshold=0.05
)
modules.to_csv(Path(out_dir, "modules.csv"))

## not using the hotspot module score calculations, so don't export them
# print("Calculating and saving module scores...")
# module_scores = hs.calculate_module_scores()
# module_scores.to_csv(Path(out_dir, "module_scores.csv"))

print("Script finished.")
