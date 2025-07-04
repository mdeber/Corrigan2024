# Script to be used following Seurat_to_MatrixMarket.R.
# 
# Args: 
# 1) the MatrixMarket output directory created in that script
# 2) the output directory for h5ad files (will be created if doesn't exist). 
#
# This script is a work in progress in terms of the number of Seurat elements
# that are preseved in the output h5ad file.

print("##### Begin script MatrixMarket_to_AnnData.py #####")

import anndata as ad
import pandas as pd
from pathlib import Path
import sys
import os
from scipy.io import mmread
from scipy.sparse import csc_matrix

mm_dir = Path(sys.argv[1])
out_dir = Path(sys.argv[2])

if not mm_dir.is_dir():
    sys.exit("MatrixMarket directory (first argument) not found.")

if not out_dir.is_dir():
    print(f"Creating output directory {out_dir.name}...")
    os.mkdir(out_dir)
else:
    print(f"Found existing output directory {out_dir.name}.")

# Import cell metadata
path_cell_meta = Path(mm_dir, "cell_metadata.csv")
if not path_cell_meta.is_file():
    sys.exit("File 'cell_metadata.csv' not found in MatrixMarket directory.")
print("Importing 'cell_metadata.csv'...")
cell_meta = pd.read_csv(path_cell_meta, index_col = 0)

# When writing out, if any columns in adata.obs have dtype 'object', an error occurs
for c in cell_meta.columns:
    if cell_meta[c].dtype == "object":
        cell_meta[c] = cell_meta[c].astype("category")

# Import cell embeddings (to be added to anndata object later)
embed_dir = Path(mm_dir, "cell_embeddings")
if embed_dir.is_dir():
    print("Founding cell_embeddings directory. Importing...")
    embedding_files = [f for f in os.scandir(embed_dir) if f.is_file() and f.name.endswith('csv')]
    embeddings = {f.name.split('.csv')[0]:pd.read_csv(f.path, index_col = 0) for f in embedding_files}
    print(f"Imported embeddings: {list(embeddings.keys())}.")

# For each Seurat assay, import matrix market and feature metadata
assays = [f for f in os.scandir(mm_dir) if f.is_dir() and f.name != 'cell_embeddings']
for assay in assays:
    print(f"Processing assay '{assay.name}'...")
    
    # Feature metadata
    path_feat_meta = Path(assay.path, "feature_metadata.csv")
    if not path_feat_meta.is_file():
        print(f"Could not find 'feature_metadata.csv' file in assay '{assay.name}'. Skipping...")
        continue

    print(f"Importing 'feature_metadata.csv' for assay '{assay.name}'...")
    feat_meta = pd.read_csv(path_feat_meta, index_col = 0)
    # -> note, the feature metadata has 1-based indexing, not 0-based...

    # Raw counts (Seurat 'counts')
    path_counts = Path(assay.path, "counts.mtx")
    if not path_counts.is_file():
        print(f"Could not find 'counts.mtx' file in assay '{assay.name}'. Skipping...")
        continue

    print(f"Importing 'counts.mtx' for assay '{assay.name}'...")
    X_raw = csc_matrix(mmread(path_counts)).transpose() # csr_matrix

    # Scaled counts (Seurat 'data')
    path_data = Path(assay.path, "data.mtx")
    if path_data.is_file():
        print(f"Found 'data.mtx' file in assay '{assay.name}'. Importing...")
        X_scaled = csc_matrix(mmread(path_data)).transpose() # csr_matrix

        print(f"Creating AnnData for assay '{assay.name}'...")
        adata = ad.AnnData(X_scaled, obs = cell_meta, var = feat_meta)
        adata.raw = ad.AnnData(X_raw, obs = cell_meta, var = feat_meta)
    else:
        print(f"Creating AnnData for assay '{assay.name}' (raw counts only)...")
        adata = ad.AnnData(X_raw, obs = cell_meta, var = feat_meta)
    
    print("Adding embeddings to AnnData object...")
    for embed in embeddings:
        adata.obsm[embed] = embeddings[embed]

    print(f"Exporting h5ad for assay '{assay.name}'...")
    path_out = Path(out_dir, f"{assay.name}.h5ad")
    adata.write(path_out, compression="gzip")

print("Script finished.")
