# This script takes in an h5ad file and exports matrix market files with
# csvs for metadata. It's a work in progress in terms of the number of 
# attributes it exports.
#
# Args:
# 1) Path to h5ad file
# 2) Output directory. If it doesn't exist, it will be created
# 3) [optional] 'counts_only': Whether to export only the counts and not the normalized counts as well.
#    Should be 'True' or 'False' (default is 'False').
# 
# Note it's missing a lot, and doesn't do anything with:
#   uns
#   varm
#   obsp
#
# This script also isn't really dealing with the issue of features being different in adata.var and adata.raw.var;
# it will just export the features from adata.raw if they're different
print("##### Begin script AnnData_to_MatrixMarket.py #####")
import anndata as ad
import pandas as pd
from pathlib import Path
import sys
import os
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

path_h5ad = Path(sys.argv[1])
out_dir = Path(sys.argv[2])

if len(sys.argv) == 4: # for 3 args
    counts_only = sys.argv[3]
    if counts_only in ['True', 'False']:
        counts_only = counts_only == 'True'
    else:
        sys.exit("Third arg given ('counts_only'), but not one of 'True' or 'False'")
else:
    counts_only = False

if not path_h5ad.is_file():
    sys.exit("h5ad file (first argument) not found.")

if not out_dir.is_dir():
    print(f"Creating output directory {out_dir.name}...")
    os.mkdir(out_dir)
else:
    print(f"Found existing output directory {out_dir.name}.")

adata = ad.read_h5ad(path_h5ad, backed = 'r')

# if nothing in adata.raw, assume adata.X is the counts
if adata.raw is None:
    counts_only = True

# Features/metadata
print("Exporting 'feature_metadata.csv'...")
path_feature_meta = Path(out_dir, "feature_metadata.csv")
if not adata.raw is None:
    if adata.raw.X.shape[1] != adata.var.shape[0]:
        print("(Using feature metadata from adata.raw)...")
        adata.raw.var.to_csv(path_feature_meta)
else:
    adata.var.to_csv(path_feature_meta)

# Cell metadata
print("Exporting 'cell_metadata.csv'...")
path_cell_meta = Path(out_dir, "cell_metadata.csv")
adata.obs.to_csv(path_cell_meta)

# Normalized data
if not counts_only:
    print("Exporting 'data.mtx'...")
    path_data = Path(out_dir, "data.mtx")
    mmwrite(path_data, csr_matrix(adata.X))

# Raw data [and note there's a problem in R if export from unsigned integer dtype]
# (and again assuming adata.X is the counts if nothing in adata.raw)
print("Exporting 'counts.mtx'...")
path_counts = Path(out_dir, "counts.mtx")
if not adata.raw is None:
    if (adata.raw.X.dtype == 'uint32'):
        mmwrite(path_counts, adata.raw.X.value.astype('int32'))
    elif (adata.raw.X.dtype == 'uint64'):
        mmwrite(path_counts, adata.raw.X.value.astype('int64'))
    else:
        mmwrite(path_counts, adata.raw.X.value)
else:
    mmwrite(path_counts, csr_matrix(adata.X.value))

# Embeddings
embed_dir = Path(out_dir, "cell_embeddings")
if not embed_dir.is_dir():
    print("Creating export directory for cell embeddings...")
    os.mkdir(embed_dir)
else:
    print("Found existing export directory for cell embeddings.")

embeddings = [x for x in adata.obsm.keys()]
for embed in embeddings:
    print(f"Exporting embeddings from obsm: {embed}...")
    path_embed = Path(embed_dir, f"{embed}.csv")
    # adata.obsm[embed].tofile(path_embed, sep = ',')
    pd.DataFrame(adata.obsm[embed]).to_csv(path_embed)

print("Script AnnData_to_MatrixMarket.py finished.\n")
