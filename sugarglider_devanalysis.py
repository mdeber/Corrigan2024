#!/usr/bin/env python
# coding: utf-8

#Preprocessing by region
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import anndata
import os
import re
import numpy as np
import scipy
import seaborn
import bbknn
import matplotlib
import matplotlib.pyplot as plt
import scrublet as scr
import seaborn as sns

sc.settings.verbosity = 0     

import matplotlib as mpl
mpl.rcParams['figure.facecolor'] = 'white'

path = 'pathtosugarglider/extended'

striatum1 = sc.read_10x_h5(path + '/striatum1/cellbender_output_filtered.h5')
striatum2 = sc.read_10x_h5(path + '/striatum2/cellbender_output_filtered.h5')
striatum3 = sc.read_10x_h5(path + '/striatum3/cellbender_output_filtered.h5')
striatum4 = sc.read_10x_h5(path + '/striatum4/cellbender_output_filtered.h5')
striatum5 = sc.read_10x_h5(path + '/striatum5/cellbender_output_filtered.h5')
striatum6 = sc.read_10x_h5(path + '/striatum6/cellbender_output_filtered.h5')
striatum7 = sc.read_10x_h5(path + '/striatum7/cellbender_output_filtered.h5')
striatum8 = sc.read_10x_h5(path + '/striatum8/cellbender_output_filtered.h5')

def make_indices_unique(adata):
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
  
    return adata

    #Apply the function to each AnnData object

striatum1 = make_indices_unique(striatum1)
striatum2 = make_indices_unique(striatum2)
striatum3 = make_indices_unique(striatum3)
striatum4 = make_indices_unique(striatum4)
striatum5 = make_indices_unique(striatum5)
striatum6 = make_indices_unique(striatum6)
striatum7 = make_indices_unique(striatum7)
striatum8 = make_indices_unique(striatum8)

adata = anndata.AnnData.concatenate(striatum1, striatum2, striatum3, 
                                    striatum4, striatum5, striatum6, striatum7, striatum8,  join='outer', batch_categories=['striatum1', 'striatum2', 'striatum3', 
                                    'striatum4', 'striatum5', 'striatum6', 'striatum7', 'striatum8'])
adata.obs_names_make_unique()
adata.var_names_make_unique()

Results_file = ('pathtosugarglider/striatum.h5ad')

#h5ad stores whole anndata data structure
adata.write(Results_file)

#Define string to append to saved figures
save_name = "sugarglider"


# ## Preprocessing of data
#Basic filtering

sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=1)

# Annotate varibles for mitochondrial and ribosomal genes

mito_genes=[name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or 'MT-' in name]
ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]

# for each cell compute fraction of counts in mito genes vs. all genes

adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
adata.obs['percent_ribo'] = np.sum(
adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)

# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1)

sc.pl.violin(adata,keys=['percent_ribo', 'percent_mito','n_counts', 'n_genes'], groupby='batch', rotation=90, multi_panel=True)

sc.pl.scatter(adata, x='n_counts', y='n_genes')

#Histograms of the distrubution of UMI counts and gene numbers

Hist2 = seaborn.distplot(adata.obs['n_genes'], kde=False)
Hist2.set_xlabel("Number of genes", fontsize=12)
Hist2.set_ylabel("Frequency", fontsize=12)
#Hist2.axvline(1000, 0,1, color='red')

plt.show()

Hist1 = seaborn.distplot(adata.obs['n_counts'], kde=False)
Hist1.set_xlabel("Count depth", fontsize=12)
Hist1.set_ylabel("Frequency", fontsize=12)

plt.show()

Hist3 = seaborn.distplot(adata.obs['n_counts'][adata.obs['n_counts']<4000], kde=False, bins=60)
Hist3.set_xlabel("Count depth", fontsize=12)
Hist3.set_ylabel("Frequency", fontsize=12)
#Hist3.axvline(1000, 0,1, color='red')

plt.show()

#Filter cells and based on mitochondrial reads
adata = adata[adata.obs.percent_ribo < .2, :]
adata = adata[adata.obs.n_genes > 1000, :]
adata = adata[adata.obs.n_genes < 6000, :]

# ## Review data after filtering

sc.pl.violin(adata,keys=['percent_ribo', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True)

#Histograms of the distrubution of UMI counts and gene numbers

Hist2 = seaborn.distplot(adata.obs['n_genes'], kde=False)
Hist2.set_xlabel("Number of genes", fontsize=12)
Hist2.set_ylabel("Frequency", fontsize=12)
#Hist2.axvline(500, 0,1, color='red')

plt.show()

Hist1 = seaborn.distplot(adata.obs['n_counts'], kde=False)
Hist1.set_xlabel("Count depth", fontsize=12)
Hist1.set_ylabel("Frequency", fontsize=12)

plt.show()

Hist3 = seaborn.distplot(adata.obs['n_counts'][adata.obs['n_counts']<4000], kde=False, bins=60)
Hist3.set_xlabel("Count depth", fontsize=12)
Hist3.set_ylabel("Frequency", fontsize=12)
Hist3.axvline(1000, 0,1, color='red')

plt.show()

# ## Identify doublets

# Estimate and predict doublets using scrublet and the expected doublet rate in the relevant 10X Genomics protocol.
sce.pp.scrublet(adata, expected_doublet_rate = 0.08, n_prin_comps = 50, verbose = True)

sce.pl.scrublet_score_distribution(adata)

# Add column to AnnData.obs with categorical Singlet/Doublet instead of boolean True/False.
conditions = [
    (adata.obs["predicted_doublet"] == True),
    (adata.obs["predicted_doublet"] == False)]

values =['Doublet', 'Singlet']

adata.obs['doublet_info'] = np.select(conditions, values)

sc.pl.violin(adata, 'doublet_score', groupby = 'doublet_info', rotation=45)

# Remove doublets
adata = adata[adata.obs['predicted_doublet'] == False]
adata.write('pathtosugarglider/extended/striatumscRNAseq_prenorm.h5ad')


# ## Combining cortical and striatal samples from updated sugar glider genome and analyzing together
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import anndata
import os
import re
import numpy as np
import scipy
import seaborn
import bbknn
import matplotlib
import matplotlib.pyplot as plt
import scrublet as scr
import seaborn as sns

sc.settings.verbosity = 0     

import matplotlib as mpl
mpl.rcParams['figure.facecolor'] = 'white'

#set a path to your working directory
Results_file=('pathtosugarglider/extended/sugarglider_concatupdatedrevisions.h5ad')

path = 'pathtosugarglider/extended/'

Cortex = sc.read_h5ad(path + 'cortexscRNAseq_prenorm.h5ad')
Striatum = sc.read_h5ad(path + 'striatumscRNAseq_prenorm.h5ad')

adata = anndata.AnnData.concatenate(Cortex, Striatum, join='outer', batch_categories=['Cortex', 'Striatum'])
adata.write(Results_file)

sc.pl.violin(adata,keys=['percent_ribo','percent_mito', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True)
adata.raw = adata


adata.layers["counts"] = adata.X.copy()
print(adata.layers["counts"])

adata.write(Results_file)


# # Normalize data and scale  
adata = sc.read_h5ad(Results_file)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
sc.pp.scale(adata, max_value=10)
adata.write(Results_file)


# ## Calculate PCA 
adata = sc.read_h5ad(Results_file)

#Perform PCA dimensional reduction for the AnnData object.
sc.tl.pca(adata, svd_solver='arpack')

#Inspect the contribution of each PC to the variance of the data.
plt.rcParams["figure.figsize"] = (5, 5)
sc.pl.pca_variance_ratio(adata, log=True)

#Visualize the PCA loadings.
plt.rcParams["figure.figsize"] = (10, 5) 
sc.pl.pca_loadings(adata, components = [1, 2, 3, 4])

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["NSG2", 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)
adata.write(Results_file)



# ## Calculate UMAP
adata = sc.read_h5ad(Results_file)

#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

#Save the AnnData object as an .h5ad file.
adata.write(Results_file)


# # Plot UMAP and markers
adata = sc.read_h5ad(Results_file)

# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Integrated_pig" + str(resolution), frameon = False, legend_loc = "on data")

sc.pl.umap(adata, color = ['batch'], title = "by_batch")

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['NEUROD2','SLA', 'SLC17A6', 'FOXG1', 'NKX2-1', 'LHX6', 'DLX1', 'DLX2', 'GAD1', 'GAD2','COL19A1', 'CHRNA3', 'TH',  'PTHLH', 'ETV1','MAF', 'TAC3', 'STXBP6','ZIC1','NXPH2', 'LHX8', 'RBP4', 'NR2F2', 'PROX1', 'CHAT', 'SST', 'NPY'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['TBR1', 'PAX6', 'OLIG1', 'OLIG2', 'HES1'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


adata.write(Results_file)


# # Isolate only inhibitory neurons
adata = sc.read_h5ad(Results_file)

inhibmat = adata[:, ['DLX1', 'DLX2', 'DLX5', 'DLX6', 'GAD1', 'GAD2']].to_df()
inhibmat['leiden'] = adata.obs['leiden']
meanmat = inhibmat.groupby('leiden').mean()
print(meanmat>meanmat.mean(0))
boolmat = (meanmat > meanmat.mean(0)).sum(1)

print(boolmat.index[boolmat>=3])

#Filter data
adata = adata[boolmat[adata.obs['leiden']].values >= 3]

# Print information about the filtering results
print(f"Number of cells after filtering: {len(adata)}")

#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()


print(adata.layers["counts"])


#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

adata.layers["counts"]

processed_counts = adata.X
print(processed_counts)

#Inspect the contribution of each PC to the variance of the data.
plt.rcParams["figure.figsize"] = (5, 5)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

#Visualize the PCA loadings.
plt.rcParams["figure.figsize"] = (10, 5) 
sc.pl.pca_loadings(adata, components = [1, 2, 3, 4])

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["GRIA1", 'TPX2', 'CENPE', 'NRG3', 'ALCAM', 'DLX5', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

#save inhibitory neurons as a separate file
adata.write('pathtosugarglider/extended/concat_inhibneurons_updatedrevisions.h5ad')


# # Load and plot inhibitory neurons and remove any excitatory contamination

import scanpy as sc
import scanpy.external as sce
import pandas as pd
import anndata
import os
import re
import numpy as np
import scipy
import seaborn
import bbknn
import matplotlib
import matplotlib.pyplot as plt
import scrublet as scr
import seaborn as sns

adata = sc.read_h5ad('pathtosugarglider/extended/concat_inhibneurons_updatedrevisions.h5ad')

adata.obs['region'] = np.where(
    adata.obs['batch'].str.contains('cortex', case=False, na=False), 'cortex',
    np.where(
        adata.obs['batch'].str.contains('striatum', case=False, na=False), 'striatum',
        'unknown'  
    )
)

# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Sugar_glider_inhibitory neurons", frameon = False, legend_loc = "on data", save = 'leidenclusters.svg')

sc.pl.umap(adata, color= 'region', title = "SugarGlider_Inhibitory_Neurons", frameon = False)

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")


plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['SLA', 'SLC17A6', 'NEUROD2', 'NEUROG2', 'NEUROD6', 'EOMES'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

excitatorymat = adata[:, ['SLA', 'SLC17A6', 'NEUROD2', 'NEUROG2', 'NEUROD6', 'EOMES']].to_df()
excitatorymat['leiden'] = adata.obs['leiden']
meanmat2 = excitatorymat.groupby('leiden').mean()
print(meanmat2>meanmat2.mean(0))
boolmat2 = (meanmat2 > meanmat2.mean(0)).sum(1)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['ETV1','MAF','TAC3', 'LHX8','TH','CHRNA3', 'CHRNA5', 'CHRNA4', 'TRHDE', 'KIT', 'COL19A1', 'ZIC1', 'ZIC2','ZIC4', 'CHAT', 'batch'],use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

#Filter data -- removing excitatory clusters
adata = adata[boolmat2[adata.obs['leiden']].values <= 3]

# Print information about the filtering results
print(f"Number of cells after filtering: {len(adata)}")


#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()

#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["GRIA1", 'TPX2', 'CENPE', 'NRG3', 'ALCAM', 'DLX5', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

#save inhibitory neurons as a separate file
adata.write('pathtosugarglider/extended/concat_inhibneurons_updatedrevisions.h5ad')


# # Load and plot cleaned up inhibitory neurons without excitatory contamination
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import anndata
import os
import re
import numpy as np
import scipy
import seaborn
import bbknn
import matplotlib
import matplotlib.pyplot as plt
import scrublet as scr
import seaborn as sns

adata = sc.read_h5ad('pathtosugarglider/extended/concat_inhibneurons_updatedrevisions.h5ad')

# Cluster umap embeddings using leiden and save umap plots.
resolution = 0.5 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Sugar_glider" + str(resolution), frameon = False, legend_loc = "on data")

adata.obs["GEMwell"] = adata.obs.index.str.split("-").str[-2]

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes', 'batch', 'GEMwell'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")



# # Remove low quality clusters and recluster
adata = adata[adata.obs['leiden'].isin(['16', '14']) == False]

#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()

#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["GRIA1", 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

#save inhibitory neurons as a separate file
adata.write('pathtosugarglider/scRNAseq_sugarglider_concat_inhibneurons_updatedrevisions.h5ad')


# # Classify clusters and examine gene expression

import scanpy as sc
import scanpy.external as sce
import pandas as pd
import anndata
import os
import re
import numpy as np
import scipy
import seaborn
import bbknn
import matplotlib
import matplotlib.pyplot as plt
import scrublet as scr
import seaborn as sns

# === Master Figure Setup Script ===

import matplotlib
import matplotlib.pyplot as plt

# === 1. Use vector-safe font (Arial) ===
matplotlib.rcParams['font.family'] = 'Arial'

# === 2. Set all font sizes to Nature-compliant 5 pt ===
matplotlib.rcParams.update({
    'font.size': 5,
    'axes.titlesize': 5,
    'axes.labelsize': 5,
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'legend.fontsize': 5,
    'figure.titlesize': 5,
    'pdf.fonttype': 42,   # Keep fonts editable in PDF (if used)
    'svg.fonttype': 'none'  # Keep fonts as text in SVG (not outlines)
})

# === 3. Force vector output only (disable rasterization) ===
def make_figure_vector_safe(fig=None):
    """Turn off rasterization for all elements in a figure."""
    fig = fig or plt.gcf()
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(False)
        for im in ax.images:
            im.set_rasterized(False)

# === 4. Helper function to save SVG in correct format and size ===
def save_nature_figure(filename, width=3.35, height=3.35):
    """
    Save figure as fully vector SVG with font size fixed for Nature.
    Width/height in inches. Defaults to 3.35 in = 85 mm = single column width.
    """
    fig = plt.gcf()
    fig.set_size_inches(width, height)
    make_figure_vector_safe(fig)
    fig.savefig(filename, format='svg', bbox_inches='tight')
    plt.close(fig)

adata= sc.read_h5ad('pathtosugarglider/scRNAseq_sugarglider_concat_inhibneurons_updatedrevisions.h5ad')

# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Integrated_sugarglider" + str(resolution), frameon = False, legend_loc = "on data")

#subcluster 
sc.tl.leiden(adata, restrict_to=('leiden', ['18']), resolution=0.1, key_added='leiden_subcluster_CRABPY')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY'], title = "sugarglider_subcluster" + str(resolution), frameon = False, legend_loc = "on data")

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['CHRNA3', 'TAC3', 'PTHLH', 'MAF', 'ETV1', 'LHX8', 'TH', 'CHAT','NPY1R', 'ISL1', 'PENK', 'TSHZ1', 'SCGN', 'NR2F2', 'PROX1'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['NEUROD2','SLA', 'SLC17A6', 'NEUROD6', 'NEUROG2', 'EOMES', 'GLI3','FOXG1', 'NKX2-1', 'LHX6', 'DLX1', 'DLX2', 'GAD1', 'ETV1','COL19A1', 'CHRNA3', 'TH', 'PTHLH','MAF', 'TAC3', 'STXBP6','ZIC1','NXPH2', 'LHX8', 'RBP4', 'NR2F2', 'PROX1', 'CHAT', 'SST', 'NPY', 'CALB1', 'TH'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['HES1', 'PAX6', 'MEIS2', 'FOXP1', 'FOXP2', 'PENK', 'TSHZ1','ISL1', 'DRD1', 'DRD2', 'SCGN', 'NPY1R', 'PROX1', 'NR2F2', 'ZIC1', 'ZIC2', 'ZIC4', 'CHAT', 'GBX1', 'LAMP5', 'CCK', 'RELN'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

#subcluster 
sc.tl.leiden(adata, restrict_to=('leiden_subcluster_CRABPY', ['5']), resolution=0.2, key_added='leiden_subcluster_CRABPY2')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY2'], use_raw=False, frameon = False, color_map = "PuRd",legend_loc = "on data")

#subcluster 
sc.tl.leiden(adata, restrict_to=('leiden_subcluster_CRABPY2', ['10']), resolution=0.2, key_added='leiden_subcluster_CRABPY3')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY3'], use_raw=False, frameon = False, color_map = "PuRd",legend_loc = "on data")

# Define a dictionary to map old cluster names to new names
cluster_name_mapping = {
    
    '7' : 'Progenitor',
    '11' : 'Progenitor',
    '21' : 'Progenitor',
    '19' : 'Progenitor',
    
    
    '18,1' : 'MGE_CRABP1/TAC3',    
    '18,0' : 'MGE_CRABP1/MAF',


    
    '2' : 'MGE_LHX6/MAF',
    '0' : 'MGE_LHX6/MAF',
 
    
    '5,1' : 'LGE_MEIS2/PAX6',
    '13' : 'LGE_MEIS2/PAX6',
    '5,0': 'LGE_MEIS2/PAX6',
    '10,0': 'LGE_MEIS2/PAX6',
    '10,2': 'LGE_MEIS2/PAX6',
    '10,1': 'LGE_MEIS2/PAX6',

    '8' : 'LGE_FOXP1/PENK',
    '1' : 'LGE_FOXP1/PENK',
    '6' : 'LGE_FOXP1/PENK',
    '16' : 'LGE_FOXP1/PENK',
    '22' : 'LGE_FOXP1/PENK',
    
    '17' : 'LGE_FOXP1/ISL1',
    '14' : 'LGE_FOXP1/ISL1',
    '9' : 'LGE_FOXP1/ISL1',
    '3' : 'LGE_FOXP1/ISL1',
    '24' : 'LGE_FOXP1/ISL1',
    
    
    '15' : 'CGE_NR2F2/PROX1',

    
    '12' : 'LGE_FOXP2/TSHZ1',
    '4' : 'LGE_FOXP2/TSHZ1',
    '9,2' : 'LGE_FOXP2/TSHZ1',
    '20' : 'LGE_FOXP2/TSHZ1',





    
   '23': 'VMF_CRABP1/LHX8',

    
 
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['leiden_subcluster_CRABPY3'] = adata.obs['leiden_subcluster_CRABPY3'].replace(cluster_name_mapping)


paldict={'CGE_NR2F2/PROX1': 'mediumpurple',
    'LGE_FOXP1/ISL1': 'royalblue',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_FOXP2/TSHZ1': '#17344c',
    'LGE_MEIS2/PAX6': 'orangered',
    'LGE_MEIS2/PAX6/SCGN': 'orange',
    'LGE-OB_MEIS2/PAX6': 'red',
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TAC3': 'fuchsia',
    'MGE_LHX6/MAF': 'skyblue',
    'MGE_LHX6/NPY': 'teal',
    'VMF_ZIC1/ZIC2': 'green',
    'Progenitor' : 'pink',
    'VMF_CRABP1/LHX8':'mediumspringgreen',

}

# Rename 'cortex' to 'Cortex' and 'striatum' to 'Striatum' in adata.obs['region']
adata.obs['region'] = adata.obs['region'].replace({'cortex': 'Cortex', 'striatum': 'Striatum'})

# Check the updated 'region' column
print(adata.obs['region'].unique())

# Define the custom colors for Cortex and Striatum
region_colors = {'Cortex': '#1d26cf', 'Striatum': '#f76f3e'}

plt.rcParams["figure.figsize"] = (10, 10)
sc.set_figure_params(fontsize = 20)
plt.rcParams["figure.figsize"] = (10, 10)

# Plot UMAP with custom colors for the 'region' column
sc.pl.umap(adata, color='region', frameon=False, size=10, save='_region.svg', palette=region_colors)
sc._settings.settings._vector_friendly=False

#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata, groupby = 'leiden_subcluster_CRABPY3', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df

adata.obs['leiden_subcluster_CRABPY3'] = adata.obs['leiden_subcluster_CRABPY3'].cat.reorder_categories(['Progenitor', 'MGE_LHX6/MAF', 'MGE_CRABP1/MAF', 'MGE_CRABP1/TAC3',
                                                                                                       'LGE_MEIS2/PAX6', 'LGE_FOXP2/TSHZ1', 'LGE_FOXP1/ISL1', 'LGE_FOXP1/PENK', 'CGE_NR2F2/PROX1',
                                                                                                      'VMF_CRABP1/LHX8'])


#Matrix of marker genes for inhibitory neuron clusters
sc.set_figure_params(fontsize = 15)
#marker genes based on Schmitz et al. 2022
marker_genes_dict = {
    'Progenitor' : ['ASPM', 'CENPF', 'CKB'],
    'MGE_LHX6/MAF': ['LHX6', 'SOX6', 'MEF2C'],
    'MGE_CRABP1/MAF': ['MAF', 'MAFB','ETV1', 'COL19A1', 'CNR1', 'MTUS2',  'KIT', ],
    'MGE_CRABP1/TAC3': ['NRTN','TRHDE','STXBP6','TAC3', 'CHRNA3'],
    'LGE_MEIS2/PAX6': ['PAX6','MEIS2'],
    'LGE_FOXP2/TSHZ1': ['EYA2','FOXP2', 'FOXP4','TSHZ1'],
    'LGE_FOXP1/ISL1': ['PBX3', 'ISL1'],
    'LGE_FOXP1/PENK': ['FOXP1', 'SIX3', 'PENK'],
    'CGE_NR2F2/PROX1': ['PDZRN3', 'NR3C2', 'NR2F2', 'PROX1', 'NPAS1'],
    'VMF_CRABP1/LHX8': ['ZIC4','ZIC1','CHAT'],

}

sc.pl.matrixplot(adata, groupby='leiden_subcluster_CRABPY3', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes=True,save='matrix.svg')

# Compute proportions of Cortex vs. Striatum for each cluster
cluster_region_counts = adata.obs.groupby(['leiden_subcluster_CRABPY3', 'region']).size().unstack(fill_value=0)
cluster_region_props = cluster_region_counts.div(cluster_region_counts.sum(axis=1), axis=0)

# Ensure order matches matrixplot
cluster_order = adata.obs['leiden_subcluster_CRABPY3'].cat.categories
cluster_region_props = cluster_region_props.loc[cluster_order]


# Set up figure with two subplots: one for barplot, one for matrixplot
fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 10]}, figsize=(25, 10))

# Top subplot: stacked barplot for cortex vs. striatum proportions
cluster_region_props.plot(kind='bar', stacked=True, color=['#1d26cf', '#f76f3e'], ax=ax[0])
# Remove the grid lines
ax[0].grid(False)
# Format barplot
ax[0].set_xticks([])
ax[0].set_ylabel('Proportion')
ax[0].set_title('Cortex vs. Striatum Proportions per Cluster')
ax[0].legend(
    title="Region", 
    labels=['Cortex', 'Striatum'], 
    loc='center left', 
    bbox_to_anchor=(1, 0.5)  # Moves legend to the right
)
         
# Bottom subplot: matrixplot
sc.pl.matrixplot(
    adata,
    groupby='leiden_subcluster_CRABPY3',
    var_names=marker_genes_dict,
    use_raw=False,
    vmin=-2,
    vmax=2,
    cmap='PiYG_r',
    swap_axes=False,
    show=False,
    ax=ax[1]  # Draw on second subplot
)

plt.savefig("proportionmatrix.svg")
plt.tight_layout()
plt.show()

# Filter the AnnData object 
clusters_of_interest = adata.obs['leiden_subcluster_CRABPY3'].isin(['MGE_CRABP1/MAF','MGE_CRABP1/TAC3'])
adata_filtered = adata[clusters_of_interest]

# List of genes of interest
genes_of_interest = ['LHX6', 'ANGPT2', 'MAF', 'CHL1','RBP4', 'ARX', 'TAC3', 'TRH','TRHDE', 'CHRNA3','ZIC1', 'ZIC2', 'ZIC4','LHX8', 'TH', 'PTHLH',  'PARD3', 'OLFM2', 'PDE1C']


# Create a dot plot for the specified genes
sc.pl.dotplot(adata_filtered, var_names=genes_of_interest, use_raw= False, vmin=-2, vmax=2, swap_axes = True, cmap='PiYG_r',groupby='leiden_subcluster_CRABPY3', save = 'macdotplot.svg')


#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata_filtered, groupby = 'region', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata_filtered.uns['rank_genes_groups']['names']).head(20)
pd.set_option('display.max_columns', 500)
small_degs_df

# Filter the AnnData object to include only clusters that start with '0840', '0841', or '0842'
clusters_of_interest = adata.obs['leiden_subcluster_CRABPY3'].isin(['MGE_CRABP1/TAC3'])
adata_filtered = adata[clusters_of_interest]

marker_genes_dict = ['LHX6', 'NKX2-1', 'SOX6', 'MAF', 'CHL1','ETV1','GALNTL6','OLFM2','NRTN', 'KIT','TRHDE','STXBP6','TAC3', 'CHRNA3',  'PARD3',  'PDE1C', 'ZIC1', 'ZIC2', 'ZIC4', 'LHX8', 'CXCR4', 'UNC5C',  'DBI', 'ADCY8', 'PCSK5']



sc.pl.dotplot(adata_filtered, groupby='region', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes=True,  figsize=(2.5,7),save='matrix.svg')

# Get counts of cells per region
region_counts = adata_filtered.obs['region'].value_counts()

# Copy the 'leiden_subcluster_CRABPY2' observation to a new 'class' observation
adata.obs['class'] = adata.obs['leiden_subcluster_CRABPY3']



# Rename specific columns in obs
adata.obs.rename(columns={'leiden_subcluster_CRABPY': 'leiden_subcluster','leiden_subcluster_CRABPY2': 'leiden_subcluster2', 'leiden_subcluster_CRABPY3': 'leiden_subcluster3'}, inplace=True)

# Extract the cell barcodes and additional obs columns
obs = adata.obs

# Extract the gene information and additional var columns
var = adata.var

if 'counts' in adata.layers:
    # Extract the count matrix from the 'counts' layer
    matrix = adata.layers['counts']
else:
    raise ValueError("The 'counts' layer does not exist in the provided AnnData object.")

# If matrix is a sparse matrix, convert it to CSR format
if not scipy.sparse.isspmatrix_csr(matrix):
    matrix = scipy.sparse.csr_matrix(matrix)

# Save obs.tsv
obs.to_csv('pathtogeo/sugarglider_barcodes.tsv', sep='\t', index=True)

# Save var.tsv
var.to_csv('pathtogeo/sugarglider_genes.tsv', sep='\t', index=True)

# Save matrix.mtx
scipy.io.mmwrite('pathtogeo/sugarglider_matrix.mtx', matrix)


adata.write('pathtosugarglider/extended/sugarglider_processed.h5ad')
adata = sc.read_h5ad('pathtosugarglider/extended/sugarglider_processed.h5ad')


# # Subset data for integration

# Assuming 'class' is stored in adata.obs['class']
class_column = "class"
max_cells_per_class = 1000

# Create an empty list to store selected cells
selected_indices = []

# Group by class and randomly select up to 1000 cells per class
for cell_class, group in adata.obs.groupby(class_column):
    selected_indices.extend(np.random.choice(group.index, size=min(max_cells_per_class, len(group)), replace=False))

# Subset the AnnData object
adata_subset = adata[selected_indices].copy()

# Verify the subset
print(adata_subset)
print(adata_subset.raw.X)

adata = adata_subset


#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()

adata.X = np.array(adata.X)  # Convert to numpy array
adata.write('pathtosugarglider/extended/sugarglider_subsetraw1k.h5ad')

