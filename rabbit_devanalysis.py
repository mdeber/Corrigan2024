#!/usr/bin/env python
# coding: utf-8

##Preprocessing by region
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
adata=sc.read_10x_h5('pathtorabbit/newgenome/striatum/RabbitE23striatum_cellbender_outputlr5_filtered.h5')
Results_file='pathtorabbit/newgenome/striatum/scRNAseq_rabbitE23.h5ad'
adata

adata.obs_names_make_unique()
adata.var_names_make_unique()

#h5ad stores whole anndata data structure
adata.write(Results_file)



#Define string to append to saved figures
save_name = "rabbit_striatum_E23"


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



#violin plots of the computed quality measures
sc.pl.violin(adata,keys=['percent_ribo','percent_mito', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True)


sc.pl.scatter(adata, x='n_counts', y='percent_mito')
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
adata = adata[adata.obs.percent_ribo < .15, :]
#adata = adata[adata.obs.percent_mito < .1, :]
adata = adata[adata.obs.n_genes < 4500, :]
adata = adata[adata.obs.n_genes > 750, :]


# ## Review data after filtering


sc.pl.violin(adata,keys=['percent_ribo','percent_mito', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True, save= save_name + 'QC_metrics')


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
sce.pp.scrublet(adata, expected_doublet_rate = 0.20, n_prin_comps = 50, verbose = True)

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

adata.write('pathtorabbit/scRNAseq_striatum_rabbit_prenorm.h5ad')



# ## Rabbit E23 concat analysis

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


path = '/pathtorabbit/'

Cortex = sc.read_h5ad(path + 'scRNAseq_cortex_rabbit_prenorm.h5ad')
Striatum = sc.read_h5ad(path + 'scRNAseq_striatum_rabbit_prenorm.h5ad')


#set a path to your working directory
Results_file=('/pathtorabbit/scRNAseq_rabbitE23_concat.h5ad')


for gem in [Cortex, Striatum]:
    print(gem.shape[0])

adata = anndata.AnnData.concatenate(Cortex, Striatum, join='outer', batch_categories=['Cortex', 'Striatum'])


adata.obs_names_make_unique()
adata.var_names_make_unique()

#Set the .raw attribute of the AnnData object
adata.raw = adata

#h5ad stores whole anndata data structure
adata.write(Results_file)


#Define string to append to saved figures
save_name = "rabbit_concat"


# ## Normalize counts scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

#Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

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
pca_var_genes = ["VIM", "MEF2C", "MYO16", "DLX6", "DLX5",
                 "GAD2"]
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ['percent_mito', 'percent_ribo', 'n_genes', 'n_counts']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


adata.write(Results_file)

#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)


#Save the AnnData object as an .h5ad file.
adata.write(Results_file)

# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Rabbit_concat" + str(resolution), frameon = False, legend_loc = "on data")

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','n_counts', 'n_genes'], wspace=0.25, ncols = 2, save = 'umapconcat_qc.png')


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['FOXG1', 'NKX2-1',  'MEF2C', 'DLX1','LHX6', 'CRABP1','ETV1','TAC3', 'CHRNA3', 'CHRNA4', 'CHRNA7','NPY', 'SST','GAD1','ZIC1','NXPH2', 'LHX8', 'RBP4', 'CHAT', 'NR2F2', 'PROX1', 'MEIS2', 'PAX6', 'SCGN', 'TSHZ1','FOXP1', 'PENK', 'ISL1', 'NPY1R'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['PAX6', 'TBR1', 'SATB2', 'CUX2', 'EOMES', 'HES1','OLIG2', 'EGFR', 'ERBB4'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


sc.pl.umap(adata, color = ['NEUROD2', 'DLX1', 'DLX2', 'GAD1', 'GAD2', 'SLA'], use_raw = False, frameon = False, color_map = "PuRd", size= 20, ncols = 3)


# ## Isolate inhibitory neurons

inhibmat = adata[:, ['DLX1','DLX2', 'DLX5', 'DLX6', 'GAD1', 'GAD2']].to_df()
inhibmat['leiden'] = adata.obs['leiden']


meanmat = inhibmat.groupby('leiden').mean()
print(meanmat>meanmat.mean(0))
boolmat = (meanmat > meanmat.mean(0)).sum(1)

print(boolmat.index[boolmat>=3])


#Filter data
adata = adata[boolmat[adata.obs['leiden']].values >= 3]

# Print information about the filtering results
print(f"Number of cells after filtering: {len(adata)}")

#Save the AnnData object as an .h5ad file.
adata.write('pathtorabbit/scRNAseq_rabbit_concat_inhibneurons.h5ad')


# # Renormalize, Recalculate PCA, recalculate UMAP
adata = sc.read_h5ad('pathtorabbit/scRNAseq_rabbit_concat_inhibneurons.h5ad')


#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()



# Directly access the raw count matrix - calling it "processed counts" because later .X will store your processed data, while the "counts" layer will 
processed_counts = adata.X
print(processed_counts)

adata.layers["counts"] = adata.X.copy() #Copy the raw data into a "layer" in the AnnData structure - this will be used for normalization, etc.

adata.layers["counts"] #This should resemble the raw data at this point

print(adata.X.shape)

#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

if isinstance(adata.X, np.matrix):
    adata.X = np.array(adata.X)  # Convert to numpy.ndarray

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

adata.layers["counts"]

processed_counts = adata.X
print(processed_counts)

#Inspect the contribution of each PC to the variance of the data.
plt.rcParams["figure.figsize"] = (5, 5)
sc.pl.pca_variance_ratio(adata, log=True)


#Visualize the PCA loadings.
plt.rcParams["figure.figsize"] = (10, 5) 
sc.pl.pca_loadings(adata, components = [1, 2, 3, 4])

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["CENPF", "TOP2A", "NYAP2",
                 "GAD1", "DCC", "CHN2",
                 "CUX2"]
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

# Convert adata.X if it's a numpy matrix
if isinstance(adata.X, np.matrix):
    adata.X = np.asarray(adata.X)

# Check and convert any layers
for key in adata.layers.keys():
    if isinstance(adata.layers[key], np.matrix):
        adata.layers[key] = np.asarray(adata.layers[key])


#Save the AnnData object as an .h5ad file.
adata.write('pathtorabbit/scRNAseq_rabbit_concat_inhibneurons.h5ad')


# # Load inhibitory neuron dataset plot and remove excitatory contamination

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

adata = sc.read_h5ad('pathtorabbit/scRNAseq_rabbit_concat_inhibneurons.h5ad')


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Rabbit_Inhibitory_Neurons_E23" + str(resolution), frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito','n_counts', 'n_genes'], wspace=0.25, ncols = 2)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['DLX1', 'DLX2', 'GAD1', 'GAD2', 'SLA', 'NEUROD2', 'SLC17A6'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['CRABP1','MAF','TAC3', 'LHX8','TH','CHRNA3', 'CHRNA5', 'CHRNA4', 'TRHDE', 'KIT', 'COL19A1', 'CHRNA7', 'ZIC1', 'ZIC2','ZIC4', 'CHAT'],use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NKX2-1', 'MEF2C',  'ETV1', 'MAF', 'MEIS2', 'PAX6', 'FOXP1', 'FOXP2', 'NR2F2', 'LHX8', 'DLX1', 'GAD1', 'HES1', 'SPC24', 'CDC20', 'TOP2A'], 
           use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 4)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NHLH1', 'NHLH2', 'TP73', 'NPY', 'SP8', 'PENK', 'ISL1', 'TSHZ1', 'MEIS2', 'PAX6', 'RELN', 'NPY1R', 'LHX8', 'CRABP1', 'ZIC1', 'ZIC2', 'CHAT'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

excitatorymat = adata[:, ['SLA', 'SLC17A6', 'NEUROD2', 'NEUROG2', 'NEUROD6', 'EOMES']].to_df()
excitatorymat['leiden'] = adata.obs['leiden']
meanmat2 = excitatorymat.groupby('leiden').mean()
print(meanmat2>meanmat2.mean(0))
boolmat2 = (meanmat2 > meanmat2.mean(0)).sum(1)

#Filter data
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
pca_var_genes = ["GRIA1", 'NRXN3', 'TPX2', 'CENPE', 'NRG3', 'ALCAM', 'DLX5', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

adata.write('pathtorabbit/scRNAseq_rabbit_concat_inhibneurons.h5ad')


# # Load plot and examine cleaned up data
adata = sc.read_h5ad('pathtorabbit/scRNAseq_rabbit_concat_inhibneurons.h5ad')


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Rabbit_Inhibitory_Neurons_E23" + str(resolution), frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')


# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['NEUROD2','SLA', 'SLC17A6', 'NEUROD6', 'NEUROG2', 'EOMES', 'GLI3','FOXG1', 'NKX2-1', 'LHX6', 'DLX1', 'DLX2', 'GAD1','CRABP1','COL19A1', 'CHRNA3', 'TH', 'PVALB', 'PTHLH', 'ETV1','MAF', 'TAC3', 'STXBP6','ZIC1','NXPH2', 'LHX8', 'RBP4', 'NR2F2', 'PROX1', 'CHAT', 'SST', 'NPY', 'CALB1', 'TH'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['HES1', 'HES5', 'PAX6', 'MEIS2', 'FOXP1', 'FOXP2', 'PENK', 'TSHZ1','ISL1', 'DRD1', 'DRD2', 'SCGN', 'NPY1R', 'PROX1', 'NR2F2', 'ZIC1', 'ZIC2', 'ZIC4', 'CHAT', 'GBX1', 'LAMP5', 'VIP'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


# # Remove low quality clusters and recluster
adata = adata[adata.obs['leiden'].isin(['13']) == False]


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
pca_var_genes = ["GRIA1", 'NRXN3', 'TPX2', 'CENPE', 'NRG3', 'ALCAM', 'DLX5', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

adata.write('pathtorabbit/scRNAseq_rabbit_concat_inhibneurons.h5ad')


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

adata = sc.read_h5ad('pathtorabbit/scRNAseq_rabbit_concat_inhibneurons.h5ad')

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





# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Integrated_rabbit" + str(resolution), frameon = False, legend_loc = "on data")

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")

#subcluster
sc.tl.leiden(adata, restrict_to=('leiden', ['1']), resolution=0.2, key_added='leiden_subcluster_CRABPY')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY'], title = "Rabbit_subcluster" + str(resolution), frameon = False, legend_loc = "on data")

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['LHX6', 'CRABP1', 'CHRNA3', 'TAC3', 'TH', 'LHX8', 'batch'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['LHX6', 'CRABP1', 'CHRNA3', 'TAC3', 'TH', 'LHX8'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['MKI67','HES1', 'HES5', 'PAX6', 'MEIS2', 'FOXP1', 'FOXP2', 'PENK', 'TSHZ1','ISL1', 'DRD1', 'DRD2', 'SCGN', 'NPY1R', 'PROX1', 'NR2F2', 'ZIC1', 'ZIC2', 'ZIC4', 'CHAT', 'GBX1', 'LAMP5', 'VIP'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata, groupby = 'leiden', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df


#subcluster 
sc.tl.leiden(adata, restrict_to=('leiden_subcluster_CRABPY', ['4']), resolution=0.2, key_added='leiden_subcluster_CRABPY2')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY2'], title = "Rabbit_subcluster" + str(resolution), frameon = False, legend_loc = "on data")


plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['MEIS2', 'PAX6', 'SCGN', 'TSHZ1', 'FOXP1', 'FOXP2','PENK', 'ISL1', 'NPY1R','PROX1', 'NR2F2', 'ZIC1', 'ZIC2', 'ZIC4', 'CHAT', 'GBX1', 'LAMP5', 'VIP'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

# Define a dictionary to map old cluster names to new names
cluster_name_mapping = {
 
    '1,2' : 'MGE_CRABP1/MAF',
    
    
    '16' : 'Progenitor',
    '13' : 'Progenitor',
    '5' : 'Progenitor',
    '14' : 'Progenitor',
    '7' : 'Progenitor',

    
    '1,1' : 'MGE_LHX6/MAF',
    '1,0' : 'MGE_LHX6/MAF',
    '0' : 'MGE_LHX6/MAF',
    '2' : 'MGE_LHX6/MAF',
    
    '18' : 'MGE_LHX6/NPY',
        
    '6' : 'CGE_NR2F2/PROX1',
    '3' : 'CGE_NR2F2/PROX1',
    
    
    '8' : 'LGE_MEIS2/PAX6',
    '10' : 'LGE_MEIS2/PAX6',



    '9' : 'LGE_FOXP1/PENK',
    '12' : 'LGE_FOXP1/ISL1',
    '4,1' : 'LGE_FOXP1/ISL1',
    '15' : 'LGE_FOXP1/ISL1/NPY1R',
    


    '17' : 'LGE_FOXP2/TSHZ1',
    '11' : 'LGE_FOXP2/TSHZ1',
    '4,0' : 'LGE_FOXP2/TSHZ1'


    
 
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['leiden_subcluster_CRABPY2'] = adata.obs['leiden_subcluster_CRABPY2'].replace(cluster_name_mapping)

paldict={'CGE_NR2F2/PROX1': 'mediumpurple',

    'LGE_FOXP1/ISL1': 'royalblue',
    'LGE_FOXP1/ISL1/NPY1R': 'blue',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_FOXP2/TSHZ1': '#17344c',
    'LGE_MEIS2/PAX6': 'orangered',
    'LGE_MEIS2/PAX6/SCGN': 'orange',

    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TAC3': 'fuchsia',
    'MGE_LHX6/MAF': 'skyblue',
    'MGE_LHX6/NPY': 'teal',
    'VMF_ZIC1/ZIC2': 'green',
    'Progenitor' : 'pink',

}
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= ['class'], palette=paldict, title = "Rabbit_inhibitory", frameon = False, size = 10, save = '_officialclusters.svg')

mge_cells = adata[adata.obs['class'] == 'MGE_CRABP1/MAF']

# Make sure the gene name matches exactly — case sensitive
gene = 'TAC3'

# Confirm gene exists
if gene in mge_cells.var_names:
    # Count how many cells express NKX2-1 (expression > 0)
    nkx2_1_positive = (mge_cells[:, gene].X > 0).sum()
    print(f"Number of MGE_CRABP1 cells expressing {gene}: {nkx2_1_positive}")
else:
    print(f"Gene {gene} not found in adata.var_names")

# Ensure batch categories are ordered correctly
adata.obs['batch'] = adata.obs['batch'].astype('category')
adata.obs['batch'].cat.reorder_categories(['Cortex', 'Striatum'])

# Define the custom colors for Cortex and Striatum
region_colors = {'Cortex': '#1d26cf', 'Striatum': '#f76f3e'}

plt.rcParams["figure.figsize"] = (10, 10)
sc.set_figure_params(fontsize = 20)
plt.rcParams["figure.figsize"] = (10, 10)

# Plot UMAP with custom colors for the 'region' column
sc.pl.umap(adata, color='batch', frameon=False, size=10, save='_region.svg', palette=region_colors)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NKX2-1', 'CRABP1', 'MAF','CHRNA3', 'TAC3', 'TH'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 6, save = 'importantmarkers.svg')

plt.rcParams["figure.figsize"] = (10, 10)

sc.pl.umap(adata, color = ['NKX2-1', 'CRABP1', 'RBP4','MAF', 'COL19A1','CHRNA3', 'TAC3', 'TH', 'LHX8' ], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3, save = 'crabp1markers')

adata.obs['leiden_subcluster_CRABPY2'] = adata.obs['leiden_subcluster_CRABPY2'].cat.reorder_categories(['Progenitor', 'MGE_LHX6/MAF',  'MGE_LHX6/NPY', 'MGE_CRABP1/MAF',
                                                                                                       'LGE_MEIS2/PAX6','LGE_FOXP2/TSHZ1','LGE_FOXP1/PENK', 'LGE_FOXP1/ISL1','LGE_FOXP1/ISL1/NPY1R', 'CGE_NR2F2/PROX1'
                                                                                                       ])


# Copy the 'leiden_subcluster_CRABPY2' observation to a new 'class' observation
adata.obs['class'] = adata.obs['leiden_subcluster_CRABPY2']

#Matrix of marker genes for inhibitory neuron clusters

#marker genes based on Schmitz et al. 2022
marker_genes_dict = {
    'Progenitor' : ['ASPM', 'CENPF', 'CKB'],
    'MGE_LHX6/MAF': ['LHX6', 'RBFOX1', 'SOX6', 'MEF2C', 'CUX2'],
    'MGE_LHX6/NPY': ['SST','NPY'],
    'MGE_CRABP1/MAF': ['ETV1','CRABP1','MAF', 'COL19A1'],
    'LGE_MEIS2/PAX6': ['MEIS2', 'PAX6'],
    'LGE_FOXP2/TSHZ1': ['EYA2','FOXP2','FOXP4','TSHZ1'],
    'LGE_FOXP1/PENK': ['FOXP1', 'SIX3', 'PENK'],
    'LGE_FOXP1/ISL1': ['RXRG', 'RARB','PBX3', 'ISL1'],
    'LGE_FOXP1/ISL1/NPY1R': ['NPY1R'],
    'CGE_NR2F2/PROX1': ['NR2F2', 'PROX1', 'SP8','NR3C2', 'PDZRN3',],

}


sc.pl.matrixplot(adata, groupby='leiden_subcluster_CRABPY2', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', save='matrix.svg', swap_axes=True)

# Compute proportions of Cortex vs. Striatum for each cluster
cluster_region_counts = adata.obs.groupby(['leiden_subcluster_CRABPY2', 'batch']).size().unstack(fill_value=0)
cluster_region_props = cluster_region_counts.div(cluster_region_counts.sum(axis=1), axis=0)

# Ensure order matches matrixplot
cluster_order = adata.obs['leiden_subcluster_CRABPY2'].cat.categories
cluster_region_props = cluster_region_props.loc[cluster_order]


# Set up figure with two subplots: one for barplot, one for matrixplot
fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 10]}, figsize=(20, 10))

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
    groupby='leiden_subcluster_CRABPY2',
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

adata.write('pathtorabbit/rabbit_processed.h5ad')


# Rename specific columns in obs
adata.obs.rename(columns={'leiden_subcluster_CRABPY': 'leiden_subcluster','leiden_subcluster_CRABPY2': 'leiden_subcluster2'}, inplace=True)

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
obs.to_csv('pathtogeo/rabbit_barcodes.tsv', sep='\t', index=True)

# Save var.tsv
var.to_csv('pathtogeo/rabbit_genes.tsv', sep='\t', index=True)

# Save matrix.mtx
scipy.io.mmwrite('pathtogeo/rabbit_matrix.mtx', matrix)

# # subset data 


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

adata.X

#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()

adata.X = np.array(adata.X)  # Convert to numpy array
adata.write('pathtorabbit/rabbit_raw.h5ad')

adata.write('pathtorabbit/rabbit_devinhibneurons.h5ad')

