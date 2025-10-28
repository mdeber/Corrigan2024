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
adata=sc.read_10x_h5('/pathto/opossum/updatedgenome/striatum1/OpossumP20Striatum1_cellbender_outputlr5_filtered.h5')
Results_file='/pathto/opossum/scRNAseq_OpossumP20striatum1.h5ad'
adata



adata.obs_names_make_unique()
adata.var_names_make_unique()



#h5ad stores whole anndata data structure
adata.write(Results_file)

#Define string to append to saved figures
save_name = "opossum_striatum1_p20"


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
adata = adata[adata.obs.percent_ribo < .2, :]
adata = adata[adata.obs.n_genes < 6000, :]
adata = adata[adata.obs.n_genes > 1000, :]


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


adata.write('/pathto/opossum/scRNAseq_striatum1_opossum_P20_prenorm.h5ad')




# Analyzing all together


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

path = '/pathto/opossum/'

P20cortex2 = sc.read_h5ad(path + 'scRNAseq_cortex2_opossum_P20_prenorm.h5ad')
P20striatum1 = sc.read_h5ad(path + 'scRNAseq_striatum1_opossum_P20_prenorm.h5ad')
P20cortex1 = sc.read_h5ad(path + 'scRNAseq_cortex1_opossum_P20_prenorm.h5ad')
P20striatum2 = sc.read_h5ad(path + 'scRNAseq_striatum2_opossum_P20_prenorm.h5ad')


#set a path to your working directory
Results_file=('/pathto/opossum/scRNAseq_opossum_concat.h5ad')



for gem in [P20cortex2, P20striatum1, P20striatum2, P20cortex1]:
    print(gem.shape[0])



adata = anndata.AnnData.concatenate(P20cortex2, P20striatum1, P20striatum2, P20cortex1, join='outer', batch_categories=['P20cortex2', 'P20striatum1', 'P20striatum2','P20cortex1'])


# Make a new column 'region' based on 'batch'
adata.obs['region'] = np.where(
    adata.obs['batch'].str.contains('cortex', case=False, na=False),
    'Cortex',
    np.where(
        adata.obs['batch'].str.contains('striatum', case=False, na=False),
        'Striatum',
        'Other'  # or np.nan if you want to leave unmatched values empty
    )
)

adata.obs_names_make_unique()
adata.var_names_make_unique()

#Set the .raw attribute of the AnnData object
adata.raw = adata

#h5ad stores whole anndata data structure
adata.write(Results_file)


#Define string to append to saved figures
save_name = "opossum_concat"


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

#Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

adata.write(Results_file)


# ## PCA


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
pca_var_genes = ["VIM", "MEF2C", "MYO16", "DLX6", "DLX5",
                 "GAD2"]
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["batch", 'percent_mito', 'percent_ribo', 'n_genes', 'n_counts']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)



sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']


#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["batch", 'percent_mito', 'percent_ribo', 'n_genes', 'n_counts']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)



#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

#Save the AnnData object as an .h5ad file.
adata.write(Results_file)


adata = sc.read_h5ad(Results_file)


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "OpossumP20_concat" + str(resolution), frameon = False, legend_loc = "on data")


# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2, save = 'umapconcat_qc.png')

sc.pl.umap(adata, color = ['batch'], wspace=0.25, ncols = 2)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['SLA', 'NEUROD2', 'SLC17A6', 'DLX1', 'DLX2', 'DLX5', 'DLX6','ZIC1', 'ZIC4','FOXG1', 'NKX2-1', 'PTPRK', 'TNS1', 'PDE3A',  'MEF2C', 'DLX1','LHX6', 'CRABP1','ETV1','TAC3', 'CHRNA3', 'CHRNA4', 'CHRNA7','NPY', 'SST','GAD1','ZIC1','NXPH2', 'LHX8', 'RBP4', 'CHAT', 'NR2F2', 'PROX1', 'MEIS2', 'PAX6', 'SCGN', 'TSHZ1','FOXP1', 'PENK', 'ISL1', 'NPY1R'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['PAX6', 'TBR1', 'SATB2', 'CUX2', 'EOMES', 'HES1','OLIG2', 'EGFR', 'ERBB4'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


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
adata.write('/pathto/opossum/scRNAseq_opossuminhibneuronsupdated.h5ad')


# # Renormalize, Recalculate PCA, recalculate UMAP



adata = sc.read_h5ad('/pathto/opossum/scRNAseq_opossuminhibneuronsupdated.h5ad')



#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()



# Directly access the raw count matrix - calling it "processed counts" because later .X will store your processed data, while the "counts" layer will 
processed_counts = adata.X

adata.layers["counts"] = adata.X.copy() #Copy the raw data into a "layer" in the AnnData structure - this will be used for normalization, etc.

adata.layers["counts"] #This should resemble the raw data at this point

print(adata.X.shape)


#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)

#Convert adata.X to numpy.ndarray
adata.X = np.array(adata.X)

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')


#Inspect the contribution of each PC to the variance of the data.
plt.rcParams["figure.figsize"] = (5, 5)
sc.pl.pca_variance_ratio(adata, log=True)


#Visualize the PCA loadings.
plt.rcParams["figure.figsize"] = (10, 5) 
sc.pl.pca_loadings(adata, components = [1, 2, 3, 4])



#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["CENPF", 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

#batch correction with harmony 
sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["batch", 'percent_mito', 'percent_ribo', 'n_genes', 'n_counts']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")


#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

# Check layers for numpy.matrix
for layer_name, layer_data in adata.layers.items():
    if isinstance(layer_data, np.matrix):
        print(f"Layer {layer_name} is a numpy.matrix")

# Check other attributes like 'X', 'obsm', 'varm', etc.
if isinstance(adata.X, np.matrix):
    print("adata.X is a numpy.matrix")

if isinstance(adata.layers["counts"], np.matrix):
    adata.layers["counts"] = np.array(adata.layers["counts"])


#Save the AnnData object as an .h5ad file.
adata.write('/pathto/opossum/scRNAseq_opossuminhibneuronsupdated.h5ad')


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


adata = sc.read_h5ad('/pathto/opossum/scRNAseq_opossuminhibneuronsupdated.h5ad')


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Opossum_Inhibitory_Neurons_p20" + str(resolution), frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito','n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['DLX1', 'DLX2', 'GAD1', 'GAD2', 'SLA', 'NEUROD2', 'SLC17A6', 'TAC3','batch'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NKX2-1', 'LHX6', 'ETV1', 'SOX6', 'MEF2C', 'COL19A1', 'ZIC1', 'ZIC2', 'ZIC4', 'TAC3', 'CER1', 'TH', 'LHX8', 'CHRNA7', 'CHRNA3','STXBP6', 'TSHZ1', 'PENK', 'SCGN', 'MEIS2', 'PAX6', 'ISL1', 'NPY1R', 'NR2F2', 'PROX1', 'VIP'], use_raw= False,frameon = False, color_map = "PuRd", size= 25, ncols = 3)

#Remove contaminating cluster
adata = adata[adata.obs['leiden'].isin(['20']) == False]


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
pca_var_genes = ['batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

#batch correction with harmony 
sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["batch", 'percent_mito', 'percent_ribo', 'n_genes', 'n_counts']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Opossum_Inhibitory_Neurons_p20" + str(resolution), frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')


#Save the AnnData object as an .h5ad file.
adata.write('/pathto/opossum/opossum_processed.h5ad')


# # Load and plot processed data

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


adata = sc.read_h5ad('/pathto/opossum/opossum_processed.h5ad')



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


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Opossum_Inhibitory_Neurons_p20", frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito','n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['DLX1', 'DLX2', 'GAD1', 'GAD2', 'SLA', 'NEUROD2', 'SLC17A6', 'TAC3','batch'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['MKI67', 'MEIS2', 'PAX6', 'ISL1', 'NPY1R', 'TSHZ1', 'PENK', 'SCGN','NR2F2', 'PROX1', 'VIP'], use_raw= False,frameon = False, color_map = "PuRd", size= 25, ncols = 3)

#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata, groupby = 'leiden', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Opossum_Inhibitory_Neurons_p20", frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['PAX6', 'MEIS2', 'SCGN', 'TSHZ1', 'FOXP1', 'ISL1', 'NPY1R','NR2F2', 'PROX1', 'PENK', 'VIP', 'CCK'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

# Define a dictionary to map old cluster names to new names
cluster_name_mapping = {
    '17' : 'MGE_CRABP1/MAF',
    '0' : 'MGE_LHX6/MAF',
    '5' : 'MGE_LHX6/MAF',
    
    '18' : 'Progenitor',
    '15' : 'Progenitor',
    '1' : 'Progenitor',
    '12' : 'Progenitor',
    '2' : 'Progenitor',
    
    
    

    '3' : 'LGE_MEIS2/PAX6',
    '5' : 'LGE_MEIS2/PAX6',
    '16' : 'LGE_MEIS2/PAX6',

    
    '8' : 'CGE_NR2F2/PROX1',
    '19' : 'CGE_NR2F2/PROX1',
    '4' : 'CGE_NR2F2/PROX1',
    '7' : 'CGE_NR2F2/PROX1',

    '13' : 'LGE_FOXP2/TSHZ1',
    '9' : 'LGE_FOXP2/TSHZ1',
    
    '10' : 'LGE_FOXP1/PENK',
    
    '6' : 'LGE_FOXP1/ISL1',
    '11' : 'LGE_FOXP1/ISL1',
    '14' : 'LGE_FOXP1/ISL1',
    

 
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['leiden'] = adata.obs['leiden'].replace(cluster_name_mapping)


sc.pl.umap(adata, color= ['leiden'], title = "Ferret" + str(1), frameon = False)

paldict={'CGE_NR2F2/PROX1': 'mediumpurple',
    'LGE_FOXP1/ISL1': 'royalblue',
    'LGE_FOXP1/ISL1/NPY1R': 'blue',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_FOXP2/TSHZ1': '#17344c',
    'LGE_MEIS2/PAX6': 'orangered',
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TAC3': 'fuchsia',
    'MGE_LHX6/MAF': 'skyblue',
    'MGE_LHX6/NPY': 'teal',
    'VMF_ZIC1/ZIC2': 'green',
    'VMF_CRABP1/LHX8': 'mediumseagreen',
    'RMTW_ZIC1/RELN': 'lightcoral',
    'Progenitor': 'pink',

}


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= ['class'], palette=paldict, title = "Opossum_P20_inhibitory", frameon = False, size = 10,save = '_officialclusters.svg')

# Define the custom colors for Cortex and Striatum
region_colors = {'Cortex': '#1d26cf', 'Striatum': '#f76f3e'}

plt.rcParams["figure.figsize"] = (10, 10)
sc.set_figure_params(fontsize = 20)
plt.rcParams["figure.figsize"] = (10, 10)

# Plot UMAP with custom colors for the 'region' column
sc.pl.umap(adata, color='region', frameon=False, size=10, save='_region.svg', palette=region_colors)


adata.obs['leiden'] = adata.obs['leiden'].cat.reorder_categories(['Progenitor', 'MGE_LHX6/MAF',  'MGE_CRABP1/MAF',
                                                                                                       'LGE_MEIS2/PAX6', 'LGE_FOXP2/TSHZ1','LGE_FOXP1/PENK', 
                                                                                                        'LGE_FOXP1/ISL1', 'CGE_NR2F2/PROX1',
                                                                                                       ])



#Matrix of marker genes for inhibitory neuron clusters
sc.set_figure_params(fontsize = 15)
#marker genes based on Schmitz et al. 2022
marker_genes_dict = {
    'Progenitor' : ['ASPM', 'CENPF', 'MKI67'],
    'MGE_LHX6/MAF': ['LHX6', 'MAF', 'SOX6', 'MEF2C'],
    'MGE_CRABP1/MAF': ['MAFB','ETV1', 'COL19A1', 'CNR1', 'MTUS2',  'KIT'],
    'LGE_MEIS2/PAX6': ['MEIS2', 'PAX6'],
    'LGE_FOXP2/TSHZ1': ['EYA2','FOXP2', 'FOXP4','TSHZ1'],
    'LGE_FOXP1/PENK': ['FOXP1', 'SIX3', 'PENK'],
    'LGE_FOXP1/ISL1': ['ISL1', 'RXRG', 'RARB'],
    'CGE_NR2F2/PROX1': ['PDZRN3', 'NR3C2', 'NR2F2', 'PROX1', 'NPAS1'],

}

sc.pl.matrixplot(adata, groupby='leiden', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes =True, save='matrix.svg')



# Compute proportions of Cortex vs. Striatum for each cluster
cluster_region_counts = adata.obs.groupby(['leiden', 'region']).size().unstack(fill_value=0)
cluster_region_props = cluster_region_counts.div(cluster_region_counts.sum(axis=1), axis=0)

# Ensure order matches matrixplot
cluster_order = adata.obs['leiden'].cat.categories
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
    groupby='leiden',
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


adata.obs['class'] = adata.obs['leiden']

sc.set_figure_params(fontsize = 50)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NKX2-1', 'ETV1', 'MAF','CHRNA3', 'TAC3', 'TH'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3, save='markers.svg')

#Save the AnnData object as an .h5ad file.
adata.write('/pathto/opossum/scRNAseq_opossuminhibneuronsupdated.h5ad')


adata = sc.read_h5ad('/pathto/opossum/scRNAseq_opossuminhibneuronsupdated.h5ad')


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
obs.to_csv('/pathto/geo/opossum_barcodes.tsv', sep='\t', index=True)

# Save var.tsv
var.to_csv('/pathto/geo/opossum_genes.tsv', sep='\t', index=True)

# Save matrix.mtx
scipy.io.mmwrite('/pathto/geo/opossum_matrix.mtx', matrix)



# ## Subset data for integrations


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

adata = adata_subset

#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()


adata.X = np.array(adata.X)  # Convert to numpy array
adata.write('/pathto/opossum/opossum_subsetraw1k.h5ad')
