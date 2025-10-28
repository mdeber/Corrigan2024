#!/usr/bin/env python
# coding: utf-8


##Preprocessing by sample
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
adata=sc.read_10x_h5('/pathto/ferret/updated_genome/cellbender_out/ventral/FerretP1ventral_cellbender_outputv11_filtered.h5')
Results_file='/pathto/ferret/updated_genome/revisions/scRNAseq_ventralupdatedv11_ferretp1.h5ad'
adata



adata.obs_names_make_unique()
adata.var_names_make_unique()

#h5ad stores whole anndata data structure
adata.write(Results_file)


#Define string to append to saved figures
save_name = "ferret_p1_ventral"


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

sc.pl.violin(adata,keys=['percent_ribo','percent_mito','n_counts', 'n_genes'],  rotation=90, multi_panel=True, save= save_name + 'QC_metrics')

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
adata = adata[adata.obs.percent_ribo < .1, :]
adata = adata[adata.obs.percent_mito < .1, :]
adata = adata[adata.obs.n_genes > 1750, :]
adata = adata[adata.obs.n_genes < 6000, :]

sc.pl.violin(adata,keys=['percent_ribo', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True, save= save_name + 'QC_metrics')

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


# Estimate and predict doublets using scrublet and the expected doublet rate in the relevant 10X Genomics protocol.
sce.pp.scrublet(adata, expected_doublet_rate = 0.08, n_prin_comps = 50, verbose = True)


sce.pl.scrublet_score_distribution(adata, save = save_name)


# Add column to AnnData.obs with categorical Singlet/Doublet instead of boolean True/False.
conditions = [
    (adata.obs["predicted_doublet"] == True),
    (adata.obs["predicted_doublet"] == False)]

values =['Doublet', 'Singlet']

adata.obs['doublet_info'] = np.select(conditions, values)

sc.pl.violin(adata, 'doublet_score', groupby = 'doublet_info', rotation=45, save = "scrublet_doublet_score_vln.png")



# Remove doublets
adata = adata[adata.obs['predicted_doublet'] == False]


adata.write('/pathto/ferret/updated_genome/revisions/scRNAseq_ventralupdatedv11_ferretp1_prenorm.h5ad')




# ## Ferret ALL concat analysis


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

path = '/pathto/ferret/updated_genome/revisions/'

P1cortex = sc.read_h5ad(path + 'scRNAseq_dorsalupdatedv11_ferretp1_prenorm.h5ad')
P1striatum = sc.read_h5ad(path + 'scRNAseq_ventralupdatedv11_ferretp1_prenorm.h5ad')
P5striatum = sc.read_h5ad(path + 'scRNAseq_striatalupdatedv11_ferretp5_prenorm.h5ad')
P5cortex = sc.read_h5ad(path + 'scRNAseq_cortexupdatedv11_ferretp5_prenorm.h5ad')
P14striatum1 = sc.read_h5ad(path + 'scRNAseq_striatum1updatedv11_ferretp14_prenorm.h5ad')
P14striatum2 =  sc.read_h5ad(path + '/scRNAseq_striatum2updatedv11_ferretp14_prenorm.h5ad')
P14cortex1 = sc.read_h5ad(path + 'scRNAseq_cortex1updatedv11_ferretp14_prenorm.h5ad')
P14cortex2 = sc.read_h5ad(path + 'scRNAseq_cortex2updatedv11_ferretp14_prenorm.h5ad')


#set a path to your working directory
Results_file=('/pathto/ferret/updated_genome/revisions/scRNAseq_ferretupdatedv11_concat.h5ad')


adata = anndata.AnnData.concatenate(P1cortex, P1striatum, P5cortex, P5striatum, P14striatum1, P14striatum2, P14cortex1, P14cortex2, join='outer', batch_categories=['P1cortex', 'P1striatum', 'P5cortex', 'P5striatum', 'P14striatum1', 'P14striatum2', 'P14cortex1', 'P14cortex2'])


adata.obs_names_make_unique()
adata.var_names_make_unique()




#Set the .raw attribute of the AnnData object
adata.raw = adata

#h5ad stores whole anndata data structure
adata.write(Results_file)


#Define string to append to saved figures
save_name = "ferret_concat"


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


# Cluster umap embeddings using leiden and save umap plots.
resolution = 2 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Ferret_P1_concat" + str(resolution), frameon = False, legend_loc = "on data")



# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2, save = 'umapconcat_qc.png')



sc.pl.umap(adata, color = ['batch'], wspace=0.25, ncols = 2)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['FOXG1', 'NKX2-1',  'MEF2C', 'DLX1','LHX6', 'CRABP1','ETV1','TAC3', 'CHRNA3', 'CHRNA4', 'CHRNA7','NPY', 'SST','GAD1','ZIC1','NXPH2', 'LHX8', 'RBP4', 'CHAT', 'NR2F2', 'PROX1', 'MEIS2', 'PAX6', 'SCGN', 'TSHZ1','FOXP1', 'PENK', 'ISL1', 'NPY1R'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)



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


#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()


# Directly access the raw count matrix - calling it "processed counts" because later .X will store your processed data, while the "counts" layer will 
processed_counts = adata.X
print(processed_counts)


adata.layers["counts"] = adata.X.copy() #Copy the raw data into a "layer" in the AnnData structure - this will be used for normalization, etc.


#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)


# Convert adata.X to numpy.ndarray
adata.X = np.array(adata.X)

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

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
                 "CUX2", 'batch']
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



#leiden clustering
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

sc.pl.umap(adata, color= 'leiden', title = "Ferret_Inhibitory_Neurons", frameon = False, legend_loc = "on data")




# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito','n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)
#sc.pl.violin(adata, ['pct_counts_ribo','pct_counts_mito'], groupby = "leiden")




excitatorymat = adata[:, ['SLA', 'SLC17A6', 'NEUROD2', 'NEUROG2', 'NEUROD6', 'EOMES']].to_df()
excitatorymat['leiden'] = adata.obs['leiden']
meanmat2 = excitatorymat.groupby('leiden').mean()
print(meanmat2>meanmat2.mean(0))
boolmat2 = (meanmat2 > meanmat2.mean(0)).sum(1)




#Filter data -- removing excitatory clusters
adata = adata[boolmat2[adata.obs['leiden']].values <= 3]

# Print information about the filtering results
print(f"Number of cells after filtering excitatory contamination: {len(adata)}")



adata.layers["counts"]



adata.X = adata.layers["counts"].todense()


print(adata.X)



#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')



sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']




#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")
#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)



sc.tl.leiden(adata, resolution = 1)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Ferret_Inhibitory_Neurons", frameon = False, legend_loc = "on data")



# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito','n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['DLX1', 'DLX2', 'GAD1', 'GAD2', 'SLA', 'NEUROD2', 'SLC17A6','EOMES', 'TBR1', 'batch'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['TAC3', 'CRABP1', 'PLPP4', 'LHX8','COL19A1', 'CHRNA3', 'TH', 'PVALB', 'PTHLH', 'ETV1','MAF', 'TAC3', 'STXBP6','ZIC1','NXPH2', 'LHX8', 'RBP4', 'NR2F2', 'TH' ], 
           use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['CHRNA3', 'CHRNA7', 'RBP4', 'SHISA6', 'ZIC1', 'ZIC4','ETV5', 'NXPH2', 'NXPH4', 'NXPH1', 'STXBP6', 'TAC3', 'LHX8', 'PLPP4', 'COL23A1', 'TRH', 'PTPRK', 'SMOC1', 'SLIT2', 'TNS1', 'NPTX1', 'PDE3A', 'IGF1', 'PLCXD3',], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)




plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NKX2-1', 'MEF2C',  'ETV1','CRABP1', 'TAC3', 'MAF', 'MEIS2', 'PAX6', 'FOXP1', 'FOXP2', 'TSHZ1', 'SCGN', 'NPY1R', 'ISL1','NR2F2', 'PROX1','CCK', 'LHX8', 'DLX1', 'GAD1', 'HES1', 'SPC24', 'CDC20', 'TOP2A'], 
           use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


adata.write(Results_file)


# # Final clean up of doublets/low quality clusters

# Remove clusters that are doublets or low quality
clusters_to_remove = ['9', '18', '16', '9', '21', '22', '12', '0', '19', '24', '20']  




adata = adata_filtered

adata.X = adata.layers["counts"].todense()

print(adata.X)



#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4) #Specify "layer = "counts" to use the counts layer for normalization
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.X = np.array(adata.X)

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')


#batch correction with harmony 
sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']



#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")


#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)



# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Ferret_Inhibitory_Neurons" + str(resolution), frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['DLX1', 'DLX2', 'GAD1', 'GAD2', 'SLA', 'NEUROD2', 'SLC17A6', 'batch'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)
#excitatory neuron markers 'Sla', 'Neurod2', 'Slc17a6'


sc.set_figure_params(fontsize = 25)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['CRABP1','MAF','TAC3', 'LHX8','TH','CHRNA3', 'CHRNA7', 'COL19A1', 'CHRNA7', 'ZIC1', 'ZIC2','ZIC4'],use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 4)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['CHRNA3', 'TAC3', 'LHX8', 'PLPP4', 'COL23A1', 'TRH', 'PTPRK', 'SMOC1', 'SLIT2', 'TNS1', 'NPTX1', 'PDE3A', 'IGF1', 'PLCXD3','VAV3'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NKX2-1', 'MEF2C',  'ETV1', 'MAF', 'MEIS2', 'PAX6', 'FOXP1', 'FOXP2', 'NR2F2', 'LHX8', 'DLX1', 'GAD1', 'HES1', 'SPC24', 'CDC20', 'TOP2A', 'CHAT', 'PVALB', 'PTHLH'], 
           use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 4)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NPY', 'ISL1', 'NPY1R', 'TSHZ1', 'PAX6', 'MEIS2', 'SCGN', 'PENK', 'NR2F2', 'PROX1'], 
           use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Ferret_Inhibitory_Neurons", frameon = False, legend_loc = "on data")


# Remove clusters that are doublets or low quality
clusters_to_remove = ['8']  # Make sure these are strings if leiden clusters are stored as strings
adata_filtered = adata[~adata.obs['leiden'].isin(clusters_to_remove)].copy()



adata = adata_filtered


adata.X = adata.layers["counts"].todense()



#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4) #Specify "layer = "counts" to use the counts layer for normalization
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)



adata.X = np.array(adata.X)


sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')


#batch correction with harmony 
sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")




#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)


# # Final dataset analysis and classifications


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
sc.pl.umap(adata, color= 'leiden', title = "Ferret_Inhibitory_Neurons" + str(resolution), frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')




plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['DLX1', 'DLX2', 'GAD1', 'GAD2', 'SLA', 'NEUROD2', 'SLC17A6', 'batch'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['CRABP1','MAF','TAC3', 'LHX8','TH','CHRNA3', 'CHRNA7', 'COL19A1', 'CHRNA7', 'ZIC1', 'ZIC2','ZIC4', 'CHAT', 'RELN'],use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 4)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NPY', 'ISL1', 'NPY1R', 'TSHZ1', 'PAX6', 'MEIS2', 'SCGN', 'PENK', 'NR2F2', 'PROX1'], 
           use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['MEIS2', 'PAX6', 'SCGN', 'leiden'], 
           use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 2)


# Create the new 'region' column based on conditions
adata.obs['region'] = adata.obs['batch'].apply(
    lambda x: 'cortex' if 'cortex' in x else ('striatum' if 'striatum' in x else None)
)


#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata, groupby = 'leiden_subcluster_CRABPY3', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df



#subcluster TAC3 cluster so that it is only TAC3 cells
sc.tl.leiden(adata, restrict_to=('leiden', ['4']), resolution=0.5, key_added='leiden_subcluster_CRABPY')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY'], title = "subcluster" + str(resolution), frameon = False, legend_loc = 'on data')



#subcluster TAC3 cluster so that it is only TAC3 cells
sc.tl.leiden(adata, restrict_to=('leiden_subcluster_CRABPY', ['0']), resolution=0.3, key_added='leiden_subcluster_CRABPY2')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY2'], title = "Pig_subcluster" + str(resolution), frameon = False, legend_loc = 'on data')


#subcluster TAC3 cluster so that it is only TAC3 cells
sc.tl.leiden(adata, restrict_to=('leiden_subcluster_CRABPY2', ['3']), resolution=0.2, key_added='leiden_subcluster_CRABPY3')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY3'], title = "Pig_subcluster" + str(resolution), frameon = False, legend_loc = 'on data')



# Define a dictionary to map old cluster names to new names
cluster_name_mapping = {
    '4,4' : 'MGE_CRABP1/TAC3',
    
    '6' : 'Progenitor',
    '13' : 'Progenitor',
    '7' : 'Progenitor',
    '10' : 'Progenitor',
    '17' : 'Progenitor',
    '9' : 'Progenitor',
    '5' : 'Progenitor',
    '2' : 'Progenitor',
    
    
    '14' : 'LGE_FOXP1/PENK',
    
    '4,0' : 'MGE_LHX6/MAF',
    '4,1' : 'MGE_LHX6/MAF',
    '4,3' : 'MGE_LHX6/MAF',
    '4,2' : 'MGE_LHX6/MAF',
    '4,5' : 'MGE_LHX6/MAF',
    '18' : 'MGE_LHX6/MAF',
    '11' : 'MGE_LHX6/MAF',
    '4,0' : 'MGE_LHX6/MAF',
    '0,0' : 'MGE_LHX6/MAF',
    '0,1' : 'MGE_LHX6/MAF',
    '15' : 'MGE_LHX6/MAF',
    
    '0,2' : 'MGE_LHX6/NPY',
    
    
    '1' : 'CGE_NR2F2/PROX1',
    

    '8' : 'LGE_FOXP2/TSHZ1',
    
    '3,1' : 'LGE_MEIS2/PAX6',
    '3,0' : 'LGE_MEIS2/PAX6',

    
    '12' : 'LGE_FOXP1/ISL1/NPY1R',
    '16' : 'LGE_FOXP1/ISL1',
    
    '3,2' : 'VMF_ZIC1/ZIC2',
    '19' : 'VMF_ZIC1/ZIC2',

    
 
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['leiden_subcluster_CRABPY3'] = adata.obs['leiden_subcluster_CRABPY3'].replace(cluster_name_mapping)



sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY3'], title = "Ferret" + str(1), frameon = False)



paldict={'CGE_NR2F2/PROX1': 'mediumpurple',
    'LGE_FOXP1/ISL1': 'royalblue',
    'LGE_FOXP1/ISL1/NPY1R': 'blue',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_FOXP2/TSHZ1': '#17344c',
    'LGE_MEIS2/PAX6': 'orangered',
    'MGE_CRABP1/MAF': 'indigo',
    'LGE_MEIS2/PAX6/SCGN': 'orange',
    'MGE_CRABP1/TAC3': 'fuchsia',
    'MGE_LHX6/MAF': 'skyblue',
    'MGE_LHX6/NPY': 'teal',
    'VMF_ZIC1/ZIC2': 'green',
    'VMF_CRABP1/LHX8': 'mediumseagreen',
    'RMTW_ZIC1/RELN': 'lightcoral',
    'Progenitor': 'pink',

}



sc._settings.settings._vector_friendly=True


sc.set_figure_params(dpi=80, dpi_save=500)




plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= ['class'], palette=paldict, title = "Ferret_inhibitory", frameon = False, size = 10, save = '_officialclusters.svg')


sc.set_figure_params(fontsize = 50)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NKX2-1', 'CRABP1', 'MAF','CHRNA3', 'TAC3', 'TH'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3, save='markers.svg')


sc.set_figure_params(fontsize = 50)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['ZIC1', 'RELN', 'CRABP1', 'LHX8', 'ZIC2'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 2)


plt.rcParams["figure.figsize"] = (10, 10)
sc.set_figure_params(fontsize = 20)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= ['batch'], frameon = False,size = 10,  save = '_batch.svg')



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



adata.obs['leiden_subcluster_CRABPY3'] = adata.obs['leiden_subcluster_CRABPY3'].cat.reorder_categories(['Progenitor', 'MGE_LHX6/MAF', 'MGE_LHX6/NPY', 'MGE_CRABP1/TAC3',
                                                                                                       'LGE_MEIS2/PAX6','LGE_FOXP2/TSHZ1','LGE_FOXP1/PENK', 


#Matrix of marker genes for inhibitory neuron clusters
sc.set_figure_params(fontsize = 13)
#marker genes based on Schmitz et al. 2022
marker_genes_dict = {
    'Progenitor' : ['ASPM', 'CENPF', 'MKI67'],
    'MGE_LHX6/MAF': ['LHX6', 'RBFOX1', 'SOX6', 'MEF2C'],
    'MGE_LHX6/NPY': ['NPY', 'TRPC6'],
    'MGE_CRABP1/TAC3': ['CRABP1','ETV1', 'KIT','NRTN', 'TRHDE','STXBP6','TAC3', 'CHRNA3', 'CHRNA7', 'PARD3',  'PDE1C'],
    'LGE_MEIS2/PAX6': ['MEIS2', 'PAX6'],
    'LGE_FOXP2/TSHZ1': ['EYA2','FOXP2', 'FOXP4','TSHZ1' ],
    'LGE_FOXP1/PENK': ['FOXP1', 'SIX3', 'PENK'],
    'LGE_FOXP1/ISL1': ['ISL1', 'RXRG', 'RARB', 'PBX3'],
    'LGE_FOXP1/ISL1/NPY1R' : ['NPY1R'],
    'CGE_NR2F2/PROX1': ['PDZRN3', 'NR3C2', 'NR2F2', 'PROX1'],
    'VMF_ZIC1/ZIC2' : ['ZIC1', 'ZIC2', 'ZIC4'],

}



sc.pl.matrixplot(adata, groupby='leiden_subcluster_CRABPY3', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes =True, save='matrix.svg')



# Compute proportions of Cortex vs. Striatum for each cluster
cluster_region_counts = adata.obs.groupby(['leiden_subcluster_CRABPY3', 'region']).size().unstack(fill_value=0)
cluster_region_props = cluster_region_counts.div(cluster_region_counts.sum(axis=1), axis=0)

# Ensure order matches matrixplot
cluster_order = adata.obs['leiden_subcluster_CRABPY3'].cat.categories
cluster_region_props = cluster_region_props.loc[cluster_order]

# Set up figure with two subplots: one for barplot, one for matrixplot
fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 10]}, figsize=(12, 4.5))

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



# Compute proportions of Cortex vs. Striatum for each cluster
cluster_region_counts = adata.obs.groupby(['leiden_subcluster_CRABPY3', 'region']).size().unstack(fill_value=0)
cluster_region_props = cluster_region_counts.div(cluster_region_counts.sum(axis=1), axis=0)

# Ensure order matches matrixplot
cluster_order = adata.obs['leiden_subcluster_CRABPY3'].cat.categories
cluster_region_props = cluster_region_props.loc[cluster_order]

### **Plot 1: Horizontal Stacked Barplot for Cortex vs. Striatum Proportions (Swapped axes)**
fig, ax = plt.subplots(figsize=(6, 10))   # Larger figure for better spacing

# Horizontal stacked barplot with explicit color definitions
cluster_region_props.plot(kind='barh', stacked=True, color=['blue', 'orange'], ax=ax, edgecolor='none')

# Format barplot
ax.set_yticks(range(len(cluster_order)))  # Ensure correct number of ticks for y-axis
ax.set_yticklabels(cluster_order, rotation=0)  # No rotation needed for y-axis labels
ax.set_xlabel('Proportion')  # X-axis label
ax.set_title('Cortex vs. Striatum Proportions per Cluster')

# Remove vertical grid lines (optional)
ax.grid(False)

# Move the legend to the right (outside the plot)
ax.legend(
    title="Region", 
    labels=['Cortex', 'Striatum'], 
    loc='upper left',
    bbox_to_anchor=(1.05, 1)  # Moves the legend outside of the plot
)

# Adjust layout for better spacing
plt.tight_layout()

# Save the plot as a .png file
plt.savefig("cortex_vs_striatum_proportions_swapped.png", dpi=300, bbox_inches='tight')

# Show the plot
plt.show()



clusters_of_interest = adata.obs['leiden_subcluster_CRABPY3'].isin(['MGE_CRABP1/TAC3'])
adata_filtered = adata[clusters_of_interest]

marker_genes_dict = ['LHX6', 'NKX2-1', 'SOX6', 'MAF', 'CHL1','CRABP1','ETV1','GALNTL6','OLFM2','NRTN', 'KIT','TRHDE','STXBP6','TAC3', 'CHRNA3', 'CHRNA7', 'PARD3',  'PDE1C', 'ZIC1', 'ZIC2', 'ZIC4', 'LHX8', 'CXCR4', 'UNC5C', 'HCN1', 'DBI','ADCY8', 'PCSK5']

sc.pl.dotplot(adata_filtered, groupby='region', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes=False,  figsize=(10,0.75),save='matrix.svg')



adata.obs['class'] = adata.obs['leiden_subcluster_CRABPY3']


adata.obs['class'].value_counts().get('MGE_CRABP1/TAC3', 0)



# Filter 
clusters_of_interest = adata.obs['class'].isin(['MGE_CRABP1/MAF', 'MGE_CRABP1/TAC3'])
adata_filtered = adata[clusters_of_interest]

# List of genes of interest
genes_of_interest = ['LHX6', 'CRABP1','ANGPT2', 'MAF', 'RBP4',  'ARX', 'TAC3', 'TRH', 'CHRNA3', 'CHRNA7', 'STXBP6', 'ZIC1', 'LHX8']

# Create a dot plot for the specified genes
sc.pl.dotplot(adata_filtered, var_names=genes_of_interest, use_raw= False, vmin=-2, vmax=2, swap_axes = True, cmap='PiYG_r',groupby='class', save = 'macdotplot.svg')



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
obs.to_csv('/pathto/geo/ferret_barcodes.tsv', sep='\t', index=True)

# Save var.tsv
var.to_csv('/pathto/geo/ferret_genes.tsv', sep='\t', index=True)

# Save matrix.mtx
scipy.io.mmwrite('/pathto/geo/ferret_matrix.mtx', matrix)




adata.write('/pathto/ferret/updated_genome/revisions/ferret_processed.h5ad')


# # Differential expression analysis


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

import scanpy as sc
import decoupler as dc

# Only needed for processing
import numpy as np
import pandas as pd

# Needed for some plotting
import matplotlib.pyplot as plt


adata= sc.read_h5ad('/pathto/ferret/updated_genome/revisions/ferret_processed.h5ad')



count = (adata.obs['class'] == 'MGE_CRABP1/TAC3').sum()
print(f"Number of cells with class 'MGE_CRABP1/TAC3': {count}")
# Filter cells with class 'MGE_CRABP1'
TAC3_cells = adata.obs[adata.obs['class'] == 'MGE_CRABP1/TAC3']

# Count how many are from 'cortex' and how many from 'striatum'
region_counts = TAC3_cells['region'].value_counts()
print(region_counts)




adata.obs['timepoint'] = np.where(adata.obs['batch'].str.contains('P14', na=False), 'P14',
                          np.where(adata.obs['batch'].str.contains('P5', na=False), 'P5',
                          np.where(adata.obs['batch'].str.contains('P1', na=False), 'P1', 'Unknown')))



pdata = dc.get_pseudobulk(
    adata,
    sample_col='batch',
    groups_col='class',
    layer='counts',
    mode='sum',
    min_cells=0,
    min_counts=0
)



sc.set_figure_params(fontsize = 10)

dc.plot_psbulk_samples(pdata, groupby=['batch', 'class'], figsize=(12, 4))


# Get filtered pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col='batch',
    groups_col='class',
    layer='counts',
    mode='sum',
    min_cells=5,
    min_counts=250
)
pdata




# Store raw counts in layers
pdata.layers['counts'] = pdata.X.copy()



# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)


# Return raw counts to X
dc.swap_layer(pdata, 'counts', X_layer_key=None, inplace=True)


sc.pl.pca(pdata, color=['batch', 'class', 'timepoint'], ncols=1, size=300)
sc.pl.pca_variance_ratio(pdata)


dc.get_metadata_associations(
    pdata,
    obs_keys = ['batch', 'region', 'timepoint', 'psbulk_n_cells', 'psbulk_counts'],  # Metadata columns to associate to PCs
    obsm_key='X_pca',  # Where the PCs are stored
    uns_key='pca_anova',  # Where the results are stored
    inplace=True,
)



dc.plot_associations(
    pdata,
    uns_key='pca_anova',  # Summary statistics from the anova tests
    obsm_key='X_pca',  # where the PCs are stored
    stat_col='p_adj',  # Which summary statistic to plot
    obs_annotation_cols = ['batch', 'region', 'timepoint','class'], # which sample annotations to plot
    titles=['Principle component scores', 'Adjusted p-values from ANOVA'],
    figsize=(7, 5),
    n_factors=10,
)


TAC3 = pdata[pdata.obs['class'] == 'MGE_CRABP1/TAC3'].copy()

dc.plot_filter_by_expr(TAC3, group='region', min_count=10, min_total_count=15)



#filter if needed
# Obtain genes that pass the thresholds
genes = dc.filter_by_expr(TAC3, group='region', min_count=10, min_total_count=15)

# Filter by these genes
TAC3 = TAC3[:, genes].copy()
TAC3


# Import DESeq2
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats


# Build DESeq2 object
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    adata=TAC3,
    design_factors='region',
    ref_level=['region', 'striatum'],
    refit_cooks=True,
    inference=inference,
)



# Compute LFCs
dds.deseq2()


stat_res = DeseqStats(
    dds,
    contrast=["region", 'striatum', 'cortex'],
    inference=inference,
)

# Compute Wald test
stat_res.summary()



# Extract results
results_df = stat_res.results_df
results_df

dc.plot_volcano_df(
    results_df,
    x='log2FoldChange',
    y='padj',
    top=20,
    figsize=(8, 4)
)



mat = results_df[['stat']].T.rename(index={'stat': 'TAC3'})
mat


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['LHX8', 'COL19A1','TAC3','CHRNA7','region'], 
           use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)




#Saving a subset for integration analysis

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



adata= adata_subset



print(adata_subset.raw.X)


#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()


adata.X = np.array(adata.X)  # Convert to numpy array


adata.write('/pathto/ferret/updated_genome/revisions/ferret_raw1K.h5ad')



