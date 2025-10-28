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
adata=sc.read_10x_h5('pathtopig/updated_genome/cellbender_out/E83striatum/PigE83striatum_cellbender_output_filtered.h5')
Results_file='pathtopig/scRNAseq_striatum_pig_E83.h5ad'


adata.obs_names_make_unique()
adata.var_names_make_unique()


#h5ad stores whole anndata data structure
adata.write(Results_file)

#Define string to append to saved figures
save_name = "pigs_tac3_striatum_E83"


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

#sc.pl.violin(adata,groupby='batch',keys=['percent_ribo','percent_mito', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True, 
             #save= save_name + 'QC_metrics')

sc.pl.violin(adata,keys=['percent_ribo','percent_mito', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True, save= save_name + 'QC_metrics')


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
adata = adata[adata.obs.percent_ribo < .1, :]
adata = adata[adata.obs.percent_mito < .1, :]
adata = adata[adata.obs.n_genes < 6000, :]
adata = adata[adata.obs.n_genes > 500, :]


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

#h5ad stores whole anndata data structure
adata.write(Results_file)

adata.write('pathtopig/scRNAseq_striatum_pig_E83_prenorm.h5ad')


# ## Identify doublets

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
adata.write('pathtopig/updated_genome/scRNAseq_striatumupdated_pig_E83_prenorm.h5ad')


# Process all samples together

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
Results_file=('pathtopig/updated_genome/scRNAseq_pig_concatupdatedrevisions.h5ad')

path = 'pathtopig/updated_genome/'

E73Striatum = sc.read_h5ad(path + 'scRNAseq_striatumupdated_pig_E73_prenorm.h5ad')
E73Cortex = sc.read_h5ad(path + 'scRNAseq_cortexupdated_pig_E73_prenorm.h5ad')
E83Striatum = sc.read_h5ad(path + 'scRNAseq_striatumupdated_pig_E83_prenorm.h5ad')
E83Cortex = sc.read_h5ad(path + 'scRNAseq_cortexupdated_pig_E83_prenorm.h5ad')

for gem in [E73Striatum, E73Cortex, E83Striatum, E83Cortex]:
    print(gem.shape[0])


adata = anndata.AnnData.concatenate(E73Striatum, E73Cortex, E83Striatum, E83Cortex, join='outer', batch_categories=['E73Striatum', 'E73Cortex', 'E83Striatum', 'E83Cortex'])


adata.raw = adata
adata.layers["counts"] = adata.X.copy()
adata.write(Results_file)


# # Normalize data and scale  

adata = sc.read_h5ad(Results_file)


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

sc.pp.scale(adata, max_value=10)
adata.write(Results_file)

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
pca_var_genes = ["NSG2", 'MT3', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["NSG2", 'MT3', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


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
sc.pl.umap(adata, color= 'leiden', title = "Integrated_pig" + str(resolution), frameon = False, legend_loc = "on data")

sc.pl.umap(adata, color = ['batch'], title = "by_batch")

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")


plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['NEUROD2','SLA', 'SLC17A6', 'FOXG1', 'NKX2-1', 'LHX6', 'DLX1', 'DLX2', 'GAD1', 'GAD2','CRABP1','COL19A1', 'CHRNA3', 'TH', 'PVALB', 'PTHLH', 'ETV1','MAF', 'TAC3', 'STXBP6','ZIC1','NXPH2', 'LHX8', 'RBP4', 'NR2F2', 'PROX1', 'CHAT', 'SST', 'NPY'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['TBR1', 'PAX6', 'OLIG1', 'OLIG2', 'HES1', 'HES5'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


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

#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

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
pca_var_genes = ["GRIA1", 'NRXN3', 'TPX2', 'CENPE', 'NRG3', 'ALCAM', 'DLX5', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']


#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["GRIA1", 'NRXN3', 'TPX2', 'CENPE', 'NRG3', 'ALCAM', 'DLX5', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)


#save inhibitory neurons as a separate file
adata.write('pathtopig/updated_genome/scRNAseq_pig_concat_inhibneurons_updatedrevisions.h5ad')


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




adata = sc.read_h5ad('pathtopig/updated_genome/scRNAseq_pig_concat_inhibneurons_updatedrevisions.h5ad')

adata.obs['region'] = np.where(
    adata.obs['batch'].str.contains('cortex', case=False, na=False), 'cortex',
    np.where(
        adata.obs['batch'].str.contains('striatum', case=False, na=False), 'striatum',
        'unknown'  
    )
)


adata.obs['timepoint'] = np.where(
    adata.obs['batch'].str.contains('E73', case=False, na=False), 'E73',
    np.where(
        adata.obs['batch'].str.contains('E83', case=False, na=False), 'E83',
        'unknown'  
)


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Pig_Inhibitory_Neurons", frameon = False, legend_loc = "on data", save = 'leidenclusters.svg')

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'timepoint', frameon = False)

sc.pl.umap(adata, color= 'region', title = "Pig_Inhibitory_Neurons", frameon = False)

sc.pl.umap(adata, color= 'batch', title = "Pig_Inhibitory_Neurons", frameon = False, save = 'batch.svg')


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
sc.pl.umap(adata, color = ['CRABP1','MAF','TAC3', 'LHX8','TH','CHRNA3', 'CHRNA5', 'CHRNA4', 'TRHDE', 'KIT', 'COL19A1', 'CHRNA7', 'ZIC1', 'ZIC2','ZIC4', 'CHAT', 'batch'],use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


#Filter data -- removing clusters where more than or equal to 3 genes have above average expression (marked as true)
adata = adata[boolmat2[adata.obs['leiden']].values <= 3]

# Print information about the filtering results
print(f"Number of cells after filtering: {len(adata)}")

print(adata.X)


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

#batch correction
sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']


#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["GRIA1", 'NRXN3', 'TPX2', 'CENPE', 'NRG3', 'ALCAM', 'DLX5', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)


#save inhibitory neurons as a separate file
adata.write('pathtopig/updated_genome/scRNAseq_pig_concat_inhibneurons_updatedrevisions.h5ad')


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



adata = sc.read_h5ad('pathtopig/updated_genome/scRNAseq_pig_concat_inhibneurons_updatedrevisions.h5ad')


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Integrated_pig" + str(resolution), frameon = False, legend_loc = "on data")


#subcluster TAC3 cluster so that it is only TAC3 cells
sc.tl.leiden(adata, restrict_to=('leiden', ['9']), resolution=0.2, key_added='leiden_subcluster_CRABPY')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY'], title = "Pig_subcluster" + str(resolution), frameon = False, legend_loc = "on data")


# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes', 'batch', 'timepoint'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")


plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['NEUROD2','SLA', 'SLC17A6', 'NEUROD6', 'NEUROG2', 'EOMES', 'GLI3','FOXG1', 'NKX2-1', 'LHX6', 'DLX1', 'DLX2', 'GAD1','CRABP1','COL19A1', 'CHRNA3', 'TH', 'PVALB', 'PTHLH', 'ETV1','MAF', 'TAC3', 'STXBP6','ZIC1','NXPH2', 'LHX8', 'RBP4', 'NR2F2', 'PROX1', 'CHAT', 'SST', 'NPY', 'CALB1', 'TH'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['HES1', 'HES5', 'PAX6', 'MEIS2', 'FOXP1', 'FOXP2', 'PENK', 'TSHZ1','ISL1', 'DRD1', 'DRD2', 'SCGN', 'NPY1R', 'PROX1', 'NR2F2', 'ZIC1', 'ZIC2', 'ZIC4', 'CHAT', 'GBX1', 'LAMP5', 'VIP'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


#sc._settings.settings._vector_friendly=True
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['LHX6', 'CRABP1', 'MAF', 'MEF2C','CHRNA3', 'ZIC1', 'ZIC4', 'TH', 'ZIC2','TAC3', 'NPY', 'SST', 'batch'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)



#Marmoset markers from 241213_TAC3_Update from Mike -- matching marmoset TAC3+/LHX8- population
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['CHRNA3', 'TAC3', 'LHX8', 'PLPP4', 'COL23A1', 'TRH', 'BMP7', 'PTPRK', 'SMOC1', 'SLIT2', 'TNS1', 'NPTX1', 'PDE3A', 'IGF1', 'PLCXD3',], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

sc.pl.umap(adata, color = ['LHX8', 'PTPRK','PTHLH', 'TNS1', 'ANGPT2', 'MAF','PDE3A', 'PLPP4', 'CHRNA3','COL19A1', 'TAC3', 'PTHLH','LHX8', 'CCK'], use_raw = False, frameon = False, color_map = "PuRd", size= 20, ncols = 3)

sc.pl.umap(adata, color = ['ENOX1', 'DDC', 'LHX8', 'HS6ST3', 'PLCB1', 'SLIT3', 'KLHL29', 'TRPC7', 'EPHA6', 'DPP10', 'PEG10'], use_raw = False, frameon = False, color_map = "PuRd", size= 20, ncols = 2)


#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata, groupby = 'leiden', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df


#Remove technical/low quality clusters
adata = adata[adata.obs['leiden_subcluster_CRABPY'].isin(['9,1', '15', '16', '17']) == False]

adata.layers["counts"]


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


#batch correction
sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["GRIA1", 'NRXN3', 'TPX2', 'CENPE', 'NRG3', 'ALCAM', 'DLX5', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)


#save inhibitory neurons as a separate file
adata.write('pathtopig/updated_genome/scRNAseq_pig_concat_inhibneurons_updatedrevisions.h5ad')


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

adata = sc.read_h5ad('pathtopig/updated_genome/scRNAseq_pig_concat_inhibneurons_updatedrevisions.h5ad')


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Integrated_pig" + str(resolution), frameon = False, legend_loc = "on data")



#subcluster 
sc.tl.leiden(adata, restrict_to=('leiden', ['7']), resolution=0.3, key_added='leiden_subcluster_CRABPY')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY'], title = "Pig_subcluster" + str(resolution), frameon = False, legend_loc = "on data")



#subcluster 
sc.tl.leiden(adata, restrict_to=('leiden_subcluster_CRABPY', ['2']), resolution=0.4, key_added='leiden_subcluster_CRABPY2')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY2'], title = "Pig_subcluster" + str(resolution), frameon = False, legend_loc = "on data")


plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['CHRNA3', 'TAC3', 'PTHLH', 'MAF', 'CRABP1', 'LHX8', 'TH', 'NPY1R', 'ISL1', 'PENK', 'VIP', 'TSHZ1', 'SCGN', 'NR2F2', 'PROX1'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['CHRNA3', 'TAC3', 'PTHLH', 'MAF', 'CRABP1', 'LHX8', 'TH', 'NPY1R', 'ISL1', 'PENK', 'VIP', 'TSHZ1', 'SCGN', 'NR2F2', 'PROX1'], use_raw = True, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes', 'batch', 'timepoint'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")


plt.rcParams["figure.figsize"] = (10, 10) 
sc.pl.umap(adata, color = ['NEUROD2','SLA', 'SLC17A6', 'NEUROD6', 'NEUROG2', 'EOMES', 'GLI3','FOXG1', 'NKX2-1', 'LHX6', 'DLX1', 'DLX2', 'GAD1','CRABP1','COL19A1', 'CHRNA3', 'TH', 'PVALB', 'PTHLH', 'ETV1','MAF', 'TAC3', 'STXBP6','ZIC1','NXPH2', 'LHX8', 'RBP4', 'NR2F2', 'PROX1', 'CHAT', 'SST', 'NPY', 'CALB1', 'TH'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


# Define a dictionary to map old cluster names to new names
cluster_name_mapping = {
    '7,1' : 'MGE_CRABP1/TAC3',
    '7,2' : 'MGE_CRABP1/TAC3',
    '7,3': 'MGE_CRABP1/TAC3',
    '7,5' : 'MGE_CRABP1/TAC3',


    
   
    '7,0': 'MGE_CRABP1/MAF',
    '2,2' : 'MGE_CRABP1/MAF',
    '2,3' : 'MGE_CRABP1/MAF',
    '5' : 'Progenitor',
    '13' : 'Progenitor',
    '16' : 'Progenitor',
    '6' : 'Progenitor',
    
    '11' : 'MGE_LHX6/MAF',
    '3' : 'MGE_LHX6/MAF',
    '2,0' : 'MGE_LHX6/MAF',

    '2,1' : 'MGE_LHX6/MAF',
    '2,5' : 'MGE_LHX6/MAF',
    '2,4' : 'MGE_LHX6/MAF',

    '14' : 'MGE_LHX6/MAF',
    '7,4' : 'MGE_LHX6/MAF',
 
    
    '0' : 'LGE_MEIS2/PAX6',
    '8' : 'LGE_MEIS2/PAX6',
    
    '12' : 'LGE_MEIS2/PAX6',
    '10': 'LGE_MEIS2/PAX6',
    '1': 'LGE_MEIS2/PAX6',


    '15' : 'LGE_FOXP1/PENK',
    
    
    '4' : 'CGE_NR2F2/PROX1',
    '17' : 'CGE_NR2F2/PROX1',
    
    '9' : 'LGE_FOXP2/TSHZ1',



    
 
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['leiden_subcluster_CRABPY2'] = adata.obs['leiden_subcluster_CRABPY2'].replace(cluster_name_mapping)

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
}



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['NKX2-1', 'CRABP1', 'MAF','CHRNA3', 'TAC3', 'TH'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 2, save='markers.svg')

plt.rcParams["figure.figsize"] = (10, 10)
sc.set_figure_params(fontsize = 20)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= ['batch'], frameon = False,size = 10,  save = '_batch.svg'

plt.rcParams["figure.figsize"] = (10, 10)
sc.set_figure_params(fontsize = 20)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= ['timepoint'], frameon = False,size = 10,  save = '_timepoint.svg')


# Rename 'cortex' to 'Cortex' and 'striatum' to 'Striatum' in adata.obs['region']
adata.obs['region'] = adata.obs['region'].replace({'cortex': 'Cortex', 'striatum': 'Striatum'})


# Define the custom colors for Cortex and Striatum
region_colors = {'Cortex': '#1d26cf', 'Striatum': '#f76f3e'}

plt.rcParams["figure.figsize"] = (10, 10)
sc.set_figure_params(fontsize = 20)
plt.rcParams["figure.figsize"] = (10, 10)

# Plot UMAP with custom colors for the 'region' column
sc.pl.umap(adata, color='region', frameon=False, size=10, save='_region.svg', palette=region_colors)


#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata, groupby = 'leiden_subcluster_CRABPY2', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df


adata.obs['leiden_subcluster_CRABPY2'] = adata.obs['leiden_subcluster_CRABPY2'].cat.reorder_categories(['Progenitor', 'MGE_LHX6/MAF', 'MGE_CRABP1/MAF', 'MGE_CRABP1/TAC3',
                                                                                                       'LGE_MEIS2/PAX6', 'LGE_FOXP2/TSHZ1', 'LGE_FOXP1/PENK', 'CGE_NR2F2/PROX1',
                                                                                                       ])



#Matrix of marker genes for inhibitory neuron clusters

#marker genes based on Schmitz et al. 2022
marker_genes_dict = {
    'Progenitor' : ['ASPM', 'CENPF', 'CKB'],
    'MGE_LHX6/MAF': ['LHX6', 'RBFOX1', 'SOX6', 'MEF2C'],
    'MGE_CRABP1/MAF': ['MAF','CRABP1','ETV1','GALNTL6','OLFM2',],
    'MGE_CRABP1/TAC3': ['NRTN', 'KIT','TRHDE','STXBP6','TAC3', 'CHRNA3', 'CHRNA7', 'PARD3',  'PDE1C', ],
    'LGE_MEIS2/PAX6': ['PAX6','MEIS2'],
    'LGE_FOXP2/TSHZ1': ['EYA2','FOXP2', 'FOXP4','TSHZ1'],
    'LGE_FOXP1/PENK': ['FOXP1', 'SIX3', 'PENK'],
    'CGE_NR2F2/PROX1': ['PDZRN3', 'NR3C2', 'NR2F2', 'PROX1'],
    
}

sc.pl.matrixplot(adata, groupby='leiden_subcluster_CRABPY2', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes=True,save='matrix.svg')
# swap_axes=True
# 
#
# , 

# Compute proportions of Cortex vs. Striatum for each cluster
cluster_region_counts = adata.obs.groupby(['leiden_subcluster_CRABPY2', 'region']).size().unstack(fill_value=0)
cluster_region_props = cluster_region_counts.div(cluster_region_counts.sum(axis=1), axis=0)

# Ensure order matches matrixplot
cluster_order = adata.obs['leiden_subcluster_CRABPY2'].cat.categories
cluster_region_props = cluster_region_props.loc[cluster_order]

# Set up figure with two subplots: one for barplot, one for matrixplot
fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 10]}, figsize=(4.5, 14))

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
    swap_axes=True,
    show=False,
    ax=ax[1]  # Draw on second subplot
)

plt.savefig("proportionmatrix.svg")
plt.tight_layout()
plt.show()


#Matrix of marker genes for inhibitory neuron clusters
sc.set_figure_params(fontsize = 12)
# Compute proportions of Cortex vs. Striatum for each cluster
cluster_region_counts = adata.obs.groupby(['leiden_subcluster_CRABPY2', 'region']).size().unstack(fill_value=0)
cluster_region_props = cluster_region_counts.div(cluster_region_counts.sum(axis=1), axis=0)

# Ensure order matches matrixplot
cluster_order = adata.obs['leiden_subcluster_CRABPY2'].cat.categories
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


# Filter the AnnData object to include only clusters that start with '0840', '0841', or '0842'
clusters_of_interest = adata.obs['leiden_subcluster_CRABPY2'].isin(['MGE_CRABP1/MAF','MGE_CRABP1/TAC3'])
adata_filtered = adata[clusters_of_interest]

# List of genes of interest
genes_of_interest = ['LHX6', 'CRABP1','ANGPT2', 'MAF','CHL1','RBP4', 'ARX', 'TAC3', 'TRH', 'CHRNA3', 'CHRNA7','ZIC1', 'LHX8', 'TH', 'PTHLH',  'PARD3', 'OLFM2', 'PDE1C']




# Create a dot plot for the specified genes
sc.pl.dotplot(adata_filtered, var_names=genes_of_interest, use_raw= False, vmin=-2, vmax=2, swap_axes = True, cmap='PiYG_r',groupby='leiden_subcluster_CRABPY2', save = 'macdotplot.svg')


#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata_filtered, groupby = 'region', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata_filtered.uns['rank_genes_groups']['names']).head(20)
pd.set_option('display.max_columns', 500)
small_degs_df

# Filter the AnnData object to include only clusters 
clusters_of_interest = adata.obs['leiden_subcluster_CRABPY2'].isin(['MGE_CRABP1/TAC3'])
adata_filtered = adata[clusters_of_interest]

marker_genes_dict = ['LHX6', 'NKX2-1', 'SOX6', 'MAF', 'CHL1','CRABP1','ETV1','GALNTL6','OLFM2','NRTN', 'KIT','TRHDE','STXBP6','TAC3', 'CHRNA3', 'CHRNA7', 'PARD3',  'PDE1C', 'ZIC1', 'ZIC2', 'ZIC4', 'LHX8', 'CXCR4', 'UNC5C', 'HCN1', 'DBI', 'ADCY8', 'PCSK5']


sc.pl.dotplot(adata_filtered, groupby='region', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes=True,  figsize=(2.5,7),save='matrix.svg')
adata_filtered.obs['class_region'] = (
    adata_filtered.obs['class'].astype(str) + "_" + adata_filtered.obs['region'].astype(str)
)



sc.set_figure_params(fontsize = 14)
# First create the combined class_region column
adata.obs['class_region'] = adata.obs['class'].astype(str) + "_" + adata.obs['region'].astype(str)

# Now you can filter
clusters_of_interest = adata.obs['class'].isin(['MGE_CRABP1/TAC3', 'MGE_CRABP1/MAF'])
adata_filtered = adata[clusters_of_interest]

marker_genes_dict = {
    'MGE_CRABP1/MAF': ['LHX6', 'NKX2-1', 'SOX6', 'CRABP1','ETV1','GALNTL6','OLFM2', 'MAF', 'MAFB', 'COL19A1', 'RBP4'],
    'MGE_CRABP1/TAC3': ['NRTN', 'KIT','TRHDE','STXBP6','TAC3', 'CHRNA3','PARD3',  'PDE1C', 'LHX8', 'ZIC1', 'ZIC2', 'ZIC4'],
    'Cortical': ['CXCR4', 'HCN1', 'UNC5C','DACH1','DBI', 'TENM3'],
    'Striatal': ['ADCY8', 'PCSK5', 'SPATS2L']}



sc.pl.matrixplot(adata_filtered, groupby='class_region', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes=False,  save='matrixcortexstriatum.svg')


# Get counts of cells per region
region_counts = adata_filtered.obs['region'].value_counts()

# Display counts
print(region_counts)


adata.obs['class'].value_counts().get('MGE_CRABP1/TAC3', 0)

plt.rcParams["figure.figsize"] = (10, 10)
sc.set_figure_params(fontsize = 20)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= ['timepoint'], frameon = False,size = 10)


# Copy the 'leiden_subcluster_CRABPY2' observation to a new 'class' observation
adata.obs['class'] = adata.obs['leiden_subcluster_CRABPY2']


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
obs.to_csv('pathtogeo/pig_barcodes.tsv', sep='\t', index=True)

# Save var.tsv
var.to_csv('pathtogeo/pig_genes.tsv', sep='\t', index=True)

# Save matrix.mtx
scipy.io.mmwrite('pathtogeo/pig_matrix.mtx', matrix)




# # Differential gene expression analysis of cortical and striatal TAC3 populations


import scanpy as sc
import decoupler as dc
import pandas as pd

# Needed for some plotting
import matplotlib.pyplot as plt


pdata = dc.get_pseudobulk(
    adata,
    sample_col='batch',
    groups_col='class',
    layer='counts',
    mode='sum',
    min_cells=0,
    min_counts=0
)

dc.plot_psbulk_samples(pdata, groupby=['batch', 'class'], figsize=(12, 4))



# Get filtered pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col='batch',
    groups_col='class',
    layer='counts',
    mode='sum',
    min_cells=10,
    min_counts=1000
)
pdata


pdata.X
# Store raw counts in layers
pdata.layers['counts'] = pdata.X.copy()

# Normalize, scale and compute pca
sc.pp.normalize_total(pdata, target_sum=1e4)
sc.pp.log1p(pdata)
sc.pp.scale(pdata, max_value=10)
sc.tl.pca(pdata)

# Return raw counts to X
dc.swap_layer(pdata, 'counts', X_layer_key=None, inplace=True)



sc.pl.pca(pdata, color=['batch', 'class'], ncols=1, size=300)
sc.pl.pca_variance_ratio(pdata)


dc.get_metadata_associations(
    pdata,
    obs_keys = ['timepoint', 'region', 'psbulk_n_cells', 'psbulk_counts'],  # Metadata columns to associate to PCs
    obsm_key='X_pca',  # Where the PCs are stored
    uns_key='pca_anova',  # Where the results are stored
    inplace=True,
)



dc.plot_associations(
    pdata,
    uns_key='pca_anova',  # Summary statistics from the anova tests
    obsm_key='X_pca',  # where the PCs are stored
    stat_col='p_adj',  # Which summary statistic to plot
    obs_annotation_cols = ['timepoint', 'region', 'class'], # which sample annotations to plot
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

import pydeseq2
print(pydeseq2.__version__)

# Build DESeq2 object
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    adata=TAC3,
    design_factors='region',
    ref_level=['region', 'Cortex'],
    refit_cooks=True,
    inference=inference,
)


# Compute LFCs
dds.deseq2()


# Extract contrast 
stat_res = DeseqStats(
    dds,
    contrast=["region", 'Cortex', 'Striatum'],
    inference=inference,
)


# Compute Wald test
stat_res.summary()


# Extract results
results_df = stat_res.results_df
results_df


sc.set_figure_params(fontsize = 20)

# Generate the plot without saving
dc.plot_volcano_df(
    results_df,
    x='log2FoldChange',
    y='padj',
    top=20,
    figsize=(8, 8),
    color_pos='#1d26cf', 
    color_neg='#f76f3e', 
    dpi=500,
    save=None
)

# Get current Axes
ax = plt.gca()

# Make the dots bigger
for coll in ax.collections:
    if hasattr(coll, "set_sizes"):
        coll.set_sizes([100])  # Try 60–150 for noticeable effect

# Optional: also tweak linewidths if you still want
for line in ax.lines:
    line.set_linewidth(0.25)

# Save manually
plt.savefig('volcanotac3.svg')



# Select cell profiles
MAF = pdata[pdata.obs['class'] == 'MGE_CRABP1/MAF'].copy()


dc.plot_filter_by_expr(MAF, group='region', min_count=10, min_total_count=15)


#filter if needed
# Obtain genes that pass the thresholds
genes = dc.filter_by_expr(MAF, group='region', min_count=10, min_total_count=15)


# Filter by these genes
MAF = MAF[:, genes].copy()
MAF

# Import DESeq2
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats

# Build DESeq2 object
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(
    adata=MAF,
    design_factors='region',
    ref_level=['region', 'Cortex'],
    refit_cooks=True,
    inference=inference,
) 
# Compute LFCs
dds.deseq2()


# Extract contrast between COVID-19 vs normal
stat_res = DeseqStats(
    dds,
    contrast=["region", 'Cortex', 'Striatum'],
    inference=inference,
)

# Compute Wald test
stat_res.summary()

# Extract results
results_dfmaf = stat_res.results_df
results_dfmaf


dc.plot_volcano_df(
    results_dfmaf,
    x='log2FoldChange',
    y='padj',
    top=20,
    figsize=(15, 15),
    color_pos='#1d26cf', 
    color_neg='#f76f3e', 
    save = 'volcanomaf.svg'
    
)


# Generate the plot without saving
dc.plot_volcano_df(
    results_dfmaf,
    x='log2FoldChange',
    y='padj',
    top=10,
    figsize=(12, 12),
    color_pos='#1d26cf', 
    color_neg='#f76f3e', 
    dpi=500,
    save=None
)

# Get current Axes
ax = plt.gca()

#  Make the dots bigger
for coll in ax.collections:
    if hasattr(coll, "set_sizes"):
        coll.set_sizes([100])  # Try 60–150 for noticeable effect

# Optional: also tweak linewidths if you still want
for line in ax.lines:
    line.set_linewidth(0.25)

# Save manually
plt.savefig('volcanomaf.svg')



# # Subset data for trying to integrate with other species


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



adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()


adata.X = np.array(adata.X)  # Convert to numpy array
adata.write('pathtopig/updated_genome/pig_subsetraw1k.h5ad')

