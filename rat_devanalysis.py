#!/usr/bin/env python
# coding: utf-8

### Preprocessing by region
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
adata=sc.read_10x_h5('pathtorat/cellbender/RatE18striatum_cellbender_output_filtered.h5')
Results_file='pathtorabbit/scRNAseq_RatE18striatum.h5ad'
adata

adata.obs_names_make_unique()
adata.var_names_make_unique()

#h5ad stores whole anndata data structure
adata.write(Results_file)

#Define string to append to saved figures
save_name = "rat_striatum_E18"


# ## Preprocessing of data
#Basic filtering

sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=1)

# Annotate varibles for mitochondrial and ribosomal genes

mito_genes=[name for name in adata.var_names if name in ['Nd1','Nd2','Nd4l','Nd4','Nd5','Nd6','Atp6','Atp8','Cytb','Cox1','Cox2','Cox3'] or 'Mt-' in name]
ribo_genes=[name for name in adata.var_names if name.startswith('Rps') or name.startswith('Rpl') ]

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
adata = adata[adata.obs.percent_mito < .2, :]
adata = adata[adata.obs.n_genes < 6000, :]
adata = adata[adata.obs.n_genes > 1500, :]


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

adata.write('pathtorat/scRNAseq_striatum_rat_E18_prenorm.h5ad')




# ## Rat with cellbender prefiltered concat analysis
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

path = '/pathto/rat/'

Striatum = sc.read_h5ad(path + 'scRNAseq_striatum_rat_E18_prenorm.h5ad')
Cortex = sc.read_h5ad(path + 'scRNAseq_cortex_rat_E18_prenorm.h5ad')

for gem in [Striatum, Cortex]:
    print(gem.shape[0])

adata = anndata.AnnData.concatenate(Striatum, Cortex, join='outer', batch_categories=['Striatum', 'Cortex'])

adata.obs_names_make_unique()
adata.var_names_make_unique()

Results_file='/pathto/rat/rat_concat.h5ad'
adata.write(Results_file)


#Define string to append to saved figures
save_name = "rat_concat"


#violin plots of the computed quality measures

sc.pl.violin(adata,groupby='batch',keys=['percent_ribo','percent_mito', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True, )

#Set the .raw attribute of the AnnData object
adata.raw = adata


# ## Normalize counts scanpy tutorial
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
pca_var_genes = ["Vim", "Mef2c", "Myo16", "Dlx6", "Dlx5",
                 "Gad2"]
sc.pl.pca(adata, color = pca_var_genes, frameon = True)

adata.write(Results_file)

#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

#Save the AnnData object as an .h5ad file.
adata.write(Results_file)


# ## Plot UMAP
adata = sc.read_h5ad(Results_file)

# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams ["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Rat_concat" + str(resolution), frameon = False, legend_loc = "on data")

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2, save = 'umapconcat_qc.png')

#Save the AnnData object as an .h5ad file.
adata.write(Results_file)

sc.pl.umap(adata, color = ['batch'], wspace=0.25, ncols = 2)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Foxg1', 'Gad1', 'Dlx1', 'Dlx2', 'Mki67','Nkx2-1', 'Lhx6', 'Crabp1','Etv1','Tac3', 'Lhx8','Tacr3','Maf','Chrna3', 'Chrna4', 'Chrna7','Rbp4', 'Nxph1', 'Nxph2', 'Cer1', 'Th','Mef2c'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3, save = 'umapgenes_concat1')

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Nkx2-1', 'Lhx6', 'Crabp1','Tac3', 'Chrna3', 'Th', 'Rbp4', 'Maf','Sst', 'Npy'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3, save = 'umapgenes_concat1')

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Mki67','Nr2f2','Meis2','St18','Pax6', 'Tbr1', 'Satb2', 'Cux2', 'Eomes', 'Hes1','Olig2', 'Gfap', 'Egfr', 'Erbb4', 'Zic1', 'Zic2'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3, save= 'moremarkers')


# ## Isolate inhibitory neurons
inhibmat = adata[:, ['Dlx1','Dlx2', 'Dlx5', 'Dlx6', 'Gad1', 'Gad2']].to_df()
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
adata.write('/pathto/rat/rat_concatinhib.h5ad')


# # Renormalize, Recalculate PCA, recalculate UMAP


#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()

#Save the AnnData object as an .h5ad file.
adata.write('/pathto/rat/rat_concatinhib_raw.h5ad')

adata = sc.read_h5ad('/pathto/rat/rat_concatinhib_raw.h5ad')

# Directly access the raw count matrix - calling it "processed counts" because later .X will store your processed data, while the "counts" layer will 
processed_counts = adata.X
print(processed_counts)


adata.layers["counts"] = adata.X.copy() #Copy the raw data into a "layer" in the AnnData structure - this will be used for normalization, etc.

adata.layers["counts"] #This should resemble the raw data at this point

#scanpy tutorial
sc.pp.normalize_total(adata, target_sum=1e4) #Specify "layer = "counts" to use the counts layer for normalization
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='batch', min_mean=0.0125, max_mean=3, min_disp=0.5)
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
pca_var_genes = ["Cenpf", "Top2a", "Nyap2",
                 "Gad1", "Dcc", "Chn2",
                 "Cux2"]
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")


#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

#Save the AnnData object as an .h5ad file.
adata.write('/pathto/rat/rat_concatinhib.h5ad')


# # Load inhibitory neuron dataset and plot
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

sc._settings.settings._vector_friendly=True

adata = sc.read_h5ad('/pathto/rat/rat_concatinhib.h5ad')

# Cluster umap embeddings using leiden and save umap plots.
resolution = 2 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Rat_Inhibitory_Neurons", frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')

excitatorymat = adata[:, ['Sla', 'Slc17a6', 'Neurod2', 'Neurog2', 'Neurod6', 'Eomes']].to_df()
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


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)
# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Rat_Inhibitory_Neurons", frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2, save = 'umapconcat_qc.png')
#sc.pl.violin(adata, ['pct_counts_ribo','pct_counts_mito'], groupby = "leiden")


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Dlx1', 'Dlx2', 'Gad1', 'Gad2', 'Sla', 'Neurod2', 'Slc17a6'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)
#excitatory neuron markers 'Sla', 'Neurod2', 'Slc17a6'


# # Remove low quality clusters
adata = adata[adata.obs['leiden'].isin(['15']) == False]

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


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")


#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

adata.write('/pathto/rat/rat_concatinhib.h5ad')


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

adata = sc.read_h5ad('/pathto/rat/rat_concatinhib.h5ad')

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
sc.pl.umap(adata, color= 'leiden', title = "Integrated_rat" + str(resolution), frameon = False, legend_loc = "on data")

#subcluster 
sc.tl.leiden(adata, restrict_to=('leiden', ['6']), resolution=0.2, key_added='leiden_subcluster_CRABPY')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY'], title = "Rat_subcluster" + str(resolution), frameon = False, legend_loc = "on data")

#subcluster
sc.tl.leiden(adata, restrict_to=('leiden_subcluster_CRABPY', ['4']), resolution=0.2, key_added='leiden_subcluster_CRABPY2')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY2'], title = "Rat_subcluster" + str(resolution), frameon = False, legend_loc = "on data")

#subcluster
sc.tl.leiden(adata, restrict_to=('leiden_subcluster_CRABPY2', ['9']), resolution=0.2, key_added='leiden_subcluster_CRABPY3')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY3'], title = "Rat_subcluster", frameon = False, legend_loc = "on data")

sc._settings.settings._vector_friendly=True

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','percent_mito', 'n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2)
sc.pl.violin(adata, ['percent_ribo','percent_mito'], groupby = "leiden")

# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['batch'], wspace=0.25, ncols = 2, frameon = False, size= 10, save = 'batch.svg')

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Nkx2-1', 'Mef2c', 'Etv1', 'Maf', 'Meis2', 'Pax6', 'Foxp1', 'Foxp2', 'Nr2f2', 'Lhx8', 'Dlx1', 'Gad1', 'Hes1', 'Spc24', 'Cdc20' ,'Top2a'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 4)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Foxg1', 'Gad1', 'Dlx1', 'Dlx2', 'Mki67','Nkx2-1', 'Lhx6', 'Crabp1','Etv1', 'Zic1', 'Zic2', 'Zic4', 'Gbx1','Tac3','Th', 'Chrna3', 'Lhx8', 'Chat','Maf','Chrna3','Rbp4'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Chat', 'Meis2', 'Pax6', 'Foxp1', 'Foxp2', 'Isl1', 'Penk', 'Npy1r', 'Tshz1', 'Scgn', 'Prox1', 'Nr2f2' ], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Npy', 'Mki67','Nr2f2', 'Prox1','Meis2','Pax6','St18','Pax6', 'Cux2', 'Hes1','Olig2', 'Egfr', 'Erbb4', 'Zic1', 'Zic2', 'Zic4', 'Nxph1', 'Isl1', 'Penk', 'Reln'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Satb2', 'Neurod6', 'Neurod2', 'Tfap2d', 'Tbr1', 'Nr2f2', 'Ascl1', 'Gadd45g', 'Foxg1', 'Emx2', 'Lhx2', 'Gsx1', 'Gsx2', 'Npy' ,'Nhlh1', 'Tp73','Lhx1', 'Lamp5'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3)

sc.set_figure_params(fontsize = 50)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Nkx2-1', 'Crabp1', 'Maf','Chrna3', 'Tac3', 'Th' ], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 6, save='keymarkersumap.svg')

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Nkx2-1', 'Crabp1', 'Rbp4','Maf', 'Col19a1','Chrna3', 'Tac3', 'Th', 'Lhx8' ], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3, save = 'crabp1markers')

#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata, groupby = 'leiden_subcluster_CRABPY3', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df

sc.pl.umap(adata, color= 'leiden_subcluster_CRABPY3', title = "Rat_Inhibitory_Neurons", frameon = False, legend_loc = "on data")


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


# Define a dictionary to map old cluster names to new names
cluster_name_mapping = {
 
    '19' : 'MGE_CRABP1/MAF',
    
    
    '14' : 'Progenitor',
    '1' : 'Progenitor',
    '8' : 'Progenitor',
    '18' : 'Progenitor',
    '2' : 'Progenitor',
    '20' : 'Progenitor',
    '7' : 'Progenitor',
    '15' : 'Progenitor',


    
    '12' : 'MGE_LHX6/MAF',
    '6,0' : 'MGE_LHX6/MAF',
    '3' : 'MGE_LHX6/MAF',
    '22' : 'MGE_LHX6/MAF',
    
    '6,1' : 'MGE_LHX6/NPY',
        
    '5' : 'CGE_NR2F2/PROX1',
    


    '4,0' : 'LGE_FOXP1/PENK',
    
    '10' : 'LGE_FOXP1/ISL1',
    '13' : 'LGE_FOXP1/ISL1',
    '17' : 'LGE_FOXP1/ISL1',
    '0' : 'LGE_FOXP1/ISL1',

    '11' : 'LGE_FOXP1/ISL1/NPY1R',
    

    '9,1' : 'LGE_FOXP2/TSHZ1',
    '16' : 'LGE_FOXP2/TSHZ1',
    
    '4,2' : 'LGE_MEIS2/PAX6',
    '9,0' : 'LGE_MEIS2/PAX6',
    '4,1' : 'LGE_MEIS2/PAX6',
    
    
    '21' : 'VMF_CRABP1/LHX8'

    
 
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['leiden_subcluster_CRABPY3'] = adata.obs['leiden_subcluster_CRABPY3'].replace(cluster_name_mapping)

paldict={'CGE_NR2F2/PROX1': 'mediumpurple',

    'LGE_FOXP1/ISL1': 'royalblue',
    'LGE_FOXP1/ISL1/NPY1R': 'blue',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_MEIS2/PAX6/SCGN': 'orange',
    'LGE_MEIS2/PAX6': 'orangered',
    'LGE_FOXP2/TSHZ1': '#17344c',
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/CHRNA3': 'fuchsia',
    'MGE_LHX6/MAF': 'skyblue',
    'MGE_LHX6/NPY': 'teal',
    'Progenitor': 'pink',

    'Transition': '#3c7632',
    'VMF_ZIC1/ZIC2': 'green',
    'VMF_CRABP1/LHX8': 'mediumspringgreen',

}

adata.obs['leiden_subcluster_CRABPY3'] = adata.obs['leiden_subcluster_CRABPY3'].cat.reorder_categories(['Progenitor', 'MGE_LHX6/MAF', 'MGE_LHX6/NPY', 'MGE_CRABP1/MAF',
                                                                                                       'LGE_MEIS2/PAX6','LGE_FOXP2/TSHZ1', 'LGE_FOXP1/PENK', 
                                                                                                        'LGE_FOXP1/ISL1', 'LGE_FOXP1/ISL1/NPY1R', 'CGE_NR2F2/PROX1',
                                                                                                         'VMF_CRABP1/LHX8'])

# Copy the 'leiden_subcluster_CRABPY2' observation to a new 'class' observation
adata.obs['class'] = adata.obs['leiden_subcluster_CRABPY3']


#Matrix of marker genes for inhibitory neuron clusters

#marker genes based on Schmitz et al. 2022
marker_genes_dict = {
    'Progenitor' : ['Aspm', 'Cenph', 'Mki67'],
    'MGE_LHX6/MAF': ['Lhx6', 'Sox6', 'Mef2c','Cux2'],
    'MGE_LHX6/NPY': ['Sst','Npy'],
    'MGE_CRABP1/MAF': ['Crabp1', 'Etv1','Maf', 'Rbp4', 'Col19a1'],
    'LGE_MEIS2/PAX6': ['Pax6', 'Meis2'],
    'LGE_FOXP2/TSHZ1': ['Foxp2', 'Tshz1'],
    'LGE_FOXP1/PENK': ['Foxp1', 'Six3', 'Penk'],
    'LGE_FOXP1/ISL1': ['Rxrg', 'Rarb', 'Pbx3','Isl1'],
    'LGE_FOXP1/ISL1/NPY1R': [ 'Npy1r'],
    'CGE_NR2F2/PROX1' : ['Nr2f2', 'Prox1', 'Sp8', 'Nr3c2', 'Pdzrn3'],
    'VMF_CRABP1/LHX8' : ['Lhx8', 'Gbx1', 'Chat']
    
}


sc.pl.matrixplot(adata, groupby='class', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes=True,save='matrix.svg')


#make a region column
adata.obs['region'] = adata.obs['batch']

# Compute proportions of Cortex vs. Striatum for each cluster
cluster_region_counts = adata.obs.groupby(['leiden_subcluster_CRABPY3', 'region']).size().unstack(fill_value=0)
cluster_region_props = cluster_region_counts.div(cluster_region_counts.sum(axis=1), axis=0)

# Ensure order matches matrixplot
cluster_order = adata.obs['leiden_subcluster_CRABPY3'].cat.categories
cluster_region_props = cluster_region_props.loc[cluster_order]

#Reorder columns so 'Cortex' comes first
cluster_region_props = cluster_region_props[['Cortex', 'Striatum']]
# Set up figure with two subplots: one for barplot, one for matrixplot
fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 10]}, figsize=(20, 10))

# Top subplot: stacked barplot for cortex vs. striatum proportions
cluster_region_props.plot(kind='bar', stacked=True, color=['#1d26cf','#f76f3e', ], ax=ax[0])
# Remove the grid lines
ax[0].grid(False)
# Format barplot
ax[0].set_xticks([])
ax[0].set_ylabel('Proportion')
ax[0].set_title('Cortex vs. Striatum Proportions per Cluster')
ax[0].legend(
    title="Region", 
    labels=['Cortex','Striatum', ], 
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


adata.write('/pathto/rat/rat_devinhibneurons.h5ad')

adata.write('/pathto/rat/rat_processed.h5ad')

adata = sc.read_h5ad('/pathto/rat/rat_processed.h5ad')

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
obs.to_csv('/pathto/geo/rat_barcodes.tsv', sep='\t', index=True)

# Save var.tsv
var.to_csv('/pathto/geo/rat_genes.tsv', sep='\t', index=True)

# Save matrix.mtx
scipy.io.mmwrite('/pathto/geo/rat_matrix.mtx', matrix)

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
adata.write('/pathto/rat/rat_subsetraw1k.h5ad')

