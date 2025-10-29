#!/usr/bin/env python
# coding: utf-8

### Preprocessing by timepoint
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
adata=sc.read_10x_h5('/pathto/mouse/Cellbender_out/E15/mouseE15_cellbender_output_filtered.h5')
Results_file='/pathto/mouse/mouse_tac2E15cellbender.h5ad'
adata


adata.obs_names_make_unique()
adata.var_names_make_unique()

#h5ad stores whole anndata data structure
adata.write(Results_file)


#Define string to append to saved figures
save_name = "E15mouse_striatum"


# ## Preprocessing of data

#Basic filtering

sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=1)


# Annotate varibles for mitochondrial and ribosomal genes

mito_genes=[name for name in adata.var_names if name in ['Nd1','Nd2','Nd4l','Nd4','Nd5','Nd6','Atp6','Atp8','Cytb','Cox1','Cox2','Cox3'] or 'mt-' in name]
ribo_genes=[name for name in adata.var_names if name.startswith('Rps') or name.startswith('Rpl') ]

# for each cell compute fraction of counts in mito genes vs. all genes

adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
adata.obs['percent_ribo'] = np.sum(
adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)

# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1)



#violin plots of the computed quality measures

sc.pl.violin(adata,keys=['percent_ribo','percent_mito', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True, 
             save= save_name + 'QC_metrics')


sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes')



#Histograms of the distrubution of UMI counts and gene numbers

Hist2 = seaborn.distplot(adata.obs['n_genes'], kde=False)
Hist2.set_xlabel("Number of genes", fontsize=12)
Hist2.set_ylabel("Frequency", fontsize=12)
#Hist2.axvline(1000, 0,1, color='red')
plt.show()

Hist3 = seaborn.distplot(adata.obs['n_genes'][adata.obs['n_genes']<4000], kde=False, bins=60)
Hist3.set_xlabel("n_genes", fontsize=12)
Hist3.set_ylabel("Frequency", fontsize=12)
#Hist3.axvline(1000, 0,1, color='red')

plt.show()

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
adata = adata[adata.obs.n_genes > 2250, :]
#adata = adata[adata.obs.n_genes < 6000, :]


# ## Review data after filtering


sc.pl.violin(adata, keys=['percent_ribo','percent_mito', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True, save= save_name + 'QC_metrics')



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


adata.write('/pathto/mouse/E15_prenorm.h5ad')



# ## Analyzing timepoints together
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

path = '/pathto/'

E15 = sc.read_h5ad(path + 'E15_prenorm.h5ad')
E17 = sc.read_h5ad(path + 'E17_prenorm.h5ad')
E18 = sc.read_h5ad(path + 'E18_prenorm.h5ad')



for gem in [E15, E17, E18]:
    print(gem.shape[0])


adata = anndata.AnnData.concatenate(E15, E17, E18, join='outer', batch_categories=['E15', 'E17', 'E18'])

adata.obs_names_make_unique()
adata.var_names_make_unique()


#add in timepoint
adata.obs['timepoint'] = 'E15'
E17 = adata.obs['batch'].str.contains('E17')
adata.obs.loc[E17, 'timepoint'] = 'E17'
E18 = adata.obs['batch'].str.contains('E18')
adata.obs.loc[E18, 'timepoint'] = 'E18'



Results_file='/pathto/mouse/mouse_concat.h5ad'


#saving all the datasets together
adata.write(Results_file)


#Define string to append to saved figures
save_name = "mouse_concat"


#violin plots of the computed quality measures

sc.pl.violin(adata,groupby='batch',keys=['percent_ribo','percent_mito', 'n_counts', 'n_genes'],  rotation=90, multi_panel=True, )


adata = adata[adata.obs.n_genes < 7000, :]


#Set the .raw attribute of the AnnData object
adata.raw = adata


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
pca_var_genes = ["Vim", "Mef2c", "Myo16", "Dlx6", "Dlx5",
                 "Gad2"]
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']



#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")


#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams ["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Mouse_concat" + str(resolution), frameon = False, legend_loc = "on data")


# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo','n_counts', 'n_genes', 'timepoint'], wspace=0.25, ncols = 2, save = 'umapconcat_qc.png')


sc.pl.umap(adata, color = ['batch'], wspace=0.25, ncols = 2)



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Foxg1', 'Gad1', 'Dlx1', 'Dlx2', 'Mki67','Nkx2-1', 'Lhx6', 'Crabp1','Etv1','Tac2', 'Lhx8','Tacr3','Maf','Chrna3', 'Chrna4', 'Chrna7','Rbp4', 'Nxph1', 'Nxph2', 'Cer1', 'Th','Mef2c'], frameon = False, color_map = "PuRd", size= 25, ncols = 3, save = 'umapgenes_concat1')


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Nkx2-1', 'Lhx6', 'Crabp1','Tac2', 'Chrna3', 'Th', 'Rbp4', 'Maf','Sst', 'Npy'], frameon = False, color_map = "PuRd", size= 25, ncols = 3, save = 'umapgenes_concat1')

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Mki67','Nr2f2','Meis2','St18','Pax6', 'Tbr1', 'Satb2', 'Cux2', 'Eomes', 'Hes1','Olig2', 'Gfap', 'Egfr', 'Erbb4', 'Zic1', 'Zic2'], frameon = False, color_map = "PuRd", size= 25, ncols = 3, save= 'moremarkers')


# ## Isolate inhibitory neurons

inhibmat = adata[:, ['Dlx1','Dlx2', 'Dlx5', 'Dlx6', 'Gad1', 'Gad2']].to_df()
inhibmat['leiden'] = adata.obs['leiden']



meanmat = inhibmat.groupby('leiden').mean()
print(meanmat>meanmat.mean(0))
boolmat = (meanmat > meanmat.mean(0)).sum(1)


print(boolmat.index[boolmat>=3])



adata.obs

#Filter data
adata = adata[boolmat[adata.obs['leiden']].values >= 3]

# Print information about the filtering results
print(f"Number of cells after filtering: {len(adata)}")


#Save the AnnData object as an .h5ad file.
adata.write('/pathto/mouse/mouse_concatinhib.h5ad')


# # Renormalize, Recalculate PCA, recalculate UMAP


#reset data to raw
print(adata,flush=True)
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()


#Save the AnnData object as an .h5ad file.
adata.write('/pathto/mouse/mouse_concatinhib_raw.h5ad')

adata = sc.read_h5ad('/pathto/mouse/mouse_concatinhib_raw.h5ad')

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
pca_var_genes = ["Cenpf", "Top2a", "Nyap2",
                 "Gad1", "Dcc", "Chn2",
                 "Cux2"]
sc.pl.pca(adata, color = pca_var_genes, frameon = True)



sc.external.pp.harmony_integrate(adata, 'batch')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")

#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)

#Save the AnnData object as an .h5ad file.
adata.write('/pathto/mouse/mouse_concatinhib.h5ad')


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

adata = sc.read_h5ad('/pathto/mouse/mouse_concatinhib.h5ad')

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
resolution = 0.5 #Adjust resolution
#sc.tl.leiden(adata, resolution = resolution)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Mouse_Inhibitory_Neurons", frameon = False, legend_loc = "on data", save = 'leiden_clusters.svg')

#subcluster TAC3 cluster so that it is only TAC3 cells
sc.tl.leiden(adata, restrict_to=('leiden', ['5']), resolution=0.1, key_added='leiden_subcluster_CRABPY')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY'], title = "Mouse_Inhibitory_Neuorns_subcluster" + str(resolution), frameon = False, legend_loc = "on data")


#subcluster TAC3 cluster so that it is only TAC3 cells
sc.tl.leiden(adata, restrict_to=('leiden_subcluster_CRABPY', ['10']), resolution=0.1, key_added='leiden_subcluster_CRABPY2')
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY2'], title = "Mouse_Inhibitory_Neuorns_subcluster" + str(resolution), frameon = False, legend_loc = "on data")


# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['percent_ribo', 'percent_mito','n_counts', 'n_genes', 'batch'], wspace=0.25, ncols = 2, save ='inhib_qcumap')
#sc.pl.violin(adata, ['pct_counts_ribo','pct_counts_mito'], groupby = "leiden")

plt.rcParams["figure.figsize"] = (10, 10)
sc.set_figure_params(fontsize = 20)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'batch', title = "Mouse_Inhibitory_Neurons" , size=10, frameon = False, save = 'batch.svg')

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Nkx2-1', 'Mef2c', 'Etv1', 'Maf', 'Meis2', 'Pax6', 'Foxp1', 'Foxp2', 'Nr2f2', 'Lhx8', 'Dlx1', 'Gad1', 'Hes1', 'Spc24', 'Cdc20' ,'Top2a'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 4, save = 'umapgenes_concat1')


sc.set_figure_params(fontsize = 50)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Crabp1', 'Maf', 'Mafb', 'Col19a1', 'Rbp4',  'Stxbp6', 'Trhde','Chrna3', 'Th', 'Tac2'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 5, save = 'umaprevisions.svg')


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Nkx2-1', 'Crabp1', 'Maf','Chrna3', 'Th', 'Tac2'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 2, save = 'umaprevisions')



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Lhx8','Zic1','Lhx6', 'Crabp1', 'Chrna3', 'Th', 'Tac2', 'Ddc', 'Dbh', 'Ddc', 'Slc18a2', 'Slc6a3','Drd2', 'Gch1', 'Gabrq'], use_raw = False, frameon = False, color_map = "PuRd", size= 50, ncols = 3)


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['Npy', 'Mki67','Nr2f2','Meis2','St18','Pax6', 'Cux2', 'Hes1','Olig2', 'Egfr', 'Erbb4', 'Zic1', 'Zic2', 'Isl1', 'Penk', 'Reln', 'Col19a1'], use_raw = False, frameon = False, color_map = "PuRd", size= 25, ncols = 3, save= 'moremarkers')

#Remove unwanted clusters.
CRABP1cluster= adata[adata.obs['leiden_subcluster_CRABPY'].isin(['5,1']),:]
CRABP1cluster

# Get the Tac4 expression data for CRABP1 cluster
tac2_expression = CRABP1cluster[:, CRABP1cluster.var_names == 'Tac2'].X

# Count how many cells are expressing Tac4
# Here, let's consider a cell to be expressing Tac4 if the expression level is greater than 0
expressing_cells = (tac2_expression > 0).sum()

print(f'Number of cells in CRABP1 cluster expressing Tac2: {expressing_cells}')


#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata, groupby = 'class', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df

# Export to CSV 
csv_filename = "small_degs_df.csv"
small_degs_df.to_csv(csv_filename, index=False)


sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY2'], title = "Mouse_Inhibitory_Neuorns_subcluster" + str(resolution), frameon = False, legend_loc = "on data")


# Define a dictionary to map old cluster names to new names
cluster_name_mapping = {
    '4' : 'Progenitor',
    '12' : 'Progenitor',
    '3' : 'MGE_LHX6/MAF',
    '0' : 'MGE_LHX6/MAF',
    '8' : 'MGE_LHX6/MAF',
    '1' : 'MGE_LHX6/MAF',
    '2' : 'MGE_LHX6/MAF',
    '9' : 'MGE_LHX6/MAF',
    '6' : 'MGE_LHX6/NPY',
    '5,1' : 'MGE_CRABP1/CHRNA3',
    '5,0' : 'MGE_CRABP1/MAF',
    '10,1' : 'CGE_NR2F2/PROX1',
    '10,0' : 'LGE_MEIS2/PAX6',
    '10,2' : 'LGE_MEIS2/PAX6',
    '11' : 'LGE_FOXP1/ISL1',
    '7' : 'LGE_FOXP1/PENK',
    
    
 
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['leiden_subcluster_CRABPY2'] = adata.obs['leiden_subcluster_CRABPY2'].replace(cluster_name_mapping)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY2'], title = "Mouse_Dev", frameon = False,  size =10,save = 'speduotimeclusters.svg')

# Copy the 'leiden_subcluster_CRABPY2' observation to a new 'class' observation
adata.obs['class'] = adata.obs['leiden_subcluster_CRABPY2']


adata.obs['class'] = adata.obs['class'].replace('MGE_CRABP1/CHRNA3', 'MGE_CRABP1/TH')


paldict={'CGE_NR2F2/PROX1': 'mediumpurple',
    'LGE_FOXP1/ISL1': 'royalblue',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_MEIS2/PAX6': 'orangered',
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TH': 'fuchsia',
    'MGE_LHX6/MAF': 'skyblue',
    'MGE_LHX6/NPY': 'teal',
    'Progenitor': 'pink',


}



plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color=  ['class'], palette=paldict, title = "Mouse_inhibitory", size =10,frameon = False, save = '_officialclustersupdated.svg')




adata.obs['class'] = adata.obs['class'].cat.reorder_categories(['Progenitor', 'MGE_LHX6/MAF', 'MGE_LHX6/NPY', 'MGE_CRABP1/MAF',
                                                                                                      'MGE_CRABP1/TH',
                                                                                                       'LGE_MEIS2/PAX6','LGE_FOXP1/PENK', 
                                                                                                        'LGE_FOXP1/ISL1', 'CGE_NR2F2/PROX1',
                                                                                                         ])



#Matrix of marker genes for inhibitory neuron clusters
sc.set_figure_params(fontsize = 13)
#marker genes based on Schmitz et al. 2022
marker_genes_dict = {
    'Progenitor' : ['Aspm', 'Cenph', 'Mki67'],
    'MGE_LHX6/MAF': ['Nkx2-1','Lhx6','Sox6','Cux2', 'Rbfox1', 'Mef2c'],
    'MGE_LHX6/NPY': ['Sst','Npy'],
    'MGE_CRABP1/MAF': ['Mafb', 'Rbp4', 'Col19a1', 'Chl1','Galntl6', 'Olfm2', 'Crabp1', 'Etv1','Kit'],
    'MGE_CRABP1/TH': ['Pde1c','Trh','Trhde', 'Lhx8', 'Chrna3', 'Stxbp6', 'Th', 'Tac2'],
    'LGE_MEIS2/PAX6': ['Pax6', 'Meis2'],
    'LGE_FOXP1/PENK': ['Foxp1', 'Six3', 'Penk'],
    'LGE_FOXP1/ISL1': ['Rxrg', 'Rarb', 'Isl1', 'Npy1r'],
    'CGE_NR2F2/PROX1' : ['Nr2f2', 'Prox1', 'Sp8', 'Pdzrn3', 'Nr3c2',]
    
}


sc.pl.matrixplot(adata, groupby='class', var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes=False, save='matrixVERTICAL.svg')


adata.write('/pathto/mouse/mouse_processed.h5ad')

#Naming for pseudotime analysis
# Define a dictionary to map old cluster names to new names
cluster_name_mapping = {
    '4' : 'Progenitor',
    '12' : 'Progenitor',
    '3' : 'Transition1',
    '0' : 'Transition1',
    '8' : 'Transition1',
    '1' : 'MGE_LHX6/MAF',
    '2' : 'MGE_LHX6/MAF',
    '9' : 'MGE_LHX6/MAF',
    '6' : 'MGE_LHX6/NPY',
    '5,1' : 'MGE_CRABP1/TH',
    '5,0' : 'MGE_CRABP1/MAF',
    '10,1' : 'CGE_NR2F2/PROX1',
    '10,0' : 'Transition2',
    '10,2' : 'LGE_MEIS2/PAX6',
    '11' : 'LGE_FOXP1/ISL1',
    '7' : 'LGE_FOXP1/PENK',
    
    
    
 
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['leiden_subcluster_CRABPY2'] = adata.obs['leiden_subcluster_CRABPY2'].replace(cluster_name_mapping)




paldict2={'CGE_NR2F2/PROX1': 'mediumpurple',
    'LGE_FOXP1/ISL1': 'royalblue',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_MEIS2/PAX6': 'orangered',
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TH': 'fuchsia',
    'MGE_LHX6/MAF': 'skyblue',
    'MGE_LHX6/NPY': 'teal',
    'Progenitor': 'pink',
    'Transition1' : 'gray',
    'Transition2' : 'lightgray'}


sc.pl.umap(adata, color= ['leiden_subcluster_CRABPY2'], palette=paldict2,frameon = False, legend_loc = "on data")



Clusters = adata[adata.obs['class'].isin(['MGE_CRABP1/MAF', 'MGE_CRABP1/CHRNA3']) == True]





plt.rcParams["figure.figsize"] = (10, 10)
#marker genes based on Schmitz et al. 2022
marker_genes_dict = {

    'MGE_CRABP1/MAF': ['Crabp1','Etv1', 'Kit','Maf', 'Mafb','Col19a1','Rbp4'],
    'MGE_CRABP1/CHRNA3': ['Trhde', 'Tac2', 'Th','Chrna3','Chrna7','Stxbp6', 'Lhx8', 'Zic1', 'Zic2', 'Zic4']

}
#'DCC', 'GAP43', 'PPP2R2B', 'RALYL', 'CACNA2D3',
# 'SYT1', 'MEF2C', 'BCL11B', 'SOX1', 'ARID1B'
# 'NFIA', 'CNTN5', 'LINGO2', 'ZIC1', 'ZIC2', 'ZIC4', 'KITLG'


sc.pl.matrixplot(Clusters, groupby=['class'] ,var_names=marker_genes_dict, vmin=-2, vmax=2, cmap='PiYG_r',   figsize=(9, 1), swap_axes = False, save='matrixtac3.svg')



# Rename specific columns in obs
adata.obs.rename(columns={'leiden_subcluster_CRABPY': 'leiden_subcluster','leiden_subcluster_CRABPY2': 'leiden_subcluster2' }, inplace=True)

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
obs.to_csv('mouse_barcodes.tsv', sep='\t', index=True)

# Save var.tsv
var.to_csv('mouse_genes.tsv', sep='\t', index=True)

# Save matrix.mtx
scipy.io.mmwrite('mouse_matrix.mtx', matrix)



adata.write('/pathto/mouse/mouse_devinhibneurons.h5ad')


# # subset data for integration



adata = sc.read_h5ad('/pathto/mouse/mouse_processed.h5ad')


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
adata.write('/pathto/mouse/mouse_raw1k.h5ad')

