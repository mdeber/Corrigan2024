#!/usr/bin/env python
# coding: utf-8

# ## Integrated cross species analysis

# In[1]:


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


# In[2]:


sc.settings.verbosity = 0     


# In[3]:


import matplotlib as mpl
mpl.rcParams['figure.facecolor'] = 'white'


# In[4]:


#Uploading raw, downsampled data (maximum of 1k cells per initial class for each species)
Mouse = sc.read_h5ad('/pathto/mouse/mouse_raw1k.h5ad')
Rat = sc.read_h5ad('/pathto/rat/rat_subsetraw1k.h5ad')
Macaque = sc.read_h5ad('/pathto/macaque_dev/mac_subsetraw1k.h5ad')
Rabbit = sc.read_h5ad('/pathto/rabbit/rabbit_raw1k.h5ad')
Pig = sc.read_h5ad('/pathto/pig/updated_genome/pig_subsetraw1k.h5ad')
Ferret = sc.read_h5ad('/pathto/ferret/updated_genome/revisions/ferret_raw1k.h5ad')
Opossum = sc.read_h5ad('/pathto/opossum/opossum_subsetraw1k.h5ad')
Sugarglider = sc.read_h5ad('/pathto/sugarglider/extended/sugarglider_subsetraw1k.h5ad')


# In[5]:


for gem in [Mouse, Macaque, Rabbit, Pig, Ferret, Opossum, Sugarglider, Rat]:
    print(gem.shape[0])


# In[6]:


#set a path to your working directory
Results_file=('/pathto/concat.h5ad')


# In[7]:


#convert mouse to human gene names
orthos=pd.read_csv('/pathto/mouse/HOM_AllOrganism.rpt',sep='\t')
#cut out other species in the reference
orthos=orthos.loc[orthos['NCBI Taxon ID'].isin([10090,9606]),:]
classcounts=orthos['DB Class Key'].value_counts()
one2one=classcounts.index[list(classcounts==2)]
orthos=orthos.loc[orthos['DB Class Key'].isin(one2one),:]
#tell it what is human and what is mouse
htab=orthos.loc[orthos['NCBI Taxon ID']==9606,:]
mtab=orthos.loc[orthos['NCBI Taxon ID']==10090,:]
genemapping = dict(zip(mtab['Symbol'], htab['Symbol']))
Mouse=Mouse[:,Mouse.var.index.isin(genemapping.keys())]
Mouse.var.index=[genemapping[x] for x in Mouse.var.index]


# In[8]:


Mouse.var_names_make_unique()
Mouse.obs_names_make_unique()


# In[9]:


Mouse.var


# In[10]:


# Load homolog file again for rat
orthos2 = pd.read_csv('/Users/emilycorrigan/Downloads/HOM_AllOrganism.rpt.txt', sep='\t')

# Keep only rat (10116) and human (9606) genes
orthos2 = orthos2.loc[orthos2['NCBI Taxon ID'].isin([10116, 9606]), :]

# Filter for one-to-one mappings
classcounts = orthos2['DB Class Key'].value_counts()
one2one = classcounts.index[classcounts == 2]
orthos2 = orthos2.loc[orthos2['DB Class Key'].isin(one2one), :]

# Separate human and rat genes
htab = orthos2.loc[orthos2['NCBI Taxon ID'] == 9606, ['DB Class Key', 'Symbol']]
rtab = orthos2.loc[orthos2['NCBI Taxon ID'] == 10116, ['DB Class Key', 'Symbol']]

# Merge on DB Class Key to ensure correct alignment
merged = rtab.merge(htab, on='DB Class Key', suffixes=('_rat', '_human'))

# Create mapping dictionary where rat gene names are mapped to human gene names (capitalized)
genemapping2 = {r: h.upper() for r, h in zip(merged['Symbol_rat'], merged['Symbol_human'])}

# Function to map rat gene names to human gene names (capitalized)
def map_gene_name(gene):
    return genemapping2.get(gene, gene)  # if gene is not in mapping, return original name

# Apply mapping to Rat.var.index, ensuring rat genes are mapped to human names (capitalized)
mapped_genes = Rat.var.index.map(map_gene_name)

# Update Rat.var.index with mapped gene names
Rat.var.index = mapped_genes


# In[11]:


Rat.var


# In[12]:


Rat.var_names_make_unique()
Rat.obs_names_make_unique()


# In[13]:


Rat.X


# In[14]:


Rat.var


# In[15]:


#concat data
adata = anndata.AnnData.concatenate(Mouse, Rat, Rabbit, Macaque, Pig, Ferret, Sugarglider, Opossum, join='outer', batch_key='species', batch_categories=['Mouse', 'Rat', 'Rabbit', 'Macaque', 'Pig', 'Ferret', 'Sugarglider', 'Opossum'])


# In[16]:


adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype(str)


# In[17]:


adata.obs['timepoint'] = adata.obs['timepoint'].astype(str)


# In[18]:


#Convert so it can save
def convert_highly_variable_to_string(adata):
    """Converts all columns containing 'highly_variable' to strings in `adata.var`."""
    for col in adata.var.columns:
        if any(key in col for key in ["highly_variable", "highly_variable_nbatches", "highly_variable_intersection"]):
            adata.var[col] = adata.var[col].astype(str)
    return adata

# Apply before saving
adata = convert_highly_variable_to_string(adata)
adata.write(Results_file)


# In[19]:


adata.write(Results_file)


# In[20]:


from scipy.sparse import issparse
# Get indices for each species
species_idx = adata.obs["species"].values  # This assumes species is stored in `adata.obs`

# Get CRABP1 column index
crabp1_idx = np.where(adata.var_names == "CRABP1")[0][0]

# Convert to dense temporarily for checking
adata_dense = adata.X.toarray() if issparse(adata.X) else adata.X

# Check which species have NaNs in CRABP1
crabp1_counts = pd.DataFrame({"species": species_idx, "CRABP1_counts": adata_dense[:, crabp1_idx]})
nan_species = crabp1_counts[crabp1_counts["CRABP1_counts"].isna()]["species"].unique()
print("Species with NaNs in CRABP1:", nan_species)


# In[21]:


# Mask for sugar glider 
sugarglider_mask = adata.obs["species"] == "Sugarglider"

# Replace only NaNs 
adata.X[sugarglider_mask, crabp1_idx] = np.nan_to_num(adata.X[sugarglider_mask, crabp1_idx])

# Confirm the fix
print(f"NaNs in CRABP1 after fix: {np.isnan(adata.X[:, crabp1_idx]).sum()}")
print(f"CRABP1 total counts after fix: {adata[:, crabp1_idx].X.sum()}")


# In[22]:


# Find common genes across all species
common_genes = set(Mouse.var_names)  & set(Ferret.var_names) &set(Rat.var_names) & set(Macaque.var_names) & set(Rabbit.var_names) & set(Pig.var_names) & set(Sugarglider.var_names) &set(Opossum.var_names)

# Ensure 'TAC3' and 'CRABP1' are included
final_genes = common_genes | {'TAC3', 'CRABP1'}

# Filter directly using Pandas index intersection (faster than list comprehension)
filtered_gene_mask = adata.var_names.isin(final_genes)
filtered_adata = adata[:, filtered_gene_mask].copy()


# In[23]:


filtered_adata.var


# In[24]:


print("TAC3 present:", "TAC3" in filtered_adata.var_names)
print("CRABP1 present:", "CRABP1" in filtered_adata.var_names)


# In[25]:


adata = filtered_adata


# In[27]:


adata.obs_names_make_unique()
adata.var_names_make_unique()


# In[28]:


#Set the .raw attribute of the AnnData object
adata.raw = adata


# In[29]:


from scipy.sparse import issparse
if issparse(adata.X):
    adata.X = adata.X.toarray()  # Convert sparse to dense


# In[30]:


adata.var = adata.var.astype(str)
adata.obs = adata.obs.astype(str)


# In[32]:


print(adata.X)


# In[33]:


#Define string to append to saved figures
save_name = "concat"


# ## Normalize counts scanpy tutorial

# In[34]:


sc.pp.filter_genes(adata, min_counts=1)


# In[35]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)


# In[36]:


#Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)


# In[37]:


adata.write(Results_file)


# ## PCA

# In[38]:


adata = sc.read_h5ad(Results_file)


# In[39]:


#Perform PCA dimensional reduction for the AnnData object.
sc.tl.pca(adata, svd_solver='arpack')


# In[40]:


#Inspect the contribution of each PC to the variance of the data.
plt.rcParams["figure.figsize"] = (5, 5)
sc.pl.pca_variance_ratio(adata, log=True)


# In[41]:


#Visualize the PCA loadings.
plt.rcParams["figure.figsize"] = (10, 5) 
sc.pl.pca_loadings(adata, components = [1, 2, 3, 4])


# In[42]:


#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ["VIM", "MEF2C", "DLX5",
                 "GAD2", 'species', 'batch']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


# # Batch correction with Harmony

# In[46]:


sc.external.pp.harmony_integrate(adata, ['batch', 'species'])


# In[47]:


adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']


# In[48]:


#Visualize genes contributing to most variance.
plt.rcParams["figure.figsize"] = (4, 4) 
pca_var_genes = ['batch', 'species']
sc.pl.pca(adata, color = pca_var_genes, frameon = True)


# In[49]:


adata.write(Results_file)


# ## Calculate UMAP

# In[51]:


adata = sc.read_h5ad(Results_file)


# In[50]:


#Compute the nearest neighbors using the knn algorithm.
sc.pp.neighbors(adata, knn = True, method = "umap", metric = "euclidean")


# In[51]:


#Compute the umap embedding based on knn-computed neightbors.
sc.tl.umap(adata)


# In[52]:


#Save the AnnData object as an .h5ad file.
adata.write(Results_file)


# # Plot UMAP

# In[5]:


adata.obs


# In[10]:


print(adata.raw.X)


# In[9]:


#save raw data
import pandas as pd
import numpy as np
from scipy import sparse

X = adata.raw.X

if sparse.issparse(X):
    # If sparse, convert to dense or use sparse-aware DataFrame
    df = pd.DataFrame.sparse.from_spmatrix(
        X,
        index=adata.obs_names,
        columns=adata.raw.var_names
    )
else:
    # If already dense, use as-is
    df = pd.DataFrame(
        X,
        index=adata.obs_names,
        columns=adata.raw.var_names
    )

# Save to file
df.to_csv("adata_raw_expression.csv", sep='\t')


# In[4]:


adata = sc.read_h5ad(Results_file)


# In[6]:


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


# In[33]:


import PyComplexHeatmap as pch
from PyComplexHeatmap import *
from PyComplexHeatmap import clustermap


# In[7]:


# Define a dictionary to map old macaque cluster names to new names
cluster_name_mapping = {
    'G1-phase_SLC1A3/ATP1A1' : 'Progenitor',
    'G2-M_UBE2C/ASPM' : 'Progenitor',
    'S-phase_MCM4/H43C' : 'Progenitor',
    'Transition' : 'Progenitor',
    'MGE_CRABP1/CHRNA3' : 'MGE_CRABP1/TH'
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['class'] = adata.obs['class'].replace(cluster_name_mapping)


# In[26]:


# Cluster umap embeddings using leiden and save umap plots.
resolution = 1 #Adjust resolution
sc.tl.leiden(adata, resolution = resolution)

plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= 'leiden', title = "Concat" + str(resolution), frameon = False, legend_loc = "on data")


# In[25]:


# Choose the colorblind palette
palette = sns.color_palette("colorblind", n_colors=10)

# Manually assign colors to labels
manual_color_assignment = {
    "Macaque": palette[1],     
    "Mouse": palette[0],     
    "Rat": palette[9], 
    "Rabbit": palette[3],         
    "Ferret": palette[2],
    "Pig": palette[4],
    "Opossum": palette[7], 
    "Sugarglider": palette[8],
    "not": "grey"
}
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['species'], palette = manual_color_assignment, frameon = False, wspace=0.25,ncols = 1,save = 'speciesumap.svg', size = 10)
save_nature_figure('naturespeciesumap.svg')


# In[7]:


adata.obs['class']


# In[8]:


#renaming during revisions
adata.obs['class'] = adata.obs['class'].replace({'MGE_CRABP1/CHRNA3': 'MGE_CRABP1/TH'})


# In[7]:


paldict={'CGE_NR2F2/PROX1': 'mediumpurple',
    'LGE_FOXP1/ISL1': 'royalblue',
    'LGE_FOXP1/ISL1/NPY1R': 'blue',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_MEIS2/PAX6/SCGN': 'orange',
    'LGE_FOXP2/TSHZ1': '#17344c',
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TAC3': 'fuchsia',
    'MGE_CRABP1/TH': 'hotpink',
    'MGE_LHX6/MAF': 'skyblue',
    'MGE_LHX6/NPY': 'teal',
    'Progenitor': 'pink',
    'LGE_MEIS2/PAX6': 'orangered',
    'RMTW_ZIC1/RELN': 'lightcoral',
    'Transition': '#3c7632',
    'VMF_ZIC1/ZIC2': 'green',
    'VMF_CRABP1/LHX8': 'mediumspringgreen',
}


# In[12]:


# Check for clusters due to technical issues.
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['class'], palette = paldict,wspace=0.25,ncols = 1, frameon= False, size = 10, save = 'original_classumap.svg')
#save_nature_figure('natureclassumap.svg')
#sc.pl.violin(adata, ['species'], groupby = "leiden_subcluster")


# In[22]:


#adding in missing info
adata.obs.loc[adata.obs['species'] == 'Mouse', 'region'] = 'Striatum'
# Add 'P20' to categories if it's not there
if 'P20' not in adata.obs['timepoint'].cat.categories:
    adata.obs['timepoint'] = adata.obs['timepoint'].cat.add_categories('P20')

# Now assign 'P20' to Opossum rows
adata.obs.loc[adata.obs['species'] == 'Opossum', 'timepoint'] = 'P20'


# In[24]:


pd.set_option('display.max_columns', 10)


# In[13]:


adata.obs['species'] = adata.obs['species'].cat.reorder_categories(['Macaque', 'Rabbit', 'Rat', 'Mouse', 'Pig', 'Ferret', 'Sugarglider', 'Opossum'])

plot_keys = []
for specie in adata.obs['species'].unique():
    is_flag = specie
    adata.obs[is_flag] = "not"
    adata.obs[is_flag][adata.obs['species'] == specie] = specie
    plot_keys.append(is_flag)


# In[14]:


plot_keys = ['Macaque', 'Rabbit', 'Rat', 'Mouse', 'Pig', 'Ferret',  'Opossum','Sugarglider',]


# In[15]:


sc.set_figure_params(fontsize = 65)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = plot_keys, groups = plot_keys, palette = manual_color_assignment,
           wspace=0,ncols = 8, frameon= False, size = 25, legend_loc = None, save = 'specieswarhol.svg')


# In[16]:


adata.obs["original_classxspecies"] = adata.obs["species"].astype(str) + "_" + adata.obs["class"].astype(str)


# In[17]:


adata.obs["original_classxspecies"] = adata.obs["original_classxspecies"].astype("category")

sc.tl.dendrogram(adata, "original_classxspecies")


# In[18]:


sc.tl.dendrogram(adata, "original_classxspecies")
sc.set_figure_params(fontsize = 10)
sc.pl.correlation_matrix(adata, "original_classxspecies",figsize=(30, 30), cmap= 'RdYlBu_r', save = 'correlationoriginal_classxspecies.svg' )


# In[9]:


#sc.set_figure_params(dpi_save= 500)


# In[68]:


plt.rcParams['axes.grid'] = False


# In[20]:


adata.obs['class']


# In[21]:


paldict1={
    'Progenitor': 'pink', 
    'MGE_LHX6/MAF': 'skyblue',  
    'MGE_LHX6/NPY': 'teal',
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TAC3': 'fuchsia',
    'MGE_CRABP1/TH': 'hotpink',
    'LGE_MEIS2/PAX6': 'orangered',

    'LGE_FOXP2/TSHZ1': '#17344c',
     'LGE_FOXP1/PENK': 'navy',
    'LGE_FOXP1/ISL1': 'royalblue',  
    'LGE_FOXP1/ISL1/NPY1R': 'blue',
    'CGE_NR2F2/PROX1': 'mediumpurple',
    'VMF_CRABP1/LHX8': 'mediumspringgreen',
    'VMF_ZIC1/ZIC2': 'green',
    'RMTW_ZIC1/RELN': 'lightcoral',
}


# In[22]:


# Step 1: Get PCA representation
X_pca_class = adata.obsm['X_pca']  # same PCA

# Step 2: Get group labels based on 'class' instead of 'integrated_class'
groups_class = adata.obs['original_classxspecies']  # assuming you have this column like 'class_species'

# Step 3: Group PCA matrix by label and compute mean per group
grouped_class = pd.DataFrame(X_pca_class).groupby(groups_class.values).mean()

# Step 4: Pearson correlation between group means
corr_matrix_class = grouped_class.T.corr(method='pearson')


# In[23]:


corr_matrix_class


# In[24]:


# Step 5: Split label into species + class
split_df_class = corr_matrix_class.index.to_series().str.split("_", n=1, expand=True)
split_df_class.columns = ['species', 'class']
split_df_class.index = corr_matrix_class.index


# In[25]:


split_df_class.index


# In[28]:


# Step 5: Split label into species + class
split_df_class = corr_matrix_class.index.to_series().str.split("_", n=1, expand=True)
split_df_class.columns = ['species', 'class']
split_df_class.index = corr_matrix_class.index

# Step 6: Create annotation dataframe
annot_df_class = pd.DataFrame({
    'Species': split_df_class['species'],
    'Class': split_df_class['class']
})

# Step 7: Define annotation colors (reuse your color dicts)
annotation_colors_class = {
    'Species': manual_color_assignment,
    'Class': paldict1  # or paldict2, whichever you're using
}

# Step 8: Annotations
plt.rcParams.update({'font.size': 16})
row_ha_class = HeatmapAnnotation(
    Class=anno_simple(df=annot_df_class['Class'], colors=paldict1, height=5),
    Species=anno_simple(df=annot_df_class['Species'], colors=manual_color_assignment),
    axis=0
)
col_ha_class = HeatmapAnnotation(
    Class=anno_simple(df=annot_df_class['Class'], colors=paldict1, height=5),
    Species=anno_simple(df=annot_df_class['Species'], colors=manual_color_assignment),
    axis=1
)

# Step 9: Plot
plt.figure(figsize=(20, 16))
cm_class = ClusterMapPlotter(
    data=corr_matrix_class,
    left_annotation=row_ha_class,
    top_annotation=col_ha_class,
    cmap='RdYlBu_r',
    vmin=-1,
    vmax=1,
    center=0,
    row_dendrogram=True,
    col_dendrogram=True,
)

# Step 10: Save figure
plt.grid(False)
plt.rcParams['savefig.dpi'] = 300
plt.savefig('correlation_heatmap_originalclassxspecies.svg', bbox_inches='tight')
plt.show()


# In[ ]:





# In[67]:


adata.obs


# In[209]:


#subcluster 
sc.tl.leiden(adata, restrict_to=('leiden', ['1']), resolution=0.4, key_added='leiden_subcluster')
sc.pl.umap(adata, color= ['leiden_subcluster'], title = "subcluster" + str(resolution), frameon = False, legend_loc = "on data")


# In[210]:


#subcluster
plt.rcParams["figure.figsize"] = (10, 10)
sc.tl.leiden(adata, restrict_to=('leiden_subcluster', ['0']), resolution=0.5, key_added='leiden_subcluster2')
sc.pl.umap(adata, color= ['leiden_subcluster2'], title = "subcluster" + str(resolution), frameon = False, legend_loc = "on data")


# In[217]:


#subcluster 
sc.tl.leiden(adata, restrict_to=('leiden_subcluster2', ['6']), resolution=0.3, key_added='leiden_subcluster3')
sc.pl.umap(adata, color= ['leiden_subcluster3'],title = "subcluster" + str(resolution), frameon = False, legend_loc = "on data")


# In[218]:


#subcluster TAC3 cluster so that it is only TAC3 cells
sc.tl.leiden(adata, restrict_to=('leiden_subcluster3', ['7']), resolution=0.4, key_added='leiden_subcluster4')
sc.pl.umap(adata, color= ['leiden_subcluster4'], title = "subcluster" , frameon = False, legend_loc = "on data")


# In[7]:


# Define a dictionary to map old cluster names to new names
cluster_name_mapping = {
 
    '1,3' : 'MGE_CRABP1/MAF',
    '1,0' : 'MGE_CRABP1/MAF',
    '1,4' : 'MGE_CRABP1/MAF',
    '14' : 'MGE_CRABP1/MAF',

    '1,5' : 'MGE_CRABP1/MAF',
    '1,1' : 'MGE_CRABP1/MAF',
    
    
    '1,2' : 'MGE_CRABP1/TAC3',
    
    '10' : 'Progenitor',
    '8' : 'Progenitor',
    '9' : 'Progenitor',


    
    '0,2' : 'MGE_LHX6/MAF',
    '0,5' : 'MGE_LHX6/MAF',
    '0,3' : 'MGE_LHX6/MAF',
    '0,0' : 'MGE_LHX6/MAF',
    '0,4' : 'MGE_LHX6/MAF',
    '0,6' : 'MGE_LHX6/MAF',
    '16' : 'MGE_LHX6/MAF',

    
    '0,1' : 'MGE_LHX6/NPY',

        
    '4' : 'CGE_NR2F2/PROX1',
    
    
    '2' : 'LGE_MEIS2/PAX6',



    '3' : 'LGE_FOXP1/PENK',
    '7,1':'LGE_FOXP1/PENK',

    
    
    '7,2' : 'LGE_FOXP1/ISL1',
    '6,2' : 'LGE_FOXP1/ISL1',
    '6,0' : 'LGE_FOXP1/ISL1',
        
    '6,1' : 'LGE_FOXP1/ISL1/NPY1R',
    '7,0' : 'LGE_FOXP1/ISL1/NPY1R',
    


    '5' : 'LGE_FOXP2/TSHZ1',
    '17' : 'LGE_FOXP2/TSHZ1',
    '12' : 'LGE_FOXP2/TSHZ1',

    '15' :'VMF_CRABP1/LHX8',
    '13' :'VMF_CRABP1/LHX8',
    
     '11' :'RMTW_ZIC1/RELN',

 
    # Add more mappings as needed
}

# Use the replace function to rename the clusters
adata.obs['leiden_subcluster4'] = adata.obs['leiden_subcluster4'].replace(cluster_name_mapping)


# In[6]:


sc.pl.umap(adata, color= ['leiden_subcluster4'], title = "subcluster", frameon = False, legend_loc = "on data")


# In[38]:


paldict={
    'LGE_FOXP2/TSHZ1': '#17344c',
    'LGE_FOXP1/ISL1/NPY1R': 'blue',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_FOXP1/ISL1': 'royalblue',  
    'MGE_LHX6/NPY': 'teal',
    'MGE_LHX6/MAF': 'skyblue',    
    'RMTW_ZIC1/RELN': 'lightcoral',
    'VMF_CRABP1/LHX8': 'mediumspringgreen',
   
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TAC3': 'fuchsia',

    'Progenitor': 'pink', 
    
    
    
    'CGE_NR2F2/PROX1': 'mediumpurple',
  #  'G1-phase_SLC1A3/ATP1A1': '#17344c',
   # 'G2-M_UBE2C/ASPM': '#19122b',




    'LGE_MEIS2/PAX6': 'orangered',

   # 'Speckle' : 'gray',
   # 'nr2f2/prox1?': 'mediumpurple',
#    'RMTW_ZIC1/RELN': 'yellow',
 #   'S-phase_MCM4/H43C': 'lawngreen',
    'VMF_ZIC1/ZIC2': 'green',

#    'VMF_LHX1/POU6F2':'',
  #  'VMF_PEG10/DLK1':'steelblue',
 #   'LGE_FOXP1/ISL1/NPY1R':'mediumpurple',
   # 'nan':'black'
}


# In[41]:


sc.pl.umap(adata, color= ['leiden_subcluster4'], title = "subcluster", palette = paldict,frameon = False, legend_loc = "on data")


# In[42]:


adata.obs['integrated_class'] = adata.obs['leiden_subcluster4'].copy()
adata.obs['integrated_class'] = adata.obs['integrated_class'].astype('category')


# In[43]:


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color= ['integrated_class'], title = "subcluster",  palette=paldict, frameon = False, size =10, save = 'integrated_clusters.svg')


# In[44]:


# Cross-tabulate counts of species per integrated_class
ct = pd.crosstab(adata.obs['integrated_class'], adata.obs['species'])

# Normalize to proportions (percent per integrated_class)
ct_prop = ct.div(ct.sum(axis=1), axis=0)
ct_prop_reversed = ct_prop[species].iloc[::-1]

# Ensure species order matches color assignment
species = [s for s in manual_color_assignment.keys() if s in ct_prop.columns]
colors = [manual_color_assignment[s] for s in species]

# Plot
ax = ct_prop_reversed.plot(
    kind='barh',
    stacked=True,
    color=colors,
    figsize=(8, 8)
)

plt.ylabel("Proportion of cells")
plt.title("Species in integrated_class")
plt.legend(title="Species", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('bargraphspecies.svg')
plt.show()


# In[ ]:





# In[45]:


sc.set_figure_params(fontsize = 10)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['integrated_class'], groups = ['MGE_CRABP1/TAC3','MGE_CRABP1/MAF'], frameon=False,size= 10,wspace=0.25, ncols = 2, save = 'integratedclassCRABP1cluster.svg')


# In[46]:


#sc.set_figure_params(fontsize = 10)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(adata, color = ['species'], frameon=False,size= 10,wspace=0.25, ncols = 2, save = 'integratedclassCRABP1clusterspecies.svg')


# In[48]:


adata.obs["integrated_classxspecies"] = adata.obs["species"].astype(str) + "_" + adata.obs["integrated_class"].astype(str)


# In[49]:


adata.obs["integrated_classxspecies"] = adata.obs["integrated_classxspecies"].astype("category")


# In[50]:


adata.obs["integrated_classxspecies"]


# In[70]:


sc.tl.dendrogram(adata, "integrated_classxspecies")
ax = sc.pl.dendrogram(adata, "integrated_classxspecies", show=False)
fig = ax.figure
fig.set_size_inches(20, 5)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=10)
for line in ax.get_lines():
    line.set_linewidth(0.5)  # Adjust the line width here
fig.tight_layout()
plt.savefig('integratedclass_dendrogram.svg')
plt.show()


# In[53]:


# Step 1: Get PCA representation
X_pca = adata.obsm['X_pca']

# Step 2: Get group labels
groups = adata.obs['integrated_classxspecies']

# Step 3: Group the PCA matrix by your label and compute the mean PCA per group
grouped = pd.DataFrame(X_pca).groupby(groups.values).mean()

# Now `grouped` is a (n_groups x n_pcs) matrix

# Step 4: Compute Pearson correlation matrix between group means
corr_matrix = grouped.T.corr(method='pearson')


# In[69]:


# Step 1: Get PCA representation
X_pca_class = adata.obsm['X_pca']  # same PCA

# Step 2: Get group labels based on 'class' instead of 'integrated_class'
groups_class = adata.obs['integrated_classxspecies']  # assuming you have this column like 'class_species'

# Step 3: Group PCA matrix by label and compute mean per group
grouped_class = pd.DataFrame(X_pca_class).groupby(groups_class.values).mean()

# Step 4: Pearson correlation between group means
corr_matrix_class = grouped_class.T.corr(method='pearson')

# Step 5: Split label into species + class
split_df_class = corr_matrix_class.index.to_series().str.split("_", n=1, expand=True)
split_df_class.columns = ['species', 'class']
split_df_class.index = corr_matrix_class.index

# Step 6: Create annotation dataframe
annot_df_class = pd.DataFrame({
    'Species': split_df_class['species'],
    'Class': split_df_class['class']
})

# Step 7: Define annotation colors (reuse your color dicts)
annotation_colors_class = {
    'Species': manual_color_assignment,
    'Class': paldict1  # or paldict2, whichever you're using
}

# Step 8: Annotations
plt.rcParams.update({'font.size': 16})
row_ha_class = HeatmapAnnotation(
    Class=anno_simple(df=annot_df_class['Class'], colors=paldict1, height=5),
    Species=anno_simple(df=annot_df_class['Species'], colors=manual_color_assignment),
    axis=0
)
col_ha_class = HeatmapAnnotation(
    Class=anno_simple(df=annot_df_class['Class'], colors=paldict1, height=5),
    Species=anno_simple(df=annot_df_class['Species'], colors=manual_color_assignment),
    axis=1
)

# Step 9: Plot
plt.figure(figsize=(20, 16))
cm_class = ClusterMapPlotter(
    data=corr_matrix_class,
    left_annotation=row_ha_class,
    top_annotation=col_ha_class,
    cmap='RdYlBu_r',
    vmin=-1,
    vmax=1,
    center=0,
    row_dendrogram=True,
    col_dendrogram=True,
)

# Step 10: Save figure
plt.grid(False)
plt.savefig('correlation_heatmap_originalclassxspecies.svg', bbox_inches='tight')
plt.show()


# In[56]:


split_df = corr_matrix.index.to_series().str.split("_", n=1, expand=True)
split_df.columns = ['species', 'integrated_class']
split_df.index = corr_matrix.index


# In[57]:


# Annotation dataframe
annotations = pd.DataFrame({
    'Species': split_df['species'],
    'Integrated Class': split_df['integrated_class']
})

# Colors (already defined)
annotation_colors = {
    'Species': manual_color_assignment,
    'Integrated Class': paldict
}


# In[29]:


plt.rcParams['axes.grid'] = False
plt.rcParams.update({'font.size': 12})


# In[37]:


paldict2={
      'Progenitor': 'pink', 
    'MGE_LHX6/MAF': 'skyblue',  
    'MGE_LHX6/NPY': 'teal',
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TAC3': 'fuchsia',
    'LGE_MEIS2/PAX6': 'orangered',

    'LGE_FOXP2/TSHZ1': '#17344c',
     'LGE_FOXP1/PENK': 'navy',
    'LGE_FOXP1/ISL1': 'royalblue',  
    'LGE_FOXP1/ISL1/NPY1R': 'blue',
'CGE_NR2F2/PROX1': 'mediumpurple',
    'VMF_CRABP1/LHX8': 'mediumspringgreen',
    'RMTW_ZIC1/RELN': 'lightcoral',

}


# In[243]:


adata.obs


# In[22]:


# Tabulate counts
counts = pd.crosstab(adata.obs['integrated_class'], adata.obs['species'])

# Normalize by row to get proportions
pivot_counts = counts.div(counts.sum(axis=1), axis=0)


# In[31]:


# Desired species order
species_order = [
    "Macaque", "Rabbit", "Rat", "Mouse", 
    "Pig", "Ferret", "Opossum", "Sugarglider"
]

# Reorder the columns in pivot_counts to match your desired order
pivot_counts = pivot_counts[species_order]
# Check the order and make sure all species are in the manual color dictionary
species_in_plot = pivot_counts.columns.tolist()
print("Species in pivot_counts:", species_in_plot)

# This will throw an error if any species is missing from the color assignment
colors = [manual_color_assignment[species] for species in species_in_plot]

# Plot
ax = pivot_counts.plot(
    kind='barh',
    stacked=True,
    figsize=(6, 10),
    color=colors
)

# Formatting
ax.set_xlabel('Integrated Class')
ax.set_ylabel('Proportion of Cells')
ax.set_title('Proportion of Each Species in Integrated Class')
ax.legend(title='Species', bbox_to_anchor=(1.05, 1), loc='upper left')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

plt.tight_layout()
plt.show
plt.savefig('/pathto/figures/proportionboxplot.svg')


# In[27]:


# Ensure the species order in pivot_counts.columns matches your palette
colors = [manual_color_assignment[species] for species in pivot_counts.columns]

# Plot with custom colors
pivot_counts.plot(kind='bar', stacked=True, figsize=(10, 6), color=colors)

# Formatting
plt.xlabel('Integrated Class')
plt.ylabel('Proportion of Cells')
plt.title('Proportion of Species in Each Integrated Class')
plt.legend(title='Species', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.show()


# In[51]:


adata.obs['integrated_class'] = adata.obs['integrated_class'].cat.reorder_categories(['Progenitor', 'MGE_LHX6/MAF', 'MGE_LHX6/NPY', 'MGE_CRABP1/MAF','MGE_CRABP1/TAC3', 'LGE_MEIS2/PAX6',
                                                                                                       'LGE_FOXP2/TSHZ1', 'LGE_FOXP1/PENK', 
                                                                                                        'LGE_FOXP1/ISL1', 'LGE_FOXP1/ISL1/NPY1R', 'CGE_NR2F2/PROX1',
                                                                                                         'VMF_CRABP1/LHX8', 'RMTW_ZIC1/RELN'])


# In[34]:


print(adata.obs['class'])


# In[50]:


adata.obs['class'] = adata.obs['class'].cat.reorder_categories(['Progenitor', 'MGE_LHX6/MAF', 'MGE_LHX6/NPY', 'MGE_CRABP1/MAF','MGE_CRABP1/TAC3', 'MGE_CRABP1/TH','LGE_MEIS2/PAX6',
                                                                                                       'LGE_FOXP2/TSHZ1', 'LGE_FOXP1/PENK', 
                                                                                                        'LGE_FOXP1/ISL1', 'LGE_FOXP1/ISL1/NPY1R', 'CGE_NR2F2/PROX1',
                                                                                                         'VMF_CRABP1/LHX8','VMF_ZIC1/ZIC2', 'RMTW_ZIC1/RELN'])


# In[63]:


adata.obs['species'] = adata.obs['species'].cat.reorder_categories(['Macaque', 'Rabbit', 'Rat', 'Mouse', 'Pig', 'Ferret', 'Sugarglider', 'Opossum'])


# In[62]:


# Combine the two columns, ensuring the order is respected for both
adata.obs['original_classxspecies_ordered'] = (
    adata.obs['class'].astype(str) + ' | ' + adata.obs['species'].astype(str)
)

# Reorder the combined column by keeping the correct order
# This will ensure both parts (class and species) stay in the correct order
combined_order = [
    f"{ic} | {sp}" for ic in adata.obs['integrated_class'].cat.categories for sp in adata.obs['species'].cat.categories
]

# Set the combined column to be ordered
adata.obs['original_classxspecies_ordered'] = pd.Categorical(
    adata.obs['original_classxspecies_ordered'],
    categories=combined_order,
    ordered=True)


# In[64]:


adata.obs


# In[78]:


#Matrix of marker genes for inhibitory neuron clusters
sc.set_figure_params(fontsize = 17)
#marker genes based on Schmitz et al. 2022
marker_genes_dict = {
    'Progenitor' : ['ASPM', 'CENPF'],
    'MGE_LHX6/MAF': ['LHX6', 'SOX6', 'MEF2C'],
    'MGE_LHX6/NPY': ['NPY'],
    'MGE_CRABP1/MAF': ['MAF','CRABP1','ETV1','KIT','COL19A1', 'MAFB','RBP4'],
    'MGE_CRABP1/TAC3': ['TRHDE','TAC3', 'TH','CHRNA3','STXBP6', 'LHX8'],
    'LGE_MEIS2/PAX6': ['MEIS2', 'PAX6'],
    'LGE_FOXP2/TSHZ1': ['FOXP2','TSHZ1', 'EYA2'],
    'LGE_FOXP1/PENK': ['FOXP1', 'SIX3', 'PENK'],
    'LGE_FOXP1/ISL1': ['ISL1', 'RXRG', 'RARB',],
    'LGE_FOXP1/ISL1/NPY1R' : ['NPY1R'],
    'CGE_NR2F2/PROX1': ['PDZRN3', 'NR3C2', 'NR2F2', 'PROX1'],
    'VMF_CRABP1/LHX8' : ['ZIC1', 'GBX1','ZIC4'],
    'RMTW_ZIC1/RELN' : ['RELN', 'TBR1'],

}

#'DCC', 'GAP43', 'PPP2R2B', 'RALYL', 'CACNA2D3',
sc.pl.matrixplot(adata, groupby=['integrated_classxspecies_ordered'], var_names=marker_genes_dict, use_raw=False, swap_axes = True, vmin=-2, vmax=2, cmap='PiYG_r', save='matrix.svg')



# In[79]:


manual_color_assignment


# In[80]:


plt.rcParams['axes.grid'] = False
plt.rcParams.update({'font.size': 12})


# In[84]:


# Flatten gene list from marker_genes_dict
genes = [gene for sublist in marker_genes_dict.values() for gene in sublist]
# Ensure only valid genes
valid_genes = [gene for gene in genes if gene in adata.var_names]
# Extract expression matrix and obs group labels
groupby_col = 'integrated_classxspecies_ordered'
df = pd.DataFrame(adata[:, valid_genes].X.toarray(), columns=valid_genes, index=adata.obs[groupby_col])
# Compute mean expression per group
matrix_data = df.groupby(df.index).mean()
# Clip to match Scanpy plot
matrix_data = matrix_data.clip(lower=-2, upper=2)

# ----------------------------------------
# Step 2: Prepare row annotations
# ----------------------------------------
# Split the index to get species and integrated class
species_list = []
integrated_class_list = []
for idx in matrix_data.index:
    idx_str = str(idx)
    if " | " in idx_str:
        class_part, species = idx_str.split(" | ", 1)
        integrated_class = class_part
    else:
        species = "unknown"
        integrated_class = idx_str
    species_list.append(species)
    integrated_class_list.append(integrated_class)

annot_df = pd.DataFrame({
    'Species': species_list,
    'Integrated Class': integrated_class_list
}, index=matrix_data.index)

# Ensure all species have colors
if 'not' in manual_color_assignment:
    del manual_color_assignment['not']
for species in annot_df['Species'].unique():
    if species not in manual_color_assignment:
        print(f"Warning: No color assigned for '{species}', using default")
        manual_color_assignment[species] = "#999999"

row_ha = HeatmapAnnotation(
    Initial_Class=anno_simple(df=annot_df['Integrated Class'], colors=paldict2, height=10),
    Species=anno_simple(df=annot_df['Species'], colors=manual_color_assignment, height=6),
    axis=0
)

# ----------------------------------------
# Step 3: Prepare column annotations
# ----------------------------------------
if isinstance(marker_genes_dict, dict):
    gene_categories = {}
    for category, gene_list in marker_genes_dict.items():
        for gene in gene_list:
            gene_categories[gene] = category
    
    gene_cats = pd.Series({gene: gene_categories.get(gene, 'Other') 
                           for gene in matrix_data.columns})
    
    gene_colors = {cat: f"C{i}" for i, cat in enumerate(set(gene_categories.values()))}
    
    col_ha = HeatmapAnnotation(
        Category=anno_simple(df=gene_cats, colors=paldict2, height=10),
        axis=1
    )
else:
    col_ha = None

# ----------------------------------------
# Step 4: Directly modify how PyComplexHeatmap handles NaN values
# ----------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

# Import necessary tools for PyComplexHeatmap and custom colormaps
from PyComplexHeatmap import ClusterMapPlotter

# Define the gray color for NaN values
nan_color = (0.5, 0.5, 0.5, 1.0)  # 50% gray

# First approach: create a custom colormap that displays NaNs as gray
def create_custom_cmap_with_nan(base_cmap_name='PiYG_r', nan_color=nan_color):
    """Create a colormap that will properly handle NaN values"""
    base_cmap = cm.get_cmap(base_cmap_name)
    custom_cmap = base_cmap.copy()
    # This is the key step - setting the "bad" color for NaN values
    custom_cmap.set_bad(color=nan_color)
    return custom_cmap

# Create our custom colormap
custom_cmap = create_custom_cmap_with_nan(nan_color=nan_color)

# ----------------------------------------
# Step 5: Plot with PyComplexHeatmap and fix NaN rendering
# ----------------------------------------
plt.figure(figsize=(20, 50))

# Make the plot with our custom colormap
cm = ClusterMapPlotter(
    data=matrix_data,
    left_annotation=row_ha,
    top_annotation=col_ha,
    cmap=custom_cmap,
    vmin=-2,
    vmax=2,
    center=0,
    row_dendrogram=False,
    col_dendrogram=False,
    row_cluster=False,
    col_cluster=False,
    show_colnames=True,
    show_rownames=False,
    subplot_gap=0,
    xticklabels_kws=dict(
        labelrotation=90,
        labelsize=25,
        labelcolor='black',
        bottom=True
    ),
    linewidths=0,
)

# Patch the generated heatmap figure to fix NaN rendering
# This is needed because some plotting libraries don't properly handle set_bad
for ax in plt.gcf().get_axes():
    # Look for the heatmap mesh in the axes collections
    for collection in ax.collections:
        # Get the colormap from the collection
        if hasattr(collection, 'get_cmap'):
            cmap = collection.get_cmap()
            if cmap is not None:
                # Force the colormap to have gray for NaN values
                collection.get_cmap().set_bad(color=nan_color)
            
            # Get the array data and check if we have NaNs
            if hasattr(collection, 'get_array'):
                array_data = collection.get_array()
                if array_data is not None and np.any(np.isnan(array_data)):
                    # Use masked array to properly handle NaNs in the rendering
                    masked_data = np.ma.masked_invalid(array_data)
                    collection.set_array(masked_data)

# ----------------------------------------
# Step 6: Add a legend for NaN values
# ----------------------------------------
# Find the heatmap axis (usually the last or second-to-last)
main_ax = None
for ax in plt.gcf().get_axes():
    if hasattr(ax, 'collections') and len(ax.collections) > 0:
        if isinstance(ax.collections[0], plt.matplotlib.collections.QuadMesh):
            main_ax = ax
            break

if main_ax is None:
    main_ax = plt.gcf().get_axes()[-1]  # Fallback to the last axis

# Create a custom legend for NaN values
handles = []
labels = []

# Add a patch for NaN values
nan_patch = mpatches.Patch(color=nan_color, label='NaN (No expression data)')
handles.append(nan_patch)
labels.append('NaN (No expression data)')

# Add the legend to the bottom-right of the main heatmap
# Adjust positioning as needed for your specific layout
main_ax.legend(handles=handles, labels=labels, loc='right', frameon=True, 
               fontsize=14, bbox_to_anchor=(1.0, 0.0))

plt.savefig('heatmapannoLEGENDNANs.svg', bbox_inches='tight')


# In[70]:


adata


# In[20]:


# Flatten the marker_genes_dict to get the genes in order
ordered_genes = [gene for gene_list in marker_genes_dict.values() for gene in gene_list]

# Sanity check that they match the columns in matrix_data
print("Matrix columns:", matrix_data.columns.tolist())
print("Ordered genes:", ordered_genes)


# In[22]:


# Combine the two columns, ensuring the order is respected for both
adata.obs['integrated_classxspecies_ordered'] = (
    adata.obs['integrated_class'].astype(str) + ' | ' + adata.obs['species'].astype(str)
)

# Reorder the combined column by keeping the correct order
# This will ensure both parts (class and species) stay in the correct order
combined_order = [
    f"{ic} | {sp}" for ic in adata.obs['integrated_class'].cat.categories for sp in adata.obs['species'].cat.categories
]

# Set the combined column to be ordered
adata.obs['integrated_classxspecies_ordered'] = pd.Categorical(
    adata.obs['integrated_classxspecies_ordered'],
    categories=combined_order,
    ordered=True)
# Plot the matrixplot
sc.pl.matrixplot(
    adata,
    groupby='integrated_classxspecies_ordered',
    use_raw=False,
    var_names=marker_genes_dict,
    vmin=-2,
    vmax=2,
    cmap='PiYG_r',
    dendrogram=False,  # Turn off dendrogram since we're using manual ordering
    swap_axes=True,
    figsize=(20, len(marker_genes_dict) * 0.5),  # Adjust based on gene count
)


# In[14]:


#Matrix of marker genes for inhibitory neuron clusters

#marker genes based on Schmitz et al. 2022
marker_genes_dict = {
    'Progenitor' : ['ASPM', 'CENPF'],
    'MGE_LHX6/MAF': ['LHX6', 'MAF', 'SOX6', 'MEF2C'],
    'MGE_LHX6/NPY': ['NPY'],
    'MGE_CRABP1/MAF': ['CRABP1','ETV1','KIT','COL19A1', 'RBP4' ],
    'MGE_CRABP1/TAC3': ['TRHDE','TAC3', 'TH','CHRNA3','STXBP6', 'LHX8'],
    'LGE_MEIS2/PAX6': ['MEIS2', 'PAX6'],
    'LGE_FOXP2/TSHZ1': ['FOXP2','TSHZ1', 'EYA2'],
    'LGE_FOXP1/PENK': ['FOXP1', 'SIX3', 'PENK'],
    'LGE_FOXP1/ISL1': ['ISL1', 'RXRG', 'RARB'],
    'LGE_FOXP1/ISL1/NPY1R' : ['NPY1R'],
    'CGE_NR2F2/PROX1': ['PDZRN3', 'NR3C2', 'NR2F2', 'PROX1'],
    'VMF_CRABP1/LHX8' : ['ZIC1', 'GBX1','ZIC4'],
    'RMTW_ZIC1/RELN' : ['RELN', 'TBR1'],

}

#'DCC', 'GAP43', 'PPP2R2B', 'RALYL', 'CACNA2D3',
sc.pl.matrixplot(adata, groupby=['integrated_class'], var_names=marker_genes_dict, use_raw=False, swap_axes = True, vmin=-2, vmax=2, cmap='PiYG_r', save='matrixconcatenated.svg')

#put tin supplement this with the box plot, aslo the original label correlation matrix


# In[15]:


adata.obs


# In[16]:


#Jaccard with class names
# 1. Reorder class, integrated_class, and species
adata.obs['class'] = adata.obs['class'].astype('category').cat.reorder_categories([
   'Progenitor', 'MGE_LHX6/MAF', 'MGE_LHX6/NPY', 'MGE_CRABP1/MAF','MGE_CRABP1/TAC3', 'MGE_CRABP1/TH','LGE_MEIS2/PAX6',
                                                                                                       'LGE_FOXP2/TSHZ1', 'LGE_FOXP1/PENK', 
                                                                                                        'LGE_FOXP1/ISL1', 'LGE_FOXP1/ISL1/NPY1R', 'CGE_NR2F2/PROX1',
                                                                                                         'VMF_CRABP1/LHX8','VMF_ZIC1/ZIC2', 'RMTW_ZIC1/RELN']
, ordered=True)

adata.obs['integrated_class'] = adata.obs['integrated_class'].cat.reorder_categories([
    'Progenitor', 'MGE_LHX6/MAF', 'MGE_LHX6/NPY', 'MGE_CRABP1/MAF', 'MGE_CRABP1/TAC3',
    'LGE_MEIS2/PAX6', 'LGE_FOXP2/TSHZ1', 'LGE_FOXP1/PENK', 'LGE_FOXP1/ISL1',
    'LGE_FOXP1/ISL1/NPY1R', 'CGE_NR2F2/PROX1', 'VMF_CRABP1/LHX8', 'RMTW_ZIC1/RELN'
], ordered=True)

adata.obs['species'] = adata.obs['species'].cat.reorder_categories([
    'Macaque', 'Rabbit', 'Rat', 'Mouse', 'Pig', 'Ferret', 'Sugarglider', 'Opossum'
], ordered=True)

# 2. Combine class/species after reordering
adata.obs['original_classxspecies'] = (
    adata.obs['class'].astype(str) + ' | ' + adata.obs['species'].astype(str)
)
adata.obs['integrated_classxspecies'] = (
    adata.obs['integrated_class'].astype(str) + ' | ' + adata.obs['species'].astype(str)
)

# 3. Build the full ordered axis labels (species-major)
species_order = list(adata.obs['species'].cat.categories)
class_order = list(adata.obs['class'].cat.categories)
integrated_class_order = list(adata.obs['integrated_class'].cat.categories)

original_order = [f"{cl} | {sp}" for sp in species_order for cl in class_order]
integrated_order = [f"{cl} | {sp}" for sp in species_order for cl in integrated_class_order]

# 4. Build sets
original_groups = adata.obs['original_classxspecies'].unique()
integrated_groups = adata.obs['integrated_classxspecies'].unique()

original_sets = {
    label: set(adata.obs_names[adata.obs['original_classxspecies'] == label])
    for label in original_groups
}
integrated_sets = {
    label: set(adata.obs_names[adata.obs['integrated_classxspecies'] == label])
    for label in integrated_groups
}

# 5. Jaccard Index
from itertools import product

jaccard_matrix = pd.DataFrame(index=original_groups, columns=integrated_groups)

for o_label, i_label in product(original_groups, integrated_groups):
    intersection = original_sets[o_label] & integrated_sets[i_label]
    union = original_sets[o_label] | integrated_sets[i_label]
    jaccard_matrix.loc[o_label, i_label] = len(intersection) / len(union) if union else np.nan

jaccard_matrix = jaccard_matrix.astype(float)

# 6. Reorder rows and columns based on species-major order
ordered_rows = [x for x in original_order if x in jaccard_matrix.index]
ordered_cols = [x for x in integrated_order if x in jaccard_matrix.columns]
jaccard_matrix = jaccard_matrix.loc[ordered_rows, ordered_cols]

# 7. Plot
plt.figure(figsize=(40, 40))
ax = sns.heatmap(
    jaccard_matrix,
    cmap="YlOrRd",
    linewidths=0,
    linecolor=None,
    square=True,
    cbar=True,
    xticklabels=True,
    yticklabels=True
)

# Rotate ticks
plt.xticks(rotation=90)
plt.yticks(rotation=0)

# Adjust title
plt.title("Jaccard Index: Class×Species vs IntegratedClass×Species (Species-major order)", fontsize=20)

# Resize colorbar and set font size
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)  # set tick font size
cbar.ax.set_ylabel("Jaccard Index", fontsize=20)  # optional: label for colorbar

# Tidy up layout
plt.tight_layout()
plt.savefig('jaccard.svg')
plt.show()
#plt.savefig('jaccard.svg')
plt.show()


# In[70]:


import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# Get the original reversed colormap
cmap_orig = plt.cm.get_cmap('RdYlBu_r')

# Function to truncate colormap
def truncate_colormap(cmap, minval=0.5, maxval=1.0, n=256):
    colors = cmap(np.linspace(minval, maxval, n))
    return mcolors.LinearSegmentedColormap.from_list(f'trunc({cmap.name},{minval},{maxval})', colors)

# Create the lighter blue colormap (skip the darkest part)
cmap_lighter_blue = truncate_colormap(cmap_orig, 0.5, 1.0)


# In[71]:


#Jaccard with color categories
# 1. Reorder categorical variables in adata.obs
adata.obs['class'] = adata.obs['class'].astype('category').cat.reorder_categories([
   'Progenitor', 'MGE_LHX6/MAF', 'MGE_LHX6/NPY', 'MGE_CRABP1/MAF','MGE_CRABP1/TAC3', 'MGE_CRABP1/TH','LGE_MEIS2/PAX6',
   'LGE_FOXP2/TSHZ1', 'LGE_FOXP1/PENK', 'LGE_FOXP1/ISL1', 'LGE_FOXP1/ISL1/NPY1R', 'CGE_NR2F2/PROX1',
   'VMF_CRABP1/LHX8','VMF_ZIC1/ZIC2', 'RMTW_ZIC1/RELN'
], ordered=True)

adata.obs['integrated_class'] = adata.obs['integrated_class'].cat.reorder_categories([
    'Progenitor', 'MGE_LHX6/MAF', 'MGE_LHX6/NPY', 'MGE_CRABP1/MAF', 'MGE_CRABP1/TAC3',
    'LGE_MEIS2/PAX6', 'LGE_FOXP2/TSHZ1', 'LGE_FOXP1/PENK', 'LGE_FOXP1/ISL1',
    'LGE_FOXP1/ISL1/NPY1R', 'CGE_NR2F2/PROX1', 'VMF_CRABP1/LHX8', 'RMTW_ZIC1/RELN'
], ordered=True)

adata.obs['species'] = adata.obs['species'].cat.reorder_categories([
    'Macaque', 'Rabbit', 'Rat', 'Mouse', 'Pig', 'Ferret', 'Sugarglider', 'Opossum'
], ordered=True)

# 2. Combine class/species strings
adata.obs['original_classxspecies'] = (
    adata.obs['class'].astype(str) + ' | ' + adata.obs['species'].astype(str)
)
adata.obs['integrated_classxspecies'] = (
    adata.obs['integrated_class'].astype(str) + ' | ' + adata.obs['species'].astype(str)
)

# 3. Create full ordered lists (species-major)
species_order = list(adata.obs['species'].cat.categories)
class_order = list(adata.obs['class'].cat.categories)
integrated_class_order = list(adata.obs['integrated_class'].cat.categories)

original_order = [f"{cl} | {sp}" for sp in species_order for cl in class_order]
integrated_order = [f"{cl} | {sp}" for sp in species_order for cl in integrated_class_order]

# 4. Build sets of cell barcodes per group
original_groups = adata.obs['original_classxspecies'].unique()
integrated_groups = adata.obs['integrated_classxspecies'].unique()

original_sets = {
    label: set(adata.obs_names[adata.obs['original_classxspecies'] == label])
    for label in original_groups
}
integrated_sets = {
    label: set(adata.obs_names[adata.obs['integrated_classxspecies'] == label])
    for label in integrated_groups
}

# 5. Compute Jaccard matrix
jaccard_matrix = pd.DataFrame(index=original_groups, columns=integrated_groups)

for o_label, i_label in product(original_groups, integrated_groups):
    intersection = original_sets[o_label] & integrated_sets[i_label]
    union = original_sets[o_label] | integrated_sets[i_label]
    jaccard_matrix.loc[o_label, i_label] = len(intersection) / len(union) if union else np.nan

jaccard_matrix = jaccard_matrix.astype(float)

# 6. Reorder matrix rows and columns by your specified orders
ordered_rows = [x for x in original_order if x in jaccard_matrix.index]
ordered_cols = [x for x in integrated_order if x in jaccard_matrix.columns]

jaccard_matrix = jaccard_matrix.loc[ordered_rows, ordered_cols]

# 7. Prepare annotation data for rows and columns
row_df = pd.Series(jaccard_matrix.index).str.split(" \| ", n=1, expand=True)
row_df.columns = ['Class', 'Species']
col_df = pd.Series(jaccard_matrix.columns).str.split(" \| ", n=1, expand=True)
col_df.columns = ['Class', 'Species']

annot_rows = pd.DataFrame({
    'Class': row_df['Class'].values,
    'Species': row_df['Species'].values
}, index=jaccard_matrix.index)

annot_cols = pd.DataFrame({
    'Class': col_df['Class'].values,
    'Species': col_df['Species'].values
}, index=jaccard_matrix.columns)

# 8. Create annotation colors dict (use your paldict and manual_color_assignment)
annotation_colors = {
    'Species': manual_color_assignment,
    'Class': paldict
}

# 9. Create HeatmapAnnotations
row_ha = HeatmapAnnotation(
    Class=anno_simple(df=annot_rows['Class'], colors=paldict, height=5),
    Species=anno_simple(df=annot_rows['Species'], colors=manual_color_assignment),
    axis=0
)

col_ha = HeatmapAnnotation(
    Class=anno_simple(df=annot_cols['Class'], colors=paldict, height=5),
    Species=anno_simple(df=annot_cols['Species'], colors=manual_color_assignment),
    axis=1
)

# Reorder matrix rows and columns explicitly (you already have this)
ordered_rows = [x for x in original_order if x in jaccard_matrix.index]
ordered_cols = [x for x in integrated_order if x in jaccard_matrix.columns]

jaccard_matrix = jaccard_matrix.loc[ordered_rows, ordered_cols]

# Then plot WITHOUT row_order or col_order arguments:
plt.figure(figsize=(20, 20))

cm = ClusterMapPlotter(
    data=jaccard_matrix,
    left_annotation=row_ha,
    top_annotation=col_ha,
    cmap=cmap_lighter_blue, #YlOrRd
    vmin=0,
    vmax=1,
    row_dendrogram=False,
    col_dendrogram=False,
    row_cluster=False,
    col_cluster=False,
    rasterized=False
)

plt.title("Jaccard Index: Class×Species vs IntegratedClass×Species (Species-major order)", fontsize=20)
plt.tight_layout()
plt.savefig('jaccard_index_annotated.svg', bbox_inches='tight')
plt.show()


# In[8]:


Cluster = adata[adata.obs['integrated_class'].isin(['MGE_CRABP1/MAF', 'MGE_CRABP1/TAC3']) == True]


# In[73]:


TAC3Cluster = adata[adata.obs['integrated_class'].isin(['MGE_CRABP1/TAC3']) == True]


# In[106]:


TAC3Cluster.obs['species_region'].value_counts()


# In[73]:


sc.tl.dendrogram(Cluster, "integrated_classxspecies")
sc.pl.dendrogram(Cluster, "integrated_classxspecies")


# In[275]:


#marker genes based on Schmitz et al. 2022
marker_genes_dict = {

    'MGE_CRABP1/MAF': ['CRABP1','ETV1','COL19A1', 'RBP4','MAF', 'MAFB','PTHLH'],
    'MGE_CRABP1/TAC3': ['KIT','TRHDE','TAC3', 'TH','CHRNA3','CHRNA5','STXBP6', 'LHX8', 'CER1'],

}


sc.pl.matrixplot(Cluster, groupby=['integrated_class'], var_names=marker_genes_dict, use_raw=False, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes =False, save='matrix.svg')



# In[273]:


#marker genes based on Schmitz et al. 2022
marker_genes_dict = {

    'MGE_CRABP1/MAF': ['CRABP1','ETV1', 'KIT','MAF', 'MAFB','COL19A1','RBP4', 'SP9', 'TIAM2'],
    'MGE_CRABP1/TAC3': ['PTPRK',  'RGS7', 'TRHDE', 'TAC3', 'TH','CHRNA3','STXBP6', 'LHX8', 'CER1', 'ZIC1', 'ZIC2', 'ZIC4', 'CACNA2D3']

}
#'DCC', 'GAP43', 'PPP2R2B', 'RALYL', 'CACNA2D3',
# 'SYT1', 'MEF2C', 'BCL11B', 'SOX1', 'ARID1B'
# 'NFIA', 'CNTN5', 'LINGO2', 'ZIC1', 'ZIC2', 'ZIC4', 'KITLG'

sc.pl.dotplot(Cluster, groupby=['integrated_class'] ,var_names=marker_genes_dict, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes =True, save='matrix.svg')



# In[75]:


#marker genes based on Schmitz et al. 2022
marker_genes_dict = {

    'MGE_CRABP1/MAF': ['CRABP1','ETV1', 'KIT','MAF', 'MAFB','COL19A1','RBP4', 'SP9', 'TIAM2'],
    'MGE_CRABP1/TAC3': ['PTPRK',  'RGS7', 'TRHDE', 'TAC3', 'TH','CHRNA3','STXBP6', 'LHX8', 'CER1', 'ZIC1', 'ZIC2', 'ZIC4', 'CACNA2D3']

}
#'DCC', 'GAP43', 'PPP2R2B', 'RALYL', 'CACNA2D3',
# 'SYT1', 'MEF2C', 'BCL11B', 'SOX1', 'ARID1B'
# 'NFIA', 'CNTN5', 'LINGO2', 'ZIC1', 'ZIC2', 'ZIC4', 'KITLG'

sc.pl.dotplot(Cluster, groupby=['integrated_class','species'] ,var_names=marker_genes_dict, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes =True, save='matrix.svg')



# In[77]:


# Define the species you want to include
species_subset = ['Macaque', 'Pig', 'Ferret', 'Sugarglider', 'Opossum']  # replace with your desired species

# Subset the adata to only include those species
subset_Cluster = Cluster[Cluster.obs['species'].isin(species_subset), :]

#marker genes based on Schmitz et al. 2022
marker_genes_dict = {

    'MGE_CRABP1/MAF': ['CRABP1','ETV1', 'KIT','MAF', 'MAFB','COL19A1','RBP4', 'SP9', 'TIAM2'],
    'MGE_CRABP1/TAC3': ['PTPRK',  'RGS7', 'TRHDE', 'TAC3', 'TH','CHRNA3','STXBP6', 'LHX8', 'CER1', 'ZIC1', 'ZIC2', 'ZIC4', 'CACNA2D3']

}
#'DCC', 'GAP43', 'PPP2R2B', 'RALYL', 'CACNA2D3',
# 'SYT1', 'MEF2C', 'BCL11B', 'SOX1', 'ARID1B'
# 'NFIA', 'CNTN5', 'LINGO2', 'ZIC1', 'ZIC2', 'ZIC4', 'KITLG'

sc.pl.dotplot(subset_Cluster, groupby=['integrated_class','species'] ,var_names=marker_genes_dict, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes =False, save='matrixCRABP1.svg')



# In[9]:


sc.set_figure_params(fontsize = 20)
# Define the species you want to include
species_subset = ['Ferret']  # replace with your desired species

# Subset the adata to only include those species
subset_Cluster = Cluster[Cluster.obs['species'].isin(species_subset), :]

#marker genes based on Schmitz et al. 2022
marker_genes_dict = {

    'MGE_CRABP1/MAF': ['NKX2-1' , 'LHX6','CRABP1','ETV1', 'KIT','MAF', 'MAFB','COL19A1','RBP4', ],
    'MGE_CRABP1/TAC3': ['TRHDE', 'TAC3', 'TH','CHRNA3','STXBP6', 'LHX8', 'ZIC1', 'ZIC2', 'ZIC4']

}
#'DCC', 'GAP43', 'PPP2R2B', 'RALYL', 'CACNA2D3',
# 'SYT1', 'MEF2C', 'BCL11B', 'SOX1', 'ARID1B'
# 'NFIA', 'CNTN5', 'LINGO2', 'ZIC1', 'ZIC2', 'ZIC4', 'KITLG'

sc.pl.matrixplot(subset_Cluster, groupby=['integrated_class','species'] ,var_names=marker_genes_dict, vmin=-2, vmax=2, cmap='PiYG_r',  figsize=(12, 1.5),swap_axes =False, save='matrixferretCRABP1.svg')



# In[ ]:





# In[80]:


# Define the species you want to include
species_subset = ['Macaque', 'Pig', 'Ferret', 'Sugarglider', 'Opossum']  # replace with your desired species

# Subset the adata to only include those species
subset_TAC3Cluster = TAC3Cluster[TAC3Cluster.obs['species'].isin(species_subset), :]

# Plot as before
sc.pl.dotplot(
    subset_TAC3Cluster, 
    groupby=['integrated_class','species'], 
    var_names=marker_genes_dict3, 
    vmin=-2, 
    vmax=2, 
    cmap='PiYG_r', 
    swap_axes=False, 
    save='matrix_Tac3subset.svg'
)


# In[84]:


# Define the species you want to include
sc.set_figure_params(fontsize = 10)
species_subset = ['Macaque', 'Rabbit', 'Rat', 'Mouse']  # replace with your desired species

# Subset the adata to only include those species
subset_TAC3Cluster = TAC3Cluster[TAC3Cluster.obs['species'].isin(species_subset), :]

# Plot as before
sc.pl.dotplot(
    subset_TAC3Cluster, 
    groupby=['integrated_class','species'], 
    var_names=marker_genes_dict3, 
    vmin=-2, 
    vmax=2, 
    cmap='PiYG_r', 
    swap_axes= False, 
    save='matrix_Tac3subsetrodents.svg'
)


# In[79]:


marker_genes_dict3 = {

    'MGE_CRABP1/TAC3': ['NKX2-1', 'LHX6', 'CRABP1', 'ETV1','TRHDE', 'CHRNA3', 'STXBP6', 'TAC3', 'TH',, 'LHX8', 'ZIC1', 'ZIC2', 'ZIC4']

}

# Define the species you want to include
sc.set_figure_params(fontsize = 10)
species_subset = ['Macaque', 'Rabbit', 'Rat', 'Mouse','Pig', 'Ferret', 'Sugarglider', 'Opossum']  # replace with your desired species

# Subset the adata to only include those species
subset_TAC3Cluster = TAC3Cluster[TAC3Cluster.obs['species'].isin(species_subset), :]

# Plot as before
sc.pl.dotplot(
    subset_TAC3Cluster, 
    groupby=['integrated_class','species'], 
    var_names=marker_genes_dict3, 
    vmin=-2, 
    vmax=2, 
    cmap='PiYG_r', 
    swap_axes= False, 
    save='matrix_Tac3subsetreviewerfig.svg'
)


# In[92]:


#marker genes based on Schmitz et al. 2022
marker_genes_dict3 = {

    'MGE_CRABP1/TAC3': ['NKX2-1', 'LHX6', 'CRABP1', 'ETV1', 'TRHDE','CHRNA3','STXBP6', 'TAC3', 'TH', 'LHX8', 'CER1', 'ZIC1', 'ZIC2', 'ZIC4']

}
#'DCC', 'GAP43', 'PPP2R2B', 'RALYL', 'CACNA2D3',
# 'SYT1', 'MEF2C', 'BCL11B', 'SOX1', 'ARID1B'
# 'NFIA', 'CNTN5', 'LINGO2', 'ZIC1', 'ZIC2', 'ZIC4', 'KITLG'

sc.pl.dotplot(TAC3Cluster, groupby=['integrated_class','species'] ,var_names=marker_genes_dict3, vmin=-2, vmax=2, cmap='PiYG_r', swap_axes =True, save='matrix.svg')



# In[93]:


#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(Cluster, groupby = 'integrated_class', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(Cluster.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df


# In[81]:


sc.set_figure_params(fontsize = 10)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(Cluster, color = ['species'], groups = ['Macaque','Pig', 'Ferret', 'Sugarglider', 'Opossum'],frameon=True,size= 50,wspace=0.25, ncols = 2, save = 'integratedspeciesCRABP1cluster.svg')


# In[276]:


# Check for clusters due to technical issues.

sc.set_figure_params(fontsize = 65)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(Cluster, color = plot_keys, groups = plot_keys, palette = manual_color_assignment,
           wspace=0,ncols = 3, frameon= False, size = 25, legend_loc = None, save = 'clusterspecieswarhol.svg')
#sc.pl.violin(adata, ['species'], groupby = "leiden_subcluster")


# In[263]:


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(Cluster, color = ['integrated_class'], size= 50,wspace=0.25, ncols = 2)


# In[207]:


sc.set_figure_params(fontsize = 50)
plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(Cluster, color = ['NKX2-1','CRABP1','COL19A1', 'MAF', 'TAC3', 'CHRNA3'], use_raw = False,frameon = False, size= 50, wspace=0.25, ncols = 2, cmap ='PuRd', save = 'TAC3CLUSTER_MARKERS.svg')


# In[83]:


# --- Setup plotting ---
sc.set_figure_params(fontsize=30)
plt.rcParams['axes.grid'] = False

# --- Setup ---
species_name = "Mouse"
class_subset = ["MGE_CRABP1/TH", "MGE_CRABP1/MAF"]

# Add combined labels if not already present
adata.obs['original_classxspecies'] = adata.obs['class'].astype(str) + "_" + adata.obs['species'].astype(str)

# --- Step 1: Subset adata for the selected original_classxspecies groups ---
subset_labels = [f"{cls}_{species_name}" for cls in class_subset]
adata_subset = adata[adata.obs['original_classxspecies'].isin(subset_labels)].copy()

# --- Step 2: Extract PCA ---
X_pca_subset = adata_subset.obsm['X_pca']

# Exclude Mouse from integrated class grouping
non_mouse_mask = adata.obs['species'] != "Mouse"
adata_non_mouse = adata[non_mouse_mask]
X_pca_non_mouse = adata_non_mouse.obsm['X_pca']

# --- Step 3: Group and average PCA values ---
grouped_original_subset = pd.DataFrame(
    X_pca_subset, index=adata_subset.obs['original_classxspecies']
).groupby(level=0).mean()

grouped_integrated_all = pd.DataFrame(
    X_pca_non_mouse, index=adata_non_mouse.obs['integrated_class']
).groupby(level=0).mean()

# --- Step 4: Compute Pearson correlation between group means ---
corr_values = np.corrcoef(
    grouped_original_subset.values, grouped_integrated_all.values
)

n_orig = grouped_original_subset.shape[0]
n_integrated = grouped_integrated_all.shape[0]

corr_matrix_final = pd.DataFrame(
    corr_values[:n_orig, n_orig:],
    index=grouped_original_subset.index,
    columns=grouped_integrated_all.index
)

# --- Step 5: Plot heatmap ---
plt.figure(figsize=(20, 7))
sns.heatmap(
    corr_matrix_final,
    cmap='RdYlBu_r',
    center=0,
    vmin=-1,
    vmax=1,
    linewidths=0.25,
    linecolor='gray',
)

plt.xticks(rotation=90, fontsize=20)
plt.yticks(rotation=0, fontsize=30)
plt.tight_layout()
plt.savefig("subset_vs_nonmouse_integrated_pca_correlationnmouse.svg")
plt.show()


# In[48]:


#Find marker genes for each leiden cluster.
sc.tl.rank_genes_groups(adata, groupby = 'integrated_class', method = 'wilcoxon', use_raw=False)
# Show the top ranked genes per cluster in a dataframe.
small_degs_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)
pd.set_option('display.max_columns', 500)
small_degs_df


# In[ ]:


plt.rcParams["figure.figsize"] = (10, 10)
sc.pl.umap(Cluster, color = ['TH', 'TAC3', 'CHRNA3', 'LHX8', 'ZIC1', 'ZIC4', 'STXBP6', 'RBP4','COL19A1',  'PTHLH'],  use_raw =False,  frameon = False, color_map = "PuRd", size=50, wspace=0.25, ncols = 3)


# In[349]:


adata.write('/pathto/concat_labeled.h5ad')


# In[5]:

#cleaning up for browsable data on ucsc 
for col in ['leiden_subcluster_CRABPY', 'leiden_subcluster_CRABPY2','leiden_subcluster_CRABPY3', 'file_name', 'hires_leiden', 'latent_cell_probability', 'phase', 'latent_time', 'region_fraction','GEMwell', 'leiden_subcluster', 'leiden_subcluster2', 'leiden_subcluster3', 'leiden_subcluster4' ]:
    if col in adata.obs:
        del adata.obs[col]

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
obs.to_csv('/pathto/geo/integrated_barcodes.tsv', sep='\t', index=True)

# Save var.tsv
var.to_csv('/pathto/geo/integrated_genes.tsv', sep='\t', index=True)

# Save matrix.mtx
scipy.io.mmwrite('/pathto/geo/integrate_matrix.mtx', matrix)

