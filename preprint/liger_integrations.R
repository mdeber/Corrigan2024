

# Packages, functions -----------------------------------------------------

suppressPackageStartupMessages({
    library(dplyr)
    library(reshape2)
    library(ggplot2)
    library(patchwork)
    library(Seurat) # v5.0.3
    library(SeuratWrappers)
    library(rliger) # v1.0.1
    library(funr) # for finding local paths
})

if (!"readxl" %in% installed.packages())
    install.packages("readxl")

### Utility functions

# for converting from primate to mouse gene symbols
toCapCase <- function(x) {
    .cap <- function(s) paste0(
        toupper(substr(s, 1, 1)),
        tolower(substr(s, 2, nchar(s)))
    )
    unname(sapply(x, .cap))
}

theme_md_classic <- function() {
    theme_classic(
        base_size = 10, 
        base_family = "Helvetica",
        base_line_size = 0.25,
        base_rect_size = 0.25
    ) %+replace% theme(
        axis.text = element_text(color = "black", size = 8),
        axis.ticks = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        complete = TRUE
    )
}

theme_md_bw <- function() {
    theme_bw(
        base_size = 10, 
        base_family = "Helvetica",
        base_line_size = 0.25,
        base_rect_size = 0.25
    ) %+replace% theme(
        axis.text = element_text(color = "black", size = 8),
        axis.ticks = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.text.y = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        complete = TRUE
    )
}


# Paths -------------------------------------------------------------------

# 
# Note: 
# funr::get_script_path() works when running script from the command line
# (`Rscript <PATH_TO_THIS_SCRIPT>`) or if sourcing in R
# (`source(<PATH_TO_THIS_SCRIPT>)`), but it won't work if you are running this
# line-by-line. In that case, manually substitute your local path to this repo
# 

path_parent <- dirname(funr::get_script_path())
path_out <- file.path(path_parent, "preprint/liger_integrations")

if (!dir.exists(path_out))
    dir.create(path_out)

# Data ---
path_abc_cl08_str_drop057 <- file.path(
    path_parent, "abc/data/abc_seurat_cl08_str_dropSubc057_0p8",
    "abc_seurat_cl08_str_dropSubc057_0p8.rds"
)

path_marm_str_gad_drop2_dropGlut <- file.path(
    path_parent, "marm_census/data/marm_str_gad_drop2_dropGlut_0p8",
    "marm_str_gad_drop2_dropGlut_0p8.rds"
)

path_siletti_recl_dropAllMeis2 <- file.path(
    path_parent, "siletti/data/siletti_neur_roi_sclust_dropAllMeis2_0p8",
    "siletti_neur_roi_sclust_dropAllMeis2_0p8.rds"
)

# Siletti atlas metadata ---
path_siletti_clust_terms <- file.path(
    path_parent, "siletti/data/cluster_annotation.xlsx"
)

path_siletti_subclust_terms <- file.path(
    path_parent, "siletti/data/subcluster_annotation.xlsx"
)

# Ortholog tables ---
# 
# In each, 'Gene.name' is marmoset gene name (with other columns for 
# 'Mouse.gene.name' or 'Human.gene.name')
path_marm_human_ortho <- file.path(
    path_parent, "ortholog_tables/marm_human_orthologs.csv"
)

path_marm_mouse_ortho <- file.path(
    path_parent, "ortholog_tables/marm_mouse_orthologs.csv"
)


# Manual cluster annotations (with mappings to the Seurat clusters) ---
path_abc_clust_anno <- file.path(
    path_parent, "abc/data/abc_seurat_cl08_str_dropSubc057_0p8",
    "abc_seurat_cl08_str_recl_drop057_clust_anno.csv"
)
path_marm_clust_anno <- file.path(
    path_parent, "marm_census/data/marm_str_gad_drop2_dropGlut_0p8",
    "marm_str_gad_drop2_dropGlut_clust_anno.csv"
)
path_siletti_clust_anno <- file.path(
    path_parent, "siletti/data/siletti_neur_roi_sclust_dropAllMeis2_0p8",
    "siletti_neur_roi_sclust_dropAllMeis2_0p8_clust_anno.csv"
)


# Imports -----------------------------------------------------------------

abc_seurat_cl08_str_recl_drop057 <- readRDS(path_abc_cl08_str_drop057)
marm_str_gad_drop2_dropGlut <- readRDS(path_marm_str_gad_drop2_dropGlut)
siletti_neur_roi_sclust_0p8_dropAllMeis2 <- readRDS(path_siletti_recl_dropAllMeis2)

siletti_clust_terms <- readxl::read_excel(path_siletti_clust_terms)
siletti_subclust_terms <- readxl::read_excel(path_siletti_subclust_terms)

df_ortho_marm_human <- read.csv(path_marm_human_ortho)
df_ortho_marm_mouse <- read.csv(path_marm_mouse_ortho)

abc_clust_anno <- read.csv(path_abc_clust_anno)
marm_clust_anno <- read.csv(path_marm_clust_anno)
siletti_clust_anno <- read.csv(path_siletti_clust_anno)

# (make cluster numbers factors to match the data objects)
abc_clust_anno$seurat_clusters <- as.factor(
    abc_clust_anno$seurat_clusters
)
marm_clust_anno$seurat_clusters <- as.factor(
    marm_clust_anno$seurat_clusters
)
siletti_clust_anno$seurat_clusters <- as.factor(
    siletti_clust_anno$seurat_clusters
)

#
# Note the siletti taxonomy labels (`siletti_clust_terms` and 
# `siletti_subclust_terms`) aren't being added to the data. They are joined 
# lazily as needed in analysis.
#


# Prepare manual cluster annotations for join -----------------------------

# These are being prepared to later join into the combined multi-species
# integrated objects
abc_clust_anno$species <- "mouse"
marm_clust_anno$species <- "marmoset"
siletti_clust_anno$species <- "human"
clust_anno <- rbind(abc_clust_anno, marm_clust_anno, siletti_clust_anno)
names(clust_anno) <- sub("^seurat_clusters$", "RNA_snn_res.0.8", names(clust_anno))

# > clust_anno
# RNA_snn_res.0.8 seurat_clusters_anno  species
# 1                0              Sst/Npy    mouse
# 2                1           Lhx6/Prox1    mouse
# 3                2          Pthlh/Pvalb    mouse
# 4                3                   Th    mouse
# 5                4          Pthlh/Pvalb    mouse
# 6                5          Pthlh/Pvalb    mouse
# 7                6              Sst/Npy    mouse
# 8                7                   Th    mouse
# 9                8                 Chat    mouse
# 10               9                   Th    mouse
# 11              10              Sst/Npy    mouse
# 12              11              Sst/Npy    mouse
# 13              12              Sst/Npy    mouse
# 14              13              Sst/Npy    mouse
# 15              14                 Chat    mouse
# 16              15                 Chat    mouse
# 17               0              SST/NPY marmoset
# 18               1          PTHLH/PVALB marmoset
# 19               2                 TAC3 marmoset
# 20               3                 TAC3 marmoset
# 21               4          PTHLH/PVALB marmoset
# 22               5          PTHLH/PVALB marmoset
# 23               6                 CHAT marmoset
# 24               7          PTHLH/PVALB marmoset
# 25               8                 TAC3 marmoset
# 26               9                 TAC3 marmoset
# 27              10              SST/NPY marmoset
# 28              11                  CCK marmoset
# 29              12                Mixed marmoset
# 30              13                  CCK marmoset
# 31              14                Mixed marmoset
# 32              15                 CHAT marmoset
# 33               0                 TAC3    human
# 34               1                 TAC3    human
# 35               2          PTHLH/PVALB    human
# 36               3          PTHLH/PVALB    human
# 37               4          PTHLH/PVALB    human
# 38               5                 TAC3    human
# 39               6                 TAC3    human
# 40               7          PTHLH/PVALB    human
# 41               8          PTHLH/PVALB    human
# 42               9              SST/NPY    human
# 43              10         Hypothalamus    human
# 44              11                Mixed    human
# 45              12              CCK/VIP    human
# 46              13                 CHAT    human
# 47              14                  CCK    human
# 48              15              SST/NPY    human

# Find gene list for integration ------------------------------------------

# need 3-way homology, will center on marmoset; let's get set of genes found in 
# the datasets and both homology tables

# dim(marm_str_gad_drop2_recl_0p8)[1] # n=19762 marm genes at start

marm_gene_in_df_ortho_human <- intersect(
    rownames(marm_str_gad_drop2_dropGlut),
    df_ortho_marm_human$Gene.name
) # n=11859 marm genes in the human homology table

marm_gene_in_df_ortho_mouse <- intersect(
    rownames(marm_str_gad_drop2_dropGlut),
    df_ortho_marm_mouse$Gene.name
) # n=10024 marm genes in mouse homology table

marm_gene_in_df_ortho_both <- intersect(
    marm_gene_in_df_ortho_human,
    marm_gene_in_df_ortho_mouse
) # n=10020 marm genes in both homology tables


# Now let's go back and find, from those genes, the ones also found in the human 
# and mouse datasets as well.

human_gene_in_marm_df_ortho_both <- df_ortho_marm_human %>% 
    subset(Gene.name %in% marm_gene_in_df_ortho_both) %>% 
    .$Human.gene.name %>% 
    intersect(
        rownames(siletti_neur_roi_sclust_0p8_dropAllMeis2)
    ) # n=10000 human genes in dataset that match those marmoset orthologs

mouse_gene_in_marm_df_ortho_both <- df_ortho_marm_mouse %>% 
    subset(Gene.name %in% marm_gene_in_df_ortho_both) %>% 
    .$Mouse.gene.name %>% 
    intersect(
        rownames(abc_seurat_cl08_str_recl_drop057)
    ) # n=9930 mouse genes in dataset that match those marmoset orthologs


# Now start subsetting the ortholog table

# n=9876 final genes found in all datasets with homologs in the homology table
df_ortho_3way <- df_ortho_marm_human %>% 
    subset(Gene.name %in% marm_gene_in_df_ortho_both) %>% 
    subset(Human.gene.name %in% human_gene_in_marm_df_ortho_both) %>% 
    subset(Mouse.gene.name %in% mouse_gene_in_marm_df_ortho_both)

# only object needed here is `df_ortho_3way`
rm(
    marm_gene_in_df_ortho_human,
    marm_gene_in_df_ortho_mouse,
    marm_gene_in_df_ortho_both,
    human_gene_in_marm_df_ortho_both,
    mouse_gene_in_marm_df_ortho_both
)

# Add back TAC3 (mouse Tac2 not called a homolog in current ensembl release)
df_ortho_3way <- df_ortho_3way %>% 
    rbind(
        data.frame(
            "Gene.name" = "TAC3",
            "Gene.stable.ID.version" = NA,
            "Human.gene.stable.ID" = NA,
            "Human.gene.name" = "TAC3",
            "Human.homology.type" = NA,
            "X.id..target.Human.gene.identical.to.query.gene" = NA,
            "X.id..query.gene.identical.to.target.Human.gene" = NA,
            "Human.Gene.order.conservation.score" = NA,
            "Human.orthology.confidence..0.low..1.high." = NA,
            "Mouse.gene.stable.ID" = NA,
            "Mouse.gene.name" = "Tac2",
            "Mouse.homology.type" = NA,
            "X.id..target.Mouse.gene.identical.to.query.gene" = NA,
            "X.id..query.gene.identical.to.target.Mouse.gene" = NA,
            "Mouse.Gene.order.conservation.score" = NA,
            "Mouse.orthology.confidence..0.low..1.high." = NA
        )
    ) %>% 
    .[order(.$Gene.name), ] 

# > nrow(df_ortho_3way)
# [1] 9877

# Drop a duplicated gene

# > df_ortho_3way$Mouse.gene.name[which(duplicated(df_ortho_3way$Mouse.gene.name))]
# [1] "Zcchc12"

df_ortho_3way <- subset(df_ortho_3way, Mouse.gene.name != "Zcchc12")

# > nrow(df_ortho_3way)
# [1] 9875



# Subset datasets for 3-way orthologs -------------------------------------

marm_subset <- subset(
    marm_str_gad_drop2_dropGlut,
    features = df_ortho_3way$Gene.name
)

abc_subset <- subset(
    abc_seurat_cl08_str_recl_drop057,
    features = df_ortho_3way$Mouse.gene.name
)

siletti_subset <- subset(
    siletti_neur_roi_sclust_0p8_dropAllMeis2,
    features = df_ortho_3way$Human.gene.name
)



# Recreate objects with renamed genes -------------------------------------

# They have to have the same gene symbols, and in the same order (the latter 
# of which Seurat::merge will do)

# ABC
mouse_rename_list <- setNames(
    df_ortho_3way$Gene.name,
    df_ortho_3way$Mouse.gene.name
)

mat <- abc_subset[["RNA"]]$counts
rownames(mat) <- unname(mouse_rename_list[rownames(mat)])

abc_subset <- CreateSeuratObject(
    counts = mat,
    meta.data = abc_subset@meta.data
)

# Siletti
human_rename_list <- setNames(
    df_ortho_3way$Gene.name,
    df_ortho_3way$Human.gene.name
)

mat <- siletti_subset[["RNA"]]$counts
rownames(mat) <- unname(human_rename_list[rownames(mat)])

siletti_subset <- CreateSeuratObject(
    counts = mat,
    meta.data = siletti_subset@meta.data
)

# For consistency, also remake the marmoset
marm_subset <- CreateSeuratObject(
    counts = marm_subset[["RNA"]]$counts,
    meta.data = marm_subset@meta.data
)

stopifnot(all(Features(abc_subset) %in% Features(marm_subset)))
stopifnot(all(Features(siletti_subset) %in% Features(marm_subset)))
rm(mat)


# Add species info and merge ----------------------------------------------

marm_subset$species <- "marmoset"
abc_subset$species <- "mouse"
siletti_subset$species <- "human"

seurat_merge <- merge(
    marm_subset,
    list(siletti_subset, abc_subset)
)

# > seurat_merge
# An object of class Seurat 
# 9875 features across 13582 samples within 1 assay 
# Active assay: RNA (9875 features, 0 variable features)
# 3 layers present: counts.1, counts.2, counts.3

seurat_merge <- JoinLayers(seurat_merge)


# Requisite Seurat pipeline -----------------------------------------------

# What's required before LIGER integration
seurat_merge <- NormalizeData(seurat_merge, verbose = FALSE)
seurat_merge <- FindVariableFeatures(seurat_merge, verbose = FALSE)
seurat_merge <- ScaleData(
    seurat_merge, split.by = "species", do.center = FALSE, verbose = FALSE
)


# Find mouse Tac2+ barcodes -----------------------------------------------

tac2_barcodes <- abc_seurat_cl08_str_recl_drop057 %>% 
    (function(x) {
        Cells(x)[
            (x$subclass == "055 STR Lhx8 Gaba") & (x[["RNA"]]$data["Tac2", ] > 0)
        ]
    })

# Those looked suitable in plots like this:
# 
# abc_seurat_cl08_str_recl_drop057 %>% 
#     (function(x) {
#         x$is_Tac2 = (x$subclass == "055 STR Lhx8 Gaba") & 
#             (x[["RNA"]]$data["Tac2", ] > 0)
#         x$is_Tac2 <- ifelse(x$is_Tac2, "Tac2+", "Tac2-")
#         x
#     }) %>% 
#     subset(subclass == "055 STR Lhx8 Gaba") %>% 
#     DimPlot(group.by = "is_Tac2") + 
#     labs(title = "Tac2+ cells in ABC subclass 055")


# Marker genes for joint embedding ----------------------------------------

# Marker genes that are still found in the 3-way ortholog map
marker_genes <- c(
    "ANGPT2",
    "ARX",
    "CHAT",
    "CHODL",
    "CHRNA3",
    "CHRNA4",
    "CHRNA7",
    "CKB",
    "COL19A1",
    "CRABP1",
    "ETV1",
    "EYA2",
    "FOXG1",
    "FOXP1",
    "FOXP2",
    "GAD1",
    "GBX1",
    "GRIK1",
    "ISL1",
    "KIT",
    "LHX6",
    "LHX8",
    "MAF",
    "MAF1",
    "MEF2C",
    "MEIS2",
    "NKX2-1",
    "NPY",
    "NPY1R",
    "NR2F2",
    "NRTN",
    "PARM1",
    "PAX6",
    "PDZRN3",
    "PENK",
    "PRLR",
    "PROX1",
    "PTHLH",
    "PTPRZ1",
    "PVALB",
    "RARB",
    "RBP4",
    "RXRG",
    "SIX3",
    "SLC17A6",
    "SLC17A7",
    "SLC35D3",
    "SOX6",
    "SST",
    "SYT1",
    "TAC3",
    "TACR3",
    "TH",
    "TRHDE",
    "TSHZ1",
    "VIP"
)


# Find optimized k range for LIGER ----------------------------------------

# This section can be skipped entirely. It takes a while to run and isn't 
# required to run anything else below.

if (FALSE) {
    liger_obj <- seuratToLiger(list(
        marmoset = subset(seurat_merge, species == "marmoset"),
        human    = subset(seurat_merge, species == "human"),
        mouse    = subset(seurat_merge, species == "mouse")
    ))
    liger_obj <- scaleNotCenter(liger_obj)
    ggsave(
        file.path(path_out, "liger_suggestK_lambda5.pdf"),
        plot = suggestK(liger_obj), # with lambda=5
        width = 100, height = 90, units = "mm"
    )
    rm(liger_obj)
}


# 500 features: Parameter sweep -------------------------------------------

seur500 <- seurat_merge %>% 
    FindVariableFeatures(nfeatures = 500) %>% 
    ScaleData(split.by = "species", do.center = FALSE)

# n=42 param combos
liger_param_sweep <- expand.grid(
    k = seq(20, 45, 5),
    lambda = 2^seq(-1, 5, 1) 
) %>% 
    as.list %>% 
    (function(x) Map(c, x$k, x$lambda)) %>% 
    (function(x) setNames(
        x, sapply(x, function(x) paste0("k", x[1], "_lambda", x[2]))
    )) %>% 
    lapply(function(x) {
        seur <- seur500
        seur@reductions <- list()
        seur <- RunOptimizeALS(
            seur, k = x[1], lambda = x[2], split.by = "species"
        )
        seur <- RunQuantileNorm(seur, split.by = "species")
        seur <- FindNeighbors(
            seur, reduction = "iNMF", dims = 1:20, verbose = FALSE
        )
        seur <- RunUMAP(
            seur, dims = 1:ncol(seur[["iNMF"]]), reduction = "iNMF", 
            verbose = FALSE
        )
        seur
    })
saveRDS(
    liger_param_sweep, 
    file.path(path_out, "liger_param_sweep_500features.rds")
)
rm(seur500)


# 500 features: Plots -----------------------------------------------------

path_sweep <- file.path(path_out, "plots_500features")

if (!dir.exists(path_sweep))
    dir.create(path_sweep)

invisible(Map(
    function(x, nm) {
        
        # --- plot species on UMAP
        path_umap_species <- file.path(path_sweep, "umap_species")
        if (!dir.exists(path_umap_species))
            dir.create(path_umap_species)
        
        DimPlot(x, group.by = "species")
        ggsave(
            file.path(path_umap_species, paste0(nm, ".pdf")),
            width = 160, height = 120, units = "mm"
        )
        
        # --- umap plot some markers together in a single file
        path_umap_combo <- file.path(path_sweep, "umap_3species_markers")
        if (!dir.exists(path_umap_combo))
            dir.create(path_umap_combo)
        
        p1 <- FeaturePlot(
            x,
            c("LHX6", "TAC3", "CHRNA3", "TH", "CRABP1"),
            split.by = "species",
            order = TRUE
        ) 
        p2 <- FeaturePlot(
            x,
            c("SST", "CHAT", "PVALB", "PTHLH", "VIP"),
            split.by = "species",
            order = TRUE
        )
        {p1 | p2} & 
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                strip.text.x = element_text(angle = 0)
            )
        ggsave(
            file.path(path_umap_combo, paste0(nm, ".pdf")),
            width = 12, height = 10, units = "in"
        )
        
        # --- plot separate umaps for all marker genes
        path_umap_markers <- file.path(path_sweep, "umap_markers")
        if (!dir.exists(path_umap_markers))
            dir.create(path_umap_markers)
        
        path_umap_param <- file.path(path_umap_markers, nm)
        if (!dir.exists(path_umap_param))
            dir.create(path_umap_param)
        
        for (gn in marker_genes) {
            FeaturePlot(x, gn, split.by = "species", order = TRUE)
            ggsave(
                file.path(path_umap_param, paste0(gn, ".pdf")),
                width = 270, height = 90, units = "mm"
            )
        }
    },
    liger_param_sweep,
    names(liger_param_sweep)
))


# 500 features: Cluster integrations --------------------------------------

# 
# note the previous clusters are present in `RNA_snn_res.0.8`;
# `seurat_clusters` is always overwritten so it's best not to use that
# 

liger_param_sweep <- lapply(
    liger_param_sweep, 
    FindClusters,
    resolution = 0.5,
    verbose = FALSE
)


# 500 features: Check number of clusters ----------------------------------

# > liger_param_sweep %>% 
#     lapply(. %>% .$RNA_snn_res.0.5 %>% nlevels) %>% 
#     # make a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
# 
# lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 19        23      23      17      18      19       18      
# k25 20        20      22      21      20      20       19      
# k30 21        25      19      23      20      24       22      
# k35 28        25      23      21      21      24       20      
# k40 22        27      23      26      22      23       26      
# k45 26        23      22      23      24      21       21   



# 500 features: Check co-clustering of Tac2 cells with themselves ---------

# First, how coherent is co-clustering of Tac2 cells?
# [out of 51 cells, how many in the top cluster?]

liger_param_sweep %>% 
    lapply(function(x) x[, Cells(x) %in% tac2_barcodes]) %>% 
    lapply(. %>% .$RNA_snn_res.0.5 %>% table %>% max) %>% 
    # make a matrix
    split(., sub("_.*", "", names(.))) %>% 
    lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
    lapply(unlist) %>% 
    do.call(rbind, .) ->
    sweep500_tac2_self_coclust

# > sweep500_tac2_self_coclust
# lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20        19      23      27      22      37       29       35
# k25        27      31      45      25      26       27       26
# k30        18      38      26      30      38       17       33
# k35        36      34      36      31      45       45       46
# k40        31      34      44      39      44       42       43
# k45        15      18      15      38      43       37       28

# As a fraction of Tac2+ cells:
# 
# > sweep500_tac2_self_coclust %>% 
#     "/"(length(tac2_barcodes)) %>% 
#     round(2)
# 
# lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20      0.37    0.45    0.53    0.43    0.73     0.57     0.69
# k25      0.53    0.61    0.88    0.49    0.51     0.53     0.51
# k30      0.35    0.75    0.51    0.59    0.75     0.33     0.65
# k35      0.71    0.67    0.71    0.61    0.88     0.88     0.90
# k40      0.61    0.67    0.86    0.76    0.86     0.82     0.84
# k45      0.29    0.35    0.29    0.75    0.84     0.73     0.55


# 500 features: Check top marm co-cluster ---------------------------------

# (Published marm cluster 7 ("07-00") has the TAC3+ striatal interneurons)

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "marmoset")) {
            subset(dsub, species == "marmoset")$CLUSTER %>% 
                table %>% 
                (function(y) names(y)[which.max(y)]) %>% 
                as.numeric
        } else {
            NA
        }
    }) -> sweep500_tac2_coclust_marm_clust

# > sweep500_tac2_coclust_marm_clust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
# lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 5         5       7       7       7       7        7       
# k25 9         NA      NA      7       7       7        7       
# k30 7         7       7       7       7       7        7       
# k35 NA        NA      7       3       NA      8        7       
# k40 4         NA      8       8       7       8        7       
# k45 9         7       7       7       7       8        8      

# > sweep500_tac2_coclust_marm_clust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
# 
# cluster freq
# 3    1
# 4    1
# 5    2
# 7   24
# 8    6
# 9    2
# <NA>    6


# 500 features: Check top marm subcluster co-cluster ----------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "marmoset")) {
            subset(dsub, species == "marmoset")$CLUSTER.SUBCLUSTER %>% 
                table %>% 
                (function(y) names(y)[which.max(y)])
        } else {
            NA
        }
    }) -> sweep500_tac2_coclust_marm_subclust


# > sweep500_tac2_coclust_marm_subclust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
# 
# lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 "05-00"   "05-00" "07-00" "07-00" "07-00" "07-00"  "07-00" 
# k25 "09-00"   NA      NA      "07-00" "07-00" "07-00"  "07-00" 
# k30 "07-00"   "07-00" "07-00" "07-00" "07-00" "07-00"  "07-00" 
# k35 NA        NA      "07-00" "03-00" NA      "08-01"  "07-00" 
# k40 "04-00"   NA      "08-01" "07-00" "07-00" "07-00"  "07-00" 
# k45 "09-00"   "07-00" "07-00" "07-00" "07-00" "08-02"  "08-02" 

# > sweep500_tac2_coclust_marm_subclust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
# 
# cluster freq
# 03-00    1
# 04-00    1
# 05-00    2
# 07-00   26
# 08-01    2
# 08-02    2
# 09-00    2
# <NA>    6



# 500 features: Check top human co-cluster --------------------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "human")) {
            subset(dsub, species == "human")$cluster_id %>% 
                table %>% 
                (function(y) names(y)[which.max(y)]) %>% 
                as.numeric
        } else {
            NA
        }
    }) -> sweep500_tac2_coclust_siletti_clust


# > sweep500_tac2_coclust_siletti_clust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
# lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 409       409     429     409     252     409      235     
# k25 288       NA      235     409     409     409      409     
# k30 258       258     409     409     409     409      409     
# k35 NA        NA      238     409     NA      NA       409     
# k40 429       392     238     409     409     409      409     
# k45 429       238     238     409     409     409      409     
 
# > sweep500_tac2_coclust_siletti_clust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
#
# cluster freq
# 235    2
# 238    4
# 252    1
# 258    2
# 288    1
# 392    1
# 409   23
# 429    3
# <NA>    5


# By far the leading `cluster_id` is 409:
# 
# > t(subset(siletti_clust_terms, `Cluster ID` == "409"))
#
# Cluster ID                       "409"
# Cluster name                     "Splat_409"
# Supercluster                     "Splatter"
# Class auto-annotation            "NEUR"
# Neurotransmitter auto-annotation "GABA"
# Neuropeptide auto-annotation     "CALCB CBLN CHGA CHGB NAMPT NUCB NXPH SCG SST TAC TRH UBL proSAAS"
# Subtype auto-annotation          NA
# Transferred MTG Label            "N/A"
# Top three regions                "Basal forebrain: 66.5%, Amygdala: 14.2%, Thalamus: 9.9%"
# Top three dissections            "Human CaB: 17.8%, Human SI: 13.0%, Human NAC: 11.5%"
# Top Enriched Genes               "GLP1R, CRABP1, CER1, TAC3, AL096799.1, TRH, LINC02203, NANOS1, 
#                                   NPR3, CHRNA3"
# Number of cells                  "6395"
# DoubletFinder score              "0.0300207"
# Total UMI                        "14166.57"
# Fraction unspliced               "0.6586517"
# Fraction mitochondrial           "0.01125453"
# H19.30.002                       "1928"
# H19.30.001                       "2386"
# H18.30.002                       "2076"
# H18.30.001                       "5"
# Fraction cells from top donor    "0.373104"
# Number of donors                 "4" 

# Of the n=6395 cells in cluster 409, our subset contains:
# 
# > sum(siletti_neur_roi_sclust_0p8_dropAllMeis2$cluster_id == "409")
# [1] 2448

# By Roi:
# > table(subset(
#     siletti_neur_roi_sclust_0p8_dropAllMeis2@meta.data,
#     cluster_id == "409"
# )$roi)
# 
# Human CaB Human NAC  Human Pu 
# 1113       620       715 


# 500 features: Check top human subcluster co-cluster ---------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "human")) {
            subset(dsub, species == "human")$subcluster_id %>% 
                table %>% 
                (function(y) names(y)[which.max(y)]) %>% 
                as.numeric
        } else {
            NA
        }
    }) -> sweep500_tac2_coclust_siletti_subclust


# > sweep500_tac2_coclust_siletti_subclust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
# lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 166       166     1809    166     771     166      198     
# k25 433       NA      198     166     166     166      166     
# k30 402       402     166     166     166     166      166     
# k35 NA        NA      175     166     NA      NA       166     
# k40 1809      771     166     166     166     166      166     
# k45 1812      173     173     166     166     166      166   

# > sweep500_tac2_coclust_siletti_subclust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
#
# cluster freq
# 166   24
# 173    2
# 175    1
# 198    2
# 402    2
# 433    1
# 771    2
# 1809    2
# 1812    1
# <NA>    5

# By far the leading `subcluster_id` is 166:
# 
# > t(subset(siletti_subclust_terms, Subcluster == "166"))
# Subcluster                                             "166"
# Cluster                                                "409"
# Supercluster                                           "Splatter"
# Transferred MTG Label (Transferred from cluster level) NA
# Class                                                  "NEUR"
# Neurotransmitter                                       "NT-GABA NT-GABA"
# Neuropeptide                                           "NP-CALCB NP-CBLN NP-CHGA NP-CHGB 
#                                                         NP-NAMPT NP-NUCB NP-NUCB NP-NXPH 
#                                                         NP-NXPH NP-SCG NP-SST NP-TAC NP-TRH 
#                                                         NP-UBL NP-proSAAS"
# Top ROIGroupFine                                       "BasalForebrain: 67.0%, 
#                                                         Amygdala: 13.0%, 
#                                                         Thalamus: 11.5%"
# Top ROI                                                "CaB: 20.3%, Pu: 12.6%, NAC: 12.2%"
# Top enriched genes                                     "['GLP1R' 'CER1' 'CRABP1' 'AL096799.1' 
#                                                          'TAC3' 'LINC02203' 'SLC22A3' 'TRH' 
#                                                          'CHRNA3' 'NANOS1']"
# Number of cells                                        "4595"
# Number of donors                                       "4"
# DoubletFinderScore                                     "0.02443994"
# TotalUMI                                               "13009.59"
# unspliced_ratio                                        "0.6554633"
# MT_ratio                                               "0.01073824"
# H18.30.001                                             "2"
# H18.30.002                                             "1476"
# H19.30.001                                             "1708"
# H19.30.002                                             "1409"
# Fraction cells from top donor                          "1708"

# Of the n=4595 cells in subcluster 166, our subset contains:
#
# > sum(siletti_neur_roi_sclust_0p8_dropAllMeis2$subcluster_id == "166")
# [1] 1979

# By Roi:
# > table(subset(
#     siletti_neur_roi_sclust_0p8_dropAllMeis2@meta.data,
#     subcluster_id == "166"
# )$roi)
# 
# Human CaB Human NAC  Human Pu 
# 931       467       581 




# 500 features: Add is_Tac2 to parm_sweep metadata ------------------------

# Dividing up at subclass level:
liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$abc_subclass_Tac2 <- x$subclass
    x$abc_subclass_Tac2[x$subclass == "055 STR Lhx8 Gaba"] <- "055 STR Lhx8 Gaba: Tac2-"
    x$abc_subclass_Tac2[Cells(x) %in% tac2_barcodes] <- "055 STR Lhx8 Gaba: Tac2+"
    x
})

# Dividing up a cluster level:
liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$abc_cluster_Tac2 <- ifelse(
        x$species == "mouse",
        ifelse(
            Cells(x) %in% tac2_barcodes,
            paste0(x$cluster, ": Tac2+"),
            paste0(x$cluster, ": Tac2-")
        ),
        NA
    )
    x
})


# Add the manual cluster annotations (prepared near start of script)
liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$seurat_clusters_anno <- left_join(
        x@meta.data, 
        clust_anno, 
        by = c("RNA_snn_res.0.8", "species")
    )$seurat_clusters_anno
    x
})


# And divide the Th cluster up by Tac2 expression as well
# (first, verifying that all Tac2 barcodes are in the Th cluster)
liger_param_sweep$k20_lambda0.5 %>% 
    subset(species == "mouse") %>% 
    .[, tac2_barcodes] %>% 
    .$seurat_clusters_anno %>% 
    "=="("Th") %>% 
    all %>% 
    stopifnot

liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$abc_seurat_anno_Tac2 <- ifelse(
        x$species == "mouse",
        as.character(x$seurat_clusters_anno),
        NA
    )
    x$abc_seurat_anno_Tac2 <- ifelse(
        x$abc_seurat_anno_Tac2 == "Th",
        ifelse(
            Cells(x) %in% tac2_barcodes,
            "Th: Tac2+",
            "Th: Tac2-"
        ),
        x$abc_seurat_anno_Tac2
    )
    x
})


# Functions to assess cross-species co-clustering -------------------------

# 
# for each species' clustering (like clusters, subclusters, etc.), find,
# for each liger integration, the top integrated cluster into which those cells
# map
# 
# depends on `liger_param_sweep` in global environment
# 
get_top_liger_clusts <- function(species, species_grouping, new_nm_col) {
    spec <- species
    for_each_integration <- lapply(liger_param_sweep, {
        . %>% 
            .@meta.data %>% 
            subset(species == spec) %>% 
            .[, c("RNA_snn_res.0.5", species_grouping)] %>% 
            table %>% 
            as.matrix %>% # rows are int_cluster, columns are species clusters
            .[, colSums(.) > 0] %>% # ignore missing groupings
            (function(x) {
                sapply(
                    apply(x, 2, which.max), # for each species cluster, top int_cluster
                    function(i) rownames(x)[i] # and get the name (basically 0-index)
                )
            }) %>% 
            data.frame(
                species_grouping = names(.),
                top_liger_clust = .,
                row.names = NULL
            ) %>% 
            setNames(., c(new_nm_col, "top_liger_clust"))
    })
    
    # combine into a single dataframe (with a column for the integration)
    for_each_integration %>% 
        Map(
            function(x, nm) mutate(x, integration = nm),
            .,
            names(.)
        ) %>% 
        unname %>% 
        do.call(rbind, .)
}


# 
# for a given subject-query pair, find the top integrated cluster for the 
# subject. then, for query cells in that same integrated cluster, find out the
# most common cluster-of-origin. 
#
# depends on `liger_param_sweep` in global environment
# 
get_top_query_coclust <- function(
        species_1, 
        species_grouping_1, 
        species_2,
        species_grouping_2
) {
    nm_1 <- paste0(species_1, "_", species_grouping_1)
    nm_2 <- paste0(species_2, "_", species_grouping_2)
    
    df_coclust <- get_top_liger_clusts(
        species = species_1, 
        species_grouping = species_grouping_1,
        new_nm_col = nm_1
    )
    
    # for each `integration`+`top_liger_clust` found, get the most common
    # species_grouping_2
    df_coclust[[nm_2]] <- unlist(Map(
        function(int_clust, int) {
            dsub <- subset(
                liger_param_sweep[[int]]@meta.data, RNA_snn_res.0.5 == int_clust
            )
            if (any(dsub$species == species_2)) {
                subset(dsub, species == species_2)[[species_grouping_2]] %>% 
                    table %>% 
                    (function(y) names(y)[which.max(y)])
            } else {
                NA
            }
        },
        df_coclust$top_liger_clust,
        df_coclust$integration
    ))
    
    as.data.frame(table(df_coclust[, c(nm_1, nm_2)], useNA = "always"))
}


# 500 features: Assess cross-species co-clustering ------------------------

# running co-clustering analysis for all combinations of our categories of 
# interest
sweep500_all_clusts_top_query_coclust <- list(
    c("mouse-subclass"),
    c("mouse-supertype"),
    c("mouse-cluster"),
    c("mouse-abc_subclass_Tac2"),
    c("mouse-abc_cluster_Tac2"),
    c("mouse-seurat_clusters_anno"),
    c("mouse-abc_seurat_anno_Tac2"),
    c("human-cluster_id"),
    c("human-subcluster_id"),
    c("human-seurat_clusters_anno"),
    c("marmoset-CLUSTER"),
    c("marmoset-CLUSTER.SUBCLUSTER"),
    c("marmoset-seurat_clusters_anno")
) %>% 
    unlist %>% 
    expand.grid(., .) %>% 
    # remove any same-species comparisons
    .[-intersect(grep("mouse", .$Var1), grep("mouse", .$Var2)), ] %>%
    .[-intersect(grep("human", .$Var1), grep("human", .$Var2)), ] %>%
    .[-intersect(grep("marmoset", .$Var1), grep("marmoset", .$Var2)), ] %>% 
    unname %>% 
    as.list %>% 
    lapply(as.character) %>% 
    (function(x) Map(c, x[[1]], x[[2]])) %>% 
    setNames(., sapply(., function(x) paste0(x[1], ".", x[2]))) %>% 
    lapply(function(x) {
        get_top_query_coclust(
            species_1 = sub("-.*", "", x[1]), 
            species_grouping_1 = sub(".*-", "", x[1]), 
            species_2 = sub("-.*", "", x[2]),
            species_grouping_2 = sub(".*-", "", x[2])
        )
    })


# 500 features: Co-clustering confusion matrix plots ----------------------

path_sweep <- file.path(path_out, "plots_500features") # in case skipped section
if (!dir.exists(path_sweep))
    dir.create(path_sweep)

path_confusion <- file.path(path_sweep, "confusion_coclustering")
if (!dir.exists(path_confusion))
    dir.create(path_confusion)

sweep500_all_clusts_top_query_coclust %>% 
    lapply(function(x) {
        out_name <- paste0(names(x)[1], "_vs_", names(x)[2], ".pdf")
        
        suppressWarnings({
            x %>% 
                ggplot(
                    aes_string(names(x)[1], names(x)[2], fill = "Freq")
                ) +
                geom_tile() + 
                scale_fill_viridis_c(
                    guide = guide_colorbar(
                        draw.ulim = FALSE,
                        draw.llim = FALSE
                    )
                ) +
                coord_cartesian(
                    expand = FALSE,
                    clip = "off"
                ) +
                labs(
                    x = paste0("subject: ", names(x)[1]),
                    y = paste0("query: ", names(x)[2])
                ) + 
                theme_md_bw() + 
                theme(
                    axis.text.x = element_text(
                        angle = 90, hjust = 1, vjust = 0.5
                    )
                )
        })
        
        ggsave(
            file.path(path_confusion, out_name),
            width = 7 * 1.3,
            height = 7,
            limitsize = FALSE
        )
    }) %>% 
    invisible


# Remake selected plots with better optimized dimensions
path_confusion_select <- file.path(path_sweep, "confusion_coclustering_select")
if (!dir.exists(path_confusion_select))
    dir.create(path_confusion_select)

sweep500_all_clusts_top_query_coclust[c(
    "human-seurat_clusters_anno.marmoset-seurat_clusters_anno",
    "mouse-abc_seurat_anno_Tac2.human-seurat_clusters_anno",
    "mouse-abc_seurat_anno_Tac2.marmoset-seurat_clusters_anno"
)] %>% 
    
    lapply(function(x) {
        out_name <- paste0(names(x)[1], "_vs_", names(x)[2], ".pdf")
        x %>% 
            ggplot(
                aes_string(names(x)[1], names(x)[2], fill = "Freq")
            ) +
            geom_tile() + 
            scale_fill_viridis_c(
                guide = guide_colorbar(
                    # draw.ulim = FALSE,
                    draw.llim = FALSE
                )
            ) +
            coord_cartesian(
                expand = FALSE,
                clip = "off"
            ) +
            labs(
                x = paste0("subject: ", sub("_.*", "", names(x)[1])),
                y = paste0("query: ", sub("_.*", "", names(x)[2]))
            ) + 
            theme_md_bw() + 
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            )
        ggsave(
            file.path(path_confusion_select, out_name),
            width = 5, 
            height = 4,
            limitsize = FALSE
        )
    }) %>% 
    invisible


# 1000 features: Parameter sweep ------------------------------------------

rm(liger_param_sweep)
invisible(gc())

seur1000 <- seurat_merge %>% 
    FindVariableFeatures(nfeatures = 1000) %>% 
    ScaleData(split.by = "species", do.center = FALSE)

# n=42 param combos
liger_param_sweep <- expand.grid(
    k = seq(20, 45, 5),
    lambda = 2^seq(-1, 5, 1) 
) %>% 
    as.list %>% 
    (function(x) Map(c, x$k, x$lambda)) %>% 
    (function(x) setNames(
        x, sapply(x, function(x) paste0("k", x[1], "_lambda", x[2]))
    )) %>% 
    lapply(function(x) {
        seur <- seur1000
        seur@reductions <- list()
        seur <- RunOptimizeALS(
            seur, k = x[1], lambda = x[2], split.by = "species"
        )
        seur <- RunQuantileNorm(seur, split.by = "species")
        seur <- FindNeighbors(
            seur, reduction = "iNMF", dims = 1:20, verbose = FALSE
        )
        seur <- RunUMAP(
            seur, dims = 1:ncol(seur[["iNMF"]]), reduction = "iNMF", 
            verbose = FALSE
        )
        seur
    })
saveRDS(
    liger_param_sweep, 
    file.path(path_out, "liger_param_sweep_1000features.rds")
)
rm(seur1000)


# 1000 features: Plots ----------------------------------------------------

path_sweep <- file.path(path_out, "plots_1000features")

if (!dir.exists(path_sweep))
    dir.create(path_sweep)

invisible(Map(
    function(x, nm) {
        
        # --- plot species on UMAP
        path_umap_species <- file.path(path_sweep, "umap_species")
        if (!dir.exists(path_umap_species))
            dir.create(path_umap_species)
        
        DimPlot(x, group.by = "species")
        ggsave(
            file.path(path_umap_species, paste0(nm, ".pdf")),
            width = 160, height = 120, units = "mm"
        )
        
        # --- umap plot some markers together in a single file
        path_umap_combo <- file.path(path_sweep, "umap_3species_markers")
        if (!dir.exists(path_umap_combo))
            dir.create(path_umap_combo)
        
        p1 <- FeaturePlot(
            x,
            c("LHX6", "TAC3", "CHRNA3", "TH", "CRABP1"),
            split.by = "species",
            order = TRUE
        ) 
        p2 <- FeaturePlot(
            x,
            c("SST", "CHAT", "PVALB", "PTHLH", "VIP"),
            split.by = "species",
            order = TRUE
        )
        {p1 | p2} & 
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                strip.text.x = element_text(angle = 0)
            )
        ggsave(
            file.path(path_umap_combo, paste0(nm, ".pdf")),
            width = 12, height = 10, units = "in"
        )
        
        # --- plot separate umaps for all marker genes
        path_umap_markers <- file.path(path_sweep, "umap_markers")
        if (!dir.exists(path_umap_markers))
            dir.create(path_umap_markers)
        
        path_umap_param <- file.path(path_umap_markers, nm)
        if (!dir.exists(path_umap_param))
            dir.create(path_umap_param)
        
        for (gn in marker_genes) {
            FeaturePlot(x, gn, split.by = "species", order = TRUE)
            ggsave(
                file.path(path_umap_param, paste0(gn, ".pdf")),
                width = 270, height = 90, units = "mm"
            )
        }
    },
    liger_param_sweep,
    names(liger_param_sweep)
))


# 1000 features: Cluster integrations -------------------------------------

# 
# note the previous clusters are present in `RNA_snn_res.0.8`;
# `seurat_clusters` is always overwritten so it's best not to use that
# 

liger_param_sweep <- lapply(
    liger_param_sweep, 
    FindClusters,
    resolution = 0.5,
    verbose = FALSE
)


# 1000 features: Check number of clusters ---------------------------------

# > liger_param_sweep %>% 
#     lapply(. %>% .$RNA_snn_res.0.5 %>% nlevels) %>% 
#     # make a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 23        22      20      25      24      21       22      
# k25 24        23      29      21      21      19       22      
# k30 20        25      21      18      20      19       21      
# k35 25        26      26      26      26      23       25      
# k40 24        24      22      19      21      20       21      
# k45 24        22      27      23      23      24       26     


# 1000 features: Check co-clustering of Tac2 cells with themselves --------

liger_param_sweep %>% 
    lapply(function(x) x[, Cells(x) %in% tac2_barcodes]) %>% 
    lapply(. %>% .$RNA_snn_res.0.5 %>% table %>% max) %>% 
    # make a matrix
    split(., sub("_.*", "", names(.))) %>% 
    lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
    lapply(unlist) %>% 
    do.call(rbind, .) ->
    sweep1000_tac2_self_coclust

# > sweep1000_tac2_self_coclust
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20        15      47      37      19      31       25       32
# k25        24      25      19      24      21       35       20
# k30        16      42      30      35      32       37       35
# k35        12      40      50      39      31       44       46
# k40        27      23      23      13      23       16       27
# k45        14      34      48      25      16       24       36

# As a fraction of Tac2+ cells:
# 
# > sweep1000_tac2_self_coclust %>% 
#     "/"(length(tac2_barcodes)) %>% 
#     round(2)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20      0.29    0.92    0.73    0.37    0.61     0.49     0.63
# k25      0.47    0.49    0.37    0.47    0.41     0.69     0.39
# k30      0.31    0.82    0.59    0.69    0.63     0.73     0.69
# k35      0.24    0.78    0.98    0.76    0.61     0.86     0.90
# k40      0.53    0.45    0.45    0.25    0.45     0.31     0.53
# k45      0.27    0.67    0.94    0.49    0.31     0.47     0.71


# 1000 features: Check top marm co-cluster --------------------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "marmoset")) {
            subset(dsub, species == "marmoset")$CLUSTER %>% 
                table %>% 
                (function(y) names(y)[which.max(y)]) %>% 
                as.numeric
        } else {
            NA
        }
    }) -> sweep1000_tac2_coclust_marm_clust

# > sweep1000_tac2_coclust_marm_clust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 4         7       7       7       7       7        7       
# k25 NA        8       7       4       7       7        7       
# k30 8         NA      NA      9       7       7        7       
# k35 8         5       NA      NA      3       NA       NA      
# k40 4         7       9       9       9       8        4       
# k45 8         7       NA      NA      7       7        NA  

# > sweep1000_tac2_coclust_marm_clust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
#
# cluster freq
#       3    1
#       4    4
#       5    1
#       7   17
#       8    5
#       9    4
#    <NA>   10



# 1000 features: Check top marm subcluster co-cluster ---------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "marmoset")) {
            subset(dsub, species == "marmoset")$CLUSTER.SUBCLUSTER %>% 
                table %>% 
                (function(y) names(y)[which.max(y)])
        } else {
            NA
        }
    }) -> sweep1000_tac2_coclust_marm_subclust


# > sweep1000_tac2_coclust_marm_subclust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
# 
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 "04-00"   "07-00" "07-00" "07-00" "07-00" "07-00"  "07-00" 
# k25 NA        "08-01" "07-00" "04-00" "07-00" "07-00"  "07-00" 
# k30 "08-02"   NA      NA      "09-00" "07-00" "07-00"  "07-00" 
# k35 "08-00"   "05-00" NA      NA      "03-00" NA       NA      
# k40 "04-00"   "07-00" "09-00" "09-00" "09-00" "09-00"  "04-00" 
# k45 "07-00"   "07-00" NA      NA      "07-00" "07-00"  NA      

# > sweep1000_tac2_coclust_marm_subclust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
# 
# cluster freq
#   03-00    1
#   04-00    4
#   05-00    1
#   07-00   18
#   08-00    1
#   08-01    1
#   08-02    1
#   09-00    5
#    <NA>   10



# 1000 features: Check top human co-cluster -------------------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "human")) {
            subset(dsub, species == "human")$cluster_id %>% 
                table %>% 
                (function(y) names(y)[which.max(y)]) %>% 
                as.numeric
        } else {
            NA
        }
    }) -> sweep1000_tac2_coclust_siletti_clust


# > sweep1000_tac2_coclust_siletti_clust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 409       409     392     409     NA      NA       NA      
# k25 NA        409     409     409     409     409      409     
# k30 409       NA      423     409     409     409      409     
# k35 409       409     NA      NA      428     NA       NA      
# k40 409       409     235     235     235     409      409     
# k45 238       409     NA      NA      409     409      NA    

# > sweep1000_tac2_coclust_siletti_clust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
#
# cluster freq
#     235    3
#     238    1
#     392    1
#     409   23
#     423    1
#     428    1
#    <NA>   12


# 1000 features: Check top human subcluster co-cluster --------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "human")) {
            subset(dsub, species == "human")$subcluster_id %>% 
                table %>% 
                (function(y) names(y)[which.max(y)]) %>% 
                as.numeric
        } else {
            NA
        }
    }) -> sweep1000_tac2_coclust_siletti_subclust

# > sweep1000_tac2_coclust_siletti_subclust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 166       165     771     167     NA      NA       NA      
# k25 NA        166     166     166     166     166      166     
# k30 166       NA      777     166     166     166      166     
# k35 166       165     NA      NA      1730    NA       NA      
# k40 166       166     193     193     193     166      166     
# k45 168       166     NA      NA      166     166      NA      

# > sweep1000_tac2_coclust_siletti_subclust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
#
# cluster freq
#     165    2
#     166   20
#     167    1
#     168    1
#     193    3
#     771    1
#     777    1
#    1730    1
#    <NA>   12


# 1000 features: Add is_Tac2 to parm_sweep metadata -----------------------

# Dividing up at subclass level:
liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$abc_subclass_Tac2 <- x$subclass
    x$abc_subclass_Tac2[x$subclass == "055 STR Lhx8 Gaba"] <- "055 STR Lhx8 Gaba: Tac2-"
    x$abc_subclass_Tac2[Cells(x) %in% tac2_barcodes] <- "055 STR Lhx8 Gaba: Tac2+"
    x
})

# Dividing up a cluster level:
liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$abc_cluster_Tac2 <- ifelse(
        x$species == "mouse",
        ifelse(
            Cells(x) %in% tac2_barcodes,
            paste0(x$cluster, ": Tac2+"),
            paste0(x$cluster, ": Tac2-")
        ),
        NA
    )
    x
})


# Add the manual cluster annotations (prepared near start of script)
liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$seurat_clusters_anno <- left_join(
        x@meta.data, 
        clust_anno, 
        by = c("RNA_snn_res.0.8", "species")
    )$seurat_clusters_anno
    x
})


# And divide the Th cluster up by Tac2 expression as well
# (first, verifying that all Tac2 barcodes are in the Th cluster)
liger_param_sweep$k20_lambda0.5 %>% 
    subset(species == "mouse") %>% 
    .[, tac2_barcodes] %>% 
    .$seurat_clusters_anno %>% 
    "=="("Th") %>% 
    all %>% 
    stopifnot

liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$abc_seurat_anno_Tac2 <- ifelse(
        x$species == "mouse",
        as.character(x$seurat_clusters_anno),
        NA
    )
    x$abc_seurat_anno_Tac2 <- ifelse(
        x$abc_seurat_anno_Tac2 == "Th",
        ifelse(
            Cells(x) %in% tac2_barcodes,
            "Th: Tac2+",
            "Th: Tac2-"
        ),
        x$abc_seurat_anno_Tac2
    )
    x
})


# 1000 features: Assess cross-species co-clustering -----------------------

# running co-clustering analysis for all combinations of our categories of 
# interest
sweep1000_all_clusts_top_query_coclust <- list(
    c("mouse-subclass"),
    c("mouse-supertype"),
    c("mouse-cluster"),
    c("mouse-abc_subclass_Tac2"),
    c("mouse-abc_cluster_Tac2"),
    c("mouse-seurat_clusters_anno"),
    c("mouse-abc_seurat_anno_Tac2"),
    c("human-cluster_id"),
    c("human-subcluster_id"),
    c("human-seurat_clusters_anno"),
    c("marmoset-CLUSTER"),
    c("marmoset-CLUSTER.SUBCLUSTER"),
    c("marmoset-seurat_clusters_anno")
) %>% 
    unlist %>% 
    expand.grid(., .) %>% 
    # remove any same-species comparisons
    .[-intersect(grep("mouse", .$Var1), grep("mouse", .$Var2)), ] %>%
    .[-intersect(grep("human", .$Var1), grep("human", .$Var2)), ] %>%
    .[-intersect(grep("marmoset", .$Var1), grep("marmoset", .$Var2)), ] %>% 
    unname %>% 
    as.list %>% 
    lapply(as.character) %>% 
    (function(x) Map(c, x[[1]], x[[2]])) %>% 
    setNames(., sapply(., function(x) paste0(x[1], ".", x[2]))) %>% 
    lapply(function(x) {
        get_top_query_coclust(
            species_1 = sub("-.*", "", x[1]), 
            species_grouping_1 = sub(".*-", "", x[1]), 
            species_2 = sub("-.*", "", x[2]),
            species_grouping_2 = sub(".*-", "", x[2])
        )
    })


# 1000 features: Co-clustering confusion matrix plots ---------------------

path_sweep <- file.path(path_out, "plots_1000features") # in case skipped section
if (!dir.exists(path_sweep))
    dir.create(path_sweep)

path_confusion <- file.path(path_sweep, "confusion_coclustering")
if (!dir.exists(path_confusion))
    dir.create(path_confusion)

sweep1000_all_clusts_top_query_coclust %>% 
    lapply(function(x) {
        out_name <- paste0(names(x)[1], "_vs_", names(x)[2], ".pdf")
        
        suppressWarnings({
            x %>% 
                ggplot(
                    aes_string(names(x)[1], names(x)[2], fill = "Freq")
                ) +
                geom_tile() + 
                scale_fill_viridis_c(
                    guide = guide_colorbar(
                        draw.ulim = FALSE,
                        draw.llim = FALSE
                    )
                ) +
                coord_cartesian(
                    expand = FALSE,
                    clip = "off"
                ) +
                labs(
                    x = paste0("subject: ", names(x)[1]),
                    y = paste0("query: ", names(x)[2])
                ) + 
                theme_md_bw() + 
                theme(
                    axis.text.x = element_text(
                        angle = 90, hjust = 1, vjust = 0.5
                    )
                )
        })
        
        ggsave(
            file.path(path_confusion, out_name),
            width = 7 * 1.3,
            height = 7,
            limitsize = FALSE
        )
    }) %>% 
    invisible


# Remake selected plots with better optimized dimensions
path_confusion_select <- file.path(path_sweep, "confusion_coclustering_select")
if (!dir.exists(path_confusion_select))
    dir.create(path_confusion_select)

sweep1000_all_clusts_top_query_coclust[c(
    "human-seurat_clusters_anno.marmoset-seurat_clusters_anno",
    "mouse-abc_seurat_anno_Tac2.human-seurat_clusters_anno",
    "mouse-abc_seurat_anno_Tac2.marmoset-seurat_clusters_anno"
)] %>% 
    
    lapply(function(x) {
        out_name <- paste0(names(x)[1], "_vs_", names(x)[2], ".pdf")
        x %>% 
            ggplot(
                aes_string(names(x)[1], names(x)[2], fill = "Freq")
            ) +
            geom_tile() + 
            scale_fill_viridis_c(
                guide = guide_colorbar(
                    # draw.ulim = FALSE,
                    draw.llim = FALSE
                )
            ) +
            coord_cartesian(
                expand = FALSE,
                clip = "off"
            ) +
            labs(
                x = paste0("subject: ", sub("_.*", "", names(x)[1])),
                y = paste0("query: ", sub("_.*", "", names(x)[2]))
            ) + 
            theme_md_bw() + 
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            )
        ggsave(
            file.path(path_confusion_select, out_name),
            width = 5, 
            height = 4,
            limitsize = FALSE
        )
    }) %>% 
    invisible


# 2000 features: Parameter sweep ------------------------------------------

rm(liger_param_sweep)
invisible(gc())

seur2000 <- seurat_merge %>% 
    FindVariableFeatures(nfeatures = 2000) %>% 
    ScaleData(split.by = "species", do.center = FALSE)

# n=42 param combos
liger_param_sweep <- expand.grid(
    k = seq(20, 45, 5),
    lambda = 2^seq(-1, 5, 1) 
) %>% 
    as.list %>% 
    (function(x) Map(c, x$k, x$lambda)) %>% 
    (function(x) setNames(
        x, sapply(x, function(x) paste0("k", x[1], "_lambda", x[2]))
    )) %>% 
    lapply(function(x) {
        seur <- seur2000
        seur@reductions <- list()
        seur <- RunOptimizeALS(
            seur, k = x[1], lambda = x[2], split.by = "species"
        )
        seur <- RunQuantileNorm(seur, split.by = "species")
        seur <- FindNeighbors(
            seur, reduction = "iNMF", dims = 1:20, verbose = FALSE
        )
        seur <- RunUMAP(
            seur, dims = 1:ncol(seur[["iNMF"]]), reduction = "iNMF", 
            verbose = FALSE
        )
        seur
    })
saveRDS(
    liger_param_sweep, 
    file.path(path_out, "liger_param_sweep_2000features.rds")
)
rm(seur2000)


# 2000 features: Plots ----------------------------------------------------

path_sweep <- file.path(path_out, "plots_2000features")

if (!dir.exists(path_sweep))
    dir.create(path_sweep)

invisible(Map(
    function(x, nm) {
        
        # --- plot species on UMAP
        path_umap_species <- file.path(path_sweep, "umap_species")
        if (!dir.exists(path_umap_species))
            dir.create(path_umap_species)
        
        DimPlot(x, group.by = "species")
        ggsave(
            file.path(path_umap_species, paste0(nm, ".pdf")),
            width = 160, height = 120, units = "mm"
        )
        
        # --- umap plot some markers together in a single file
        path_umap_combo <- file.path(path_sweep, "umap_3species_markers")
        if (!dir.exists(path_umap_combo))
            dir.create(path_umap_combo)
        
        p1 <- FeaturePlot(
            x,
            c("LHX6", "TAC3", "CHRNA3", "TH", "CRABP1"),
            split.by = "species",
            order = TRUE
        ) 
        p2 <- FeaturePlot(
            x,
            c("SST", "CHAT", "PVALB", "PTHLH", "VIP"),
            split.by = "species",
            order = TRUE
        )
        {p1 | p2} & 
            theme(
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                strip.text.x = element_text(angle = 0)
            )
        ggsave(
            file.path(path_umap_combo, paste0(nm, ".pdf")),
            width = 12, height = 10, units = "in"
        )
        
        # --- plot separate umaps for all marker genes
        path_umap_markers <- file.path(path_sweep, "umap_markers")
        if (!dir.exists(path_umap_markers))
            dir.create(path_umap_markers)
        
        path_umap_param <- file.path(path_umap_markers, nm)
        if (!dir.exists(path_umap_param))
            dir.create(path_umap_param)
        
        for (gn in marker_genes) {
            FeaturePlot(x, gn, split.by = "species", order = TRUE)
            ggsave(
                file.path(path_umap_param, paste0(gn, ".pdf")),
                width = 270, height = 90, units = "mm"
            )
        }
    },
    liger_param_sweep,
    names(liger_param_sweep)
))


# 2000 features: Cluster integrations -------------------------------------

# 
# note the previous clusters are present in `RNA_snn_res.0.8`;
# `seurat_clusters` is always overwritten so it's best not to use that
# 

liger_param_sweep <- lapply(
    liger_param_sweep, 
    FindClusters,
    resolution = 0.5,
    verbose = FALSE
)


# 2000 features: Check number of clusters ---------------------------------

# > liger_param_sweep %>% 
#     lapply(. %>% .$RNA_snn_res.0.5 %>% nlevels) %>% 
#     # make a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 23        22      25      20      23      21       25      
# k25 30        26      23      26      23      23       20      
# k30 32        26      23      22      21      22       21      
# k35 26        27      27      21      24      21       17      
# k40 27        23      18      21      20      19       19      
# k45 27        30      24      23      25      20       23    


# 2000 features: Check co-clustering of Tac2 cells with themselves --------

liger_param_sweep %>% 
    lapply(function(x) x[, Cells(x) %in% tac2_barcodes]) %>% 
    lapply(. %>% .$RNA_snn_res.0.5 %>% table %>% max) %>% 
    # make a matrix
    split(., sub("_.*", "", names(.))) %>% 
    lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
    lapply(unlist) %>% 
    do.call(rbind, .) ->
    sweep2000_tac2_self_coclust

# > sweep2000_tac2_self_coclust
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20        46      28      44      18      44       16       27
# k25        28      41      38      21      36       35       37
# k30        32      45      33      17      20       33       21
# k35        32      48      24      32      46       38       50
# k40        13      29      40      42      20       22       17
# k45        36      31      46      47      48       48       34

# As a fraction of Tac2+ cells:
# 
# > sweep2000_tac2_self_coclust %>% 
#     "/"(length(tac2_barcodes)) %>% 
#     round(2)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20      0.90    0.55    0.86    0.35    0.86     0.31     0.53
# k25      0.55    0.80    0.75    0.41    0.71     0.69     0.73
# k30      0.63    0.88    0.65    0.33    0.39     0.65     0.41
# k35      0.63    0.94    0.47    0.63    0.90     0.75     0.98
# k40      0.25    0.57    0.78    0.82    0.39     0.43     0.33
# k45      0.71    0.61    0.90    0.92    0.94     0.94     0.67


# 2000 features: Check top marm co-cluster --------------------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "marmoset")) {
            subset(dsub, species == "marmoset")$CLUSTER %>% 
                table %>% 
                (function(y) names(y)[which.max(y)]) %>% 
                as.numeric
        } else {
            NA
        }
    }) -> sweep2000_tac2_coclust_marm_clust

# > sweep2000_tac2_coclust_marm_clust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 NA        4       8       7       7       4        4       
# k25 NA        8       NA      NA      NA      NA       NA      
# k30 NA        3       4       7       7       7        7       
# k35 NA        NA      4       7       NA      7        NA      
# k40 8         7       4       8       8       4        4       
# k45 9         NA      NA      NA      7       7        NA 

# > sweep2000_tac2_coclust_marm_clust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
#
# cluster freq
#       3    1
#       4    8
#       7   11
#       8    5
#       9    1
#    <NA>   16



# 2000 features: Check top marm subcluster co-cluster ---------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "marmoset")) {
            subset(dsub, species == "marmoset")$CLUSTER.SUBCLUSTER %>% 
                table %>% 
                (function(y) names(y)[which.max(y)])
        } else {
            NA
        }
    }) -> sweep2000_tac2_coclust_marm_subclust


# > sweep2000_tac2_coclust_marm_subclust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
# 
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 NA        "04-00" "08-01" "07-00" "07-00" "04-00"  "04-00" 
# k25 NA        "08-00" NA      NA      NA      NA       NA      
# k30 NA        "03-00" "04-00" "07-00" "07-00" "07-00"  "07-00" 
# k35 NA        NA      "04-00" "07-00" NA      "07-00"  NA      
# k40 "08-00"   "07-00" "04-00" "08-02" "08-01" "04-00"  "04-00" 
# k45 "09-00"   NA      NA      NA      "07-00" "07-00"  NA     

# > sweep2000_tac2_coclust_marm_subclust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
# 
# cluster freq
#   03-00    1
#   04-00    8
#   07-00   11
#   08-00    2
#   08-01    2
#   08-02    1
#   09-00    1
#    <NA>   16


# 2000 features: Check top human co-cluster -------------------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "human")) {
            subset(dsub, species == "human")$cluster_id %>% 
                table %>% 
                (function(y) names(y)[which.max(y)]) %>% 
                as.numeric
        } else {
            NA
        }
    }) -> sweep2000_tac2_coclust_siletti_clust


# > sweep2000_tac2_coclust_siletti_clust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 NA        392     290     409     429     392      NA      
# k25 NA        235     235     235     NA      NA       NA      
# k30 NA        378     409     409     409     409      409     
# k35 NA        NA      392     423     NA      409      NA      
# k40 409       409     238     235     238     238      238     
# k45 409       NA      NA      NA      409     409      409  

# > sweep2000_tac2_coclust_siletti_clust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
#
# cluster freq
#     235    4
#     238    4
#     290    1
#     378    1
#     392    3
#     409   13
#     423    1
#     429    1
#    <NA>   14


# 2000 features: Check top human subcluster co-cluster --------------------

liger_param_sweep %>% 
    lapply(function(x) {
        top_tac2_clust <- x[, Cells(x) %in% tac2_barcodes] %>% 
            .$RNA_snn_res.0.5 %>% 
            table %>% 
            (function(y) names(y)[which.max(y)])
        
        dsub <- subset(x, RNA_snn_res.0.5 == top_tac2_clust)
        
        if (any(dsub$species == "human")) {
            subset(dsub, species == "human")$subcluster_id %>% 
                table %>% 
                (function(y) names(y)[which.max(y)]) %>% 
                as.numeric
        } else {
            NA
        }
    }) -> sweep2000_tac2_coclust_siletti_subclust

# > sweep2000_tac2_coclust_siletti_subclust %>% 
#     # make into a matrix
#     split(., sub("_.*", "", names(.))) %>% 
#     lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
#     do.call(rbind, .)
#
#     lambda0.5 lambda1 lambda2 lambda4 lambda8 lambda16 lambda32
# k20 NA        777     555     165     1809    1502     NA      
# k25 NA        198     198     198     NA      NA       NA      
# k30 NA        854     166     166     166     166      166     
# k35 NA        NA      777     777     NA      166      NA      
# k40 166       166     173     193     173     173      173     
# k45 166       NA      NA      NA      166     166      166       

# > sweep2000_tac2_coclust_siletti_subclust %>% 
#     unlist %>% 
#     table(useNA = "always") %>% 
#     as.data.frame %>% 
#     setNames(c("cluster", "freq")) %>% 
#     print.data.frame(row.names = FALSE)
#
# cluster freq
#     165    1
#     166   12
#     173    4
#     193    1
#     198    3
#     555    1
#     777    3
#     854    1
#    1502    1
#    1809    1
#    <NA>   14


# 2000 features: Add is_Tac2 to parm_sweep metadata -----------------------

# Dividing up at subclass level:
liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$abc_subclass_Tac2 <- x$subclass
    x$abc_subclass_Tac2[x$subclass == "055 STR Lhx8 Gaba"] <- "055 STR Lhx8 Gaba: Tac2-"
    x$abc_subclass_Tac2[Cells(x) %in% tac2_barcodes] <- "055 STR Lhx8 Gaba: Tac2+"
    x
})

# Dividing up a cluster level:
liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$abc_cluster_Tac2 <- ifelse(
        x$species == "mouse",
        ifelse(
            Cells(x) %in% tac2_barcodes,
            paste0(x$cluster, ": Tac2+"),
            paste0(x$cluster, ": Tac2-")
        ),
        NA
    )
    x
})


# Add the manual cluster annotations (prepared near start of script)
liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$seurat_clusters_anno <- left_join(
        x@meta.data, 
        clust_anno, 
        by = c("RNA_snn_res.0.8", "species")
    )$seurat_clusters_anno
    x
})


# And divide the Th cluster up by Tac2 expression as well
# (first, verifying that all Tac2 barcodes are in the Th cluster)
liger_param_sweep$k20_lambda0.5 %>% 
    subset(species == "mouse") %>% 
    .[, tac2_barcodes] %>% 
    .$seurat_clusters_anno %>% 
    "=="("Th") %>% 
    all %>% 
    stopifnot

liger_param_sweep <- lapply(liger_param_sweep, function(x) {
    x$abc_seurat_anno_Tac2 <- ifelse(
        x$species == "mouse",
        as.character(x$seurat_clusters_anno),
        NA
    )
    x$abc_seurat_anno_Tac2 <- ifelse(
        x$abc_seurat_anno_Tac2 == "Th",
        ifelse(
            Cells(x) %in% tac2_barcodes,
            "Th: Tac2+",
            "Th: Tac2-"
        ),
        x$abc_seurat_anno_Tac2
    )
    x
})


# 2000 features: Assess cross-species co-clustering -----------------------

# running co-clustering analysis for all combinations of our categories of 
# interest
sweep2000_all_clusts_top_query_coclust <- list(
    c("mouse-subclass"),
    c("mouse-supertype"),
    c("mouse-cluster"),
    c("mouse-abc_subclass_Tac2"),
    c("mouse-abc_cluster_Tac2"),
    c("mouse-seurat_clusters_anno"),
    c("mouse-abc_seurat_anno_Tac2"),
    c("human-cluster_id"),
    c("human-subcluster_id"),
    c("human-seurat_clusters_anno"),
    c("marmoset-CLUSTER"),
    c("marmoset-CLUSTER.SUBCLUSTER"),
    c("marmoset-seurat_clusters_anno")
) %>% 
    unlist %>% 
    expand.grid(., .) %>% 
    # remove any same-species comparisons
    .[-intersect(grep("mouse", .$Var1), grep("mouse", .$Var2)), ] %>%
    .[-intersect(grep("human", .$Var1), grep("human", .$Var2)), ] %>%
    .[-intersect(grep("marmoset", .$Var1), grep("marmoset", .$Var2)), ] %>% 
    unname %>% 
    as.list %>% 
    lapply(as.character) %>% 
    (function(x) Map(c, x[[1]], x[[2]])) %>% 
    setNames(., sapply(., function(x) paste0(x[1], ".", x[2]))) %>% 
    lapply(function(x) {
        get_top_query_coclust(
            species_1 = sub("-.*", "", x[1]), 
            species_grouping_1 = sub(".*-", "", x[1]), 
            species_2 = sub("-.*", "", x[2]),
            species_grouping_2 = sub(".*-", "", x[2])
        )
    })


# 2000 features: Co-clustering confusion matrix plots ---------------------

path_sweep <- file.path(path_out, "plots_2000features") # in case skipped section
if (!dir.exists(path_sweep))
    dir.create(path_sweep)

path_confusion <- file.path(path_sweep, "confusion_coclustering")
if (!dir.exists(path_confusion))
    dir.create(path_confusion)

sweep2000_all_clusts_top_query_coclust %>% 
    lapply(function(x) {
        out_name <- paste0(names(x)[1], "_vs_", names(x)[2], ".pdf")
        
        suppressWarnings({
            x %>% 
                ggplot(
                    aes_string(names(x)[1], names(x)[2], fill = "Freq")
                ) +
                geom_tile() + 
                scale_fill_viridis_c(
                    guide = guide_colorbar(
                        draw.ulim = FALSE,
                        draw.llim = FALSE
                    )
                ) +
                coord_cartesian(
                    expand = FALSE,
                    clip = "off"
                ) +
                labs(
                    x = paste0("subject: ", names(x)[1]),
                    y = paste0("query: ", names(x)[2])
                ) + 
                theme_md_bw() + 
                theme(
                    axis.text.x = element_text(
                        angle = 90, hjust = 1, vjust = 0.5
                    )
                )
        })
        
        ggsave(
            file.path(path_confusion, out_name),
            width = 7 * 1.3,
            height = 7,
            limitsize = FALSE
        )
    }) %>% 
    invisible


# Remake selected plots with better optimized dimensions
path_confusion_select <- file.path(path_sweep, "confusion_coclustering_select")
if (!dir.exists(path_confusion_select))
    dir.create(path_confusion_select)

sweep2000_all_clusts_top_query_coclust[c(
    "human-seurat_clusters_anno.marmoset-seurat_clusters_anno",
    "mouse-abc_seurat_anno_Tac2.human-seurat_clusters_anno",
    "mouse-abc_seurat_anno_Tac2.marmoset-seurat_clusters_anno"
)] %>% 
    
    lapply(function(x) {
        out_name <- paste0(names(x)[1], "_vs_", names(x)[2], ".pdf")
        x %>% 
            ggplot(
                aes_string(names(x)[1], names(x)[2], fill = "Freq")
            ) +
            geom_tile() + 
            scale_fill_viridis_c(
                guide = guide_colorbar(
                    # draw.ulim = FALSE,
                    draw.llim = FALSE
                )
            ) +
            coord_cartesian(
                expand = FALSE,
                clip = "off"
            ) +
            labs(
                x = paste0("subject: ", sub("_.*", "", names(x)[1])),
                y = paste0("query: ", sub("_.*", "", names(x)[2]))
            ) + 
            theme_md_bw() + 
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            )
        ggsave(
            file.path(path_confusion_select, out_name),
            width = 5, 
            height = 4,
            limitsize = FALSE
        )
    }) %>% 
    invisible


# Aggregate summaries of 500, 1000, 2000 features -------------------------

rm(liger_param_sweep)
invisible(gc())

# max number of abc tac2 cells in a single cluster
tac2_coclust <- list(
    f500  = sweep500_tac2_self_coclust,
    f1000 = sweep1000_tac2_self_coclust,
    f2000 = sweep2000_tac2_self_coclust
)

# most frequent marmoset cluster with which abc tac2 cells co-cluster
tac2_marm_clust_coclust <- list(
    f500 = sweep500_tac2_coclust_marm_clust,
    f1000 = sweep1000_tac2_coclust_marm_clust,
    f2000 = sweep2000_tac2_coclust_marm_clust
) %>% 
    # make into matrices
    lapply({
        . %>% 
            split(., sub("_.*", "", names(.))) %>% 
            lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
            do.call(rbind, .)
    })

# most frequent marmoset subcluster with which abc tac2 cells co-cluster
tac2_marm_subclust_coclust <- list(
    f500 = sweep500_tac2_coclust_marm_subclust,
    f1000 = sweep1000_tac2_coclust_marm_subclust,
    f2000 = sweep2000_tac2_coclust_marm_subclust
) %>% 
    # make into matrices
    lapply({
        . %>% 
            split(., sub("_.*", "", names(.))) %>% 
            lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
            do.call(rbind, .)
    })

# most frequent siletti cluster with which abc tac2 cells co-cluster
tac2_siletti_clust_coclust <- list(
    f500 = sweep500_tac2_coclust_siletti_clust,
    f1000 = sweep1000_tac2_coclust_siletti_clust,
    f2000 = sweep2000_tac2_coclust_siletti_clust
) %>% 
    # make into matrices
    lapply({
        . %>% 
            split(., sub("_.*", "", names(.))) %>% 
            lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
            do.call(rbind, .)
    })

# most frequent siletti subcluster with which abc tac2 cells co-cluster
tac2_siletti_subclust_coclust <- list(
    f500 = sweep500_tac2_coclust_siletti_subclust,
    f1000 = sweep1000_tac2_coclust_siletti_subclust,
    f2000 = sweep2000_tac2_coclust_siletti_subclust
) %>% 
    # make into matrices
    lapply({
        . %>% 
            split(., sub("_.*", "", names(.))) %>% 
            lapply(function(x) setNames(x, sub(".*_", "", names(x)))) %>% 
            do.call(rbind, .)
    })


# ---
# Make into dataframes
# ---

df_tac2_coclust <- reshape2::melt(tac2_coclust) %>% 
    setNames(c("k", "lambda", "tac2_cells_in_clust", "features")) %>% 
    mutate(
        k = sub("k", "", k),
        lambda = sub("lambda", "", lambda),
        features = sub("f", "", features)
    ) %>% 
    .[, c(1, 2, 4, 3)]

df_tac2_marm_clust_coclust <- reshape2::melt(tac2_marm_clust_coclust) %>% 
    setNames(c("k", "lambda", "top_marm_clust", "features")) %>% 
    mutate(
        k = sub("k", "", k),
        lambda = sub("lambda", "", lambda),
        features = sub("f", "", features),
        top_marm_clust = as.character(top_marm_clust)
    ) %>% 
    .[, c(1, 2, 4, 3)]

df_tac2_siletti_clust_coclust <- reshape2::melt(tac2_siletti_clust_coclust) %>% 
    setNames(c("k", "lambda", "top_siletti_clust", "features")) %>% 
    mutate(
        k = sub("k", "", k),
        lambda = sub("lambda", "", lambda),
        features = sub("f", "", features),
        top_siletti_clust = as.character(top_siletti_clust)
    ) %>% 
    .[, c(1, 2, 4, 3)]

df_tac2_siletti_subclust_coclust <- reshape2::melt(tac2_siletti_subclust_coclust) %>% 
    setNames(c("k", "lambda", "top_siletti_subclust", "features")) %>% 
    mutate(
        k = sub("k", "", k),
        lambda = sub("lambda", "", lambda),
        features = sub("f", "", features),
        top_siletti_subclust = as.character(top_siletti_subclust)
    ) %>% 
    .[, c(1, 2, 4, 3)]

# ---
# Join together
# ---

df_tac2_marm_siletti_coclust <- left_join(
    df_tac2_coclust,
    df_tac2_marm_clust_coclust,
    by = c("k", "lambda", "features")
) %>% left_join(
    df_tac2_siletti_clust_coclust,
    by = c("k", "lambda", "features")
) %>% left_join(
    df_tac2_siletti_subclust_coclust,
    by = c("k", "lambda", "features")
)


# Print co-clustering summaries -------------------------------------------

# Marmoset cluster co-clustering frequency:
#
# > cbind(table(df_tac2_marm_siletti_coclust$top_marm_clust, useNA = "always"))
#       .
# 3     3
# 4    13
# 5     3
# 7    52
# 8    16
# 9     7
# <NA> 32

# Human cluster co-clustering:
#
# > cbind(table(df_tac2_marm_siletti_coclust$top_siletti_clust, useNA = "always"))
#       .
# 235   9
# 238   9
# 252   1
# 258   2
# 288   1
# 290   1
# 378   1
# 392   5
# 409  59
# 423   2
# 428   1
# 429   4
# <NA> 31

# Human subcluster co-clustering:
# 
# > cbind(table(df_tac2_marm_siletti_coclust$top_siletti_subclust, useNA = "always"))
#       .
# 165   3
# 166  56
# 167   1
# 168   1
# 173   6
# 175   1
# 193   4
# 198   5
# 402   2
# 433   1
# 555   1
# 771   3
# 777   4
# 854   1
# 1502  1
# 1730  1
# 1809  3
# 1812  1
# <NA> 31


# Plot co-clustering histograms -------------------------------------------

# Marmoset clusters
df_tac2_marm_siletti_coclust %>% 
    ggplot(aes(top_marm_clust)) + 
    geom_histogram(
        stat = "count",
        color = "black",
        linewidth = 0.25,
        fill = "gray75"
    ) + 
    coord_cartesian(
        expand = FALSE,
        clip = "off",
        ylim = c(0, 60)
    ) +
    labs(
        x = "Marmoset (Krienen atlas) cluster\nco-clustering with ABC Tac2+ cells",
        y = "Occurrences across 126 integrations"
    ) +
    theme_md_classic()

ggsave(
    file.path(path_out, "freq_abc_tac2_coclust_marm_clust.pdf"),
    width = 90, height = 80, units = "mm"
)

# Human clusters
df_tac2_marm_siletti_coclust %>% 
    ggplot(aes(top_siletti_clust)) + 
    geom_histogram(
        stat = "count",
        color = "black",
        linewidth = 0.25,
        fill = "gray75"
    ) + 
    coord_cartesian(
        expand = FALSE,
        clip = "off",
        ylim = c(0, 60)
    ) +
    labs(
        x = "Human (Siletti atlas) cluster\nco-clustering with ABC Tac2+ cells",
        y = "Occurrences across 126 integrations"
    ) +
    theme_md_classic()

ggsave(
    file.path(path_out, "freq_abc_tac2_coclust_human_clust.pdf"),
    width = 90, height = 80, units = "mm"
)


# Human subclusters
df_tac2_marm_siletti_coclust %>% 
    ggplot(aes(top_siletti_subclust)) + 
    geom_histogram(
        stat = "count",
        color = "black",
        linewidth = 0.25,
        fill = "gray75"
    ) + 
    coord_cartesian(
        expand = FALSE,
        clip = "off",
        ylim = c(0, 60)
    ) +
    labs(
        x = "Human (Siletti atlas) subcluster\nco-clustering with ABC Tac2+ cells",
        y = "Occurrences across 126 integrations"
    ) +
    theme_md_classic() + 
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

ggsave(
    file.path(path_out, "freq_abc_tac2_coclust_human_subclust.pdf"),
    width = 90, height = 80, units = "mm"
)

