# This script generates the main ABC dataset used for integrations and hotspot:
#   'abc_seurat_cl08_str_dropSubc057_0p8'

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(Seurat) # v5.0.3
    library(funr) # for finding local paths
})
# also package 'tibble' is used without importing its namespace


# Paths -------------------------------------------------------------------

# 
# Note: 
# funr::get_script_path() works when running script from the command line
# (`Rscript <PATH_TO_THIS_SCRIPT>`) or if sourcing in R
# (`source(<PATH_TO_THIS_SCRIPT>)`), but it won't work if you are running this
# line-by-line. In that case, manually substitute your local path to 
# '.../abc/data'
# 

path_abc_data <- file.path(dirname(funr::get_script_path()), "data")
path_abc_seurat_str <- file.path(path_abc_data, "WMB-10Xv3-STR-raw.rds")
path_abc_cell_meta_cluster_anno <- file.path(
    path_abc_data, "cell_metadata_with_cluster_annotation.csv"
)


# Import ABC atlas cell metadata ------------------------------------------

# the main table that we need:
df_cell_meta_cluster_anno <- read.csv(path_abc_cell_meta_cluster_anno)


# Import striatal cells ---------------------------------------------------

abc_str <- readRDS(path_abc_seurat_str)

# > abc_str
# An object of class Seurat 
# 32285 features across 285167 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts


# Add metadata ------------------------------------------------------------

abc_str$feature_matrix_label <- "WMB-10Xv3-STR" # to match metadata
abc_str$chemistry <- "10Xv3"

# > names(abc_str@meta.data) %>%
#     (function(x) x[x %in% names(df_cell_meta_cluster_anno)])
# [1] "cell_barcode"         "library_label"        "feature_matrix_label"

abc_str@meta.data <- left_join(
    abc_str@meta.data,
    df_cell_meta_cluster_anno[, c(
        "cell_barcode",
        "library_label",
        "feature_matrix_label",
        "donor_label",
        "donor_genotype",
        "donor_sex",
        "cluster_alias",
        "neurotransmitter",
        "class",
        "subclass",
        "supertype",
        "cluster"
    )],
    by = join_by(cell_barcode, library_label, feature_matrix_label)
)
rownames(abc_str@meta.data) <- Cells(abc_str)



# Subset to class 08 ------------------------------------------------------

abc_seurat_cl08_str_recl <- subset(abc_str, class == "08 CNU-MGE GABA")

# > abc_seurat_cl08_str_recl
# An object of class Seurat 
# 32285 features across 4785 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts


# Change gene names to use symbol as identifier ---------------------------

# Must subset and rename genes
# 
# Can't have multiple genes with same symbol:
# > sum(duplicated(abc_seurat_cl08_str_recl[["RNA"]][[]]$gene_symbol))
# [1] 40

abc_seurat_cl08_str_recl[["RNA"]][[]] <- abc_seurat_cl08_str_recl[["RNA"]][[]] %>% 
    tibble::rownames_to_column("ENSEMBL")

# can't seem to remove duplicates within a column; so will remove them by ENS
# ID first
gn_table <- abc_seurat_cl08_str_recl[["RNA"]][[]]
idx_drop <- which(duplicated(gn_table$gene_symbol))

gn_table <- gn_table[-idx_drop, ]
stopifnot(!any(duplicated(gn_table$gene_symbol)))

gn_table <- gn_table %>% 
    tibble::remove_rownames() %>% 
    tibble::column_to_rownames("gene_symbol")

counts_mat <- abc_seurat_cl08_str_recl[["RNA"]]$counts[-idx_drop, ]
# dim(counts_mat)
# [1] 32245  4785
rownames(counts_mat) <- rownames(gn_table)

abc_seurat_cl08_str_recl <- CreateSeuratObject(
    counts_mat,
    meta.data = abc_seurat_cl08_str_recl@meta.data
)


# Re-normalize ------------------------------------------------------------

abc_seurat_cl08_str_recl <- NormalizeData(abc_seurat_cl08_str_recl)


# Remove subclass057 ------------------------------------------------------

abc_seurat_cl08_str_dropSubc057 <- subset(
    abc_seurat_cl08_str_recl,
    subclass != "057 NDB-SI-MA-STRv Lhx8 Gaba"
)


# Re-cluster, export (with subc057) ---------------------------------------

# abc_seurat_cl08_str_recl <- FindVariableFeatures(abc_seurat_cl08_str_recl)
# abc_seurat_cl08_str_recl <- ScaleData(abc_seurat_cl08_str_recl)
# abc_seurat_cl08_str_recl <- RunPCA(abc_seurat_cl08_str_recl)
# 
# abc_seurat_cl08_str_recl <- FindNeighbors(
#     abc_seurat_cl08_str_recl, dims = 1:30
# )
# abc_seurat_cl08_str_recl_0p8 <- FindClusters(
#     abc_seurat_cl08_str_recl, res = 0.8
# )
# abc_seurat_cl08_str_recl_0p8 <- RunUMAP(
#     abc_seurat_cl08_str_recl_0p8, dims = 1:30
# )
# 
# 
# # export RDS, cluster counts
# clust_dir <- file.path(path_abc_data, "abc_seurat_cl08_str_recl_0p8")
# if (!dir.exists(clust_dir))
#     dir.create(clust_dir)
# 
# SeuratObject::SaveSeuratRds(
#     abc_seurat_cl08_str_recl_0p8,
#     file.path(clust_dir, "abc_seurat_cl08_str_recl_0p8.rds")
# )
# 
# write.csv(
#     setNames(
#         as.data.frame(table(abc_seurat_cl08_str_recl_0p8$seurat_clusters)),
#         c("cluster", "cells")
#     ),
#     file.path(clust_dir, "cell_counts_by_cluster.csv"),
#     row.names = FALSE
# )



# Re-cluster, export (without subc057) ------------------------------------

abc_seurat_cl08_str_dropSubc057 <- FindVariableFeatures(
    abc_seurat_cl08_str_dropSubc057
)
abc_seurat_cl08_str_dropSubc057 <- ScaleData(abc_seurat_cl08_str_dropSubc057)
abc_seurat_cl08_str_dropSubc057 <- RunPCA(abc_seurat_cl08_str_dropSubc057)

abc_seurat_cl08_str_dropSubc057 <- FindNeighbors(
    abc_seurat_cl08_str_dropSubc057, dims = 1:30
)
abc_seurat_cl08_str_dropSubc057_0p8 <- FindClusters(
    abc_seurat_cl08_str_dropSubc057, res = 0.8
)
abc_seurat_cl08_str_dropSubc057_0p8 <- RunUMAP(
    abc_seurat_cl08_str_dropSubc057_0p8, dims = 1:30
)

# export RDS, cluster counts
clust_dir <- file.path(path_abc_data, "abc_seurat_cl08_str_dropSubc057_0p8")
if (!dir.exists(clust_dir))
    dir.create(clust_dir)

SeuratObject::SaveSeuratRds(
    abc_seurat_cl08_str_dropSubc057_0p8,
    file.path(clust_dir, "abc_seurat_cl08_str_dropSubc057_0p8.rds")
)

write.csv(
    setNames(
        as.data.frame(table(abc_seurat_cl08_str_dropSubc057_0p8$seurat_clusters)),
        c("cluster", "cells")
    ),
    file.path(clust_dir, "cell_counts_by_cluster.csv"),
    row.names = FALSE
)


# Marker genes ------------------------------------------------------------

# Using 2 lists for some different plots; `marker_genes` is a subset of 
# `marker_genes2`
marker_genes <- c(
    "Foxg1",
    "Nkx2-1",
    "Lhx6",
    "Lhx8",
    "Stxbp6",
    "Angpt2",
    "Tac2",
    "Trh",
    "Chrna3",
    "Chrna4",
    "Chrna7",
    "Chrnb4",
    "Gfra2",
    "Prlr",
    "Syt1",
    "Grik1",
    "Maf1",
    "Etv1",
    "Crabp1",
    "Pvalb",
    "Pthlh",
    "Th",
    "Chat",
    "Npy", 
    "Sst",
    "Mia",
    "Gad1",
    "Rbp4",
    "Arx",
    "Nxph1",
    "Nxph2",
    "Tacr3",
    "Zic1",
    "Gbx1",
    "Isl1",
    "Ptprz1",
    "Slc35d3",
    "Sorcs1",
    "Parm1",
    "Chodl",
    "Col19a1",
    "Vip",
    "Slc17a6",
    "Slc17a7",
    "Avp",
    "Penk",
    "Htr1f",
    "Cck",
    "Meis2"
)

marker_genes2 <- c(
    "Foxg1",
    "Nkx2-1",
    "Lhx6",
    "Lhx8",
    "Stxbp6",
    "Angpt2",
    "Tac2",
    "Trh",
    "Chrna3",
    "Chrna4",
    "Chrna7",
    "Chrnb4",
    "Gfra2",
    "Prlr",
    "Syt1",
    "Grik1",
    "Maf1",
    "Etv1",
    "Crabp1",
    "Pvalb",
    "Pthlh",
    "Th",
    "Chat",
    "Npy",
    "Sst",
    "Mia",
    "Gad1",
    "Rbp4",
    "Arx",
    "Nxph1",
    "Nxph2",
    "Tacr3",
    "Zic1",
    "Gbx1",
    "Isl1",
    "Ptprz1",
    "Slc35d3",
    "Sorcs1",
    "Parm1",
    "Chodl",
    "Col19a1",
    "Vip",
    "Slc17a6",
    "Slc17a7",
    "Avp",
    "Penk",
    "Htr1f",
    "Cck",
    "Meis2",
    "Aspn",
    "Cenpf",
    "Ckb",
    "Rbfox1",
    "Sox6",
    "Mef2c",
    "Maf",
    "Nrtn",
    "Kit",
    "Trhde",
    "Pax6",
    "Scgn",
    "Eya2",
    "Foxp2",
    "Tshz1",
    "Foxp1",
    "Six3",
    "Rxrg",
    "Rarb",
    "Elmod1",
    "Npy1r",
    "Pdzrn3",
    "Nr3c2",
    "Nr2f2",
    "Prox1",
    "Zic4",
    "Lamp5"
)

# > length(marker_genes)
# [1] 49

# length(marker_genes2)
# [1] 76

stopifnot(all(marker_genes %in% marker_genes2))
stopifnot(
    all(marker_genes2 %in% rownames(abc_seurat_cl08_str_dropSubc057_0p8))
)


# Plots (without subc057) -------------------------------------------------

# Using `marker_genes` for the dotplots, but `marker_genes2` for the individual
# UMAP & violin plots

plot_folder <- file.path(clust_dir, "plots")

if (!dir.exists(plot_folder))
    dir.create(plot_folder)

# Seurat clusters UMAP
pdf(file.path(plot_folder, "umap.pdf"))
print(DimPlot(abc_seurat_cl08_str_dropSubc057_0p8))
invisible(dev.off())

pdf(file.path(plot_folder, "umap_labelled.pdf"))
print(DimPlot(abc_seurat_cl08_str_dropSubc057_0p8, label = TRUE))
invisible(dev.off())


# Marker gene dot plots
pdf(file.path(plot_folder, "marker_dotplots.pdf"), width = 12)
print({
    DotPlot(abc_seurat_cl08_str_dropSubc057_0p8, marker_genes) + 
        theme(
            axis.text.x = element_text(
                angle = 90, hjust = 1, vjust = 0.5
            )
        )
})
invisible(dev.off())

# Marker gene dot plots (clustered)
pdf(file.path(plot_folder, "marker_dotplots_clustered.pdf"), width = 12)
print({
    DotPlot(abc_seurat_cl08_str_dropSubc057_0p8, marker_genes, cluster.idents = TRUE) + 
        theme(
            axis.text.x = element_text(
                angle = 90, hjust = 1, vjust = 0.5
            )
        )
})
invisible(dev.off())


# Marker gene violin plots
marker_violin_dir <- file.path(plot_folder, "marker_violin_plots")
if (!dir.exists(marker_violin_dir))
    dir.create(marker_violin_dir)

for (j in marker_genes2) {
    out <- file.path(marker_violin_dir, paste0(j, ".pdf"))
    if (file.exists(out))
        next
    pdf(out, width = 12)
    print(VlnPlot(abc_seurat_cl08_str_dropSubc057_0p8, j, raster=FALSE) + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ))
    invisible(dev.off())
}


# Marker gene UMAP feature plots
umap_marker_folder <- file.path(plot_folder, "marker_umap_plots")
if (!dir.exists(umap_marker_folder))
    dir.create(umap_marker_folder)

for (j in marker_genes2) {
    out <- file.path(umap_marker_folder, paste0(j, ".pdf"))
    if (file.exists(out))
        next
    pdf(out)
    print(FeaturePlot(abc_seurat_cl08_str_dropSubc057_0p8, j))
    invisible(dev.off())
}

# Metadata UMAP feature plots
# (try to plot all metadata columns, and make a pdf if it works)
umap_metadata_folder <- file.path(plot_folder, "metadata_umap_plots")
if (!dir.exists(umap_metadata_folder))
    dir.create(umap_metadata_folder)

for (nm in names(abc_seurat_cl08_str_dropSubc057_0p8@meta.data)) {
    if (nm == "cell_barcode")
        next
    
    pdf(NULL)
    suppressWarnings(suppressMessages(
        plt <- try(
            print(FeaturePlot(abc_seurat_cl08_str_dropSubc057_0p8, nm)), 
            silent = TRUE
        )
    ))
    if (identical(class(plt), "try-error")) {
        suppressWarnings(suppressMessages(
            plt <- try(
                print(DimPlot(abc_seurat_cl08_str_dropSubc057_0p8, group.by = nm)), 
                silent = TRUE
            )
        ))
    }
    invisible(dev.off())
    if (identical(class(plt), "try-error"))
        next
    pdf(file.path(umap_metadata_folder, paste0(nm, ".pdf")))
    print(plt)
    invisible(dev.off())
}
