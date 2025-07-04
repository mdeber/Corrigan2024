# This script generates the main Siletti dataset used for integrations and hotspot:
#   'siletti_neur_roi_sclust_dropAllMeis2_0p8'

suppressPackageStartupMessages({
    library(Matrix)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(Seurat) # v5.0.3
    library(funr) # for finding local paths
})


# Paths -------------------------------------------------------------------

# 
# Note: 
# funr::get_script_path() works when running script from the command line
# (`Rscript <PATH_TO_THIS_SCRIPT>`) or if sourcing in R
# (`source(<PATH_TO_THIS_SCRIPT>)`), but it won't work if you are running this
# line-by-line. In that case, manually substitute your local path to 
# '.../siletti/data'
# 

path_siletti_data <- file.path(dirname(funr::get_script_path()), "data")
path_mtx     <- file.path(path_siletti_data, 'Neurons_select_roi_X.mtx')
path_obs_csv <- file.path(path_siletti_data, 'Neurons_select_roi_obs.csv')
path_var_csv <- file.path(path_siletti_data, 'Neurons_var.csv')


# Import function for Siletti data ----------------------------------------

import_siletti_mm <- function(path_mtx, path_obs_csv, path_var_csv) {
    # This script is to be run after creating csv and MatrixMarket files in
    # python
    
    # --- Import gene (var) metadata ---
    gene_meta <- read.csv(path_var_csv) # n=59480
    
    # drop genes where multiple ens id's have same gene symbol, and also drop
    # the non-standard genes (includes transgenes & other non-std chromosomes)
    gene_dupes <- unique(gene_meta$Gene[duplicated(gene_meta$Gene)])
    idx_dupe <- which(gene_meta$Gene %in% gene_dupes) # n=1175
    idx_nonstd <- which(
        gene_meta$Chromosome == "chrEXTRA" | grepl("_", gene_meta$Chromosome)
    ) # n=68
    idx_drop <- union(idx_dupe, idx_nonstd) # n=1228
    gene_meta <- gene_meta[-idx_drop, ] # n=58252
    
    # index by gene symbol, and get rid of underscores (which seurat would do)
    rownames(gene_meta) <- gene_meta$Gene
    gene_meta$Gene <- NULL
    rownames(gene_meta) <- gsub("_", "-", rownames(gene_meta))
    
    # --- Import cell (obs) metadata ---
    cell_meta <- read.csv(path_obs_csv, row.names = "CellID")
    
    # --- Import MatrixMarket (mtx) ---
    counts <- as(t(Matrix::readMM(path_mtx)), "CsparseMatrix")[-idx_drop, ]
    stopifnot(identical(
        dim(counts), 
        c(nrow(gene_meta), nrow(cell_meta)))
    )
    colnames(counts) <- rownames(cell_meta)
    rownames(counts) <- rownames(gene_meta)
    
    ### Make Seurat
    seur <- SeuratObject::CreateSeuratObject(
        counts = counts,
        assay = "RNA",
        meta.data = cell_meta
    )
    seur[["RNA"]]@meta.data <- gene_meta
    invisible(gc())
    return(seur)
}


# Import ROI-selected neurons ---------------------------------------------

siletti_neur_roi <- import_siletti_mm(
    path_mtx = path_mtx,
    path_obs_csv = path_obs_csv,
    path_var_csv = path_var_csv
)

# > siletti_neur_roi
# An object of class Seurat 
# 58252 features across 88044 samples within 1 assay 
# Active assay: RNA (58252 features, 0 variable features)
# 1 layer present: counts

# > as.matrix(table(siletti_neur_roi$supercluster_term))
# [,1]
# Amygdala excitatory                  11
# CGE interneuron                     711
# Deep-layer corticothalamic and 6b   130
# Deep-layer intratelencephalic       682
# Deep-layer near-projecting            7
# Eccentric medium spiny neuron      5542
# LAMP5-LHX6 and Chandelier           106
# Medium spiny neuron               73861
# MGE interneuron                     238
# Miscellaneous                        27
# Splatter                           6716
# Upper rhombic lip                     1
# Upper-layer intratelencephalic       12


# Subset to supercluster tersm --------------------------------------------

siletti_neur_roi_sclust <- subset(siletti_neur_roi, supercluster_term %in% c(
    "CGE interneuron", "MGE interneuron", "Splatter"
))
rm(siletti_neur_roi)
invisible(gc())

# > siletti_neur_roi_sclust
# An object of class Seurat 
# 58252 features across 7665 samples within 1 assay 
# Active assay: RNA (58252 features, 0 variable features)
# 1 layer present: counts


# Recluster ---------------------------------------------------------------

siletti_neur_roi_sclust <- NormalizeData(
    siletti_neur_roi_sclust, verbose = FALSE
)
siletti_neur_roi_sclust <- FindVariableFeatures(
    siletti_neur_roi_sclust, verbose = FALSE
)
siletti_neur_roi_sclust <- ScaleData(siletti_neur_roi_sclust, verbose = FALSE)
siletti_neur_roi_sclust <- RunPCA(siletti_neur_roi_sclust, verbose = FALSE)

siletti_neur_roi_sclust_0p8 <- FindNeighbors(
    siletti_neur_roi_sclust, dims = 1:30, verbose = FALSE
)
siletti_neur_roi_sclust_0p8 <- FindClusters(
    siletti_neur_roi_sclust_0p8, res = 0.8, verbose = FALSE
)
siletti_neur_roi_sclust_0p8 <- suppressMessages(RunUMAP(
    siletti_neur_roi_sclust_0p8, dims = 1:30, verbose = FALSE
))


# Export RDS, cluster counts ----------------------------------------------

# this isn't used anywhere downstream
if (FALSE) {
    clust_dir <- file.path(path_siletti_data, "siletti_neur_roi_sclust_0p8")
    if (!dir.exists(clust_dir))
        dir.create(clust_dir)
    
    SeuratObject::SaveSeuratRds(
        siletti_neur_roi_sclust_0p8,
        file.path(clust_dir, "siletti_neur_roi_sclust_0p8.rds")
    )
    
    write.csv(
        setNames(
            as.data.frame(table(siletti_neur_roi_sclust_0p8$seurat_clusters)),
            c("cluster", "cells")
        ),
        file.path(clust_dir, "cell_counts_by_cluster.csv"),
        row.names = FALSE
    )
}


# Remove MEIS2+ cells -----------------------------------------------------

# Subsetting based on my initial reclustering

siletti_neur_roi_sclust_dropAllMeis2_0p8 <- subset(
    siletti_neur_roi_sclust_0p8,
    seurat_clusters %in% setdiff(
        levels(siletti_neur_roi_sclust_0p8$seurat_clusters),
        c("10", "13", "14", "17", "18", "19", "21")
    )
)

# > siletti_neur_roi_sclust_dropAllMeis2_0p8
# An object of class Seurat 
# 58252 features across 6617 samples within 1 assay 
# Active assay: RNA (58252 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap


# Recluster again ---------------------------------------------------------

siletti_neur_roi_sclust_dropAllMeis2_0p8 <- NormalizeData(
    siletti_neur_roi_sclust_dropAllMeis2_0p8, verbose = FALSE
)
siletti_neur_roi_sclust_dropAllMeis2_0p8 <- FindVariableFeatures(
    siletti_neur_roi_sclust_dropAllMeis2_0p8, verbose = FALSE
)
siletti_neur_roi_sclust_dropAllMeis2_0p8 <- suppressWarnings(ScaleData(
    siletti_neur_roi_sclust_dropAllMeis2_0p8, verbose = FALSE
))
siletti_neur_roi_sclust_dropAllMeis2_0p8 <- RunPCA(
    siletti_neur_roi_sclust_dropAllMeis2_0p8, verbose = FALSE
)

siletti_neur_roi_sclust_dropAllMeis2_0p8 <- FindNeighbors(
    siletti_neur_roi_sclust_dropAllMeis2_0p8, dims = 1:30, verbose = FALSE
)
siletti_neur_roi_sclust_dropAllMeis2_0p8 <- FindClusters(
    siletti_neur_roi_sclust_dropAllMeis2_0p8, res = 0.8, verbose = FALSE
)
siletti_neur_roi_sclust_dropAllMeis2_0p8 <- RunUMAP(
    siletti_neur_roi_sclust_dropAllMeis2_0p8, dims = 1:30, verbose = FALSE
)


# Fix metadata numeric-to-factor ------------------------------------------

# > class(siletti_neur_roi_sclust_dropAllMeis2_0p8$cluster_id)
# [1] "integer"
# > class(siletti_neur_roi_sclust_dropAllMeis2_0p8$subcluster_id)
# [1] "integer"

siletti_neur_roi_sclust_dropAllMeis2_0p8$cluster_id <- as.factor(
    siletti_neur_roi_sclust_dropAllMeis2_0p8$cluster_id
)
siletti_neur_roi_sclust_dropAllMeis2_0p8$subcluster_id <- as.factor(
    siletti_neur_roi_sclust_dropAllMeis2_0p8$subcluster_id
)


# Export RDS, cluster counts ----------------------------------------------

clust_dir <- file.path(
    path_siletti_data, "siletti_neur_roi_sclust_dropAllMeis2_0p8"
)
if (!dir.exists(clust_dir))
    dir.create(clust_dir)

SeuratObject::SaveSeuratRds(
    siletti_neur_roi_sclust_dropAllMeis2_0p8,
    file.path(clust_dir, "siletti_neur_roi_sclust_dropAllMeis2_0p8.rds")
)

write.csv(
    setNames(
        as.data.frame(
            table(siletti_neur_roi_sclust_dropAllMeis2_0p8$seurat_clusters)
        ),
        c("cluster", "cells")
    ),
    file.path(clust_dir, "cell_counts_by_cluster.csv"),
    row.names = FALSE
)


# Marker genes ------------------------------------------------------------

# Note: marker genes AVP and MIA were missing from marmoset data, but are
# present in the Siletti data

# Using 2 lists for some different plots; marker_genes is a subset of 
# marker_genes2

marker_genes <- c(
    "FOXG1",     
    "NKX2-1",    
    "LHX6",      
    "LHX8",      
    "STXBP6",    
    "ANGPT2",    
    "TAC3",      
    "TRH",       
    "CHRNA3",    
    "CHRNA4",    
    "CHRNA7",    
    "CHRNB4",    
    "GFRA2",     
    "PRLR",      
    "SYT1",      
    "GRIK1",     
    "MAF1",      
    "ETV1",      
    "CRABP1",    
    "PVALB",     
    "PTHLH",     
    "TH",        
    "CHAT",      
    "NPY",       
    "SST",       
    "MIA",
    "GAD1",      
    "RBP4",      
    "ARX",       
    "NXPH1",     
    "NXPH2",     
    "TACR3",     
    "ZIC1",      
    "GBX1",      
    "ISL1",      
    "PTPRZ1",    
    "SLC35D3",   
    "SORCS1",    
    "PARM1",     
    "CHODL",     
    "COL19A1",   
    "VIP",       
    "SLC17A6",   
    "SLC17A7",   
    "AVP",
    "PENK",      
    "HTR1F",     
    "CCK",       
    "MEIS2",
    "LAMP5"
)

marker_genes2 <- c(
    "FOXG1",
    "NKX2-1",
    "LHX6",
    "LHX8",
    "STXBP6",
    "ANGPT2",
    "TAC3",
    "TRH",
    "CHRNA3",
    "CHRNA4",
    "CHRNA7",
    "CHRNB4",
    "GFRA2",
    "PRLR",
    "SYT1",
    "GRIK1",
    "MAF1",
    "ETV1",
    "CRABP1",
    "PVALB",
    "PTHLH",
    "TH",
    "CHAT",
    "NPY",
    "SST",
    "MIA",
    "GAD1",
    "RBP4",
    "ARX",
    "NXPH1",
    "NXPH2",
    "TACR3",
    "ZIC1",
    "GBX1",
    "ISL1",
    "PTPRZ1",
    "SLC35D3",
    "SORCS1",
    "PARM1",
    "CHODL",
    "COL19A1",
    "VIP",
    "SLC17A6",
    "SLC17A7",
    "AVP",
    "PENK",
    "HTR1F",
    "CCK",
    "MEIS2",
    "ASPN",
    "CENPF",
    "CKB",
    "RBFOX1",
    "SOX6",
    "MEF2C",
    "MAF",
    "NRTN",
    "KIT",
    "TRHDE",
    "PAX6",
    "SCGN",
    "EYA2",
    "FOXP2",
    "TSHZ1",
    "FOXP1",
    "SIX3",
    "RXRG",
    "RARB",
    "ELMOD1",
    "NPY1R",
    "PDZRN3",
    "NR3C2",
    "NR2F2",
    "PROX1",
    "ZIC4",
    "LAMP5"
)

# > length(marker_genes)
# [1] 50

# > length(marker_genes2)
# [1] 76

stopifnot(all(marker_genes %in% marker_genes2))
stopifnot(
    all(marker_genes2 %in% rownames(siletti_neur_roi_sclust_0p8))
)



# Plots -------------------------------------------------------------------

# Using `marker_genes` for the dotplots, but `marker_genes2` for the individual
# UMAP & violin plots

plotdir <- file.path(clust_dir, "plots")

if (!dir.exists(plotdir))
    dir.create(plotdir)

# Seurat clusters UMAP
pdf(file.path(plotdir, "umap.pdf"))
print(DimPlot(siletti_neur_roi_sclust_dropAllMeis2_0p8))
invisible(dev.off())

pdf(file.path(plotdir, "umap_labelled.pdf"))
print(DimPlot(siletti_neur_roi_sclust_dropAllMeis2_0p8, label = TRUE))
invisible(dev.off())


# Marker gene dot plots
pdf(file.path(plotdir, "marker_dotplots.pdf"), width = 12)
print({
    DotPlot(siletti_neur_roi_sclust_dropAllMeis2_0p8, marker_genes) + 
        theme(
            axis.text.x = element_text(
                angle = 90, hjust = 1, vjust = 0.5
            )
        )
})
invisible(dev.off())

# Marker gene dot plots (clustered)
pdf(file.path(plotdir, "marker_dotplots_clustered.pdf"), width = 12)
print({
    DotPlot(siletti_neur_roi_sclust_dropAllMeis2_0p8, marker_genes, cluster.idents = TRUE) + 
        theme(
            axis.text.x = element_text(
                angle = 90, hjust = 1, vjust = 0.5
            )
        )
})
invisible(dev.off())


# Marker gene violin plots
marker_violin_dir <- file.path(plotdir, "marker_violin_plots")
if (!dir.exists(marker_violin_dir))
    dir.create(marker_violin_dir)

for (j in marker_genes2) {
    out <- file.path(marker_violin_dir, paste0(j, ".pdf"))
    if (file.exists(out))
        next
    pdf(out, width = 12)
    print(VlnPlot(siletti_neur_roi_sclust_dropAllMeis2_0p8, j, raster=FALSE) + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ))
    invisible(dev.off())
}


# Marker gene UMAP feature plots
umap_marker_folder <- file.path(plotdir, "marker_umap_plots")
if (!dir.exists(umap_marker_folder))
    dir.create(umap_marker_folder)

for (j in marker_genes2) {
    out <- file.path(umap_marker_folder, paste0(j, ".pdf"))
    if (file.exists(out))
        next
    pdf(out)
    print(FeaturePlot(siletti_neur_roi_sclust_dropAllMeis2_0p8, j))
    invisible(dev.off())
}


# Metadata violin plots
# (try to plot all metadata columns, and make a pdf if it works)
metadata_violin_dir <- file.path(plotdir, "metadata_violin_plots")
if (!dir.exists(metadata_violin_dir))
    dir.create(metadata_violin_dir)

for (j in names(siletti_neur_roi_sclust_dropAllMeis2_0p8@meta.data)) {
    suppressWarnings(suppressMessages({
        pdf(NULL)
        plt <- try(
            print(VlnPlot(siletti_neur_roi_sclust_dropAllMeis2_0p8, j) + 
                      theme(axis.text.x = element_text(
                          angle = 90, hjust = 1, vjust = 0.5
                      ))), 
            silent = TRUE
        )
        invisible(dev.off())
    }))
    if (identical(class(plt), "try-error"))
        next
    pdf(file.path(metadata_violin_dir, paste0(j, ".pdf")), width = 12)
    print(plt)
    invisible(dev.off())
}


# Metadata UMAP feature plots
# (try to plot all metadata columns, and make a pdf if it works)
metadata_umap_dir <- file.path(plotdir, "metadata_umap_plots")
if (!dir.exists(metadata_umap_dir))
    dir.create(metadata_umap_dir)

for (nm in names(siletti_neur_roi_sclust_dropAllMeis2_0p8@meta.data)) {
    pdf(NULL)
    suppressWarnings(suppressMessages(
        plt <- try(print(FeaturePlot(
            siletti_neur_roi_sclust_dropAllMeis2_0p8, nm
        )), silent = TRUE)
    ))
    if (identical(class(plt), "try-error")) {
        suppressWarnings(suppressMessages(
            plt <- try(print(DimPlot(
                siletti_neur_roi_sclust_dropAllMeis2_0p8, group.by = nm
            )), silent = TRUE)
        ))
    }
    invisible(dev.off())
    if (identical(class(plt), "try-error"))
        next
    pdf(file.path(metadata_umap_dir, paste0(nm, ".pdf")))
    print(plt)
    invisible(dev.off())
}


