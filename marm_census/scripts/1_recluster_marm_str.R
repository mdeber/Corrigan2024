# This script generates the main marmoset dataset used for integrations and hotspot:
#   'marm_str_gad_drop2_dropGlut_0p8'

suppressPackageStartupMessages({
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
# '.../marm_census/data'
# 

path_marm_data <- file.path(dirname(funr::get_script_path()), "data")
path_marm_seurat_str <- file.path(path_marm_data, "striatum.GAD.rds")
path_str_gad_cell_meta <- file.path(path_marm_data, "marm_str_gad_cell_meta.csv")


# Import striatal cells ---------------------------------------------------

str_gad_seurat <- readRDS(path_marm_seurat_str)

# > str_gad_seurat
# An object of class Seurat 
# 19762 features across 6249 samples within 1 assay 
# Active assay: RNA (19762 features, 0 variable features)
# 1 layer present: counts


# Import & join metadata --------------------------------------------------

str_gad_cell_meta <- read.csv(path_str_gad_cell_meta)
rownames(str_gad_cell_meta) <- str_gad_cell_meta[, 1]
str_gad_cell_meta <- str_gad_cell_meta[, -1]

str_gad_seurat <- AddMetaData(str_gad_seurat, str_gad_cell_meta)

# Fix CLUSTER metadata
str_gad_seurat$CLUSTER <- as.factor(
    str_gad_seurat$CLUSTER
)


# Drop unassigned cells ---------------------------------------------------

# > sum(is.na(str_gad_seurat$CLUSTER))
# [1] 342

str_gad_seurat <- str_gad_seurat[, !is.na(str_gad_seurat$CLUSTER)]

# > str_gad_seurat
# An object of class Seurat 
# 19762 features across 5907 samples within 1 assay 
# Active assay: RNA (19762 features, 0 variable features)
# 1 layer present: counts


# Remove cluster 02 (MEIS2+) ----------------------------------------------

str_gad_drop2 <- subset(str_gad_seurat, CLUSTER != 2)

# > str_gad_drop2
# An object of class Seurat 
# 19762 features across 3930 samples within 1 assay 
# Active assay: RNA (19762 features, 0 variable features)
# 1 layer present: counts


# Recluster ---------------------------------------------------------------

str_gad_drop2 <- NormalizeData(str_gad_drop2, verbose = FALSE)
str_gad_drop2 <- FindVariableFeatures(str_gad_drop2, verbose = FALSE)
str_gad_drop2 <- ScaleData(str_gad_drop2, verbose = FALSE)
str_gad_drop2 <- RunPCA(str_gad_drop2, verbose = FALSE)

marm_str_gad_drop2_recl_0p8 <- FindNeighbors(
    str_gad_drop2, dims = 1:30, verbose = FALSE
)
marm_str_gad_drop2_recl_0p8 <- FindClusters(
    marm_str_gad_drop2_recl_0p8, res = 0.8, verbose = FALSE
)
marm_str_gad_drop2_recl_0p8 <- RunUMAP(
    marm_str_gad_drop2_recl_0p8, dims = 1:30, verbose = FALSE
)


# Export RDS, cluster counts ----------------------------------------------

# this isn't used anywhere downstream
if (FALSE) {
    clustdir <- file.path(path_marm_data, "marm_str_gad_drop2_recl_0p8")
    if (!dir.exists(clustdir))
        dir.create(clustdir)
    
    SeuratObject::SaveSeuratRds(
        marm_str_gad_drop2_recl_0p8, 
        file.path(clustdir, "marm_str_gad_drop2_recl_0p8.rds")
    )
    
    write.csv(
        setNames(
            as.data.frame(table(marm_str_gad_drop2_recl_0p8$seurat_clusters)),
            c("cluster", "cells")
        ),
        file.path(clustdir, "cell_counts_by_cluster.csv"),
        row.names = FALSE
    )
}


# Drop excitatory cells ---------------------------------------------------

# Some small clusters appear excitatory; and they look like they might be 
# doublets that weren't removed

# > as.matrix(t(AverageExpression(
#     marm_str_gad_drop2_recl_0p8,
#     features = c("SLC17A6", "SLC17A7"),
#     group.by = "seurat_clusters"
# )[[1]]))
# 
# SLC17A6      SLC17A7
# g0  0.000000000 0.0091909906
# g1  0.001710245 0.0014552859
# g2  0.001284669 0.0025837911
# g3  0.000000000 0.0009428861
# g4  0.000000000 0.0000000000
# g5  0.000661197 0.0039418081
# g6  0.001624780 0.0311567592
# g7  0.000000000 0.0014177745
# g8  0.004838973 0.0207768691
# g9  0.222108826 2.8241632816
# g10 0.000000000 0.0000000000
# g11 0.000000000 0.0187992894
# g12 0.005416097 0.8351321635
# g13 0.007815431 1.6834824317
# g14 0.009264380 0.8360616205
# g15 0.000000000 0.0000000000
# g16 0.000000000 0.0000000000
# g17 0.026912143 0.0000000000
# g18 1.172019331 0.0000000000
# g19 0.213487770 2.4543052440
# g20 0.000000000 0.0000000000
# g21 0.000000000 0.0000000000

marm_drop_clust <- c("9", "12", "13", "14", "18", "19")

# How many cells is that?
# > sum(marm_str_gad_drop2_recl_0p8$seurat_clusters %in% marm_drop_clust)
# [1] 551

# unlike with dataframes, subset() of seurat obj doesn't work with negation
marm_keep_clust <- setdiff(
    levels(marm_str_gad_drop2_recl_0p8$seurat_clusters),
    marm_drop_clust
)
marm_str_gad_drop2_dropGlut <- subset(
    marm_str_gad_drop2_recl_0p8, 
    seurat_clusters %in% marm_keep_clust
)

# > marm_str_gad_drop2_dropGlut
# An object of class Seurat 
# 19762 features across 3379 samples within 1 assay 
# Active assay: RNA (19762 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap


# Recluster again ---------------------------------------------------------

# So note that we're no longer using the cluster names that were found 
# previously

marm_str_gad_drop2_dropGlut <- NormalizeData(
    marm_str_gad_drop2_dropGlut, verbose = FALSE
)
marm_str_gad_drop2_dropGlut <- FindVariableFeatures(
    marm_str_gad_drop2_dropGlut, verbose = FALSE
)
marm_str_gad_drop2_dropGlut <- suppressWarnings(ScaleData(
    marm_str_gad_drop2_dropGlut, verbose = FALSE
))
marm_str_gad_drop2_dropGlut <- RunPCA(
    marm_str_gad_drop2_dropGlut, verbose = FALSE
)

marm_str_gad_drop2_dropGlut_0p8 <- FindNeighbors(
    marm_str_gad_drop2_dropGlut, dims = 1:30, verbose = FALSE
)
marm_str_gad_drop2_dropGlut_0p8 <- FindClusters(
    marm_str_gad_drop2_dropGlut_0p8, res = 0.8, verbose = FALSE
)
marm_str_gad_drop2_dropGlut_0p8 <- RunUMAP(
    marm_str_gad_drop2_dropGlut_0p8, dims = 1:30, verbose = FALSE
)


# Export RDS, cluster counts ----------------------------------------------

clustdir <- file.path(path_marm_data, "marm_str_gad_drop2_dropGlut_0p8")
if (!dir.exists(clustdir))
    dir.create(clustdir)

SeuratObject::SaveSeuratRds(
    marm_str_gad_drop2_dropGlut_0p8,
    file.path(clustdir, "marm_str_gad_drop2_dropGlut_0p8.rds")
)

write.csv(
    setNames(
        as.data.frame(table(marm_str_gad_drop2_dropGlut_0p8$seurat_clusters)),
        c("cluster", "cells")
    ),
    file.path(clustdir, "cell_counts_by_cluster.csv"),
    row.names = FALSE
)


# Marker genes ------------------------------------------------------------

# Using 2 lists for some different plots; marker_genes is a subset of 
# marker_genes2
# 
# Note: marker genes AVP, MIA, and CENPF are missing from marmoset data, but are
# present in the Siletti data

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
    # "MIA",       
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
    # "AVP",       
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
    # "MIA",
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
    # "AVP",
    "PENK",
    "HTR1F",
    "CCK",
    "MEIS2",
    "ASPN",
    # "CENPF",
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


# Plots -------------------------------------------------------------------

plotdir <- file.path(clustdir, "plots")

if (!dir.exists(plotdir))
    dir.create(plotdir)

# Seurat clusters UMAP
pdf(file.path(plotdir, "umap.pdf"))
print(DimPlot(marm_str_gad_drop2_dropGlut_0p8))
invisible(dev.off())

pdf(file.path(plotdir, "umap_labelled.pdf"))
print(DimPlot(marm_str_gad_drop2_dropGlut_0p8, label = TRUE))
invisible(dev.off())


# Marker gene dot plots
pdf(file.path(plotdir, "marker_dotplots.pdf"), width = 12)
print({
    DotPlot(marm_str_gad_drop2_dropGlut_0p8, marker_genes) + 
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
    DotPlot(marm_str_gad_drop2_dropGlut_0p8, marker_genes, cluster.idents = TRUE) + 
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
    print(VlnPlot(marm_str_gad_drop2_dropGlut_0p8, j, raster=FALSE) + theme(
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
    print(FeaturePlot(marm_str_gad_drop2_dropGlut_0p8, j))
    invisible(dev.off())
}


# Metadata violin plots
# (try to plot all metadata columns, and make a pdf if it works)
metadata_violin_dir <- file.path(plotdir, "metadata_violin_plots")
if (!dir.exists(metadata_violin_dir))
    dir.create(metadata_violin_dir)

for (j in names(marm_str_gad_drop2_dropGlut_0p8@meta.data)) {
    suppressWarnings(suppressMessages({
        pdf(NULL)
        plt <- try(
            print(VlnPlot(marm_str_gad_drop2_dropGlut_0p8, j) + theme(axis.text.x = element_text(
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

for (nm in names(marm_str_gad_drop2_dropGlut_0p8@meta.data)) {
    pdf(NULL)
    suppressWarnings(suppressMessages(
        plt <- try(print(FeaturePlot(marm_str_gad_drop2_dropGlut_0p8, nm)), silent = TRUE)
    ))
    if (identical(class(plt), "try-error")) {
        suppressWarnings(suppressMessages(
            plt <- try(print(DimPlot(marm_str_gad_drop2_dropGlut_0p8, group.by = nm)), silent = TRUE)
        ))
    }
    invisible(dev.off())
    if (identical(class(plt), "try-error"))
        next
    pdf(file.path(metadata_umap_dir, paste0(nm, ".pdf")))
    print(plt)
    invisible(dev.off())
}
