# Script to be used following AnnData_to_MatrixMarket.py.
# 
# Args: 
# 1) the MatrixMarket output directory created in that script
# 2) path for the output Seurat RDS file
#
# This script is very basic as-is and doesn't do everything. Everything is 
# for assay "RNA".

cat(paste0("##### Begin script MatrixMarket_to_Seurat.R #####\n"))

suppressPackageStartupMessages({
    library(Matrix)
    library(Seurat)
})
args <- commandArgs(trailingOnly = TRUE)
mm_dir <- args[1]
rds_path <- args[2]

if (!dir.exists(mm_dir))
    stop("MatrixMarket directory (first argument) not found.")

# Import cell metadata
path_cell_meta <- file.path(mm_dir, "cell_metadata.csv")
if (!file.exists(path_cell_meta))
    stop("File 'cell_metadata.csv' not found in MatrixMarket directory.")
cat("Importing 'cell_metadata.csv'...\n")
cell_meta <- read.csv(path_cell_meta)
rownames(cell_meta) <- cell_meta[,1]
cell_meta <- cell_meta[, -1, drop = FALSE]


# Import feature metadata  [no guarantee row.names is correct]
path_feat_meta = file.path(mm_dir, "feature_metadata.csv")
if (!file.exists(path_feat_meta))
    stop("File 'feature_metadata.csv' not found in MatrixMarket directory.")
cat("Importing 'feature_metadata.csv'...\n")
feat_meta <- read.csv(path_feat_meta)
rownames(feat_meta) <- feat_meta[,1]
feat_meta <- feat_meta[, -1, drop = FALSE]


# Import embeddings
embed_dir <- file.path(mm_dir, "cell_embeddings")
if (dir.exists(embed_dir) && length(dir(embed_dir)) > 0) {
    cat("Founding cell_embeddings directory. Importing...\n")
    embeds <- setNames(
        dir(embed_dir, full.names = TRUE), 
        sub(".csv", "", dir(embed_dir))
    )
    embeds <- lapply(embeds, read.csv)
    embeds <- lapply(embeds, function(x) as.matrix(x[,-1]))
    stopifnot(all(sapply(embeds, nrow) == nrow(cell_meta)))
    embeds <- lapply(embeds, `rownames<-`, rownames(cell_meta))
    embeds <- lapply(embeds, `colnames<-`, NULL)
    names(embeds) <- sub("^X_", "", names(embeds))
    embeds <- Map(
        function(mat, nm) {
            embed_key <- paste0(gsub("_", "", nm), "_")
            suppressWarnings(SeuratObject::CreateDimReducObject(
                embeddings = mat,
                assay = "RNA",
                key = embed_key
            ))
        },
        embeds,
        names(embeds)
    )
}

# Import raw data
path_counts <- file.path(mm_dir, "counts.mtx")
if (!file.exists(path_counts))
    stop("File 'counts.mtx' not found in MatrixMarket directory.")
cat("Importing 'counts.mtx'...\n")
counts <- as(t(readMM(path_counts)), "CsparseMatrix")


# Import normalized data
path_data <- file.path(mm_dir, "data.mtx")
if (file.exists(path_data)) {
    cat("Found 'data.mtx' file in MatrixMarket directory. Importing...\n")
    data <- as(t(readMM(path_data)), "CsparseMatrix")
}

# Create Seurat object
cat("Creating Seurat object...\n")
seur <- SeuratObject::CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    meta.data = cell_meta
)
rownames(seur) <- rownames(feat_meta)
colnames(seur) <- rownames(cell_meta)
seur[["RNA"]]@meta.data <- feat_meta
if ("data" %in% names(.GlobalEnv))
    seur[["RNA"]]$data <- data
if ("embeds" %in% names(.GlobalEnv))
    seur@reductions <- embeds

# Export
cat("Exporting Seurat RDS...\n")
SeuratObject::SaveSeuratRds(seur, rds_path)

cat("Script MatrixMarket_to_Seurat.R finished.\n\n")
