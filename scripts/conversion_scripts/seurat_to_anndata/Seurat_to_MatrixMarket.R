# This script takes in a Seurat RDS file and exports matrix market files with
# csvs for metadata. It's a work in progress in terms of the number of 
# seurat attributes it exports.
#
# Args:
# 1) Path to seurat rds
# 2) Output directory. If it doesn't exist, it will be created
#
# Currently embeddings are just exported separately from the assays, although
# it can be checked e.g., seur@reductions[[reduction]]@assay.used, so they
# could be grouped by assay...

cat(paste0("##### Begin script Seurat_to_MatrixMarket.R #####\n"))

suppressPackageStartupMessages({library(Seurat)})
args <- commandArgs(trailingOnly = TRUE)
rds_path <- args[1]
out_dir <- args[2]

cat(paste0("Loading ", rds_path, "...\n"))
seur <- readRDS(rds_path)

if (dir.exists(out_dir)) {
    cat("Export directory already exists.\n")
} else {
    cat("Creating export directory...\n")
    dir.create(out_dir)
}

# Cell Metadata
cat("Exporting cell metadata...\n")
write.csv(
    seur@meta.data,
    file.path(out_dir, "cell_metadata.csv")
)

# Assays
cat("Processing assays...\n")
for (assay in names(seur@assays)) {
    assay_dir <- file.path(out_dir, assay)
    if (dir.exists(assay_dir)) {
        cat(paste0("Found existing export directory for '", assay, "'.\n"))
    } else {
        cat(paste0("Creating export directory for '", assay, "'...\n"))
        dir.create(assay_dir)
    }
    
    cat(paste0("Exporting feature metadata for assay '", assay, "'...\n"))
    fmeta <- seur[[assay]]@meta.data
    rownames(fmeta) <- rownames(seur)
    write.csv(fmeta, file.path(assay_dir, "feature_metadata.csv"))
    
    for (layer in Layers(seur[[assay]])) {
        if (layer == "scale.data") {
            cat("Skipping layer 'scale.data'...\n")
        } else {
            cat(paste0(
                "Exporting MatrixMarket for layer '", 
                layer, "' in assay '", assay, "'...\n"
            ))
            Matrix::writeMM(
                seur[[assay]][layer],
                file.path(assay_dir, paste0(layer, ".mtx"))
            )
        }
    }
}
cat("Finished processing assays.\n")

# Reductions
embed_dir <- file.path(out_dir, "cell_embeddings")
if (dir.exists(embed_dir)) {
    cat("Found existing export directory for cell embeddings.\n")
} else {
    cat("Creating export directory for cell embeddings...\n")
    dir.create(embed_dir)
}
for (reduction in Reductions(seur)) {
    cat(paste0("Exporting embeddings from reduction '", reduction, "'...\n"))
    write.csv(
        seur@reductions[[reduction]]@cell.embeddings,
        file.path(embed_dir, paste0(reduction, ".csv"))
    )
}

cat("Script Seurat_to_MatrixMarket.R finished.\n")
