suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(Seurat) # v5.0.3
    library(funr) # for finding local paths
})

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


# Paths, imports ----------------------------------------------------------

# 
# Note: 
# funr::get_script_path() works when running script from the command line
# (`Rscript <PATH_TO_THIS_SCRIPT>`) or if sourcing in R
# (`source(<PATH_TO_THIS_SCRIPT>)`), but it won't work if you are running this
# line-by-line. In that case, manually substitute your local path to 
# '.../marm_census/data'
# 

path_marm_data <- file.path(dirname(funr::get_script_path()), "data")
clustdir <- file.path(path_marm_data, "marm_str_gad_drop2_dropGlut_0p8")
marm_str_gad_drop2_dropGlut_0p8 <- readRDS(
    file.path(clustdir, "marm_str_gad_drop2_dropGlut_0p8.rds")
)
hsdir <- file.path(clustdir, "hotspot")
str_drop2_dropGlut_hs_modules <- read.csv(file.path(hsdir, "modules.csv"))


# Gene module functions ---------------------------------------------------

RunGeneModuleScore <- function(object, features) {
    # Similar to hotspot module score
    # (didn't do their pre-smoothing across neighboring cells)
    # 
    # The module score is PC1 after doing PCA on only the genes in the module;
    # scores are returned with cell barcodes retained as rownames
    suppressWarnings({
        object <- ScaleData(object, features = features, verbose = FALSE)
        object <- RunPCA(object, features = features, npcs = 2, verbose = FALSE)
    }) # error if try npcs=1
    pc1 <- FetchData(object, vars = 'PC_1') # rownames are cell barcodes
    
    # i think this is theoretically possible
    if (mean(object@reductions$pca@feature.loadings[,1]) < 0)
        pc1[,1] <- -1*pc1[,1]
    
    # rescale min-to-max score from 0 to 1
    pc1[,1] <- (pc1[,1] - min(pc1[,1])) / (max(pc1[,1]) - min(pc1[,1]))
    setNames(pc1, 'module_score')
}

GeneModuleScorePlot <- function(object, features, ...) {
    object$module_score <- RunGeneModuleScore(object, features)$module_score
    FeaturePlot(object, 'module_score', ...)
}


# Module score UMAP plots -------------------------------------------------

dir_hs_mscore_umap <- file.path(hsdir, "md_module_score_umap_plots")

if (!dir.exists(dir_hs_mscore_umap))
    dir.create(dir_hs_mscore_umap)

for (i in seq_len(max(str_drop2_dropGlut_hs_modules$Module))) {
    GeneModuleScorePlot(
        marm_str_gad_drop2_dropGlut_0p8,
        subset(str_drop2_dropGlut_hs_modules, Module == i)$X
    ) + 
        labs(
            color = paste0("module ", i, "\nscore"),
            title = NULL
        ) + 
        theme_md_classic() +
        theme(legend.title = element_text(hjust = 0.5))
    
    ggsave(
        file.path(dir_hs_mscore_umap, paste0(i, ".pdf")),
        width = 150, height = 150, units = "mm"
    )
}


# Export genes by module --------------------------------------------------

dir_hs_genes_by_module <- file.path(hsdir, "genes_by_module")

if (!dir.exists(dir_hs_genes_by_module))
    dir.create(dir_hs_genes_by_module)

str_drop2_dropGlut_hs_modules %>% 
    (function(modules_csv) {
        split(modules_csv[[1]], modules_csv[[2]])
    }) %>% 
    setNames(., paste0("module_", names(.))) %>% 
    setNames(., sub("module_-1", "no_module", names(.))) %>% 
    Map(function(x, nm) {
        out <- file.path(dir_hs_genes_by_module, paste0(nm, ".txt"))
        write(x, out)
    }, ., names(.)) %>% 
    invisible()
