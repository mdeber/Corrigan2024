#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(Seurat) # v5.0.3
    library(ggpubr)
    library(ComplexHeatmap)
})

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


# Paths -------------------------------------------------------------------

# 
# Note: 
# funr::get_script_path() works when running script from the command line
# (`Rscript <PATH_TO_THIS_SCRIPT>`) or if sourcing in R
# (`source(<PATH_TO_THIS_SCRIPT>)`), but it won't work if you are running this
# line-by-line. In that case, manually substitute your local path to this repo.
# 
pdir   <- file.path(funr::get_script_path())
outdir <- file.path(pdir, "correlation_analysis")

if (!dir.exists(outdir))
    dir.create(outdir)

pdirs_data <- list(
    marm    = file.path(pdir, "marm_census/data/marm_str_gad_drop2_dropGlut_0p8"),
    abc     = file.path(pdir, "abc/data/abc_seurat_cl08_str_dropSubc057_0p8"),
    siletti = file.path(pdir, "siletti/data/siletti_neur_roi_sclust_dropAllMeis2_0p8")
)
stopifnot(all(sapply(pdirs_data, file.exists)))


paths_rds <- list(
    marm    = file.path(pdirs_data$marm, "marm_str_gad_drop2_dropGlut_0p8.rds"),
    abc     = file.path(pdirs_data$abc, "abc_seurat_cl08_str_dropSubc057_0p8.rds"),
    siletti = file.path(pdirs_data$siletti, "siletti_neur_roi_sclust_dropAllMeis2_0p8.rds")
)
stopifnot(all(sapply(paths_rds, file.exists)))


paths_clust_anno <- lapply(pdirs_data, file.path, "clust_anno.csv")
stopifnot(all(sapply(paths_clust_anno, file.exists)))


# query is human genes (GRCh38.p14)
path_df_ortho <- file.path(
    pdir, "ortholog_tables/ens_112_ortho_human_ferret_pig_mouse_marmoset.csv"
)
df_ortho <- read.csv(path_df_ortho)



# Import Seurat objects ---------------------------------------------------

seur_list <- lapply(paths_rds, readRDS)


# Import, join cluster annotations ----------------------------------------

clust_anno <- lapply(paths_clust_anno, read.csv)

# for safety (actually not necessary)
clust_anno <- lapply(clust_anno, function(x) {
    names(x) <- sub("^seurat_clusters$", "RNA_snn_res.0.8", names(x))
    x$RNA_snn_res.0.8 <- as.factor(x$RNA_snn_res.0.8)
    x
})

stopifnot(identical(names(seur_list), names(clust_anno)))
seur_list <- Map(
    function(x, df) {
        x@meta.data <- left_join(x@meta.data, df, by = "RNA_snn_res.0.8")
        rownames(x@meta.data) <- Cells(x)
        x
    },
    seur_list,
    clust_anno
)


# Get single ortholog table -----------------------------------------------

# Same thing done in the scVI integration script
map_dat_colname = list(
    marm    = 'White.tufted.ear.marmoset.gene.name',
    abc     = 'Mouse.gene.name',
    siletti = 'Gene.name'
)

# Fix TAC3
# 
# Specifically marmoset is the column missing the value (everything else has 
# it, including mouse as Tac2)
df_ortho$White.tufted.ear.marmoset.gene.name[df_ortho$Gene.name == "TAC3"] <- "TAC3"

# Get single ortholog table
df_ortho_3way <- df_ortho
cols_keep <- unique(c(
    "Gene.name", "Gene.stable.ID", "Gene.stable.ID.version",
    unname(unlist(map_dat_colname))
))
# > cols_keep
# [1] "Gene.name"                           "Gene.stable.ID"                     
# [3] "Gene.stable.ID.version"              "White.tufted.ear.marmoset.gene.name"
# [5] "Mouse.gene.name" 

df_ortho_3way <- unique(df_ortho_3way[, cols_keep])

for (i in seq_along(map_dat_colname)) {
    dname <- names(map_dat_colname)[i]
    cname <- map_dat_colname[[i]]
    df_ortho_3way <- df_ortho_3way[
        df_ortho_3way[[cname]] %in% rownames(seur_list[[dname]]), 
    ]
    dupes <- df_ortho_3way[[cname]][duplicated(df_ortho_3way[[cname]])]
    df_ortho_3way <- df_ortho_3way[
        !df_ortho_3way[[cname]] %in% dupes,
    ]
}


# Subset Seurat & rename genes --------------------------------------------

seur_list <- Map(
    function(x, cname) subset(x, features = df_ortho_3way[[cname]]),
    seur_list,
    map_dat_colname
)

# Note these are not in consistent order, so we want to change the gene names 
# with an explicit rename mapping
#
# this only keeps cell metadata and counts
seur_list <- Map(
    function(x, cname) {
        rename_list <- setNames(df_ortho_3way$Gene.name, df_ortho_3way[[cname]])
        # keep counts
        counts <- x[["RNA"]]$counts
        rownames(counts) <- unname(rename_list[rownames(counts)])
        counts <- counts[order(rownames(counts)), ]
        # also keep data (to keep normalization before subsetting)
        data <- x[["RNA"]]$data
        rownames(data) <- unname(rename_list[rownames(data)])
        data <- data[order(rownames(data)), ]
        # new seurat object
        out <- CreateSeuratObject(counts = counts, meta.data = x@meta.data)
        out[["RNA"]]$data <- data
        out
    }, 
    seur_list, 
    map_dat_colname
)


# Modify cluster annotations ----------------------------------------------

# Marmoset: Define TAC3 LHX8+ vs. LHX8-
seur_list$marm$seurat_clusters_anno_orig <- seur_list$marm$seurat_clusters_anno
seur_list$marm$seurat_clusters_anno[seur_list$marm$seurat_clusters == "8"] <- "TAC3/LHX8+"
seur_list$marm$seurat_clusters_anno[seur_list$marm$seurat_clusters_anno == "TAC3"] <- "TAC3/LHX8-"

# Mouse: Use clust_anno, but ABC clusters for SUBC055 Th cells
seur_list$abc$seurat_clusters_anno_orig <- seur_list$abc$seurat_clusters_anno
seur_list$abc$seurat_clusters_anno <- ifelse(
    seur_list$abc$seurat_clusters_anno == "Th" & seur_list$abc$subclass == "055 STR Lhx8 Gaba",
    paste("ABC", seur_list$abc$cluster),
    seur_list$abc$seurat_clusters_anno
)


# Get metacells -----------------------------------------------------------

avg_exp_by_clust <- seur_list %>% 
    lapply(AverageExpression, group.by = "seurat_clusters_anno") %>% 
    lapply(unname) %>% 
    unlist %>% 
    lapply(as.matrix) %>% 
    Map(
        function(x, nm) {
            nm <- switch(
                nm, 
                "marm"    = "marmoset", 
                "siletti" = "human",
                "abc"     = "mouse"
            )
            colnames(x) <- paste0(nm, ".", colnames(x))
            x
        }, 
        ., 
        names(.)
    )


# Couple quick checks
avg_exp_by_clust %>%
    lapply(rownames) %>%
    .[1:2] %>%
    unname %>%
    do.call(identical, .) %>% 
    stopifnot

avg_exp_by_clust %>%
    lapply(rownames) %>%
    .[2:3] %>%
    unname %>%
    do.call(identical, .) %>% 
    stopifnot

# > colnames(do.call(cbind, avg_exp_by_clust))
#  [1] "marmoset.CCK"                   "marmoset.CHAT"                  "marmoset.Mixed"                
#  [4] "marmoset.PTHLH/PVALB"           "marmoset.SST/NPY"               "marmoset.TAC3/LHX8-"           
#  [7] "marmoset.TAC3/LHX8+"            "mouse.ABC 0839 STR Lhx8 Gaba-1" "mouse.ABC 0840 STR Lhx8 Gaba-1"
# [10] "mouse.ABC 0841 STR Lhx8 Gaba-1" "mouse.ABC 0842 STR Lhx8 Gaba-1" "mouse.ABC 0843 STR Lhx8 Gaba-2"
# [13] "mouse.ABC 0844 STR Lhx8 Gaba-2" "mouse.Chat"                     "mouse.Lhx6/Prox1"              
# [16] "mouse.Pthlh/Pvalb"              "mouse.Sst/Npy"                  "mouse.Th"                      
# [19] "human.CCK"                      "human.CCK/VIP"                  "human.CHAT"                    
# [22] "human.Hypothalamus"             "human.Mixed"                    "human.PTHLH/PVALB"             
# [25] "human.SST/NPY"                  "human.TAC3" 



# Prune non-target clusters -----------------------------------------------

# Dropping these non-target clusters as in other analyses
# Plus the non-subc55 Th cells (as it's only a very few cells)
avg_exp_by_clust %>% 
    do.call(cbind, .) %>% 
    cor() %>% 
    rownames %>% 
    (function(x) {
        stopifnot(all(rownames(x) == colnames(x)))
        which(!x %in% c("human.Hypothalamus", "human.Mixed", "marmoset.Mixed", "mouse.Th"))
    }) -> idx_keep


# Heatmap Pearson correlation (all genes) ---------------------------------

avg_exp_by_clust %>% 
    do.call(cbind, .) %>% 
    cor() %>% 
    .[idx_keep, idx_keep] %>% 
    rownames %>% 
    (function(x) {
        ifelse(
            grepl("mouse", x), 
            "mouse", 
            ifelse(grepl("TAC3", x), "tac3", "other")
        )
    }) ->
    row_groups

# Match colors: Colorbrewer RdYlBu_r
col_fun <- circlize::colorRamp2(
    seq(0.4, 1, length.out = 11), 
    rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
)

pdf(
    file.path(outdir, "heatmap_pearson_cor_cluster2_avg_for_all_genes_splitmouse_splitTAC3_prune2_RdYlBu.pdf"),
    width = 9, height = 7
)
avg_exp_by_clust %>% 
    do.call(cbind, .) %>% 
    cor() %>% 
    .[idx_keep, idx_keep] %>% 
    ### ALSO RENAME SOME ROWS
    (function(mat) {
        x <- rownames(mat)
        x <- sub("mouse.ABC ", "mouse.ABC_", x)
        x <- sub(" .*", "", x)
        rownames(mat) <- colnames(mat) <- x
        mat
    }) %>% 
    Heatmap(
        .,
        heatmap_legend_param = list(title = "expression\ncorrelation"),
        row_split = row_groups,
        column_split = row_groups,
        col = col_fun
    )
dev.off()



# Inset: point plot all mouse-primate cor for selected lines --------------

# For the main mouse clusters, want to show all the mouse-primate correlations
# (this is for the figure inset)

avg_exp_by_clust %>% 
    do.call(cbind, .) %>% 
    cor() %>% 
    .[idx_keep, idx_keep] %>% 
    
    # rename some rows
    (function(mat) {
        x <- rownames(mat)
        x <- sub("mouse.ABC ", "mouse.ABC_", x)
        x <- sub(" .*", "", x)
        rownames(mat) <- colnames(mat) <- x
        mat
    }) %>% 
    
    # get a dataframe for all non-redundant value pairs
    reshape2::melt(
        varnames = c("mouse", "primate"), 
        value.name = "pearson"
    ) %>% 
    subset(grepl("^mouse", mouse) & !grepl("^mouse", primate)) %>% 
    subset(mouse %in% c(
        "mouse.Chat",
        "mouse.Sst/Npy",
        "mouse.Pthlh/Pvalb",
        "mouse.ABC_0840",
        "mouse.ABC_0841",
        "mouse.ABC_0842"
    )) %>% 
    mutate(
        species = sub("\\..*", "", primate),
        homology = "nonhomolog",
        homology = ifelse(
            toupper(sub(".*\\.", "", mouse)) == sub(".*\\.", "", primate),
            "homolog",
            homology
        ),
        homology = ifelse(
            mouse %in% c("mouse.ABC_0840", "mouse.ABC_0841", "mouse.ABC_0842") & grepl("marmoset.TAC3", primate),
            "homolog",
            homology
        ),
        homology = ifelse(
            mouse %in% c("mouse.ABC_0840", "mouse.ABC_0841", "mouse.ABC_0842") & primate == "human.TAC3",
            "homolog",
            homology
        )
    ) %>% 
    mutate(
        query = sub(".*\\.", "", mouse),
        primate_clust = sub(".*\\.", "", primate)
    ) -> df


clust_cols <- c(
    "CCK"         = "#9370DB",
    "CHAT"        = "#008000",
    "PTHLH/PVALB" = "#4B0082",
    "SST/NPY"     = "#008080",
    "TAC3/LHX8-"  = "#FF00FF",
    "TAC3"        = "#FF00FF",
    "CCK/VIP"     = "#8A2BE2",
    "TAC3/LHX8+"  = "#F5BF82"
)


p <- ggplot() + 
    geom_point(
        data = df,
        aes(query, pearson, color = primate_clust, shape = species, group = species),
        position = position_dodge(width = 0.5),
    ) + 
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_color_manual(values = clust_cols) +
    scale_shape_manual(values = c(19, 1)) +
    labs(
        x = NULL, 
        y = "Mouse-Primate Pearson Correlation"
    ) +
    theme_md_classic() + 
    theme(
        axis.text.x  = element_text(angle = 30, hjust = 1),
        axis.title   = element_text(size = 5),
        axis.text    = element_text(size = 5),
        legend.text  = element_text(size = 5),
        legend.title = element_text(size = 5)
    )
ggsave(
    file.path(outdir, "pointplot_mouse_v_primate_cor_cluster2_avg_for_all_genes_all_clust.pdf"),
    p, width = 70, height = 60, units = "mm"
)

