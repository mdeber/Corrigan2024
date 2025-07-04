#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
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


# Paths -------------------------------------------------------------------

# 
# Note: 
# funr::get_script_path() works when running script from the command line
# (`Rscript <PATH_TO_THIS_SCRIPT>`) or if sourcing in R
# (`source(<PATH_TO_THIS_SCRIPT>)`), but it won't work if you are running this
# line-by-line. In that case, manually substitute your local path to this repo.
# 
pdir   <- file.path(funr::get_script_path())
outdir <- file.path(pdir, "barplots_type_proportions")

if (!dir.exists(outdir))
    dir.create(outdir)

pdirs_data <- list(
    marm    = file.path(pdir, "marm_census/data/marm_str_gad_drop2_dropGlut_0p8"),
    abc     = file.path(pdir, "abc/data/abc_seurat_cl08_str_dropSubc057_0p8"),
    siletti = file.path(pdir, "siletti/data/siletti_neur_roi_sclust_dropAllMeis2_0p8")
)
stopifnot(all(sapply(pdirs_data, file.exists)))


paths_rds <- list(
    marm    = file.path(pdirs_data$marm,    "marm_str_gad_drop2_dropGlut_0p8.rds"),
    abc     = file.path(pdirs_data$abc,     "abc_seurat_cl08_str_dropSubc057_0p8.rds"),
    siletti = file.path(pdirs_data$siletti, "siletti_neur_roi_sclust_dropAllMeis2_0p8.rds")
)
stopifnot(all(sapply(paths_rds, file.exists)))


paths_clust_anno <- lapply(pdirs_data, file.path, "clust_anno.csv")
stopifnot(all(sapply(paths_clust_anno, file.exists)))


path_hmba_cell_meta <- file.path(pdir, "marm_hmba/data/cells_annotated.csv")


experiment_colors <- c(
    marm      = "#1F77B4",
    abc       = "#FF7F0E",
    siletti   = "#2CA02C",
    marm_hmba = "#D3D3D3"
)


# to match new labels used for figure outputs below
experiment_colors_rename <- experiment_colors
names(experiment_colors_rename) <- sapply(
    names(experiment_colors_rename),
    switch,
    "abc"       = "Yao2023",
    "marm_hmba" = "Marmoset HMBA",
    "marm"      = "Krienen2023",
    "siletti"   = "Siletti2023"
)


# Imports (except HMBA) ---------------------------------------------------

cells_anno <- lapply(paths_rds, function(x) readRDS(x)@meta.data)
names(cells_anno) <- names(paths_rds)

clust_anno <- lapply(
    paths_clust_anno, 
    function(x) setNames(read.csv(x), c("seurat_clusters", "clust_anno"))
)
names(clust_anno) <- names(paths_clust_anno)


# Join in cluster names (clust_anno)
stopifnot(identical(names(cells_anno), names(clust_anno)))
cells_anno <- lapply(cells_anno, mutate, seurat_clusters = as.factor(seurat_clusters))
clust_anno <- lapply(clust_anno, mutate, seurat_clusters = as.factor(seurat_clusters))
cells_anno <- Map(left_join, cells_anno, clust_anno, by = list("seurat_clusters"))



# Import HMBA -------------------------------------------------------------

cells_anno$marm_hmba <- read.csv(path_hmba_cell_meta)

# Remove LHX8 subtype distinction from marm_hmba
cells_anno$marm_hmba$clust_anno_orig <- cells_anno$marm_hmba$clust_anno
cells_anno$marm_hmba$clust_anno[cells_anno$marm_hmba$clust_anno == "TAC3/LHX8"] <- "TAC3"


# Define putative ABC TAC3 homologs ---------------------------------------

# Define TAC3 as ABC clusters 0840+0841+0842

# Save original abc annotations
cells_anno$abc$clust_anno_orig <- cells_anno$abc$clust_anno

# overwrite
cells_anno$abc$clust_anno <- ifelse(
    cells_anno$abc$cluster %in% c("0840 STR Lhx8 Gaba_1", "0841 STR Lhx8 Gaba_1", "0842 STR Lhx8 Gaba_1"),
    "TAC3", 
    cells_anno$abc$clust_anno_orig
)


# Cell type ratios (TAC3 to each major class) -----------------------------

cells_anno %>% 
    lapply(mutate, clust_anno = toupper(clust_anno)) %>% 
    Map(function(x, nm) mutate(x, experiment = nm), ., names(.)) %>% 
    lapply(. %>% .[, c("experiment", "clust_anno")]) %>% 
    unname %>% 
    do.call(rbind, .) %>% 
    table  %>% 
    (function(x) data.frame(
        experiment = rownames(x),
        tac3_pvalb = apply(x, 1, function(y) y["TAC3"] / y["PTHLH/PVALB"]),
        tac3_sst = apply(x, 1, function(y) y["TAC3"] / y["SST/NPY"]),
        tac3_chat = apply(x, 1, function(y) y["TAC3"] / y["CHAT"])
    )) ->
    df_tac3_ratios2

# # the initial results of table():
#            clust_anno
# experiment   CCK CCK/VIP CHAT HYPOTHALAMUS LAMP5/LHX6 LHX6/PROX1 MIXED PTHLH/PVALB SST/NPY TAC3   TH
#   abc          0       0  285            0          0        434     0         969    1222  321  355
#   marm       135       0  279            0          0          0    98        1234     522 1111    0
#   marm_hmba 1111    3846  645            0       2143          0     0        1880    1028 1529    0
#   siletti    190     230  196          339          0          0   236        2564     399 2463    0
# 
# > print.data.frame(df_tac3_ratios2, row.names=FALSE)
#  experiment tac3_pvalb  tac3_sst  tac3_chat
#         abc  0.2012384 0.1595745  0.6842105
#        marm  0.9003241 2.1283525  3.9820789
#   marm_hmba  0.8132979 1.4873541  2.3705426
#     siletti  0.9606084 6.1729323 12.5663265


df_tac3_ratios2 %>% 
    reshape2::melt(.) %>% 
    mutate(
        variable = sub("tac3_pvalb", "TAC3 : PTHLH/PVALB", variable),
        variable = sub("tac3_sst",   "TAC3 : SST/NPY", variable),
        variable = sub("tac3_chat",  "TAC3 : CHAT", variable),
        experiment = sub("abc", "Yao2023", experiment),
        experiment = sub("marm_hmba", "Marmoset HMBA", experiment),
        experiment = sub("marm", "Krienen2023", experiment),
        experiment = sub("siletti", "Siletti2023", experiment)
    ) %>% 
    ggplot(aes(experiment, value, fill = experiment)) + 
    facet_wrap(~variable, scales = "free") +
    geom_col(position = "dodge") +
    scale_fill_manual(values = experiment_colors_rename) +
    scale_y_continuous(
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.1))
    ) +
    labs(
        x = NULL,
        y = "cell count ratio",
        fill = NULL
    ) +
    theme_md_classic() + 
    theme(
        legend.position = "none",
        axis.text.x  = element_text(angle = 30, hjust = 1),
        strip.text.x = element_text(size = 5),
        axis.title   = element_text(size = 5),
        axis.text    = element_text(size = 5)
    ) -> 
    p

ggsave(
    file.path(outdir, "barplot_cellratio_tac3_to_each_pvalb_sst_chat.pdf"),
    p, width = 80, height = 40, units = "mm"
)


# Proportion of inhibitory interneurons -----------------------------------

# only keep common types in the denominator

cells_anno %>% 
    lapply(mutate, clust_anno = toupper(clust_anno)) %>% 
    Map(function(x, nm) mutate(x, experiment = nm), ., names(.)) %>% 
    lapply(. %>% .[, c("experiment", "clust_anno")]) %>% 
    unname %>% 
    do.call(rbind, .) %>% 
    table %>% 
    .[, !colnames(.) %in% c("HYPOTHALAMUS", "CCK/VIP", "LAMP5/LHX6")] %>% 
    (function(x) data.frame(
        experiment = rownames(x),
        ratio = apply(x, 1, function(y) y["TAC3"] / sum(y))
    )) %>%  
        mutate(
        experiment = sub("abc", "Yao2023", experiment),
        experiment = sub("marm_hmba", "Marmoset HMBA", experiment),
        experiment = sub("marm", "Krienen2023", experiment),
        experiment = sub("siletti", "Siletti2023", experiment)
    ) %>% 
    ggplot(aes(experiment, ratio, fill = experiment)) + 
    geom_col(position = "dodge") +
    scale_fill_manual(values = experiment_colors_rename) +
    scale_y_continuous(
        labels = . %>% scales::percent(., big.mark = ""),
        limits = c(0, 1),
        expand = c(0, 0)
    ) +
    labs(
        x = NULL,
        y = "TAC3 cell type abundance",
        fill = NULL
    ) +
    theme_md_classic() + 
    theme(
        legend.position = "none",
        axis.text.x  = element_text(angle = 30, hjust = 1),
        strip.text.x = element_text(size = 5),
        axis.title   = element_text(size = 5),
        axis.text    = element_text(size = 5)
    ) -> p

ggsave(
    file.path(outdir, "barplot_cellratio_tac3_to_combined_drop_HYPO_VIP_LAMP5.pdf"),
    p, width = 40, height = 40, units = "mm"
)
