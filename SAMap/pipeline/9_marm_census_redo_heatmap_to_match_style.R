# Quick match SAMap output heatmap to correlation analysis heatmap

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(Seurat) # v5.0.3
    library(ggpubr)
    library(ComplexHeatmap)
    library(viridis)
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

# 
# Note: 
# funr::get_script_path() works when running script from the command line
# (`Rscript <PATH_TO_THIS_SCRIPT>`) or if sourcing in R
# (`source(<PATH_TO_THIS_SCRIPT>)`), but it won't work if you are running this
# line-by-line. In that case, manually substitute your local path
# 
dir_samap     <- dirname(funr::get_script_path())
dir_samap_res <- file.path(dir_samap, "results_marm_census")

samap_res <- read.csv(
    file.path(dir_samap_res, "mapping_table_clust_anno_splitsubc55.csv")
)
rownames(samap_res) <- samap_res[, 1]
samap_res <- as.matrix(samap_res[, -1])

x <- rownames(samap_res)
x <- sub("hg_", "human.", x)
x <- sub("cj_", "marmoset.", x)
x <- sub("mm_", "mouse.", x)
x <- ifelse(
    grepl("mouse.08", x),
    sub("mouse.", "mouse.ABC_", x),
    x
)
x <- sub(" .*", "", x)
x <- sub("mouse.Th", "mouse.Th (non sc55)", x)

rownames(samap_res) <- colnames(samap_res) <- x

# potentially drop:
# "marmoset.LAMP5/LHX6"


# Fix diagonal ------------------------------------------------------------

samap_res_fix <- samap_res
stopifnot(identical(rownames(samap_res), colnames(samap_res)))
for (i in 1:nrow(samap_res_fix)) {
    samap_res_fix[i, i] <- 1
}


# All clusters ------------------------------------------------------------

row_groups <- ifelse(grepl("TAC3", rownames(samap_res_fix)), "tac3", "other")
row_groups <- factor(row_groups, levels = c("other", "tac3"))

col_fun <- circlize::colorRamp2(seq(0, 1, length.out = 100), viridis(100))

p <- Heatmap(
    samap_res_fix,
    heatmap_legend_param = list(title = "SAMap"),
    row_split = row_groups,
    column_split = row_groups,
    col = col_fun
)


pdf(
    file.path(dir_samap_res, "mapping_table_clust_anno_splitsubc55_complexheatmap.pdf"),
    width = 9, height = 7
)
p
dev.off()


# Fix ordering ------------------------------------------------------------

# Define group labels
row_groups <- ifelse(grepl("TAC3", rownames(samap_res_fix)), "tac3", "other")
col_groups <- ifelse(grepl("TAC3", colnames(samap_res_fix)), "tac3", "other")

# Set the desired order
desired_group_order <- c("other", "tac3")  # will be shown top to bottom / left to right

# Convert to factors with correct order
row_groups <- factor(row_groups, levels = desired_group_order)
col_groups <- factor(col_groups, levels = desired_group_order)

# Get custom row/column order by clustering within each group
get_clustered_order <- function(mat, axis = c("row", "column")) {
    axis <- match.arg(axis)
    split_labels <- if (axis == "row") row_groups else col_groups
    full_order <- c()
    
    for (grp in levels(split_labels)) {
        idx <- which(split_labels == grp)
        submat <- if (axis == "row") mat[idx, , drop = FALSE] else mat[, idx, drop = FALSE]
        hc <- hclust(dist(if (axis == "row") submat else t(submat)))
        ordered_idx <- idx[hc$order]
        full_order <- c(full_order, ordered_idx)
    }
    
    return(full_order)
}

custom_row_order <- get_clustered_order(samap_res_fix, "row")
custom_col_order <- get_clustered_order(samap_res_fix, "column")

# Define color scale
col_fun <- circlize::colorRamp2(seq(0, 1, length.out = 100), viridis(100))

# Final heatmap with full control over split order
p <- Heatmap(
    samap_res_fix,
    name = "SAMap",
    col = col_fun,
    row_split = row_groups,
    column_split = col_groups,
    row_order = custom_row_order,
    column_order = custom_col_order,
    # show_row_names = FALSE,
    # show_column_names = FALSE,
    heatmap_legend_param = list(title = "SAMap")
)

pdf(
    file.path(dir_samap_res, "mapping_table_clust_anno_splitsubc55_complexheatmap_reorder_nodend.pdf"),
    width = 9, height = 7
)
p
dev.off()


# drop Th non-sc55 --------------------------------------------------------

stopifnot(all(rownames(samap_res_fix) == colnames(samap_res_fix)))
idx_drop <- which(rownames(samap_res_fix) == "mouse.Th (non sc55)")
samap_res_fix_drop <- samap_res_fix[-idx_drop, -idx_drop]
stopifnot(all(rownames(samap_res_fix_drop) == colnames(samap_res_fix_drop)))


row_groups <- ifelse(grepl("TAC3", rownames(samap_res_fix_drop)), "tac3", "other")
row_groups <- factor(row_groups, levels = c("other", "tac3"))

col_fun <- circlize::colorRamp2(seq(0, 1, length.out = 100), viridis(100))

p <- Heatmap(
    samap_res_fix_drop,
    heatmap_legend_param = list(title = "SAMap"),
    row_split = row_groups,
    column_split = row_groups,
    col = col_fun
)

pdf(
    file.path(dir_samap_res, "mapping_table_clust_anno_splitsubc55_dropOtherTh_complexheatmap.pdf"),
    width = 9, height = 7
)
p
dev.off()


# Fix ordering ------------------------------------------------------------

# Define group labels
row_groups <- ifelse(grepl("TAC3", rownames(samap_res_fix_drop)), "tac3", "other")
col_groups <- ifelse(grepl("TAC3", colnames(samap_res_fix_drop)), "tac3", "other")

# Set the desired order
desired_group_order <- c("other", "tac3")  # will be shown top to bottom / left to right

# Convert to factors with correct order
row_groups <- factor(row_groups, levels = desired_group_order)
col_groups <- factor(col_groups, levels = desired_group_order)

# Get custom row/column order by clustering within each group
get_clustered_order <- function(mat, axis = c("row", "column")) {
    axis <- match.arg(axis)
    split_labels <- if (axis == "row") row_groups else col_groups
    full_order <- c()
    
    for (grp in levels(split_labels)) {
        idx <- which(split_labels == grp)
        submat <- if (axis == "row") mat[idx, , drop = FALSE] else mat[, idx, drop = FALSE]
        hc <- hclust(dist(if (axis == "row") submat else t(submat)))
        ordered_idx <- idx[hc$order]
        full_order <- c(full_order, ordered_idx)
    }
    
    return(full_order)
}

custom_row_order <- get_clustered_order(samap_res_fix_drop, "row")
custom_col_order <- get_clustered_order(samap_res_fix_drop, "column")

# Define color scale
col_fun <- circlize::colorRamp2(seq(0, 1, length.out = 100), viridis(100))

# Final heatmap with full control over split order
p <- Heatmap(
    samap_res_fix_drop,
    name = "SAMap",
    col = col_fun,
    row_split = row_groups,
    column_split = col_groups,
    row_order = custom_row_order,
    column_order = custom_col_order,
    # show_row_names = FALSE,
    # show_column_names = FALSE,
    heatmap_legend_param = list(title = "SAMap")
)


pdf(
    file.path(dir_samap_res, "mapping_table_clust_anno_splitsubc55_dropOtherTh_complexheatmap_reorder_nodend.pdf"),
    width = 9, height = 7
)
p
dev.off()
