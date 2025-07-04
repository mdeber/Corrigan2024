

# Packages, utility functions ---------------------------------------------

suppressPackageStartupMessages({
    library(Matrix)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(Seurat) # v5.0.3
    library(parallel)
    library(funr) # for finding local paths
})

# for converting from primate to mouse gene symbols
toCapCase <- function(x) {
    .cap <- function(s) paste0(
        toupper(substr(s, 1, 1)),
        tolower(substr(s, 2, nchar(s)))
    )
    unname(sapply(x, .cap))
}

# Hotspot functions -------------------------------------------------------

RunGeneModuleScore <- function(object, features) {
    # Similar to hotspot module score
    # (didn't do their pre-smoothing across neighboring cells)
    # 
    # The module score is PC1 after doing PCA on only the genes in the module;
    # scores are returned with cell barcodes retained as rownames.
    #
    # Correlates very, very highly with seurat's module scores, but way less 
    # convoluted and runs a lot faster, but not compatible with integration
    suppressWarnings({
        object <- ScaleData(object, features = features, verbose = FALSE)
        object <- RunPCA(object, features = features, npcs = 2, verbose = FALSE)
    }) # error if try npcs=1
    pc1 <- FetchData(object, vars = 'PC_1') # rownames are cell barcodes
    
    # I think this is theoretically possible, so:
    if (mean(object@reductions$pca@feature.loadings[,1]) < 0)
        pc1[,1] <- -1*pc1[,1]
    
    # rescale min-to-max score from 0 to 1
    pc1[,1] <- (pc1[,1] - min(pc1[,1])) / (max(pc1[,1]) - min(pc1[,1]))
    setNames(pc1, 'module_score')
}

GeneModuleScorePlot <- function(object, features, ...) {
    object$module_score <- RunGeneModuleScore(object, features)$module_score
    FeaturePlot(object, 'module_score', ..., order = TRUE)
}

hsModuleTableToList <- function(df) {
    # Input: 'run_hotspot.py' outputs a file 'modules.csv'; and the input for 
    # this function is `df <- read.csv('modules.csv')` 
    #
    # Output: A list, elements are named "module_1", "module_2", etc., each of
    # which contains vector of all the genes in that module
    modules <- seq_len(max(as.numeric(df$Module))) # drops module '-1' (no module)
    setNames(
        lapply(modules, function(i) subset(df, Module == i)$X) ,
        paste0("module_", modules)
    )
}

getAllModuleScores <- function(object, module_list, ncores = NCORES) {
    # Inputs: 
    #   seurat object
    #   module_list from hsModuleTableToList()
    # Output: 
    #   dataframe, cell barcodes for rownames, a column for each module
    #   (named 'module_1', 'module_2', etc.)
    out <- mclapply(
        module_list, 
        function(x) RunGeneModuleScore(object, x), 
        mc.cores = ncores
    )
    out <- Map(setNames, out, names(out))
    do.call(cbind, out)
}

getMarkerGeneModuleCorrelations <- function(
        object, all_module_scores, markers, ncores=NCORES
) {
    # Inputs:
    #   object: seurat object
    #   all_module_scores: output of `getAllModuleScores()`
    #   markers: vector of all marker genes of interest
    # Output:
    #   A list with an element for each marker gene. Each element contains a 
    #   vector with length=length(module_list), each element the spearman 
    #   correlation of the marker gene to that module's score
    if (!all(markers %in% rownames(object))) {
        warning("Not all markers found in object. Subsetting markers...")
        markers <- markers[markers %in% rownames(object)]
    }
    dsub <- object[["RNA"]]$data[markers, ] # avoid caching seurat obj
    cor_fun <- function(marker) {
        sapply(
            seq_len(ncol(all_module_scores)), 
            function(i) suppressWarnings(cor(
                all_module_scores[, i], 
                dsub[marker, ], 
                method = "spearman"
            ))
        )
    }
    markers <- setNames(markers, markers)
    mscore_cors_by_marker <- mclapply(markers, cor_fun, mc.cores = ncores)
    return(mscore_cors_by_marker)
}


# Paths -------------------------------------------------------------------

# 
# Note: 
# funr::get_script_path() works when running script from the command line
# (`Rscript <PATH_TO_THIS_SCRIPT>`) or if sourcing in R
# (`source(<PATH_TO_THIS_SCRIPT>)`), but it won't work if you are running this
# line-by-line. In that case, manually substitute your local path to this repo
# 

pdir <- file.path(funr::get_script_path())
path_out <- file.path(path_parent, "reciprocal_hotspot_projections")

pdirs_data <- list(
    marm = file.path(pdir, "marm_census/data/marm_str_gad_drop2_dropGlut_0p8"),
    abc = file.path(pdir, "abc/data/abc_seurat_cl08_str_dropSubc057_0p8"),
    siletti = file.path(pdir, "siletti/data/siletti_neur_roi_sclust_dropAllMeis2_0p8"),
    dev_ferret = file.path(pdir, "dev_ferret/data"),
    dev_mouse = file.path(pdir, "dev_mouse/data"),
    dev_pig = file.path(pdir, "dev_pig/data"),
    dev_macaque = file.path(pdir, "dev_macaque/data")
)
stopifnot(all(sapply(pdirs_data, dir.exists)))

paths_rds <- list(
    marm = file.path(pdirs_data$marm, "marm_str_gad_drop2_dropGlut_0p8.rds"),
    abc = file.path(pdirs_data$abc, "abc_seurat_cl08_str_dropSubc057_0p8.rds"),
    siletti = file.path(pdirs_data$siletti, "siletti_neur_roi_sclust_dropAllMeis2_0p8.rds"),
    dev_ferret = file.path(pdirs_data$dev_ferret, "ferret_processed.rds"),
    dev_mouse = file.path(pdirs_data$dev_mouse, "mouse_processed.rds"),
    dev_pig = file.path(pdirs_data$dev_pig, "pig_processed.rds"),
    dev_macaque = file.path(pdirs_data$dev_macaque, "MacaqueDevInhibitoryNeurons.rds")
)
stopifnot(all(sapply(paths_rds, file.exists)))

paths_hs <- lapply(pdirs_data, file.path, "hotspot")
stopifnot(all(sapply(paths_hs, dir.exists)))

path_df_ortho <- file.path(
    pdir, "ortholog_tables/ens_112_ortho_human_ferret_pig_mouse_marmoset.csv"
)
df_ortho <- read.csv(path_df_ortho)

if (!dir.exists(path_out))
    dir.create(path_out)


# Import Seurat objects ---------------------------------------------------

seur_list <- lapply(paths_rds, readRDS)

# dev macaque data lacks normalization; but actually not relevant in this script
seur_list$dev_macaque <- NormalizeData(seur_list$dev_macaque)

# > seur_list
# $marm
# An object of class Seurat 
# 19762 features across 3379 samples within 1 assay 
# Active assay: RNA (19762 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap
# 
# $abc
# An object of class Seurat 
# 32245 features across 3586 samples within 1 assay 
# Active assay: RNA (32245 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap
# 
# $siletti
# An object of class Seurat 
# 58252 features across 6617 samples within 1 assay 
# Active assay: RNA (58252 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap
# 
# $dev_ferret
# An object of class Seurat 
# 26253 features across 10540 samples within 1 assay 
# Active assay: RNA (26253 features, 0 variable features)
# 2 layers present: counts, data
# 3 dimensional reductions calculated: pca_harmony, pca, umap
# 
# $dev_mouse
# An object of class Seurat 
# 28153 features across 71023 samples within 1 assay 
# Active assay: RNA (28153 features, 0 variable features)
# 2 layers present: counts, data
# 3 dimensional reductions calculated: pca_harmony, pca, umap
# 
# $dev_pig
# An object of class Seurat 
# 30284 features across 8140 samples within 1 assay 
# Active assay: RNA (30284 features, 0 variable features)
# 2 layers present: counts, data
# 3 dimensional reductions calculated: pca_harmony, pca, umap
# 
# $dev_macaque
# An object of class Seurat 
# 58628 features across 109111 samples within 1 assay 
# Active assay: RNA (58628 features, 0 variable features)
# 2 layers present: counts, data
# 2 dimensional reductions calculated: pca, umap


# Import hotspot modules --------------------------------------------------

hs_module_csv_paths <- lapply(paths_hs, file.path, "modules.csv")
hs_module_tables <- lapply(hs_module_csv_paths, read.csv)

# and make into lists
hs_module_lists <- lapply(hs_module_tables, hsModuleTableToList)

# > lengths(hs_module_lists) # number of modules per dataset
# marm      abc      siletti    dev_ferret    dev_mouse    dev_pig  dev_macaque 
#  111      128          155            87           98        106          142 

# 250429, siletti now has 133


# Marker genes for each species -------------------------------------------

#
# note that all of these lists (`seur_list`, `hs_module_lists`, `markers`) are
# fully parallel with one another and match up (and which is why there are two
# elements in `markers` with all the mouse markers, because there are two
# mouse datasets).
# 
# the outputs generated in sections below will also be parallel
# 

markers <- list(
    marm = NULL,
    mouse = NULL,
    human = NULL,
    ferret = NULL,
    mouse_again = NULL, # to make this parallel with the data
    pig = NULL,
    macaque = NULL
)

markers$human <- c(
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

stopifnot(all(markers$human %in% rownames(seur_list$siletti)))

# marmoset data lacks a few of those markers
markers$marm <- markers$human[!markers$human %in% c("MIA", "AVP", "CENPF")]
stopifnot(all(markers$marm %in% rownames(seur_list$marm)))

# mouse only needs to fix Tac2; dev_mouse data also has all the genes
markers$mouse <- toCapCase(markers$human)
markers$mouse <- sub("Tac3", "Tac2", markers$mouse)
stopifnot(all(markers$mouse %in% rownames(seur_list$abc)))
stopifnot(all(markers$mouse %in% rownames(seur_list$dev_mouse)))
markers$mouse_again <- markers$mouse

# genes missing from the ferret data
markers$ferret <- markers$human[
    !markers$human %in% c("PVALB", "TH", "MIA", "NRTN")
]
stopifnot(all(markers$ferret %in% rownames(seur_list$dev_ferret)))

# genes missing from the pig data
markers$pig <- markers$human[!markers$human == "CHRNB4"]
stopifnot(all(markers$pig %in% rownames(seur_list$dev_pig)))

# genes missing from macaque data
markers$macaque <- markers$human[!markers$human == "MIA"]


# Calculate all module scores ---------------------------------------------

# takes around 10 minutes; could be more optimized
all_module_scores <- Map(
    getAllModuleScores,
    object      = seur_list,
    module_list = hs_module_lists,
    ncores      = list(1)
)


# Calculate all marker gene correlations ----------------------------------

# for every marker gene, its correlations with all hotspot modules
# (see function documentation above)

# also takes on the order of 10 minutes... also could be more optimized
all_marker_module_correlations <- Map(
    getMarkerGeneModuleCorrelations,
    object            = seur_list,
    all_module_scores = all_module_scores,
    markers           = markers,
    ncores            = list(1)
)


# Get most correlated modules for each marker -----------------------------

# list (of species) with lists (of markers) of integers (max module)
max_marker_module_correlations <- lapply(
    all_marker_module_correlations, lapply, which.max
)


# Write out those max correlated modules ----------------------------------

max_cor_dir <- file.path(path_out, "most_cor_modules_for_markers")
if (!dir.exists(max_cor_dir))
    dir.create(max_cor_dir)

max_marker_module_correlations %>% 
    Map(function(x, nm) {
        # some missing values in pig data
        x[lengths(x) == 0] <- NA
        df <- data.frame(
            gene = names(x),
            max_cor_module = unlist(x),
            row.names = NULL
        )
        write.csv(
            df, 
            file.path(max_cor_dir, paste0(nm, ".csv")),
            row.names = FALSE
        )
    }, ., names(.)) %>% 
    invisible


# Reciprocal module projection function -----------------------------------

# this is the function that takes modules from one species, and makes the 
# module score plots in the other species as well. 
# 
# this function is specifically written to make use of the objects as they are
# made in this script
checkModulesAcrossSpecies <- function(
        ref_species, # options are elements of names(seur_list)
        gene,
        seur_list = seur_list,
        hs_module_lists = hs_module_lists,
        max_marker_module_correlations = max_marker_module_correlations
) {
    max_mod <- max_marker_module_correlations[[ref_species]][[gene]]
    ref_mod <- hs_module_lists[[ref_species]][[max_mod]] # genes in the module
    
    # want to plot in a certain order depending on the reference species;
    # just popping out the name of the ref_species and putting it first:
    species_list <- unique(c(ref_species, names(seur_list)))
    
    # gene symbols to match up with species_list 
    # (this should work regardless of whether mouse or non-mouse gene names
    # are used)
    is_species_mouse <- species_list %in% c("abc", "dev_mouse")
    mk_gene_symbols <- ifelse(
        is_species_mouse,
        toCapCase(gene), 
        toupper(gene)
    )
    mk_gene_symbols <- sub("TAC2", "TAC3", mk_gene_symbols)
    mk_gene_symbols <- sub("Tac3", "Tac2", mk_gene_symbols)
    
    # for plot titles
    proper_names <- sapply(
        species_list,
        switch,
        siletti = "Human (Siletti 2023)", 
        marm = "Marmoset (Krienen 2023)", 
        abc = "Mouse (Yao 2023)",
        dev_ferret = "Developing ferret",
        dev_mouse = "Developing mouse",
        dev_pig = "Developing pig",
        dev_macaque = "Developing macaque"
    )
    
    # [For the marker genes, this info was already saved to csv files]
    cat(paste0(
        "For ", proper_names[1], " gene ", gene, ", using module ", max_mod, "\n"
    ))
    
    # get gene symbols for all genes in module, for all species
    module_genes <- rep(list(ref_mod), length(species_list))
    names(module_genes) <- species_list
    module_genes[is_species_mouse] <- lapply(
        module_genes[is_species_mouse], toCapCase
    )
    module_genes[!is_species_mouse] <- lapply(
        module_genes[!is_species_mouse], toupper
    )
    module_genes <- lapply(module_genes, function(x) {
        sub("TAC2", "TAC3", x)
        sub("Tac3", "Tac2", x)
        x
    })
    
    # subset for any missing genes
    module_genes <- Map(
        function(genes, object) genes[genes %in% rownames(object)],
        module_genes,
        seur_list[species_list]
    )
    
    plots <- Map(
        function(species, gene, proper_name, module) {
            if (gene %in% rownames(seur_list[[species]])) {
                p1 <- FeaturePlot(seur_list[[species]], gene, order = TRUE) +
                    labs(title = paste0(proper_name, "\n\n", gene))
            } else {
                p1 <- ggplot() + 
                    geom_blank() + 
                    annotate(
                        geom = "text", 
                        label = "Gene not in data", 
                        x = 0, 
                        y = 0
                    ) +
                    labs(title = paste0(proper_name, "\n\n", gene)) + 
                    cowplot::theme_cowplot() + # used by Seurat::FeaturePlot
                    theme(
                        plot.title = element_text(hjust = 0.5),
                        axis.line = element_blank(),
                        axis.text = element_blank(),
                        axis.title = element_blank(),
                        axis.ticks = element_blank()
                    )
            }
            p2 <- GeneModuleScorePlot(seur_list[[species]], module)
            p1 / p2
        },
        species_list,
        mk_gene_symbols,
        proper_names,
        module_genes
    )
    do.call(wrap_plots, c(plots, nrow = 1))
}


# Reciprocal module projections -------------------------------------------


### For adult datasets ---

reciprocal_module_dir <- file.path(path_out, "reciprocal_projections_3species")
if (!dir.exists(reciprocal_module_dir))
    dir.create(reciprocal_module_dir)

adult_idx <- 1:3

mclapply(names(seur_list)[adult_idx], function(ref) {
    save_dir <- file.path(reciprocal_module_dir, paste0("modules_from_", ref))
    if (!dir.exists(save_dir))
        dir.create(save_dir)
    
    list_idx <- which(names(seur_list) == ref)
    for (gn in markers[[list_idx]]) {
        
        outpath <- file.path(save_dir, paste0(gn, ".pdf"))
        if (file.exists(outpath))
            next
        
        # necessary for the pig
        if (length(max_marker_module_correlations[[ref]][[gn]]) == 0)
            next
        
        # Error on pig "GFRA2" (pig module_47); module lacks genes with real names
        p <- try(checkModulesAcrossSpecies(
            ref_species = ref,
            gene = gn,
            seur_list = seur_list[adult_idx],
            hs_module_lists = hs_module_lists[adult_idx],
            max_marker_module_correlations = max_marker_module_correlations[adult_idx]
        ), silent = TRUE)
        
        if (any(class(p) == "patchwork")) {
            ggsave(outpath, p, height = 220, width = 440, units = "mm")
        }
    }
}, mc.cores = 3) %>% 
    invisible


### For all modules in all datasets ---

reciprocal_module_dir <- file.path(path_out, "reciprocal_projections_6species")
if (!dir.exists(reciprocal_module_dir))
    dir.create(reciprocal_module_dir)

mclapply(names(seur_list), function(ref) {
    save_dir <- file.path(reciprocal_module_dir, paste0("modules_from_", ref))
    if (!dir.exists(save_dir))
        dir.create(save_dir)
    
    list_idx <- which(names(seur_list) == ref)
    for (gn in markers[[list_idx]]) {
        
        outpath <- file.path(save_dir, paste0(gn, ".pdf"))
        if (file.exists(outpath))
            next
        
        # necessary for the pig
        if (length(max_marker_module_correlations[[ref]][[gn]]) == 0)
            next
        
        # Error on pig "GFRA2" (pig module_47); module lacks genes with real names
        p <- try(checkModulesAcrossSpecies(
            ref_species = ref,
            gene = gn,
            seur_list = seur_list,
            hs_module_lists = hs_module_lists,
            max_marker_module_correlations = max_marker_module_correlations
        ), silent = TRUE)
        
        if (any(class(p) == "patchwork")) {
            ggsave(outpath, p, height = 220, width = 440, units = "mm")
        }
    }
}, mc.cores = 2) %>% 
    invisible


### For all modules in developing datasets ---

reciprocal_module_dir <- file.path(path_out, "reciprocal_projections_dev_data")
if (!dir.exists(reciprocal_module_dir))
    dir.create(reciprocal_module_dir)

dev_idx <- 4:7
# > names(seur_list)[dev_idx]
# [1] "dev_ferret"  "dev_mouse"  "dev_pig"   "dev_macaque" 

mclapply(names(seur_list)[dev_idx], function(ref) {
    save_dir <- file.path(reciprocal_module_dir, paste0("modules_from_", ref))
    if (!dir.exists(save_dir))
        dir.create(save_dir)
    
    list_idx <- which(names(seur_list) == ref)
    for (gn in markers[[list_idx]]) {
        
        outpath <- file.path(save_dir, paste0(gn, ".pdf"))
        if (file.exists(outpath))
            next
        
        # necessary for the pig
        if (length(max_marker_module_correlations[[ref]][[gn]]) == 0)
            next
        
        # Error on pig "GFRA2" (pig module_47); module lacks genes with real names
        p <- try(checkModulesAcrossSpecies(
            ref_species = ref,
            gene = gn,
            seur_list = seur_list[dev_idx],
            hs_module_lists = hs_module_lists[dev_idx],
            max_marker_module_correlations = max_marker_module_correlations[dev_idx]
        ), silent = TRUE)
        
        if (any(class(p) == "patchwork")) {
            ggsave(outpath, p, height = 220, width = 440, units = "mm")
        }
    }
}, mc.cores = 2) %>% 
    invisible

# (Export genes by modules for dev datasets) ------------------------------

# Wasn't done in any previous scripts for these datasets

lapply(paths_hs[dev_idx], function(path) {
    if (!dir.exists(file.path(path, "genes_by_module")))
        dir.create(file.path(path, "genes_by_module"))
    NULL
})

Map(
    function(mod_list, path) {
        mod_list %>% 
            Map(function(x, nm) {
                out <- file.path(path, "genes_by_module", paste0(nm, ".txt"))
                write(x, out)
            }, ., names(.)) %>% 
            invisible()
    },
    hs_module_lists[dev_idx],
    paths_hs[dev_idx]
) %>% invisible()

