# Mike DeBerardine
# Princeton Neuroscience Institute
# February 2025
#
# [!] Throughout, this script expects a `gene_id` column in the GTF file. If not
#     present, you should edit that in, or replace `gene_id` everywhere in this
#     script. 
#
# [!] There is one more place at the end of this script where all the fields
#     much match what is in your GTF. If your GTF doesn't match that, you either
#     need to edit your GTF or that section to ensure the new elements being
#     added fit with your annotations.
# 
# [!] Chrom sizes are UCSC-style files. On fa.fai index file, run to generate:
#     `cut -f1,2 input.fa.fai > chrom.sizes`.
#     The output is 2 columns, tab-separated, no header, first column is chrom
#     name, second is the length of the chromosome.
# 
path_gtf_in  <- '/PATH/TO/Pbrev_HiC/Tvul_vs_Pbrev_a05s05d5_withmito_hotfix.gtf'
path_gtf_out <- '/PATH/TO/Pbrev_HiC/Tvul_vs_Pbrev_a05s05d5_withmito_hotfix_3kbExt.gtf'
path_chrom_sizes <- '/PATH/TO/Pbrev_HiC/Pbrev_HiC_chrom.sizes'
ext_length <- 3e3    # the length of the 3' extension to be applied to all genes
min_ext_length <- 20 # the minimum length of a clipped 3' extension to be kept


# Begin -------------------------------------------------------------------

suppressPackageStartupMessages({
    library(dplyr)
    library(rtracklayer)
})

cat("Importing, checking GTF...\n")

gr <- import(path_gtf_in)

if (!"gene_id" %in% names(mcols(gr)))
    stop("GTF file missing the `gene_id` field.")

# double-check no genes have elements on both strands
if (any(duplicated(unique(data.frame(gr$gene_id, strand(gr)))$gene_id)))
    stop("Some genes have elements on both strands. This must be fixed.")


# Get chromosome lengths --------------------------------------------------

cat("Importing, adding chromosome sizes...\n")
chrom_sizes <- read.csv(path_chrom_sizes, sep = "\t", header = FALSE)
chrom_sizes <- setNames(chrom_sizes[[2]], chrom_sizes[[1]])
seqlengths(gr) <- chrom_sizes[names(seqlengths(gr))]


# Create initial 3' end extensions ----------------------------------------

cat("Generating 3' end extension segments...\n")

# separate transcript elements by their genes
grl <- split(gr, gr$gene_id)

# find the 3' most position of each gene
gene_ends <- unlist(resize(range(grl), 1, fix = "end"))
gene_ends$gene_id <- names(gene_ends)

# shift them downstream by 1bp (so they don't overlap the existing annotation)
gene_ends <- suppressWarnings(
    shift(gene_ends, ifelse(strand(gene_ends) == "+", 1, -1))
)

# get segments for the 3' end extensions (will be trimmed if they exceed the chromosome)
gene_ext <- trim(suppressWarnings(
    promoters(gene_ends, 0, ext_length)
))


# Find conflicts ----------------------------------------------------------

# Reduce-by-gene: for each gene annotation, a set of non-overlapping segments
# that cover all its annotations. (Note that doesn't mean different genes won't
# overlap one another, and it's possible in rare cases for a gene to end up with
# multiple non-contiguous segments).
reduced_genes <- unlist(reduce(grl))
reduced_genes$gene_id <- names(reduced_genes)
names(reduced_genes) <- NULL

cat("Finding conflicts...\n")

# Find all cases where a 3' extension overlaps another gene
hits <- findOverlaps(gene_ext, reduced_genes)
conflicts <- data.frame(
    extended_gene = gene_ext$gene_id[queryHits(hits)],
    downstream_gene = reduced_genes$gene_id[subjectHits(hits)]
)

cat("Out of a total number of genes:\n")
cat(paste0(length(unique(gr$gene_id)), "\n"))

cat("Total number of conflicts (intersections between 3' extensions and downstream genes:\n")
cat(paste0(nrow(conflicts), "\n"))

cat("Unique genes whose 3' end extensions are affected:\n")
cat(paste0(length(unique(conflicts$extended_gene)), "\n"))

cat("Unique downstream genes affected::\n")
cat(paste0(length(unique(conflicts$downstream_gene)), "\n"))

cat("Total unique genes implicated:\n")
cat(paste0(length(unique(c(conflicts$extended_gene, conflicts$downstream_gene))), "\n"))

cat("Separately, number of genes whose extensions are clipped by the end of a chromosome:\n")
cat(paste0(sum(width(gene_ext) < ext_length), "\n"))


# Clip or remove conflicting 3' end extensions ----------------------------

cat("\nClipping 3' extensions for those affected genes...\n")

# 
# For all cases where a 3' extension overlaps another gene, determine the 
# maximum extent (here, the genomic position of that maximum extent) of a 3' 
# extension that won't overlap another gene.
# 
# Note: this code is written to be robust to multiple overlapping conflicts, but
# it is still possible (but very uncommon) for the 3' extensions themselves to
# to overlap one another. This isn't addressed here and will just be left to the
# aligner e.g., in cellranger/STAR, reads in those overlaps won't be counted.
# 
max_extension_site <- split(conflicts$downstream_gene, conflicts$extended_gene) %>% 
    lapply(function(x) subset(reduced_genes, gene_id %in% x)) %>% 
    as("GRangesList") %>% 
    range %>% 
    resize(1)

# For the affected extensions, make new ones. 
# First, try to bisect each extension into keepable-vs-unkeepable segments.
gene_ext_clipped <- as(gene_ext, "GRangesList")[names(max_extension_site)]
gene_ext_clipped <- disjoin(pc(gene_ext_clipped, max_extension_site))

# If we were able to bisect the the extension, there will be 3 ranges for each.
# It not, there are only 2. So drop any that we were unable to bisect.
extensions_drop <- names(gene_ext_clipped)[lengths(gene_ext_clipped) < 3]
gene_ext_clipped <- gene_ext_clipped[lengths(gene_ext_clipped) == 3]

cat("3' extensions that fully overlap another gene (likely a pre-existing gene conflict):\n")
cat(paste0(length(extensions_drop), "\n"))

# After bisecting the 3' extensions, take the upstream one
gene_ext_clipped <- endoapply(gene_ext_clipped, function(x) {
    if (unique(strand(x)) == "+") x[1] else x[3]
})
gene_ext_clipped <- unlist(gene_ext_clipped)
gene_ext_clipped$gene_id <- names(gene_ext_clipped)

# Also drop any genes whose extension would be shorter than the minimum 
# extension length.
#
# For this, also want to include extensions that may be off the chromosome
gene_ext[names(gene_ext_clipped)] <- gene_ext_clipped
idx_short <- which(width(gene_ext) < min_ext_length)
extensions_drop <- union(extensions_drop, names(gene_ext)[idx_short])
gene_ext <- gene_ext[!names(gene_ext) %in% extensions_drop]

cat("3' extensions dropped for being less than `min_ext_length` after clipping:\n")
cat(paste0(length(idx_short), "\n"))

cat("Total number of genes not getting extensions:\n")
cat(paste0(length(extensions_drop), "\n"))

cat("Total number of genes that are getting clipped extensions:\n")
cat(paste0(sum(width(gene_ext) < ext_length), "\n"))

cat("Total number of genes that are getting extensions:\n")
cat(paste0(length(gene_ext), "\n"))


# Create new transcripts of the 3' extensions -----------------------------

# Each gene with an extension will get a new transcript (w/ exon) representing
# that extension

#
# [!] This is the section where the fields must match your GTF. 
# 
# -> This code is based on a GTF file for which a gene entry looks like this. 
#    If yours is different, match that! (Note the example gene has one 
#    transcript, which is coding, others might have many transcripts which may
#    not be coding. We're not annotating CDS).
#
# > grl[[2]]
# GRanges object with 3 ranges and 7 metadata columns:
#           seqnames            ranges strand |   source       type     score
#              <Rle>         <IRanges>  <Rle> | <factor>   <factor> <numeric>
# [1] HiC_scaffold_3 23742352-23743487      - |  Liftoff transcript        NA
# [2] HiC_scaffold_3 23742352-23743487      - |  Liftoff exon              NA
# [3] HiC_scaffold_3 23742420-23743487      - |  Liftoff CDS               NA
# 
#         phase      transcript_id     gene_id   gene_name
#     <integer>        <character> <character> <character>
# [1]      <NA> rna-XM_036762030.1 gene-A4GALT      A4GALT
# [2]      <NA> rna-XM_036762030.1 gene-A4GALT      A4GALT
# [3]         0 rna-XM_036762030.1 gene-A4GALT      A4GALT
# -------
#     seqinfo: 474 sequences from an unspecified genome; no seqlengths

cat("Adding in 3' extensions...\n")

# get a template for making the new annotations, start with existing transcript
# annotations, then make them match up with the `gene_ext`
tx_ext <- subset(gr, type == "transcript")
names(tx_ext) <- tx_ext$gene_id
tx_ext <- tx_ext[!duplicated(names(tx_ext))]
tx_ext <- tx_ext[names(gene_ext)] 

# now start editing the specific metadata fields to match the GTF
tx_ext$source <- "ThreePrimeExtension"
tx_ext$score  <- NA
tx_ext$phase  <- NA
tx_ext$transcript_id <- paste0(tx_ext$gene_id, "_Extension")

# now apply that metadata to the final extensions; but then we want to follow
# that by making new annotations for a "transcript" and an "exon" for each 
# extension
mcols(gene_ext) <- mcols(tx_ext)
names(gene_ext) <- NULL 
tx_ext <- gene_ext
ex_ext <- gene_ext
ex_ext$type <- "exon"
to_add <- c(tx_ext, ex_ext)
to_add <- split(to_add, to_add$gene_id)

# now add those to the transcripts
grl[names(to_add)] <- pc(grl[names(to_add)], unname(to_add))


# Export modified GTF -----------------------------------------------------

cat("Exporting modified GTF...\n")
out <- unname(unlist(grl))
export(out, path_gtf_out)

cat("Script finished.\n")
