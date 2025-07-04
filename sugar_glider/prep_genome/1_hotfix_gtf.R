# Error message on cellranger mkfref:
# ===
# mkref has FAILED
# Error building reference package
# Invalid gene annotation input: in GTF
# /mnt/cup/labs/krienen/genomic_annotations/cellranger_ref/Pbrev_HiC/Pbrev_HiC/genes/genes.gtf.gz
# records for gene_id gene-LOC118843072 are not contiguous in the file
# ===
#
# The issue is that there are individual genes with multiple transcripts, and
# those transcripts don't appear next to each other in the gtf.
#
# The solution is to reorganize the lines in the gtf such that all elements for
# a single gene are next to each other.
#
# The code below works because each transcript for a given gene is organized
# correctly e.g., a set of 3 lines might be "transcript"-"exon"-"CDS". It's just 
# that there are multiple transcripts for some genes that don't appear next to
# each other.
# 
suppressPackageStartupMessages({
    library(rtracklayer)
})
path_gtf <- '/PATH/TO/Pbrev_HiC/Tvul_vs_Pbrev_a05s05d5_withmito.gtf'
path_out <- '/PATH/TO/Pbrev_HiC/Tvul_vs_Pbrev_a05s05d5_withmito_hotfix.gtf'
gtf <- import(path_gtf)
gr  <- unlist(split(gtf, gtf$gene_id))
names(gr) <- NULL
# export(gr, path_out) # not so fast...

# ===
# Problem number 2: Some genes have transcripts on multiple strands. This fails
# when cellranger runs.
# ===

strands <- split(as.character(strand(gr)), gr$gene_id)
strands <- lapply(strands, unique)
gn_fix  <- names(strands)[lengths(strands) > 1]

# > gn_fix
# [1] "mt-gene-tRNA-Ser"

# Just one gene

idx <- which(gr$gene_id == gn_fix)

# > gr[idx]
# GRanges object with 4 ranges and 7 metadata columns:
#     seqnames        ranges strand |   source       type     score
# <Rle>     <IRanges>  <Rle> | <factor>   <factor> <numeric>
# [1] HiC_scaffold_16348 141970-142027      - |   RefSeq transcript        NA
# [2] HiC_scaffold_16348 141970-142027      - |   RefSeq exon              NA
# [3] HiC_scaffold_16348 146694-146762      + |   RefSeq transcript        NA
# [4] HiC_scaffold_16348 146694-146762      + |   RefSeq exon              NA
# 
# phase   transcript_id          gene_id   gene_name
# <integer>     <character>      <character> <character>
# [1]      <NA> mt-rna-tRNA-Ser mt-gene-tRNA-Ser mt-tRNA-Ser
# [2]      <NA> mt-rna-tRNA-Ser mt-gene-tRNA-Ser mt-tRNA-Ser
# [3]      <NA> mt-rna-tRNA-Ser mt-gene-tRNA-Ser mt-tRNA-Ser
# [4]      <NA> mt-rna-tRNA-Ser mt-gene-tRNA-Ser mt-tRNA-Ser

# Will just append a 1 and 2 to the two different transcripts
gr$transcript_id[idx[1:2]] <- paste0(gr$transcript_id[idx[1:2]], "1")
gr$gene_id[idx[1:2]]       <- paste0(gr$gene_id[idx[1:2]], "1")
gr$gene_name[idx[1:2]]     <- paste0(gr$gene_name[idx[1:2]], "1")

gr$transcript_id[idx[3:4]] <- paste0(gr$transcript_id[idx[3:4]], "2")
gr$gene_id[idx[3:4]]       <- paste0(gr$gene_id[idx[3:4]], "2")
gr$gene_name[idx[3:4]]     <- paste0(gr$gene_name[idx[3:4]], "2")

# export(gr, path_out)


# ===
# Problem 3: The mitochondrial transcripts don't have exons. From cellranger:
# 
# Transcripts with no exons were detected in the GTF.
# Please ensure there is at least one `exon` entry for each `transcript` in the GTF.
# The following transcripts are lacking an exon entry:
# Transcript ID: mt-gene-ATP6, Gene Name: mt-ATP6
# Transcript ID: mt-gene-ATP8, Gene Name: mt-ATP8
# Transcript ID: mt-gene-COX1, Gene Name: mt-COX1
# Transcript ID: mt-gene-COX2, Gene Name: mt-COX2
# Transcript ID: mt-gene-COX3, Gene Name: mt-COX3
# Transcript ID: mt-gene-CYTB, Gene Name: mt-CYTB
# Transcript ID: mt-gene-ND1, Gene Name: mt-ND1
# Transcript ID: mt-gene-ND2, Gene Name: mt-ND2
# Transcript ID: mt-gene-ND3, Gene Name: mt-ND3
# Transcript ID: mt-gene-ND4, Gene Name: mt-ND4
# Transcript ID: mt-gene-ND4L, Gene Name: mt-ND4L
# Transcript ID: mt-gene-ND5, Gene Name: mt-ND5
# Transcript ID: mt-gene-ND6, Gene Name: mt-ND6
# ===

grl <- split(gr, gr$transcript_id)
chk <- sapply(grl, function(x) "exon" %in% x$type)

# Verify the cellranger error
# 
# > names(grl)[!chk]
# [1] "mt-gene-ATP6" "mt-gene-ATP8" "mt-gene-COX1" "mt-gene-COX2" "mt-gene-COX3"
# [6] "mt-gene-CYTB" "mt-gene-ND1"  "mt-gene-ND2"  "mt-gene-ND3"  "mt-gene-ND4" 
# [11] "mt-gene-ND4L" "mt-gene-ND5"  "mt-gene-ND6"
# 

# They all have a CDS but lack an exon:
# 
# > unique(lapply(grl[!chk], function(x) as.character(x$type)))
# [[1]]
# [1] "transcript" "CDS"

# So will copy the CDS and call that the exon.

grl[!chk] <- endoapply(grl[!chk], function(x) {
    stopifnot(all(x$type == c("transcript", "CDS")))
    x[3] <- x[2]
    x$type[2] <- "exon"
    x$phase[2] <- NA
    x
})
gr <- unlist(grl)
names(gr) <- NULL

# and make sure that the genes are still contiguous before exporting
gr <- unlist(split(gr, gr$gene_id))
names(gr) <- NULL
export(gr, path_out)

