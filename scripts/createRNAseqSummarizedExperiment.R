#'#################################################################################
#'#################################################################################
#' Create RangedSummarizedExperiment
#' This code generates a SummarizedExperiment from a matrix of RNAseq counts.
#' It can be reused for other projects but modifying some parts, highlighted with ##
#' Things to take into account:
#' - Genes annotation is loaded from gff file. 
#' - Counts data is expected to come from gene quantification.
#' - Sample names is counts is corrected to match pheno names.
#' - Keep genes with more than 0.5 cpm in more than 5 samples 
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
countFile <- args[1]
annotFile <- args[2]
phenoFile <- args[3]

## Load libraries ####
library(SummarizedExperiment)
library(edgeR)

## Load data
counts <- read.table(countFile, header = TRUE, row.names = 1)
annot <- read.delim(gzfile(annotFile), comment.char = "#", header = FALSE, as.is = TRUE)
load(phenoFile)

## Adapt data
### Remove rows corresponding to unmapped reads
counts <- counts[-grep("N_", rownames(counts)), ]

#'#################################################################################
#'#################################################################################
## Map sample names to pheno ids
colnames(counts) <- gsub(".", "", fixed = TRUE, colnames(counts))
colnames(counts) <- gsub("X", "", fixed = TRUE, colnames(counts))

#'#################################################################################
#'#################################################################################

## Convert counts to matrix
counts <- data.matrix(counts)

### Add colnames to annotation
colnames(annot) <- c("seqname", "source", "feature", "start", "end", "score", 
                     "strand", "frame", "attribute")

### Select only exons
gene_annot <- subset(annot, feature == "gene")

### Add gene name column
gene_annot$geneID <- gsub("gene_id ", "", sapply(strsplit(gene_annot$attribute, ";"), `[`, 1))
rownames(gene_annot) <- gene_annot$geneID

## Create GRanges
annotGR <- makeGRangesFromDataFrame(gene_annot, keep.extra.columns = TRUE)
names(annotGR) <- annotGR$geneID

## Change row.names pheno to sampleID
rownames(pheno) <- pheno$SampleID

### Create SummarizedExperiment
rse <- SummarizedExperiment(assays = list(counts = counts), 
                            rowRanges = annotGR[rownames(counts)], 
                            colData = pheno[colnames(counts), ])
save(rse, file = "RNAseq_RangedSE_allGenes.Rdata")

### Filter genes with low counts (criterium in header)
rse <- rse[rowSums(cpm(assay(rse)) > 0.5) > 5, ]
save(rse, file = "RNAseq_RangedSE_expressedGenes.Rdata")

### Filter genes with low counts (criterium in header)
rse <- rse[!seqnames(rowRanges(rse)) %in% c("chrX", "chrY"), ]
save(rse, file = "RNAseq_RangedSE_autosomicGenes.Rdata")
