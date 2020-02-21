#'#################################################################################
#'#################################################################################
#' Create GenomicRatioSets
#' This code prepares a GenomicRatioSet for analysis, by removing different types of
#' probes:
#' - Probes not measuring methylation (SNPs, CH probes)
#' - Crosshibridizing probes
#' - Probes with SNPs
#' - Probes in sexual chromosomes
#' 
#' For crosshibridizing probes and probes with SNPs, we used annotation from PMID: 27924034
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gsetfile <- args[1]
manifest <- args[2]

## Load libraries ####
library(minfi)

## Load dataset ####
load(gsetfile)
grAnnot <- readRDS(manifest)

## Remove technical duplicate
gset <- gset[, colnames(gset) != "G2CS12COBIS"]
save(gset, file = "gset.allProbes.Rdata")

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = "gset.allCpGs.Rdata")

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file = "gset.filterAnnotatedProbes.Rdata")

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = "gset.autosomic.Rdata")
