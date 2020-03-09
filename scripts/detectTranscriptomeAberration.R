#'#################################################################################
#'#################################################################################
#' Detect transcriptome aberrations
#' This code detects transcriptome aberrations with OUTRIDER
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
rseFile <- args[1]
gtfFile <- args[2]

## Load libraries ####
library(OUTRIDER)

load(rseFile)

ods <- OutriderDataSet(rse)

## Filter non-expressed genes
ods.filt <- filterExpression(ods, gtfFile = gtfFile, addExpressedGenes=TRUE)

ods.cor <- estimateSizeFactors(ods.filt)
ods.cor <- findEncodingDim(ods.cor)
ods.cor <- controlForConfounders(ods.cor)

# run full OUTRIDER pipeline
ods.res <- OUTRIDER(ods.cor)
save(ods.cor, file = "OUTRIDER_QC.Rdata")
save(ods.res, file = "OUTRIDER_res.Rdata")