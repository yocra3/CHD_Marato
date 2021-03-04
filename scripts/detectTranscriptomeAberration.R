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
ods.filt <- filterExpression(ods, gtfFile = gtfFile, addExpressedGenes = TRUE)

register(SerialParam())

message("Size Factors")
ods.cor <- estimateSizeFactors(ods.filt)
message("Encoding dimentions")
ods.cor <- findEncodingDim(ods.cor)
message("Confounders")
ods.cor <- controlForConfounders(ods.cor)

# run full OUTRIDER pipeline
ods.res <- OUTRIDER(ods.cor)
save(ods.filt, ods.cor, file = "OUTRIDER_QC.Rdata")
save(ods.res, file = "OUTRIDER_res.Rdata")