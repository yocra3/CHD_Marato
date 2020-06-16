#!/usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Run test for MonoAllelic Expression (MAE) with tMAE
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
ASEtab <- args[1]

## Load libraries
library(tMAE)

## Run tMAE
mae <- fread(ASEtab, header = FALSE)
colnames(mae) <- c("contig",	"position",	"variantID",	"refAllele",	"altAllele",
	"refCount",	"altCount",	"totalCount",	"lowMAPQDepth",	"lowBaseQDepth",	"rawDepth",
  "otherBases",	"improperPairs",	"MAE_ID")
resMAE <- DESeq4MAE(mae, minCoverage = 10)

save(resMAE, file = "tMAE_res.Rdata")
