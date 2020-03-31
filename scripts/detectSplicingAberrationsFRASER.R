#!/usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Run FRASER analysis from FRASER object with psi values computed.
#' It can be reused for other projects but modifying some parts, highlighted with ##
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
fraser <- args[1]

## Load libraries
library(FRASER)

## Load object
load(fraser)

## Filter lowly expressed junctions
fdsQC1 <- filterExpression(fds, filter = FALSE)
fdsfilt <- fdsQC1[mcols(fdsQC1, type="j")[,"passed"],]

## Compute Hyperparameter
register(SerialParam())
fdsfilt <- optimHyperParams(fdsfilt, correction = "PCA", type = "psi5", 
                            plot = FALSE, BPPARAM = SerialParam(),
                            internalThreads = 1)
fdsfilt <- optimHyperParams(fdsfilt, correction = "PCA", type = "psi3", 
                            plot = FALSE, BPPARAM = SerialParam(),
                            internalThreads = 1)
fdsfilt <- optimHyperParams(fdsfilt, correction = "PCA", type = "psiSite", 
                            plot = FALSE, BPPARAM = SerialParam(),
                            internalThreads = 1)
## Run FRASER
types <- c("psi5", "psi3", "psiSite")
names(types) <- types
qs <- sapply(types, bestQ, fds = fdsfilt)
fdsres <- FraseR(fdsfilt, q = qs, correction = "PCA")
save(fdsQC1, file = "FRASER_QC.Rdata")
save(fdsres, file = "FRASER_results.Rdata")