#'#################################################################################
#'#################################################################################
#' Apply ComBat to a GenomicRatioSet
#' ComBat will be applied to M-values
#' The result will be stored as a GenomicRatioSet
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(minfi)
library(sva)

load("gset.Rdata")

pheno <- colData(gset)
##### Add pheno for technical replicate
pheno["G2CS12COBIS", ] <- pheno["G2CS12CO", ]
batch <- pheno[["SampleBatch"]]

## Create ComBat model
modcombat <- model.matrix(~ Sex + pathClass, data = pheno)

## Run comBat
m <- getM(gset)
m[m == -Inf] <- -10
m[m == Inf] <- 10

combat_M <- ComBat(dat = m, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots = FALSE)

beta <- ilogit2(combat_M)
assay(gset) <- beta
save(gset, file = "gset.normalizedComBat.GenomicRatioSet.Rdata")
