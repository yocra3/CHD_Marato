#'#################################################################################
#'#################################################################################
#' Generate files for CIBERSORTX
#'#################################################################################
#'#################################################################################

# Load libraries
library(recount)
library(SummarizedExperiment)

load("results/RNAseq/Bioc_objects/v3/RNAseq_RangedSE_expressedGenes.Rdata")

expmat <- getTPM(rse, length_var = NULL)
rownames(expmat) <- gsub("\\..*", "", rownames(expmat))
expdf <- data.frame(Gene = rownames(expmat), expmat)
write.table(expdf, file = "results/RNAseq/CIBERSORT/RNAseq_mat.tab", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
