#'#################################################################################
#'#################################################################################
#' Generate QC reports from meffil
#' Run from base folder
#'#################################################################################
#'#################################################################################

## Load libraries
library(meffil)

## Define names
vers <- "v1"
QC_folder <- paste0("results/methylation/QC_intermediate/", vers)
reportsName <- paste("reports/meffil", vers, Sys.Date(), sep = "_")

## QC report raw
load(paste0(QC_folder, "/qcsummary.Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(reportsName, "_raw.html"))

load(paste0(QC_folder, "/qcsummary.clean.Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(reportsName, "_clean.html"))

load(paste0(QC_folder, "/norm.summary.Rdata"))
meffil.normalization.report(norm.summary, output.file = paste0(reportsName, "_normalization.html"))



# ## Check correlation in duplicated sample
# cpgMeans <- rowMeans(norm.beta)
# correlsglobal <- cor(norm.beta - cpgMeans)
# 
# ## Add to report
# hist(correlsglobal[upper.tri(correlsglobal)], main = "Mean-centered correlation", xlab = "Pearson r")
# median(correlsglobal[upper.tri(correlsglobal)])
# correlsglobal["G2CS12COBIS", "G2CS12CO"]
# ## Median centered correlation is close to 0, as expected
# ## Highest correlation is found for duplicated samples (G2CS12CO)
