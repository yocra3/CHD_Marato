#'#################################################################################
#'#################################################################################
#' Generate QC reports from meffil
#' Run from base folder
#'#################################################################################
#'#################################################################################

## Load libraries
library(meffil)

## Define names
vers <- "v4"
QC_folder <- paste0("results/methylation/QC_intermediate/", vers)
reportsName <- paste("reports/meffil", vers, Sys.Date(), sep = "_")

## QC report raw
load(paste0(QC_folder, "/qcsummary.Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(reportsName, "_raw.html"))

load(paste0(QC_folder, "/qcsummary.clean.Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(reportsName, "_clean.html"))

load(paste0(QC_folder, "/norm.summary.Rdata"))
meffil.normalization.report(norm.summary, output.file = paste0(reportsName, "_normalization.html"))
