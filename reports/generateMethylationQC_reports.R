#'#################################################################################
#'#################################################################################
#' Generate QC reports from meffil
#' Run from base folder
#'#################################################################################
#'#################################################################################

## Load libraries
library(meffil)
library(minfi)

## Define names
vers <- "v5"
QC_folder <- paste0("results/methylation/QC_intermediate/", vers)
reportsName <- paste("reports/meffil", vers, Sys.Date(), sep = "_")

## QC report raw
load(paste0(QC_folder, "/qcsummary.Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(reportsName, "_raw.html"))

load(paste0(QC_folder, "/qcsummary.clean.Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(reportsName, "_clean.html"))

load(paste0(QC_folder, "/norm.summary.Rdata"))
meffil.normalization.report(norm.summary, output.file = paste0(reportsName, "_normalization.html"))

load(paste0(QC_folder, "/gset.normalizedComBat.GenomicRatioSet.Rdata"))
load(paste0(QC_folder, "/norm.obj.pc.Rdata"))

beta.pcs <- meffil.methylation.pcs(getBeta(gset), probe.range = 40000)

batch_var <- c("Slide", "Array", "Sex", "GestAge", "Status", "ExtractionBatch", "pathGroup", "pathClass", "MolecularCause", "SampleBatch", "Rearrangements", "Del22q11", "N_genes")

## Define normalization parameters
norm.parameters <- meffil.normalization.parameters(
	norm.objects,
	variables = batch_var,
	control.pcs = seq_len(5),
	batch.pcs = seq_len(5),
	batch.threshold = 0.01
)
norm.summary <- meffil.normalization.summary(norm.objects, pcs = beta.pcs, parameters = norm.parameters)
meffil.normalization.report(norm.summary, output.file = paste0(reportsName, "_normalization_postComBat.html"))
