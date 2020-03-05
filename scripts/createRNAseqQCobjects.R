#'#################################################################################
#'#################################################################################
#' Create QC objects for RNAseq analysis
#' This code generates a SummarizedExperiment from a matrix of RNAseq counts.
#' It can be reused for other projects but modifying some parts, highlighted with ##
#' Things to take into account:
#' - Genes annotation is loaded from gff file. 
#' - Counts data is expected to come from gene quantification.
#' - Sample names is counts is corrected to match pheno names.
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
rseFile <- args[1]

## Load libraries ####
library(DESeq2)

## Load dataset
load(rseFile)

## Create DEseqData
ded <- DESeqDataSet(rse, design = ~ Status)

## Vst normalization
vsd <- vst(ded)
distVSD <- dist(t(assay(vsd)))

## rld normalization
rld <- rlog(ded)
distRLD <- dist(t(assay(rld)))

save(vsd, distVSD, rld, distRLD, file = "RNAseq_QCobject.Rdata")
