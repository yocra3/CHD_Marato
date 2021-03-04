#!/usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Create FRASER object
#' This script quantifies splice junctions and creates FRASER object for downstream
#' analysis
#' It can be reused for other projects but modifying some parts, highlighted with ##
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
phenoPath <- args[1]
bamPath <- args[2]
cores <- as.numeric(args[3])

## Load libraries
library(FRASER)
library(dplyr)

## Create FRASER settings
load(phenoPath)
bams <- dir(bamPath, pattern = "bam$", full.names = TRUE)
bamsdf <- data.frame(bamFile = bams, stringsAsFactors = FALSE) %>%
  mutate(sampleT = strsplit(bamFile, "/"),
         sampleT = sapply(sampleT, function(x) x[[length(x)]]),
         sampleT = gsub(".bam", "", sampleT),
         sampleID = gsub("-", "", sampleT)) %>%
  select(bamFile, sampleID)

colData <- pheno %>%
  mutate(sampleID = SampleID) %>%
  right_join(bamsdf, by = "sampleID")

settings <- FraseRDataSet(colData = DataFrame(colData), 
                          workingDir = ".")


## Generate FRASER object
fds <- countRNAData(settings, NcpuPerSample = cores, countDir = ".")

## Compute PSI values
fds <- calculatePSIValues(fds)
save(fds, file = "FRASER_obj.Rdata")