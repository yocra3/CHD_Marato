#!/usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Prepare files needed for sashimi plot
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
fraserObj <- args[1]

## Load libraries
library(FRASER)

## Load files
load(fraserObj)

## Create data.frame per sample collapsed ranges
res <- results(fdsres)

## Remove aberrations in chromosome MT
res <- res[seqnames(res) != "MT"]

## Add aberration due to a deletion in G2CS6CO
resg2c6 <- results(fdsres, padjCutoff = 0.1, sampleIDs = "G2CS6CO")
resg2c6 <- resg2c6[seqnames(resg2c6) == "14"]

## Add aberration due to a mutation in G2CS12CO
resg2c12 <-  results(fdsres, padjCutoff = 1, sampleIDs = "G2CS12CO")
resg2c12 <- subsetByOverlaps(resg2c12, GRanges("11:66052051-66056638"))
res <- c(res, resg2c6, resg2c12)

mergeRanges <- function(id){
  range <- reduce(subset(res, sampleID == id)) + 100
  df <- data.frame(range = as.character(range), 
                   bamFile = paste0(id, ".txt"),
                   rangeName = paste(id, seq_len(length(range)), sep = "."))
  df
}
map <- Reduce(rbind, lapply(unique(res$sampleID), mergeRanges ))
write.csv(map, file = "ggsashimi_vars.csv", quote = FALSE, row.names = FALSE, 
          col.names = TRUE )


## Create data.frame with sample info for each sample in FRASER
pheno <- colData(fdsres)[, c("SampleID", "bamFile", "Status")]
pheno <- subset(pheno, SampleID != "G3CS5CO")
cases <- unique(gsub(".txt", "", map$bamFile))

lapply(cases, function(id) {
  
  pheno$Status <- as.character(pheno$Status)

  pheno[pheno$SampleID == id, "Status"] <- id
  if (id == "16B12251") {
    pheno[pheno$SampleID == id, "Status"] <- "M16B12251"
  }
  write.table(pheno, file = paste0(id, ".txt"), quote = FALSE, row.names = FALSE, 
              col.names = FALSE, sep = "\t")
}
)
