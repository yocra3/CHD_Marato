#'#################################################################################
#'#################################################################################
#' Detect candidate regions of harboring an epimutation
#' Epimutations are groups of probes having a very different methylation pattern in 
#' a samples than in the overall group. To detect epimutation, we will use the 
#' approach described in PMID: 30929737
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
gsetPath <- args[1]
sample <- args[2]

## Load libraries ####
library(minfi)

## Load data
load(gsetPath)

## Run bumphunter ####
### Define model
pdata <- colData(gset)
pdata$samp <- pdata$SampleID == sample
model <- model.matrix(~ samp, pdata)

## Run bumphunter
bumps <- bumphunter(gset, model, cutoff = 0.1)$table

if (!is.na(bumps)) {

  bumps$sample <- sample
  
  ## Select bumps with at least 5 cpgs
  bumps.sel <- subset(bumps, L >= 5)
  
  ## Compute anova
  getFstatAnova <-  function(i){
    row <- bumps.sel[i, ]
    mini <- getBeta(gset[row$indexStart:row$indexEnd, ])
    modSum <- summary(manova(t(mini) ~ model))
    modSum$stats[, 3][1]
  }
  bumps.sel$Anova <- sapply(seq_len(nrow(bumps.sel)), getFstatAnova)
  
  write.table(bumps.sel, file = "all_bumps.txt", col.names = FALSE, row.names = FALSE,
              sep = "\t", quote = FALSE)

  ## Select bumps with anova higher than 40
  bumps.final <- subset(bumps.sel, Anova > 40)
  write.table(bumps.final, file = "sel_bumps.txt", col.names = FALSE, row.names = FALSE,
              sep = "\t", quote = FALSE)
}