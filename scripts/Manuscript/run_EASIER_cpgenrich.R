#'#################################################################################
#'#################################################################################
#' Run enrichment of CpGs Episignature with EASIER
#' - Use docker: yocra3/rstudio_easier:1.0
#'#################################################################################
#'#################################################################################

## Load libraries
library(EASIER)

## Read CpGs ####
load("results/methylation/Episignatures/cohort.training.Rdata")

file <- "CpGs.Episignatures.txt"
outdir <- "results/methylation/Episignatures/"

cpgs <- get_annotattions( svm_comb$feats, artype = "EPIC", filename = file,
                          outdir = outdir)
cpgs$chromosome <- substr(cpgs$chr, 4, length(cpgs$chr))
cpgs$rs_number <- cpgs$CpGs

## Enrichment with missmethyl ####
## Passar argumentos?
miss_enrich <- missMethyl_enrichment(cpgs, outdir, file,  "EPIC", all = TRUE, plots = FALSE )

# Molecular Signatures Database enrichment ####
msd_enrich <- MSigDB_enrichment(cpgs, outdir, file,  "EPIC", all = TRUE)

# Enrichment with ConsensusPathDB ####
geneUniv <- list(episignature = getUniqueGenes(cpgs$UCSC_RefGene_Name))


acFSet <- c('C', 'P', 'G2', 'G3')
acType <- 'entrez-gene'

CPDB_enrich <- lapply(names(geneUniv), function( data, accFSet, genes ) {
  print(data)
  lapply(accFSet,
         get_consensusPdb_OverRepresentation,
         entityType='genes',
         accNumbers=na.omit(as.character(eval(parse(text = paste0("genes$",data))))),
         accType=acType,
         outputdir = NULL)},
  accFSet = acFSet, genes = geneUniv)