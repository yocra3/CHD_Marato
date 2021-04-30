#'#################################################################################
#'#################################################################################
#' Run enrichment of CpGs Episignature with EASIER
#' - Use docker: yocra3/rstudio_easier:1.0
#'#################################################################################
#'#################################################################################

## Load libraries
library(EASIER)
library(tidyverse)

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

# Regions enrichment ####
unsignif_df <- get_annotation_unlisted_CpGs(cpgs$rs_number, "EPIC" )
cpgs$signif <- 'yes'
unsignif_df$signif <- 'no'
cpgs <- rbind(cpgs[,c(2:47, 50)], as.data.frame(unsignif_df) )
cpgs$rs_number <- cpgs$Name

# Get descriptives
get_descriptives_GenePosition(cpgs$UCSC_RefGene_Group, cpgs$signif , "CpGlist", 
                              outputdir = paste0(outdir, "GenePosition/Fisher_CpGlist/"), 
                              outputfile =file)
GenePosition <- getAllFisherTest(cpgs$signif, cpgs$UCSC_RefGene_Group, 
                                 outputdir = paste0(outdir, "GenePosition/Fisher_CpGlist/"), 
                                 outputfile = file,
                                 plots = TRUE )
GenePosition %>%
  filter(Data != "ExonBnd") %>%
  mutate(OR = as.numeric(OR),
         OR.inf = as.numeric(OR.inf),
         OR.sup = as.numeric(OR.sup),
         Positions = factor(Data, levels = c("TSS1500", "TSS200","5'UTR", "1stExon", "Body", "3'UTR"))) %>%
  ggplot(aes(x = Positions, y = OR)) + 
  geom_bar(stat = "identity", fill = "steelblue1", width = 0.5) +
  geom_errorbar(aes(ymin = OR.inf, ymax = OR.sup), width = 0.2) + 
  scale_y_continuous(trans = "log2", 
                   breaks = scales::trans_breaks("log2", function(x) round(2^x, 2))) +
  theme_classic(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_hline(yintercept = 1) + xlab("Gene position")


plot_GenePosition(GenePosition, 
                  outputdir = paste0(outdir, "GenePosition/Fisher_CpGlist/"), 
                  outputfile = paste0("CpGlist_",file), main = )

