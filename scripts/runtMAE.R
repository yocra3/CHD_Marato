#!/usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Run test for MonoAllelic Expression (MAE) with tMAE
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
ASEtab <- args[1]

## Load libraries
library(tMAE)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)

## Run tMAE 
### Restrict to position with at least 20 SNPs
mae <- fread(ASEtab, header = FALSE)
colnames(mae) <- c("contig",	"position",	"variantID",	"refAllele",	"altAllele",
	"refCount",	"altCount",	"totalCount",	"lowMAPQDepth",	"lowBaseQDepth",	"rawDepth",
  "otherBases",	"improperPairs",	"MAE_ID")
resMAE <- DESeq4MAE(mae, minCoverage = 20)

## Annotate variants to genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
rng2 <- GRanges(seqnames = paste0("chr", resMAE$contig),
                ranges = IRanges(start = resMAE$position, end = resMAE$position),
                strand = '*')
loc <- locateVariants(rng2, txdb, AllVariants())
loc.gen <- loc[!is.na(loc$GENEID)]
mcols(loc.gen) <- mcols(loc.gen)[, c("LOCATION", "GENEID")]
loc.gen.u <- unique(loc.gen)
loc.gen.df <- data.frame(loc.gen.u)
loc.gen.df$posID <- paste0(loc.gen.df$seqnames, loc.gen.df$start)

## Discard positions in introns
resMAE <- resMAE %>%
  mutate(posID = paste0("chr", contig, position),
         MAE = ifelse(resMAE$altRatio > 0.1 & resMAE$altRatio < 0.9, "biallelic", "monoallelic")) %>%
  full_join(dplyr::select(loc.gen.df, c(posID, LOCATION, GENEID)), by = "posID")


genes <- unique(resMAE$GENEID)
genesMAE <- sapply(genes, function(g){
  resmae_gene <- subset(resMAE, GENEID == g & LOCATION == "coding")
  if (nrow(resmae_gene) == 0){
    return(NA)
  } else{
    return(mean(resmae_gene$MAE == "biallelic"))
  }
})
genesDF <- data.frame(GENEID = names(genesMAE), 
                      gene_MAE = ifelse(is.na(genesMAE), "non-coding", 
                                              ifelse(genesMAE > 0.5, "biallelic", "monoallelic")),
                      row.names = NULL)
                      
resMAE <- resMAE %>%
  left_join(genesDF, by = "GENEID")


save(resMAE, file = "tMAE_res.Rdata")
