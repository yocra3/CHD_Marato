#!/usr/local/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Prioritize CNVs are reformat data from CNV
#' Remove CNVs that fulfill these criteria:
#' - Common CNVs: reference CNV contains >20% sample CNV. Common calls are from nstd186 (AF > 0.01)
#' - ClinVar Bening CNVs: overlap with ClinVar benign or likely benign calls from nstd102. Reference CNV contains >20% sample CNV.
#' - Segmental Duplications: overlap with segmental duplications from UCSC. Segmental duplication contains >50% sample CNV.
#' Create three tables with selected CNVs:
#' - Pathogenic: overlap with ClinVar pathogenic or likely pathogenic calls from nstd102. Sample CNV contains >80% reference CNV.
#' - OMIM: Overlap with exons of OMIM genes
#' - GENCODE: Overlap with exons of any GENCODE gene
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
cnvFile <- args[1]
outFile <- args[2]

# Load libraries
library(VariantAnnotation)
library(openxlsx)

## Load data
cnv <- readVcf(cnvFile)

vars <- info(cnv)
ini <- nrow(vars)

# Remove CNVs ####
## Remove common CNVS
com <- sum(vars$commonCNV20 != "")
vars.com <- vars[vars$commonCNV20 == "", ]

## Remove benign CNVS
ben <- sum(vars.com$beningCNV20 != "")
vars.ben <- vars.com[vars.com$beningCNV20 == "", ]

## Remove CNVS in segmental Duplications
seg <- sum(vars.ben$segDups50 != "")
vars.filt <- vars.ben[vars.ben$segDups50 == "", ]

# Create sets ####
## Create pathogenic set
vars.path <- vars.filt[vars.filt$pathoCNV80 != "", ]
path <- nrow(vars.path)

## Create OMIM set
vars.omim <- vars.filt[vars.filt$omimGENES != "" & vars.filt$gencodeEXONS != "", ]
omim <- nrow(vars.omim)

## Create GENCODE set (Exclude CNVs already in OMIM)
vars.gencode <- vars.filt[vars.filt$omimGENES == "" & vars.filt$gencodeGENES != "" & vars.filt$gencodeEXONS != "", ]
gencode <- nrow(vars.gencode)

# Create Excel sheets ####
createGRfromVCF <- function(vcf){
  
  cnvGR <- rowRanges(vcf)
  mcols(cnvGR) <- info(vcf)
  cnvGR$POS <- start(cnvGR)
  start(cnvGR) <- pmin(cnvGR$POS, cnvGR$END)
  end(cnvGR) <- pmax(cnvGR$POS, cnvGR$END)
  cnvGR
}

createReportTable <- function(filtTab, cnv){
  
  ## Select columns
  selCols <- c("SVLEN", "SVTYPE", "commonCNV", "clinCNV", "gencodeGENES", "omimGENES")
  df <- filtTab[, selCols]
  cnv.filt <- cnv[rownames(df), ]
  
  df$ReadDepth <- vapply(geno(cnv.filt)$NRD, function(x) 
    mean(as.numeric(strsplit(x, "|", fixed = TRUE)[[1]])), numeric(1))
  df$E <- geno(cnv.filt)$E[, 1]
  df$q0 <- geno(cnv.filt)$q0[, 1]
  df$Coordinates <- as.character(createGRfromVCF(cnv.filt))
  df$SVLEN <- width(GRanges(df$Coordinates))
  df[, c("Coordinates", "ReadDepth", "E", "q0", selCols)]
}

pathdf <- createReportTable(vars.path, cnv)
omimdf <- createReportTable(vars.omim, cnv)
gencodedf <- createReportTable(vars.gencode, cnv)

sumTable <- data.frame(Description = c(
  "Initial Number CNVs",
  "Common CNVs",
  "ClinVar Benign CNVs",
  "CNVs in Segmental Duplications",
  "Total ClinVar Pathogenic CNVs",
  "Total CNVs in OMIM genes",
  "Total CNVs in other GENCODE genes"),
  Number = c(ini, com, ben, seg, path, omim, gencode))

write.xlsx(list(pathdf, omimdf, gencodedf, sumTable), 
           file = outFile,
           rowNames = FALSE,
           colNames = TRUE, 
           sheetName = c("ClinVar Pathogenic CNVs",
                         "CNVs in OMIM genes", 
                         "CNVs in GENCODE genes",
                         "Prioritization summary"))