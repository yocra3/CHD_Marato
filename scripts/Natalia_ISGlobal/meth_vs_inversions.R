#'#################################################################################
#'#################################################################################
#' Explore correlations between DNA methylation and inversions
#'#################################################################################
#'#################################################################################

# Load libraries ####
library(scoreInvHap)
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(cowplot)
library(limma)
# Load methylation ####
load("results/methylation/finalQC_files/2021-03-08/gset.autosomic.Rdata")


# Define inversions ####
vcf_paths <- dir("results/VariantCalling/SNV/", pattern = ".gz$", full.names = TRUE)

## inv17q21.31 ####
range <- inversionGR["inv17_007"]
seqlevelsStyle(range) <- "NCBI"
## Subset VCFs to inv17 region
a <- lapply(vcf_paths, function(x) {
  vcf <- readVcf(x, genome = "hg19", 
                 param = ScanVcfParam(which = range))
  id <- strsplit(basename(x), "_")[[1]][2]
  writeVcf(vcf, paste0("scripts/Natalia_ISGlobal/tempVCFs/", id, ".inv17.vcf.gz"), index = TRUE)
})
## Merge inv17 vcfs
system("bcftools merge --force-samples scripts/Natalia_ISGlobal/tempVCFs/*inv17.vcf.bgz -m none -o scripts/Natalia_ISGlobal/tempVCFs/merged.inv17.vcf.gz -O z")
vcf <- readVcf("scripts/Natalia_ISGlobal/tempVCFs/merged.inv17.vcf.gz", genome = "hg19")

check <- checkSNPs(vcf)
vcf.filt <- check$genos
inv17sc <- scoreInvHap(SNPlist = vcf.filt, inv = "inv17_007")
save(inv17sc, file = "scripts/Natalia_ISGlobal/inv17Genos.scoreInvHap.Rdata")

## inv16 ####
range <- inversionGR["inv16_009"]
seqlevelsStyle(range) <- "NCBI"
## Subset VCFs to inv17 region
a <- lapply(vcf_paths, function(x) {
  vcf <- readVcf(x, genome = "hg19", 
                 param = ScanVcfParam(which = range))
  id <- strsplit(basename(x), "_")[[1]][2]
  writeVcf(vcf, paste0("scripts/Natalia_ISGlobal/tempVCFs/", id, ".inv16.vcf.gz"), index = TRUE)
})
## Merge inv17 vcfs
system("bcftools merge --force-samples scripts/Natalia_ISGlobal/tempVCFs/*inv16.vcf.bgz -m none -o scripts/Natalia_ISGlobal/tempVCFs/merged.inv16.vcf.gz -O z")
vcf <- readVcf("scripts/Natalia_ISGlobal/tempVCFs/merged.inv16.vcf.gz", genome = "hg19")

check <- checkSNPs(vcf)
vcf.filt <- check$genos
inv16sc <- scoreInvHap(SNPlist = vcf.filt, inv = "inv16_009")
save(inv16sc, file = "scripts/Natalia_ISGlobal/inv16Genos.scoreInvHap.Rdata")


# Test difference in methylation due to inversions in cases ####
## inv8 ####
cases <- gset[, gset$Status != "Control"]

load("scripts/Natalia_ISGlobal/inv8Genos.scoreInvHap.Rdata")

inv8class <- classification(inv8sc)
cases$inv8 <- factor(inv8class[colnames(cases)], levels = c("NN", "NI", "II"))
cases$inv8add <- as.numeric(cases$inv8) - 1
modelcases <- model.matrix(~ inv8add + Sex, colData(cases))
lmfitcas <- lmFit(getBeta(cases[, rownames(modelcases)]), design = modelcases)
lmFitecas <- eBayes(lmfitcas)
tabcas <- topTable(lmFitecas, n = Inf, coef = 2)
tabcas <- cbind(tabcas, data.frame(rowRanges(cases)[rownames(tabcas)]))
tabcas8 <- subset(tabcas, seqnames == "chr8")
save(tabcas8, file = "scripts/Natalia_ISGlobal/inv8_cpgAssoc.Rdata")

## inv17 ####
load("scripts/Natalia_ISGlobal/inv17Genos.scoreInvHap.Rdata")

inv17class <- classification(inv17sc)
cases$inv17 <- factor(inv17class[colnames(cases)], levels = c("NN", "NI", "II"))
cases$inv17add <- as.numeric(cases$inv17) - 1
modelcases <- model.matrix(~ inv17add + Sex, colData(cases))
lmfitcas <- lmFit(getBeta(cases[, rownames(modelcases)]), design = modelcases)
lmFitecas <- eBayes(lmfitcas)
tabcas <- topTable(lmFitecas, n = Inf, coef = 2)
tabcas <- cbind(tabcas, data.frame(rowRanges(cases)[rownames(tabcas)]))
tabcas17 <- subset(tabcas, seqnames == "chr17")
save(tabcas17, file = "scripts/Natalia_ISGlobal/inv17_cpgAssoc.Rdata")

## inv16 ####
load("scripts/Natalia_ISGlobal/inv16Genos.scoreInvHap.Rdata")

inv16class <- classification(inv16sc)
cases$inv16 <- factor(inv16class[colnames(cases)], levels = c("NN", "NI", "II"))
cases$inv16add <- as.numeric(cases$inv16) - 1
modelcases <- model.matrix(~ inv16add + Sex, colData(cases))
lmfitcas <- lmFit(getBeta(cases[, rownames(modelcases)]), design = modelcases)
lmFitecas <- eBayes(lmfitcas)
tabcas <- topTable(lmFitecas, n = Inf, coef = 2)
tabcas <- cbind(tabcas, data.frame(rowRanges(cases)[rownames(tabcas)]))
tabcas16 <- subset(tabcas, seqnames == "chr16")
save(tabcas16, file = "scripts/Natalia_ISGlobal/inv16_cpgAssoc.Rdata")
