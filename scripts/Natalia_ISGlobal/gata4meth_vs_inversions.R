#'#################################################################################
#'#################################################################################
#' Explore correlations between CpGs in GATA4 and inv8p23.1
#'#################################################################################
#'#################################################################################

# Load libraries ####
library(scoreInvHap)
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(cowplot)

# Load methylation ####
load("results/methylation/finalQC_files/2021-03-08/gset.autosomic.Rdata")
gata4_cpgs <- read.table("data/Natalia_cpgs_gata4.txt")$V1


# Load genotypes ####
vcf_paths <- dir("results/VariantCalling/SNV/", pattern = ".gz$", full.names = TRUE)
range <- inversionGR["inv8_001"]
seqlevelsStyle(range) <- "NCBI"
## Subset VCFs to inv8 region
a <- lapply(vcf_paths, function(x) {
  vcf <- readVcf(x, genome = "hg19", 
               param = ScanVcfParam(which = range))
  id <- strsplit(basename(x), "_")[[1]][2]
  writeVcf(vcf, paste0("scripts/Natalia_ISGlobal/tempVCFs/", id, ".inv8.vcf.gz"), index = TRUE)
  })


## Merge inv8 vcfs
system("bcftools merge --force-samples scripts/Natalia_ISGlobal/tempVCFs/*.vcf.bgz -m none -o scripts/Natalia_ISGlobal/tempVCFs/merged.inv8.vcf.gz -O z")
vcf <- readVcf("scripts/Natalia_ISGlobal/tempVCFs/merged.inv8.vcf.gz", genome = "hg19")

check <- checkSNPs(vcf)
vcf.filt <- check$genos
inv8sc <- scoreInvHap(SNPlist = vcf.filt, inv = "inv8_001")
save(inv8sc, file = "scripts/Natalia_ISGlobal/inv8Genos.scoreInvHap.Rdata")


## Test difference in methylation in CpGs detected by Natalia ####
inv8class <- classification(inv8sc)
gset$inv8 <- factor(inv8class[colnames(gset)], levels = c("NN", "NI", "II"))
gset$inv8add <- as.numeric(gset$inv8) - 1
model <- model.matrix(~ inv8add + Sex + Status, colData(gset))
lmfit <- lmFit(getBeta(gset[, rownames(model)]), design = model)
lmFite <- eBayes(lmfit)
tab <- topTable(lmFite, n = Inf, coef = 2)
tab[rownames(tab) %in% gata4_cpgs, ]
mapCpgs <- gata4_cpgs[gata4_cpgs %in% rownames(gset)]

a <- lapply(mapCpgs, function(cpg){
  df <- data.frame(cpg = as.vector(getBeta(gset[cpg, ])), inv =  gset$inv8, row.names = NULL)
  df <- df[!is.na(df$inv), ]
  ggplot(df, aes(x = inv, y = methy)) + 
    geom_boxplot() +
    ggtitle(cpg)
})
plot_grid(plotlist = a)
## Test interaction between inversion and status in methylation in CpGs detected by Natalia ####
model2 <- model.matrix(~ inv8add*Status + Sex, colData(gset))
lmfit2 <- lmFit(getBeta(gset[, rownames(model2)]), design = model2)
lmFite2 <- eBayes(lmfit2)
tab2 <- topTable(lmFite2, n = Inf, coef = 5)
tab2[rownames(tab2) %in% gata4_cpgs, ]

b <- lapply(mapCpgs, function(cpg){
  df <- data.frame(methy = as.vector(getBeta(gset[cpg, ])), inv =  gset$inv8, 
                   status = gset$Status, row.names = NULL)
  df <- df[!is.na(df$inv), ]
  ggplot(df, aes(x = inv, y = methy)) + 
    geom_boxplot() +
    ggtitle(cpg) +
    facet_grid(~ status)
})
plot_grid(plotlist = b)
