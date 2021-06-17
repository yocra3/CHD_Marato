#'#################################################################################
#'#################################################################################
#' Explore correlations between CpGs in GATA4 and CHD
#'#################################################################################
#'#################################################################################

# Load libraries and data ####
library(minfi)
library(limma)
library(dplyr)

load("results/methylation/finalQC_files/2021-03-08/gset.autosomic.Rdata")
gata4_cpgs <- read.table("data/Natalia_cpgs_gata4.txt")$V1

## Test difference in methylation in CpGs detected by Natalia ####
pheno <- colData(gset)
pheno$Status <- factor(pheno$Status , levels = c("Control", "Case"))
model <- model.matrix(~ Status + Sex, pheno)
lmfit <- lmFit(getBeta(gset), design = model)
lmFite <- eBayes(lmfit)
tab <- topTable(lmFite, n = Inf, coef = 2)
tab[rownames(tab) %in% gata4_cpgs, ]

gset.gata4 <- rowRanges(gset)[rownames(gset) %in% gata4_cpgs]

mat <- getBeta(gset.gata4)
cols <- c("#000000", "#CC6666", "#9999CC", "#66CC99")
names(cols) <- c("Control", "Conotruncal Malformations", 
                 "Left heart hypoplasia", "Other malformations")

colPlot <- cols[gset$pathClass]
heatmap(mat, scale = "none", 
        col = gray.colors(100, start = 0, end = 1),
        ColSideColors = colPlot, cexCol = 2)
par(mfrow = c(4, 5))
lapply(rownames(mat), function(cpg){
  df <- data.frame(cpg = as.vector(getBeta(gset[cpg, ])), status =  gset$Status)
  boxplot(cpg ~ status, main = cpg, df)
})
hist(tab[rownames(mat) , ] $logFC)


## Test difference in methylation in CpGs in GATA4 region ####
gset.reg <- subsetByOverlaps(gset, range(rowRanges(gset.gata4)))
tab[rownames(gset.reg) , ] %>% arrange(P.Value) %>% head()
plot(density(tab[rownames(gset.reg) , ] $logFC))


heatmap(getBeta(gset.reg), scale = "none", 
        col = gray.colors(100, start = 0, end = 1),
        ColSideColors = colPlot, cexCol = 2)

par(mfrow = c(4, 5))
lapply(rownames(gset.reg), function(cpg){
  df <- data.frame(cpg = as.vector(getBeta(gset[cpg, ])), status =  gset$Status)
  boxplot(cpg ~ status, main = cpg, df)
})

## Test CpGs in Clara's paper ####
lit_cpgs <- c("cg09626984", "cg20279283", "cg01546563", "cg19172575", "cg13434842")
tab[lit_cpgs , ] %>% arrange(P.Value) %>% head()
heatmap(getBeta(gset[lit_cpgs, ]), scale = "none", 
        col = gray.colors(100, start = 0, end = 1),
        ColSideColors = colPlot, cexCol = 2)
plot(density(tab[lit_cpgs , ] $logFC))
hist(tab[lit_cpgs , ] $logFC)
par(mfrow = c(2, 3))
lapply(lit_cpgs, function(cpg){
  df <- data.frame(cpg = as.vector(getBeta(gset[cpg, ])), status =  gset$Status)
  boxplot(cpg ~ status, main = cpg, df)
})
