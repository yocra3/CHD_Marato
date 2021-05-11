#'#################################################################################
#'#################################################################################
#' Select CpGs for validation with another technique
#'#################################################################################
#'#################################################################################

# Load libraries ####
library(tidyverse)
library(cowplot)
library(minfi)
library(limma)

# Load data ####
load("results/methylation/Episignatures/cohort.training.Rdata")
load("results/methylation/finalQC_files/v5/gset.filterAnnotatedProbes.Rdata")
load("results/G2CS2CO/epimutations_chrX.Rdata")
bumps <- read.delim("results/methylation/Epimutations/v5/bumpsSel.txt", header = FALSE)
colnames(bumps) <- c("chromosome", "start", "end", "value", "area", "cluster", "indexStart",
                     "indexEnd", "L", "clusterL", "Sample", "Anova")

# Heatmap CpGs for validation ####
gset$pathClassSimp <- ifelse(gset$pathClass %in% c("Conotruncal Malformations", "Left heart hypoplasia"),
                             "Conotruncal Malformations or Left heart hypoplasia", gset$pathClass)

model <- model.matrix(~ pathClassSimp + Sex, colData(gset))
lmfit <- lmFit(getBeta(gset), design = model)
lmFite <- eBayes(lmfit)
tab <- topTable(lmFite, n = Inf)

tab.con <-  topTable(lmFite, n = Inf, coef = 2)
tab.con.f <- tab.con[rownames(tab.con) %in% svm_comb$feats, ]
tab.comp <-  topTable(lmFite, n = Inf, coef = 3)
tab.comp.f <- tab.comp[rownames(tab.comp) %in% svm_comb$feats, ]

## Select top 5 CpGs more associated with each comparison
cpgs_episig <- c(rownames(tab.con.f)[1:5], rownames(tab.comp.f)[1:5])
mat <- getBeta(gset[cpgs_episig, ])

cols <- c("#000000", "#CC6666", "#9999CC", "#66CC99")
names(cols) <- c("Control", "Conotruncal Malformations", 
                 "Left heart hypoplasia", "Other malformations")

colPlot <- cols[gset$pathClass]


png("figures/episignature_cohort_methylation_minimalCpGs.png", height = 800)
heatmap(mat, scale = "none", 
        col = gray.colors(100, start = 0, end = 1),
        ColSideColors = colPlot, cexCol = 2)
dev.off()

cpgs_episigGR <- rowRanges(gset)[cpgs_episig]
cpgs_episig.df <- data.frame(CpG = cpgs_episig, Chromosome = seqnames(cpgs_episigGR),
                             Position = start(cpgs_episigGR), Type = "Episignature", 
                             Sample = "All - 31 samples")


## Make list of Cpgs in epimutations
samps <- c("G2CS8CO", "G4CS8CO", "G4CS6CO", "G2CS4CO", "G2CS6CO", "G1CS2CO")
bumps.manuscript <- subset(bumps, Sample %in% samps)
bumps.manuscript <- subset(bumps.manuscript, !(Sample == "G4CS8CO" & chromosome == "chr11"))
bumps.manuscript <- subset(bumps.manuscript, !(Sample == "G2CS4CO" & start == 34460298))
colnames(bumps.sel) <- c("chromosome", "start", "end", "value", "area", "cluster", "indexStart",
                     "indexEnd", "L", "clusterL", "Sample", "Anova")
bumps.manuscript <- rbind(bumps.manuscript, bumps.sel)
bumpsGR <- makeGRangesFromDataFrame(bumps.manuscript, keep.extra.columns = TRUE)
cpgs_epimutGR <- subsetByOverlaps(rowRanges(gset), bumpsGR)

cpgs_epimut.df <- data.frame(CpG = names(cpgs_epimutGR), Chromosome = seqnames(cpgs_epimutGR),
                             Position = start(cpgs_epimutGR), Type = "Epimutation", 
                             Sample = rep(bumpsGR$Sample, bumpsGR$L))
cpgs_val.df <- rbind(cpgs_epimut.df, cpgs_episig.df)
write.table(cpgs_val.df, file = "scripts/Manuscript/cpgs_validation.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
