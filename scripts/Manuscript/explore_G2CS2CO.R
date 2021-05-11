#'#################################################################################
#'#################################################################################
#' Explore results for G2CS2CO
#'#################################################################################
#'#################################################################################

## Wordking directory in project's folder
## Load libraries and data ####
library(DESeq2)
library(ggplot2)
library(SummarizedExperiment)
library(cowplot)
library(dplyr)
library(ggplot2)
library(tidyr)

load("results/RNAseq/Bioc_objects/v3/RNAseq_RangedSE_expressedGenes.Rdata")
load("results/RNAseq/OUTRIDER_aberrations/v4/OUTRIDER_res.Rdata")
load("results/methylation/finalQC_files/v5/gset.filterAnnotatedProbes.Rdata")

## Create DEseqData
rse$Sample <- ifelse(rse$SampleID == "G2CS2CO", "G2CS2CO", "Other")
ded <- DESeqDataSet(rse, design = ~ Status)

## Vst normalization
vst <- vst(ded)

## Select females samples in chromosome X
vst.fem_chrX <- vst[seqnames(vst) == "chrX", vst$Sex == "Female"]
png("figures/G2CS2CO_PCA_gexp_chrX.png", height = 200, width = 300)
plotPCA(vst.fem_chrX, intgroup = "Sample") +
  theme_bw() +
  ggtitle("Normalized RNAseq values") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) 
dev.off()

outres <- results(ods.res, all = TRUE)
fems <- vst.fem_chrX$SampleID

out_g2cs2 <- subset(outres, sampleID == "G2CS2CO")
out_g2cs2$chrX <- ifelse(out_g2cs2$geneID %in% rownames(vst.fem_chrX), "chrX", "Autosomic")

png("figures/G2CS2CO_aberrantgexp_chrX.png", height = 200, width = 300)
ggplot(out_g2cs2, aes(x = l2fc, color = chrX)) + geom_density() +
  xlim(c(-2, 2)) + geom_vline(xintercept = 0) +
  ggtitle("G2CS2CO chrX expression")
dev.off()


## Explore epimutations in chrX in females ####
### Subset dataset
gset.fem <- gset[seqnames(gset) == "chrX", gset$Sex == "F"]

### Define model
pdata <- colData(gset.fem)
pdata$samp <- pdata$SampleID == "G2CS2CO"
model <- model.matrix(~ samp, pdata)

## Run bumphunter
bumps <- bumphunter(gset.fem, model, cutoff = 0.1)$table
bumps$sample <- "G2CS2CO"

## Select bumps with at least 5 cpgs
bumps.sel <- subset(bumps, L >= 5)
  
## Compute anova
getFstatAnova <-  function(i){
  row <- bumps.sel[i, ]
  mini <- getBeta(gset.fem[row$indexStart:row$indexEnd, ])
  modSum <- summary(manova(t(mini) ~ model))
  modSum$stats[, 3][1]
}
bumps.sel$Anova <- sapply(seq_len(nrow(bumps.sel)), getFstatAnova)
save(bumps.sel, file = "results/G2CS2CO/epimutations_chrX.Rdata")

## Select bumps with anova higher than 40
bumps.final <- subset(bumps.sel, Anova > 40)
## No epimutations found

makePlot <- function(row, bumps, gset){
  
  sample <- bumps[row, "sample"]
  idx <- bumps[row, "indexStart"]:bumps[row, "indexEnd"]
  
  beta <- data.frame(getBeta(gset[idx, ]))
  beta$pos <- start(rowRanges(gset[idx, ]))
  beta %>%
    gather("Sample", "Methylation", seq_len(ncol(gset))) %>%
    mutate(Sample = ifelse(Sample == "X16A66", "16A66", Sample)) %>%
    ggplot(aes(x = pos, y = Methylation, color = Sample == sample, group = Sample)) +
    geom_point() + 
    geom_line() +
    ggtitle(paste(sample, bumps[row, "chromosome"], bumps[row, "start"], bumps[row, "end"])) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
}
png("results/G2CS2CO/epimutation.png")
makePlot(1, bumps.sel, gset.fem)
dev.off()

png("results/G2CS2CO/methylation_heatmap_females_chrX.png", height = 800)
heatmap(getBeta(gset.fem), scale = "none", 
        col = gray.colors(100, start = 0, end = 1))
dev.off()