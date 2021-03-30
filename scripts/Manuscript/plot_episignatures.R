#'#################################################################################
#'#################################################################################
#' Episignaturesplot
#' Figure with general PCA and heatmaps of episignatures
#'#################################################################################
#'#################################################################################

## Working directory in project's folder

# Load libraries ####
library(tidyverse)
library(cowplot)
library(minfi)
library(DESeq2)
library(dendextend)

# Load data ####
# Methylation 
load("results/methylation/Episignatures/cohort.training.Rdata")
load("results/methylation/finalQC_files/v5/gset.autosomic.Rdata")
load("results/methylation/Episignatures/GSE62629.GenomicRatioSet.EPICcpgs.Rdata")

# RNAseq 
load("results/RNAseq/QC_files/v3/RNAseq_QCobject.Rdata")
load("results/methylation/Episignatures/transcriptomicSignature.Rdata")

# PCA methylation ####
meth.pcs <- prcomp(t(getBeta(gset)), rank. = 6)
varprop <- round(meth.pcs$sdev**2/sum(meth.pcs$sdev**2)*100, 2)

pc <- meth.pcs$x %>%
  data.frame() %>%
  mutate(group = gset$pathClass,
         group = factor(group, levels = c("Control", "Conotruncal Malformations", 
                                          "Left heart hypoplasia", "Other malformations"))) %>%
  ggplot(aes(x = PC4, y = PC6, col = group)) +
  geom_point() +
  scale_color_manual(name = "", values = c("#000000", "#CC6666", "#9999CC", "#66CC99")) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC4: ", varprop[4], "% variance")) +
  scale_y_continuous(name = paste0("PC6: ", varprop[6], "% variance")) +
  ggtitle("PCA of whole genome methylation") +
  theme(plot.title = element_text(hjust = 0.5))

png("figures/pc_all_methy.png", height = 400, width = 1120)
pc
dev.off()
pc_p <- ggdraw() + draw_image("figures/pc_all_methy.png")


# heatmap selected markers ####
mat <- getBeta(gset[svm_comb$feats, ])

cols <- c("#000000", "#CC6666", "#9999CC", "#66CC99")
names(cols) <- c("Control", "Conotruncal Malformations", 
                 "Left heart hypoplasia", "Other malformations")

colPlot <- cols[gset$pathClass]


png("figures/episignature_cohort_methylation.png", height = 800)
heatmap(mat, scale = "none", 
        col = gray.colors(100, start = 0, end = 1),
        ColSideColors = colPlot, cexCol = 2)
dev.off()

heat_markers <- ggdraw() + draw_image("figures/episignature_cohort_methylation.png")


# heatmap external dataset ####
mat.GSE <- data.matrix(getBeta(gset.GSE62629.filt[svm_comb$feats, ]))
cols2 <- c("#000000", "#999999", "blue", "green")
names(cols2) <- c("right atrium normal", "right ventricle normal", "Tetralogy of Fallot", 
                  "ventricular septal defect")

colPlot2 <- cols2[gset.GSE62629.filt$Group]


png("figures/episignature_GEO_methylation.png", height = 800)
heatmap(mat.GSE, scale = "none", 
        col = gray.colors(100, start = 0, end = 1),
        ColSideColors = colPlot2, cexCol = 2)

dev.off()

heat_GEO <- ggdraw() + draw_image("figures/episignature_GEO_methylation.png")

# heatmap RNAseq ####
selGenes <- unique(rownames(subset(eQTMs.df, padj.comb < 0.05)))
selGenes.case <- unique(rownames(subset(eQTMs.case.df, padj.comb < 0.05)))

colPlotRNA <- cols[vsd$pathClass]

png("figures/episignature_RNA.png", height = 800)
heatmap(assay(vsd[unique(selGenes), ]), scale = "none", 
        ColSideColors = colPlotRNA, cexCol = 2)
dev.off()

heat_RNA <- ggdraw() + draw_image("figures/episignature_RNA.png")

heats <- plot_grid(heat_markers, heat_GEO, heat_RNA, ncol = 3, labels = LETTERS[2:4])


png("figures/episignature_panel.png", height = 800, width = 1000)
plot_grid(pc, heats, nrow = 2, rel_heights = c(3, 5), labels = c("A", ""))
dev.off()


# heatmap external dataset full ####
hcGSE <- as.dendrogram(hclust(dist(t(getBeta(gset.GSE62629.filt)))))

ord1 <- order.dendrogram(hcGSE)
labels_colors(hcGSE) <- cols2[gset.GSE62629.filt[, ord1]$Group]

png("figures/hclust_all_GEO_methylation.png")
plot(hcGSE, main = "EPIC CpGs GSE62629")
dev.off()

# Check cpgs from  PMID: 26429603 ####
cpgsFiles <- dir("data", pattern = "26429603", full.names = TRUE)
cpgsList <- lapply(cpgsFiles, read.csv, skip = 1, as.is = TRUE)
lit_cpgs <- unique(unlist(lapply(cpgsList, function(x) x$TargetID)))

mat_lit <- getBeta(gset[rownames(gset) %in% lit_cpgs, ])

cols <- c("#000000", "#CC6666", "#9999CC", "#66CC99")
names(cols) <- c("Control", "Conotruncal Malformations", 
                 "Left heart hypoplasia", "Other malformations")

colPlot <- cols[gset$pathClass]
heatmap(mat_lit, scale = "none", 
        col = gray.colors(100, start = 0, end = 1),
        ColSideColors = colPlot, cexCol = 2)                   
