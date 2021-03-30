#'#################################################################################
#'#################################################################################
#' Plot skewed transcription of chromosome X
#'#################################################################################
#'#################################################################################

# Load libraries and data ####
library(ggplot2)
library(SummarizedExperiment)
library(dplyr)
library(cowplot)

load("results/methylation/finalQC_files/v5/gset.filterAnnotatedProbes.Rdata")
load("results/RNAseq/ASE/v1/tMAE_res.Rdata")

## Rename samples to match methylation names
resMAE.filt <- resMAE[ refAllele != "N"]
females <- colnames(gset)[gset$Sex == "F"]
resMAE.filt$MAE_ID <- gsub("-", "", resMAE.filt$MAE_ID)
females.sel <- females[females %in% unique(resMAE.filt$MAE_ID)]
  
## Plot Chromosome X MAE in females ####
mae.plot <- resMAE.filt %>%
    filter(MAE_ID %in% females.sel & contig == "X") %>%
    ggplot(aes(x = position/1e6, y = altRatio)) +
    geom_point() +
    ggtitle(paste("Allelic expression of chromosome X")) +
    facet_wrap(~ MAE_ID, nrow = 3) +
    theme_bw() +
    scale_x_continuous(name = "Coordinates (in Mb)") +
    scale_y_continuous(name = "Alternative allele expression") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = c(0.9, 0.1), linetype = "dashed")


## Plot heatmap chromosome X in females ####
beta_X_fem <- getBeta(gset[seqnames(rowRanges(gset)) == "chrX", gset$Sex == "F"])

png("figures/chrX_methylation.png", height = 700)
heatmap(beta_X_fem, Rowv = NA, scale = "none", 
        col = gray.colors(100, start = 0, end = 1))
dev.off()

heat_chrX <- ggdraw() + draw_image("figures/chrX_methylation.png")
plot_grid(mae.plot, heat_chrX, ncol = 2)

## Plot Distribution Methylation in chrX females ####
pc <- prcomp(t(beta_X_fem))
vars <- round((pc$sdev)^2/sum((pc$sdev)^2)*100, 2)
meth_pc <- pc$x %>%
  data.frame() %>%
  mutate(Sample = colnames(beta_X_fem),
         Skewing = ifelse(!Sample %in% c("G2CS2CO", "G2CS13CO"), "Cohort", Sample)) %>%
  ggplot(aes(x = PC1, y = PC2, color = Skewing)) +
  ggtitle("PCA of chrX methylation in females") + 
  geom_point() +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1: ", vars[1], "% variance")) +
  scale_y_continuous(name = paste0("PC2: ", vars[2], "% variance")) +
  scale_color_manual(name = "", values = c("black", "red", "blue")) +
  theme(plot.title = element_text(hjust = 0.5))

bot <- plot_grid(meth_pc, heat_chrX, nrow = 2)

png("figures/chrX_skewing.png", width = 700)
plot_grid(mae.plot, meth_pc, ncol = 2, labels = c("A", "B"))
dev.off()
