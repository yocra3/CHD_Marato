#'#################################################################################
#'#################################################################################
#' Figure RNAseq clustering
#'#################################################################################
#'#################################################################################

## Wordking directory in project's folder
## Load libraries and data ####
library(DESeq2)
library(ggplot2)

load("results/RNAseq/QC_files/v3/RNAseq_QCobject.Rdata")

vp <- plotPCA(vsd, intgroup = "Status") +
  theme_bw() +
  ggtitle("Normalized RNAseq values") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)) 
png("figures/RNAseq_clustering.png", height = 350, width = 700)
vp
dev.off()
