#'#################################################################################
#'#################################################################################
#' Figures of OUTLIER expression
#' Important!! This file also contains the test for differential expression of CNVs
#'#################################################################################
#'#################################################################################

## Wordking directory in project's folder

# Load libraries ####
library(dplyr)
library(OUTRIDER)
library(cowplot)
library(GenomicRanges)
library(Gviz)
library(biomaRt)
library(ggplot2)

## Define functions ####
varOutliers <- function(index, varRange, outrider, geneRanges, window){
  
  selGR <- varRange[index]
  selGenes <- subsetByOverlaps(geneRanges, selGR + window)
  subset(outrider, geneID %in% names(selGenes) & 
           sampleID == selGR$SampleID)
  
}
gvizPlot <- function(region, cnv){
  
  
  ideo <- IdeogramTrack(genome = "hg19", chromosome = as.character(seqnames(region)))
  basetracks <- GenomeAxisTrack()
  
  bm <- useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")
  biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = as.character(seqnames(region)), 
                                      start = start(region), 
                                      end = end(region), name = "ENSEMBL", biomart = bm)
  tracks_Highlight <- HighlightTrack(trackList = list(basetracks, biomTrack),
                                           start = start(cnv), end = end(cnv),
                                           chromosome = as.character(seqnames(region)),
                                           col = "#7EA577", fill = "#C6D7C3",
                                           alpha = 0.4,
                                           inBackground = FALSE)
  plotTracks(c(ideo, tracks_Highlight),
             transcriptAnnotation = "SYMBOL")
  
}

getOutlierDist <- function(GR){
  
  outs <- varOutliers(index = 1, varRange = GR, 
                      outrider = outres, geneRanges = GRgenes, window = 1)
  all <- subset(outres, sampleID == GR$SampleID)
  all$class <- ifelse(all$geneID %in% outs$geneID, "Genes in SV", "Genes outside SV")
  all
}

plotDistribution <- function(df){
  ggplot(df, aes(x = l2fc, fill = class)) + 
    geom_density(alpha = 0.5) +
    theme_bw() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(name = "log2FC", limits = c(-2.5, 2.5)) +
    scale_y_continuous(name = "Density")
}


# Load data ####
load("results/RNAseq/OUTRIDER_aberrations/v3/OUTRIDER_res.Rdata")
load("results/RNAseq/Bioc_objects/v3/RNAseq_RangedSE_autosomicGenes.Rdata")
GRgenes <- rowRanges(rse)

## Load variants 
vars <- read.delim("results/SNVs/prioritizedSNVsTable.tab", header = TRUE)
varsGR <- GRanges(vars$"Coordinates..GRCh37.")
varsGR$SampleID <- vars$ID
varsGR$Variant.type <- vars$Variant.type

## Select samples with gene expression and CNVs
varsGR.gexp <- varsGR[varsGR$SampleID %in% ods.res$SampleID & varsGR$Variant.type %in% c("Duplication", "Deletion")]

## Divide in large and short CNVs
cnvsGR <- varsGR.gexp[width(varsGR.gexp) < 2e5]
svGR <- varsGR.gexp[width(varsGR.gexp) > 2e5]

outres <- results(ods.res, all = TRUE)

## CNVs ####
cnvResL <- lapply(seq_len(length(cnvsGR)), varOutliers, varRange = cnvsGR,
                  outrider = outres, geneRanges = GRgenes, window = 100e3)
cnvRes <- Reduce(rbind, cnvResL)
cnvRes$Variant.type <- rep(cnvsGR$Variant.type, sapply(cnvResL, nrow))
cnvRessig <- subset(cnvRes, pValue < 0.05) %>% unique() 

## G2CS2CO ####
png("figures/G2CS2CO_gviz.png", height = 400, width = 960)
gvizPlot(cnvsGR[1] + 15e3, cnvsGR[1])
dev.off()

lama4 <- plotExpectedVsObservedCounts(ods.res, "ENSG00000112769.20_7", basePlot = TRUE,
                                     groups = factor(ifelse(ods.res$SampleID == "G2CS2CO", 
                                                            "G2CS2CO", "Population")), 
                                     main = "LAMA4 Expression", pch = 19) +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))
g2cs2igv <- ggdraw() + draw_image("figures/IGV_figures/G2CS2CO_chr6del.png")
g2cs2gviz <- ggdraw() + draw_image("figures/G2CS2CO_gviz.png")

g2cs2genes <- plot_grid(g2cs2igv, g2cs2gviz, nrow = 2, labels = c("A", "B"))
g2cs2 <- plot_grid(g2cs2genes, lama4, ncol = 2, labels = c("", "C"))

png("figures/G2CS2CO_panel.png", height = 800, width = 1800)
g2cs2
dev.off()

## G2CS6CO ####
png("figures/G2CS6CO_gviz.png", height = 400, width = 960)
gvizPlot(cnvsGR[2], cnvsGR[2])
dev.off()

mnat1 <- plotExpectedVsObservedCounts(ods.res, "ENSG00000020426.11_4", basePlot = TRUE,
                                      groups = factor(ifelse(ods.res$SampleID == "G2CS6CO", 
                                                             "G2CS6CO", "Population")), 
                                      main = "MNAT1 Expression") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))
g2cs6igv <- ggdraw() + draw_image("figures/IGV_figures/G2CS6CO_chr14del.png")
g2cs6gviz <- ggdraw() + draw_image("figures/G2CS6CO_gviz.png")
g2cs6sp <- ggdraw() + draw_image("results/RNAseq/FRASER_splicing_aberrations_plots/v2/G2CS6CO.1.png")

g2cs6genes <- plot_grid(g2cs6igv, g2cs6gviz, nrow = 2, labels = c("A", "B"))
g2cs6 <- plot_grid(g2cs6genes, mnat1, ncol = 2, labels = c("", "C"))
g2cs6 <-  plot_grid(g2cs6gviz, g2cs6igv, mnat1, g2cs6sp, labels = LETTERS[1:4],
                    rel_heights = c(2, 3))
png("figures/G2CS6CO_panel.png", height = 1100, width = 1800)
g2cs6
dev.off()


## Non-significant gene
png("figures/G2CS6COb_gviz.png", height = 400, width = 960)
gvizPlot(cnvsGR[3], cnvsGR[3])
dev.off()

obscn <- plotExpectedVsObservedCounts(ods.res, "ENSG00000154358.21_5", basePlot = TRUE,
                                      groups = factor(ifelse(ods.res$SampleID == "G2CS6CO", 
                                                             "G2CS6CO", "Population")), 
                                      main = "OBSCN Expression") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))
g2cs6igvb <- ggdraw() + draw_image("figures/IGV_figures/G2CS6CO_chr1dup.png")
g2cs6gvizb <- ggdraw() + draw_image("figures/G2CS6COb_gviz.png")

g2cs6bgenes <- plot_grid(g2cs6igvb, g2cs6gvizb, nrow = 2, labels = c("A", "B"))
g2cs6b <- plot_grid(g2cs6bgenes, obscn, ncol = 2, labels = c("", "C"))

png("figures/G2CS6COb_panel.png", height = 800, width = 1800)
g2cs6b
dev.off()


## G2CS13CO ####
png("figures/G2CS13CO_gviz.png", height = 600, width = 1480)
gvizPlot(cnvsGR[6], cnvsGR[6])
dev.off()

gpx1 <- plotExpectedVsObservedCounts(ods.res, "ENSG00000233276.5_7", basePlot = TRUE,
                                      groups = factor(ifelse(ods.res$SampleID == "G2CS13CO", 
                                                             "G2CS13CO", "Population")), 
                                      main = "GPX1 Expression") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))
g2cs13igv <- ggdraw() + draw_image("figures/IGV_figures/G2CS13CO_chr3dup.png")
g2cs13gviz <- ggdraw() + draw_image("figures/G2CS13CO_gviz.png")

g2cs13genes <- plot_grid(g2cs13igv, g2cs13gviz, nrow = 2, labels = c("A", "B"))
g2cs13 <- plot_grid(g2cs13genes, gpx1, ncol = 2, labels = c("", "C"))

png("figures/G2CS13CO_panel.png", height = 800, width = 1800)
g2cs13
dev.off()

## G4CS13CO ####
png("figures/G4CS13CO_gviz.png", height = 400, width = 960)
gvizPlot(cnvsGR[9], cnvsGR[9])
dev.off()

lims2 <- plotExpectedVsObservedCounts(ods.res, "ENSG00000072163.19_5", basePlot = TRUE,
                                     groups = factor(ifelse(ods.res$SampleID == "G4CS13CO", 
                                                            "G4CS13CO", "Population")), 
                                     main = "LIMS2 Expression") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))
g4cs13igv <- ggdraw() + draw_image("figures/IGV_figures/G4CS13CO.chr2dup.png")
g4cs13gviz <- ggdraw() + draw_image("figures/G4CS13CO_gviz.png")

g4cs13genes <- plot_grid(g4cs13igv, g4cs13gviz, nrow = 2, labels = c("A", "B"))
g4cs13 <- plot_grid(g4cs13genes, lims2, ncol = 2, labels = c("", "C"))

png("figures/G4CS13CO_panel.png", height = 800, width = 1800)
g4cs13
dev.off()


## G4CS3CO ####
png("figures/G4CS3CO_gviz.png", height = 600, width = 1480)
gvizPlot(cnvsGR[7], cnvsGR[7])
dev.off()

setd1b <- plotExpectedVsObservedCounts(ods.res, "ENSG00000139718.10_4", basePlot = TRUE,
                                     groups = factor(ifelse(ods.res$SampleID == "G4CS3CO", 
                                                            "G4CS3CO", "Population")), 
                                     main = "SETD1B Expression") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))
g4cs3igv <- ggdraw() + draw_image("figures/IGV_figures/G4CS3CO.chr12dup.png")
g4cs3gviz <- ggdraw() + draw_image("figures/G4CS3CO_gviz.png")

g4cs3genes <- plot_grid(g4cs3igv, g4cs3gviz, nrow = 2, labels = c("A", "B"))
g4cs3 <- plot_grid(g4cs3genes, setd1b, ncol = 2, labels = c("", "C"))

png("figures/G4CS3CO_panel.png", height = 800, width = 1800)
g4cs3
dev.off()



## G2CS10CO ####
png("figures/G2CS10CO_gviz.png", height = 400, width = 960)
gvizPlot(cnvsGR[5], cnvsGR[5])
dev.off()

maff <- plotExpectedVsObservedCounts(ods.res, "ENSG00000185022.12_7", basePlot = TRUE,
                                     groups = factor(ifelse(ods.res$SampleID == "G2CS10CO", 
                                                            "G2CS10CO", "Population")), 
                                     main = "MAFF Expression") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))
tmem <- plotExpectedVsObservedCounts(ods.res, "ENSG00000198792.13_7", basePlot = TRUE,
                                     groups = factor(ifelse(ods.res$SampleID == "G2CS10CO", 
                                                            "G2CS10CO", "Population")), 
                                     main = "TMEM184B Expression") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))

csnk1e <- plotExpectedVsObservedCounts(ods.res, "ENSG00000213923.13_8", basePlot = TRUE,
                                     groups = factor(ifelse(ods.res$SampleID == "G2CS10CO", 
                                                            "G2CS10CO", "Population")), 
                                     main = "CSNK1E Expression") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))

g2cs10igv <- ggdraw() + draw_image("figures/IGV_figures/G2CS10CO.chr22dup.png")
g2cs10gviz <- ggdraw() + draw_image("figures/G2CS10CO_gviz.png")

g2cs10genes <- plot_grid(g2cs10igv, g2cs10gviz, nrow = 2, labels = c("A", "B"))
g2cs10 <- plot_grid(g2cs10genes, maff, tmem, csnk1e, nrow = 2, ncol = 2,  labels = c("", "C", "D", "E"))

png("figures/G2CS10CO_panel.png", height = 1600, width = 1800)
g2cs10
dev.off()

## Not significant
png("figures/G2CS10COb_gviz.png", height = 400, width = 960)
gvizPlot(cnvsGR[4], cnvsGR[4])
dev.off()

gne <- plotExpectedVsObservedCounts(ods.res, "ENSG00000159921.16_5", basePlot = TRUE,
                                     groups = factor(ifelse(ods.res$SampleID == "G2CS10CO", 
                                                            "G2CS10CO", "Population")), 
                                     main = "GNE Expression") +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 20),
        axis.text  = element_text(size = 16),
        legend.text = element_text(size = 20))

g2cs10bigv <- ggdraw() + draw_image("figures/IGV_figures/G2CS10CO.chr10dup.png")
g2cs10bgviz <- ggdraw() + draw_image("figures/G2CS10COb_gviz.png")

g2cs10bgenes <- plot_grid(g2cs10bigv, g2cs10bgviz, nrow = 2, labels = c("A", "B"))
g2cs10b <- plot_grid(g2cs10bgenes, gne, ncol = 2, labels = c("", "C"))


png("figures/G2CS10COb_panel.png", height = 800, width = 1800)
g2cs10b
dev.off()


## SVs ####
svResL <- lapply(seq_len(length(svGR)), function(x) getOutlierDist(svGR[x]))
svPlots <- lapply(svResL, plotDistribution)

g2cs4 <- svResL[[2]]
delGenes <- svResL[[3]]$geneID[svResL[[3]]$class == "Genes in SV"]
g2cs4$class <- ifelse(g2cs4$geneID %in% delGenes, "Genes in chr7pter deletion", 
                  ifelse(g2cs4$class == "Genes in SV", "Genes in chr3pter duplication", "Genes outside SVs"))
g2cs4p <- plotDistribution(g2cs4) +
  scale_fill_manual(values = c("blue", "red", "grey")) +
  ggtitle("Expression of G2CS4CO")

g2cs2p <- svPlots[[1]] +
  scale_fill_manual(values = c("red", "grey"), labels = c("Genes in chr9qter deletion", "Genes outside SV")) +
  ggtitle("Expression of G2CS2CO")

g2cs6p <- svPlots[[4]] +
  scale_fill_manual(values = c("red", "grey"), labels = c("Genes in chr17 deletion", "Genes outside SV")) +
  ggtitle("Expression of G2CS6CO")

g2cs8p <- svPlots[[5]] +
  scale_fill_manual(values = c("blue", "grey"), labels = c("Genes in chr15 mosaic trisomy", "Genes outside SV")) +
  ggtitle("Expression of G2CS8CO")

svsplot_panel <- plot_grid(g2cs4p, g2cs8p, g2cs2p, g2cs6p, ncol = 2, nrow = 2, labels = LETTERS[1:4])
png("figures/SV_expression_panel.png", height = 450, width = 900)
svsplot_panel
dev.off()


lapply(svResL[c(1, 4:5)], function(x) {
  x$class <- factor(x$class, levels = c("Genes outside SV", "Genes in SV"))
  summary(lm(l2fc ~ class, x))
})

g2cs4$class <- factor(g2cs4$class, levels = c("Genes outside SVs", "Genes in chr7pter deletion",  "Genes in chr3pter duplication"))
summary(lm(l2fc ~ class, g2cs4))
