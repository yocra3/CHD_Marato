#'#################################################################################
#'#################################################################################
#' Figure G4CS8CO
#' Make a plot of Epimutations, RNAseq outliers and Genetic variant
#'#################################################################################
#'#################################################################################

## Working directory in project's folder

# Load libraries ####
library(tidyverse)
library(OUTRIDER)
library(cowplot)
library(minfi)
library(Gviz)
library(biomaRt)

# Load files ####
bumps <- read.table("results/methylation/Epimutations/v5/bumpsSel.txt")
load("results/RNAseq/OUTRIDER_aberrations/v3/OUTRIDER_res.Rdata")
load("results/methylation/finalQC_files/v5/gset.autosomic.Rdata")

## Get IGV figure 
igv <- ggdraw() + draw_image("./figures/G4CS8CO_PRDM16.png")

## Epimutation values ####
bump <- subset(bumps, V11 == "G4CS8CO" & V1 == "chr21")
bumpGR <- makeGRangesFromDataFrame(bump, seqnames.field = "V1", start.field = "V2",
                                   end.field = "V3")
mat <- getBeta(subsetByOverlaps(gset, bumpGR))
cpgsGR <- rowRanges(subsetByOverlaps(gset, bumpGR))

## Prepare df for plotting values
betadf <- mat %>%
  data.frame() %>%
  mutate(CpG = rownames(.),
         pos = start(cpgsGR)) %>%
  gather(ID, Beta, -c(45, 46)) %>%
  mutate(Sample = ifelse(ID == "G4CS8CO", "G4CS8CO", "Cohort")) 

means <- rowMeans(mat[, colnames(mat) != "G4CS8CO"])
sd <- rowSds(mat[, colnames(mat) != "G4CS8CO"])
statsdf <- data.frame(means = means, sd = sd, CpG = names(means), pos = start(cpgsGR))


epi_plot <- ggplot(data = statsdf, aes(x = pos)) + 
  geom_line(data = betadf, aes(x = pos, y = Beta, group = ID, color = Sample), linetype = "longdash") + 
  geom_point(data = betadf, aes(x = pos, y = Beta, group = ID, color = Sample)) +
  scale_colour_manual( values = c("black","blue")) +
  geom_ribbon(aes(ymin = means - 2*sd, ymax =  means + 2*sd), fill = "gray39", alpha = 0.4) +
  geom_ribbon(aes(ymin = means - sd, ymax =  means + sd), fill = "gray98", alpha = 0.4) +
  geom_line(aes(y = means), color = "black") +
  geom_point(aes(y = means), show.legend = TRUE) +
  lims(y = c(0,1)) +  
  theme_bw() + 
  labs(x = "Coordinates chr21 (hg19)") + 
  labs(y = "DNA methylation level")

## Epimutation - genetic context ####
ideo <- IdeogramTrack(genome = "hg19", chromosome = "chr21")
basetracks <- GenomeAxisTrack()

bm <- useMart(host = "grch37.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL",
              dataset = "hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr21", 
                                    start = start(bumpGR) - 10e3, 
                                    end = end(bumpGR) + 10e3, name = "ENSEMBL", biomart = bm)
tracks_Highlight <- Gviz::HighlightTrack(trackList = list(basetracks, biomTrack),
                                         start = start(bumpGR), end = end(bumpGR),
                                         chromosome = "chr21",
                                         col = "#7EA577", fill = "#C6D7C3",
                                         alpha = 0.4,
                                         inBackground = FALSE)
png("figures/G4CS8CO_epimutation_genes.png", height = 200)
plotTracks(c(ideo, tracks_Highlight),
           transcriptAnnotation = "SYMBOL")
dev.off()

epi_context <- ggdraw() + draw_image("figures/G4CS8CO_epimutation_genes.png")


## Expression outliers ####
# COL6A1
col6 <- plotExpectedVsObservedCounts(ods.res, "ENSG00000142156.15_3", basePlot = TRUE,
                   groups = factor(ifelse(ods.res$SampleID == "G4CS8CO", 
                                          "G4CS8CO", "Population")), 
                   main = "COL6A1 Expression") +
  theme(legend.title = element_blank())
                             
prdm16 <- plotExpectedVsObservedCounts(ods.res, "ENSG00000142611.17_4", basePlot = TRUE,
                                     groups = factor(ifelse(ods.res$SampleID == "G4CS8CO", 
                                                            "G4CS8CO", "Population")), 
                                     main = "PRDM16 Expression") +
  theme(legend.title = element_blank())


upper <- plot_grid(igv, prdm16, ncol = 2, labels = c("A", "B"))
epi <- plot_grid(epi_context, epi_plot, nrow = 2)
down <- plot_grid(epi, col6, ncol = 2, labels = c("C", "D"))

png("figures/G4CS8CO_panel.png", height = 600, width = 900)
plot_grid(upper, down, nrow = 2)
dev.off()
