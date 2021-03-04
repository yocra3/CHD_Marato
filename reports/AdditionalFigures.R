plotEpi <- function(range, gset, sample){
  gset.f <- subsetByOverlaps(gset, range)
  beta <- data.frame(getBeta(gset.f))
  
  beta$pos <- start(rowRanges(gset.f))
  beta %>%
    gather("Sample", "Methylation", seq_len(ncol(gset))) %>%
    mutate(Sample = ifelse(Sample == "X16A66", "16A66", Sample)) %>%
    ggplot(aes(x = pos, y = Methylation, color = Sample == sample, group = Sample)) +
    geom_point() + 
    geom_line() +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
}

## Plot Chromosome X MAE in females
plotMAE <- function(sample, chrom, maeDF){
  tab <- subset(maeDF, MAE_ID == sample & contig == chrom)
  plot(tab$position, tab$altRatio, main = sample, ylim = c(0, 1), xlab = "Position",
       ylab = "AB")
}
females <- c("G2-CS12-CO", "G1-CS3-CO", "G4-CS3-CO", "G2-CS2-CO", "G2-CS11-CO", "G2-CS13-CO",
             "G3-CS1-CO")
par(mfrow = c(2, 4))
lapply(females, plotMAE, chrom = "X", maeDF = resMAE.filt)

## Plot heatmap chromosome X in females
load("../results/methylation/finalQC_files/v3/gset.filterAnnotatedProbes.Rdata")
beta_X_fem <- getBeta(gset[seqnames(rowRanges(gset)) == "chrX", gset$Sex == "F"])

heatmap(beta_X_fem, scale = "none",col = cm.colors(256))
## Plot Distribution Methylation in females
methX <- gather(data.frame(b), sample, methy) %>%
  mutate(group = ifelse(sample == "X16A66", "16A66", 
                        ifelse(sample == "G2CS2CO", "G2CS2CO", "Rest"))) 
  ggplot(methX, aes(x = methy, color = sample)) +
  geom_density(data = methX[methX$group == "Rest",], color = "grey")
  +
  facet_grid(~  group)
