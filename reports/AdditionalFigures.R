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

