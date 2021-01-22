#'#################################################################################
#'#################################################################################
#' Get QC stats and plots of methylation data for manuscript
#' - Use autosomic probes after filtering
#' - Use 100k probes to compute the PCs
#'#################################################################################
#'#################################################################################

# Load libraries
library(minfi)
library(meffil)
library(cowplot)
library(ggplot2)
library(dplyr)

# Load data
load("results/methylation/finalQC_files/v3/gset.autosomic.Rdata")

## Compute PCs
meth.pcs <- meffil.methylation.pcs(getBeta(gset), probe.range = 100e3)

### Assess individual features
batch.vars <- c("Sex", "Status", "pathGroup", "SampleBatch")
names(batch.vars) <- batch.vars

df.pheno <- meth.pcs %>%
	data.frame() %>%
	mutate(Sex = ifelse(gset$Sex == "M", "Male", "Female"),
		   Status = gset$Status,
		   pathGroup = ifelse(is.na(gset$pathGroup) | gset$pathGroup %in% c("G1", "G2"), "Conotruncal",
		   				   ifelse(gset$pathGroup == "G4", "LHH", 
		   				   	   ifelse(gset$pathGroup == "G5", "Complex", "Control"))),
		   SampleBatch = gset$SampleBatch)



plots <- lapply(batch.vars, function(x){
	ggplot(df.merge, aes_string(x = "PC1", y = "PC2", color = x)) + 
		geom_point() +
		theme_bw() +
		ggtitle(x)
})
plot_grid(plotlist = plots, nrow = 2, ncol = 2)

PC1.test <-  lapply(batch.vars, function(x) summary(lm(formula(paste("PC1 ~",  x)), df.pheno)))

computelm <- function(x, y, df) summary(lm(formula(paste(y, "~",  x)), df))

getPvals <- function(lm) {
	fstats <- lm$fstatistic
	pf(fstats[1], fstats[2], fstats[3], lower.tail = FALSE)
}
	
pcs <- paste0("PC", 1:6)
names(pcs) <- pcs
PCsTests <- lapply(batch.vars, function(x) lapply(pcs, function(y) computelm(x, y, df.pheno)))
pvals.df <- sapply(PCsTests, sapply, getPvals)