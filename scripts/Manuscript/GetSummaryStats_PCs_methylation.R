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
load("results/methylation/finalQC_files/v4/gset.autosomic.Rdata")

## Compute PCs
meth.pcs <- prcomp(t(getBeta(gset)), rank. = 10)

varprop <- round(meth.pcs$sdev**2/sum(meth.pcs$sdev**2)*100, 2)

### Assess individual features
batch.vars <- c("Sex", "Status", "pathClass")
names(batch.vars) <- batch.vars

df.pheno <- meth.pcs$x %>%
	data.frame() %>%
	mutate(Sex = ifelse(gset$Sex == "M", "Male", "Female"),
		   Status = gset$Status,
		   pathClass = gset$pathClass)

computelm <- function(x, y, df) summary(lm(formula(paste(y, "~",  x)), df))

getPvals <- function(lm) {
	fstats <- lm$fstatistic
	pf(fstats[1], fstats[2], fstats[3], lower.tail = FALSE)
}
	
pcs <- paste0("PC", 1:10)
names(pcs) <- pcs
PCsTests <- lapply(batch.vars, function(x) lapply(pcs, function(y) computelm(x, y, df.pheno)))
pvals.df <- sapply(PCsTests, sapply, getPvals)

save(pvals.df, meth.pcs, file = "results/methylation/finalQC_files/v4/meth_10PCs.Rdata")

## Check all PCS 1-8
plots <- lapply(seq(1, 8, 2), function(x){
	ggplot(df.pheno, aes_string(x = paste0("PC", x), y = paste0("PC", x + 1), color = "pathClass")) +
		geom_point() +
		theme_bw() +
		scale_x_continuous(name =  paste0("PC", x, " (", varprop[x], "%)")) +
		scale_y_continuous(name =  paste0("PC", x+1, " (", varprop[x+1], "%)"))
	
})
plot_grid(plotlist = plots, nrow = 2, ncol = 2)

## Plot for manuscript - PCs 3 and 4 (run out of docker)
plot <- df.pheno %>%
	mutate(pathClass = factor(pathClass, 
							  levels = c("Control", "Conotruncal Malformations", "Left heart hypoplasia", "Other malformations"))) %>%
ggplot(aes(x = PC3, y = PC4, color = pathClass)) +
	geom_point() +
	theme_bw() +
	scale_x_continuous(name = paste0("PC3 (", varprop[3], "%)")) +
	scale_y_continuous(name = paste0("PC4 (", varprop[4], "%)")) +
	scale_color_manual(name = "Pathologic Class", 
					   values = c("#000000", "#CC6666", "#9999CC", "#66CC99"))
png(filename = "reports/methylation.pathClass.PCs.png",  width = 1500, height = 1300, res = 300)
plot
dev.off()
