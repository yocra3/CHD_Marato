#'#################################################################################
#'#################################################################################
#' Plot aberrant splicing 
#'#################################################################################
#'#################################################################################

# Load libraries
library(cowplot)
makePlot <- function(path) ggdraw() + draw_image(path)

## G2CS12CO ####
g2cs12igv <- makePlot("figures/IGV_figures/G2CS12CO_YIFsplicing.png")
g2cs12sp <- makePlot("results/RNAseq/FRASER_splicing_aberrations_plots/v2/G2CS12CO.1.png")

g2cs12 <- plot_grid(g2cs12igv, g2cs12sp, nrow = 2, labels = c("A", "B"),
                    rel_heights = c(2, 3))

png("figures/G2CS12CO_panel.png", height = 800, width = 1000)
g2cs12
dev.off()

## G2CS6CO ####
g2cs6igv <- makePlot("figures/IGV_figures/G2CS6CO.ZFAS1splicing.png")
g2cs6sp <- makePlot("results/RNAseq/FRASER_splicing_aberrations_plots/v2/G2CS6CO.2.png")

g2cs6 <- plot_grid(g2cs6igv, g2cs6sp, nrow = 2, labels = c("A", "B"),
                    rel_heights = c(2, 3))

png("figures/G2CS6CO_panel_ZFAS1.png", height = 800, width = 1000)
g2cs6
dev.off()

## G4CS13CO ####
g4cs13igv <- makePlot("figures/IGV_figures/G4CS13CO.RPS6KB2splicing.png")
g4cs13sp <- makePlot("results/RNAseq/FRASER_splicing_aberrations_plots/v2/G4CS13CO.1.png")

g4cs13 <- plot_grid(g4cs13igv, g4cs13sp, nrow = 2, labels = c("A", "B"),
                   rel_heights = c(2, 3))

png("figures/G4CS13CO_panel_RPS6KB2.png", height = 800, width = 1000)
g4cs13
dev.off()