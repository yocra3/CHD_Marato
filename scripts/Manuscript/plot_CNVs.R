#'#################################################################################
#'#################################################################################
#' Plot CNVs without gene expression information
#'#################################################################################
#'#################################################################################

# Load libraries
library(cowplot)

## G1CS2CO ####
g1cs2top <- ggdraw() + draw_image("figures/IGV_figures/G1CS2CO.chr16del.png")
g1cs2bot <- ggdraw() + draw_image("figures/IGV_figures/G1CS2CO.chr22del.png")

g1cs2 <- plot_grid(g1cs2top, g1cs2bot, nrow = 2, labels = c("A", "B"))

png("figures/G1CS2CO_panel.png", height = 800, width = 1000)
g1cs2
dev.off()

## G2CS2CO ####
g2cs2top <- ggdraw() + draw_image("figures/IGV_figures/G2CS2CO.chr9Del.png")
g2cs2bot <- ggdraw() + draw_image("figures/IGV_figures/G2CS2CO.chr9DelBP.png")

g2cs2 <- plot_grid(g2cs2top, g2cs2bot, nrow = 2, labels = c("A", "B"))

png("figures/G2CS2CO_chr9_panel.png", height = 800, width = 1000)
g2cs2
dev.off()

## G2CS3CO ####
g2cs3top <- ggdraw() + draw_image("figures/IGV_figures/G2CS3CO.chr1dup.png")
g2cs3bot <- ggdraw() + draw_image("figures/IGV_figures/G2CS3CO_chr11dup.png")

g2cs3 <- plot_grid(g2cs3top, g2cs3bot, nrow = 2, labels = c("A", "B"))

png("figures/G2CS3CO_panel.png", height = 800, width = 1000)
g2cs3
dev.off()

## G2CS4CO ####
makePlot <- function(path) ggdraw() + draw_image(path)
g2cs4ch3 <- makePlot("figures/IGV_figures/G2CS4CO.chr3Translocation.png")
g2cs4ch3bp <- makePlot("figures/IGV_figures/G2CS4CO.chr3TranslocationBP.png")
g2cs4ch7 <- makePlot("figures/IGV_figures/G2CS4CO.chr7Translocation.png")
g2cs4ch7bp <- makePlot("figures/IGV_figures/G2CS4CO.chr7TranslocationBP.png")

g2cs4 <- plot_grid(g2cs4ch3, g2cs4ch3bp, g2cs4ch7, g2cs4ch7bp, nrow = 2, 
                   labels = LETTERS[1:4])

png("figures/G2CS4CO_panel.png", height = 800, width = 2000)
g2cs4
dev.off()

## G4CS7CO ####
g4cs7top <- makePlot("figures/IGV_figures/G4CS7CO.chr3dup.png")
g4cs7bot <- makePlot("figures/IGV_figures/G4CS7CO.chr17del.png")

g4cs7 <- plot_grid(g4cs7top, g4cs7bot, nrow = 2, labels = LETTERS[1:2])

png("figures/G4CS7CO_panel.png", height = 800, width = 1000)
g4cs7
dev.off()
