qc.objects <- meffil.qc(combSheet, verbose = TRUE)
pvalsmat <- meffil.load.detection.pvalues(qc.objects)

bumps <- read.delim(paste0("results/methylation/Epimutations/v2/bumpsSel.txt"), 
                    header = FALSE)

load("results/methylation/finalQC_files/v1/gset.autosomic.Rdata")

pvalsFilt <- pvalsmat[rownames(gset), ]

colnames(bumps) <- c("chromosome", "start", "end", "value", "area", "cluster", "indexStart",
                     "indexEnd", "L", "clusterL", "Sample", "Anova")

getPvals <- function(bumps, row){
  idx <- bumps[row, "indexStart"]:bumps[row, "indexEnd"]
  pvalsmat[idx, bumps[row, "Sample"]]
  
}

lapply(seq_len(nrow(bumps)), getPvals, bumps = bumps)

range <- GRanges("chr22:19744226-19754855")
