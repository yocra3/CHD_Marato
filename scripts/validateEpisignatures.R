#'#################################################################################
#'#################################################################################
#' Validate episignatures with external data
#' 1. Run SVM in external dataset (GSE62629)
#' 2. Check overlap between features selected and differential features from the literature
#'#################################################################################
#'#################################################################################

## Load libraries
library(minfi)
library(e1071)
library(tidyr)
library(dplyr)
library(sva)

## Load files
load("results/methylation/finalQC_files/v4/gset.autosomic.Rdata")
load("results/methylation/Episignatures/cohort.training.Rdata")
load("data/GEO/GSE62629/GenomicRatioSet.Rdata")

## Test SVM in GSE62629 ####
### Match GSE62629 probes with EPIC probes
over <- findOverlaps(rowRanges(gset), rowRanges(gset.GSE62629))

gset.GSE62629.filt <- gset.GSE62629[to(over), ]
rownames(gset.GSE62629.filt) <- rownames(gset)[from(over)]
save(gset.GSE62629.filt, file = "results/methylation/Episignatures/GSE62629.GenomicRatioSet.EPICcpgs.Rdata")
pred.GSE62629 <-  predict(svm_all$svm, t(getBeta(gset.GSE62629.filt[svm_all$feats, ])))


mat <- cbind(getBeta(gset.GSE62629.filt[svm_all$feats, ]),
			 getBeta(gset[svm_all$feats, ]))

pc_comb <- prcomp(t(mat))
pc_GSE <- prcomp(t(getBeta(gset.GSE62629.filt)))
pc_GSE_svm <- prcomp(t(getBeta(gset.GSE62629.filt[svm_all$feats, ])))
save(pred.GSE62629, pc_comb, pc_GSE, pc_GSE_svm, file = "results/methylation/Episignatures/GSE62629.results.Rdata")
