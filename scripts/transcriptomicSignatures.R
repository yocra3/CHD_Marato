#'#################################################################################
#'#################################################################################
#' Explore transcriptomic signatures 
#' Check whether episignatures correlate wit different transcriptomic signatures
#'#################################################################################
#'#################################################################################

## Load libraries and data ####
library(DESeq2)
library(minfi)

load("results/methylation/finalQC_files/v5/gset.autosomic.Rdata")
load("results/methylation/Episignatures/cohort.training.Rdata")

load("results/RNAseq/QC_files/v3/RNAseq_QCobject.Rdata")
load("results/RNAseq/Bioc_objects/v3/RNAseq_RangedSE_expressedGenes.Rdata")

## Subset GRset to CpGs in episignature and samples with RNAseq
comSamps <- intersect(colnames(gset), colnames(rse))
gset.epi <- gset[svm_comb$feats, comSamps]
rse.com <- rse[, comSamps]

## Select genes around CpGs ####
# vsd.epi <- subsetByOverlaps(vsd, rowRanges(gset.epi))
# pc <- prcomp(t(assay(vsd.epi)))
# apply(pc$x, 2, function(x) summary(lm(x ~ vsd.epi$pathClass) ))
# plot(pc$x[, 5:6], col = factor(vsd.epi$pathClass), pch = 16)
# plotPCA(vsd.epi, intgroup = "pathClass")
## We cannot distinguish groups with transcriptomics

## Associate CpGs to gene expression ###
eQTMs <- lapply(rownames(gset.epi), function(cpg){
  methy.filt <- gset.epi[cpg, ]
  rse2 <- rse.com 
  rse2$methy <- getBeta(methy.filt)
  dds <- DESeqDataSet(rse2, ~ Sex + methy)
  dds <- DESeq(dds, BPPARAM = BiocParallel::MulticoreParam(10))
  res <- results(dds)
  res <- res[names(subsetByOverlaps(rowRanges(rse), rowRanges(methy.filt) + 500e3)), ]
})
names(eQTMs) <- rownames(gset.epi)
eQTMs.df <- Reduce(rbind, eQTMs)
eQTMs.df$padj.comb <- p.adjust(eQTMs.df$pvalue, method = "BH")
eQTMs.df$cpg <- rep(names(eQTMs), sapply(eQTMs, nrow))

## Output genes list for enrichment
selGenes <- unique(rownames(subset(eQTMs.df, padj.comb < 0.05)))
selGenes <- gsub("\\..*", "", selGenes)
write.table(selGenes, file = "results/methylation/Episignatures/eQTM_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

background <- rownames(rse)
background <- gsub("\\..*", "", background)
write.table(background, file = "results/methylation/Episignatures/background_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


## Associate CpGs to gene expression  without controls ###
gset.case <- gset.epi[, gset.epi$pathClass != "Control"]
rse.case <- rse.com[, rse.com$pathClass != "Control"]
eQTMs.case <- lapply(rownames(gset.case), function(cpg){
  methy.filt <- gset.case[cpg, ]
  rse2 <- rse.case 
  rse2$methy <- getBeta(methy.filt)
  dds <- DESeqDataSet(rse2, ~ Sex + methy)
  dds <- DESeq(dds, BPPARAM = BiocParallel::MulticoreParam(10))
  res <- results(dds)
  res <- res[names(subsetByOverlaps(rowRanges(rse), rowRanges(methy.filt) + 500e3)), ]
})
names(eQTMs.case) <- rownames(gset.case)
eQTMs.case.df <- Reduce(rbind, eQTMs.case)
eQTMs.case.df$padj.comb <- p.adjust(eQTMs.case.df$pvalue, method = "BH")
eQTMs.case.df$cpg <- rep(names(eQTMs.case), sapply(eQTMs.case, nrow))


save(eQTMs.case, eQTMs, file = "results/methylation/Episignatures/transcriptomicSignature.list.Rdata")
save(eQTMs.case.df, eQTMs.df, file = "results/methylation/Episignatures/transcriptomicSignature.Rdata")
