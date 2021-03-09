#'#################################################################################
#'#################################################################################
#' Download data from GEO using GEOQuery
#'#################################################################################
#'#################################################################################

## Load library
library(GEOquery)
library(minfi)

# GSE62629
pheno <- getGEO("GSE62629", GSEMatrix = TRUE)
#files <- getGEOSuppFiles("GSE62629", baseDir = "data/GEO") ## Download data (manually download and decompress)

gzfiles <- dir("data/GEO/GSE62629", pattern = "gz", full.names = TRUE)

methy <- lapply(gzfiles, function(x) read.delim(gzfile(x), header = FALSE, as.is = TRUE))
methmat <- Reduce(cbind, lapply(methy, `[`, 3))
sampNames <- gsub("data/GEO/GSE62629/", "", gzfiles)
colnames(methmat) <- gsub("_.*", "", sampNames)

annot <- methy[[1]][, 1:2]
annot$start <- annot$V2 + 1 ## Add 1 to ensure matching with EPIC annotation
annotGR <- makeGRangesFromDataFrame(annot, seqnames.field = "V1",
									end.field = "start")

gset.GSE62629 <- GenomicRatioSet(gr = annotGR, Beta = methmat, colData = pData(pheno[[1]]))
gset.GSE62629$Group <- ifelse(gset.GSE62629$"disease state:ch1" == "normal", 
                              paste(gset.GSE62629$source_name_ch1, gset.GSE62629$"disease state:ch1"), 
                              gset.GSE62629$"disease state:ch1")
save(gset.GSE62629, file = "data/GEO/GSE62629/GenomicRatioSet.Rdata")