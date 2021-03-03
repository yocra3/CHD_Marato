#'#################################################################################
#'#################################################################################
#' Compute CNV calling descriptives
#' - Compute IQR
#'#################################################################################
#'#################################################################################

logs <- dir("results/CNVs/log")
names(logs) <- sapply(strsplit(logs, ".", fixed = TRUE), `[`, 1)

ncnvs <- sapply(logs, function(x) {
	tab <- read.delim(paste0("results/CNVs/log/", x))
	tab$Number[3]
})
ncnvs <- ncnvs[names(ncnvs) != "G3CS5CO"]
q <- signif(quantile(ncnvs, c(0.25, 0.5, 0.75)), 4)
paste0(q[2], " (", q[1], "-", q[3], ")")

Status <- ifelse(grepl("G3", names(ncnvs)), "Control", "Case")
wilcox.test(ncnvs ~ Status)

tab <- data.frame(ncnvs)