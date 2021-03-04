#'#################################################################################
#'#################################################################################
#' Compute CNV calling descriptives
#' - Compute IQR
#'#################################################################################
#'#################################################################################

logs <- dir("results/Mosaics/SNVs/logs/prioritization/")
names(logs) <- sapply(strsplit(logs, ".", fixed = TRUE), `[`, 1)

nmosaic <- sapply(logs, function(x) {
	tab <- read.delim(paste0("results/Mosaics/SNVs/logs/prioritization/", x))
	as.numeric(strsplit(tab[1,], " ")[[1]][4])
})
nmosaic <- nmosaic[names(nmosaic) != "G3CS5CO"]
q <- signif(quantile(nmosaic, c(0.25, 0.5, 0.75)), 4)
paste0(q[2], " (", q[1], "-", q[3], ")")

Status <- ifelse(grepl("G3", names(nmosaic)), "Control", "Case")
wilcox.test(nmosaic ~ Status)
