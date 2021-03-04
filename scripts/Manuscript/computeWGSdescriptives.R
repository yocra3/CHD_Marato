#'#################################################################################
#'#################################################################################
#' Compute WGS descriptives
#' - Use median and mean coverage 
#' - Compute IQR
#'#################################################################################
#'#################################################################################

multiQCpath <- "reports/multiQC/multiqc_data/multiqc_general_stats.txt"

tab <- read.table(multiQCpath, sep = "\t", header = TRUE, as.is = TRUE)

## Select bam rows and qualimap columns
tab.bam <- tab[grepl("recal", tab$Sample), grepl("QualiMap", colnames(tab)) | colnames(tab) == "Sample"]

## Add status column
tab.bam$Status <- ifelse(grepl("G3", tab.bam$Sample), "Control", "Case")

## Remove sample with T21
tab.bam <- subset(tab.bam, Sample != "G3CS9CO.recal")


computeIQR <- function(x)  {
	q <- signif(quantile(x, c(0.25, 0.5, 0.75)), 4)
	paste0(q[2], " (", q[1], "-", q[3], ")")
}	

sel.cols <- paste0("QualiMap_mqc.generalstats.qualimap.", 
				   c("mapped_reads", "median_coverage", "mean_coverage",
				     "10_x_pc", "30_x_pc", "50_x_pc", "general_error_rate"))



sumTab <- data.frame(Cases = sapply(sel.cols, function(col) computeIQR(subset(tab.bam, Status == "Case")[[col]])),
					 Controls = sapply(sel.cols, function(col) computeIQR(subset(tab.bam, Status == "Control")[[col]]))
)
rownames(sumTab) <- c("Number of mapped reads (Million)",
					  "Median coverage at each targeted base (X)",
					  "Mean coverage at each targeted base (X)",
					  "Percentage of bases read at least 10X",
					  "Percentage of bases read at least 30X",
					  "Percentage of bases read at least 50X",
					  "Mean Error Rate")

write.table(sumTab, file = "reports/paper/WGS_sumStats.tsv", sep = "\t", quote = FALSE)
tests <- lapply(sel.cols, function(x) wilcox.test(formula(paste(x, "~ Status")), tab.bam))
	   
## Select variants and bcftools columns
tab.vars <- tab[grepl("Strelka", tab$Sample), grepl("Bcftools", colnames(tab)) | colnames(tab) == "Sample"]
## Add status column
tab.vars$Status <- ifelse(grepl("G3", tab.vars$Sample), "Control", "Case")

## Remove sample with T21
tab.vars <- subset(tab.vars, Sample != "Strelka_G3CS5CO_variants")

computeIQR(tab.vars$Bcftools.Stats_mqc.generalstats.bcftools_stats.number_of_records)
computeIQR(tab.vars$Bcftools.Stats_mqc.generalstats.bcftools_stats.number_of_SNPs)

var.cols <- paste0("Bcftools.Stats_mqc.generalstats.bcftools_stats.number_of_", 
				   c("records", "SNPs"))
lapply(var.cols, function(x) wilcox.test(formula(paste(x, "~ Status")), tab.vars))



