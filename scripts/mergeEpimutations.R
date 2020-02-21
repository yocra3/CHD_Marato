#'#################################################################################
#'#################################################################################
#' Merge results from detectEpimutation.R
#'#################################################################################
#'#################################################################################

args <- commandArgs(trailingOnly=TRUE)
l <- lapply(args, function(x){
  load(x)
  bumps.final
})
bumps.combined <- Reduce(rbind, l)
save(bumps.combined, file = "bumps.combined.Rdata")
write.table(bumps.combined, sep = "\t", row.names = FALSE, quote = FALSE, 
            file = "bumps.combined.txt")
