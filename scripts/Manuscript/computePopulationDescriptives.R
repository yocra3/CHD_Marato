#'#################################################################################
#'#################################################################################
#' Compute population descriptives
#' - Rely on code from preparePhenotypesv2.R
#' - Test difference in gestational age between cases and controls
#'#################################################################################
#'#################################################################################

## Load libraries ####
library(readxl)
library(dplyr)

inFold <- "."
### Load and adapt samples data
raw <- read_excel(paste0(inFold, "/data/mostresMARATO_genomes_transcriptoma_Metiloma2.xlsx"), 
				  sheet = "Metiloma_plantilla")

## Remove rows without sample id
raw <- raw[!is.na(raw[, "Name...1"]), ]

## Select desired variables
selvars <- c("Name...1", "SG", "Sex")

## Correct gestional age
pheno <- raw[, selvars]
pheno$GestAge <- pheno$SG
pheno$GestAge[pheno$GestAge == "?"] <- NA
pheno$GestAge <- as.numeric(pheno$GestAge)

## Correct sex
pheno$Sex[!pheno$Sex %in% c("Male", "Female")] <- NA

## Simplify SampleID and adhere to conventions
### Change spaces with a dash
pheno$SampleID <- pheno$`Name...1`
pheno$SampleID <- gsub(" |-|Mallorca|_", "", pheno$SampleID)

## Add Pathologic groups
pheno$pathGroup <- sapply(strsplit(pheno$SampleID, "CS"), `[`, 1)
pheno$pathGroup[!pheno$pathGroup %in% paste0("G", 1:5)] <- NA
pheno$pathGroup[pheno$pathGroup == "G3"] <- "Control"

## Add Status
pheno$Status <- factor(ifelse(is.na(pheno$pathGroup) | pheno$pathGroup != "Control", "Case", "Control"))

## Remove sample with T21
pheno <- subset(pheno, SampleID != "G3CS5CO")

## Compare gestational age between the groups
wilcox.test(GestAge ~ Status, pheno)
