#'#################################################################################
#'#################################################################################
#' Prepare phenotype data - version 2
#' This code process original phenotype data to create required variables for the 
#' analysis
#' Important decisions:
#' - Use new molecular diagnosis file. This file was generated after running
#' all other analysis
#' - Create a batch variable
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
inFold <- args[1]


## Load libraries ####
library(readxl)
library(dplyr)

### Load and adapt samples data
raw <- read_excel(paste0(inFold, "/data/mostresMARATO_genomes_transcriptoma_Metiloma2.xlsx"), 
                      sheet = "Metiloma_plantilla")
### Load genetic data
genetics <- read_excel(paste0(inFold, "/data/Datos_Pacientes.xlsx"), 
                      sheet = "SampleSummary")

## Load samples sheets
batch1 <- read_excel(paste0(inFold, "/data/Muestras2-2016_Datosentregar.xlsx"), 
                     sheet = "Datos entregados", skip = 10)
batch2 <- read_excel(paste0(inFold, "/data/Muestras10-2018_datosentregar.xlsx"), skip = 8)


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
pheno$pathGroup[!pheno$pathGroup %in% paste0("G", 1:5)] <- "G5"
pheno$pathGroup[pheno$pathGroup == "G3"] <- "Control"

## Add Pathologic groups
pheno$pathClass <- ifelse(pheno$pathGroup == "G2", "Conotruncal Malformations",
						  ifelse(pheno$pathGroup == "Control", "Control", 
						  	   ifelse(pheno$pathGroup == "G4", "Left heart hypoplasia", "Other malformations")))


## Add Status
pheno$Status <- factor(ifelse(is.na(pheno$pathGroup) | pheno$pathGroup != "Control", "Case", "Control"))

# Complete phenotypes with genetic info
## Correct IDs
genetics$SampleID <- gsub(" |-|_Mallorca", "", genetics$ID)
genetics$Del22q11 <- genetics$`22q11.2del`

## Add missing sex to phenotypes
genSex <- data.frame(genetics)
rownames(genSex) <- genSex$SampleID
pheno$Sex[is.na(pheno$Sex)] <- genSex[pheno$SampleID[is.na(pheno$Sex)], "Sex"]

## Merge pheno with genetic summary
pheno <- genetics %>% 
  select(Rearrangements, MolecularCause, N_genes, SampleID, Del22q11) %>%
  right_join(pheno, by = "SampleID")

## Add recruitment bacths
batch1$SampleID <- gsub("-| ", "", batch1$Etiqueta)
batch2$SampleID <- gsub("-| ", "", batch2$Etiqueta)

pheno$SampleBatch <- ifelse(pheno$SampleID %in% batch1$SampleID, "VHIR1",
                            ifelse(pheno$SampleID %in% batch2$SampleID, "VHIR2",
                                   ifelse(!grepl("CO$", pheno$SampleID), "Mallorca", "Lab")))

## Create final object with created variables
finalVars <- c("SampleID", "Sex", "GestAge", "Status", "pathGroup", "pathClass", "MolecularCause", 
               "Rearrangements", "Del22q11", "SampleBatch", "N_genes")
pheno <- pheno[, finalVars]

save(pheno, file = "pheno.Rdata")
write.table(pheno, file = "pheno.tab", sep = "\t", quote = FALSE, row.names = FALSE)

## Create variables codebook
codebook <- data.frame(Variable = finalVars, 
                       Definition = c("Sample identifier", "Sex of the sample obtained from biobank or from genetics",
                                   "Gestional age in weeks", "Case (fetus with CHD) or Control",
                                   "CHD classification from pathologists", 
                                   "If YES, the sample has a potential molecular cause",
                                   "If YES, the sample has a big CNV (>1Mb) or a rearrangement",
                                   "If YES, the sample has del22q11.1",
                                   "Batch when the sample was introduced in the study. VHIR1: first batch from VHIR. VHIR2: second batch from VHIR. Lab: Samples present in the lab. Mallorca: Samples from Mallorca.",
                                   "Number of genes affected. Monogenic: variants in one gene, Oligogenic: several genes affected, Rearrangement: big translocations/CNVs. None: no variants detected."))
write.table(codebook, file = "codebook.tab", sep = "\t", quote = FALSE, row.names = FALSE)
