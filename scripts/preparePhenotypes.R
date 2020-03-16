#'#################################################################################
#'#################################################################################
#' Prepare phenotype data
#' This code process original phenotype data to create required variables for the 
#' analysis
#' Important decisions:
#' - Ignore case/control and del22q variables and extract them from identifiers
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
raw_gen <- read_excel(paste0(inFold, "/data/Resumen_Pacientes_modv2.xlsx"), skip = 1)

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
pheno$pathGroup[!pheno$pathGroup %in% paste0("G", 1:5)] <- NA
pheno$pathGroup[pheno$pathGroup == "G3"] <- "Control"

## Add Status
pheno$Status <- factor(ifelse(is.na(pheno$pathGroup) | pheno$pathGroup != "Control", "Case", "Control"))

## Complete phenotypes with genetic info
### Select relevant columns
genvars <- c("ID", "Sexo", "Gen", "Variante", "Priorization")
genetics <- raw_gen[, genvars]

### Propagate id
genetics$SampleID <- genetics$ID[!is.na(genetics$ID)][cumsum(!is.na(genetics$ID))]
genetics$SampleID <- gsub(" |-|Mallorca|_", "", genetics$SampleID)

### Propagate sex
genetics$Sex <- genetics$Sexo[!is.na(genetics$Sexo)][cumsum(!is.na(genetics$Sexo))]

## Add missing sex to phenotypes
genSex <- unique(genetics[, c("SampleID", "Sex")]) %>%
  mutate(Sex = gsub("*", "", Sex, fixed = TRUE)) %>%
  data.frame()
rownames(genSex) <- genSex$SampleID

pheno$Sex[is.na(pheno$Sex)] <- genSex[pheno$SampleID[is.na(pheno$Sex)], "Sex"]

## Add variants identified by Luis
## Genetic cause is TRUE for samples with a detected gene or samples from Group 1 (Di George deletion)
genetics <- subset(genetics, !is.na(Priorization))
genSum <- function(x) {if (length(unique(x)) == 1) unique(x) else paste(x, collapse = "/")}
pheno <- genetics %>% group_by(SampleID) %>%
  summarize(Gene = genSum(Gen)) %>%
  right_join(pheno, by ="SampleID") %>%
  mutate(GeneticCause = ifelse(is.na(Gene) & pathGroup != "G1", FALSE, TRUE))


## Add recruitment bacths
batch1$SampleID <- gsub("-| ", "", batch1$Etiqueta)
batch2$SampleID <- gsub("-| ", "", batch2$Etiqueta)

pheno$SampleBatch <- ifelse(pheno$SampleID %in% batch1$SampleID, "VHIR1",
                            ifelse(pheno$SampleID %in% batch2$SampleID, "VHIR2",
                                   ifelse(!grepl("CO$", pheno$SampleID), "Mallorca", "Lab")))

## Create final object with created variables
finalVars <- c("SampleID", "Sex", "GestAge", "Status", "pathGroup", "GeneticCause", "Pos", "Gene", "SampleBatch")
pheno <- pheno[, finalVars]

save(pheno, file = "pheno.Rdata")
write.table(pheno, file = "pheno.tab", sep = "\t", quote = FALSE, row.names = FALSE)

## Create variables codebook
codebook <- data.frame(Variable = finalVars, 
                       Definition = c("Sample identifier", "Sex of the sample obtained from biobank or from genetics",
                                   "Gestional age in weeks", "Case (fetus with CHD) or Control",
                                   "CHD classification from pathologists", 
                                   "If TRUE, fetus with Di George syndrom or with a pathogenic candidate variant",
                                   "Gene containing the pathogenic variant",
                                   "Batch when the sample was introduced in the study. VHIR1: first batch from VHIR. VHIR2: second batch from VHIR. Lab: Samples present in the lab. Mallorca: Samples from Mallorca."))
write.table(codebook, file = "codebook.tab", sep = "\t", quote = FALSE, row.names = FALSE)
