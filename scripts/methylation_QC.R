#'#################################################################################
#'#################################################################################
#' Run methylation QC data
#' This code performs quality control to methylation data. It can be reused for other
#' projects but modifying some parts, highlighted with ##
#' Important decisions:
#' - Remove samples not belonging to the project.
#' - Remove samples based on QC
#' - Use values from meffil vignette in all parameters
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
inFold <- args[1]
csvPattern <- args[2]
phenoPath <- args[3]
cores <- args[4]
genosPath <- args[5]

## Load libraries ####
library(meffil)
library(minfi)
library(readxl)
library(ggplot2)
library(dplyr)

## Set number of cores
options(mc.cores = cores)


#'#################################################################################
#'#################################################################################
## Parameters definition
qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.1,
  detectionp.cpgs.threshold             = 0.1, 
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5,
  snp.concordance.threshold             = 0.95,
  sample.genotype.concordance.threshold = 0.8
)
pcs <- 5 

## Batch variables
batch_var <- c("Slide", "Array", "Sex", "GestAge", "Status", "pathGroup", "pathClass", "MolecularCause", "SampleBatch", "Rearrangements", "Del22q11", "N_genes")

#'#################################################################################
#'#################################################################################


## Prepare sample sheet ####
### Load predefined sample sheet
samplesheet <- meffil.read.samplesheet(base = inFold, pattern = csvPattern)
### Load and adapt samples data
load(phenoPath)

## Merge both
samplesheet <- mutate(samplesheet, SampleID = Sample_Name)
combSheet <- left_join(select(samplesheet, -Sex), pheno, by = "SampleID")

## Change sex to M / F
combSheet <- mutate(combSheet, Sex = substring(Sex, 1, 1))

#'#################################################################################
#'#################################################################################
## To be changed in other projects
## Remove IVD sample - not belongs to the project
combSheet <- subset(combSheet, Sample_Name != "IVD")

## Remove G3CS5CO sample - control with T21 
combSheet <- subset(combSheet, Sample_Name != "G3CS5CO")

#'#################################################################################
#'#################################################################################

## Generate QC report
### Load genotypes
genos <- meffil.extract.genotypes(genosPath)

#'#################################################################################
#'#################################################################################
## To be changed in other projects
## Change sample ids from Mallorca to match names
colnames(genos) <- gsub("Mallorca", "", colnames(genos))
#'#################################################################################
#'#################################################################################

genotypes <- genos[, match(combSheet$SampleID, colnames(genos))]

## Load methylation data and QC ####
qc.objects <- meffil.qc(combSheet, verbose = TRUE)
qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters = qc.parameters,
  genotypes = genotypes
  )
save(qc.objects, file = "qc.objects.Rdata")
save(qc.summary, file = "qcsummary.Rdata")


#'#################################################################################
#'#################################################################################
## To be changed in other projects
## Check that with one cycle removing samples is enough
## Remove bad samples based on QC report and rerun QC
outlier <- qc.summary$bad.samples
qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
save(qc.objects, file = "qc.objects.clean.Rdata")

qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
save(qc.summary, file = "qcsummary.clean.Rdata")

#'#################################################################################
#'#################################################################################

## Report filtered samples and probes
write.table(outlier, file = "removed.samples.txt", quote = FALSE, row.names = FALSE,
            sep = "\t")
write.table(qc.summary$bad.cpgs, file = "removed.probes.txt", quote = FALSE, row.names = FALSE,
            sep = "\t")

## Run functional normalization ####
#'#################################################################################
#'#################################################################################
## To be changed in other projects
### Select number PCs (run, see plot and adapt pcs number)
y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot, filename = "pc.fit.pdf", height = 6, width = 6)
#'#################################################################################
#'#################################################################################

norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs = pcs)
save(norm.objects, file = "norm.obj.pc.Rdata")

## Define normalization parameters
norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables = batch_var,
  control.pcs = seq_len(pcs),
  batch.pcs = seq_len(pcs),
  batch.threshold = 0.01
)
norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove = qc.summary$bad.cpgs$name, verbose = TRUE)
save(norm.beta, file = "norm.beta.Rdata")

## Check covariables
beta.pcs <- meffil.methylation.pcs(norm.beta, probe.range = 8e5)
norm.summary <- meffil.normalization.summary(norm.objects, pcs = beta.pcs, parameters = norm.parameters)
save(norm.summary, file = "norm.summary.Rdata")

## Create GenomicRatioSet
rownames(combSheet) <- combSheet$SampleID
gset <- makeGenomicRatioSetFromMatrix(norm.beta, pData = combSheet[colnames(norm.beta), ],
                                      array = "IlluminaHumanMethylationEPIC",
                                      annotation = "ilm10b4.hg19")
save(gset, file = "gset.Rdata")

