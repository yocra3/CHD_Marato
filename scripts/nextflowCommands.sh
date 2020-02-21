#'#################################################################################
#'#################################################################################
#' Nextflow commands used in the project
#' These can be modified in older versions. Check log file in each results folder.
#'#################################################################################
#'#################################################################################

## Prepare phenotypes
nextflow run yocra3/CHD_Marato/workflows/preparePhenotype.nf -with-docker yocra3/rsession_chd_marato:release-1.2.3 --version v1

## Get SNPs from WGS present in methylation data
nextflow run  yocra3/CHD_Marato/workflows/prepareSNPsMethylation.nf --inFold data/ExomeVCFs \
--sampleAnnot data/CHD_marato_sampleSummary.tab --with-docker --version v1

## QC and normalization of methylation data
nextflow run workflows/methylation_QC_normalization.nf --inFold "./data/methylation/" \
--sampleSheet sampleSheet_Roser_EPIC_concatenados.csv \
--phenoPath results/phenotypes/v1/pheno.Rdata --cores 16 \
--genosPath results/methylation/SNPs/v1/genos.raw --methyAnnot data/EPIC.hg19.manifest.rds \
--version v1 -with-docker yocra3/rsession_chd_marato:release-1.2.3
