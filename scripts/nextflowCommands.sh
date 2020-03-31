#'#################################################################################
#'#################################################################################
#' Nextflow commands used in the project
#' These can be modified in older versions. Check log file in each results folder.
#'#################################################################################
#'#################################################################################

## Prepare phenotypes
nextflow run yocra3/CHD_Marato/workflows/preparePhenotype.nf \
-with-docker yocra3/rsession_chd_marato:release-1.2.4 --version v2

## Get SNPs from WGS present in methylation data
nextflow run  yocra3/CHD_Marato/workflows/prepareSNPsMethylation.nf --inFold data/ExomeVCFs \
--sampleAnnot data/CHD_marato_sampleSummary.tab -with-docker --version v1

## QC and normalization of methylation data
nextflow run yocra3/CHD_Marato/workflows/methylation_QC_normalization.nf --inFold "./data/methylation/" \
--sampleSheet sampleSheet_Roser_EPIC_concatenados.csv \
--phenoPath results/phenotypes/v2/pheno.Rdata --cores 9 \
--genosPath results/methylation/SNPs/v1/genos.raw --methyAnnot data/EPIC.hg19.manifest.rds \
--version v2 -with-docker yocra3/rsession_chd_marato:release-1.2.4

## Epimutations detection
nextflow run yocra3/CHD_Marato/workflows/detectEpiMutations.nf --sampleAnnot results/phenotypes/v2/pheno.tab \
--gset results/methylation/finalQC_files/v2/gset.autosomic.Rdata --version v3 \
-with-docker yocra3/rsession_chd_marato:release-1.2.4

## RNAseq quantification
nextflow run yocra3/CHD_Marato/workflows/RNAseqQuantification.nf --inFold data/RNAseq_fastq/ \
-with-docker yocra3/ubuntu_genomicutils:release-0.99.3 --version v2 --cores 15

## Create RNAseq objects and QC
nextflow run yocra3/CHD_Marato/workflows/RNAseq_QC.nf --countsPath results/RNAseq/quantification/v2/geneQuantification.txt \
--phenoPath results/phenotypes/v2/pheno.Rdata -with-docker yocra3/rsession_chd_marato:release-1.2.4  \
--version v2

## Run OUTRIDER
nextflow run yocra3/CHD_Marato/workflows/RNAseq_aberrations_OUTRIDER.nf \
--rangedSE results/RNAseq/Bioc_objects/v2/RNAseq_RangedSE_autosomicGenes.Rdata \
-with-docker yocra3/rsession_chd_marato:release-1.2.4  \
--version v2

## RNAseq variant calling
nextflow run yocra3/CHD_Marato/workflows/RNAseq_VariantCalling.nf \
--inFold results/RNAseq/alignment/v2/ --version v1

## Run splicing aberrations FRASER
nextflow run yocra3/CHD_Marato/workflows/RNAseq_splice_aberrations_FRASER.nf \
--bamFold results/RNAseq/alignment/v2/ --phenoPath results/phenotypes/v2/pheno.Rdata \
--version v1
