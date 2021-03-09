#'#################################################################################
#'#################################################################################
#' Nextflow commands used in the project
#' These can be modified in older versions. Check log file in each results folder.
#'#################################################################################
#'#################################################################################

## Prepare phenotypes
nextflow run yocra3/CHD_Marato/workflows/preparePhenotype.nf \
-with-docker yocra3/rsession_chd_marato:release-1.2.4 --version v2

## Prepare phenotypes (v3 - after classifying samples with other analysis)
nextflow run workflows/preparePhenotype.nf \
-with-docker yocra3/rsession_chd_marato:release-1.2.5 --version v4


## Get SNPs from WGS present in methylation data
nextflow run  yocra3/CHD_Marato/workflows/prepareSNPsMethylation.nf --inFold data/ExomeVCFs \
--sampleAnnot data/CHD_marato_sampleSummary.tab -with-docker --version v1

## QC and normalization of methylation data
nextflow workflows/methylation_QC_normalization.nf --inFold "./data/methylation/" \
--sampleSheet sampleSheet_Roser_EPIC_concatenados.csv \
--phenoPath results/phenotypes/v4/pheno.Rdata --cores 9 \
--genosPath results/methylation/SNPs/v1/genos.raw --methyAnnot data/EPIC.hg19.manifest.rds \
--version v5 -with-docker yocra3/rsession_chd_marato:release-1.2.5 -resume

## Epimutations detection
nextflow run yocra3/CHD_Marato/workflows/detectEpiMutations.nf --sampleAnnot results/phenotypes/v3/pheno.tab \
--gset results/methylation/finalQC_files/v5/gset.autosomic.Rdata --version v5 \
-with-docker yocra3/rsession_chd_marato:release-1.2.5 -resume

## RNAseq quantification
nextflow run yocra3/CHD_Marato/workflows/RNAseqQuantification.nf --inFold data/RNAseq_fastq/ \
-with-docker yocra3/ubuntu_genomicutils:release-0.99.3 --version v2 --cores 15

## Create RNAseq objects and QC
nextflow run workflows/RNAseq_QC.nf --countsPath results/RNAseq/quantification/v2/geneQuantification.txt \
--phenoPath results/phenotypes/v3/pheno.Rdata -with-docker yocra3/rsession_chd_marato:release-1.2.5  \
--version v3

## Run OUTRIDER
nextflow run workflows/RNAseq_aberrations_OUTRIDER.nf \
--rangedSE results/RNAseq/Bioc_objects/v3/RNAseq_RangedSE_autosomicGenes.Rdata \
-with-docker yocra3/rsession_chd_marato:release-1.2.5  \
--version v3

## RNAseq variant calling
nextflow run yocra3/CHD_Marato/workflows/RNAseq_VariantCalling.nf \
--inFold results/RNAseq/alignment/v2/ --version v1

## Run splicing aberrations FRASER
nextflow run yocra3/CHD_Marato/workflows/RNAseq_splice_aberrations_FRASER.nf \
--bamFold results/RNAseq/alignment/v2/ --phenoPath results/phenotypes/v2/pheno.Rdata \
--version v1 --cores 15

## CNV calling
nextflow run yocra3/CHD_Marato/workflows/callCNVs.nf --bams "data/WGS/BAMS/*.bam" \
--genome "hg19" --sampleAnnot data/CHD_marato_sampleSummary.tab --version v1

## ASE analysis
nextflow run yocra3/CHD_Marato/workflows/RNAseq_AllelicExpression.nf --bamsDir 'results/RNAseq/sortAlignment/v1/*.bam' \
--vcfDir 'data/WGS/VCFs/*.vcf.gz' --mapFile data/sample_WGSbatch_map.tab -resume --version v1


## Repeat MuliQC with all samples in a report
multiqc results/Reports/  --ignore *L00* --outdir reports/multiQC

## Run docker in project folder
docker run -v $PWD:$PWD -it yocra3/rsession_chd_marato:release-1.2.5 bash
