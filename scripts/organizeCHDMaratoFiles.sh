# Script to reorganize files from CHD Marató generated in batches

## Generate final folders
mkdir -p results/Alignments/
mkdir -p results/VariantCalling/SNV/
mkdir results/VariantCalling/Manta/
mkdir results/Reports

## Move files to folders
batches=(batch1 batch2 batch3 batch4 batch5)

### Move bam files
for batch in ${batches[@]}
do
	for fold in $batch/Preprocessing/*
	do
		mv $fold/Recalibrated/* results/Alignments/
	done
done

### Move VCF files
for batch in ${batches[@]}
do
	for fold in $batch/VariantCalling/*
	do
		mv $fold/Strelka/*variants* results/VariantCalling/SNV/
		mv $fold/Manta/* results/VariantCalling/Manta/
		ls $fold/Strelka/*variants* 
		ls $fold/Manta/* 
	done
done

### Move Reports
for batch in ${batches[@]}
do
	mkdir results/Reports/$batch/
	mv $batch/pipeline_info/*.* results/Reports/$batch/
	mv $batch/Reports/* results/Reports/$batch/
done

