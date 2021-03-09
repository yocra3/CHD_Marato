#!/bin/bash
#SBATCH -J runNextflow
#SBATCH -N 1
#SBATCH -o nextflow.out # STDOUT
#SBATCH -e nextflow.err # STDERR
#SBATCH --mem=1g
#SBATCH --tasks-per-node=1

module load nextflow
module load Singularity

## Run nextflow in batches and copy results to projects
### Results and work folder should be removed after each batch
### Batch 1
nextflow run nf-core/sarek -profile singularity -resume --tools 'Manta,Strelka' --genome GRCh37 --igenomes_ignore False --igenomes_base /homes/users/cruizg/projects/references/ --species 'homo_sapiens' --input ids00.tsv --max_memory '60.GB' --no_gatk_spark
rsync -azvh results/ ~/projects/CHD_MARATO/batch1
rm -rf work/ results/

### Batch 2
nextflow run nf-core/sarek -profile singularity -resume --tools 'Manta,Strelka' --genome GRCh37 --igenomes_ignore False --igenomes_base /homes/users/cruizg/projects/references/ --species 'homo_sapiens' --input ids01.tsv --max_memory '60.GB' --no_gatk_spark
rsync -azvh results/ ~/projects/CHD_MARATO/batch2
rm -rf work/ results/

### Batch 3
nextflow run nf-core/sarek -profile singularity -resume --tools 'Manta,Strelka' --genome GRCh37 --igenomes_ignore False --igenomes_base /homes/users/cruizg/projects/references/ --species 'homo_sapiens' --input ids02.tsv --max_memory '60.GB' --no_gatk_spark
rsync -azvh results/ ~/projects/CHD_MARATO/batch3
rm -rf work/ results/

### Batch 4
nextflow run nf-core/sarek -profile singularity -resume --tools 'Manta,Strelka' --genome GRCh37 --igenomes_ignore False --igenomes_base /homes/users/cruizg/projects/references/ --species 'homo_sapiens' --input ids03.tsv --max_memory '60.GB' --no_gatk_spark
rsync -azvh results/ ~/projects/CHD_MARATO/batch4
rm -rf work/ results/

### Batch 5
nextflow run nf-core/sarek -profile singularity -resume --tools 'Manta,Strelka' --genome GRCh37 --igenomes_ignore False --igenomes_base /homes/users/cruizg/projects/references/ --species 'homo_sapiens' --input ids04_filt.tsv --max_memory '60.GB' --no_gatk_spark
rsync -azvh results/ ~/projects/CHD_MARATO/batch5
rm -rf work/ results/

