#'#################################################################################
#'#################################################################################
#' Set up server for CHD_MARATO project
#'#################################################################################
#'#################################################################################

## Project folder: /home/SHARED/PROJECTS/CHD_MARATO/

# Create folders

mkdir data
ln -s /media/Lacie_1/DATA/Methylation_Marato/ data/methylation

ln -s /media/Lacie_1/DATA/MARATÓ/RNASEQ/ data/RNAseq_fastq

mkdir data/RNAseq

## Add links to WGS bams
mkdir data/WGS/
mkdir data/WGS/BAMS/
ln -s /media/carlos/PORSCHE_2/RAW_DATA/Marato_Fetus/BAMS/*S?.bam data/WGS/BAMS/
ln -s /media/carlos/PORSCHE_1/RAW_DATA/Marato_Fetus/BAMS/*S?.bam data/WGS/BAMS/
ln -s /media/carlos/PORSCHE_1/RAW_DATA/Marato_Fetus/BAMS/*S??.bam data/WGS/BAMS/
ln -s /media/carlos/CANVIO_2/RAW_DATA/Marato_Fetus/BAMS/*S?.bam data/WGS/BAMS/
ln -s /media/carlos/CANVIO_1/RAW_DATA/Marato_Fetus/BAMS/*S?.bam data/WGS/BAMS/

ln -s /media/carlos/PORSCHE_2/RAW_DATA/Marato_Fetus/BAMS/*S?.bam.bai data/WGS/BAMS/
ln -s /media/carlos/PORSCHE_1/RAW_DATA/Marato_Fetus/BAMS/*S?.bam.bai data/WGS/BAMS/
ln -s /media/carlos/PORSCHE_1/RAW_DATA/Marato_Fetus/BAMS/*S??.bam.bai data/WGS/BAMS/
ln -s /media/carlos/CANVIO_2/RAW_DATA/Marato_Fetus/BAMS/*S?.bam.bai data/WGS/BAMS/
ln -s /media/carlos/CANVIO_1/RAW_DATA/Marato_Fetus/BAMS/*S?.bam.bai data/WGS/BAMS/


## Remove bams from other projects
rm data/WGS/BAMS/624603*
rm data/WGS/BAMS/624604*
rm data/WGS/BAMS/62460506_S1.bam
rm data/WGS/BAMS/62460517_S2.bam

## Add links to WGS fastqs
mkdir data/WGS/FASTQ/
ln -s /media/carlos/PORSCHE_1/RAW_DATA/Marato_Fetus/20180104_aeav/*.fastq.gz data/WGS/FASTQ/

### Copy quantification data from dropbox
### Download general sample data to data folder

mkdir results
mkdir results/methylation
mkdir results/methylation/QC_intermediate
mkdir results/methylation/finalQC_files
mkdir results/methylation/SNPs

mkdir scripts
mkdir workflows

## Function to extract VCFs from zip to folder (run in /home/SHARED/PROJECTS/CHD_MARATO/data/WGS/VCFs/)
function extractVCF {
        
    unzip $2/$1.zip -d /home/SHARED/PROJECTS/CHD_MARATO/data/WGS/VCFs/ -x *fastq*
    rm *.xls*
    rm *bam*
    rm *summary*
    rm downloadzip.log 
    vcf=$(basename "`echo $1*.vcf.gz`" .gz)
    gunzip ${vcf}.gz
    bgzip $vcf
    tabix -p vcf ${vcf}.gz
}



