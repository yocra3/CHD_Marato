#'#################################################################################
#'#################################################################################
#' Set up server for CHD_MARATO project
#'#################################################################################
#'#################################################################################

## Project folder: /home/SHARED/PROJECTS/CHD_MARATO/

# Create folders

mkdir data
ln -s /media/Lacie_1/DATA/Methylation_Marato/ data/methylation

mkdir data/RNAseq
### Copy quantification data from dropbox
### Download general sample data to data folder

mkdir results
mkdir results/methylation
mkdir results/methylation/QC_intermediate
mkdir results/methylation/finalQC_files
mkdir results/methylation/SNPs

mkdir scripts