## Script to generate input file pipeline CNVs

## Create file matching fastq and WGS sample
for f in `ls results/Alignments/*.recal.bam`
do
p1=$(basename $f)
samp=$(echo $p1 | cut -d'.' -f 1)
echo -e "$samp\t$f\tresults/VariantCalling/SNV/Strelka_${samp}_variants.vcf.gz"
done > temp.txt

## Add Sex
join -t $'\t' -j 1 -o 1.1,1.2,1.3,2.2  <(sort -k1 temp.txt)  <(sort -k1 data/pheno.tab) > cnvCallingPaths.tsv

## Set sex as F and M
sed -i 's/Female/F/g' cnvCallingPaths.tsv
sed -i 's/Male/M/g' cnvCallingPaths.tsv
