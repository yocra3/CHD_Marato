## Script to generate input file for sarek

## Create file matching fastq and WGS sample
for f in `ls WGS/fastq/*R1_001.fastq.gz`
do
p1=$(basename $f)
samp=$(echo $p1 | cut -d'_' -f 1)
group=$(echo $p1 | cut -d'_' -f 2)
lane=$(echo $p1 | cut -d'_' -f 3)
p2="WGS/fastq/${samp}_${group}_${lane}_R2_001.fastq.gz"

echo -e "$samp\tXX\t0\t$lane\t$f\t$p2"
done > fastqPaths.tsv

## Create mapping file
join -j 1 -o 1.2,2.2,2.3,1.4,2.4,2.5,2.6 -t $'\t' <(paste wgsBatch_MaratoID.map wgsBatch_MaratoID.map | grep -Ev 'AUT|COHO' | sort -k1) <(sort -k1 fastqPaths.tsv) > inputWGS.tsv

cut -f1 inputWGS.tsv | sort -u | split - ids -l 10 -d 
for id in ids*
do
 grep -f $id inputWGS.tsv > ${id}.tsv
 rm $id
done
