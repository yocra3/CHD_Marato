#!/bin/bash
#SBATCH -J extractFastq
#SBATCH -N 1
#SBATCH -o extractFastq.out # STDOUT
#SBATCH -e extractFastq.err # STDERR
#SBATCH --mem=1g
#SBATCH --tasks-per-node=1

dir=WGS/fastq/

for f in $dir/*.zip
do
  unzip $f *.fastq.g* -d $dir
done
md5sum -c *.md5 > md5.check
