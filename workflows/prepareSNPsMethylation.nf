/*
 * Extract SNPs present in Illumina 450K from WGS data
 */

// Override baseDir from nextflow repository to current directory
baseDir = "$PWD"
date = java.time.LocalDate.now()
params.outdir = "results/methylation/SNPs/${date}"
params.sampleAnnot = baseDir
params.inFold = baseDir
params.vcfref = "/home/SHARED/DATA/REFERENCES/hg19/hg19/All_20180418.vcf.gz" //Downloaded from ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz

// Select containers
container_ubuntu = 'yocra3/ubuntu_genomicutils:0.99.1'
container_R = 'yocra3/rsession_chd_marato:release-1.2.2'

/*
 * Convert parameters with strings to files
 */
vcfref = file(params.vcfref)
sampleAnnot = file(params.sampleAnnot)

// Store workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $nextflow.timestamp
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

vcfs = Channel.fromPath("${params.inFold}/*.vcf")


// Divide VCF variants in SNPs
process convertToSNP {

  container container_ubuntu

  wordEndsVcf = ~/vcf$\w+/

  input:
  file (vcf) from vcfs

  output:
  file("decomposed.vcf") into decom

  """
  /home/SHARED/Software/vt/vt  decompose_blocksub $vcf -o decomposed.vcf
  """
}

// Get list of SNPs from meffil
process getSNPsIDs {

  container container_R

  output:
  file "snp-names.txt" into snpIDs

  """
  Rscript -e 'writeLines(meffil::meffil.snp.names(), con = "snp-names.txt")'
  """
}

// List positions of Illumina SNPs based on 1000 Genomes
process getSNPsPos {

  container container_ubuntu

  input:
  file vcfref from vcfref
  file snpNames from snpIDs

  output:
  file "snps_pos.txt" into snpPos
  file "1000KG.pos.recode.vcf" into vcfRefPos

  """
  vcftools --gzvcf ${vcfref} --snps ${snpNames} --recode --out 1000KG.pos
  zgrep -v "#" 1000KG.pos.recode.vcf | cut -f1,2 | awk '{print "chr"\$0}' > snps_pos.txt
  """
}

/*
* Select SNPS present in Illumina by position
*/
process extractSNPs {

  container container_ubuntu

  input:
  file (vcfdec) from decom
  file snpPos from snpPos

  output:
  file("vcffilt.gz") into filtvcfs
  file("vcffilt.gz.tbi") into filtvcfsidx

  """
  bgzip $vcfdec
  tabix -p vcf ${vcfdec}.gz
  bcftools view ${vcfdec}.gz -R $snpPos -Oz -o vcffilt.gz
  tabix -p vcf vcffilt.gz
  """
}

/*
* Merge VCFs in one file
*/
process mergeVCFs {

  container container_ubuntu

  input:
  file("vcffilt*.vcf.gz") from filtvcfs.toList()
  file("vcffilt*.vcf.gz.tbi") from filtvcfsidx.toList()

  output:
  file "merged.vcf.gz" into mergedVcf
  file "merged.vcf.gz.tbi" into mergedVcfidx

  """
  bcftools merge --missing-to-ref vcffilt*.vcf.gz -Oz -o "merged.vcf.gz"
  tabix -p vcf "merged.vcf.gz"
  """
}

/*
* Annotate VCFs with rs id
*/
process annotateVCF {

  container container_ubuntu

  input:
  file(mergedVcf) from mergedVcf
  file(mergedVcfidx) from mergedVcfidx
  file(vcfRefPos) from vcfRefPos

  output:
  file("mergedVcf.rs.gz") into annotVcf

  """
  awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}' ${vcfRefPos} > 1000KG.pos.chr.vcf ## Add chr to name to match exome vcfs
  bgzip 1000KG.pos.chr.vcf
  tabix -p vcf 1000KG.pos.chr.vcf.gz
  bcftools annotate --annotations 1000KG.pos.chr.vcf.gz \
    --columns ID --threads 20 --output-type z\
    --output mergedVcf.rs.gz ${mergedVcf}
  """
}

/*
* Generate genotype file suitable for meffil
*/
process processPhenotypes {

    container container_ubuntu

   	publishDir params.outdir, mode: 'copy'

   	if ( params.version != null ){
   		publishDir "results/methylation/SNPs/${params.version}", mode: 'copy'
   	}

   	input:
    file(annotVcf) from annotVcf
    file (sampleTab) from sampleAnnot
   	val logText from "$workflowInfo"

    output:
    file 'genos.raw' into genos
    file 'log.txt' into logCh

    """
    cut -f23,30 ${sampleTab} | tail -n +2 | awk '{print \$1, \$1, \$2, \$2}' > sample.map
    sed -i 's/_//g' sample.map
    plink --vcf ${annotVcf} --make-bed --out genos
    plink -bfile genos --update-ids sample.map --recodeA --out genos
    echo "$logText" > log.txt
    """
}
