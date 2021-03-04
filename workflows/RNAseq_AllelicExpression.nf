/*
* Run variant calling on RNAseq data
* This workflow starts with bam files generated with STAR (two steps)
*/

// Override baseDir from nextflow repository to current directory
baseDir = "$PWD"
date = java.time.LocalDate.now()

// Default parameters -- Adapt to experiment!!!
params.outDir = "results/RNAseq/ASE/"
params.bamsDir = ""
params.vcfDir = ""
params.fastaRef = "/home/SHARED/DATA/REFERENCES/GRCh37/Sequence/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz" //Compressed file
params.varsRefexome = "/home/SHARED/DATA/REFERENCES/GRCh37/SNP_annotation/gnomad.exomes.r2.1.1.sites.vcf.bgz"
params.varsRefgenome = "/home/SHARED/DATA/REFERENCES/GRCh37/SNP_annotation/gnomad.genomes.r2.1.1.sites.vcf.bgz"
params.version = null
params.mapFile = ""

/*
* Convert parameters with strings to files
*/
fastaRef = file("${params.fastaRef}")
varsRefexome = file("${params.varsRefexome}")
varsRefexomeidx = file("${params.varsRefexome}.tbi")
varsRefgenome = file("${params.varsRefgenome}")
varsRefgenomeidx = file("${params.varsRefgenome}.tbi")
mapFile = file("${params.mapFile}")

// Select containers
container_ubuntu = 'yocra3/ubuntu_genomicutils:release-0.99.4'
container_gatk = 'broadinstitute/gatk:4.1.5.0'
container_R = 'yocra3/rsession_chd_marato:release-1.2.5'


// Store workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

bams = Channel.fromPath("${params.bamsDir}").map { file -> tuple(file.simpleName, file) }
vcfs = Channel.fromPath("${params.vcfDir}").map { name -> tuple(name.simpleName.replaceFirst(/_S[0-9]/, ""), file(name)) }
mapChan = Channel.value(mapFile.text).splitCsv(sep: '\t')


// Prepare Fasta reference
process prepareFasta {

  container container_ubuntu

  input:
  file(fasta) from fastaRef

  output:
  file("ref.fasta") into fasta
  file("ref.fasta.fai") into fastaidx

  """
  gzip -cd $fasta > ref.fasta
  samtools faidx ref.fasta
  """
}

process createFastDict {

  container container_gatk

  input:
  file("ref.fasta") from fasta

  output:
  file("ref.dict") into fastadict

  """
  gatk CreateSequenceDictionary -R ref.fasta
  """

}

//
// Preprocess BAM for GATK analysis. Add read group, sample name and sort bam.
//
process correctBam{

  container container_ubuntu

  input:
  set val(sample), file(bam) from bams

  output:
  set val(sample), file("sort.bam") into sortbams

  """
  java -jar ~/picard.jar AddOrReplaceReadGroups \
        I=$bam \
        O=sort.bam \
        RGID=1 \
        RGLB=HW5CWBBXX \
        RGPL=ILLUMINA \
        RGPU=K00171 \
        RGSM=$sample \
        SORT_ORDER=coordinate
    """
}


// Filter VCF
process filterVCF {

  container container_ubuntu

  input:
  set val(sample), file(vcf) from vcfs

  output:
  set val(sample), file("filt.vcf.gz"), file("filt.vcf.gz.tbi") into filtVCF

  """
  bcftools view -m2 -M2 -v snps $vcf | sed 's/chr//' | bgzip > filt.vcf.gz
  tabix -p vcf filt.vcf.gz
  """
}

merge = sortbams.join(mapChan).map{ it -> tuple(it[2], it[0], it[1])}.join(filtVCF)


// Run Allelic Expression
process aseInference {

  container container_gatk

  publishDir "${params.outDir}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.outDir}/${params.version}", mode: 'copy'
  }

  input:
  set val(batch), val(sample), file(bam), file(vcf), file(vcftbi) from merge
  file(fasta) from fasta
  file(fastaidx) from fastaidx
  file(fastadic) from fastadict

  output:
  file("${sample}.table") into outTable

  """
  gatk ASEReadCounter \
  -R $fasta \
  -I $bam \
  -V $vcf \
  -O ${sample}.table
  sed -i "s/\$/\t$sample/" ${sample}.table
  sed -i "s/improperPairs\t$sample/improperPairs\tSampID/" ${sample}.table
  """
}

// Merge ASE table
process mergeTable {

  container container_ubuntu

  input:
  file("*.table") from outTable.collect().ifEmpty([])

  output:
  file("ASE_all.tab") into aseTab

  """
  for t in *.table
  do
    tail -n +2 \$t >> ASE_all.tab
  done
  """
}

// Run MAE analysis
process runMAE{

  container container_R

  publishDir "${params.outDir}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.outDir}/${params.version}", mode: 'copy'
  }

  input:
  file(tab) from aseTab
  val logText from "$workflowInfo"

  output:
  file("tMAE_res.Rdata")
  file("log.txt") into logF

  """
  runtMAE.R $tab
  echo "$logText" > log.txt
  """

}
