/*
* Run variant calling on RNAseq data
* This workflow starts with bam files generated with STAR (two steps)
*/

// Override baseDir from nextflow repository to current directory
baseDir = "$PWD"
date = java.time.LocalDate.now()

// Default parameters -- Adapt to experiment!!!
params.variantsDir = "results/RNAseq/Variants/"
params.inFold = baseDir
params.fastaRef = "/home/SHARED/DATA/REFERENCES/GRCh37/Sequence/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz" //Compressed file
params.varsRefexome = "/home/SHARED/DATA/REFERENCES/GRCh37/SNP_annotation/gnomad.exomes.r2.1.1.sites.vcf.bgz"
params.varsRefgenome = "/home/SHARED/DATA/REFERENCES/GRCh37/SNP_annotation/gnomad.genomes.r2.1.1.sites.vcf.bgz"
params.version = null

/*
* Convert parameters with strings to files
*/
fastaRef = file("${params.fastaRef}")
varsRefexome = file("${params.varsRefexome}")
varsRefexomeidx = file("${params.varsRefexome}.tbi")
varsRefgenome = file("${params.varsRefgenome}")
varsRefgenomeidx = file("${params.varsRefgenome}.tbi")


// Select containers
container_ubuntu = 'yocra3/ubuntu_genomicutils:release-0.99.4'
container_gatk = 'broadinstitute/gatk:4.1.5.0'


// Store workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

bams = Channel.fromPath("${params.inFold}/*.bam")


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


/*
* Preprocess BAM for GATK analysis. Add read group, sample name and sort bam.
*/
process correctBam{

  container container_ubuntu

  input:
  file(bam) from bams


  output:
  file("sort.bam") into sortbams

  """
  bam=$bam
  samp=\${bam/Aligned.out.bam/}
  java -jar ~/picard.jar AddOrReplaceReadGroups \
        I=$bam \
        O=sort.bam \
        RGID=1 \
        RGLB=HW5CWBBXX \
        RGPL=ILLUMINA \
        RGPU=K00171 \
        RGSM=\$samp \
        SORT_ORDER=coordinate
    """
}

// Mark duplicate reads
process markDuplicates{

  container container_ubuntu

  input:
  file(bam) from sortbams

  output:
  file("dups.bam") into dupbams

  """
  java -jar ~/picard.jar MarkDuplicates \
  I=$bam \
  O=dups.bam \
  M=dup_metrics.txt
  """
}


// Convert CIGAR to DNA convention
process convertCIGAR{

  container container_gatk

  input:
  file(fasta) from fasta
  file(fastaidx) from fastaidx
  file(fastadic) from fastadict
  file(bam) from dupbams

  output:
  file("outbam.bam") into cigarbam

  """
  gatk SplitNCigarReads \
     -R $fasta \
     -I $bam \
     -O outbam.bam
  """

}

cigarbam.into { recalBam; baseBam }

// Recalibrate base score
process recalibrateBase {

  container container_gatk

  input:
  file(bam) from baseBam
  file(varsRefexome)  from varsRefexome
  file(varsRefexomeidx)  from varsRefexomeidx
  file(varsRefgenome)  from varsRefgenome
  file(varsRefgenomeidx)  from varsRefgenomeidx
  file(fasta) from fasta
  file(fastaidx) from fastaidx
  file(fastadic) from fastadict


  output:
  file("recal_data.table") into recalTable

  """
  gatk BaseRecalibrator \
   -I $bam \
   -R $fasta \
   --known-sites $varsRefexome \
   --known-sites $varsRefgenome \
   -O recal_data.table
  """
}

// Apply recalibration from recalibrateBase
process applyRecalibration {

  container container_gatk

  input:
  file(bam) from recalBam
  file(baseTable) from recalTable
  file(fasta) from fasta
  file(fastaidx) from fastaidx
  file(fastadic) from fastadict

  output:
  file("output.bam") into calBams

  """
  gatk ApplyBQSR \
    -R $fasta \
    -I $bam \
    --bqsr-recal-file $baseTable \
    -O output.bam
  """

}

// Run variant calling
process variantCalling {

  container container_gatk

  input:
  file(bam) from calBams
  file(fasta) from fasta
  file(fastaidx) from fastaidx
  file(fastadic) from fastadict

  output:
  file("output.vcf.gz") into vcfs
  file("output.vcf.gz.tbi") into vcfsidx

  """
  gatk --java-options "-Xmx4g" HaplotypeCaller  \
  -R $fasta \
  -I $bam \
  -O output.vcf.gz
  """
}

/*
* Merge VCFs in one file
*/
process mergeVCFs {

  container container_ubuntu

  publishDir "${params.variantsDir}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.variantsDir}/${params.version}", mode: 'copy'
  }

  input:
  file("output*.vcf.gz") from vcfs.toList()
  file("output*.vcf.gz.tbi") from vcfsidx.toList()
  val logText from "$workflowInfo"

  output:
  file "RNAseq.vcf.gz" into mergedVcf
  file "RNAseq.vcf.gz.tbi" into mergedVcfidx

  """
  bcftools merge --missing-to-ref output*.vcf.gz -Oz -o "RNAseq.vcf.gz"
  tabix -p vcf "RNAseq.vcf.gz"
  echo "$logText" > log.txt
  """
}
