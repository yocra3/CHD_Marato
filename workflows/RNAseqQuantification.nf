/*
* Run RNAseq quantification in genes, transcripts and splice sites
*/

// Override baseDir from nextflow repository to current directory
baseDir = "$PWD"
date = java.time.LocalDate.now()

// Default parameters -- Adapt to experiment!!!
params.QCdir = "results/RNAseq/QC_sequencing/"
params.mapDir = "results/RNAseq/alignment/"
params.quantDir = "results/RNAseq/quantification/"
params.inFold = baseDir
params.fastqcScreen_refs = "/home/SHARED/DATA/REFERENCES/FastQ_Screen_Genomes"
params.fastqcScreen_config = "$baseDir/scripts/fastq_screen.conf"
params.version = null
params.fastaRef = "/home/SHARED/DATA/REFERENCES/GRCh37/Sequence/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz" //Compressed file
params.gtfRef = "/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz" //Compressed file
params.cores = 1

/*
* Convert parameters with strings to files
*/
fastqcScreen_refs = file("${params.fastqcScreen_refs}")
fastqcScreen_config = file("${params.fastqcScreen_config}")
fastaRef = file("${params.fastaRef}")
gtfRef = file("${params.gtfRef}")
cpus = params.cores

// Store workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

// Create channels with paired fastq files
fastqs = Channel.fromFilePairs("${params.inFold}/*{1,2}.fastq.gz", flat: true)
fastqs.into { fastq_QC; fastq_QCscreen; fastq_STAR }

// Run FastQC on reads
process runFastQC {

  input:
  set pair_id, file(read1), file(read2) from fastq_QC

  output:
  file("*_fastqc.{zip,html}") into fastqc_results

  """
  fastqc ${read1} ${read2}
  """

}

// Run FastQC_screen on reads
process runFastQCscreen {

  input:
  file fastqcS_dir from fastqcScreen_refs
  file config from fastqcScreen_config
  set pair_id, file(read1), file(read2) from fastq_QCscreen

  output:
  file("${pair_id}_1_screen.txt") into fastqcS1
  file("${pair_id}_2_screen.txt") into fastqcS2

  """
  fastq_screen --conf ${config} ${read1} --aligner bowtie2
  fastq_screen --conf ${config} ${read2} --aligner bowtie2
  """

}

// Index genome
process indexGenome {

  input:
  file(fasta) from fastaRef
  file(gtf) from gtfRef

  output:
  file("indexes/") into indexes
  file("genes.gtf") into gtfFile

  """
  gzip -cd $fasta > genome.fa
  gzip -cd $gtf > genes.gtf
  sed 's/^chr//' genes.gtf > genesChr.gtf

  mkdir indexes
  STAR --runMode genomeGenerate --genomeDir indexes \
            --genomeFastaFiles genome.fa \
            --sjdbGTFfile genesChr.gtf \
            --runThreadN $cpus
  """
}


// Align reads with STAR
process AlignReads {

  label 'RNAseq_align'

  publishDir "${params.mapDir}/$date", pattern: '*.{bam,txt}', mode: 'copy'
  publishDir "${params.quantDir}/$date", pattern: '*.{tab,txt}', mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.mapDir}/${params.version}", pattern: '*.{bam,txt}', mode: 'copy'
    publishDir "${params.quantDir}/${params.version}", pattern: '*.{tab,txt}', mode: 'copy'
  }

  input:
  file(indexes) from indexes
  set pair_id, file(read1), file(read2) from fastq_STAR
  val logText from "$workflowInfo"

  output:
  file ("${pair_id}SJ.out.tab") into junction
  file ("${pair_id}ReadsPerGene.out.tab") into geneQuanti
  file ("${pair_id}Aligned.out.bam") into align
  file ("${pair_id}Log.final.out") into logSTAR
  file 'log.txt'

  """
  STAR --genomeDir $indexes \
      --readFilesIn ${read1} ${read2} \
      --readFilesCommand zcat \
      --outSAMtype BAM Unsorted \
      --quantMode GeneCounts \
      --outFileNamePrefix ${pair_id}  \
      --runThreadN $cpus
  echo "$logText" > log.txt
  """
}

process multiqc {

  publishDir "${params.QCdir}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.QCdir}/${params.version}", mode: 'copy'
  }

  input:
  file ('FastQC_screen/*') from fastqcS1.collect().ifEmpty([])
  file ('FastQC_screen/*') from fastqcS2.collect().ifEmpty([])
  file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
  file ('align/*') from geneQuanti.collect().ifEmpty([])
  file ('align/*') from logSTAR.collect().ifEmpty([])
  val logText from "$workflowInfo"

  output:
  file "multiqc_report.html" into multiqc_report
  file "multiqc_data"
  file 'log.txt'

  script:

  """
  multiqc .
  echo "$logText" > log.txt
  """
}
