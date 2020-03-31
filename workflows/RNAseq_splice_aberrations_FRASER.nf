/*
* Detect RNAseq splicing aberrations with FRASER
* This workflow starts with bam files generated with STAR (two steps)
*/

// Override baseDir from nextflow repository to current directory
baseDir = "$PWD"
date = java.time.LocalDate.now()

// Default parameters -- Adapt to experiment!!!
params.resDir = "results/RNAseq/FRASER_splicing_aberrations/"
params.bamFold = baseDir
params.phenoPath = baseDir
params.version = null
params.cores = 1

/*
* Convert parameters with strings to files
*/
pheno = file("${params.phenoPath}")
cpus = params.cores

// Select containers
container_ubuntu = 'yocra3/ubuntu_genomicutils:release-0.99.4'
container_R = 'yocra3/rsession_chd_marato:release-1.2.4'


// Store workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

bams = Channel.fromPath("${params.bamFold}/G3*.bam")
              .map { file -> tuple(file.simpleName.replaceFirst(/Aligned/, ""), file) }

// Sort bam
process sortBam{

  container container_ubuntu

  input:
  set sampID, file(bam) from bams

  output:
  set sampID, file("sort.bam") into sortbams

  """
  samtools sort -o sort.bam $bam
  """
}

// Index bam
process indexBams {

  container container_ubuntu

  input:
  set sampID, file(bam) from sortbams

  output:
  set sampID, file(bam), file("${bam}.bai") into indexbams

  """
  samtools index -b $bam
  """
}


process filterBams {

  container container_ubuntu

  input:
  set sampID, file(bam), file(bamidx) from indexbams

  output:
  file("${sampID}.bam") into bamsfilt
  file("${sampID}.bam.bai") into bamsidxfilt

  """
  samtools view -b $bam "22" > ${sampID}.bam
  samtools index -b ${sampID}.bam
  """
}


//Create FRASER object
process createFRASER{

  container container_R

  publishDir "${params.resDir}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.resDir}/${params.version}", mode: 'copy'
  }

  input:
  file(bam) from bamsfilt.toList()
  file(bamidx) from bamsidxfilt.toList()
  file(pheno) from pheno
  val(cores) from cpus

  output:
  file("FRASER_obj.Rdata")
  file("FRASER_QC.Rdata")
  file("FRASER_results.Rdata")
  file("nonSplitCounts.tsv.gz") into nonSplit
  file("savedObjects/") into hdf5fold
  file("splitCounts.tsv.gz") into split

  """
  createFRASERobject.R $pheno ./ $cores
  detectSplicingAberrationsFRASER.R FRASER_obj.Rdata
  """

}
