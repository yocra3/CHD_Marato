/*
 * Create Bioc Objects and run Quality control of RNAseq data
 */

baseDir = "$PWD"
date = java.time.LocalDate.now()
params.version = null

params.qcdir = "results/RNAseq/QC_files/"
params.outdir = "results/RNAseq/Bioc_objects/"
params.countsPath = ""
params.phenoPath = ""
params.gtfRef = "/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz"
//params.gtfRef should map to the same file used in mapping

/*
 * Convert parameters with strings to files
 */
countsPath = file(params.countsPath)
phenoPath = file(params.phenoPath)
gtfRef = file(params.gtfRef)

// Store workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

process createRangedSummarizedExperiment {

  publishDir "${params.outdir}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.outdir}/${params.version}",  mode: 'copy'
  }

  input:
  file (counts) from countsPath
  file (pheno) from phenoPath
  file (annot) from gtfRef
  val logText from "$workflowInfo"

  output:
  file("RNAseq_RangedSE_allGenes.Rdata") into rseAll
  file("RNAseq_RangedSE_expressedGenes.Rdata") into rseFilt
  file("RNAseq_RangedSE_autosomicGenes.Rdata") into rse

  file 'log.txt'

  """
  Rscript $baseDir/scripts/createRNAseqSummarizedExperiment.R $counts $annot $pheno
  echo "$logText" > log.txt
  """
}

process createQC {

 publishDir "${params.qcdir}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.qcdir}/${params.version}",  mode: 'copy'
  }

  input:
  file (rseF) from rse
  val logText from "$workflowInfo"
  file (counts) from countsPath //Añado esta linea para que encuentre el script

  output:
  file("RNAseq_QCobject.Rdata") into objQC
  file 'log.txt'

  """
  Rscript $baseDir/scripts/createRNAseqQCobjects.R $rseF
  echo "$logText" > log.txt
  """
}
