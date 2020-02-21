/*
 * Extract SNPs present in Illumina 450K from WGS data
 */

// Override baseDir from nextflow repository to current directory
baseDir = "$PWD"
date = java.time.LocalDate.now()

params.outdir = "results/methylation/Epimutations/"
params.sampleAnnot = baseDir
params.gset = baseDir

sampleAnnot = file(params.sampleAnnot)
gset = file(params.gset)


workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

process getCases {

  input:
  file(pheno) from sampleAnnot

  output:
  file("cases.txt") into casesFile

  """
  grep Case $pheno | cut -f1 > cases.txt
  """

}


cases = casesFile.splitText()

process runEpimutations {

  input:
  val sample from cases
  file(gset) from gset

  output:
  file("bumps.Rdata") into bumps

  """
  Rscript $baseDir/scripts/detectEpimutations.R $gset $sample
  """

}

process mergeBumps {

  publishDir "${params.outdir}/${date}", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.outdir}/${params.version}", mode: 'copy'
  }

  input:
  file("bumps*.Rdata") from bumps.toList()
  val logText from "$workflowInfo"

  output:
  file("bumps.combined.Rdata") into bumpsR
  file("bumps.combined.txt") into bumpsText
  file 'log.txt' into logCh

  """
  Rscript $baseDir/scripts/mergeEpimutations.R bumps*.Rdata
  echo "$logText" > log.txt
  """
}
