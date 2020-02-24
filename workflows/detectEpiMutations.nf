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
  file("all_bumps.txt") into allbumps
  file("sel_bumps.txt") into selbumps

  """
  touch all_bumps.txt
  touch sel_bumps.txt
  Rscript $baseDir/scripts/detectEpimutations.R $gset $sample
  """

}

process mergeBumps {

  publishDir "${params.outdir}/${date}", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.outdir}/${params.version}", mode: 'copy'
  }

  input:
  file("all_bumps*.txt") from allbumps.toList()
  file("sel_bumps*.txt") from selbumps.toList()
  val logText from "$workflowInfo"

  output:
  file("bumpsAll.txt") into bumpsall
  file("bumpsSel.txt") into bumpsfilt
  file 'log.txt' into logCh

  """
  cat all_bumps*.txt >> bumpsAll.txt
  cat sel_bumps*.txt >> bumpsSel.txt
  echo "$logText" > log.txt
  """
}
