/*
* Run RNAseq quantification in genes, transcripts and splice sites
*/

// Override baseDir from nextflow repository to current directory
baseDir = "$PWD"
date = java.time.LocalDate.now()

// Default parameters -- Adapt to experiment!!!
params.resDir = "results/RNAseq/OUTRIDER_aberrations/"
params.rangedSE = ""
params.gtfRef = "/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz" //Compressed file
params.version = null

/*
* Convert parameters with strings to files
*/
rseFile = file("${params.rangedSE}")
gtfRef = file("${params.gtfRef}")

// Store workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

// Run epimutation algorithm to each sample
process runOUTRIDER {

  publishDir "${params.resDir}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.resDir}/${params.version}", mode: 'copy'
  }

  input:
  file(rseFile) from rseFile
  file(gtfRef) from gtfRef
  val logText from "$workflowInfo"

  output:
  file("OUTRIDER_QC.Rdata")
  file("OUTRIDER_res.Rdata")
  file 'log.txt'

  """
  Rscript $baseDir/scripts/detectTranscriptomeAberration.R $rseFile $gtfRef
  echo "$logText" > log.txt
  """

}
