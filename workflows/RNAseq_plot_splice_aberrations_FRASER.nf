/*
* Make plots of splicing aberrations detected with FRASER
*/

// Override baseDir from nextflow repository to current directory
baseDir = "$PWD"
date = java.time.LocalDate.now()

// Default parameters -- Adapt to experiment!!!
params.plotsDir = "results/RNAseq/FRASER_splicing_aberrations_plots/"
params.bams = ""
params.fraserFold = ""
params.gtfRef = "/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz" //Compressed file
params.version = null

/*
* Convert parameters with strings to files
*/
gtfRef = file("${params.gtfRef}")

/*
* These lines assumes that FRASER results are stored in FRASER_results.Rdata
* and that FRASER objects points to current working directory.
*/
fraserRes = file("${params.fraserFold}FRASER_results.Rdata")
fraserFold = file("${params.fraserFold}/savedObjects/")

// Create bams channnel
bams = Channel.fromPath("${params.bams}/*.bam")
bais = Channel.fromPath("${params.bams}/*.bam.bai")

// Select containers
container_R = 'yocra3/rsession_chd_marato:release-1.2.5'
container_ubuntu = 'yocra3/ubuntu_genomicutils:release-0.99.5'
container_ggSashimi = 'guigolab/ggsashimi:0.5.0'

// Store workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

process prepareTargetFiles {

  container container_R

  input:
  file("res.Rdata") from fraserRes
  file(fold) from fraserFold

  output:
  file("*.txt") into sampFiles
  file("ggsashimi_vars.csv") into ggsashimilines

  """
  prepareSashimiFiles.R res.Rdata
  """

}

process decompressGTF {

  container container_ubuntu

  input:
  file(gtfRef) from gtfRef

  output:
  file("genesChr.gtf") into geneRef

  """
  gzip -cd $gtfRef | awk -F "\t" '\$3=="exon"||\$3=="transcript"' > genes.gtf
  sed 's/^chr//' genes.gtf > genesChr.gtf
  """

}

ggsashimilines
    .splitCsv(header:true)
    .map{ row -> tuple(row.range, row.bamFile, row.rangeName) }
    .set{ sashimi }

process makeSashimiPlots  {

  container container_ggSashimi
  containerOptions = "--entrypoint bash"

  publishDir "${params.plotsDir}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.plotsDir}/${params.version}", mode: 'copy'
  }

  input:
  set range, bamLink, plotName from sashimi
  file("*") from bams.toList()
  file("*") from bais.toList()
  file("*") from sampFiles
  file(gtfRef) from geneRef
  val logText from "$workflowInfo"

  output:
  file("${plotName}.png") into plots
  file("log.txt")

  """
  /sashimi-plot.py -b $bamLink -c $range -M 10 -C 3 -O 3 -A median --alpha 1 -g $gtfRef \
  -F png -o $plotName
  echo "$logText" > log.txt
  """


}
