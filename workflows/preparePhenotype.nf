/*
 * Performs cleaning of phenotype data
 */

baseDir = "$PWD"
date = java.time.LocalDate.now()

params.outdir = "results/phenotypes/${date}"
params.inFold = baseDir

workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""


process processPhenotypes {
   	publishDir params.outdir, mode: 'copy'

   	if ( params.version != null ){
   		publishDir "results/phenotypes/${params.version}", mode: 'copy'
   	}

   	input:
   	path inFold from "${params.inFold}"
   	val logText from "$workflowInfo"

    output:
    file 'pheno.tab' into pheno_text
    file 'pheno.Rdata' into pheno_R
    file 'codebook.tab' into codebook
    file 'log.txt' into logCh

    """
    Rscript $baseDir/scripts/preparePhenotypesv2.R $inFold
    echo "$logText" > log.txt
    """
}
