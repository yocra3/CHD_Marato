/*
 * Performs cleaning of phenotype data
 */

baseDir = "$PWD"

workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $nextflow.timestamp
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""


process processPhenotypes {
   	publishDir $baseDir, mode: 'copy'  
    
  
   	input:
   	val logText from "$workflowInfo"

    output:
    file 'log.txt' into logCh

    """
    echo "$logText" > log.txt
    """
}

