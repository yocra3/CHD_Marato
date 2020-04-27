/*
* Run RNAseq quantification in genes, transcripts and splice sites
*/

// Override baseDir from nextflow repository to current directory
baseDir = "$PWD"
date = java.time.LocalDate.now()

// Default parameters -- Adapt to experiment!!!
params.bams = ""
params.genome = ""
params.binSize = 100
params.sampleAnnot = ""
params.outFolder = "results/CNVs/"
params.version = null
params.fastaRef = "/home/SHARED/DATA/REFERENCES/GRCh37/Sequence/Homo_sapiens.GRCh37.dna_rm.primary_assembly.fa.gz" //Compressed file
params.gapsRef = "/home/SHARED/DATA/REFERENCES/Archive/RLCR/RLCRs.txt"

// Select containers
container_ubuntu = 'yocra3/ubuntu_genomicutils:release-0.99.5'

/*
* Convert parameters with strings to files
*/
fastaRef = file("${params.fastaRef}")
gapsRef = file("${params.gapsRef}")
sampleAnnot = file(params.sampleAnnot)

bams = Channel.fromPath("${params.bams}")
              .map { file -> tuple(file.simpleName, file) }

// Store workflow info
workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

// Extract read mapping
process extractReadMapping {

  container container_ubuntu

  label 'process_medium'

  input:
  set sampID, file(bam) from bams
  val(genome) from params.genome

  output:
  set sampID, file("${sampID}.root") into root

  """
  cnvnator -genome $genome -root ${sampID}.root -tree $bam
  """

}

// Split fasta in chromosomes
process splitFasta {

  container container_ubuntu

  input:
  file(fasta) from fastaRef

  output:
  file("*.fa") into fastaChr

  """
  gzip -cd $fasta > genome.fasta

  ## Split fasta in chromosomes
  csplit -s -z genome.fasta '/>/' '{*}'
  for i in xx*
  do
    n=\$(sed 's/>// ; s/ .*// ; 1q' "\$i") ; \
    mv "\$i" "\$n.fa" ; \
  done
  """
}

// Generate histograms
process generateHistogram{

  container container_ubuntu

  input:
  set sampID, file("file.root") from root
  val(genome) from params.genome
  val(bin_size) from params.binSize
  file(fasta) from fastaChr

  output:
  set sampID, file("${sampID}.root") into root_hist

  """
  cp file.root "${sampID}.root"
  cnvnator -genome $genome -root "${sampID}.root" -his $bin_size
  """
}


// Compute statistics
process computeStats {

  container container_ubuntu

  label 'process_low'

  input:
  set sampID, file("file.root") from root_hist
  val(bin_size) from params.binSize

  output:
  set sampID, file("${sampID}.root") into root_stats

  """
  cp file.root "${sampID}.root"
  cnvnator -root "${sampID}.root" -stat $bin_size
  """
}

// Partition signal
process partSignal {

  container container_ubuntu

  label 'process_low'

  input:
  set sampID, file("file.root") from root_stats
  val(bin_size) from params.binSize

  output:
  set sampID, file("${sampID}.root") into root_part

  """
  cp file.root "${sampID}.root"
  cnvnator -root "${sampID}.root" -partition $bin_size -ngc
  """
}


// Call CNVs
process callCNVs {

  container container_ubuntu

  input:
  set sampID, file("file.root") from root_part
  val(bin_size) from params.binSize

  output:
  set sampID, file("${sampID}.txt") into calls

  """
  cp file.root "${sampID}.root"
  cnvnator -root "${sampID}.root" -call $bin_size -ngc > ${sampID}.txt
  """
}

// Format results
process formatResults {

  container container_ubuntu

  input:
  set sampID, file(call) from calls

  output:
  set sampID, file("${sampID}.filtered.txt") into filterCalls

  """
  format_cnvnator_results.py $call $sampID
  """
}

// Merge CNVs and remove CNVs in gap locations
process mergeCNVs {

  container container_ubuntu

  input:
  set sampID, file(call) from filterCalls
  file(gaps) from gapsRef

  output:
  set sampID, file("merged/${call}.cluster.txt") into mergedCNVs

  """
  echo -e $call"\t"$sampID >> ALT_CNVN.txt

  mkdir calls
  mkdir merged
  mv $call calls

  merge_cnvnator_results.py -i ./calls/ -a ALT_CNVN.txt -o merged -g $gaps
  """
}

// Convert CNVnator output to VCF
process convert2VCF {

  container container_ubuntu


  input:
  set sampID, file(merged) from mergedCNVs
  file (sampleTab) from sampleAnnot

  output:
  file "*.vcf" into vcfs

  """
  id="$sampID"
  samp=`(grep \${id/_S[0-9]*/} ${sampleTab} | cut -f30 | sed 's/-//g')`

  echo "##fileformat=VCFv4.3
##fileDate=$date
##source=CNVnator
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=NRD,Number=1,Type=String,Description=\"Normalized Read depth\">
##FORMAT=<ID=E,Number=1,Type=Float,Description=\"e-val2\">
##FORMAT=<ID=q0,Number=1,Type=Float,Description=\"q0\">
##FORMAT=<ID=NCNV,Number=1,Type=Integer,Description=\"Number of CNVs\">
##FORMAT=<ID=LCNVS,Number=1,Type=Integer,Description=\"Length of CNVS\">
##FORMAT=<ID=LG,Number=1,Type=Integer,Description=\"Length of gaps\">
##FORMAT=<ID=PG,Number=1,Type=Float,Description=\"Proportion of Gaps\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\${samp}" > \${samp}.vcf
  awk  '{OFS = "\t"}!/^#/{print \$2, \$3, ".", "N", ".", "<"\$5">", "PASS", "SVTYPE="\$5";END="\$4";SVLEN="\$6,\
   "GT:NRD:E:q0:NCNV:LCNVS:LG:PG", "1/0:"\$7":"\$8":"\$9":"\$10":"\$11":"\$12":"\$13  }' \
   $merged >> \${samp}.vcf
   sed -i 's/:-/:\\./g' \${samp}.vcf ## Change empty values to .

  """
}

process createFastaIndex {
  container container_ubuntu

  input:
  file(fasta) from fastaRef

  output:
  file("hg19.fasta.fai") into fastaIDX

  """
  gzip -cd $fasta > hg19.fasta
  samtools faidx hg19.fasta
  """

}

// Compress, sort and index vcf
process prepareVCF {

  container container_ubuntu

  publishDir "${params.outFolder}/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.outFolder}/${params.version}", mode: 'copy'
  }

  input:
  file(vcf) from vcfs
  val logText from "$workflowInfo"
  file(fastaidx) from fastaIDX

  output:
  file("${vcf}.gz") into sortVCF
  file("${vcf}.gz.tbi") into tbiVCF
  file 'log.txt'

  """


  bcftools reheader -f $fastaidx -o header.vcf $vcf
  bcftools sort -o ${vcf}.gz -O z  header.vcf
  tabix -p vcf ${vcf}.gz
  echo "$logText" > log.txt
  """
}
