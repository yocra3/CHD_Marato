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
params.commonCNV = "/home/SHARED/DATA/REFERENCES/GRCh37/CNVs/nstd186.GRCh37.variant_call.vcf.gz"
params.clinvarCNV = "/home/SHARED/DATA/REFERENCES/GRCh37/CNVs/nstd102.GRCh37.variant_call.vcf.gz"
params.segDups = "/home/SHARED/DATA/REFERENCES/GRCh37/Repeats/segDups_hg19.txt.gz"
params.gtfRef = "/home/SHARED/DATA/REFERENCES/GRCh37/GenesAnnotation/gencode.v33lift37.annotation.gtf.gz" //Compressed file
params.omim = "/home/SHARED/DATA/REFERENCES/GRCh38/OMIM/genemap2.txt.gz"

// Select containers
container_ubuntu = 'yocra3/ubuntu_genomicutils:release-0.99.5'
container_R = 'yocra3/rsession_chd_marato:release-1.2.4'

/*
* Convert parameters with strings to files
*/
fastaRef = file("${params.fastaRef}")
gapsRef = file("${params.gapsRef}")
sampleAnnot = file(params.sampleAnnot)
commonCNV = file("${params.commonCNV}")
clinvarCNV = file("${params.clinvarCNV}")
segDups = file("${params.segDups}")
gtfRef = file("${params.gtfRef}")
omim = file("${params.omim}")

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
##INFO=<ID=END,Number=1,Type=Integer,Description='End position of the variant described in this record'>
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description='Difference in length between REF and ALT alleles'>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description='Type of structural variant'>
##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>
##FORMAT=<ID=NRD,Number=1,Type=String,Description='Normalized Read depth'>
##FORMAT=<ID=E,Number=1,Type=Float,Description='e-val2'>
##FORMAT=<ID=q0,Number=1,Type=Float,Description='q0'>
##FORMAT=<ID=NCNV,Number=1,Type=Integer,Description='Number of CNVs'>
##FORMAT=<ID=LCNVS,Number=1,Type=Integer,Description='Length of CNVS'>
##FORMAT=<ID=LG,Number=1,Type=Integer,Description='Length of gaps'>
##FORMAT=<ID=PG,Number=1,Type=Float,Description='Proportion of Gaps'>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\${samp}" > \${samp}.CNVs.vcf
  awk  '{OFS = "\t"}!/^#/{print \$2, \$3, ".", "N", ".", "<"\$5">", "PASS", "SVTYPE="\$5";END="\$4";SVLEN="\$6,\
   "GT:NRD:E:q0:NCNV:LCNVS:LG:PG", "1/0:"\$7":"\$8":"\$9":"\$10":"\$11":"\$12":"\$13  }' \
   $merged >> \${samp}.CNVs.vcf
   sed -i 's/:-/:\\./g' \${samp}.CNVs.vcf ## Change empty values to .
   sed -i 's/\\x27/\\x22/g' \${samp}.CNVs.vcf ## Change ' for "

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


  input:
  file(vcf) from vcfs
  file(fastaidx) from fastaIDX

  output:
  set file("${vcf}.gz"), file("${vcf}.gz.tbi")  into sortVCF

  """
  bcftools reheader -f $fastaidx -o header.vcf $vcf
  bcftools sort -o ${vcf}.gz -O z  header.vcf
  tabix -p vcf ${vcf}.gz
  """
}

// Add annotation to CNVs
process annotateVCF {

  container container_R

  publishDir "${params.outFolder}/VCFs/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.outFolder}/VCFs/${params.version}", mode: 'copy'
  }

  input:
  set file(vcf), file(vcftbi) from sortVCF
  file(commonCNV)
  file(clinvarCNV)
  file(segDups)
  file(gtfRef)
  file(omim)

  output:
  set samp, file("${samp}_annotated.vcf.gz"), file("${samp}_annotated.vcf.gz.tbi")  into annotatedVCF

  script:
  samp = vcf.toString() - '.vcf.gz'
  """
  annotateCNVs.R $vcf $commonCNV $clinvarCNV $segDups $gtfRef $omim ${samp}_annotated.vcf

  bgzip ${samp}_annotated.vcf
  tabix -p vcf ${samp}_annotated.vcf.gz
  """

}

/*
* Prioritize variants:
* - Remove CNVs with overlap > 20% with commonCNVs
* - Remove CNVs with overlap > 50% with segmental duplications
* - Select CNVs with overlap > 80% pathogenic variants (subset1)
* - Select CNVs overlapping exons in OMIM genes (subset2)
* - Select CNVs overlapping exons in GENCODE genes (subset3)
*/
process filterVCF {

  container container_R

  publishDir "${params.outFolder}/CNVpriorizedTables/$date", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.outFolder}/CNVpriorizedTables/${params.version}", mode: 'copy'
  }

  input:
  set samp, file(vcf), file(vcftbi) from annotatedVCF

  output:
  file("${samp}.Prioritization.xlsx") into priorCNV

  script:
  """
  filterCNVs.R $vcf ${samp}.Prioritization.xlsx
  """

}

// Write log files
workflow.onComplete = {
  new File("${params.outFolder}/VCFs/$date/log.txt").text = "$workflowInfo"
  if ( params.version != null ){
    new File("${params.outFolder}/VCFs/${params.version}/log.txt").text = "$workflowInfo"
  }

  new File("${params.outFolder}/CNVpriorizedTables/$date/log.txt").text = "$workflowInfo"
  if ( params.version != null ){
    new File("${params.outFolder}/CNVpriorizedTables/${params.version}/log.txt").text = "$workflowInfo"
  }
}
