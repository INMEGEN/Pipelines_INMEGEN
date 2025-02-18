process featureCounts {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir}:/ref"
  publishDir params.out + "/featureCounts" , mode: 'copy'

  input:
  tuple val(sample), path(bam)

  output:
  tuple val(sample), path("${sample}_fcounts.tsv"), emit: fcounts
  tuple val(sample), path("quants/${sample}_counts.tsv"), emit: counts
  path("*.jcounts")
  path("*.summary")

  script:
  """
  featureCounts -a /ref/${params.gtfname} -T ${params.ncrs} -B -C -g gene_id -J -p --countReadPairs ${params.ss} -t exon -o ${sample}_fcounts.tsv ${bam}

  awk 'NR == 2 {gsub(/_Aligned\\.sortedByCoord\\.out\\.bam/, "", \$7); print \$1,\$6,\$7; next} NR > 2 {print \$1,\$6,\$7}' OFS='\t' ${sample}_fcounts.tsv > ${sample}_counts.tsv

  mkdir quants
  mv ${sample}_counts.tsv quants/
  """  
}
