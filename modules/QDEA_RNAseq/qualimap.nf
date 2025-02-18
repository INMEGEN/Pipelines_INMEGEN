process qualimap {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir_star}:/ref"
  publishDir params.out +"/qualimap", mode: 'copy'

  input:
  tuple val(sample), path(bam)

  output:
  path("${sample}/*"), emit: qcmap

  script:
  """
  mkdir ${sample}
  qualimap rnaseq -bam ${bam} -gtf /ref/${params.gtfname} -outdir ${sample}/ -outfile ${sample}_report.pdf -p ${params.QMstranded} -s --java-mem-size=12G
  """
}
