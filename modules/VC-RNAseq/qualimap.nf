process qualimap {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir_star}:/ref"  
  publishDir params.out +"/qualimap", mode: 'copy'

  input:
  tuple val(sample), path(bam), path(index)

  output:
  path("${sample}/*"), emit: qcmap

  script:
  """
  mkdir ${sample}
  qualimap rnaseq -a proportional -bam ${bam} -gtf /ref/${params.gtfname} -outdir ${sample}/ -outfile ${sample}_report.pdf  -p strand-specific-reverse -pe -s --java-mem-size=12G
  """
}
