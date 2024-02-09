process multiqc {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out, mode: 'copy'

  input:
  val(sample_1)
  path(dir_all)

  output:
  path("multiqc/*")   , emit: multiqc_fq_data

  script:
  """
    multiqc -o multiqc/ ${dir_all}
  """
}
