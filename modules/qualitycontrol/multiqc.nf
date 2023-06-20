process multiqc {
  container 'pipelinesinmegen/pipelines_inmegen:latest'
  cache 'lenient'
  publishDir params.out, mode: 'copy'

  input:
  val(sample_1)
  path(dir_all)

  output:
  path("multiqc/*")   , emit: multiqc_fq_data

  script:
  """
    mkdir -p multiqc

    multiqc -o multiqc/ ${dir_all}
  """
}

