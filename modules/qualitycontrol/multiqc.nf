process multiqc {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out, mode: 'copy'

  input:
  val(sample_1)
  path(dir_all)
  val(name)

  output:
  path("multiqc/*")   , emit: multiqc_fq_data

  script:
  """
    mkdir -p multiqc/${name}

    multiqc -o multiqc/${name} ${dir_all}
  """
}

