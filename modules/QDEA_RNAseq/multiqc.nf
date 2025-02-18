process multiqc {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out, mode: 'copy'

  input:
  val(dir_1)
  file(config)
  path(dir_all)

  output:
  path("multiqc/*")   , emit: multiqc_fq_data

  script:
  """
  multiqc -c ${config} -o multiqc/ ${dir_all}

  echo "Done"  
  """
}

