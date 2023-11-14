process trim_Galore {
  container 'pipelinesinmegen/pipelines_inmegen:public'
  cache 'lenient'
  publishDir params.out + "/trimming_files", mode: 'copy'

  input:
  tuple val(sample), path(read_1), path(read_2)

  output:
  tuple val(sample), path("${sample}/*_trimmed.fq.gz")      , emit: trim_fq
  tuple val(sample), path("${sample}/*report.txt")          , emit: trim_report
  tuple val(sample), path("${sample}/*_trimmed_fastqc.html"), emit: trim_html
  tuple val(sample), path("${sample}/*_trimmed_fastqc.zip") , emit: trim_zip

  script:
  """
  trim_galore -o ${sample} --fastqc ${read_1} ${read_2}
  """
}
