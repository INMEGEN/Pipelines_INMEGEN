process trim_Galore {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out + "/trimming_files", mode: 'copy'

  input:
  tuple val(sample_id), val(sample), val(RG), val(PU), path(read_1), path(read_2)

  output:
  tuple val(sample_id), val(sample), val(RG), val(PU), path("*.fq.gz")  , emit: trim_fq
  tuple val(sample), path("*report.txt")                                , emit: trim_report
  tuple val(sample), path("*_fastqc.*")                                 , emit: trim_fastqc

  script:
  """
  trim_galore -j ${params.ncrs} --basename ${sample_id} --fastqc --paired ${read_1} ${read_2}

  mv ${sample_id}_val_1.fq.gz ${sample_id}_R1_trimmed.fq.gz
  mv ${sample_id}_val_2.fq.gz ${sample_id}_R2_trimmed.fq.gz
  """
}
