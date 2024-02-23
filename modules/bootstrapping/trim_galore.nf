process trim_Galore {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out + "/trimming_files", mode: 'copy'

  input:
  tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path(R1), path(R2)

  output:
  tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path("${sample_id}_R1_trimmed.fq.gz") ,path("${sample_id}_R2_trimmed.fq.gz")  , emit: trim_fq
  tuple val(sample), path("*report.txt")                                                 , emit: trim_report
  tuple val(sample), path("*_fastqc.*")                                                  , emit: trim_fastqc

  script:
  """
  trim_galore -j ${params.ncrs} --basename ${sample_id} --gzip --fastqc --paired ${R1} ${R2}

  mv ${sample_id}_val_1.fq.gz ${sample_id}_R1_trimmed.fq.gz
  mv ${sample_id}_val_2.fq.gz ${sample_id}_R2_trimmed.fq.gz
  """
}
