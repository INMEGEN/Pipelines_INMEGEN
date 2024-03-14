process trim_Galore {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out + "/trimming_files", mode: 'copy'

  input:
  tuple val(sample), path(R1), path(R2)

  output:
  tuple val(sample), path("data/${sample}_R1_trimmed.fq.gz") ,path("data/${sample}_R2_trimmed.fq.gz")  , emit: trim_fq
  tuple val(sample), path("*report.txt")                                                 , emit: trim_report
  tuple val(sample), path("*_fastqc.*")                                                  , emit: trim_fastqc

  script:
  """
  trim_galore -j ${params.ncrs} --basename ${sample} --gzip --fastqc --paired ${R1} ${R2}
 
  mkdir data
 
  mv ${sample}_val_1.fq.gz data/${sample}_R1_trimmed.fq.gz
  mv ${sample}_val_2.fq.gz data/${sample}_R2_trimmed.fq.gz
  """
}
