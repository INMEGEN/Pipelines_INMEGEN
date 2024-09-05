process trimmomatic {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out + "/trimmomatic", mode: 'symlink'

  input:
  tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path(R1), path(R2)
  file(adapters)

  output:
  tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: trim_fq
  path("*un.trimmed.fastq.gz")

  script:
  """
  trimmomatic PE -threads ${params.ncrs} ${R1} ${R2} \
  ${sample_id}_R1.trimmed.fastq.gz ${sample_id}_R1un.trimmed.fastq.gz \
  ${sample_id}_R2.trimmed.fastq.gz ${sample_id}_R2un.trimmed.fastq.gz \
  ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30;
  """
}
