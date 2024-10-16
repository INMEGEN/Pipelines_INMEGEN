process fastp {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out + "/trimmed_files", mode: 'symlink'

  input:
  tuple val(sample), path(R1), path(R2)

  output:
  tuple val(sample), path("${sample}_R1.trimmed.fq.gz"), path("${sample}_R2.trimmed.fq.gz"), emit: trim_fq
  tuple path("${sample}_fastp.html"), path("${sample}_fastp.json")

  script:
  """
  fastp --in1 ${R1} --in2 ${R2} --out1 ${sample}_R1.trimmed.fq.gz --out2 ${sample}_R2.trimmed.fq.gz -g 10 -q 20 -l 50 -h ${sample}_fastp.html

  mv fastp.json ${sample}_fastp.json
  """
}
