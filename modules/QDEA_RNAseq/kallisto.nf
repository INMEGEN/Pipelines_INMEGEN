process kallisto {
  container 'pipelinesinmegen/pipelines_inmegen:latest'
  containerOptions "-v ${params.refdir}:/db"
  cache 'lenient'
  publishDir params.out + "/kallisto_quants", mode: 'copy'

  input:
  tuple val(sample), path(reads)
  
  output:
  tuple val(sample), path("${sample}/abundance.h5") , emit: abundance_h5
  tuple val(sample), path("${sample}/abundance.tsv"), emit: abundance_tsv
  tuple val(sample), path("${sample}/run_info.json"), emit: run_info  

  script:
  """
  kallisto quant -i /db/${params.refname} -o ${sample} --rf-stranded ${reads[0]} ${reads[1]}
  """
}

