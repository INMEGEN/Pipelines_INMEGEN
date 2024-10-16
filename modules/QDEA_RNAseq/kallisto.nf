process kallisto {
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir}:/ref"
  cache 'lenient'
  publishDir params.out + "/kallisto_quants", mode: 'copy'

  input:
  tuple val(sample), path(R1), path(R2)
  
  output:
  tuple val(sample), path("${sample}/abundance.h5") , emit: abundance_h5
  tuple val(sample), path("${sample}/abundance.tsv"), emit: abundance_tsv
  tuple path("${sample}/${sample}_kallisto.log"), path("${sample}/run_info.json"), emit: run_info  

  script:
  """
  kallisto quant -i /ref/${params.refname} -o ${sample} ${R1} ${R2} 2> ${sample}_kallisto.log

  mv ${sample}_kallisto.log ${sample}/
  """
}

