process fastqc {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out + "/fastqc", mode:'copy'
  
  input: 
  tuple val(sample), path(reads)
  
  output:
  path("${sample}/*"), emit: fq_files

  script:
  """
    mkdir -p ${sample}
    fastqc -o ${sample} -f fastq -q ${reads[0]} ${reads[1]}
  """   
}

