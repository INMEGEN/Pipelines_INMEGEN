process fastqc {
  container 'pipelinesinmegen/pipelines_inmegen:latest'
  cache 'lenient'
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

