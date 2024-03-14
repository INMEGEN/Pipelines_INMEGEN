process fastqc {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out + "/fastqc", mode:'copy'
 
  input: 
  tuple val(sample), path(R1), path(R2)  
  
  output:
  path("${sample}/*"), emit: fq_files

  script:
  """
   mkdir -p ${sample}
 
  fastqc -o ${sample} -f fastq -q ${R1} ${R2}
  """   
}

