process salmonML {
  cache 'lenient'
  container 'combinelab/salmon:latest'
  containerOptions "-v ${params.refdir}:/ref"
  publishDir params.out +"/salmon", mode: 'copy'

  input:
  tuple val(sample), path(R1), path(R2)
  file(gene_ids)

  output:
  tuple val(sample), path("${sample}/quant.sf")            , emit: quant  
  tuple val(sample), path("${sample}/quant.genes.sf")      , emit: quant_genes
  path("${sample}/*")

  script:
  """
  salmon quant -p ${params.ncrs} -i /ref -l ${params.slib_type} -g ${gene_ids}  -1 ${R1[0]} ${R1[1]}  -2 ${R2[0]} ${R2[1]} --validateMappings --gcBias --seqBias -o ${sample}
  """
}
