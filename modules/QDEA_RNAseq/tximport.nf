process tximport {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir}:/ref"
  publishDir params.out + "/smcounts", mode: 'copy'

  input:
  val(sample_k)
  file(sample_info)
  path(salmon_dir)
  file(script)

  output:
  path("${params.mcounts}")         , emit: mcounts
  path("${params.mcounts_tpm}")     , emit: mcounts_tpm
  path("*.tsv") 
  path("*.log")                     

  script:
  """
  mkdir -p /wdir/salmon  
  cp -r ${salmon_dir}/* /wdir/salmon

   Rscript ${script} \
   --working_dir /wdir \
   --sample_info ${sample_info} \
   --dir_quants "salmon" \
   --gtf_file /ref/${params.gtfname} \
   --countmat ${params.mcounts} \
   --countpm ${params.mcounts_tpm}

  rm -r ${salmon_dir}
  """
}
