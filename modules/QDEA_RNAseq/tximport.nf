process tximport_q {
  container 'pipelines_inmegen:latest'
  containerOptions "-v ${params.refdir}:/ref"
  cache 'lenient'
  publishDir params.out + "/resultados_q", mode: 'copy'

  input:
  val(sample_k)
  file(sample_info)
  path(klx_dir)
  file(script)

  output:
  path("${params.mcounts}")         , emit: mcounts
  path("${params.mcounts_tpm}")     , emit: mcounts_tpm
  path("*.log")                     , emit: R_sesion_info

  script:
  """
  mkdir -p /wdir/kallisto_quants  
  cp -r ${klx_dir}/* /wdir/kallisto_quants

   Rscript ${script} \
   --working_dir /wdir \
   --sample_info ${sample_info} \
   --dir_quants "kallisto_quants" \
   --gtf_file /ref/${params.gtfname} \
   --countsmat ${params.mcounts} \
   --countpm ${params.mcounts_tpm}

  rm -r ${klx_dir}
  """
}
