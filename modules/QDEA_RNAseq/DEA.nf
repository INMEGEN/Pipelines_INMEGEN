process tximport_deseq2 {
  container 'pipelinesinmegen/pipelines_inmegen:latest'
  containerOptions "-v ${params.refdir}:/ref"
  cache 'lenient'
  publishDir params.out + "/resultados", mode: 'copy'

  input:
  val(sample_k)
  file(sample_info)
  path(klx_dir)
  file(script)

  output:
  path("${params.pca_plot_name}")      , emit: pca_plot
  path("${params.heatmap_name}")       , emit: heatmap_f
  path("${params.volcano_plot_name}")  , emit: volcano_plot
  path("${params.results_c_file}")     , emit: results_c
  path("${params.results_file_name}")  , emit: results
  path("${params.resultsf_file_name}") , emit: results_f
  path("${params.mcounts}")            , emit: mcounts
  path("${params.mcounts_tpm}")        , emit: mcounts_tpm
  //path("${params.list_gt}")          , emit: gtl_n
  path("*.log")                        , emit: R_sesion_info

  script:
  """ 
  mkdir -p /wdir/kallisto_quants  
  cp -r ${klx_dir}/* /wdir/kallisto_quants 

  Rscript ${script} \
           --working_dir /wdir \
           --sample_info ${sample_info} \
           --dir_quants "kallisto_quants" \
           --gtf_file /ref/${params.gtfname} \
           --condition1 ${params.condition_1} \
           --condition2 ${params.condition_2} \
           --Log2FC_th ${params.th_l2fc} \
           --p_adj_th ${params.th_padj} \
           --outdir_pca ${params.pca_plot_name} \
           --out_p_hm ${params.heatmap_name} \
           --outdir_vp ${params.volcano_plot_name} \
           --outdir_cvs ${params.results_c_file} \
           --outres_cvs ${params.results_file_name} \
           --outfilt_cvs ${params.resultsf_file_name} \
           --countsmat ${params.mcounts} \
           --countpm ${params.mcounts_tpm} \
           --data_b ${params.data_set} \
           --onto ${params.onto_type} \
           --gtl_o ${params.list_gt}

   rm -r ${klx_dir}
  """
}
