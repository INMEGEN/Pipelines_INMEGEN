process tximport_deseq2 {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir}:/ref"
  publishDir params.out + "/DEA_" + params.DEAname , mode: 'copy'

  input:
  val(sample_k)
  file(sample_info)
  path(klx_dir)
  file(script)

  output:
  path("${params.pca_plot_name}")      , emit: pca_plot
  path("${params.heatmap_name}")       , emit: heatmap_f
  path("${params.volcano_plot_name}")  , emit: volcano_plot
  path("${params.results_name}")       , emit: results
  path("${params.deg_name}")           , emit: results_f
  path("${params.mcounts}")            , emit: mcounts
  path("${params.mcounts_tpm}")        , emit: mcounts_tpm
  path("enrichr_*")                    , emit: erichR
  path("*.rds")                        , emit: rds
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
           --cond_colum ${params.cond_colum} \
           --condition1 ${params.condition_1} \
           --condition2 ${params.condition_2} \
           --Log2FC_th ${params.th_l2fc} \
           --p_adj_th ${params.th_padj} \
           --outdir_pca ${params.pca_plot_name} \
           --out_p_hm ${params.heatmap_name} \
           --outdir_vp ${params.volcano_plot_name} \
           --out_res ${params.results_name} \
           --out_deg ${params.deg_name} \
           --countmat ${params.mcounts} \
           --countpm ${params.mcounts_tpm} \
  """
}

