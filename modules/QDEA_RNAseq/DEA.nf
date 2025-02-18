process DEA {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir}:/ref"
  publishDir params.out + "/DEA", mode: 'copy'

  input:
  path(mcounts)
  file(metadata)
  file(script)

  output:
  path("*_mqc.png")                                      , emit: clustering_plot
  path("${params.DEAname}/${params.volcano_plot_name}")  , emit: volcano_plot
  path("${params.DEAname}/${params.results_name}")       , emit: results
  path("${params.DEAname}/${params.deg_name}")           , emit: results_f
  path("${params.DEAname}/heatmap_*.png")                , emit: plots
  path("${params.DEAname}/DESeq2_dds.rds")               , emit: rds
  path("${params.DEAname}/*.log")                        , emit: R_sesion_info

  script:
  """  
  mkdir -p ${params.DEAname}
  
  Rscript ${script} \
           --meta_data ${metadata} \
           --gtf_file /ref/${params.gtfname} \
           --condition1 ${params.condition_1} \
           --condition2 ${params.condition_2} \
           --Log2FC_th ${params.th_l2fc} \
           --p_adj_th ${params.th_padj} \
           --outdir_pca ${params.pca_plot_name} \
           --out_p_hm ${params.heatmap_name} \
           --outdir_vp ${params.volcano_plot_name} \
           --out_res ${params.results_name} \
           --out_deg ${params.deg_name} \
           --nsamples ${params.nsamples}
  
  mv ${params.volcano_plot_name} ${params.results_name} ${params.deg_name} heatmap_zscore.png heatmap_log2.png DESeq2_dds.rds ${params.DEAname}/
  mv *.log ${params.DEAname}/
  """
}
