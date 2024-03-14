process mothur_classify_OTUS {
    cache 'lenient'
    conda params.conda_env
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    //containerOptions "-v ${params.ref_green}:/ref"
    publishDir params.out + "/mothur_OTUS", mode:'copy'

    input:
    tuple val(project_id), path(fasta_file), path(count_table)

    output:
    tuple val(project_id), path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.list"), \
    path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.cons.taxonomy")                           , emit: mth_otu
    path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.shared")                                       , emit: mth_otu_shared
    tuple val(project_id), path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.cons.tax.summary") , emit: mth_tax_summary    
    path("${project_id}*")

    script:
    """
   ${params.mothur}/mothur "#dist.seqs(fasta=${fasta_file})"

   ${params.mothur}/mothur "#cluster(column=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.dist, \
   count=${count_table}, method=opti ,cutoff=0.03 )"

   ${params.mothur}/mothur "#make.shared(list=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.list, \
   count=${count_table}, label=0.03)"

   ${params.mothur}/mothur "#classify.seqs(fasta=${fasta_file}, count=${count_table}, \
   template=${params.ref_green}/gg_13_8_99.fasta, taxonomy=${params.ref_green}/gg_13_8_99.gg.tax, cutoff=80, processors=${params.ncrs})"

   ${params.mothur}/mothur "#classify.otu(list=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.list,\
   taxonomy=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.gg.wang.taxonomy, \
   count=${count_table}, label=0.03, cutoff=80, basis=otu, probs=F)"

   ${params.mothur}/mothur "#rarefaction.single(shared=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.shared)"

   ${params.mothur}/mothur  "#summary.single(shared=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.shared, \
   calc=nseqs-coverage-sobs-chao-shannon-invsimpson, subsample=T)"

   ${params.mothur}/mothur "#make.biom(shared=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.shared, label=0.03, \
   reftaxonomy=${params.ref_green}/gg_13_8_99.gg.tax, constaxonomy=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.cons.taxonomy)" 
    """
}
