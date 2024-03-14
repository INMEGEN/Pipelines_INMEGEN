process mothur_QC2 {
    cache 'lenient'
    conda '/scratch/home/dperez/programas/miniconda3/envs/nf-meta'
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    //containerOptions "-v ${params.ref_silva}:/ref_silva"
    publishDir params.out + "/mothur_QC2", mode:'copy'

    input:
    tuple val(project_id), path(fasta), path(count_table)

    output:
    tuple val(project_id), path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta"), \
    path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table"), emit: mth_qc2
    path("${project_id}.*")

    script:
    """
   ${params.mothur}/mothur "#classify.seqs(fasta=${fasta}, count=${count_table}, reference=${params.ref_silva}/silva.nr_v132.align, \
   taxonomy=${params.ref_silva}/silva.nr_v132.tax, cutoff=80, processors=${params.ncrs})"

   ${params.mothur}/mothur "#remove.lineage(fasta=${fasta}, count=${count_table}, \
   taxonomy=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.nr_v132.wang.taxonomy, \
   taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

   ${params.mothur}/mothur "#summary.seqs(fasta=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.fasta, \
   count=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, processors=${params.ncrs})"
    """
}
