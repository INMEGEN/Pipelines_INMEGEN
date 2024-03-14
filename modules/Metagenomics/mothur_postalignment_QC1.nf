process mothur_QC1 {
    cache 'lenient'
    conda '/scratch/home/dperez/programas/miniconda3/envs/nf-meta'
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    //containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/mothur_QC1", mode:'copy'

    input:
    tuple val(project_id), path(align)
    tuple val(project_id), path(fasta), path(count_table)

    output:
    tuple val(project_id), path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.fasta"), \
    path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.count_table"), emit: mth_qc1
    path("${project_id}.*")

    script:
    """
   ${params.mothur}/mothur "#summary.seqs(fasta=${align}, count=${count_table}, processors=${params.ncrs})"

   ${params.mothur}/mothur "#screen.seqs(fasta=${align}, count=${count_table}, summary=${project_id}.trim.contigs.good.unique.summary, start=13862, \
   end=23444, maxhomop=6, processors=${params.ncrs})"

   ${params.mothur}/mothur "#summary.seqs(fasta=${project_id}.trim.contigs.good.unique.good.align, \
   count=${project_id}.trim.contigs.good.good.count_table, processors=${params.ncrs})"

   ${params.mothur}/mothur "#filter.seqs(fasta=${project_id}.trim.contigs.good.unique.good.align, \
   processors=${params.ncrs}, vertical=T, trump=.)"

   ${params.mothur}/mothur "#unique.seqs(fasta=${project_id}.trim.contigs.good.unique.good.filter.fasta, \
   count=${project_id}.trim.contigs.good.good.count_table)"

   ${params.mothur}/mothur "#pre.cluster(fasta=${project_id}.trim.contigs.good.unique.good.filter.unique.fasta, \
   count=${project_id}.trim.contigs.good.unique.good.filter.count_table, diffs=2, processors=${params.ncrs})"
    """
}
