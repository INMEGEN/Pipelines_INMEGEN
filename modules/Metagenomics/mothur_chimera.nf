process mothur_chimera {
    cache 'lenient'
    conda params.conda_env
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/mothur_chimera", mode:'copy'

    input:
    tuple val(project_id), path(fasta_file), path(count_table)    

    output:
    tuple val(project_id), path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta"), \
    path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table"), emit: mth_chm
    path("${project_id}.*")

    script:
    """
    ${params.mothur}/mothur "#chimera.vsearch(fasta=${fasta_file}, count=${count_table}, dereplicate=t)"
    
    ${params.mothur}/mothur "#summary.seqs(fasta=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, \
    count=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table)"
    """
}
