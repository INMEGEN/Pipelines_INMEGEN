process mothur_alignment {
    cache 'lenient'
    conda params.conda_env
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    //containerOptions "-v ${params.ref_silva}:/ref_silva"
    publishDir params.out + "/mothur_align", mode:'copy'

    input:
    tuple val(project_id), path(fasta_file), path(count_table)    
    
    output:
    tuple val(project_id), path("${project_id}.trim.contigs.good.unique.align"), emit: mth_aling
    path("${project_id}.*")

    script:
    """
   ${params.mothur}/mothur "#align.seqs(fasta=${fasta_file}, reference=${params.ref_silva}/silva.nr_v132.align, processors=${params.ncrs})"
    """
}
