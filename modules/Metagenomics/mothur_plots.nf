process mothur_plots {
    cache 'lenient'
    conda params.conda_env
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    //containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/mothur_plots", mode:'copy'

    input:
    tuple val(project_id), path(list_file), path(tax_file)
    path(shared_file)
    tuple val(project_id), path(phylip_file)
    file(Rplot)

    output:
    path("*.png")

    script:
    """
    Rscript ${Rplot} --l ${list_file} --tx ${tax_file} --t ${phylip_file} --s ${shared_file} --m ${params.metadata} --o ./ 
    """
}
