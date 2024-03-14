process krona {
    cache 'lenient'
    conda params.conda_env
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    //containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/mothur_plots", mode:'copy'

    input:
    tuple val(project_id), path(summary_file)
    file(Pyscript)

    output:
    path("${project_id}*"), emit: krona_ch

    script:
    """
    python ${Pyscript} ${summary_file} > ${project_id}_krona_output.xml

    ktImportXML -o ${project_id}_krona.html.html ${project_id}_krona_output.xml
    """
}
