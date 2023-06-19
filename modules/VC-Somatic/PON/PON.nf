process createSomaticPanelofNormals {
    cache 'lenient'
    container 'pipelines_inmegen:latest'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/panelofNormals", mode:'copy'
       
    input:
    tuple val(project_id), path(database)

    output:
    path("${project_id}_PON.vcf.gz")      , emit: PON_out
    path("${project_id}_PON.vcf.gz.tbi")

    script:
    """
    gatk CreateSomaticPanelOfNormals -R /ref/${params.refname} --germline-resource /ref/${params.onlygnomad} --output ${project_id}_PON.vcf.gz -V gendb://${database}
    """
}
