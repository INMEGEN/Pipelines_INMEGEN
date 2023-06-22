process genotypeGVCFs {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/join_vcfs", mode:'copy'

    input:
    tuple val(project_id), path(database)

    output:
    tuple val(project_id), path("${project_id}_raw_variants.vcf.gz"), emit: gvcfs_out
    path("${project_id}_raw_variants.vcf.gz.tbi"), emit: gvcfs_index

    script:
    """
    gatk GenotypeGVCFs \
     -R /ref/${params.refname} \
     -V gendb://${database} \
     -O ${project_id}_raw_variants.vcf.gz 
    """
}

