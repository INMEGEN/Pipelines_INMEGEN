process genotypeGVCFs {
    containerOptions "-v ${params.refdir}:/ref"
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/raw_vcfs", mode:'copy'

    input:
    tuple val(project_id), path(database)

    output:
    tuple val(project_id), path("${project_id}_raw_variants.vcf.gz"), path("${project_id}_raw_variants.vcf.gz.tbi"), emit: gvcfs_out

    script:
    """
    gatk GenotypeGVCFs \
     -R /ref/${params.refname} \
     -V gendb://${database} \
     -O ${project_id}_raw_variants.vcf.gz 
    """
}

