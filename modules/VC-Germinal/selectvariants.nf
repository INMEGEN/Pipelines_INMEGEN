process selectVariants {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/raw_vcfs", mode:'symlink'

    input:
    tuple val(project_id), path(variants)
    path(variants_index)

    output:
    tuple val(project_id), path("${project_id}_raw_snps.vcf"),    emit: snps_ch
    tuple val(project_id), path("${project_id}_raw_indels.vcf"),  emit: indels_ch

    script:
    """
    gatk SelectVariants \
        -R /ref/${params.refname} \
        -V ${variants} \
        -select-type SNP \
        -O ${project_id}_raw_snps.vcf 

    gatk SelectVariants \
        -R /ref/${params.refname} \
        -V ${variants} \
        -select-type INDEL \
        -O ${project_id}_raw_indels.vcf 
    """
}
