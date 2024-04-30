process joinvcfs {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/filtered_vcfs", mode:'copy'

    input:
    tuple val(project_id), path(vcf_snps)  , path(vcf_snps_index)
    tuple val(project_id), path(vcf_indels), path(vcf_indels_index)

    output:
    tuple val(project_id), path("${project_id}_filtered.vcf.gz"), path("${project_id}_filtered.vcf.gz.tbi"),  emit: join_vars_filt

    script:
    """
    bcftools concat -a -Oz -o ${project_id}_filtered.vcf.gz ${vcf_snps} ${vcf_indels}
    tabix  ${project_id}_filtered.vcf.gz
    """
}
