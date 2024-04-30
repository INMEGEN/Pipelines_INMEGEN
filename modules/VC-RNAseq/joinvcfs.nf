process joinvcfs {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/filtered_vcfs", mode:'copy'

    input:
    tuple val(sample), path(vcf_snps), path(vcf_snps_index), path(vcf_indels), path(vcf_indels_index)

    output:
    tuple val(sample), path("${sample}_filtered.vcf.gz"), path("${sample}_filtered.vcf.gz.tbi"),  emit: join_vars_filt

    script:
    """
    bcftools concat -a -Oz -o ${sample}_filtered.vcf.gz ${vcf_snps} ${vcf_indels}
    tabix  ${sample}_filtered.vcf.gz
    """
}
