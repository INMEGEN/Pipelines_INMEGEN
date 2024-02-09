process filterIndels {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir_star}:/ref"
    publishDir params.out + "/filtered_vcfs", mode:'symlink'

    input:
    tuple val(sample), path(raw_indels), path(raw_indels_idx)

    output:
    tuple val(sample), path("${sample}_filtered_indels.vcf.gz"), path("${sample}_filtered_indels.vcf.gz.tbi"), emit: filtered_indels

    script:
    """
    gatk VariantFiltration \
        -R /ref/${params.refname_star} \
        -V ${raw_indels} \
        -O ${sample}_filtered_indels.vcf.gz \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}
