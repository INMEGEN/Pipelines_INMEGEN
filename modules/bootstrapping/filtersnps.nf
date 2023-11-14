process filterSnps {
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    cache 'lenient'
    publishDir params.out + "/filtered_vcfs", mode:'copy'

    input:
    tuple val(pair_id), path(raw_snps)

    output:
    tuple val(pair_id), \
    path("${pair_id}_filtered_snps.vcf"), \
    path("${pair_id}_filtered_snps.vcf.idx"),     emit: filtered_snps

    script:
    """
    gatk VariantFiltration \
        -R /ref/${params.refname} \
        -V ${raw_snps} \
        -O ${pair_id}_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    """
}
