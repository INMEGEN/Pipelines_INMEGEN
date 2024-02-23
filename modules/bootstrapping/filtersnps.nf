process filterSnps {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir_star}:/ref"
    publishDir params.out + "/filtered_vcfs", mode:'copy'

    input:
    tuple val(sample), path(raw_snps), path(raw_snps_idx)

    output:
    tuple val(sample), path("${sample}_filtered_snps.vcf.gz"), path("${sample}_filtered_snps.vcf.gz.tbi"),  emit: filtered_snps

    script:
    """
    gatk VariantFiltration \
        -R /ref/${params.refname_star} \
        -V ${raw_snps} \
        -O ${sample}_filtered_snps.vcf.gz \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    """
}
