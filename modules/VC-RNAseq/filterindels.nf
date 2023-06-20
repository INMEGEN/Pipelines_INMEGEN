process filterIndels {
    container 'pipelinesinmegen/pipelines_inmegen:latest'
    containerOptions "-v ${params.refdir}:/ref"
    cache 'lenient'
    publishDir params.out + "/filtered_vcfs", mode:'symlink'

    input:
    tuple val(pair_id), path(raw_indels)

    output:
    tuple val(pair_id), \
          path("${pair_id}_filtered_indels.vcf"), \
          path("${pair_id}_filtered_indels.vcf.idx"), emit: filtered_indels

    script:
    """
    gatk VariantFiltration \
        -R /ref/${params.refname} \
        -V ${raw_indels} \
        -O ${pair_id}_filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0"
    """
}
