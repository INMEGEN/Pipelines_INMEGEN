process haplotypeCallerERC {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/RAW_gvcfs", mode:'copy'

    input:
    tuple val(pair_id), path(input_bam)

    output:
    tuple val(pair_id), path("${pair_id}_raw_variants.g.vcf.gz"), emit: hc_erc_out
    path("${pair_id}_raw_variants.g.vcf.gz.tbi"), emit: hc_erc_index

    script:
    """
    gatk HaplotypeCaller \
        -R /ref/${params.refname} \
        -I $input_bam \
        -O ${pair_id}_raw_variants.g.vcf.gz \
        -max-mnp-distance 0 \
        -ERC GVCF \
        -L /ref/${params.intervalname}  
    """
}
