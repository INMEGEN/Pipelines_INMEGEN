process haplotypeCaller {
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    cache 'lenient'
    publishDir params.out + "/raw_vcfs", mode:'copy'
    
    input:
    tuple val(pair_id), path(input_bam)

    output:
    tuple val(pair_id), path("${pair_id}_raw_variants.vcf"), emit: hc_output
    path("${pair_id}_raw_variants.vcf.idx"), emit: hc_output_index

    script:
    """
    gatk HaplotypeCaller \
        -R /ref/${params.refname} \
        -I ${input_bam} \
        -O ${pair_id}_raw_variants.vcf 
    """
}
