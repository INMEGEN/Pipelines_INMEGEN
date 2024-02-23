process haplotypeCaller {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/raw_vcfs", mode:'copy'

    input:
    tuple val(sample), path(input_bam)

    output:
    tuple val(sample), path("${sample}_raw_variants.vcf.gz"), path("${sample}_raw_variants.vcf.gz.tbi"), emit: hc_output

    script:
    """
    gatk HaplotypeCaller \
        -R /ref/${params.refname} \
        -I ${input_bam} \
        -O ${sample}_raw_variants.vcf.gz
    """
}
