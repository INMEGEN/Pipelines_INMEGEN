process haplotypeCaller {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir_star}:/ref"
    publishDir params.out + "/raw_vcfs", mode:'copy'
    
    input:
    tuple val(sample), path(input_bam), path(input_bam_idx)

    output:
    tuple val(sample), path("${sample}_raw_variants.vcf.gz"), path("${sample}_raw_variants.vcf.gz.tbi"), emit: hc_output

    script:
    """
    gatk HaplotypeCaller \
        -R /ref/${params.refname_star} \
        -I ${input_bam} \
        -O ${sample}_raw_variants.vcf.gz 
    """
}
