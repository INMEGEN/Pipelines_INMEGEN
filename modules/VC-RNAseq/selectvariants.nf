process selectVariants {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir_star}:/ref"
    publishDir params.out + "/raw_vcfs", mode:'copy'

    input:
    tuple val(sample), path(variants), path(variants_index)

    output:
    tuple val(sample), path("${sample}_raw_snps.vcf.gz"),   path("${sample}_raw_snps.vcf.gz.tbi"),    emit: snps_ch
    tuple val(sample), path("${sample}_raw_indels.vcf.gz"), path("${sample}_raw_indels.vcf.gz.tbi"),  emit: indels_ch

    script:
    """
    gatk SelectVariants \
        -R /ref/${params.refname_star} \
        -V ${variants} \
        -select-type SNP \
        -O ${sample}_raw_snps.vcf.gz

    gatk SelectVariants \
        -R /ref/${params.refname_star} \
        -V ${variants} \
        -select-type INDEL \
        -O ${sample}_raw_indels.vcf.gz 
    """
}
