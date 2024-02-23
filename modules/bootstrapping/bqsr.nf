process bqsr {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/bqsr", mode:'symlink'

    input:
    tuple val(sample), path(input_bam), path(input_bam_idx) 
    tuple val(sample), path(filtered_snps), path(filtered_snps_index) 
    tuple val(sample), path(filtered_indels), path(filtered_indels_index)  

    output:    
    tuple val(sample), path("${sample}_recalibration_data.table"), path("${sample}_post_recal_data.table"), emit: analyze_covariates
    tuple val(sample), path("${sample}_recalibrated.bam"),                                                  emit: recalibrated_bam
    path("${sample}_recalibrated.bai")

    script:
    """
    mkdir -p tmp/bqsr/${sample}

    gatk SelectVariants \
        --exclude-filtered \
        -V ${filtered_snps} \
        -O ${sample}_bqsr_snps.vcf

    gatk SelectVariants \
        --exclude-filtered \
        -V ${filtered_indels} \
        -O ${sample}_bqsr_indels.vcf

    gatk BaseRecalibrator \
        -R /ref/${params.refname} \
        -I ${input_bam} \
        --known-sites ${sample}_bqsr_snps.vcf \
        --known-sites ${sample}_bqsr_indels.vcf \
        -O ${sample}_recalibration_data.table \
        --tmp-dir tmp/bqsr/${sample} 

    gatk ApplyBQSR \
        -R /ref/${params.refname} \
        -I ${input_bam} \
        -bqsr ${sample}_recalibration_data.table \
        -O ${sample}_recalibrated.bam \
        --tmp-dir tmp/bqsr/${sample} 

    gatk BaseRecalibrator \
        -R /ref/${params.refname} \
        -I ${sample}_recalibrated.bam \
        --known-sites ${sample}_bqsr_snps.vcf \
        --known-sites ${sample}_bqsr_indels.vcf \
        -O ${sample}_post_recal_data.table \
        --tmp-dir tmp/bqsr/${sample}
 
    rm -r tmp/bqsr/${sample}
    """ 
}
