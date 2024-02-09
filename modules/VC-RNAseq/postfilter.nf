process postfiltervcf {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/post_filtered_vcfs" , mode:'copy'

    input:
    tuple val(sample), path(vcf_file), path(vcf_index)

    output:
    tuple val("${sample}"), path("${sample}_postfilter.vcf.gz"), path("${sample}_postfilter.vcf.gz.tbi"), emit: filt_pass_vcf

    script:
    """
    bcftools view -f "PASS" -e "FORMAT/DP < 10" -o ${sample}_postfilter.vcf.gz ${vcf_file}
    tabix ${sample}_postfilter.vcf.gz
    """
}
