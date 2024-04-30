process postfilter {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'    
    publishDir params.out + "/postfiltered_vcfs" , mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file), path(index_file)

    output:
    tuple val(sample_id), path("${sample_id}_postfilter.vcf.gz"),    emit: filt_pass_vcf

    script:
    """
    bcftools view -f "PASS" -e "FORMAT/DP < 10" -e "FORMAT/AD[1:1] < 10" -Oz --output ${sample_id}_postfilter.vcf.gz ${vcf_file}
    """
}
