process postfilter {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:local'    
    publishDir params.out + "/filtered_vcfs" , mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file)

    output:
    tuple val(sample_id), path("${sample_id}_postfilter.vcf.gz"),    emit: filt_pass_vcf

    script:
    """
    bcftools view -f "PASS" -e "FORMAT/DP < 10"  --output ${sample_id}_tmp_postfilter.vcf.gz ${vcf_file}

    bcftools view -e "FORMAT/AD[*:1] < 10"  --output ${sample_id}_postfilter.vcf.gz ${sample_id}_tmp_postfilter.vcf.gz
    """
}
