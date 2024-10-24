process postfilter {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'    
    publishDir params.out + "/postfiltered_vcfs" , mode:'copy'

    input:
    tuple val(sample_id), path(vcf_file), path(index_file)

    output:
    tuple val(sample_id), path("${sample_id}_postfilter.vcf.gz"), path("${sample_id}_postfilter.vcf.gz.tbi"),    emit: filt_pass_vcf

    script:
    """
    bcftools view -f "PASS" -e "FORMAT/DP < 10"  --output ${sample_id}_tmp_postfilter.vcf ${vcf_file}

    bgzip -c ${sample_id}_tmp_postfilter.vcf > ${sample_id}_tmp_postfilter.vcf.gz
    tabix ${sample_id}_tmp_postfilter.vcf.gz

    bcftools view -e "FORMAT/AD[1:1] < 10"  --output ${sample_id}_postfilter.vcf ${sample_id}_tmp_postfilter.vcf.gz

    bgzip -c ${sample_id}_postfilter.vcf > ${sample_id}_postfilter.vcf.gz
    tabix ${sample_id}_postfilter.vcf.gz
    """
}
