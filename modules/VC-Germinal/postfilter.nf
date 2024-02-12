process postfiltervcf {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/filtered_vcfs" , mode:'copy'

    input:
    tuple val(project_id), path(vcf_file), path(vcf_index)

    output:
    tuple val("${project_id}"), path("${project_id}_postfilter.vcf.gz"), path("${project_id}_postfilter.vcf.gz.tbi"), emit: filt_pass_vcf

    script:
    """
    bcftools view -f "PASS" -o ${project_id}_postfilter.vcf.gz ${vcf_file}
    tabix ${project_id}_postfilter.vcf.gz
    """
}
