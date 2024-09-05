process snpEff {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:An1'
    publishDir params.out + "/snpEff", mode:'copy'

    input:
    tuple val(sample_id), path(filtered_vcfs), path(filtered_vcfs_idx)
    
    output:
    tuple val("${sample_id}"), path("${sample_id}_snpEff.vcf.gz"), emit: snpeff_ch_vcf
    tuple path("*.csv"), path("*.txt"), emit: snpeff_ch_txt

    script:
    """
    snpEff -v GRCh38.99 -csvStats ${sample_id}_snpEff_stats.csv ${filtered_vcfs} > ${sample_id}_snpEff.vcf

    bgzip -c ${sample_id}_snpEff.vcf > ${sample_id}_snpEff.vcf.gz 
    """
}
