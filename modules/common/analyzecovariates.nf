process analyzeCovariates{
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/bqsr", mode:'copy'

    input:
    tuple val(pair_id), path(recal_table), path(post_recal_table)

    output:
    tuple val(pair_id), path("${pair_id}_recalibration_plots.pdf"), emit: analyzed_covariates_ch

    script:
    """
    gatk AnalyzeCovariates \
        -before ${recal_table} \
        -after ${post_recal_table} \
        -plots ${pair_id}_recalibration_plots.pdf 
     """
}
