process getMetrics {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir_star}:/ref"
    publishDir params.out + "/metrics", mode:'copy'

    input:
    tuple val(sample), path(bam), path(bam_idx)

    output:
    tuple path("${sample}_alignment_metrics.txt"), \
          path("${sample}_insert_metrics.txt"), \
          path("${sample}_insert_size_histogram.pdf"),  emit: metrics_qc_ch

    script:
    """
    picard CollectAlignmentSummaryMetrics \
        -R /ref/${params.refname_star} \
        -I ${bam} \
        -O ${sample}_alignment_metrics.txt

    picard CollectInsertSizeMetrics \
        -I ${bam} \
        -O ${sample}_insert_metrics.txt \
        -H ${sample}_insert_size_histogram.pdf
    """
}
