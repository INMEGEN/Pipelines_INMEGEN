process getMetrics {
    cache 'lenient'
    publishDir params.out, mode:'copy'

    input:
    tuple val(pair_id), path(sorted_dedup_reads)

    output:
    tuple val(pair_id),
          path("metrics/${pair_id}_alignment_metrics.txt"), \
          path("metrics/${pair_id}_insert_metrics.txt"), \
          path("metrics/${pair_id}_insert_size_histogram.pdf"),  emit: metrics_qc_ch

    script:
    """
   mkdir -p metrics
   cp ${sorted_dedup_reads} metrics/

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/metrics:/data -v "${params.refdir}":/ref  pipelinesinmegen/pipelines_inmegen:public \
   java -jar /usr/bin/picard.jar CollectAlignmentSummaryMetrics \
        -R /ref/${params.refname} \
        -I /data/${sorted_dedup_reads} \
        -O /data/${pair_id}_alignment_metrics.txt

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/metrics:/data pipelinesinmegen/pipelines_inmegen:public \
   java -jar /usr/bin/picard.jar CollectInsertSizeMetrics \
        -I /data/${sorted_dedup_reads} \
        -O /data/${pair_id}_insert_metrics.txt \
        -H /data/${pair_id}_insert_size_histogram.pdf
   
   rm metrics/${sorted_dedup_reads}
    """
}
