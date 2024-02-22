process getMetrics {
    cache 'lenient'
    publishDir params.out, mode:'copy'

    input:
    tuple val(sample), path(bam), path(bam_idx)

    output:
    tuple val(sample),
          path("metrics/${sample}_alignment_metrics.txt"), \
          path("metrics/${sample}_insert_metrics.txt"), \
          path("metrics/${sample}_insert_size_histogram.pdf"),  emit: metrics_qc_ch

    script:
    """
   mkdir -p metrics
   cp ${bam} metrics/

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/metrics:/data -v "${params.refdir}":/ref  pipelinesinmegen/pipelines_inmegen:public \
   java -jar /usr/bin/picard.jar CollectAlignmentSummaryMetrics \
        -R /ref/${params.refname} \
        -I /data/${bam} \
        -O /data/${sample}_alignment_metrics.txt

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/metrics:/data pipelinesinmegen/pipelines_inmegen:public \
   java -jar /usr/bin/picard.jar CollectInsertSizeMetrics \
        -I /data/${bam} \
        -O /data/${sample}_insert_metrics.txt \
        -H /data/${sample}_insert_size_histogram.pdf
   
   rm metrics/${bam}
    """
}
