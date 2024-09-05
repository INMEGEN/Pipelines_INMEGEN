process mark_duplicates {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/markedDuplicates", mode:'symlink'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_markilluminaadapters.bam"),          emit: mark_duplicates_bam
    tuple val(sample), path("${sample}_markilluminaadapters_metrics.txt"),  emit: mark_duplicates_metrics

    script:
    """
    picard MarkIlluminaAdapters \
           -I ${reads} \
           -O ${sample}_markilluminaadapters.bam \
           -M ${sample}_markilluminaadapters_metrics.txt
    """
}
