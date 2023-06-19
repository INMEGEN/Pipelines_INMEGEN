process mark_duplicates {
    cache 'lenient'
    publishDir params.out, mode:'symlink'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("markedDuplicates/${sample}_markilluminaadapters.bam"),          emit: mark_duplicates_bam
    tuple val(sample), path("markedDuplicates/${sample}_markilluminaadapters_metrics.txt"),  emit: mark_duplicates_metrics

    script:
    """
    mkdir -p markedDuplicates
    cp ${reads} markedDuplicates/

    docker run --cpus ${params.ncrs} -v \$PWD/markedDuplicates:/data pipelines_inmegen:latest \
    java -jar /usr/bin/picard.jar MarkIlluminaAdapters \
      -I /data/${reads} \
      -O /data/${sample}_markilluminaadapters.bam \
      -M /data/${sample}_markilluminaadapters_metrics.txt

    rm markedDuplicates/${reads}
    """
}
