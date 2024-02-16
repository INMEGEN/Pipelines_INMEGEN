process markDuplicatesSpark {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/dedup_sorted", mode:'copy'

    input:
    tuple val(pair_id), path(aligned_reads)

    output:
    tuple val(pair_id), path("${pair_id}_sorted_dedup.bam"), path("${pair_id}_sorted_dedup.bam.bai"),   emit: bam_for_variant_calling
    tuple val(pair_id), path("${pair_id}_dedup_metrics.txt"),   emit: dedup_qc_ch

    script:
    """
    mkdir -p markduplicates/${pair_id}

    gatk MarkDuplicatesSpark \
        -I ${aligned_reads} \
        -M ${pair_id}_dedup_metrics.txt \
        -O ${pair_id}_sorted_dedup.bam \
        --tmp-dir markduplicates/${pair_id}
    
    rm -r markduplicates/${pair_id}
    """
}
