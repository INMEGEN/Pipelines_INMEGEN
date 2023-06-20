process calculateContamination {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:latest'    
    publishDir params.out + "/contamination_tables", mode:'symlink'

    input:
    tuple val(tumor_id), path(input_bam)
    file(interval_list)
    file(common_biallelic)
    file(common_biallelic_idx)

    output:
    tuple val(tumor_id), path("${tumor_id}_segments.tsv"), path("${tumor_id}_calculatecontamination.table"),  emit: cont_tables

    script:
    """
    gatk GetPileupSummaries \
         -I ${input_bam} \
         -V ${common_biallelic} \
         -L ${interval_list} \
         -O ${tumor_id}_getpileupsummaries.table

    gatk CalculateContamination \
        -I ${tumor_id}_getpileupsummaries.table \
        -tumor-segmentation ${tumor_id}_segments.tsv \
        -O ${tumor_id}_calculatecontamination.table
    """
}
