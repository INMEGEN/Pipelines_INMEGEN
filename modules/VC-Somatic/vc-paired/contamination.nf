process calculateContamination {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'    
    publishDir params.out + "/contamination_tables", mode:'symlink'

    input:
    tuple val(tumor_id), path(tumor_bam), val(normal_id), path(normal_bam)
    file(interval_list)
    file(common_biallelic)
    file(common_biallelic_idx)

    output:
    tuple val(tumor_id), path("${tumor_id}_segments.tsv"), path("${tumor_id}_calculatecontamination.table"),  emit: cont_tables

    script:
    """
    gatk GetPileupSummaries \
         -I ${tumor_bam} \
         -V ${common_biallelic} \
         -L ${interval_list} \
         -O ${tumor_id}_getpileupsummaries.table

    gatk GetPileupSummaries \
         -I ${normal_bam} \
         -V ${common_biallelic} \
         -L ${interval_list} \
         -O ${normal_id}_getpileupsummaries.table

    gatk CalculateContamination \
        -I ${tumor_id}_getpileupsummaries.table \
        -matched ${normal_id}_getpileupsummaries.table \
        -segments ${tumor_id}_segments.tsv \
        -O ${tumor_id}_calculatecontamination.table
    """
}
