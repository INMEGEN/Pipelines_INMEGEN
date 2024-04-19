process filterMutectCalls {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/filteredVCF" , mode:'copy'

    input:
    tuple val(tumor_id), val(input_vcf), path(orient_model), path(seg_table), path(cont_table)

    output:
    tuple val(tumor_id), path("${tumor_id}_filtered.vcf.gz"), path("${tumor_id}_filtered.vcf.*"), emit: filt_vcf

    script:
    """
    gatk FilterMutectCalls \
        -R /ref/${params.refname} \
        -V ${input_vcf} \
        -O ${tumor_id}_filtered.vcf.gz \
        --tumor-segmentation ${seg_table} \
        --contamination-table ${cont_table} \
        --ob-priors ${orient_model}
    """
}
