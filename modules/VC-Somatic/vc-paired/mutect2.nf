process mutect2 {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/unfilteredVCFs", mode:'copy'

    input:
    tuple val(tumor_id), path(tumor_bam), val(normal_id), path(normal_bam)
    file(interval_list)
    file(panel_normales)
    file(panel_normales_idx)

    output:
    tuple val(tumor_id), path("${tumor_id}_unfiltered.vcf.gz"), path("${tumor_id}_read-orientation-model.tar.gz"), emit: unfilt_vcf
    path("${tumor_id}_unfiltered.vcf.*")

    script:
    """
   gatk Mutect2 \
   -R /ref/${params.refname} \
   -I ${tumor_bam} \
   -I ${normal_bam} \
   -normal ${normal_id} \
   -L ${interval_list} \
   -imr ALL \
   --germline-resource /ref/${params.onlygnomad} \
   -pon ${panel_normales} \
   --f1r2-tar-gz ${tumor_id}_f1r2.tar.gz \
   -O ${tumor_id}_unfiltered.vcf.gz

   gatk LearnReadOrientationModel -I ${tumor_id}_f1r2.tar.gz -O ${tumor_id}_read-orientation-model.tar.gz
    """
}
