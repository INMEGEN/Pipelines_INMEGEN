process mutect2 {
    cache 'lenient'
    container 'pipelines_inmegen:latest'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/unfilteredVCFs", mode:'copy'

    input:
    tuple val(tumor_id), path(tumor_bam), val(normal_id), path(normal_bam)
    file(interval_list)
    file(panel_normales)
    file(panel_normales_idx)

    output:
    tuple val(tumor_id), path("${tumor_id}_unfiltered.vcf"), path("${tumor_id}_read-orientation-model.tar.gz"), emit: unfilt_vcf
    path("${tumor_id}_unfiltered.vcf.stats")

    script:
    """
   gatk Mutect2 \
   -R /ref/${params.refname} \
   -I ${tumor_bam} \
   -I ${normal_bam} \
   -normal ${normal_id} \
   -L ${interval_list} \
   -imr ALL \
   --interval-padding ${params.pading} \
   --germline-resource /ref/${params.onlygnomad} \
   -pon ${panel_normales} \
   --genotype-germline-sites ${params.germline_sites} \
   --genotype-pon-sites ${params.pon_sites} \
   --f1r2-tar-gz ${tumor_id}_f1r2.tar.gz \
   -O ${tumor_id}_unfiltered.vcf

   gatk LearnReadOrientationModel -I ${tumor_id}_f1r2.tar.gz -O ${tumor_id}_read-orientation-model.tar.gz
    """
}
