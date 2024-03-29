process mutect2forPanelofNormals {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref" 
    publishDir params.out + "/vcfsforPON", mode:'copy'
    
    input:
    tuple val(sample), path(input_bam)

    output:
    tuple val(sample), path("${sample}_for_pon.vcf.gz"), emit: mtf_PON_out
    path("${sample}_for_pon.vcf.gz.tbi")
    path("${sample}_for_pon.vcf.gz.stats")

    script:
    """
   gatk Mutect2 \
   -R /ref/${params.refname} \
   -I ${input_bam} \
   -O ${sample}_for_pon.vcf.gz \
   -tumor ${sample} \
   -max-mnp-distance 0 
    """
}
