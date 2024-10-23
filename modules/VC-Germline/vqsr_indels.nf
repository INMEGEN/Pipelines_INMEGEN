process vqsrindels {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/filtered_vcfs", mode:'copy'
    
    input:
    tuple val(project_id), path(raw_indels), path(raw_indels_idx)

    output:
    tuple val(project_id), path("${project_id}_filtered_indels.vcf.gz"), path("${project_id}_filtered_indels.vcf.gz.tbi"),   emit: indels_filt_ch

    script:
    """   
   mkdir -p vqsr/tmp

   gatk VariantRecalibrator \
   -R /ref/${params.refname} \
   -V ${raw_indels} \
   -resource:mills,known=false,training=true,truth=true,prior=12.0 /ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /ref/Homo_sapiens_assembly38.dbsnp138.vcf \
   -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
   -mode INDEL \
   -O ${project_id}_indels_output.recal \
   --tranches-file ${project_id}_indels_output.tranches \
   --tmp-dir vqsr/tmp
   
   gatk ApplyVQSR \
   -R /ref/${params.refname} \
   -V ${raw_indels} \
   -O ${project_id}_filtered_indels.vcf.gz \
   -ts-filter-level 99.0 \
   --tranches-file ${project_id}_indels_output.tranches \
   --recal-file ${project_id}_indels_output.recal \
   -mode INDEL \
   --create-output-variant-index true \
   --tmp-dir vqsr/tmp

  rm -r vqsr/tmp
    """
}
