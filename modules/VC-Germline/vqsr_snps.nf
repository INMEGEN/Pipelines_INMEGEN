process vqsrsnps {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/filtered_vcfs", mode:'copy'

    input:
    tuple val(project_id), path(raw_snps), path(raw_snps_idx)

    output:
    tuple val(project_id), path("${project_id}_filtered_snps.vcf.gz"), path("${project_id}_filtered_snps.vcf.gz.tbi"),    emit: snps_filt_ch

    script:
    """
   mkdir -p vqsr/tmp

   gatk VariantRecalibrator \
   -R /ref/${params.refname} \
   -V ${raw_snps} \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /ref/hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 /ref/1000G_omni2.5.hg38.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /ref/Homo_sapiens_assembly38.dbsnp138.vcf \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
   -mode SNP \
   -O ${project_id}_snps_output.recal \
   --tranches-file ${project_id}_snps_output.tranches \
   --tmp-dir vqsr/tmp

   gatk ApplyVQSR \
   -R /ref/${params.refname} \
   -V ${raw_snps} \
   -O ${project_id}_filtered_snps.vcf.gz \
   -ts-filter-level 99.5 \
   --tranches-file ${project_id}_snps_output.tranches \
   --recal-file ${project_id}_snps_output.recal \
   -mode SNP \
   --create-output-variant-index true \
   --tmp-dir vqsr/tmp
   
   rm -r vqsr/tmp
    """
}
