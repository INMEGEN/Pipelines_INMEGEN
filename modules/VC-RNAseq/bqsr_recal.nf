process bqsr {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.ref_dir_bqsr}:/db -v ${params.refdir_star}:/ref"
    //publishDir params.out + "/bqsr", mode:'symlink'

    input:
    tuple val(sample), path(reads)    

    output:
    tuple val(sample), path("${sample}_recalibration_data.table"), path("${sample}_post_recal_data.table"), emit: analyze_covariates
    tuple val(sample), path("${sample}_recalibrated.bam"), path("${sample}_recalibrated.bai"),              emit: recalibrated_bam

    script:
    """
    mkdir -p tmp/bqsr/${sample}

    gatk BaseRecalibrator \
    -I ${reads} \
    -R /ref/${params.refname_star} \
    --known-sites /db/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites /db/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites /db/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${sample}_recalibration_data.table \
    --tmp-dir tmp/bqsr/${sample}

    gatk ApplyBQSR \
    -R /ref/${params.refname_star} \
    -I ${reads} \
    -bqsr ${sample}_recalibration_data.table \
    -O ${sample}_recalibrated.bam \
    --tmp-dir tmp/bqsr/${sample}

    gatk BaseRecalibrator \
    -R /ref/${params.refname_star} \
    -I ${sample}_recalibrated.bam \
    --known-sites /db/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites /db/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites /db/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${sample}_post_recal_data.table \
    --tmp-dir tmp/bqsr/${sample}

    rm -r tmp/bqsr/${sample}
    """
}
