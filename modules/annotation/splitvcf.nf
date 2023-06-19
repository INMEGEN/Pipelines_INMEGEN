process splitVCFs {
    container 'pipelines_inmegen:latest'
    containerOptions "-v ${params.refdir}:/ref"
    cache 'lenient'
    publishDir params.out + "/vcfs_persample" , mode:'copy'

    input:
    tuple val(project_id), path(project_vcf_snps), path(project_vcf_indels)

    output:
    path("*.ann.vcf"),    emit: ann_vcf_persample

    script:
    """ 
    mkdir -p tmp/splitvcf
    
    for sample in `bcftools query -l ${project_vcf_snps}`
    do
    gatk SelectVariants -R /ref/${params.refname} -V ${project_vcf_snps} -O \${sample}_filtered_snps_${project_id}.ann.vcf -sn \${sample} --exclude-non-variants true --tmp-dir tmp/splitvcf 
    done

    for sample in `bcftools query -l ${project_vcf_indels}`
    do
    gatk SelectVariants -R /ref/${params.refname} -V ${project_vcf_indels} -O \${sample}_filtered_indels_${project_id}.ann.vcf -sn \${sample} --exclude-non-variants true --tmp-dir tmp/splitvcf
    done

    rm -r tmp/splitvcf
    """
}
