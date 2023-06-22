process splitVCFs {
    cache 'lenient'
    container 'pipelines_inmegen:latest'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/vcfs_persample" , mode:'copy'

    input:
    tuple val(project_id), path(project_vcf_snps)

    output:
    path("*.ann.vcf"),    emit: ann_vcf_persample

    script:
    """ 
    mkdir -p tmp/splitvcf
    
    for sample in `bcftools query -l ${project_vcf_snps}`
    do
    gatk SelectVariants -R /ref/${params.refname} -V ${project_vcf_snps} -O \${sample}_filtered_${project_id}.ann.vcf -sn \${sample} --exclude-non-variants true --tmp-dir tmp/splitvcf 
    done

    rm -r tmp/splitvcf
    """
}
