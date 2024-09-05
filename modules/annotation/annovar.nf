process annovar {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:An1'
    containerOptions "-v ${params.annovar}/humandb/:/humandb"
    publishDir params.out + "/annovar" , mode:'copy'

    input:
    tuple val(sample), path(filtered_vcfs), path(filtered_vcfs_idx)


    output:
    tuple val(sample), path("${sample}_annovar.hg38_multianno.vcf.gz"), path("${sample}_annovar.hg38_multianno.vcf.gz.tbi"), emit: annovar_ch_vcf
    tuple val(sample), path("${sample}_annovar.hg38_multianno.txt"), emit: annovar_ch_txt
    path("*.avinput")

    script:
    """
    table_annovar ${filtered_vcfs} /humandb/ --buildver hg38 -out ${sample}_annovar \
    -remove -protocol refGene,ensGene,avsnp150,clinvar_20221231,gnomad312_genome,cosmic92_coding,dbnsfp33a -operation g,g,f,f,f,f,f --vcfinput --thread ${params.ncrs}
   
    diff=\$((\$(grep "Start" ${sample}_annovar.hg38_multianno.txt | wc -w)-\$(grep "QUAL" ${sample}_annovar.hg38_multianno.vcf | wc -w)))	
     
    sed -i "s/^Chr.*\$/\$(echo -e "\$(grep "Start" ${sample}_annovar.hg38_multianno.txt | cut -d\$'\t' -f1-\$diff)\t\$(grep "QUAL" ${sample}_annovar.hg38_multianno.vcf)")/" ${sample}_annovar.hg38_multianno.txt

    bgzip -c ${sample}_annovar.hg38_multianno.vcf > ${sample}_annovar.hg38_multianno.vcf.gz
    tabix ${sample}_annovar.hg38_multianno.vcf.gz
    """
}
