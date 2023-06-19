process annovar {
    cache 'lenient'
    publishDir params.out, mode:'copy'

    input:
    tuple val(project_id), path(filtered_snps)  , path(filtered_snps_index)
    tuple val(project_id), path(filtered_indels), path(filtered_indeles_index)

    output:
    tuple val("${project_id}_annovar"), path("annovar/${project_id}_annovar_snps*.vcf"), path("annovar/${project_id}_annovar_indels*.vcf"), emit: annovar_ch
    path("annovar/*_multianno.txt")
    path("annovar/*.avinput")

    script:
    """
   mkdir -p annovar/
   cp ${filtered_snps} ${filtered_snps_index} annovar/
   cp ${filtered_indels} ${filtered_indeles_index} annovar/

   docker run --cpus ${params.ncrs} -v "\$PWD/annovar/":/data -v "${params.annovar}/humandb/":/humandb pipelines_inmegen \
   perl /annovar/table_annovar.pl /data/${filtered_snps} /humandb/ --buildver hg38 -out /data/${project_id}_annovar_snps \
   -remove -protocol refGene,ensGene,avsnp150,clinvar_20221231,gnomad312_genome,cosmic92_coding,dbnsfp33a -operation g,g,f,f,f,f,f -nastring . -vcfinput

   docker run --cpus ${params.ncrs} -v "\$PWD/annovar/":/data -v "${params.annovar}/humandb/":/humandb pipelines_inmegen \
   perl /annovar/table_annovar.pl /data/${filtered_indels} /humandb/ --buildver hg38 -out /data/${project_id}_annovar_indels \
   -remove -protocol refGene,ensGene,avsnp150,clinvar_20221231,gnomad312_genome,cosmic92_coding,dbnsfp33a -operation g,g,f,f,f,f,f -nastring . -vcfinput
    """
}

