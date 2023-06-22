process annovar {
    cache 'lenient'
    publishDir params.out, mode:'copy'

    input:
    tuple val(project_id), path(filtered_vcfs)

    output:
    tuple val("${project_id}_annovar"), path("annovar/${project_id}_annovar*.vcf"), emit: annovar_ch
    path("annovar/*_multianno.txt")
    path("annovar/*.avinput")

    script:
    """
   mkdir -p annovar/
   cp ${filtered_vcfs} annovar/

   docker run --cpus ${params.ncrs} -v "\$PWD/annovar/":/data -v "${params.annovar}/humandb/":/humandb pipelines_inmegen:latest \
   perl /annovar/table_annovar.pl /data/${filtered_vcfs} /humandb/ --buildver hg38 -out /data/${project_id}_annovar \
   -remove -protocol refGene,ensGene,avsnp150,clinvar_20221231,gnomad312_genome,cosmic92_coding,dbnsfp33a -operation g,g,f,f,f,f,f -nastring . -vcfinput
    """
}

