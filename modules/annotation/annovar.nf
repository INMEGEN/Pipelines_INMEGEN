process annovar {
    cache 'lenient'
    publishDir params.out, mode:'copy'

    input:
    tuple val(sample_id), path(filtered_vcfs), path(filtered_vcfs_idx)


    output:
    tuple val("${sample_id}"), path("annovar/*_multianno.vcf.gz"), path("annovar/*_multianno.vcf.gz.tbi"), emit: annovar_ch_vcf
    tuple val("${sample_id}"), path("annovar/*_multianno.txt"), emit: annovar_ch_txt
    path("annovar/*.avinput")

    script:
    """
   mkdir -p annovar/
   cp ${filtered_vcfs} annovar/

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v "\$PWD/annovar/":/data -v "${params.annovar}/humandb/":/humandb pipelinesinmegen/pipelines_inmegen:public \
   perl /annovar/table_annovar.pl /data/${filtered_vcfs} /humandb/ --buildver hg38 -out /data/${sample_id}_annovar \
   -remove -protocol refGene,ensGene,avsnp150,clinvar_20221231,gnomad312_genome,cosmic92_coding,dbnsfp33a -operation g,g,f,f,f,f,f --vcfinput --thread ${params.ncrs}
   
   cd annovar/
 
   diff=\$((\$(grep "Start" ${sample_id}_annovar.hg38_multianno.txt | wc -w)-\$(grep "QUAL" ${sample_id}_annovar.hg38_multianno.vcf | wc -w)))	
     
   sed -i "s/^Chr.*\$/\$(echo -e "\$(grep "Start" ${sample_id}_annovar.hg38_multianno.txt | cut -d\$'\t' -f1-\$diff)\t\$(grep "QUAL" ${sample_id}_annovar.hg38_multianno.vcf)")/" ${sample_id}_annovar.hg38_multianno.txt

   bgzip -c ${sample_id}_annovar.hg38_multianno.vcf > ${sample_id}_annovar.hg38_multianno.vcf.gz
   tabix ${sample_id}_annovar.hg38_multianno.vcf.gz 
   
   rm ${filtered_vcfs}
    """
}
