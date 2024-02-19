process metricswes {
   cache 'lenient'
   container 'pipelinesinmegen/pipelines_inmegen:public'
   publishDir params.out + "/metrics", mode: 'copy'

   input:
   tuple val(sample_id), path(input_bam), path(bam_idx)
   file(bed_file)
   file(bed_file_w)

   output:
   tuple val(sample_id), path("summary/${sample_id}_QCmetrics.txt"), emit: summary_file
   path("${sample_id}*")

   script:
   """
   mkdir summary

   totalcounts=\$(samtools view -q 1 -F 3840 -c ${input_bam})
   onbedcounts=\$(samtools view  -q 1 -F 3840 -L ${bed_file} -c ${input_bam})
   onbedcounts2=\$(samtools view  -q 1 -F 3840 -L ${bed_file_w} -c ${input_bam})   

   ontarget=\$(awk "BEGIN {x=\$totalcounts;y=\$onbedcounts;print y/x}")
   ontargetp=\$(awk "BEGIN {x=\$ontarget;y=100;print x*y}")

   ontarget2=\$(awk "BEGIN {x=\$totalcounts;y=\$onbedcounts2;print y/x}")
   ontargetp2=\$(awk "BEGIN {x=\$ontarget2;y=100;print x*y}")

   mosdepth -t ${params.ncrs} -b ${bed_file} --thresholds 0,1 ${sample_id} ${input_bam}
   zcat ${sample_id}.thresholds.bed.gz | awk '{print \$(NF-1),"\\t",\$NF}' > ${sample_id}_thresholds.txt
   
   dpth=\$(grep "total_region" ${sample_id}.mosdepth.summary.txt | awk -F'\t' '{ print \$4"x" }')

   totalbases=\$(awk '{split(\$0,a,"\\t"); sum += a[1]} END {print sum}' ${sample_id}_thresholds.txt)
   ontargetbases=\$(awk '{split(\$0,a,"\\t"); sum += a[2]} END {print sum}' ${sample_id}_thresholds.txt)
    
   Bontarget=\$(awk "BEGIN {x=\$ontargetbases;y=\$totalbases;print x/y}")
   Bontargetp=\$(awk "BEGIN {x=\$Bontarget;y=100;print x*y}")   

   cd summary
   echo -e "${sample_id}	\$dpth	\$totalcounts	\$onbedcounts	\$onbedcounts2	\$ontargetp	\$ontargetp2	\$totalbases	\$ontargetbases	\$Bontargetp" >> ${sample_id}_QCmetrics.txt    
   """
}
