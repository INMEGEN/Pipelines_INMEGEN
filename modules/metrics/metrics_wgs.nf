process metricswgs {
   cache 'lenient'
   container 'pipelinesinmegen/pipelines_inmegen:public'
   publishDir params.out + "/metrics", mode: 'copy'

   input:
   tuple val(sample_id), path(input_bam), path(bam_idx)

   output:
   tuple val(sample_id), path("summary/${sample_id}_QCmetrics.txt"), emit: summary_file
   path("${sample_id}*")

   script:
   """
   mkdir summary

   samtools coverage -w 32 -o ${sample_id}_cov_hist.txt ${input_bam}   

   mosdepth -t ${params.ncrs} -n --fast-mode -b 500 --thresholds 0,1 ${sample_id} ${input_bam}
   zcat ${sample_id}.thresholds.bed.gz | awk '{print \$(NF-1),"\\t",\$NF}' > ${sample_id}_thresholds.txt

   dpth=\$(grep "total_region" ${sample_id}.mosdepth.summary.txt | awk -F'\t' '{ print \$4"x" }')

   totalbases=\$(awk '{split(\$0,a,"\\t"); sum += a[1]} END {print sum}' ${sample_id}_thresholds.txt)
   ontargetbases=\$(awk '{split(\$0,a,"\\t"); sum += a[2]} END {print sum}' ${sample_id}_thresholds.txt)

   Bontarget=\$(awk "BEGIN {x=\$ontargetbases;y=\$totalbases;print x/y}")
   Bontargetp=\$(awk "BEGIN {x=\$Bontarget;y=100;print x*y}")

   cd summary
   echo -e "${sample_id}	\$dpth	\$totalbases	\$ontargetbases	\$Bontargetp" >> ${sample_id}_QCmetrics.txt
   """
}
