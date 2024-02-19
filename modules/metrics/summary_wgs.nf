process summary_wgs {
   cache 'lenient'
   container 'pipelinesinmegen/pipelines_inmegen:public'
   publishDir params.out + "/summary_metrics", mode: 'copy'

   input:
   val(list)
   val(project_id)
   path(summary_dir)

   output:
   tuple val(project_id), path("${project_id}_summary_QCmetrics.txt"), emit: summary_file

   script:
   """
   echo	"Sample id	Depth mean	N alingned reads	N aligned reads on target	On target reads percent	Total bases(in bed)	N bases on target	Percent of bases on target" >> ${project_id}_summary_QCmetrics.txt

   cat ${summary_dir}/*.txt >> ${project_id}_summary_QCmetrics.txt
   """
}
