process mergeSam {
   cache 'lenient'
   container 'pipelinesinmegen/pipelines_inmegen:public'
   publishDir params.out + "/merged_sam", mode:'symlink'

   input:
   tuple val(key), path(bam_files)

   output:
   tuple val(key), path("${key}_merged.sam"),     emit: merged_sam_ch

   script:
   """
   picard MergeSamFiles \
   ${'-INPUT '+bam_files.join(' -INPUT ')} \
   -OUTPUT ${key}_merged.sam   
   """
}
