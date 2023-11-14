process mergeSam {
   cache 'lenient'
   publishDir params.out, mode:'symlink'

   input:
   tuple val(key), path(bam_files)

   output:
   tuple val(key), path("merged_sam/${key}_merged.sam"),     emit: merged_sam_ch

   script:
   """    
   mkdir -p merged_sam
   cp ${bam_files} merged_sam/

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/merged_sam:/data pipelinesinmegen/pipelines_inmegen:public \
   java -jar /usr/bin/picard.jar MergeSamFiles \
   ${'-INPUT /data/'+bam_files.join(' -INPUT /data/')} \
   -OUTPUT /data/${key}_merged.sam -MSD true   

   rm merged_sam/${bam_files}
   """
}
