process mergeSam {
   cache 'lenient'
   publishDir params.out + "/merged_sam", mode:'symlink'

   input:
   tuple val(key), path(bam_files)

   output:
   tuple val(key), path("${key}_merged.sam"),     emit: merged_sam_ch

   script:
   """    
   docker run --cpus ${params.ncrs} -v "\$PWD":/data pipelines_inmegen \
   java -jar /usr/bin/picard.jar MergeSamFiles \
   ${'-INPUT /data/'+bam_files.join(' -INPUT /data/')} \
   -OUTPUT /data/${key}_merged.sam -MSD true
   """
}
