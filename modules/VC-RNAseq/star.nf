process star {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir}:/ref"
  publishDir params.out +"/alignments", mode: 'symlink' 

  input:
  tuple val(sample_id), val(sample), val(RG), val(PU), path(read_1), path(read_2)

  output:
    tuple val(sample_id), path("*.sam"),   emit: aligned_reads_ch
    //path("*sortedByCoord.out.bam.bai"),  emit: aligned_idx_ch

  script:
      readGroup = \
          "ID:${RG}\\tLB:${sample}.${PU}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${sample}"

      prefix = "${sample_id}" + "_"
  """
    STAR --runMode alignReads \
         --genomeDir /ref/ \
         --runThreadN ${params.ncrs} \
         --runDirPerm All_RWX \
         --readFilesIn ${read_1} ${read_2} \
         --outFileNamePrefix $prefix \
         --outReadsUnmapped None \
         --twopassMode Basic \
         --twopass1readsN -1 \
         --outSAMattrRGline ${readGroup} \
         --readFilesCommand "gunzip -c" \
         --outSAMunmapped Within \
         --outSAMtype SAM \
         --outSAMattributes NH HI AS nM

  """
}
