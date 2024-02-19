process star {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir_star}:/ref"
  publishDir params.out +"/alignments", mode: 'copy'

  input:
  tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path(R1), path(R2)

  output:
  tuple val(sample_id), path("*.sam"), emit: aligned_reads_ch  

  script:
      readGroup = \
          "ID:${PU}	PU:${PU}.${sample} PL:${PL}	LB:${LB}	SM:${sample}"

      prefix = "${sample_id}" + "_"
  """
    STAR --runMode alignReads \
         --genomeDir /ref/ \
         --runThreadN ${params.ncrs} \
         --runDirPerm All_RWX \
         --readFilesIn ${R1} ${R2} \
         --outFileNamePrefix $prefix \
         --outReadsUnmapped None \
         --twopassMode Basic \
         --twopass1readsN -1 \
         --outSAMattrRGline $readGroup \
         --readFilesCommand "gunzip -c" \
         --outSAMunmapped Within \
         --outSAMtype SAM \
         --outSAMattributes NH HI AS nM
  """
}
