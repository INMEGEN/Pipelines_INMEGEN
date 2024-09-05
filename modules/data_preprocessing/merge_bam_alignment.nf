process merge_bam_alignment {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir}:/ref"
  publishDir params.out + "/merge_bam_algn", mode:'symlink'
  
  input: 
  tuple val(sample), path(reads_1), path(reads_2)
  
  output:
  tuple val(sample), path("${sample}_merged.bam"), emit: bam_alignment_ch

  script:
  """
  picard MergeBamAlignment \
         -ALIGNED ${reads_1} \
         -UNMAPPED ${reads_2} \
         -O ${sample}_merged.bam \
         -R /ref/${params.refname} \
         --CREATE_INDEX true \
         --ADD_MATE_CIGAR true \
         --CLIP_ADAPTERS false \
         --CLIP_OVERLAPPING_READS true \
         --INCLUDE_SECONDARY_ALIGNMENTS true \
         --MAX_INSERTIONS_OR_DELETIONS -1 \
         --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
         --ATTRIBUTES_TO_RETAIN XS
  """
}
