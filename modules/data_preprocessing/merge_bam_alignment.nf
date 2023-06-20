process merge_bam_alignment {
  cache 'lenient'
  publishDir params.out, mode:'copy'
  
  input: 
  tuple val(sample), path(reads_1), path(reads_2)
  
  output:
  tuple val(sample), path("merge_bam_algn/${sample}_merged.bam"),     emit: bam_alignment_ch

  script:
  """
  mkdir -p merge_bam_algn
  cp ${reads_1} merge_bam_algn/
  cp ${reads_2} merge_bam_algn/

  docker run --cpus ${params.ncrs} -v \$PWD/merge_bam_algn:/data -v "${params.refdir}":/ref pipelinesinmegen/pipelines_inmegen:latest \
  java -jar /usr/bin/picard.jar MergeBamAlignment \
      -ALIGNED /data/${reads_1} \
      -UNMAPPED /data/${reads_2} \
      -O /data/${sample}_merged.bam \
      -R /ref/${params.refname} \
      --CREATE_INDEX true \
      --ADD_MATE_CIGAR true \
      --CLIP_ADAPTERS false \
      --CLIP_OVERLAPPING_READS true \
      --INCLUDE_SECONDARY_ALIGNMENTS true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --ATTRIBUTES_TO_RETAIN XS
  
  rm merge_bam_algn/${reads_1}
  rm merge_bam_algn/${reads_2} 
  """
}
