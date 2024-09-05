process fastq_to_sam {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/uBams", mode:'symlink'

    input:
    tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path(R1), path(R2)
         
    output:
    tuple val(sample_id), path("${sample_id}_fastqtosam.bam"),  emit: fastq_to_sam_ch
        
    script:   
    """
    picard FastqToSam \
          -FASTQ ${R1} \
          -FASTQ2 ${R2} \
          -OUTPUT ${sample_id}_fastqtosam.bam \
          -READ_GROUP_NAME ${PU} \
          -SAMPLE_NAME ${sample} \
          -LIBRARY_NAME ${LB} \
          -PLATFORM_UNIT ${PU}.${sample} \
          -PLATFORM ${PL}
    """
}
