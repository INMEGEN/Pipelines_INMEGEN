process sam_to_fastq {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "sam_to_fq", mode:'symlink'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_samtofastq.fastq"),     emit: sam_to_fastq_ch

    script:
    """
    picard SamToFastq \
           -I ${reads} \
           -FASTQ ${sample}_samtofastq.fastq \
           -CLIPPING_ATTRIBUTE XT \
           -CLIPPING_ACTION 2 \
           -INTERLEAVE true \
           -NON_PF true
    """
}
