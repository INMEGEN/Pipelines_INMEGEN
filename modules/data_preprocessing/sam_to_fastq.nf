process sam_to_fastq {
    cache 'lenient'
    publishDir params.out, mode:'symlink'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("sam_to_fq/${sample}_samtofastq.fastq"),     emit: sam_to_fastq_ch

    script:
    """
    mkdir sam_to_fq
    cp ${reads} sam_to_fq/

    docker run --cpus ${params.ncrs} -v \$PWD/sam_to_fq:/data pipelinesinmegen/pipelines_inmegen:latest \
    java -jar /usr/bin/picard.jar SamToFastq \
      -I /data/${reads} \
      -FASTQ /data/${sample}_samtofastq.fastq \
      -CLIPPING_ATTRIBUTE XT \
      -CLIPPING_ACTION 2 \
      -INTERLEAVE true \
      -NON_PF true

    rm sam_to_fq/${reads}
    """
}
