process align {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/aligned_reads", mode:'symlink'

    input:
    tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path(R1), path(R2)

    output:
    tuple val(sample_id), path("${sample_id}_aligned_reads.sam"),     emit: aligned_reads_ch

    script:
    readGroup = \
        "@RG\\tID:${PU}\\tPU:${PU}.${sample}\\tPL:${PL}\\tLB:${LB}\\tSM:${sample}"

    """
    bwa mem \
        -K 100000000 \
        -v 3 \
        -t ${params.ncrs} \
        -Y \
        -R \"${readGroup}\" \
        /ref/${params.refname} \
        ${R1} \
        ${R2} > ${pair_id}_aligned_reads.sam
    """
}
