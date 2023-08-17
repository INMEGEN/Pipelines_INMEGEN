process align {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/aligned_reads", mode:'copy'

    input:
    tuple val(pair_id), val(sample), val(RG), val(PU), path(reads)

    output:
    tuple val(pair_id), path("${pair_id}_aligned_reads.sam"),     emit: aligned_reads_ch

    script:
    readGroup = \
        "@RG\\tID:${RG}\\tLB:${sample}.${PU}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${sample}"

    """
    bwa mem \
        -K 100000000 \
        -v 3 \
        -t ${params.ncrs} \
        -Y \
        -R \"${readGroup}\" \
        /ref/${params.refname} \
        ${reads[0]} \
        ${reads[1]} > ${pair_id}_aligned_reads.sam
    """
}
