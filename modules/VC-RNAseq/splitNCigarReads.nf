process splitNCigarReads{
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/splitncigarReads", mode:'symlink'

    input:
    tuple val(pair_id), path(sorted_bam)

    output:
    tuple val(pair_id), path("${pair_id}_output.bam"), emit: split_bam
    
    script:
    """
    gatk SplitNCigarReads \
      -R /ref/${params.refname} \
      -I ${sorted_bam} \
      -O ${pair_id}_output.bam 
    """
}

