process splitNCigarReads{
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    containerOptions "-v ${params.refdir_star}:/ref"
    //publishDir params.out + "/splitncigarReads", mode:'symlink'

    input:
    tuple val(sample), path(sorted_bam), path(sorted_bam_idx)

    output:
    tuple val(sample), path("${sample}_output.bam"), emit: split_bam
    
    script:
    """
    gatk SplitNCigarReads \
      -R /ref/${params.refname_star} \
      -I ${sorted_bam} \
      -O ${sample}_output.bam 
    """
}

