process align {
    cache 'lenient'
    container 'pipelines_inmegen:latest'
    containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out +"/Alignments", mode:'copy'

    input:
    tuple val(sample), path(reads) 
         
    output:
    tuple val(sample), path("${sample}_aligned_reads.sam"),     emit: aligned_reads
        
    script:
    """   
    bwa mem -K 100000000 -M -t ${params.ncrs} -p /ref/${params.refname} ${reads} >  ${sample}_aligned_reads.sam
    """
}
