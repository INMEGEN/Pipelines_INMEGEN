process fastq_to_sam {
    cache 'lenient'
    publishDir params.out, mode:'symlink'

    input:
    tuple val(sample_id), val(sample), val(RG), val(PU), path(read1), path(read2)
         
    output:
    tuple val(sample_id), path("uBams/${sample_id}_fastqtosam.bam"),  emit: fastq_to_sam_ch
        
    script:   
    """
    mkdir -p uBams
    cp ${read1} uBams/
    cp ${read2} uBams/

    docker run --cpus ${params.ncrs} -v \$PWD/uBams:/data pipelines_inmegen:latest \
    java -jar /usr/bin/picard.jar FastqToSam \
             -FASTQ /data/${read1} \
             -FASTQ2 /data/${read2} \
             -OUTPUT /data/${sample_id}_fastqtosam.bam \
             -READ_GROUP_NAME ${RG} \
             -SAMPLE_NAME ${sample} \
             -LIBRARY_NAME ${sample}"."${PU} \
             -PLATFORM_UNIT ${params.pl} \
             -PLATFORM ${params.pl}

    rm uBams/${read1}
    rm uBams/${read2}
    """
}
