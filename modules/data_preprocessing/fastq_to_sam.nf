process fastq_to_sam {
    cache 'lenient'
    publishDir params.out, mode:'symlink'

    input:
    tuple val(sample), val(sample_id), val(PU), val(PL), val(LB), path(R1), path(R2)
         
    output:
    tuple val(sample_id), path("uBams/${sample_id}_fastqtosam.bam"),  emit: fastq_to_sam_ch
        
    script:   
    """
    mkdir -p uBams
    cp ${R1} ${R2} uBams/

    docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/uBams:/data pipelinesinmegen/pipelines_inmegen:public \
    java -jar /usr/bin/picard.jar FastqToSam \
             -FASTQ /data/${R1} \
             -FASTQ2 /data/${R2} \
             -OUTPUT /data/${sample_id}_fastqtosam.bam \
             -READ_GROUP_NAME ${PU} \
             -SAMPLE_NAME ${sample} \
             -LIBRARY_NAME ${LB} \
             -PLATFORM_UNIT ${PU}.${sample} \
             -PLATFORM ${PL}
    
    cd uBams
    rm ${R1} ${R2}
    """
}
