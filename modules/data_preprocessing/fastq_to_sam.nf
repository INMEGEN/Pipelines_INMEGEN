process fastq_to_sam {
    cache 'lenient'
    publishDir params.out, mode:'symlink'

    input:
    tuple val(sample_id), val(sample), val(RG), val(PU), path(reads)
         
    output:
    tuple val(sample_id), path("uBams/${sample_id}_fastqtosam.bam"),  emit: fastq_to_sam_ch
        
    script:   
    """
    mkdir -p uBams
    cp ${reads[0]} uBams/
    cp ${reads[1]} uBams/

    docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/uBams:/data pipelinesinmegen/pipelines_inmegen:public \
    java -jar /usr/bin/picard.jar FastqToSam \
             -FASTQ /data/${reads[0]} \
             -FASTQ2 /data/${reads[1]} \
             -OUTPUT /data/${sample_id}_fastqtosam.bam \
             -READ_GROUP_NAME ${RG} \
             -SAMPLE_NAME ${sample} \
             -LIBRARY_NAME ${sample}"."${PU} \
             -PLATFORM_UNIT ${params.pl} \
             -PLATFORM ${params.pl}

    rm uBams/${reads[0]}
    rm uBams/${reads[1]}
    """
}
