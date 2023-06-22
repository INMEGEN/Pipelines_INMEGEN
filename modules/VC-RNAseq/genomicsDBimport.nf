process genomicsDBimport {
    cache 'lenient'
    container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/Genomicsdb", mode:'symlink'

    input:
    val(sample_map)
    val(project_id)
    file(interval_list)

    output:
    tuple val(project_id), path("${project_id}_database"), emit: genomics_db

    script:
    """
    echo -e $sample_map | sed 's/^.//' | sed 's/.\$//' | sed 's/],/+/g' | tr '+' "\n" | tr ',' "\t" | tr -d "] [" >> cohort.sample_map

    mkdir -p genomicsdb/tmp

    gatk GenomicsDBImport \
       --genomicsdb-workspace-path ${project_id}_database \
       --batch-size ${params.batchsize} \
       --sample-name-map cohort.sample_map \
       --interval-merging-rule ALL \
       -L ${interval_list} \
       --merge-input-intervals ${params.interval} \
       --tmp-dir genomicsdb/tmp \
       --reader-threads ${params.GI_cores} \
       --max-num-intervals-to-import-in-parallel ${params.GI_cores}

    rm -r genomicsdb/tmp
    """
}
