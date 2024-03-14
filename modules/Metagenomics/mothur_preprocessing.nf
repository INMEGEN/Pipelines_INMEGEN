process mothur_preprocessing {
    cache 'lenient'
    conda params.conda_env
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    publishDir params.out + "/mothur_preprocessing", mode:'copy'

    input:
    val(raw_data)    
    val(trimming_dir)
    val(project_id)

    output:
    tuple val(project_id), path("${project_id}.trim.contigs.good.unique.fasta"), path("${project_id}.trim.contigs.good.count_table"),  emit: mth_pre
    path("${project_id}.*")

    script:
    """
   ${params.mothur}/mothur "#make.file(inputdir=${trimming_dir}/, outputdir=./, type=gz, prefix=${project_id})"

   ${params.mothur}/mothur "#make.contigs(file=${project_id}.files, inputdir=${trimming_dir}/, outputdir=./, processors=${params.ncrs})"

   ${params.mothur}/mothur "#summary.seqs(fasta=${project_id}.trim.contigs.fasta, processors=${params.ncrs})"

   ${params.mothur}/mothur "#screen.seqs(fasta=${project_id}.trim.contigs.fasta, count=${project_id}.contigs.count_table, summary=${project_id}.trim.contigs.summary, \
   maxambig=0, maxlength=275, maxhomop=8)"

   ${params.mothur}/mothur "#summary.seqs(fasta=${project_id}.trim.contigs.good.fasta, count=${project_id}.contigs.good.count_table, processors=${params.ncrs})"
   
   ${params.mothur}/mothur "#unique.seqs(fasta=${project_id}.trim.contigs.good.fasta, count=${project_id}.contigs.good.count_table)"
    """
}
