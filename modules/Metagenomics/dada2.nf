process dada2 {
    cache 'lenient'
    conda params.conda_env
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    //containerOptions "-v ${params.ref_dada2}:/ref_dada2"
    publishDir params.out + "/dada2", mode:'copy'

    input:
    val(raw_data)    
    file(Rdada2)

    output:
    path("results/*"), emit: dada2_out
    path("plots/*"),   emit: dada2_plots

    script:
    """
    mkdir plots
    mkdir results

   Rscript ${Rdada2} \
     --input ${raw_data}/ \
     --output results/ \
     --reference ${params.ref_dada2}/ \
     --plots plots/
    """
}
