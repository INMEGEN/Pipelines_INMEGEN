process mothur_phylogenetic {
    cache 'lenient'
    conda '/scratch/home/dperez/programas/miniconda3/envs/nf-meta'
    //container 'pipelinesinmegen/pipelines_inmegen:public'
    //containerOptions "-v ${params.refdir}:/ref"
    publishDir params.out + "/mothur_phylogenetic_analysis", mode:'copy'

    input:
    tuple val(project_id), path(fasta_file), path(count_table)
    tuple val(project_id), path(list_file),  path(tax_file)
    file(Pyscript)

    output:
    tuple val(project_id), path("${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.rep.otu_modified.phylip.tre"), emit: mt_phylo
    path("${project_id}.*")

    script:
    """
   ${params.mothur}/mothur "#phylotype(taxonomy=${tax_file})"

   ${params.mothur}/mothur "#dist.seqs(fasta=${fasta_file}, output=lt)"

   ${params.mothur}/mothur "#clearcut(phylip=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.phylip.dist)"

   ${params.mothur}/mothur "#get.oturep(phylip=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.phylip.dist, \
   list=${list_file}, count=${count_table}, fasta=${fasta_file})"

   python ${Pyscript} ${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.rep.fasta
   
   ${params.mothur}/mothur "#dist.seqs(fasta=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.rep.otu_modified.fasta, output=lt)"

   ${params.mothur}/mothur "#clearcut(phylip=${project_id}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.rep.otu_modified.phylip.dist)"
    """
}
