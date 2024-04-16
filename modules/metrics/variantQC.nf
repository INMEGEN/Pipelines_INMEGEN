process variantQC {
   cache 'lenient'
   container 'ghcr.io/bimberlab/discvrseq:latest'
   containerOptions "-v ${params.refdir}:/ref"
   publishDir params.out + "/summary_metrics", mode: 'copy'

   input:
   tuple val(project_id), path(join_vcf), path(vcf_index)

   output:
   tuple val(project_id), path("${project_id}_variantQC.html"), emit: summary_QC

   script:
   """
   VariantQC \
     -R /ref/${params.refname} \
     -V ${join_vcf} \
     -O ${project_id}_variantQC.html \
     --threads ${params.ncrs}
   """
}
