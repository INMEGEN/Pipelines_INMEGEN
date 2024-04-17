process variantQC {
   cache 'lenient'
   publishDir params.out, mode: 'copy'

   input:
   tuple val(project_id), path(join_vcf), path(vcf_index)

   output:
   tuple val(project_id), path("vcfs_metrics/${project_id}_variantQC.html"), emit: summary_QC

   script:
   """
   mkdir -p /vcfs_metrics
   cp ${join_vcf} ${vcf_index} /vcfs_metrics

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/vcfs_metrics:/data -v "${params.refdir}":/ref \
                ghcr.io/bimberlab/discvrseq:latest VariantQC \
                -R /ref/${params.refname} \
                -V ${join_vcf} \
                -O ${project_id}_variantQC.html \
                --threads ${params.ncrs}

   cd vcfs_metrics/
   rm ${join_vcf} ${vcf_index}
   """
}
