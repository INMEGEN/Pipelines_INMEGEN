process variantQC {
   cache 'lenient'
   publishDir params.out , mode: 'copy'

   input:
   tuple val(project_id), path(join_vcf), path(vcf_index)

   output:
   tuple val(project_id), path("summary_metrics/${project_id}_variantQC.html"), emit: summary_QC

   script:
   """
   mkdir -p summary_metrics

   cp ${join_vcf} ${vcf_index} summary_metrics/

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/summary_metrics:/data -v "${params.refdir}":/ref ghcr.io/bimberlab/discvrseq:latest VariantQC \
                -R /ref/${params.refname} \
                -V /data/${join_vcf} \
                -O /data/${project_id}_variantQC.html \
                --threads ${params.ncrs}

   cd summary_metrics/
   rm ${join_vcf} ${vcf_index}
   """
}
