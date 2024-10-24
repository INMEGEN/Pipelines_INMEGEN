process variantQC {
   cache 'lenient'
   publishDir params.out , mode: 'copy'

   input:
   val(files_all)
   val(project_id)
   
   output:
   tuple val(project_id), path("vcfs_stats/${project_id}_filtered_all.vcf.gz"), \
                          path("vcfs_stats/${project_id}_filtered_all.vcf.gz.tbi"), emit: merge_vcf
   tuple val(project_id), path("vcfs_stats/${project_id}_variantQC.html")         , emit: summary_vcfQC

   script:
   """   
   mkdir -p vcfs_stats

   cp ${params.out}/filtered_vcfs/*_filtered.vcf.gz* vcfs_stats/

   cd vcfs_stats/

   ls *_filtered.vcf.gz > list_tmp.txt

   awk '{print "/data/"\$0}' list_tmp.txt > vcf_list.txt

   cd ..   

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/vcfs_stats:/data pipelinesinmegen/pipelines_inmegen:public bcftools merge \
                -m none \
                -l /data/vcf_list.txt \
                -Oz \
                -o /data/${project_id}_filtered_all.vcf.gz \
                --threads ${params.ncrs}

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/vcfs_stats:/data pipelinesinmegen/pipelines_inmegen:public tabix /data/${project_id}_filtered_all.vcf.gz

   docker run --cpus ${params.ncrs} --user="\$(id -u):\$(id -g)" -v \$PWD/vcfs_stats:/data -v "${params.refdir}":/ref ghcr.io/bimberlab/discvrseq:latest VariantQC \
                -R /ref/${params.refname} \
                -V /data/${project_id}_filtered_all.vcf.gz \
                -O /data/${project_id}_variantQC.html \
                --threads ${params.ncrs}
   """
}
