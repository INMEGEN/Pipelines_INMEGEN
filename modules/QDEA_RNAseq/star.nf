process star {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  containerOptions "-v ${params.refdir_star}:/ref"
  publishDir params.out +"/alignments", mode: 'copy'

  input:
  tuple val(sample), path(R1), path(R2)

  output:
  tuple val(sample), path("${sample}_Aligned.sortedByCoord.out.bam")    , emit: aligned_ch
  tuple val(sample), path("${sample}_Aligned.toTranscriptome.out.bam")  , emit: aligned_toTrascript_ch
  tuple val(sample), path("countsU/${sample}_CountsPerGene.out.tab")    , emit: cuentas_uns
  tuple val(sample), path("countsR/${sample}_CountsPerGenerf.out.tab")  , emit: cuentas_rf
  path("counts/*")
  path("statistics/*")

  script:
  prefix = "${sample}" + "_"
  RG = "ID:${params.project_name}	SM:${sample}"
  """
    STAR --runMode alignReads \
         --genomeDir /ref/ \
         --sjdbGTFfile /ref/${params.gtfname} \
         --runThreadN ${params.ncrs} \
         --runDirPerm All_RWX \
         --readFilesIn ${R1} ${R2} \
         --outFileNamePrefix $prefix \
         --outReadsUnmapped None \
         --twopassMode Basic \
         --twopass1readsN -1 \
         --outSAMattrRGline $RG \
         --readFilesCommand "gunzip -c" \
         --quantMode TranscriptomeSAM GeneCounts \
         --outSAMunmapped Within \
         --outSAMtype BAM SortedByCoordinate

  echo "Gene_id	${sample}" > ${sample}_CountsPerGene.out.tab

  sed '1,4d' ${sample}_ReadsPerGene.out.tab | awk -F'\t' '{print \$1,"\\t",\$2}' >> ${sample}_CountsPerGene.out.tab

  echo "Gene_id	${sample}" > ${sample}_CountsPerGenerf.out.tab

  sed '1,4d' ${sample}_ReadsPerGene.out.tab | awk -F'\t' '{print \$1,"\\t",\$4}' >> ${sample}_CountsPerGenerf.out.tab

  mkdir counts countsR countsU statistics

  mv ${sample}_CountsPerGene.out.tab countsU/

  mv ${sample}_CountsPerGenerf.out.tab countsR/

  mv *.tab* counts/

  mv *_Log* statistics/
  """
}
