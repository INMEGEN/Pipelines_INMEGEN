process tools_verion {
  cache 'lenient'
  container 'pipelinesinmegen/pipelines_inmegen:public'
  publishDir params.out + "/tools_verions", mode: 'copy'
  
  input:
  path(salmon_version)

  output:
  path("tools_versions.txt"), emit: versions 

  script:
  """
  echo "Versiones de las herramientas utilizadas en el flujo de analisis QDEA" > tools_versions.txt
  echo "\$(fastp -v 2>&1)" >> tools_versions.txt 
  echo "\$(fastqc -v)" >> tools_versions.txt
  echo "star \$(STAR --version)" >> tools_versions.txt
  qualimap -h 2>&1 | grep 'QualiMap v.' >> tools_versions.txt
  featureCounts -v 2>&1 | grep 'featureCounts v' >> tools_versions.txt
  cat ${salmon_version} >> tools_versions.txt
  echo "\$(multiqc --version 2>&1)" >> tools_versions.txt 
  """
}
