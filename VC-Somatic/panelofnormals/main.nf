#!/usr/bin/env nextflow
// Workflow    : Creación de panel de normales para identificación de variantes somáticas
// Institution : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática
// Versión     : 0.1
// Docker image - pipelinesinmegen/pipelines_inmegen  -

nextflow.enable.dsl=2

include { mutect2forPanelofNormals      } from "../../modules/VC-Somatic/PON/mutect2PON.nf" 
include { genomicsDBimport              } from "../../modules/VC-Somatic/PON/genomicsDBimport.nf"
include { createSomaticPanelofNormals   } from "../../modules/VC-Somatic/PON/PON.nf"

// Print some pipeline information
println "Pipelines Inmegen"
println "Flujo de trabajo: Panel de normales para Mutect2"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Directorio de la referencia: $params.refdir"
println "Directorio de salida: $params.out"
println " "

workflow {

// Declare some parameters
    interval_list=file("${params.interval_list}")

// Data processing 
   Channel.fromPath("${params.sample_info}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sample = "${row.Sample}"
                         def bam = file("${row.Path}")
                 return [ sample, bam ]
               }.set { bam_ch }

   mutect2forPanelofNormals(bam_ch)
   
    vcf_files = mutect2forPanelofNormals.out.mtf_PON_out.toList()

   genomicsDBimport(vcf_files,"${params.project_name}",interval_list)
   
   createSomaticPanelofNormals(genomicsDBimport.out.genomics_db)
}
