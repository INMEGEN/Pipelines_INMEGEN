#!/usr/bin/env nextflow
// Workflow    : Creación de panel de normales para identificación de variantes somáticas GATK4
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdifección de bioinformatica del INMEGEN
// Versión     : 0.1 
// Docker image - pipelines_inmegen:latest -

nextflow.enable.dsl=2

include { mutect2forPanelofNormals      } from "../../modules/VC-Somatic/PON/mutect2PON.nf" 
include { genomicsDBimport              } from "../../modules/VC-Somatic/PON/genomicsDBimport.nf"
include { createSomaticPanelofNormals   } from "../../modules/VC-Somatic/PON/PON.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Panel de normales para Mutect2"
println "Imagen de docker: pipelines_inmegen:latest"
println " "
println "Nombre del proyecto: $params.project_name"
println "Archivos bam (muestras normales): $params.bams"
println "Información de las muestras: $params.sample_sheet"
println "Referencia: $params.ref"
println "Directorio de salida: $params.out"
println " "

workflow {
 
   Channel.fromPath("${params.sample_sheet}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sampleID = "${row.Sample}"
                         def bam_f = file("${row.Path}")
                 return [ sampleID, bam_f ]
               }.set { ready_bam_ch }

   interval_list=file("${params.interval_list}")

   mutect2forPanelofNormals(ready_bam_ch,interval_list)
   
    vcf_files = mutect2forPanelofNormals.out.mtf_PON_out.toList()   
    project_id="${params.project_name}"

   genomicsDBimport(vcf_files,project_id,interval_list)
   
   createSomaticPanelofNormals(genomicsDBimport.out.genomics_db)
}
