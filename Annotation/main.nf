#!/usr/bin/env nextflow
// Workflow    : Anotación de variantes
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática del INMEGEN
// Versión     : 0.1 
// Docker image - pipelines_inmegen:latest -

nextflow.enable.dsl=2

include { annovar            } from "../modules/annotation/annovar_single.nf"
include { splitVCFs          } from "../modules/annotation/splitvcf_single.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Anotación de variantes"
println "Imagen de docker: pipelines_inmegen:latest"
println " "
println "Nombre del proyecto: $params.project_name"
println "Archivos vcfs : $params.vcfs"
println "Información de las muestras: $params.sample_sheet"
println "Directorio de salida: $params.out"
println " "

workflow {
 
   Channel.fromPath("${params.sample_sheet}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sampleID = "${row.Sample}"
                         def dir_vcf = file("${row.Path}")
                 return [ sampleID, dir_vcf ]
               }.set { vcfs_ch }

// Variant annotation
// Nota: recuerda obtener tu propia liga de annovar e incluirla en el Docker file 

    if ("${params.multiple_samples}" == true){ 

   annovar(vcfs_ch)
   splitVCFs(annovar.out.annovar_ch)

    } 
    else {

    annovar(vcfs_ch)

    }

}
