#!/usr/bin/env nextflow
// Workflow    : Anotación de variantes
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática del INMEGEN
// Versión     : 0.1 
// Docker image - pipelines_inmegen:latest -

nextflow.enable.dsl=2

include { annovar            } from "../modules/annotation/annovar.nf"
include { splitVCFs          } from "../modules/annotation/splitvcf.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Anotación de variantes"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Directorio de salida: $params.out"
println " "

workflow {
 
   Channel.fromPath("${params.sample_info}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sample = "${row.Sample}"
                         def path = file("${row.Path}")
                         def index = file("${row.Index}")
                         return [ sample, path, index ]
               }.set { vcfs_ch }

// Variant annotation

    if ("${params.multiple_samples}" == true){ 

   annovar(vcfs_ch)
   splitVCFs(annovar.out.annovar_ch_vcf,"annotated")

    } 
    else {

    annovar(vcfs_ch)

    }
}
