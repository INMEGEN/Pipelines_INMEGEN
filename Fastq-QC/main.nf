#!/usr/bin/env nextflow
// Workflow    : Control de calidad de archivos FASTQ
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdifección de bioinformatica del INMEGEN
// Versión     : 0.1
// Docker image - pipelinesinmegen/pipelines_inmegen -

nextflow.enable.dsl=2

// Processes for this workflow
include { fastqc                             } from "../modules/qualitycontrol/fastqc.nf"
include { fastqScreen                        } from "../modules/qualitycontrol/fastq_screen.nf"
include { multiqc                            } from "../modules/qualitycontrol/multiqc.nf"

// Some useful information
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Control de calidad de los archivos fastq"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Directorio de salida: $params.out"
println " "


workflow {

// Declare some parameters 
   config=file("${params.fastq_screen_config}")

// Data preprocessing

   Channel.fromPath("${params.sample_info}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sample = "${row.Sample}"
                         def R1 = file("${row.R1}")
                         def R2 = file("${row.R2}")
                 return [ sample, R1, R2 ]
               }
          .set { read_pairs_ch }

    fastqc(read_pairs_ch)

    fastqScreen(read_pairs_ch,config)

// Summary

   multiqc(fastqScreen.out.screen_ch.collect(),fastqc.out.fq_files.collect(),"${params.out}")
}
