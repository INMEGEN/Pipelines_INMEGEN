#!/usr/bin/env nextflow
// Workflow    : Cuantificación y análisis de expresión diferencial de datos RNA-seq
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática (INMEGEN)
// Versión     : 0.1
// Docker image - pipelines_inmegen:public -

nextflow.enable.dsl=2

include { trim_Galore              } from "../modules/QDEA_RNAseq/trim_galore.nf"
include { kallisto                 } from "../modules/QDEA_RNAseq/kallisto.nf"
include { tximport_deseq2          } from "../modules/QDEA_RNAseq/DEA.nf"
include { tximport_q               } from "../modules/QDEA_RNAseq/tximport.nf"
include { multiqc                  } from "../modules/QDEA_RNAseq/multiqc.nf"
         
// Imprimir algunos directorios importantes
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Cuantificación y Análisis de Expresión Diferencial"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Tipo de análisis (true = cuantificación y DEG, false = sólo cuantificación): $params.QDEA"
println "Directorio de la referencia: $params.refdir"
println "Directorio de salida: $params.out"
println " "

workflow {

// Data preprocessing
   Channel.fromPath("${params.csv_info}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sample = "${row.Sample_name}"
                         def read1  = file("${row.R1}")
                         def read2  = file("${row.R2}")
                 return [ sample, read1, read2 ]
               }
          .set { read_pairs_ch}
    
    trim_Galore(read_pairs_ch)
    
     tg_dir = "${params.out}"+"/trimming_files"    
    multiqc(trim_Galore.out.trim_fq.collect(), tg_dir, "trimming_data")
    
    kallisto(trim_Galore.out.trim_fq)
    
      sample_info=file("${params.csv_info}")
      klx_files="${params.out}"+"/kallisto_quants"


    if ("${params.QDEA}" == true){
      script_r=file("${params.rscript_DEA_dir}") 
    tximport_deseq2(kallisto.out.abundance_h5.collect(),sample_info,klx_files,script_r) 
    
    multiqc(tximport_deseq2.out.R_sesion_info,"${params.out}")
    } 
    else { 
      script_r2=file("${params.rscript_q_dir}")
    tximport_q(kallisto.out.abundance_h5.collect(),sample_info,klx_files,script_r2) 

    multiqc(tximport_q.out.R_sesion_info,"${params.out}")
    }
}
