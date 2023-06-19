#!/usr/bin/env nextflow

// Workflow    : Cuantificación y análisis de expresión diferencial
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdifección de bioinformatica (INMEGEN)
// Versión     : 0.1 
// Docker image - pipelines_inmegen:latest -

nextflow.enable.dsl=2

include { fastqc                   } from "../modules/qualitycontrol/fastqc.nf"
include { multiqc                  } from "../modules/qualitycontrol/multiqc.nf"
include { trim_Galore              } from "../modules/QDEA_RNAseq/trim_galore.nf"
include { kallisto                 } from "../modules/QDEA_RNAseq/kallisto.nf"
include { tximport_deseq2          } from "../modules/QDEA_RNAseq/DEA.nf"
include { tximport_q               } from "../modules/QDEA_RNAseq/tximport.nf"
         
// Imprimir algunos directorios importantes
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Cuantificación y Análisis de Expresión"
println "Imagen de docker: pipelines_inmegen:latest"
println " "
println "Datos crudos: $params.reads"
println "Información de las muestras: $params.csv_info"
println "Realizar análsis de expresión diferencial: $params.QDEA"
println "Referencia: $params.ref"
println "Directorio de salida: $params.out"
println " "

workflow qualitycontrol {

   data_fq = Channel.fromFilePairs("${params.reads}")
                    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
                       
   fastqc(data_fq)
   
   analisis_dir = "${params.out}"+"/fastqc"
   multiqc(fastqc.out.fq_files.collect(),analisis_dir)   
}

workflow {

 // Subworkflow for quality control

   qualitycontrol()

// Data preprocessing

   Channel.fromPath("${params.csv_info}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sample = "${row.Sample_name}"
                         def read1  = file("${row.R1}")
                         def read2  = file("${row.R2}")
                 return [ sample, read1, read2 ]
               }
          .set { read_pairs_ch}
  
    //read_pairs_ch.view()
    
    trim_Galore(read_pairs_ch)
    
     tg_dir = "${params.out}"+"/trimming_files"    
       
    multiqc(trim_Galore.out.trim_fq.collect(), tg_dir)
    
    kallisto(trim_Galore.out.trim_fq)
    
      sample_info=file("${params.csv_info}")
      klx_files="${params.out}"+"/kallisto_quants"

    if ("${params.QDEA}" == true){
      script_r=file("${params.rscript_DEA_dir}") 
    tximport_deseq2(kallisto.out.abundance_h5.collect(),sample_info,klx_files,script_r) 
    } 
    else { 
      script_r2=file("${params.rscript_q_dir}")
    tximport_q(kallisto.out.abundance_h5.collect(),sample_info,klx_files,script_r2) 
    }
}
