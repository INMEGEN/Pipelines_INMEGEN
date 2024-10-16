#!/usr/bin/env nextflow
// Workflow    : Cuantificación y análisis de expresión diferencial de datos RNA-seq
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática (INMEGEN)
// Versión     : 0.1
// Docker image - pipelinesinmegen/pipelines_inmegen -

nextflow.enable.dsl=2

include { fastp                    } from "../modules/QDEA_RNAseq/fastp.nf"
include { fastqc                   } from "../modules/QDEA_RNAseq/fastqc.nf"
include { star                     } from "../modules/QDEA_RNAseq/star.nf"
include { qualimap                 } from "../modules/QDEA_RNAseq/qualimap.nf"
include { kallisto                 } from "../modules/QDEA_RNAseq/kallisto.nf"
include { starcmatrix              } from "../modules/QDEA_RNAseq/starcmatrix.nf"
include { tximport_deseq2          } from "../modules/QDEA_RNAseq/DEA.nf"
include { tximport_q               } from "../modules/QDEA_RNAseq/tximport.nf"
include { multiqc                  } from "../modules/QDEA_RNAseq/multiqc.nf"
         
// Imprimir algunos directorios importantes
println "Pipelines INMEGEN"
println "Flujo de trabajo: Cuantificación y Análisis de Expresión Diferencial"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Tipo de análisis (true = cuantificación y DEG, false = sólo cuantificación): $params.QDEA"
println "Directorio de salida: $params.out"
println " "

workflow {
// Some Necessary files 
    rDEA=file("${params.r_DEA}")
    rQ=file("${params.rQ}")
    sample_info=file("${params.metadata}")

// Data preprocessing
    Channel.fromPath("${params.sample_info}" )
           .splitCsv(sep:"\t", header: true)
           .map { row ->  def sample = "${row.Sample}"
                          def R1  = file("${row.R1}")
                          def R2  = file("${row.R2}")
                 return [ sample, R1, R2 ]
                }
           .set { read_pairs_ch}

// Trimming data
    fastp(read_pairs_ch)

    fastqc(fastp.out.trim_fq)

// Aling to genome 
    star(fastp.out.trim_fq)

// Quality control
    qualimap(star.out.aligned_ch)

// Star matrix counts   
    starcmatrix(star.out.cuentas_rf.collect(), "${params.out}"+"/alignments/countsR")

// Quantifying transcript abundance in RNA-seq data
    kallisto(fastp.out.trim_fq)
    
// Differential expression analysis based on transcript abundance
    if ("${params.QDEA}" == true){
    tximport_deseq2(kallisto.out.abundance_h5.collect(),sample_info,"${params.out}"+"/kallisto_quants",rDEA) 

    multiqc(tximport_deseq2.out.mcounts_tpm,qualimap.out.qcmap.collect(),"${params.out}")
    } 
    else { 
    tximport_q(kallisto.out.abundance_h5.collect(),sample_info,"${params.out}"+"/kallisto_quants",rQ)
 
    multiqc(tximport_q.out.mcounts_tpm,qualimap.out.qcmap.collect(),"${params.out}")
    }
}
