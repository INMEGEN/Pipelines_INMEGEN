#!/usr/bin/env nextflow
// Workflow    : Cuantificación y análisis de expresión diferencial de datos RNA-seq
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática (INMEGEN)
// Versión     : 0.1
// Docker image - pipelinesinmegen/pipelines_inmegen -

nextflow.enable.dsl=2

include { fastp            } from "../modules/QDEA_RNAseq/fastp.nf"
include { fastqc           } from "../modules/QDEA_RNAseq/fastqc.nf"
include { star             } from "../modules/QDEA_RNAseq/star.nf"
include { starML           } from "../modules/QDEA_RNAseq/starML.nf"
include { qualimap         } from "../modules/QDEA_RNAseq/qualimap.nf"
include { salmon           } from "../modules/QDEA_RNAseq/salmon.nf"
include { salmonML         } from "../modules/QDEA_RNAseq/salmonML.nf"
include { featureCounts    } from "../modules/QDEA_RNAseq/featureCounts.nf"
include { fmcounts         } from "../modules/QDEA_RNAseq/fmcounts.nf"
include { DEA              } from "../modules/QDEA_RNAseq/DEA.nf"
include { tximport         } from "../modules/QDEA_RNAseq/tximport.nf"
include { multiqc          } from "../modules/QDEA_RNAseq/multiqc.nf"
include { salmon_version   } from "../modules/QDEA_RNAseq/salmon_version.nf"
include { tools_verion     } from "../modules/QDEA_RNAseq/tools_verion.nf"

         
// Imprimir algunos directorios importantes
println "Pipelines INMEGEN"
println "Flujo de trabajo: Cuantificación y Análisis de Expresión Diferencial"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Varios lanes por muestra (true = sí, false = no): $params.multiple_lanes"
println "Tipo de análisis (true = cuantificación y DEG, false = sólo cuantificación): $params.QDEA"
println "Directorio de salida: $params.out"
println " "

workflow {
// Some Necessary files 
    rDEA=file("${params.r_DEA}")
    rQ=file("${params.rQ}")
    sample_info=file("${params.metadata}")
    mqc_config=file("${params.mqc_config}")
    genes_file=file("${params.genes_file}")

// Data preprocessing
    Channel.fromPath("${params.sample_info}" )
           .splitCsv(sep:"\t", header: true)
           .map { row ->  def sample   = "${row.Sample}"
                          def sampleid = "${row.SampleID}"
                          def R1  = file("${row.R1}")
                          def R2  = file("${row.R2}")
                 return [ sample, sampleid, R1, R2 ]
                }
           .set { input_ch }

// Trimming data
    fastp(input_ch)

    fastqc(fastp.out.trim_fq)

// Quantifing abundance and differential expresion genes 
   if ("${params.multiple_lanes}" == true){

     xa = fastp.out.trim_fq.collect().flatten().collate( 3 ).groupTuple() 

     salmonML(xa,genes_file)

     sa = fastp.out.trim_fq.collect().flatten().collate( 3 ).groupTuple() | starML

     qualimap(starML.out.aligned_ch)

     featureCounts(starML.out.aligned_ch)

     fmcounts(featureCounts.out.fcounts.collect(),"${params.out}" + "/featureCounts/quants")

     tximport(salmonML.out.quant.collect(),sample_info,"${params.out}"+"/salmon",rQ)

     salmon_version(tximport.out.mcounts)

     tools_verion(salmon_version.out.sversion)

     if ("${params.QDEA}" == true){
 
     DEA(fmcounts.out.mcounts,sample_info,rDEA)

     reports_ch = DEA.out.rds.concat(qualimap.out.qcmap.collect()).concat(salmonML.out.quant.collect()).collect()

     multiqc(reports_ch,mqc_config,"${params.out}")
     }
     else { 
    
     reports_ch = qualimap.out.qcmap.collect().concat(salmonML.out.quant.collect()).collect()
     
     multiqc(reports_ch,mqc_config,"${params.out}")
     }
    }
    else {

     salmon(fastp.out.trim_fq,genes_file)

     star(fastp.out.trim_fq)
   
     featureCounts(star.out.aligned_ch)

     qualimap(star.out.aligned_ch)

     fmcounts(featureCounts.out.fcounts.collect(),"${params.out}"+"/featureCounts/quants")

     tximport(salmon.out.quant.collect(),sample_info,"${params.out}"+"/salmon",rQ)

     salmon_version(tximport.out.mcounts)

     tools_verion(salmon_version.out.sversion)

     if ("${params.QDEA}" == true){

     DEA(fmcounts.out.mcounts,sample_info,rDEA)

     reports_ch = DEA.out.rds.concat(qualimap.out.qcmap.collect()).concat(salmon.out.quant.collect()).collect()

     multiqc(reports_ch,mqc_config,"${params.out}")
     }
     else { 

     reports_ch = qualimap.out.qcmap.collect().concat(salmon.out.quant.collect()).collect()
 
     multiqc(reports_ch,mqc_config,"${params.out}")
     }
    }
}
