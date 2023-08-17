#!/usr/bin/env nextflow
// Workflow     : Identificación de variantes de datos RNA-seq con GATK4
// Institución  : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer   : Subdirección de genómica poblacional y subdirección de bioinformática del INMEGEN
// Versión      : 0.1 
// Docker image : - pipelines_inmegen:public -

nextflow.enable.dsl=2

// Processes for this workflow
include { fastqc                        } from "../modules/qualitycontrol/fastqc.nf"
include { multiqc                       } from "../modules/qualitycontrol/multiqc.nf"
include { trim_Galore                   } from "../modules/qualitycontrol/trim_galore.nf"
include { star as align                 } from "../modules/VC-RNAseq/star.nf"
include { mergeSam                      } from "../modules/VC-RNAseq/mergesamfiles.nf"
include { markDuplicatesSpark           } from "../modules/common/markDuplicatesSpark.nf"
include { getMetrics                    } from "../modules/metricts/getmetrics.nf"
include { Metrics                       } from "../modules/metricts/metrics.nf"
include { splitNCigarReads              } from "../modules/VC-RNAseq/splitNCigarReads.nf"
include { bqsr                          } from "../modules/VC-RNAseq/bqsr_recal.nf"
include { analyzeCovariates             } from "../modules/metricts/analyzecovariates.nf"
include { haplotypeCaller               } from "../modules/VC-RNAseq/haplotypecaller.nf"
include { selectVariants                } from "../modules/VC-RNAseq/selectvariants.nf"
include { filterSnps                    } from "../modules/VC-RNAseq/filtersnps.nf"
include { filterIndels                  } from "../modules/VC-RNAseq/filterindels.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Indentificación de variantes de RNAseq (experimental)"
println "Contenedor: pipelinesinmegen/pipelines_inmegen:public "
println " "
println "Nombre del proyecto: $params.project_name"
println "Datos crudos: $params.reads"
println "Información de las muestras: $params.sample_sheet"
println "Varios lanes por muestra: $params.multiple_samples"
println "Referencia: $params.ref"
println "Directorio de salida: $params.out"
println " "

workflow qualitycontrol {

   data_fq = Channel.fromFilePairs("${params.reads}")
                    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
                       
   fastqc(data_fq)
   
   analisis_dir = "${params.out}"+"/fastqc"
   multiqc(fastqc.out.fq_files.collect(), analisis_dir, "raw_data")   
}

workflow {

// subflujo de trabajo para el análisis de calidad de las muestras
   qualitycontrol()
  
// Data preprocessing
   Channel.fromPath("${params.sample_sheet}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sampleID = "${row.SampleID}"
                         def sample = "${row.Sample_name}"
                         def RG = "${row.RG}"
                         def PU = "${row.PU}"
                         def read1 = file("${row.R1}")
                         def read2 = file("${row.R2}")
                 return [ sampleID, sample, RG, PU, read1, read2 ]
               }
          .set { read_pairs_ch }

   trim_Galore(read_pairs_ch)

     tg_dir = "${params.out}"+"/trimming_files"
   multiqc(trim_Galore.out.trim_fq.collect(), tg_dir, "trimming_data")

   align(trim_Galore.out.trim_fq)
 
    if ("${params.multiple_samples}" == true){ 
     xa = align.out.aligned_reads_ch.collect().flatten().collate( 2 )
     xa.map { a , b -> def key = a.toString().tokenize('_').get(0)
                       def keyS = a.toString().tokenize('_').get(2)
                       return tuple("${key}" + "_" + "${keyS}", b)
            }.groupTuple() | mergeSam
    markDuplicatesSpark(mergeSam.out.merged_sam_ch)
    } 
    else {
    markDuplicatesSpark(align.out.aligned_reads_ch)
    }

   getMetrics(markDuplicatesSpark.out.bam_for_variant_calling)

   Metrics(markDuplicatesSpark.out.bam_for_variant_calling)

   splitNCigarReads(markDuplicatesSpark.out.bam_for_variant_calling)

   bqsr(splitNCigarReads.out.split_bam)

   analyzeCovariates(bqsr.out.analyze_covariates)
   
// Proceso para el llamado de variantes
   haplotypeCaller(bqsr.out.recalibrated_bam)

   selectVariants(haplotypeCaller.out.hc_output,haplotypeCaller.out.hc_output_index)
  
   filterSnps(selectVariants.out.snps_ch)

   filterIndels(selectVariants.out.indels_ch)
}
