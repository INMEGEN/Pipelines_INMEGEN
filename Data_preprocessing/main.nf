#!/usr/bin/env nextflow
// Workflow: Data_preprocessing for GATK short variant calling
// Institución: Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer: Subdirección de genómica poblacional y subdirección de bioinformática (INMEGEN)
// Versión: 0.1 
// Docker image - pipelines_inmegen:public -

nextflow.enable.dsl=2

include { fastqc                        } from "../modules/qualitycontrol/fastqc.nf"
include { multiqc                       } from "../modules/qualitycontrol/multiqc.nf"
include { trim_Galore                   } from "../modules/qualitycontrol/trim_galore.nf"
include { fastq_to_sam                  } from "../modules/data_preprocessing/fastq_to_sam.nf"
include { mark_duplicates               } from "../modules/data_preprocessing/mark_duplicates.nf"
include { sam_to_fastq                  } from "../modules/data_preprocessing/sam_to_fastq.nf"
include { align                         } from "../modules/data_preprocessing/align.nf"
include { merge_bam_alignment           } from "../modules/data_preprocessing/merge_bam_alignment.nf"
include { mergeSam                      } from "../modules/common/mergesamfiles.nf"
include { markDuplicatesSpark           } from "../modules/common/markDuplicatesSpark.nf"
include { bqsr                          } from "../modules/data_preprocessing/bqsr_recal.nf"
include { getMetrics                    } from "../modules/metricts/getmetrics.nf"
include { Metrics                       } from "../modules/metricts/metrics.nf"
include { analyzeCovariates             } from "../modules/metricts/analyzecovariates.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Preprocesamiento de datos para GATK4"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen:public"
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

// Subworkflow for quality control
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
          .set { read_pairs_ch}

   trim_Galore(read_pairs_ch)
    
     tg_dir = "${params.out}"+"/trimming_files"
   multiqc(trim_Galore.out.trim_fq.collect(), tg_dir, "trimming_data")

   fastq_to_sam(trim_Galore.out.trim_fq)
       // canal que junta la salida de fastq to sam
       ch_fqtsam=fastq_to_sam.out.fastq_to_sam_ch.collect().flatten().collate( 2 )
 
   mark_duplicates(fastq_to_sam.out.fastq_to_sam_ch)
   
   sam_to_fastq(mark_duplicates.out.mark_duplicates_bam)
          
   align(sam_to_fastq.out.sam_to_fastq_ch)
       // canales que juntan las salidas de fastq to sam y align para hacer merge por nombre de muestra
       // los acomoda en un formato [sampleID, path fastqtosam, path align bam]
       ch_aln=align.out.aligned_reads.collect().flatten().collate( 2 )
       ch_aln.join(ch_fqtsam).groupTuple().flatten().collate( 3 ).set{joinfmerge}
   
   merge_bam_alignment(joinfmerge)
  
     if ("${params.multiple_samples}" == true){ 

     xa = merge_bam_alignment.out.bam_alignment_ch.collect().flatten().collate( 2 )
     xa.map { a , b -> def key = a.toString().tokenize('_').get(0)
                       def keyS = a.toString().tokenize('_').get(2)
                       return tuple("${key}" + "_" + "${keyS}", b)
            }.groupTuple() | mergeSam

    markDuplicatesSpark(mergeSam.out.merged_sam_ch)

    } 
    else {

    markDuplicatesSpark(merge_bam_alignment.out.bam_alignment_ch)

    }
 
   getMetrics(markDuplicatesSpark.out.bam_for_variant_calling)

   Metrics(markDuplicatesSpark.out.bam_for_variant_calling)

   bqsr(markDuplicatesSpark.out.bam_for_variant_calling)

   analyzeCovariates(bqsr.out.analyze_covariates)
}
