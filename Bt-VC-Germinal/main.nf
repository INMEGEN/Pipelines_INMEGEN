#!/usr/bin/env nextflow
// Workflow    : Identificación conjunta de variantes germinales utilizando bootstrapping
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática del INMEGEN
// Versión     : 0.1 
// Docker image - pipelinesinmegen/pipelines_inmegen:public -

nextflow.enable.dsl=2

// Processes for this workflow
include { fastqc                        } from "../modules/qualitycontrol/fastqc.nf"
include { multiqc                       } from "../modules/qualitycontrol/multiqc.nf"
include { trim_Galore                   } from "../modules/qualitycontrol/trim_galore.nf"
include { align                         } from "../modules/VC-Germinal/bwa_germinal.nf"
include { mergeSam                      } from "../modules/common/mergesamfiles.nf"
include { markDuplicatesSpark           } from "../modules/common/markDuplicatesSpark.nf"
include { getMetrics                    } from "../modules/metricts/getmetrics.nf"
include { Metrics                       } from "../modules/metricts/metrics.nf"
include { bqsr                          } from "../modules/bootstrapping/bqsr.nf"
include { analyzeCovariates             } from "../modules/metricts/analyzecovariates.nf"
include { haplotypeCallerERC            } from "../modules/VC-Germinal/haplotypecaller_erc.nf"
include { haplotypeCaller               } from "../modules/bootstrapping/haplotypecaller.nf"
include { genomicsDBimport              } from "../modules/VC-Germinal/genomicsDBimport.nf"
include { genotypeGVCFs                 } from "../modules/VC-Germinal/genotypegvcfs.nf"
include { selectVariants                } from "../modules/VC-Germinal/selectvariants.nf"
include { filterSnps                    } from "../modules/bootstrapping/filtersnps.nf"
include { filterIndels                  } from "../modules/bootstrapping/filterindels.nf"

// Some useful information
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Identificación conjunta de variantes germinales con GATK4"
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

workflow bootstrapping {

   take: data_1 
   
   main:
           
   haplotypeCaller(data_1)
   
   selectVariants(haplotypeCaller.out.hc_output,haplotypeCaller.out.hc_output_index)
   
   filterSnps(selectVariants.out.snps_ch)
   
   filterIndels(selectVariants.out.indels_ch)
         
   emit:
   filterSnps_out = filterSnps.out.filtered_snps
   filterIndels_out = filterIndels.out.filtered_indels     
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

// Subflujo de trabajo para generar estándar para  BQSR    
   bootstrapping(markDuplicatesSpark.out.bam_for_variant_calling)

   getMetrics(markDuplicatesSpark.out.bam_for_variant_calling)
   Metrics(markDuplicatesSpark.out.bam_for_variant_calling)

   bqsr(markDuplicatesSpark.out.bam_for_variant_calling,bootstrapping.out.filterSnps_out,bootstrapping.out.filterIndels_out)

   analyzeCovariates(bqsr.out.analyze_covariates)

// Variant calling
   haplotypeCallerERC(bqsr.out.recalibrated_bam)

    hc_files = haplotypeCallerERC.out.hc_erc_out.toList()
    project_id="${params.project_name}"
    interval_list=file("${params.interval_list}")    

   genomicsDBimport(hc_files,project_id,interval_list)

   genotypeGVCFs(genomicsDBimport.out.genomics_db)
   
   selectVariants(genotypeGVCFs.out.gvcfs_out,genotypeGVCFs.out.gvcfs_index)

   filterSnps(selectVariants.out.snps_ch)

   filterIndels(selectVariants.out.indels_ch)
}
