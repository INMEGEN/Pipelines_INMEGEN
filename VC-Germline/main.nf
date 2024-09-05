#!/usr/bin/env nextflow
// Workflow    : Identificación conjunta de variantes germinales con GATK4
// Institution : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática
// Versión     : 0.1
// Docker image - pipelinesinmegen/pipelines_inmegen -

nextflow.enable.dsl=2

// Processes for this workflow
include { fastqc                             } from "../modules/VC-Germline/fastqc.nf"
include { multiqc                            } from "../modules/VC-Germline/multiqc.nf"
include { fastp                              } from "../modules/VC-Germline/fastp.nf"
include { align                              } from "../modules/VC-Germline/bwa_germline.nf"
include { mergeSam                           } from "../modules/VC-Germline/mergesamfiles.nf"
include { markDuplicatesSpark                } from "../modules/common/markDuplicatesSpark.nf"
include { getMetrics                         } from "../modules/metrics/getmetrics.nf"
include { metricswes                         } from "../modules/metrics/metrics_wes.nf"
include { metricswgs                         } from "../modules/metrics/metrics_wgs.nf"
include { summary_wes                        } from "../modules/metrics/summary_wes.nf"
include { summary_wgs                        } from "../modules/metrics/summary_wgs.nf"
include { bqsr                               } from "../modules/VC-Germline/bqsr_recal.nf"
include { analyzeCovariates                  } from "../modules/common/analyzecovariates.nf"
include { haplotypeCallerERC                 } from "../modules/VC-Germline/haplotypecaller_erc.nf"
include { genomicsDBimport                   } from "../modules/VC-Germline/genomicsDBimport.nf"
include { genotypeGVCFs                      } from "../modules/VC-Germline/genotypegvcfs.nf"
include { selectVariants                     } from "../modules/VC-Germline/selectvariants.nf"
include { vqsrsnps                           } from "../modules/VC-Germline/vqsr_snps.nf"
include { vqsrindels                         } from "../modules/VC-Germline/vqsr_indels.nf"
include { joinvcfs                           } from "../modules/VC-Germline/joinvcfs.nf"
include { variantQC                          } from "../modules/metrics/variantQC.nf"
include { postfiltervcf                      } from "../modules/VC-Germline/postfilter.nf"
include { splitVCFs                          } from "../modules/VC-Germline/splitvcf.nf"
include { snpEff                             } from "../modules/annotation/snpEff.nf"

// Print some pipelines information
println "Pipelines Inmegen"
println "Flujo de trabajo: Identificación conjunta de variantes germinales"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen:public"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Tipo de análisis (true = WES, false = WGS): $params.wes"
println "Varios lanes por muestra (true = sí, false = no): $params.multiple_lanes"
println "Directorio de la referencia: $params.refdir"
println "Directorio de salida: $params.out"
println " "

workflow {

// Declare some parameters 
   adapters=file("${params.adapters}")
   bed_file=file("${params.bed_file}")
   bed_filew=file("${params.bed_filew}")
   interval_list=file("${params.interval_list}")

// Data preprocessing

   Channel.fromPath("${params.sample_info}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sample = "${row.Sample_name}"
                         def sample_id = "${row.SampleID}"
                         def PU = "${row.RG_PU}"
                         def PL = "${row.RG_PL}"
                         def LB = "${row.RG_LB}"
                         def R1 = file("${row.R1}")
                         def R2 = file("${row.R2}")
                 return [ sample, sample_id, PU, PL, LB , R1, R2 ]
               }
          .set { read_pairs_ch }

    fastp(read_pairs_ch,adapters)

    fastqc(fastp.out.trim_fq)

// Align and mark duplicates

    align(fastp.out.trim_fq)
 
    if ("${params.multiple_lanes}" == true){ 

     xa = align.out.aligned_reads_ch.collect().flatten().collate( 2 )
     xa.map { a , b -> def key = a.toString().tokenize('_').get(0)
                       return tuple("${key}" , b)
            }.groupTuple() | mergeSam

    markDuplicatesSpark(mergeSam.out.merged_sam_ch)

    } 
    else {

    markDuplicatesSpark(align.out.aligned_reads_ch)

    }

// Get metrics 
   
    if ("${params.wes}" == true){

    metricswes(markDuplicatesSpark.out.bam_for_variant_calling,bed_file,bed_filew)
    summary_wes(metricswes.out.summary_file.collect(),"${params.project_name}","${params.out}"+"/metrics/summary")

    }
    else {

    metricswgs(markDuplicatesSpark.out.bam_for_variant_calling)
    summary_wgs(metricswgs.out.summary_file.collect(),"${params.project_name}","${params.out}"+"/metrics/summary") 
    }

   getMetrics(markDuplicatesSpark.out.bam_for_variant_calling)

// Base quality recalibration

   bqsr(markDuplicatesSpark.out.bam_for_variant_calling)

   analyzeCovariates(bqsr.out.analyze_covariates)

// Variant calling

   haplotypeCallerERC(bqsr.out.recalibrated_bam)

    hc_files = haplotypeCallerERC.out.hc_erc_out.toList()

   genomicsDBimport(hc_files,"${params.project_name}",interval_list)

   genotypeGVCFs(genomicsDBimport.out.genomics_db)

   selectVariants(genotypeGVCFs.out.gvcfs_out)

   vqsrsnps(selectVariants.out.snps_ch)

   vqsrindels(selectVariants.out.indels_ch)

// Concatenate and potsfilter snps + indels 
 
   joinvcfs(vqsrsnps.out.snps_filt_ch,vqsrindels.out.indels_filt_ch)
   
   postfiltervcf(joinvcfs.out.join_vars_filt)
   splitVCFs(postfiltervcf.out.filt_pass_vcf,"filtered")

// Variant annotation

   snpEff(postfiltervcf.out.filt_pass_vcf)

// Variant summary

   variantQC(joinvcfs.out.join_vars_filt) 

   multiqc(snpEff.out.snpeff_ch_txt,"${params.out}")
}
