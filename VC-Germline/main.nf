#!/usr/bin/env nextflow
// Workflow    : Identificación conjunta de variantes germinales con GATK4
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática del INMEGEN
// Versión     : 0.2
// Docker image - pipelinesinmegen/pipelines_inmegen -

nextflow.enable.dsl=2

// Processes for this workflow
include { fastqc                             } from "../modules/VC-Germline/fastqc.nf"
include { multiqc                            } from "../modules/VC-Germline/multiqc.nf"
include { trimmomatic                        } from "../modules/VC-Germline/trimmomatic.nf"
include { align                              } from "../modules/VC-Germline/bwa_germline.nf"
include { mergeSam                           } from "../modules/VC-Germline/mergesamfiles.nf"
include { markDuplicatesSpark                } from "../modules/common/markDuplicatesSpark.nf"
include { getMetrics                         } from "../modules/metricts/getmetrics.nf"
include { metricswes                         } from "../modules/metricts/metrics_wes.nf"
include { metricswgs                         } from "../modules/metricts/metrics_wgs.nf"
include { summary_wes                        } from "../modules/metricts/summary_wes.nf"
include { summary_wgs                        } from "../modules/metricts/summary_wgs.nf"
include { bqsr                               } from "../modules/VC-Germinal/bqsr_recal.nf"
include { analyzeCovariates                  } from "../modules/common/analyzecovariates.nf"
include { haplotypeCallerERC                 } from "../modules/VC-Germline/haplotypecaller_erc.nf"
include { genomicsDBimport                   } from "../modules/VC-Germline/genomicsDBimport.nf"
include { genotypeGVCFs                      } from "../modules/VC-Germline/genotypegvcfs.nf"
include { selectVariants                     } from "../modules/VC-Germline/selectvariants.nf"
include { vqsrsnps                           } from "../modules/VC-Germline/vqsr_snps.nf"
include { vqsrindels                         } from "../modules/VC-Germline/vqsr_indels.nf"
include { joinvcfs                           } from "../modules/VC-Germline/joinvcfs.nf"
include { postfiltervcf                      } from "../modules/VC-Germline/postfilter.nf"
include { splitVCFs as split_filt            } from "../modules/VC-Germline/splitvcf.nf"

// Some useful information
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Identificación conjunta de variantes germinales con GATK4"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen:public"
println " "
println "Nombre del proyecto: $params.project_name"
println "Directorio con los datos crudos: $params.reads"
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

    trimmomatic(read_pairs_ch,adapters)
 
    fastqc(trimmomatic.out.trim_fq)

// Align and mark duplicates 

    align(trimmomatic.out.trim_fq)
 
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

    metricswgs(markDuplicatesSpark.out.bam_for_variant_calling,bed_file)
    summary_wgs((metricswgs.out.summary_file.collect()),"${params.project_name}","${params.out}"+"/metrics/summary") 
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
   split_filt(postfiltervcf.out.filt_pass_vcf,"filtered")

// Variant summary

   multiqc(split_filt.out.vcf_persample.collect(),"${params.out}")
}
