#!/usr/bin/env nextflow
// Workflow     : Identificación de variantes de datos RNA-seq
// Institution  : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer   : Subdirección de genómica poblacional y subdirección de bioinformática
// Versión      : 0.1
// Docker image : - pipelinesinmegen/pipelines_inmegen -

nextflow.enable.dsl=2

// Processes for this pipeline
include { fastqc                             } from "../modules/VC-RNAseq/fastqc.nf"
include { fastp                              } from "../modules/VC-RNAseq/fastp.nf"
include { star                               } from "../modules/VC-RNAseq/star.nf"
include { mergeSam                           } from "../modules/VC-RNAseq/mergesamfiles.nf"
include { markDuplicatesSpark                } from "../modules/common/markDuplicatesSpark.nf"
include { getMetrics                         } from "../modules/VC-RNAseq/getmetrics.nf"
include { splitNCigarReads                   } from "../modules/VC-RNAseq/splitNCigarReads.nf"
include { bqsr                               } from "../modules/VC-RNAseq/bqsr_recal.nf"
include { analyzeCovariates                  } from "../modules/common/analyzecovariates.nf"
include { haplotypeCaller                    } from "../modules/VC-RNAseq/haplotypecaller.nf"
include { selectVariants                     } from "../modules/VC-RNAseq/selectvariants.nf"
include { filterSnps                         } from "../modules/VC-RNAseq/filtersnps.nf"
include { filterIndels                       } from "../modules/VC-RNAseq/filterindels.nf"
include { joinvcfs                           } from "../modules/VC-RNAseq/joinvcfs.nf"
include { variantQC                          } from "../modules/VC-RNAseq/variantQC.nf"
include { postfiltervcf                      } from "../modules/VC-RNAseq/postfilter.nf"
include { snpEff                             } from "../modules/annotation/snpEff.nf"
include { multiqc                            } from "../modules/VC-RNAseq/multiqc.nf"

// Print some parameters
println "Pipelines Inmegen"
println "Flujo de trabajo: Indentificación de variantes con datos de RNA-seq"
println "Contenedor: Docker pipelinesinmegen/pipelines_inmegen "
println " "
println "Nombre del proyecto: $params.project_name"
println "Archivo con la información de las muestras: $params.sample_info"
println "Datos con varios lanes por muestra (true = sí, false = no): $params.multiple_samples"
println "Directorio con la referencia de star: $params.refdir_star"
println "Directorio de salida: $params.out"
println " "

workflow {
  
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
                 return [ sample, sample_id, PU, PL, LB, R1, R2 ]
               }
          .set { read_pairs_ch }

   fastp(read_pairs_ch)

   fastqc(fastp.out.trim_fq)

   star(fastp.out.trim_fq)

    if ("${params.multiple_samples}" == true){ 
     xa = star.out.aligned_reads_ch.collect().flatten().collate( 2 )
     xa.map { a , b -> def key  = a.toString().tokenize('_').get(0)
                       return tuple("${key}", b)
            }.groupTuple() | mergeSam
    markDuplicatesSpark(mergeSam.out.merged_sam_ch)
    } 
    else {
    markDuplicatesSpark(star.out.aligned_reads_ch)
    }

   getMetrics(markDuplicatesSpark.out.bam_for_variant_calling)

   splitNCigarReads(markDuplicatesSpark.out.bam_for_variant_calling)

   bqsr(splitNCigarReads.out.split_bam)

   analyzeCovariates(bqsr.out.analyze_covariates)

// Variant calling
  
   haplotypeCaller(bqsr.out.recalibrated_bam)
  
// Filtering vcfs by type

   selectVariants(haplotypeCaller.out.hc_output)

   filterSnps(selectVariants.out.snps_ch)

   filterIndels(selectVariants.out.indels_ch)

// Concatenate snps + indels and potsfiltering   

   ch_fsnps   = filterSnps.out.filtered_snps.collect().flatten().collate( 3 )
   ch_findels = filterIndels.out.filtered_indels.collect().flatten().collate( 3 )
   
   ch_fsnps.join(ch_findels).groupTuple().flatten().collate( 5 ).set{join_vcfs}

   joinvcfs(join_vcfs)

   variantQC(joinvcfs.out.join_vars_filt.collect(),"${params.project_name}")

   postfiltervcf(joinvcfs.out.join_vars_filt)

// Variant annotation

   snpEff(postfiltervcf.out.filt_pass_vcf)

// Summary

   multiqc(postfiltervcf.out.filt_pass_vcf.collect(),"${params.out}")
}
