#!/usr/bin/env nextflow
// Workflow    : Identificación conjunta de variantes germinales utilizando bootstrapping y GATK4
// IMPORTANTE  : Ya qué, se utiliza bootstrapping para la recalibración de la calidad de las bases este pipeline puede utilizarse para especies distintas al humano.
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática del INMEGEN
// Versión     : 0.1 
// Docker image - pipelinesinmegen/pipelines_inmegen -

nextflow.enable.dsl=2

// Processes for this workflow
include { multiqc                       } from "../modules/VC-Germline/multiqc.nf"
include { trim_Galore                   } from "../modules/bootstrapping/trim_galore.nf"
include { align                         } from "../modules/VC-Germinal/bwa_germline.nf"
include { mergeSam                      } from "../modules/VC-Germline/mergesamfiles.nf"
include { markDuplicatesSpark           } from "../modules/common/markDuplicatesSpark.nf"
include { getMetrics                    } from "../modules/metrics/getmetrics.nf"
include { metricswes                    } from "../modules/metrics/metrics_wes.nf"
include { metricswgs                    } from "../modules/metrics/metrics_wgs.nf"
include { summary_wes                   } from "../modules/metrics/summary_wes.nf"
include { summary_wgs                   } from "../modules/metrics/summary_wgs.nf"
include { bqsr                          } from "../modules/bootstrapping/bqsr.nf"
include { analyzeCovariates             } from "../modules/common/analyzecovariates.nf"
include { haplotypeCallerERC            } from "../modules/VC-Germinal/haplotypecaller_erc.nf"
include { haplotypeCaller               } from "../modules/bootstrapping/haplotypecaller.nf"
include { genomicsDBimport              } from "../modules/VC-Germinal/genomicsDBimport.nf"
include { genotypeGVCFs                 } from "../modules/VC-Germinal/genotypegvcfs.nf"
include { selectVariants                } from "../modules/VC-Germinal/selectvariants.nf"
include { filterSnps                    } from "../modules/bootstrapping/filtersnps.nf"
include { filterIndels                  } from "../modules/bootstrapping/filterindels.nf"
include { joinvcfs                      } from "../modules/VC-Germline/joinvcfs.nf"
include { variantQC                     } from "../modules/metrics/variantQC.nf"
include { postfiltervcf                 } from "../modules/VC-Germline/postfilter.nf"
include { splitVCFs as split_filt       } from "../modules/VC-Germline/splitvcf.nf"

// Some useful information
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Identificación conjunta de variantes germinales con GATK4"
println "IMPORTANTE: Este flujo utiliza bootstrapping para BQSR"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Tipo de análisis (true = WES, false = WGS): $params.wes"
println "Varios lanes por muestra (true = sí, false = no): $params.multiple_lanes"
println "Directorio de la referencia: $params.refdir"
println "Directorio de salida: $params.out"
println " "

workflow bootstrapping {

   take: data_1 
   
   main:
           
   haplotypeCaller(data_1)
   
   selectVariants(haplotypeCaller.out.hc_output)

   filterSnps(selectVariants.out.snps_ch)

   filterIndels(selectVariants.out.indels_ch)
 
   emit:
   filterSnps_out   = filterSnps.out.filtered_snps
   filterIndels_out = filterIndels.out.filtered_indels     
} 

workflow {

// Declare some parameters 

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

    trim_Galore(read_pairs_ch)

    align(trim_Galore.out.trim_fq)
 
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

// Bootstrapping for BQSR    

   bootstrapping(markDuplicatesSpark.out.bam_for_variant_calling)

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

// Concatenate snps + indels and potsfiltering   

   ch_fsnps   = filterSnps.out.filtered_snps.collect().flatten().collate( 3 )
   ch_findels = filterIndels.out.filtered_indels.collect().flatten().collate( 3 )

   ch_fsnps.join(ch_findels).groupTuple().flatten().collate( 5 ).set{join_vcfs}

   joinvcfs(join_vcfs)

   variantQC(joinvcfs.out.join_vars_filt)

   postfiltervcf(joinvcfs.out.join_vars_filt)

// Summary

   multiqc(postfiltervcf.out.filt_pass_vcf.collect(),"${params.out}")
}
