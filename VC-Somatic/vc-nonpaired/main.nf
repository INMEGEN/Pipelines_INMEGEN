#!/usr/bin/env nextflow
// Workflow    : Identificación de variantes somáticas
// Important   : Pipeline para muestras sin emparejar
// Institution : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática
// Versión     : 0.1 
// Docker image - pipelinesinmegen/pipelines_inmegen  -

nextflow.enable.dsl=2

include {  mutect2                   } from "../../modules/VC-Somatic/vc-nonpaired/mutect2.nf"
include {  calculateContamination    } from "../../modules/VC-Somatic/vc-scommon/contamination.nf"
include {  filterMutectCalls         } from "../../modules/VC-Somatic/vc-scommon/filtermutect.nf"
include {  variantQC                 } from "../../modules/VC-Somatic/vc-scommon/variantQC.nf"
include {  postfilter                } from "../../modules/VC-Somatic/vc-nonpaired/postfilter.nf"
include {  annovar                   } from "../../modules/annotation/annovar.nf"
include {  snpEff                    } from "../../modules/annotation/snpEff.nf"
include {  multiqc                   } from "../../modules/VC-Somatic/vc-scommon/multiqc.nf"

// Print some pipelines information
println "Pipelines Inmegen"
println "Flujo de trabajo: Indentificación de variantes somáticas"
println "Nota: Pipeline para muestras sin emparejar (sólo tumor) "
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Directorio de la referencia: $params.refdir"
println "Directorio de salida: $params.out"
println " "

workflow {

// Declare some parameters
   interval_list          = file("${params.interval_list}")
   panel_normales         = file("${params.panel_normales}")
   panel_normales_index   = file("${params.panel_normales_idx}")
   common_biallelic       = file("${params.common_biallelic}")
   common_biallelic_index = file("${params.common_biallelic_idx}")

// Data processing
     Channel.fromPath("${params.sample_sheet_tumor}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def tumor_id   = "${row.Tumor_ID}"
                         def tumor_bam  = file("${row.Tumor_Path}")
                 return [ tumor_id, tumor_bam]
         }.set { bam_ch } 

   mutect2(bam_ch,interval_list,panel_normales, panel_normales_index)

     unfilt=mutect2.out.unfilt_vcf.collect().flatten().collate( 3 )

   calculateContamination(bam_ch, interval_list, common_biallelic, common_biallelic_index)

      contables=calculateContamination.out.cont_tables.collect().flatten().collate( 3 )
      unfilt.join(contables).groupTuple().flatten().collate( 5 ).set{forfilter}

   filterMutectCalls(forfilter)
   
   postfilter(filterMutectCalls.out.filt_vcf)

// Variant annotation

   snpEff(postfiltervcf.out.filt_pass_vcf)

   annovar(postfiltervcf.out.filt_pass_vcf)

// Variant summary

   variantQC(filterMutectCalls.out.filt_vcf.collect()."${params.project_name}")

   multiqc(snpEff.out.snpeff_ch_txt,"${params.out}")
}
