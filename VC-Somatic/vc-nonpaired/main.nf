#!/usr/bin/env nextflow
// Workflow    : Identificación de variantes somáticas con GATK4
// Importante  : Pipeline para muestras sin emparejar
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdirección de bioinformática del INMEGEN
// Versión     : 0.1 
// Docker image - pipelinesinmegen/pipelines_inmegen  -

nextflow.enable.dsl=2

include {  mutect2                   } from "../../modules/VC-Somatic/vc-nonpaired/mutect2.nf"
include {  calculateContamination    } from "../../modules/VC-Somatic/vc-scommon/contamination.nf"
include {  filterMutectCalls         } from "../../modules/VC-Somatic/vc-scommon/filtermutect.nf"
include {  postfilter                } from "../../modules/VC-Somatic/vc-nonpaired/postfilter.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Pipelines INMEGEN"
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

     Channel.fromPath("${params.sample_sheet_tumor}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def tumor_id   = "${row.Tumor_ID}"
                         def tumor_bam  = file("${row.Tumor_Path}")
                 return [ tumor_id, tumor_bam]
         }.set { bam_ch } 

   interval_list          = file("${params.interval_list}")
   panel_normales         = file("${params.panel_normales}")
   panel_normales_index   = file("${params.panel_normales_idx}")
   common_biallelic       = file("${params.common_biallelic}")
   common_biallelic_index = file("${params.common_biallelic_idx}")
 
   mutect2(bam_ch,interval_list,panel_normales, panel_normales_index)

     unfilt=mutect2.out.unfilt_vcf.collect().flatten().collate( 3 )

   calculateContamination(bam_ch, interval_list, common_biallelic, common_biallelic_index)

      contables=calculateContamination.out.cont_tables.collect().flatten().collate( 3 )
      unfilt.join(contables).groupTuple().flatten().collate( 5 ).set{forfilter}

   filterMutectCalls(forfilter)

   postfilter(filterMutectCalls.out.filt_vcf) 
}
