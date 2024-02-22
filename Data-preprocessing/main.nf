  #!/usr/bin/env nextflow
// Workflow:    Preprocesamiento de datos para identificación de variantes con GATK4  
// Institución: Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer:  Subdirección de genómica poblacional y subdirección de bioinformática (INMEGEN)
// Versión:     0.1 
// Docker image - pipelinesinmegen/pipelines_inmegen -

nextflow.enable.dsl=2

include { fastqc                        } from "../modules/data_preprocessing/fastqc.nf"
include { trimmomatic                   } from "../modules/data_preprocessing/trimmomatic.nf"
include { fastq_to_sam                  } from "../modules/data_preprocessing/fastq_to_sam.nf"
include { mark_duplicates               } from "../modules/data_preprocessing/mark_duplicates.nf"
include { sam_to_fastq                  } from "../modules/data_preprocessing/sam_to_fastq.nf"
include { align                         } from "../modules/data_preprocessing/align.nf"
include { merge_bam_alignment           } from "../modules/data_preprocessing/merge_bam_alignment.nf"
include { mergeSam                      } from "../modules/data_preprocessing/mergesamfiles.nf"
include { markDuplicatesSpark           } from "../modules/common/markDuplicatesSpark.nf"
include { getMetrics                    } from "../modules/metrics/getmetrics.nf"
include { metricswes                    } from "../modules/metrics/metrics_wes.nf"
include { metricswgs                    } from "../modules/metrics/metrics_wgs.nf"
include { summary_wes                   } from "../modules/metrics/summary_wes.nf"
include { summary_wgs                   } from "../modules/metrics/summary_wgs.nf"
include { bqsr                          } from "../modules/data_preprocessing/bqsr_recal.nf"
include { analyzeCovariates             } from "../modules/common/analyzecovariates.nf"
include { multiqc                       } from "../modules/data_preprocessing/multiqc.nf"

// Imprimir la ruta de algunos directorios importantes
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Preprocesamiento de datos para GATK4"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
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
                 return [ sample, sample_id, PU, PL, LB, R1, R2 ]
               }
          .set { read_pairs_ch}

   trimmomatic(read_pairs_ch,adapters)

   fastqc(trimmomatic.out.trim_fq)    

   fastq_to_sam(trimmomatic.out.trim_fq)

       ch_fqtsam=fastq_to_sam.out.fastq_to_sam_ch.collect().flatten().collate( 2 )

   mark_duplicates(fastq_to_sam.out.fastq_to_sam_ch)

   sam_to_fastq(mark_duplicates.out.mark_duplicates_bam)

   align(sam_to_fastq.out.sam_to_fastq_ch)

       ch_aln=align.out.aligned_reads.collect().flatten().collate( 2 )
       ch_aln.join(ch_fqtsam).groupTuple().flatten().collate( 3 ).set{joinfmerge}

   merge_bam_alignment(joinfmerge)

   if ("${params.multiple_lanes}" == true){ 

    xa = merge_bam_alignment.out.bam_alignment_ch.collect().flatten().collate( 2 )
    xa.map { a , b -> def key = a.toString().tokenize('_').get(0)
                      return tuple("${key}" , b)
           }.groupTuple() | mergeSam

    markDuplicatesSpark(mergeSam.out.merged_sam_ch)

   } 
   else {

    markDuplicatesSpark(merge_bam_alignment.out.bam_alignment_ch)

   }

   getMetrics(markDuplicatesSpark.out.bam_for_variant_calling)

   if ("${params.wes}" == true){

    metricswes(markDuplicatesSpark.out.bam_for_variant_calling,bed_file,bed_filew)
    summary_wes(metricswes.out.summary_file.collect(),"${params.project_name}","${params.out}"+"/metrics/summary")

   }
   else {

    metricswgs(markDuplicatesSpark.out.bam_for_variant_calling,bed_file)
    summary_wgs(metricswgs.out.summary_file.collect(),"${params.project_name}","${params.out}"+"/metrics/summary") 
   }

    bqsr(markDuplicatesSpark.out.bam_for_variant_calling)

    analyzeCovariates(bqsr.out.analyze_covariates)

    multiqc(analyzeCovariates.out.analyzed_covariates_ch.collect(),"${params.out}")
}
