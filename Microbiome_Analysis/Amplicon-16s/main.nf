#!/usr/bin/env nextflow
// Workflow    : Metagenómica -Amplicon 16S- 
// Institución : Instituto Nacional de Medicina Genómica (INMEGEN)
// Maintainer  : Subdirección de genómica poblacional y subdifección de bioinformática del INMEGEN
// Versión     : 0.1 
// Docker image - pipelinesinmegen/pipelines_inmegen -

nextflow.enable.dsl=2

// Processes for this workflow

include { fastqc                       } from "../../modules/Metagenomics/fastqc.nf"
include { trim_Galore                  } from "../../modules/Metagenomics/trim_galore.nf"
include { dada2                        } from "../../modules/Metagenomics/dada2.nf"
include { mothur_preprocessing         } from "../../modules/Metagenomics/mothur_preprocessing.nf"
include { mothur_alignment             } from "../../modules/Metagenomics/mothur_alignment.nf"
include { mothur_QC1                   } from "../../modules/Metagenomics/mothur_postalignment_QC1.nf"
include { mothur_QC2                   } from "../../modules/Metagenomics/mothur_postalignment_QC2.nf"
include { mothur_chimera               } from "../../modules/Metagenomics/mothur_chimera.nf"
include { mothur_classify_OTUS         } from "../../modules/Metagenomics/mothur_classify_OTUS.nf"
include { mothur_phylogenetic          } from "../../modules/Metagenomics/mothur_phylogenetic.nf"
include { mothur_plots                 } from "../../modules/Metagenomics/mothur_plots.nf"
include { krona                        } from "../../modules/Metagenomics/krona.nf"
include { multiqc                      } from "../../modules/Metagenomics/multiqc.nf"

// Some useful information
println " "
println "Pipelines INMEGEN"
println "Flujo de trabajo: Metagenómica - Amplicon 16S -"
println "Imagen de docker: pipelinesinmegen/pipelines_inmegen"
println " "
println "Nombre del proyecto: $params.project_name"
println "Información de las muestras: $params.sample_info"
println "Directorio con las referencias: $params.ref"
println "Directorio de salida: $params.out"
println " "

workflow {
   
// Declare some parameters 

   rdada2   =file("${params.rdada2}")
   rplots   =file("${params.rplots}")
   pyphylip =file("${params.pyphylip}")
   pykrona  =file("${params.pykrona}")
   data_dir = "${params.out}" + "/trimming_files/data"
   raw_data = "${params.data_dir}"

// Data preprocessing

   Channel.fromPath("${params.sample_info}" )
          .splitCsv(sep:"\t", header: true)
          .map { row ->  def sample = "${row.Sample}"
                         def R1 = file("${row.R1}")
                         def R2 = file("${row.R2}")
                 return [ sample, R1, R2 ]
               }
          .set { read_pairs_ch }

   fastqc(read_pairs_ch)

   trim_Galore(read_pairs_ch)

// Dada2 workflow

   dada2(raw_data,rdada2)

// Mothur workflow

   mothur_preprocessing(trim_Galore.out.trim_fq.collect(),data_dir,"${params.project_name}")
   
   mothur_alignment(mothur_preprocessing.out.mth_pre)

   mothur_QC1(mothur_alignment.out.mth_aling,mothur_preprocessing.out.mth_pre)

   mothur_chimera(mothur_QC1.out.mth_qc1)

   mothur_QC2(mothur_chimera.out.mth_chm)

   mothur_classify_OTUS(mothur_QC2.out.mth_qc2)

   mothur_phylogenetic(mothur_QC2.out.mth_qc2,mothur_classify_OTUS.out.mth_otu,pyphylip)

   mothur_plots(mothur_classify_OTUS.out.mth_otu, mothur_classify_OTUS.out.mth_otu_shared,mothur_phylogenetic.out.mt_phylo,rplots)

   krona(mothur_classify_OTUS.out.mth_tax_summary,pykrona)

// Amplicon 16s summary

   multiqc(krona.out.krona_ch.collect(),"${params.out}")
}
