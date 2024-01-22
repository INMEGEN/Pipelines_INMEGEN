# Flujo de trabajo llamado de variantes de datos RNAseq utilizando NextFlow y GATK4.

Este pipeline realiza la identificación de variantes a partir de archivos de secuenciación masiva (RNA-seq).

**Nota:** Por el momento el análisis sólo está disponible para datos de lectura corta de humano (Illumina paired-end).

**Nota:** GATK no soporta  la identificación conjunta de variantes en datos de RNA-seq, consulta el siguiente [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) para más información.

## Solicitud de servicio

Para solicitar este flujo de trabajo como servicio debes de entregar al personal de INMEGEN: 

- Archivos de secuenciación FASTQ (Illumina paired-end).
- Archivo con la información experimental (identificador de la muestra, identificador de la librería, plataforma de secuenciación,  número de lane).


## Implementando este flujo por tu cuenta: Instrucciones de uso 

Los archivos que necesitas se describen en el apartando **"Solicitud de servicio"**.

### Preparar el ambiente de trabajo

1. Te debes asegurar de contar con las siguientes herramientas informaticas:
	- [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7)
	- [Docker](https://docs.docker.com/) (23.0.5)
	- Imagen de docker pipelinesinmegen/pipelines_inmegen:public, la puedes descargar con el comando: 

                docker pull pipelinesinmegen/pipelines_inmegen:public

 
2. Asegurarse de contar con los siguientes archivos, necesarios para el pipeline:
	- Genoma hg38
	- Índice del genoma de referencia (generado con SAMTOOLS faidx)
	- Índice de [STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
	- Archivos de recalibración de BQSR y VQSR **se pueden descargar del** [bundle de GATK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).

**Nota:** El directorio bin contiene un bash script para generar estos acrhivos y sólo resta descargar los archivos para BQSR y VQSR. 
**Se recomienda que todos estos archivos se encuentren en el mismo directorio.**

### Ejecutar el flujo de trabajo

Para correr este pipeline se deben clonar este repositorio y ejecutar las siguientes instrucciones:

 1. Completar el archivo sample_info.tsv con la información que se describe en la sección - Formato del archivo con la información de las muestras -
 2. Editar el archivo de nextflow.config con la siguiente información:

	- Ruta de los archivos *fastq* (params.reads)
	- Ruta del directorio de salida de nextflow (params.outdir)
	- Nombre del proyecto (params.project_name)
 	- Si son multiples lanes por muestra (params.multiple_samples)
	- Ruta (params.refdir) y nombre (params.refname) del índice de STAR del genoma de referencia
 	- Ruta de los archivos para la recalibración de las bases (params.ref_dir)
	- Ruta del archivo sample_info.tsv (params.sample_sheet)
 	- Tenología de secuenciación (params.pl) ej. ILLUMINA, SOLID, LS454, HELICOS y PACBIO
	- Tipo de secuencidor  (params.pl) ej. Nextseq, Hiseq, etc,
 	- Número de núcleos por proceso (parámetro runOptions) 
	- Número de procesos que se ejecutarán de forma simultánea (parámetro queueSize)

Para opciones de configuración especificas para tu servidor o cluster puedes consultar la siguiente [liga](https://www.nextflow.io/docs/latest/config.html) 

  3. Ejecutar el comando: 

		bash run_nextflow.sh /path/to/out/dir


#### Formato del archivo con la información experimental

Para tener un buen control de los archivos a procesar (formato fastq pareados {Read_1,Read_2}), en el archivo sample_info.tsv incluir la siguiente información por columna:

 - Sample_ID   = Nombre asignadp a los archivos fastq (no tiene que ser identico al nombre del archivo), se debe utilizar el formato [identificador - número de muestra - número de lane]
 - Sample_name = Nombre de la muestras, se recomienda el formato [identificados - número de muestra]
 - RG          = Nombre del grupo de lectura de la muestra, revisar la siguiente [liga](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) para más información
 - PU          = Plataforma + Libreria + número de lane
 - R1          = Ruta absoluta del archivo fastq R1
 - R2          = Ruta absoluta del archivo fastq R2

 Ejemplo:
 
		Sample_ID	Sample_name	RG	PU	R1	R2
  		ID_S001_L001	ID_S001		Read_gorup_id	FLOWCELL_BARCODE	Path/to/fastq_R1.fq	Path/to/fastq_R2.fq

En caso de que el nombre del archivo sólo sea ID + lane (ID_L001), agregar S1 a la columna SampleID como se muestra a continuación 

		Sample_ID	Sample_name	RG	PU	R1	R2
  		ID_S1_L001	ID_S1		Read_gorup_id	FLOWCELL_BARCODE	Path/to/fastq_R1.fq	Path/to/fastq_R2.fq

**Nota:** Recuerda que el archivo debe estar separado por tabulador (\t).

#### Las herramientas utilizadas para correr este flujo de trabajo son:

 - FastQC (0.11.9)
 - MultiQC (1.11)
 - Trim Galore (0.6.7)
 - GATK (4.2.6.1)
 - R (4.2.3)
 - BWA (0.7.17)
 - Picard Tools (2.27.5)
 - Samtools (1.12)
 - Bcftools (1.12)

## Diagrama de flujo del pipeline 

Para una mayor descripción de la información del pipeline [consulta el las buenas pŕacticas de gatk para la identificación de variantes de datos de RNAseq](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-). A continuación se anexa un diagrama de flujo del pipeline ejecutado. 


![Flujo identificación de variantes rnaseq](../flowcharts/flow_VCRNA.PNG)
