# Flujo de trabajo: Pre-procesamiento de datos para el descubrimiento de variantes con GATK

Este flujo de trabajo realiza el preprocesamiento de archivos de secuenciación masiva (WGS/WES) en formato *fastq* para realizar una indentificación de variantes con GATK. 
Con este flujo de trabajo de obtiene un archivo bam limpio como lo indica la siguiente [liga](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently) de GATK, lo que es parte del pre-procesamiento de datos para el descubrimiento de variantes cuyo tutorial de gatk se encuentra [aquí](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery) 

**Nota:** Por el momento el análisis sólo está disponible para datos de lectura corta en humano (ilummina paired-end).

## Solicitud de servicio

Para solicitar este flujo de trabajo como servicio debes de entregar al personal de INMEGEN: 

- Archivos de secuenciación FASTQ (Illumina paired-end).
- Archivo con la información experimental (identificador de la muestra, identificador de la librería, plataforma de secuenciación,  número de lane).
- En caso de WES específicar el kit utilizado.

## Implementando este flujo por tu cuenta: Instrucciones de uso 

Los archivos que necesitas se describen en el apartando **"Solicitud de servicio"**.

### Preparar el ambiente de trabajo

Si deseas utilizar este flujo de trabajo sin apoyo del personal del INMEGEN sigue las siguientes instrucciones:

### Preparar el ambiente de trabajo

1. Te debes asegurar de contar con las siguientes herramientas informaticas:
        - [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7)
        - [Docker](https://docs.docker.com/) (23.0.5)
        - Imagen de docker pipelinesinmegen/pipelines_inmegen:public, la puedes descargar con el comando: 

                docker pull pipelinesinmegen/pipelines_inmegen:public

 
2. Asegurarse de contar con los siguientes archivos, necesarios para el pipeline:
        - Genoma hg38
        - Índice del genoma de referencia (generado con SAMTOOLS faidx)
        - Índice de [BWA](https://bio-bwa.sourceforge.net/bwa.shtml)
        - Archivos de recalibración de BQSR y VQSR

**Todos estos archivos se pueden descargar del** [bundle de GATK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&for
ceOnObjectsSortingFiltering=false).
** Se recomienda que todos estos archivos se encuentreb en el mismo directorio.

### Ejecutar el flujo de trabajo

Para correr este pipeline se deben de ejecutar las siguientes instrucciones:

 1. Generar el archivo sample_info.tsv con la información que se describe en la sección - Formato del archivo con la información de las muestras -
 2. Editar el archivo de nextflow.config con la siguiente información:

        - Ruta de los archivos *fastq*
        - Ruta del directorio de salida de nextflow
        - Nombre del proyecto 
        - Ruta y nombre del genoma de referencia
        - Ruta del índice de BWA
        - Ruta del archivo sample_info.tsv
        - Número de muestras (parámetro batchsize)
        - Ruta del archivo con la lista de intervalos, en el caso de WES es el archivo BED del kit, para más información consulta la siguiente [liga](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists)
        - Número de núcleos por proceso (parámetro runOptions) 
        - Número de procesos que se ejecutarán de forma simultánea (parámetro queueSize)

Para opciones de configuración especificas para tu servidor o cluster puedes consultar la siguiente [liga](https://www.nextflow.io/docs/latest/config.html) 

  4. Ejecutar el comando: 

                bash run_nextflow.sh /path/to/out/dir

### Formato del archivo con la información de las muestras

En el archivo sample_info.tsv incluir la siguiente información por columna:
 
			Sample	RG	PU	R1	R1
Dónde: 
 - Sample = Nombre completo de los archivos, se recomienda el formato [identificador único-número de muestra-número de lane]
 - RG     = Nombre del grupo de lectura de la muestra, revisar la siguiente [liga](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) para más información
 - PU     = Plataforma + Libreria + número de lane
 - R1     = Ruta absoluta del archivo fastq R1
 - R2     = Ruta absoluta del archivo fastq R2
			
**Nota:** Recuerda que el archivo debe estar separado por tabulador (\t).

#### Este flujo de trabajo requiere de las siguientes herramientas:

 - FastQC (0.11.9)
 - MultiQC (1.11)
 - Trim Galore (0.6.7)
 - GATK (4.2.6.1)
 - R (4.2.3)
 - BWA (0.7.17)
 - Picard Tools (2.27.5)
 - Samtools (1.12)
 - Bcftools (1.12)

#### Diagrama con el flujo de trabajo 

![Flujo data_preprocessing](../flowcharts/flujo_DP.PNG)
