# Flujo de trabajo identificación de variantes somáticas utilizando NextFlow y GATK

Este pipeline realiza la identificación de variantes somáticas a partir de archivos de secuenciación masiva (WGS/WES) y se divide en 3 subflujos de trabajos correspondientes a una configuración de análisis en específico:

- Modo tumor-only [vc-nonpaired]
- Modo pareado [vc-paired]
- Panel de normales [panelofnormals]

**Nota:** Por el momento el análisis sólo está disponible para datos de lectura corta en humano (ilummina paired-end).
Para conocer más sobre la indentificación de variantes somáticas con GATK4 (Mutect2) consulta la siguiente [liga](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2).

   ### Modo tumor-only
Este flujo de trabajo identifica las variantes somáticas utilizando únicamente el panel de normales para distinguir las alteraciones no somáticas.

   ### Modo pareado
Al momento de identificar las variantes somáticas este flujo de trabajo utiliza una muestra normal correspondiente al mismo paciente pero que no sea de tejido tumoral. Que junto al panel de normales aumenta la precisión de la identificación de variantes somáticas.

### Panel de normales
Para distinguir las variantes que derivan de las muestras normales (tejido sano que no tiene alteraciones somáticas) de las muestras tumorales se utiliza esta herramienta. Para general el VCF denominado panel de normales es necesario tener al menos 40 muestras normales procesadas de la misma manera que las muestras a procesar (tumor).
Para mayor información ver el siguiente [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-).

**Importante**
En el caso de que no contar con 40 muestras normales secuenciadas con las mismas condiciones para generar el panel de normales, se utilizará el que proporciona GATK de 1000 genomas. 

## Solicitud de servicio

Para solicitar este flujo de trabajo como servicio debes de entregar al personal de INMEGEN: 

- Archivos de secuenciación FASTQ (Illumina paired-end).
- Archivo con la información experimental (los identificadores de las muestras indicando si son normales o de tumor, identificador de la librería, plataforma de secuenciación,  número de lane).
- En caso de WES específicar el kit utilizado.

## Instrucciones de uso 

### Preparar el ambiente de trabajo

1. Te debes asegurar de contar con las siguientes herramientas informaticas:
	- [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7)
	- [Docker](https://docs.docker.com/) (23.0.5)
	- Imagen de docker pipelinesinmegen/pipelines_inmegen:public, la puedes descargar con el comando: 

                docker pull pipelinesinmegen/pipelines_inmegen:public

3. Asegurarse de contar con los siguientes archivos, necesarios para el pipeline:
	- Genoma hg38
	- Índice del genoma de referencia (generado con SAMTOOLS faidx)
	- Índice de [BWA](https://bio-bwa.sourceforge.net/bwa.shtml)
	- Archivos de recalibración de BQSR y VQSR
	- Archivo [gnomAD VCF](https://gnomad.broadinstitute.org/downloads/)

**Todos estos archivos se pueden descargar del** [bundle de GATK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).
** Se recomienda que todos estos archivos se encuentreb en el mismo directorio.

**IMPORTANTE**
Estos flujos de trabajo utilizan archivos bam previamente procesados con el flujo [Data-preprocessing](https://github.com/INMEGEN/Pipelines_INMEGEN/tree/Principal/Data_preprocessing)

### Ejecutar el flujo de trabajo

Para correr este pipeline se deben de ejecutar las siguientes instrucciones:

 1. Después generar el archivo sample_*.tsv con la información que se describe en la sección - Formato del archivo con la información de las muestras -
 2. Editar el archivo de nextflow.config con la siguiente información:

	- Ruta de los archivos *BAM*
	- Ruta del directorio de salida de nextflow
	- Nombre del proyecto 
	- Ruta y nombre del genoma de referencia
	- Ruta del panel de normales (Modo paired y unpired)
	- Ruta del o los archivos sample_info.tsv (normales y tumor)
	- Ruta del archivo gnomAD VCF
	- Nombre del índice de BWA
	- Ruta del archivo con la lista de intervalos, en el caso de WES es el archivo BED del kit, para más información consulta la siguiente [liga](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists)
	- Parámetros de mutect2 (parámetros: pading, germline_sites, pon_sites)
	- Número de núcleos por proceso (parámetro runOptions) 
        - Número de procesos que se ejecutarán de forma simultánea (parámetro queueSize)

  3. Ejecutar el comando correspondiente a cada subflujo de trabajo: 

                bash run_nextflow.sh /path/to/out/dir

### Formato del archivo con la información experimental

##### Panel de normales 

Para tener un buen control de los archivos a procesar (formato bam), el archivo sample_*.tsv debe de incluir la siguiente información por columna:
 
		Sample	Path	

 - Sample   = Nombre de la muestras normales, se recomienda el formato [nombre de la muestra - número de muestra]
 - Path     = Ruta absoluta del archivo bam de la muestra normal

**Nota:** Recuerda que el archivo debe estar separado por tabulador (\t).

##### Modo paired
Para tener un buen control de los archivos a procesar (formato bam), el archivo sample_*.tsv debe de incluir la siguiente información por columna:
 
                Tumor_id	Tumor_Path	Normal_id	Normal_Path     

 - Tumor_id    = Nombre de la muestras tumor, se recomienda el formato [nombre de la muestra - número de muestra]
 - Tumor_Path  = Ruta absoluta del archivo bam de la muestra Tumor_id
 - Normal_id   = Nombre de la muestras normales, se recomienda el formato [nombre de la muestra - número de muestra]
 - Normal_Path = Ruta absoluta del archivo bam de la muestra Normal_id

**Nota:** Los identificadores por renglón deben de pertenecer a la misma muestra (paciente)
**Nota:** Recuerda que el archivo debe estar separado por tabulador (\t).
       
##### Modo unpaired       

Para tener un buen control de los archivos a procesar (formato bam), el archivo sample_*.tsv debe de incluir la siguiente información por columna:
 
                Tumor_id        Tumor_Path

 - Tumor_id    = Nombre de la muestras tumor, se recomienda el formato [nombre de la muestra - número de muestra]
 - Tumor_Path  = Ruta absoluta del archivo bam de la muestra Tumor_id

**Nota:** Recuerda que el archivo debe estar separado por tabulador (\t).

#### Las herramientas utilizadas para correr este flujo de trabajo son:

 - GATK (4.2.6.1)
 - R (4.2.3)
 - Picard Tools (2.27.5)

## Diagrama de flujo del pipeline 

Para una mayor descripción de la información del pipeline ejecutado se anexa el siguiente diagrama de flujo basado en [las buenas prácticas de GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-).

![Flujo identificación de variantes somaticas](../flowcharts/flujo_VCS.PNG)
