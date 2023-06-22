# Flujo de trabajo identificación de variantes somáticas utilizando NextFlow y GATK

Este pipeline realiza la identificación de variantes a partir de archivos de secuenciación masiva (WGS/WES).
En caso de trabajar con el genoma hg38, los archivos como el índice de [BWA](http://bio-bwa.sourceforge.net/) y los archivos de recalibración de BQSR y VQSR se pueden descargar del [bundle de GATK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).  
Lo único que necesitas es tus archivos de lectura fastq, en caso de WES el kit utilizado y los indentificadores de las muestras indicando si son normales o de tumor.
Además, si son múltiples lanes por muestra es necesario especificar a que muestras están asociadas (revisar la información solicitada por el archivo sample_info.tsv).

Este flujo de trabajo se divide en 3 subflujos de trabajos correspondientes a una configuración de análisis en especifico

	- Panel de normales [PON]
	- Modo pareado [vc-paired]
	- Modo tumor-only [vc-nonpaired]

**Nota:** En el caso de que no contar con 40 muestras normales secuenciadas en las mismas condiciones para generar el panel de normales se utilizará el que proporciona GATK. 
**Nota:** Por el momento el análisis sólo está disponible para datos de lectura corta (ilummina paired-end).

## Instrucciones de uso 

Primero se debe asegurar que se cuenta con [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7), [Docker](https://docs.docker.com/) (23.0.5) y la imagen de docker pipe
linesinmegen/pipelines_inmegen:latest.

Estos flujos de trabajo utilizan archivos bam previamente procesados con el pipeline de preprocesamiento de datos [Data_preprocessing].

 1. Seleccionar una ruta y el nombre para el directorio de salida
 2. Después generar el archivo sample_*.tsv con la información que se describe en la sección - Formato del archivo con la información de las muestras -
 3. Editar el archivo de nextflow.config con la siguiente información:

	- Ruta de los archivos *fastq*
	- Ruta del directorio de salida de nextflow
	- Nombre del proyecto 
	- Ruta del índice de BWA
	- Ruta del archivo sample_*.tsv
	- Nombre del índice de BWA
	- Ruta del archivo con la lista de intervalos
	- Ruta del directorio de annovar
	- Parámetros de mutect2
	- Condiciones del análisis (número de núcleos a utilizar por proceso, número de procesos simultáneos e información adicional)

  4. Ejecutar el comando correspondiente a cada subflujo de trabajo: 

                bash run_nextflow.sh /path/to/out/dir

### Formato del archivo con la información experimental

	##### Panel de normales
 
Para tener un buen control de los archivos a procesar (formato bam), el archivo sample_*.tsv debe de incluir la siguiente información por columna:
 
		Sample	Path	

 - Sample   = Nombre completo de los archivos e identificador de las muestras normales (debe de contener el nombre o identificador de la muestra)
 - Path     = Ruta absoluta del archivo bam de la muestra sample

	##### Modo paired

Para tener un buen control de los archivos a procesar (formato bam), el archivo sample_*.tsv debe de incluir la siguiente información por columna:
 
                Tumor_id	Tumor_Path	Normal_id	Normal_Path     

 - Tumor_id    = Nombre de los archivos e identificador de las muestras de tumor (debe de contener el nombre o identificador de la muestra)
 - Tumor_Path  = Ruta absoluta del archivo bam de la muestra Tumor_id
 - Normal_id   = Nombre de los archivos e identificador de las muestras normales (debe de contener el nombre o identificador de la muestra)
 - Normal_Path = Ruta absoluta del archivo bam de la muestra Normal_id

**Nota:** Los identificadores por renglón deben de pertenecer a la misma muestra (paciente)

	##### Modo paired
       
Para tener un buen control de los archivos a procesar (formato bam), el archivo sample_*.tsv debe de incluir la siguiente información por columna:
 
                Tumor_id        Tumor_Path

 - Tumor_id    = Nombre de los archivos e identificador de las muestras de tumor (debe de contener el nombre o identificador de la muestra)
 - Tumor_Path  = Ruta absoluta del archivo bam de la muestra Tumor_id



#### Las herramientas utilizadas para correr este flujo de trabajo son:
>
> - FastQC (0.11.9)
> - MultiQC (1.11)
> - Openjdk (11.0.13 o superior)
> - GATK (4.2.6.1)
> - BWA (0.7.17-r 1188)
> - Picard Tools (2.0.1)
> - Samtools (1.6.0)
> - Annovar
> - bcftools (1.14)
>

## Diagrama de flujo del pipeline 

Para una mayor descripción de la información del pipeline ejecutado se anexa el siguiente diagrama de flujo basado en las buenas prácticas de GATK.
[Consulta el flujo de identificación de variantes somáticas de GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-)

![Flujo identificación de variantes somaticas](../flowcharts/flujo_VCS.PNG)
