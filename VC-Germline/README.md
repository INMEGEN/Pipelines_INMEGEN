# Flujo de trabajo identificación conjunta de variantes germinales a partir de datos WGS/WES 

Este pipeline realiza la identificación conjunta de variantes germinales a partir de archivos de secuenciación masiva (WGS/WES).

**Nota:** 
 - Por el momento el análisis sólo está disponible para datos ilummina paired-end en humano.  
 - Si se desea otra especie revisar el flujo identificación conjunta de variantes germinales a partir de datos WGS/WES con bootstrapping. Este flujo se provee como parte de este repositorio pero no ha sido probado por personal del INMEGEN en otras especies.


## Solicitud de servicio

Para solicitar este flujo de trabajo como servicio debes de entregar al personal de INMEGEN: 

- Archivos de secuenciación FASTQ (Illumina paired-end).
- Archivo con la información experimental (identificador de la muestra, identificador de la librería, plataforma de secuenciación,  número de lane).
- En caso de WES específicar el kit utilizado.


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
	- Índice de [BWA](https://bio-bwa.sourceforge.net/bwa.shtml)
	- Archivos de recalibración de BQSR y VQSR

**Todos estos archivos se pueden descargar del** [bundle de GATK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).
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
 
		Sample_ID	Sample_name	RG	PU	R1	R2

 - Sample_ID   = Nombre completo de los archivos, se recomienda el formato [identificador único-número de muestra-número de lane]
 - Sample_name = Nombre de la muestras, se recomienda el formato [nombre de la muestra - número de muestra]
 - RG          = Nombre del grupo de lectura de la muestra, revisar la siguiente [liga](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) para más información 
 - PU          = Plataforma + Libreria + número de lane
 - R1          = Ruta absoluta del archivo fastq R1
 - R2          = Ruta absoluta del archivo fastq R2

**Nota:** Recuerda que el archivo debe estar separado por tabulador (\t).

### Las herramientas utilizadas para correr este flujo de trabajo son:

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

 Se anexa el siguiente diagrama de flujo con la descripción completa del pipeline ejecutado. Este pipeline está basado en las buenas prácticas de GATK: [identificación de variantes germinal de GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

![Flujo identificación de variantes germinal](../flowcharts/flujo_VCG.PNG)
