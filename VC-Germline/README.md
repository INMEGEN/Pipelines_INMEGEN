# Flujo de trabajo identificación conjunta de variantes germinales a partir de datos WGS/WES 

Este pipeline realiza la identificación conjunta de variantes germinales a partir de archivos de secuenciación masiva (WGS/WES).

**Nota:** 
 - Por el momento el análisis sólo está disponible para datos ilummina paired-end en humano.  
 - Si se desea otra especie revisar el flujo identificación conjunta de variantes germinales a partir de datos WGS/WES con bootstrapping. Este flujo se provee como parte de este repositorio pero no ha sido probado por personal del INMEGEN en otras especies.


## Solicitando este servicio

Para solicitar este flujo de trabajo como servicio debes de entregar al personal de INMEGEN 

- Archivos de secuenciación FASTQ (Illumina paired-end).
- Archivo con la información experimental (identificador de la muestra, identificador de la librería, plataforma de secuenciación,  número de lane).
- En caso de WES específicar el kit utilizado.


## Utilizando este flujo por tu cuenta: Instrucciones de uso 

Los archivos que necesitas se describen en el apartando **"Solicitando este servicio"**.

### Preparación de ambiente de trabajo

1. Te debes asegurar que contar con:
  - [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7),
  - [Docker](https://docs.docker.com/) (23.0.5) y
  - la imagen de docker pipelinesinmegen/pipelines_inmegen:public. 
2. Descargar el genoma y otros archivos necesarios:
  - Descargar el genoma hg38
  - Obtener su índice con SAMTOOLS faidx
  - Obtener su índice con [BWA](https://bio-bwa.sourceforge.net/bwa.shtml)
  - Descargar los archivos de recalibración de BQSR y VQSR del [bundle de GATK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).

### Ejecución del flujo de trabajo

Asegúrate de estar parado en el directorio de este flujo de trabajo.

 1. Seleccionar una ruta y el nombre para el directorio de salida [DANIEL ¿?¿? no entendi como que tengo que hacer?, donde la defino o como la selecciono?]
 2. Después generar el archivo sample_info.tsv con la información que se describe en la sección - Formato del archivo con la información de las muestras -
 3. Editar el archivo de nextflow.config con la siguiente información:

        - Ruta de los archivos *fastq*
        - Ruta del directorio de salida de nextflow
        - Nombre del proyecto 
        - Ruta del índice de BWA
        - Ruta del archivo sample_info.tsv
        - Nombre del índice de BWA
        - Ruta del archivo con la lista de intervalos [DANIEL ¿?¿? esto es solo para WES y es el BED del kit, no? hay que especificarlo]
        - Condiciones del análisis (número de núcleos a utilizar por proceso, número de procesos simultáneos e información adicional)[DANIEL ¿?¿? seguro esta parte tiene un formato específico, hay que especificarlo]

  4. Ejecutar el comando: 

                bash run_nextflow.sh /path/to/out/dir

#### Formato del archivo con la información de las muestras

En el archivo sample_info.tsv incluir la siguiente información por columna:
 
		Sample_ID	Sample_name	RG	PU	R1	R2

 - Sample_ID   = Nombre completo de los archivos e identificador de las muestras (debe poseer el nombre o identificador, el número de muestra y número de lane) [DANIEL ¿?¿? seguro esta parte tiene un formato específico, hay que especificarlo]
 - Sample_name = Nombre o identificador de la muestra (comprende el identificador y el número de muestra) [DANIEL ¿?¿? que es el numero de la muestra]
 - RG          = ReadGroup name - El nombre del grupo de lectura de las muestras a procesar [DANIEL ¿?¿? esto tambien tiene un formato, hay que especificarlo, y decir que esperas que venga aquí]
 - PU          = Plataforma + Libreria + número de lane
 - R1          = Ruta absoluta del archivo fastq R1
 - R2          = Ruta absoluta del archivo fastq R2

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

 Se anexa el siguiente diagrama de flujo con la descripción completa del pipeline ejecutado. Este pipeline está basado en las buenas prácticas de GATK.

[Consulta el flujo de identificación de variantes germinal de GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

![Flujo identificación de variantes germinal](../flowcharts/flujo_VCG.PNG)
