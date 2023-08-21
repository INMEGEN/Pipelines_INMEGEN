# Flujo de trabajo llamado de variantes de datos RNAseq utilizando NextFlow y GATK4.

Este pipeline realiza la identificación de variantes a partir de archivos de secuenciación masiva (RNA-seq).

**Nota:** Por el momento el análisis sólo está disponible para datos de lectura corta (Illumina paired-end).
**Nota:** GATK no soporta  la identificación conjunta de variantes en datos de RNA-seq, ver el siguiente [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) para más información.


### Para solicitar este flujo de trabajo como servicio debes de entregar al personal de INMEGEN 

- Archivos de lectura fastq (Illumina paired-end).
- Archivo con la información experimental (identificador de la muestra, plataforma y librería de secuenciación, si son múltiples lanes especificar el número).

## Instrucciones de uso 

Si deseas utilizar este flujo de trabajo sin apoyo del personal del INMEGEN sigue las siguientes instrucciones.

En caso de trabajar con el genoma humano de referencia hg38, la forma de crear el índice de STAR se encuentra en el directorio bin. Los archivos adicionales (bases de datos de snps e indeles conocidos) se pueden obtener del [bundle de GATK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).  

Primero se debe asegurar que se cuenta con [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7), [Docker](https://docs.docker.com/) (23.0.5) y la imagen de docker pipelinesinmegen/pipelines_inmegen:latest.

 1. Seleccionar una ruta y el nombre para el directorio de salida
 2. Después generar el archivo sample_info.csv con la información que se describe en la sección - Formato del archivo con la información experimental -
 3. Editar el archivo de nextflow.config con la siguiente información:

        - Ruta de los archivos *fastq*
        - Ruta del directorio de salida de nextflow
        - Nombre del proyecto 
        - Ruta del índice de STAR
        - Ruta del archivo sample_info.csv
        - Nombre del índice de STAR
        - Ruta del directorio con los archivos de snps e indeles conocidos
        - Ruta del directorio de annovar
        - Condiciones del análisis (número de núcleos a utilizar por proceso, número de procesos simultáneos e información adicional)

  4. Ejecutar el comando: 

                bash run_nextflow.sh /path/to/out/dir

#### Formato del archivo con la información experimental

Para tener un buen control de los archivos a procesar (formato fastq pareados {Read_1,Read_2}), en el archivo sample_info.tsv incluir la siguiente información por columna:
 
		Sample_ID	Sample_name	RG	PU	R1	R2

 - Sample_ID   = Nombre completo o identificador de las muestras (debe contener el nombre o identificador, el número de muestra y el número de lane)
 - Sample_name = Nombre de los archivos o identificador de la muestra (comprende el nombre o identificador y el número de muestra)
 - RG          = ReadGroup name -El nombre del grupo de lectura de las muestras a procesar
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

Para una mayor descripción de la información del pipeline ejecutado se anexa el siguiente diagrama de flujo basado en las buenas prácticas de gatk.

[Consulta el las buenas pŕacticas de gatk para la identificación de variantes de datos de RNAseq](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-) 


![Flujo identificación de variantes rnaseq](../flowcharts/flow_vc-rnaseq.png)
