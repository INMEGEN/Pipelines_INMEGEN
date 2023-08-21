# Flujo de trabajo: Pre-procesamiento de datos para el descubrimiento de variantes con GATK

Este pipeline realiza el preprocesamiento de archivos de secuenciación masiva (WGS/WES) en formato *fastq*. 
Con este flujo de trabajo de obtiene un archivo bam limpio como lo indica el siguiente [link](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently) de GATK, lo que es parte del pre-procesamiento de datos para el descubrimiento de variantes cuyo tutorial de gatk se encuentra [aquí](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
Este flujo de trabajo se utiliza para el preprocesamiento de archivos fastq para el flujo de trabajo identidicación de variantes somáticas. 

**Nota:** Por el momento el análisis sólo está disponible para datos de lectura corta (ilummina paired-end).

### Para solicitar este flujo de trabajo como servicio debes de entregar al personal de INMEGEN 

- Archivos de lectura fastq (Illumina paired-end).
- Archivo con la información experimental (identificador de la muestra, plataforma y librería de secuenciación, si son múltiples lanes especificar el número).
- En caso de WES especificar el kit utilizado y los identificadores de las muestras indicando si son normales o de tumor.

**Nota:** Esta información también es necesaria para utilizar los flujos de trabajo de sin asistencia.

## Instrucciones para ejecutar el pipeline

Si deseas utilizar este flujo de trabajo sin apoyo del personal del INMEGEN sigue las siguientes instrucciones:

Asegurarse que se cuenta con la instalación de [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7), [Docker](https://docs.docker.com/) (23.0.5) y la imagen de docker pipelinesinmegen/pipelines_inmegen:public.

 1. Seleccionar una ruta y el nombre para el directorio de salida.
 2. Para generar el archivo sample_info.tsv es necesario contar con un archivo con la información del arreglo experimental. Sin embargo, sólo es neceario que el archivo contenga la siguiente información: 
 
			Sample	RG	PU	R1	R1
Dónde: 
 - Sample = nombre de la muestra
 - RG     = nombre del grupo de lectura (necesario para el uso de GATK)
 - PU     = nombre de la unidad donde se llevó a cabo la secuenciación (aquí se agrega el lane y la librería)
 - R1     = ruta absoluta del archivo fastq R1
 - R2     = ruta absoluta del archivo fastq R2
			
 3. Editar el archivo de nexflow.config con la siguiente información:
	- Ruta del directorio de salida 
	- Ruta al archivo sample_info.tsv
	- Ruta del índice de BWA (*.fasta)
	- Nombre del índice de BWA
	- Ruta absoluta al archivo con la lista de intervalos (*.interval_list)
	- Nombre del proyecto
	- Indicar si son múltiples lanes por muestra 

 4. Ejecutar el comando: 

		bash run_nextflow.sh /path/to/out/dir

En caso de algún error en la ejecución modificar el origen del error y correr de nuevo el comando arriba descrito.

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
