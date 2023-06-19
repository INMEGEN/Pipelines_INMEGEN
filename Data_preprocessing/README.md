# Flujo de trabajo: Pre-procesamiento de datos para el descubrimiento de variantes con GATK4 (pipeline pre-GATK)

Con este pipeline se procesan archivos de secuenciación masiva (WGS/WES) en un formato *fastq*. 
Con este flujo de trabajo de obtiene un archivo *.bam limpio como lo indica el siguiente [link](https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently) de GATK, lo que es parte del pre-procesamiento de datos para el descubrimiento de variantes cuyo tutorial de gatk se encuentra [aquí](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)

#### Este flujo de trabajo requiere de las siguientes herramientas:
> 
> - NextFlow (22.04.5)
> - FastQC (0.11.9) 
> - MultiQC (1.11)
> - Openjdk (11.0.13 o superior)
> - GATK (4.2.6.1)
> - BWA (0.7.17-r 1188)
> - Picard Tools (2.0.1)
> - Samtools (1.6.0)
> 

Es necesario contar un un archivo *.txt con la información del arreglo experimental. También, es necesario el índice de [BWA](http://bio-bwa.sourceforge.net/) que corresponde al genoma de de referencia a emplear.  

## Instructivo de uso

## Para Ejecutar por primera vez el pipeline 
###Microambiente conda nf-gatk
 
 1. Verificar que el entorno de trabajo cuente con las herramientas descritas en la sección herramientas necesarias, en caso de que no se tengan instaladas, se recomienda generar un microambiente conda de NextFlow e instalar ahí las herramientas necesarias (excepto, GATK, Picard y Snpeff, de estas se especifica la ruta en el archivo de configuración). Para realizar lo anterior se recomienda seguir los siguientes pasos:
			
			conda install -c bioconda nextflow 
			conda install -c bioconda bwa
			conda install -c bioconda samtools
			conda update samtools 
			
 **Nota:** Es importante actualizar SamTools para poder instalar todas sus herramientas
 2. Descargar GATK (4.2.6.1) y Picard (2.0.1)
 3. Verificar que los archivos de nextflow estén disponibles (descargar de este repositorio), los archivos necesarios son:
	- archivo: main.nf
	- archivo: modules.nf
	- archivo: nextflow.config
	- archivo: sample_sheet.txt
 4. Verificar la ruta del directorio de Openjdk 8 (ver sección microambiente conda para GATK)
 5. Archivos a analizar en un formato *fastq*
 6. La ruta con los índices de BWA

Una vez que se cuenta con lo enunciado en las secciones anteriores se puede proceder a correr el pipeline pre-GATK siguiendo los siguientes pasos:

## Ejecutar flujo de trabajo VC_GATK

 1. Seleccionar un nombre para el directorio donde se concentrarán los archivos del análisis, algunas recomendaciones son: utilizar la fecha del análisis y una breve descripción del diseño experimental, p.ej. AAAA-M-D_DEA-nombre del ensayo. Esta carpeta se denominará directorio principal
 2. Para generar el archivo sample_sheet.txt es necesario contar con un archivo con la información del arreglo experimental. Sin embargo, sólo es neceario que el archivo contenga la siguiente información: 
 
			Sample	RG	PU	R1	R1
Dónde: 
 - Sample = nombre de la muestra
 - RG = nombre del grupo de lectura (necesario para el uso de GATK)
 - PU = nombre de la unidad donde se llevó a cabo la secueciación (aquí se agrega el lane)
 - R1 = ruta al archivo fastq R1
 - R2 = ruta al archivo fastq R2
			
 3. Editar el archivo de nexflow.config con la siguiente información:
	- Ruta al archivo sample_sheet.txt
	- Ruta del índice de BWA (*.fasta)
	- Ruta del ejecutable de GATK
	- Ruta del ejecutable de Picard
 4. Antes de ejecutar el pipeline activar el microambiente conda denominado nf-gatk utilizando:
 
			conda activate nf-gatk
			
 5. En el directorio donde están alojados los archivos de nexflow correr el siguiente comando para correr el pipeline
 
			nextflow run <pipeline name>
			
Si se quiere un informe del rendimiento de Nextflow se puede optar por:
 
			nextflow run <pipeline name> -with-report <file name>
			
Una vez ejecutado el pipeline, verificar que la carpeta de resultados dentro del directorio de salida contenga los archivos *.VCF y los informes *.html.
 
**NOTA 1:** En caso de ocurrir un error debido a una ruta mal escrita o algún evento similar, se puede editar los archivos de nextflow corrigiendo los errores, guardarlo de nuevo y reanudar el pipeline del análisis utilizando la opción -resume, como se muestra a continuación: 

			nextflow run <pipeline name> -resume

### Microambiente conda para GATK 

En este microambiente se pueden correr las diferentes herramientas de GATK que utilizan la versión 8 de Openjdk. Además, en este microambiente se puede correr el pipeline de [Mohammed Khalfan](https://github.com/gencorefacility/variant-calling-pipeline-gatk4). Sin embargo, es un ambiente alterno al que se utiliza para correr el pipeline VC_GATK

Los pasos para montar este microambiente son:
1. Instalar GATK siguiendo las instrucciones de este [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4)
2. Activar el ambiente conda GATK con el siguiente comando

			conda activate GATK
			
3. Instalar la correcta versión de Java 8, utilizando el siguiente comando
			
			conda install openjdk==8.0.332=h166bdaf_0 
			
4. Utilizar conda para instalar BWA, Samtools y Nextflow (siguiendo ese orden) con los comandos

			conda install -c bioconda bwa
			conda install -c bioconda samtools
			conda install -c bioconda nextflow

Como se puede ver en el orden de los pasos, en este microambiente NextFlow es la última herramienta que se instala (por cuestiones de compatibilidad), mientras que en el microambiente nf-gatk es la primera 

