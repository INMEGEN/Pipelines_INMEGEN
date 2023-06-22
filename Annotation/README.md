# Flujo de trabajo identificación de variantes somáticas utilizando NextFlow y GATK

Este pipeline realiza la anotación de variantes a partir de archivos vcfs filtrados.
En caso de trabajar con el genoma hg38, los archivos de referencia se pueden descargar del [bundle de GATK](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).  

**Nota:** Para poder utilizar annovar es necesario conseguir tu propia liga en el siguiente [link](https://www.openbioinformatics.org/annovar/annovar_download_form.php).

## Instrucciones de uso 

Primero se debe asegurar que se cuenta con [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7), [Docker](https://docs.docker.com/) (23.0.5) y la imagen de docker pipe
linesinmegen/pipelines_inmegen:latest.

Estos flujos de trabajo utilizan archivos vcf previamente procesados con el pipeline identificación de variantes germinal, somáticas o RNAseq.

 1. Seleccionar una ruta y el nombre para el directorio de salida
 2. Después generar el archivo sample_*.tsv con la información que se describe en la sección - Formato del archivo con la información de las muestras -
 3. Editar el archivo de nextflow.config con la siguiente información:

	- Ruta de los archivos *vcf*
	- En caso de utilizar un vcf producto del llamado conjunto de variantes [VC-Germinal] seleccionar multiple_samples=true
	- Ruta del directorio de salida de nextflow
	- Nombre del proyecto 
	- Ruta de la referencia
	- Ruta del archivo sample_sheet.tsv
	- Nombre de la referencia
	- Ruta del directorio de annovar (en caso de uso académico)
	- Condiciones del análisis (número de núcleos a utilizar por proceso, número de procesos simultáneos e información adicional)

  4. Ejecutar el comando correspondiente a cada subflujo de trabajo: 

                bash run_nextflow.sh /path/to/out/dir

### Formato del archivo con la información experimental
 
Para tener un buen control de los archivos a procesar (formato bam), el archivo sample_*.tsv debe de incluir la siguiente información por columna:
 
		Sample	Path	

 - Sample   = Nombre completo de los archivos e identificador de las muestras normales (debe de contener el nombre o identificador de la muestra)
 - Path     = Ruta absoluta del archivo vcf de la muestra sample

#### Las herramientas utilizadas para correr este flujo de trabajo son:
>
> - GATK (4.2.6.1)
> - Bcftools (1.14.0)
> - Annovar
>
