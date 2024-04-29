# Anotación de variantes (archivos vcf) con Annovar

Este pipeline realiza la anotación de variantes a partir de archivos VCF.

**Nota:** Para poder utilizar annovar es necesario conseguir tu propia liga en el siguiente [link](https://www.openbioinformatics.org/annovar/annovar_download_form.php) y descargar tus propias bases de datos para realizar la anotación.

## Instrucciones de uso 

### Preparar el ambiente de trabajo

1. Te debes asegurar de contar con las siguientes herramientas informaticas:
        - [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7)
        - [Docker](https://docs.docker.com/) (23.0.5)
        - Imagen de docker pipelinesinmegen/pipelines_inmegen:public, la puedes clonar con el comando:

		docker pull pipelinesinmegen/pipelines_inmegen:public

2. Asegurarse de contar con las siguientes bases de datos, necesarias para el pipeline:
	- refGene
 	- ensGene
   	- avsnp150
   	- clinvar
   	- gnomad312_genome
   	- cosmic92
   	- dbnsfp33a

Consulta esta [liga](https://annovar.openbioinformatics.org/en/latest/user-guide/startup/#a-useful-tutorial) para más información.

### Ejecutar el flujo de trabajo

Para correr este pipeline se deben de ejecutar las siguientes instrucciones:

 1. Generar el archivo sample_info.tsv con la información que se describe en la sección - Formato del archivo con la información de las muestras -
 2. Editar el archivo de nextflow.config con la siguiente información:
	- Ruta de los archivos *vcf*
	- En caso de utilizar un vcf producto del llamado conjunto de variantes [VC-Germinal] seleccionar multiple_samples=true
	- Ruta del directorio de salida de nextflow
	- Nombre del proyecto 
	- Ruta de la referencia
	- Ruta del archivo sample_sheet.tsv
	- Nombre de la referencia
	- Ruta del directorio de annovar
	- Número de núcleos por proceso (parámetro runOptions)
	- Número de procesos que se ejecutarán de forma simultánea (parámetro queueSize)

  3. Ejecutar el comando correspondiente a cada subflujo de trabajo: 

         bash run_nextflow.sh /path/to/out/dir

### Formato del archivo con la información experimental
 
Para tener un buen control de los archivos a procesar (formato vcf), el archivo sample_info.tsv debe de incluir la siguiente información por columna:
 
	Sample	Path	Index
  	Sample_name	/paht/to/vcf_file	/path/to/vcf_index

 - Sample   = Nombre completo de los archivos vcf, se recomienda el formato [nombre de la muestra - número de muestra]
 - Path     = Ruta absoluta del archivo vcf de la muestra sample

#### Las herramientas utilizadas para correr este flujo de trabajo son:

 - GATK (4.2.6.1)
 - Bcftools (1.14.0)
 - Annovar
