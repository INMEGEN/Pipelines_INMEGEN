[![INMEGEN](./flowcharts/inmegen_t1.png)](https://www.inmegen.gob.mx/)
#  Pipelines INMEGEN
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.7-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Flujos de trabajo automatizados 

Este repositorio contiene diversos flujos de trabajo (pipelines) desarrollados y automatizados en el Instituto Nacional de Medicina Genómica (INMEGEN).
El principal objetivo de estos pipelines es el procesamiento de datos provenientes de secuenciación masiva (Whole Genome Sequencing [WGS]/Whole Exome Sequencing [WES], RNA sequencing [RNA-seq])

Los siguientes directorios de este repositorio contienen un flujo de trabajo. En cada uno de estos directorios [marcado entre corchetes] se encuentran las instrucciones necesarias para ejecutar ese flujo de trabajo: 

 - Cuantificación y Análisis de expresión diferencial a partir de datos de RNA-seq [QDEA-RNAseq]

 - Identificación conjunta de variantes germinales a partir de datos WGS/WES [VC-Germline]
   
 - Identificación conjunta de variantes germinales a partir de datos WGS/WES con bootstrapping [BT-VC-Germinal]

 - Identificación de variantes somáticas a partir de datos WGS/WES [VC-Somatic]

 - Identificación de variantes a partir de datos de RNA-seq  [VC-RNAseq]
   
**Nota:** El pipeline [VC-Germline] requiere los archivos necesarios para realizar VQSR (variant quality score recalibration), por lo que  está diseñado para utilizar el bundle de gatk para el genoma humano hg38. En caso de no contar con estos archivos se recomienda utilizar el flujo [BT-VC-Germinal], ya que este sólo necesita la referencia, el índice de SAMTOOLS y el índice de BWA creados a partir de la referencia. 


Algunos directorios contienen pipelines de procesamiento que son requeridos por más de un flujo de trabajo
 
 - Preprocesamiento de archivos FASTQ [Data-preprocessing]

 - Anotación de variantes [Annotation]

Otros directorios contienen información general:

 - Los flujos de trabajo están divididos en procesos. Estos procesos se encuentran en el directorio **modules**.
   
 - El diagrama de flujo correspondiente a cada pipeline se encuentra en la carpeta **flowcharts**. 




##  Instrucciones para ejecutar los pipelines 

Para ejecutar los pipelines es necesario contar con [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7) y [Docker](https://docs.docker.com/) (23.0.5)

Además, es necesario clonar la imagen de docker de este repositorio con el comando 

		docker pull pipelinesinmegen/pipelines_inmegen:public

En caso de querer construir la imagen de docker con el Dockerfile que se encuentra en el directorio docker/ utilizar:

		docker build -t pipelines_inmegen:public -f Dockerfile .

y modificar el tag de la imagen con el comando:

               docker tag pipelines_inmegen:public pipelinesinmegen/pipelines_inmegen:public

Debes asegurarte de que el directorio de docker cuente con suficiente espacio para generar la imagen (~ 6 GB)

Finalmente, clonar el repositorio github de interés.

## Políticas de uso

Los flujos de este repositorio pueden ser descargados y utilizados sin restricciones para uso académico. En caso de utilizar cualquiera de estos flujos solicitamos incluir la siguiente frase en los productos académicos generados: “Agradecemos a la Subdirección de Genómica Poblacional y a la Subdirección de Bioinformática del Instituto Nacional de Medicina Genómica por proveer flujos de trabajo que han sido utilizados de forma parcial o total como parte del análisis de este trabajo (We acknowledge the Population Genomics and the Bioinformatics Departments from the National Institute of Genomic Medicine for providing workflows that were, either partially or completely, used as part of the analysis in this work )”

Si requieres el apoyo de nuestro personal para implementar alguno de estos flujos de trabajo en tus datos, este se considerará un servicio. Por lo tanto, se cobrará de acuerdo a los tabuladores existentes en la cartera de servicios INMEGEN. NOTA: por un tiempo limitado, estos servicios estarán disponibles gratuitamente para personal interno del INMEGEN.

En ningún caso nuestros servicios incluyen interpretación de resultados o generación de resultados más allá de los descritos como parte del flujo de trabajo.

En caso de querer iniciar una colaboración académica con alguno de los miembros de este proyecto favor de contactarnos directamente.

No está permitida la utilización de estos flujos de trabajo con fines comerciales por terceros.  

### Equipo de desarrollo
Daniel Pérez-Calixto [dperez@inmegen.gob.mx](dperez@inmegen.gob.mx)

Laura Gómez-Romero [lgomez@inmegen.gob.mx](lgomez@inmegen.gob.mx)

Alejandra Cervera Taboada [acerverat@inmegen.gob.mx](acerverat@inmegen.gob.mx)

## Uso de Docker
Para garantizar la reproducibilidad y repetibilidad de los flujos de trabajo, estos se automatizaron con [NextFlow](https://www.nextflow.io/docs/latest/index.html) y se contenerizaron con [Docker](https://docs.docker.com/). El Dockerfile utilizado para generar las imágenes empleadas se encuentran en el directorio docker/

## Contacto
Cualquier duda o comentario escribir a [nuestro correo de contacto](dperez@inmegen.gob.mx)
