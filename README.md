# Pipelines INMEGEN
## Flujos de trabajo automatizados con NextFlow

Este repositorio contiene diversos flujos de trabajo (pipelines) desarrollados y automatizados en el Instituto Nacional de Medicina Genómica (INMEGEN).
El principal objetivo de estos pipelines es el procesamiento de datos provenientes de secuenciación masiva (WGS/WES, RNA-seq)

Los flujos de trabajo de acuerdo a los diferentes directorios de este repositorio son:


 - Cuantificación y Análisis de expresión diferencial (RNA-seq) [QDEA_RNAseq]

 - Identificación conjunta de variantes Germinal (WGS/WES) [VC-Germline]

 - Identificación de variantes somáticas (WGS/WES) [VC-Somatic]

 - Identificación de variantes RNA-seq (experimental) [VC-RNAseq]
 
 - Preprocesamiento de archivos fastq [Data_preprocessing]

 - Anotación de variantes [Annotacion]

 - Identificación conjunta de variantes germinal con bootstrapping [Bt-VC-Germinal]
 

Las instrucciones necesarias para ejecutar cada flujo de trabajo se encuentra en su directorio correspondiente [marcado entre corchetes].

Al estar automatizados con [NextFlow](https://www.nextflow.io/docs/latest/index.html), los flujos de trabajo están divididos en procesos, que se encuentran en el directorio modules. 

El diagrama de flujo correspondiente a cada pipeline se encuentra en la carpeta flowcharts. 

**Nota:** Particularmente, el pipeline de identificación conjunta de variantes germinal está pensando para utilizar la ref hg38 del bundle de gatk, ya que, contiene los archivos necesarios para realizar VQSR (variant quality score recalibration).
En el caso de no contar con estos archivos como lo es el caso de especies diferentes al humano, se recomienda utilizar el flujo con bootstrapping [Bt-VC-Germinal], ya que, este sólo necesita la referencia, el índice de la referencia y el índice de BWA creado a partir de la referencia. 

##  Instrucciones para ejecutar los pipelines 

Para ejecutar los pipelines se debe contar con [NextFlow](https://www.nextflow.io/docs/latest/index.html) (22.10.7) y [Docker](https://docs.docker.com/) (23.0.5)

Además, es necesario clonar la imagen de docker con el comando 

		docker push pipelinesinmegen/pipelines_inmegen:public

En caso de querer construir la imagen de docker con el Dockerfile que se encuentra en el directorio docker/ utilizar:

		docker build -t pipelines_inmegen:public -f Dockerfile .

y modificar el tag de la imagen con el comando:

               docker tag pipelines_inmegen:public pipelinesinmegen/pipelines_inmegen:public

Debes asegurarte que el directorio de docker cuente con suficiente espacio para generar la imagen (~ 6 GB)

Finalmente clonar el repositorio github de interés.

## Políticas de uso

Los flujos de este repositorio pueden ser descargados y utilizados sin restricciones para uso académico. En caso de utilizar cualquiera de estos flujos solicitamos incluir la siguiente frase en los productos académicos generados: “Agradecemos a la Subdirección de Genómica Poblacional y a la Subdirección de Bioinformática del Instituto Nacional de Medicina Genómica por proveer flujos de trabajo que han sido utilizados de forma parcial o total como parte del análisis de este trabajo (We acknowledge the Population Genomics and the Bioinformatics Departments from the National Institute of Genomic Medicine for providing workflows that were, either partially or completely, used as part of the analysis in this work )”

Si requieres el apoyo de nuestro personal para implementar alguno de estos flujos de trabajo en tus datos, este se considerará un servicio. Por lo tanto, se cobrará de acuerdo a los tabuladores existentes en la cartera de servicios INMEGEN. NOTA: por un tiempo limitado, este servicio estará disponible gratuitamente para personal interno del INMEGEN.

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
