# Pipelines INMEGEN
## Flujos de trabajo automatizados con NextFlow

Este repositorio contiene diversos flujos de trabajo (pipelines) desarrollados y automatizados en el Instituto Nacional de Medicina Genómica (INMEGEN).
El principal objetivo de estos pipelines es el procesamiento de datos provenientes de secuenciación masiva (WGS/WES, RNA-seq)

Los flujos de trabajo son:


 - Cuantificación y Análisis de expresión diferencial (RNA-seq)

 - Indentificación conjunta de variantes Germinal (WGS/WES)

 - Indentificación de variantes somáticas (WGS/WES)

 - Indentificación de variantes RNA-seq (experimental)
 
 - Data_preprocessing (preprocesamiento de archivos fastq de GATK)
 

Las instrucciones necesarias para ejecutar cada flujo de trabajo se encuentra en su directorio correspondiente.

El directorio modules contiene todos los procesos que se usan en cada uno de los pipelines. 

Los diagramas de flujo correspondientes a cada flujo de trabajo se encuentran en la carpeta imagenes 

### En caso de utilizar los flujos de trabajo de este repositorio
En caso de utilizar cualquiera de estos pipelines solicitamos incluir la siguiente frase en los productos académicos generados: “Agradecemos a la Subdirección de Genómica Poblacional y a la Subdirección de Bioinformática del Instituto Nacional de Medicina Genómica por proveer flujos de trabajo utilizados de forma parcial o total como parte del análisis de este trabajo. (We acknowledge the population genomics and the bioinformatics departments from the National Institute of Genomic Medicine for providing workflows that were, either partially or completely, used as part of the analysis in this work )”

Si quieres aplicar alguno de estos flujos de trabajo con el apoyo de nuestro personal para implementarlo en tus datos este se considerará un servicio, por lo tanto, se cobrará de acuerdo a los tabuladores existentes en la cartera de servicios INMEGEN. NOTA: por un tiempo limitado, este servicio estará disponible gratuitamente para personal interno del INMEGEN.

En ningún caso se incluye interpretación de resultados o generación de resultados más allá de los descritos como parte del flujo de trabajo.

En caso de querer iniciar una colaboración académica con alguno de los miembros de este proyecto favor de contactarnos directamente.

No está permitida la utilización de estos flujos de trabajo con fines comerciales por terceros.  

### Equipo de desarrollo
Daniel Pérez-Calixto [dperez@inmegen.gob.mx](dperez@inmegen.gob.mx)

Laura Gómez-Romero [lgomez@inmegen.gob.mx](lgomez@inmegen.gob.mx)

Alejandra Cervera Taboada [acerverat@inmegen.gob.mx](acerverat@inmegen.gob.mx)

### Uso de Docker
Para garantizar la reproducibilidad y repetibilidad de los flujos de trabajo, estos se automatizaron con NexFlow y se contenizaron con Docker. El Dockerfile que se utilizó en para generar las imágenes empleadas se encuentran en el directorio docker/

### Contacto
Cualquier duda o comentario escribir a [nuestro correo de contacto](dperez@inmegen.gob.mx)
