# Cuantificación y análisis de expresión diferencial
## Descripción de los archivos de salida del flujo de trabajo 

Como parte del servicio de análisis bioinformático del INMEGEN [Q&DEA]
Se entregarán los siguientes directorios con los siguientes archivos:

- Directorio: Archivos_alineados


- Directorio: Matrices_dcuestas


- Directorio: Genes_DEA 


- Directorio: Reportes de calidad

Esta carpeta contiene dos archivos que resumen las métricas de calidad del análisis.

El archivo **[nombre del proyecto]_multiqc.html** contiene el resumen del reporte de calidad de las lecturas después de que se eliminaron los adaptadores y las lecturas de mala calidad para par de archivos **Fastq** (R1 + R2). También incluyen el número de lecturas y bases on target junto a la profundidad de cada muestra (archivo bam).y algunas métricas adicionales de la calidad del alineamiento. 


Adicionalmente, en esta carpeta se encuentra el **Subdirectorio FastQC** donde se encuentra el reporte de calidad obtenido con FastQC de cada pareja de archivos Fastq (R1 + R2). 
