# Análisis de calidad de archivos de secuenciación masiva en formato Fastq

## Descripción de los archivos de salida del flujo de trabajo [Fastq-QC]

Como parte de los servicios de análisis bioinformáticos del INMEGEN se hace una primera inspección de calidad de los archivos **Fastq** que permita en primera instancia conocer la calidad y alguna contaminación existente de las muestras a analizar.
Se entregarán los siguientes directorios con los siguientes archivos:

- Directorio: **Reportes de calidad**

Esta carpeta contiene un archivo html que resumen las métricas de calidad del análisis.

El archivo **[nombre del proyecto]_multiqc.html** contiene el resumen del reporte de calidad de las lecturas por cada pareja de archivos **Fastq** (R1 + R2) junto a los resultados de Fastq Screen 

Adicionalmente, en esta carpeta se encuentra el **Subdirectorio FastQC** donde se encuentra el reporte de calidad obtenido con FastQC de cada pareja de archivos **Fastq** (R1 + R2). 

**NOTA:** Todos los archivos VCFs se entregarán compresos en un formato bgzip.
