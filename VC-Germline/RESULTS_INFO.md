# Identificación conjunta de variantes germinales a partir de datos WGS/WES

## Descripción de los archivos de salida del flujo de trabajo [VC-Germinal]

Como parte de los servicios de análisis bioinformáticos del Inmegen, después de ejecutar el flujo de trabajo se entregarán los siguientes directorios organizados de la siguiente manera:

    Folio de Proyecto/
    ├── Alineamientos
    ├── Analisis_de_calidad
    │   ├── Inspeccion_Secuencias
    │   └── Reportes_Calidad 
    └── Resultados
        ├── Reportes_de_calidad
        ├── Variantes
        └── Variantes_anotadas

### - Directorio: **Alineamientos**

Este directorio contiene los archivos alineados a hg38 (genoma humano versión GRCh38) **por muestra** en formato [bam](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm).

**MOTA:** Regularmente los archivos alineados son de un peso aproximado que oscila entre ~1 Gb a 20 Gb por lo que se recomienda elegir un lugar con suficiente espacio para la transferencia de dichos archivos

### - Directorio: **Analisis_Calidad**

  - **Reporte_Calidad.html**

  Resumen del análisis de calidad de la secuenciación (FASTQ: R1 + R2). Informe generado con MultiQC con las salidas de FastQC y Fastq_Screen.

#### -- Subdirectorio, **Inspeccion_Secuencias**:
      
 - Archivos de salida de **Fastq_Screen** con la inspección completa del origen de las secuencias utilizando genomas de bacterias, hongos y algunas especies comunes.

#### -- Subdirectorio, Reportes_Calidad:
      
 - Archivos de salida de **FastQC** por cada archivo (FASTQ: R1 + R2).

**IMPORTANTE:** El análisis de calidad y la inspección de las secuencias realizados con **FastQC** y **Fastq_Screen** se llevaron a cabo sobre los archivos **FASTQ** sin ningún tratamiento previo.

### - Direcotrio: **Resultados**

 - **Reporte_de_Calidad_Analisis.html**
   Archivo generado con MultiQC. Este informe incluye diversas métricas de calidad y alineamiento generadas por las siguientes herramientas bioinformáticas:

     1 **FastQC** y **Fastp**: Evaluación de la calidad de las lecturas después del recorte de adaptadores y la eliminación de lecturas de baja calidad.
   
     2 **GATK:** Métricas de alineamiento al genoma de referencia.
   
     3 **Mosdepth:** Comprobación de la cobertura y profundidad de las regiones del genoma alineadas.
   
     4 **SNPeff:** Métricas de la calidad de las variantes.

 - **Reporte_variantQC.html**
   Diversos estadísticos del número y tipo de variantes encontradas.

   Este reporte clasifica a las variantes en tres categorías:
   1. RAW (número total de variantes sin filtrar).
   2. Filtered (número de variantes que NO pasaron algún filtro de VQSR).
   3. Called (número de variantes que pasaron los filtros de VQSR y fueron marcadas con la bandera PASS).

 - **Resumen_cobertura.txt**
   Resumen rápido de la profundidad y lecturas on-target de las muestras.
   NOTA: En el caso de secuenciación de genoma completo las lecturas on-target son aquellas que han sido alineadas

#### -- Subdirectorio: **Variantes**

 - **[Folio]_variantes.vcf.gz**
   Archivo VCF con las variantes identificadas de forma conjunta que pasaron los filtros de [VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR) (bandera **PASS**) de todas las muestras.

**NOTA**: Se incluye el subdirectorio *variantes_por_muestra* que contiene un archivo VCF por muestra.

#### -- Subdirectorio: **Variantes_anotadas** 

 - **[Folio]_vars_anotadas.annovar.vcf.gz**
   Variantes identificadas de forma conjunta y anotadas con annovar. Los catálogos de genes utilizados son: *refGene* y *ensGene*, así como las bases de datos *avSNP*, *CLINVAR*, *gnomAD*, *COSMIC* y *dbNSFP*.

 - **[Folio]_vars_anotadas.annovar.txt**
   Misma información que el archivo NovaGermline_vars_anotadas.annovar.vcf.gz pero en un formato tabular.

**NOTA**: En caso de existir más de un alelo alternativo, este se coloca en un renglón diferente. Entonces, para la correcta interpretación del genotipo es necesario remitirse a la columna ALT (diferente de Alt) la cual describe todos los alelos encontrados en las distintas muestras.

 - **[Folio]_vars_anotadas.snpEff.vcf.gz**
   Variantes identificadas de forma conjunta y anotadas con SnpEff con el catálogo *GRCh38.99*.

**NOTA**: Se incluye el subdirectorio *variantes_por_muestra* que contiene un archivo VCF con las variantes anotadas por muestra.

**NOTA**: Todos los archivos VCFs se entregarán compresos en un formato bgzip. Para mayor información del formato de llamado de variantes (VCF) [consulte esta liga](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/VCF-Format.htm).
  
**NOTA**: Dependiendo de la configuración de la solicitud del servicio, también puede haber un archivo por tipo de variante (SNPs o INDELs) o un archivo que contenga ambas (SNPs + INDELs).

## Recursos utilizados

   - [Bundle GATK, genoma hg38](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/)

### Versiones de las herramientas utilizadas en el flujo de analisis QDEA
  - Fastp 0.23.4
  - FastQC v0.12.1
  - BWA MEM 
  - Mosdpth 
  - GATK
  - Picard 
  - BCFTools
  - SAMTools
  - SNPeff
  - Annovar 
  - Multiqc v1.25.1
  - R 4.4.2
  - NextFlow v24.10.4.5934
  - Docker v25.0.1
