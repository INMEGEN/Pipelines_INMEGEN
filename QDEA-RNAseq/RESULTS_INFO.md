[![INMEGEN](../flowcharts/inmegen_t1.png)](https://www.inmegen.gob.mx/)

# Servicios de Análisis Bioinformáticos
## Análisis de expresión diferencial (QDEA)

Descripción de los entregables del flujo de cuantificación y análisis de expresión diferencial

Los archivos entregables del análisis se colocaron en los siguientes directorios:

    Folio de Proyecto/ 
    ├── Alineamientos 
    ├── Analisis_Calidad 
    │   ├── Inspeccion_Secuencias 
    │   └── Reportes_Calidad 
    ├─── Resultados 
    │   ├── Expresion 
    │   │   ├── Genes 
    │   │   └── Transcritos 
    │   ├── Expresion_Diferencial 
    │   └── Reportes_Calidad 
    ├── NextFlow 
    │   └── modules 
    └─── R  

## Contenido de cada uno de los directorios 

### Directorio: Alineamientos

 - Secuencias alineadas al genoma humano versión GRCh38, resultados en formato BAM por muestra.

**NOTA:** Cada archivo de alineamiento tiene un peso considerable, por lo que es necesario elegir un lugar con suficiente espacio de disco para la descarga de estos archivos.

### Directorio: Analisis_Calidad

  - **Reporte_Calidad.html***

  Resumen del análisis de calidad de la secuenciación (FASTQ: R1 + R2). Informe generado con MultiQC con las salidas de FastQC y Fastq_Screen.

#### Subdirectorio, **Inspeccion_Secuencias**:
      
 - Archivos de salida de **Fastq_Screen** con la inspección completa del origen de las secuencias utilizando genomas de bacterias, hongos y algunas especies comunes.

#### Subdirectorio, Reportes_Calidad:
      
 - Archivos de salida de **FastQC** por cada archivo (FASTQ: R1 + R2).

**IMPORTANTE:** El análisis de calidad y la inspección de las secuencias realizados con **FastQC** y **Fastq_Screen** se llevaron a cabo sobre los archivos **FASTQ** sin ningún tratamiento previo.

### Directorio: Resultados 

 - **Reporte_de_Calidad_Analisis.html**
   Archivo generado con MultiQC. Este informe incluye diversas métricas de calidad y alineamiento generadas por las siguientes herramientas bioinformáticas:

     1 **FastQC** y **Fastp**: Evaluación de la calidad de las lecturas después del recorte de adaptadores y la eliminación de lecturas de baja calidad.
   
     2 **STAR:** Métricas de alineamiento al genoma de referencia.
   
     3 **QualiMap:** Comprobación de la cobertura y distribución de las regiones del genoma alineadas.
   
     4 **Salmon:** Métricas de pseudomapeo.
   
     5 **Análisis de Expresión:** Correlación de Spearman y análisis de componentes principales (PCA) entre las diferentes muestras.

#### Subdirectorio: Expresion

Contiene los siguientes subdirectorios.

#### Subdirectorio: Expresion/Genes
Cuantificación de lecturas a nivel de gen, este subdirectorio contiene:

 - **fmcounts.tsv**

   Matriz de cuentas utilizada para realizar el análisis de expresión diferencial.
    
 - **fmcounts_tpm.tsv**

   Matriz de cuentas TPM.

   **NOTA:** La cuantificación a nivel de gen se realizó con FeactureCounts

#### Subdirectorio: Expresion/Transcritos
Cuantificación de lecturas a nivel de transcrito, este subdirectorio contiene:

 - **cuentastx.tsv**

   Matriz de cuentas a nivel de transcrito.
    
 - subdirectorio por muestra con la cuantificación de **salmon**.

   Archivo quant.sf
   
   **NOTA:** La cuantificación a nivel de gen se realizó con Salmon.

#### Subdirectorio: Expresion_Diferencial

 - PCA_mqc.png
 
 Representación gráfica de las dos primeras componentes del análisis de componentes principales (PCA) aplicado a todas las muestras.

 - Correlacion_entre_muestras_mqc.png
 
   Mapa de calor con la correlación de Spearman calculada entre todas las muestras.

#### Subdirectorio por cada comparación realizada, contiene:
 
 - **DEG_DESeq2_filtrados.tsv**

   Tabla de genes diferencialmente expresados filtrados por **|Log2FC| > 1** y ** p-ajustado < 0.05**.

 - **DEG_DESeq2.tsv**

   Tabla completa con los estadísticos del análisis de expresión diferencial generado por **DESeq2**.

 - **heatmap_log2.png**

   Mapa de calor agrupado por condición, basado en el escalamiento mediante logaritmo en base 2 (**log2**).

 - **heatmap_zscore.png**

   Mapa de calor agrupado por condición, basado en el **z-score**, una medida estadística que indica cuántas desviaciones estándar se encuentra cada valor de expresión por encima o por debajo de la media.

 - **volcano_plot.png**

   Gráfico de volcán para visualizar los genes diferencialmente expresados en función de **Log2FC** y **p-ajustado**.

 - **DESeq2_dds.rds**

   Objeto de R de **DESeq2** que almacena los datos de expresión génica y la información experimental asociada.

 - **R_sessionInfo.log**
                 
   Información de la sesión de R, contiene un resumen del número de genes con un recuento total de lecturas distinto de cero y las versiones de las librerías de R utilizadas.

#### Subdirectorio: Reportes_Calidad

 Archivos generados por **FastQC** para cada archivo (FASTQ: R1 + R2) después de la eliminación de adaptadores y secuencias de baja calidad.

### Directorio: NextFlow

 - Archivo **main.nf**

   Archivo de NextFlow con el flujo de trabajo principal del análisis.

 - Archivo **nextflow.config**

   Archivo de NextFlow con la configuración del análisis.

 - Archivos **sample_info.tsv** y **metadata.tsv**

  Archivos con la información de las muestras.

#### Subdirectorio modules:
 
 Contiene los procesos utilizados para ejecutar el flujo principal de trabajo (archivo main.nf).

### Directorio: R

 - **DEA.R**

  Script de R utilizado para realizar la expresión diferencial.

 - **Q.R**

  Script de R utilizado para realizar la cuantificación a nivel de transcrito.

## Recursos utilizados

 - Referencia DNA:

https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

 - Referencia cDNA:

https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

 - Archivo GTF:
https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz

### Versiones de las herramientas utilizadas en el flujo de analisis QDEA

  - Fastp 0.23.4
  - FastQC v0.12.1
  - STAR 2.7.11b
  - QualiMap v.2.3
  - FeatureCounts v2.0.8
  - Salmon 1.10.3
  - Multiqc v1.25.1
  - R 4.4.2
  - NextFlow v4.10.4.5934
  - Docker v25.0.1
