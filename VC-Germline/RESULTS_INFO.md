# Identificación conjunta de variantes germinales a partir de datos WGS/WES

## Descripción de los archivos de salida del flujo de trabajo [VC-Germinal]

Como parte de los servicios de análisis bioinformáticos del INMEGEN, después de ejecutar el flujo de trabajo se entregarán los siguientes directorios con los siguientes archivos:

- Directorio: **Archivos_bam**

Este directorio contiene los alineados a hg38 (genoma humano versión GRCh38) **por muestra** en formato [bam](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm).
Estos archivos está ordenados y con los duplicados ya marcados. 

- Directorio: **VCFs_filtrados**

Esta carpeta contiene diversos archivos en formato [VCF](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/VCF-Format.htm) con las variantes identificadas y marcadas con los filtros de GATK [VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR) de todas las muestras.
Esta carpeta también contiene el subdirectorio **vcfs_pormuestra** donde se proporciona un archivo VCF por muestra.

**NOTA:** Dependiendo de la configuración de la solicitud del servicio, puede haber un archivo por tipo de variante (SNPs o INDELs) o un archivo que contenga ambas (SNPs + INDELs).

- Directorio: **VCFs_anotados** 

Esta carpeta contiene diversos archivos en formato [VCF](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/VCF-Format.htm) con las variantes identificadas, que pasaron los filtros de VQSR (bandera **PASS**) y anotadas con los catálogos de genes refGene y ensGene, junto con las bases de datos avSNP, CLINVAR, gnomAD, COSMIC y dbNSFP utilizando Annovar. 
Para más información de las bases de datos utilizadas consultar la siguiente [liga](https://annovar.openbioinformatics.org/en/latest/user-guide/filter/#overview). 
Esta carpeta también contiene el subdirectorio **vcfs_pormuestra** donde se proporciona un archivo VCF por muestra.  

También, se incluye un archivo de texto separado por tabulador (/t) que contiene la información del VCF anotado.

**NOTA:** Dependiendo de la configuración de la solicitud del servicio, puede haber un archivo por tipo de variante (SNPs o INDELs) o un archivo que contenga ambas (SNPs + INDELs). 

- Directorio: **Reportes de calidad**

Esta carpeta contiene dos archivos que resumen las métricas de calidad del análisis.

El archivo **[nombre del proyecto]_multiqc.html** contiene el resumen del reporte de calidad de las lecturas después de que se eliminaron los adaptadores y las lecturas de mala calidad para cada par de archivos **Fastq** (R1 + R2). También, se incluyen el número de lecturas y bases on target junto a la profundidad de cada muestra (archivo bam) y algunas métricas adicionales de la calidad del alineamiento. 

El archivo **[nombre del proyecto]_variantqc.html** contiene un resumen del número de variantes encontradas en conjunto y por muestra. También, contiene métricas por tipo de variante (SNPs e INDELs). 
Este reporte clasifica a las variantes en tres categorías:

-RAW son las variantes sin filtro
-Filtered son las variantes que no pasaron algún filtro de VQSR
-Called son las variantes que pasaron los filtros de VQSR (bandera **PASS**) 

Adicionalmente, en esta carpeta se encuentra el **Subdirectorio FastQC** donde se encuentra el reporte de calidad obtenido con FastQC de cada pareja de archivos **Fastq** (R1 + R2). 

**NOTA:** Todos los archivos VCFs se entregarán compresos en un formato bgzip.
