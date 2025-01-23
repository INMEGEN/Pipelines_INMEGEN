# Identificación de variantes somáticas a partir de datos WES/WGS
## Descripción de los archivos de salida del flujo de trabajo [VC-Somatic]

Como parte de los servicios de análisis bioinformáticos del Inmegen, después de ejecutar el flujo de trabajo se entregarán los siguientes directorios con los siguientes archivos:

### **Los archivos entregables del análisis se colocaron en los siguientes directorios:**
Folio_Proyecto
├── Alineamientos
├── Analisis_de_calidad
│   └── Fastqc
├── Fastq
└── Resultados
    ├── Reportes_de_calidad
    ├── Pannel_de_normales
    ├── Variantes
    └── Variantes_anotadas


#### - Directorio: **Alineamientos**

Este directorio contiene las secuencias alineados a hg38 (genoma humano versión GRCh38) **por muestra** en formato [bam](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm).

**MOTA:** Regularmente los archivos alineados son de un peso aproximado que oscila entre ~ 1 Gb a 20 Gb por lo que se recomienda elegir un lugar con suficiente espacio para la transferencia de dichos archivos

## El directorio Analisis_de_calidad contiene:

  - **[Folio]**_Calidad.html**
    Resumen del análisis de calidad de cada archivo FASTQ (R1 + R2).

  - Subdirectorio **Fastqc**:
    Reportes del análisis de calidad por cada archivo FASTQ.


#### - Direcotrio: **Resultados**

Este directorio contiene las siguientes subdirectorios:

- Subdirectorio **Reportes_de_calidad**:

    - **[Folio]_CalidadAnalisis.html**
      Métricas de calidad del análisis.
      Los nombres de los renglones que aparecen en la sección de estadísticos generales tiene la siguiente estructura por cada muestra:

                            _R1 y _R2 (archivos FASTQ)
                                |
                        _L001 (archivos FASTQ p/carril p/muestra)
                                |
                        _aligned_reads (archivo BAM p/muestra)
                                |
                        _sorted_dedup (archivo BAM ordenado y marcado p/muestra)
                        _sorted_dedup (archivo BAM ordenado y marcado p/muestra)
                                |
                        _SX (reporte profundidad del archivo sorted_dedup)
                                |
                        _snpEff_stats (archivo VCF filtrado y anotado)

       Cada tipo de archivo está asociado a las siguientes estadísticas:

                  ── _R1 y _R2 reporta para cada archivo FASTQ:
                        %Dups  (porcentaje de lecturas duplicadas)
                        %GC    (porcentaje promedio del contenido de GC)
                        M Seqs (número total de secuencias)

                  ── _L001 reporta para cada combinación muestra-carril:
                        % Duplication           (porcentaje de duplicados después del filtrado de lecturas por calidad)
                        M Reads After Filtering (número total de lecturas que pasaron el filtrado)
                        GC content              (porcentaje del contenido de GC después del filtrado)
                        %PF                     (porcentaje de lecturas después del filtrado)
                        %Adapter                (porcentaje de adaptadores)

                  ── _aligned_reads reporta para cada archivo alineado por muestra:
                        Duplication     (porcentaje de lecturas duplicadas)

                  ── _sorted_dedup reporta por muestra:
                        Insert Size     (la mediana del tamaño del inserto de los fragmentos de secuenciación por muestra)

                  ── _SX reporta para cada muestra:
                        >30X            (porcentaje del genoma con al menos una profundidad de 30X)
                        Median          (mediana de la profundidad por muestra)
                        Mean Cov        (promedio de la profundidad por muestra)

                  ── _snpEff_stats reporta para el llamado conjunto de todas las muestras:
                        Change rate     (tasa de cambio)
                        Ts/Tv           (tasa de Ts/Tv)
                        M Variants      (número de variantes después del filtrado)

      - **[Folio]_variantQC.html**
        Diversos estadísticos del número y tipo de variantes encontradas.
        Este reporte clasifica a las variantes en tres categorías:
        
        1. RAW (número total de variantes sin filtrar).
        2. Filtered (número de variantes que NO pasaron algún filtro de VQSR).
        3. Called (número de variantes que pasaron los filtros de VQSR y fueron marcadas con la bandera PASS).

     - **Resumen_cobertura.txt**
       Resumen rápido de la profundidad y lecturas on-target de las muestras.
       NOTA: En el caso de secuenciación de genoma completo las lecturas on-target son aquellas que han sido alineadas.

Los siguientes directorios contienen diversos archivos en formato de llamado de variantes [VCF](https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/VCF-Format.htm) con las variantes identificadas, los archivos VCF en los directorios de **Variantes** y **Variantes_anotadas**, además de pasar por el filtro de GATK [FilterMutectCalls](https://gatk.broadinstitute.org/hc/en-us/articles/360036856831-FilterMutectCalls), se filtraron por una profundidad de 10 y que tuvieran al menos 10 lecturas en el alelo de referencia, se proporciona un archivo VCF por muestra.

  - Subdirectorio: **Pannel_de_normales**
  - 
    **[Folio]_PON.vcf.gz**
    En caso de haber proporcionado 40 muestras normales, se proporciona el panel de normales creado. 

  - Subdirectorio: **Variantes**
 
    **[Folio]**_variantes.vcf.gz con las variantes identificadas de forma conjunta que pasaron los filtros de calidad.

    **NOTA**: Se incluye el subdirectorio variantes_por_muestra que contiene un archivo VCF por muestra

  - Subdirectorio: **Variantes_anotadas**
       **[Folio]_vars_anotadas.annovar.vcf.gz**:
        Variantes identificadas por muestra y anotadas con annovar. Los catálogos de genes utilizados son: refGene y ensGene, así como las bases de datos avSNP, CLINVAR, gnomAD, COSMIC y dbNSFP.

       **[Folio]_vars_anotadas.annovar.txt**:
       Misma información que el archivo NovaGermline_vars_anotadas.annovar.vcf.gz pero en un formato tabular.

**NOTA**: En caso de existir más de un alelo alternativo, este se coloca en un renglón diferente. Entonces, para la correcta interpretación del genotipo es necesario remitirse a la columna ALT (diferente de Alt) la cual describe todos los alelos encontrados en las distintas muestras.   

  **[Folio]_vars_anotadas.snpEff.vcf.gz**:
    Variantes identificadas por muestra y anotadas con SnpEff con el catálogo GRCh38.99.

**NOTA**: Se incluye el subdirectorio variantes_por_muestra que contiene un archivo VCF con las variantes anotadas por muestra.

Para más información de las bases de datos puedes consultar la siguiente [liga](https://annovar.openbioinformatics.org/en/latest/user-guide/filter/#overview).

**NOTA**: Todos los archivos VCFs se entregan comprimidos en un formato bgzip.

Para cualquier duda o información adicional favor de escribir al correo: serviciosbioinfo@inmegen.edu.mx
(END)

