# Flujo de trabajo cuantificación y análisis de expresión diferencial (pipeline Q&DEA)

Con este pipeline se analizan archivos de secuenciación masiva (*RNA-seq*) en un formato *fastq*. Se puede obtener sólo la matriz de cuentas TPM o la matriz de cuentas más un análisis de los genes diferencialmente expresados entre dos condiciones distintas.

#### Las herramientas necesarias para correr este flujo de trabajo son:
> 
> - NextFlow (22.04.5)
> - R (4.2.0) 
> - FastQC (0.11.8) 
> - MultiQC (1.13.deb0)
> - Trim Galore (0.6.7) 
> - Kallisto (0.44.0)
> 
Se recomienda utilizar conda para instalar estas herramientas

Además de los programas arriba enunciados, son necesarias las siguientes librerías:

#### Librerías de R:
> 
> - Tximport (1.22.0)
> - readr (2.1.2)
> - BUSpaRse (1.8.0)
> - tximportData (1.22.0)
> - dplyr (1.0.9)
> - DESeq2 (1.34.0)
> - ggplot2 (3.3.6)
> - ggrepel (0.9.1)
> - PCAtools (2.6.0)
> - EnhancedVolcano (1.13.2)
> - optparse (1.7.1)
> - ComplexHeatmap (2.10.0)
> - topGO (2.46.0)
> - GeneTonic (1.6.4)
> - org.Hs.eg.db (3.14.0)
> - org.Mm.eg.db (3.15.0)
> - AnnotationDbi (1.56.2)
> - SummarizedExperiment (1.24.0)
 
También, debe haberse instalado previamente las librerías rhdf5 (2.38.1) y pcaExplorer (2.22.0)
Estas librerías se pueden instalar en R utilizando BiocManager de Bioconductor
  

Adicionalmente, es necesario contar con un archivo csv con la información del arreglo experimental, el índice de [kallisto](https://pachterlab.github.io/kallisto/manual) del transcriptoma de referencia deseado y del archivo *.gtf con las anotaciones de los genes del transcriptoma de referencia.  

## Instructivo de uso

### Para Ejecutar por primera vez el pipeline 
 
 1. Verificar que el entorno de trabajo cuente con las herramientas descritas en la sección herramientas necesarias, en caso de que no se tengan instaladas, se recomienda generar un microambiente conda de NextFlow e instalar ahí las herramientas necesarias
 2. Verificar que los archivos de nextflow estén disponibles (descargar de este repositorio), los archivos necesarios son:
	- archivo: main.nf
	- archivo: modules.nf
	- archivo: nextflow.config
 3. Verificar que los scripts de R que utilizan Tximport y DESeq2 se encuentre disponibles y en un directorio conocido (descargar de este repostorio)
 4. Contar con el archivo csv que contiene la información del arreglo experimental (ver sección **formato del archivo con la información experimental**)
 5. Contar con el índice de kallisto y el archivo gtf del genoma utilizado
 6. Archivos a analizar en un formato *fastq*

Debido a que se utiliza DESeq2 se recomienda al menos tener 3 réplicas por condición.

Una vez que se cuenta con lo enunciado en las secciones anteriores se puede proceder a correr el pipeline Q&DEA siguiendo los siguientes pasos:

### Ejecutar flujo de trabajo Q&DEA

 1. Editar el archivo de nexflow.config con la siguiente información:
        - Ruta del directorio de salida
	- Ruta de los archivos fastq
        - Ruta del índice de kallisto (*.idx)
        - Ruta del archivo *.gtf
        - Ruta del archivo *.csv con la información experimental
	- Ruta del archivo cvs con la información del arreglo experimental
 2. Verificar que el archivo csv contenga las condiciones a comparar y estos sean los mismos en el archivo de configuración
 6. En el archivo de configuración definir los umbrales de los valores p-value y Log2FC
 7. Verificar los nombres de los archivos de salida (archivo de cuentas, archivo de cuentas TPM, tabla de genes, tabla de genes diferencialmente expresados y los nombres de las gráficas). Se recomienda que los archivos con las listas de genes se nombren con la información de la comparación realizada p.ej. result-tratado-vs-control.csv
 8. En el directorio donde están alojados los archivos de nexflow ejecutar el comando:
 
			nextflow run <pipeline name>
			
Si se quiere un informe del rendimiento de Nextflow se puede optar por:
 
			nextflow run <pipeline name> -with-report <file name>
			
Una vez ejecutado el pipeline, verificar que la carpeta de resultados dentro del directorio de salida contenga los archivos *.csv y las gráficas generadas.
 
**NOTA 1:** En caso de ocurrir un error debido a una ruta mal escrita o algún evento similar, se puede editar los archivos de nextflow corrigiendo los errores, guardarlo de nuevo y resumir el pipeline del análisis utilizando la opción -resume, como se muestra a continuación: 

			nextflow run <pipeline name> -resume

### Formato del archivo con la información experimental 

El archivo csv con la descripción del diseño experimental de las muestras debe de contener al menos la siguiente información:

> 
> - name_id = nombre común de identificación de las muestras   
> - names = nombres de las muestras
> - condition = describe la condición experimental de las muestras 
> - adicional = información adicional del arreglo experimental, estas son otras columnas que agregan información sobre las muestras, por ejemplo: # de réplica, lote, etc.
> 
**Nota:** Es imprescindible conservar los nombres de las columnas name_id, names y condition 
