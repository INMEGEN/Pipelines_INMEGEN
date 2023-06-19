#!/bin/sh
# Directorio de salida de nextflow
path=$1

## Ejecutar nextflow
nextflow run main.nf -resume -with-trace trace_pon.txt -with-report report_pon.html -with-timeline timeline_pon.html

mkdir -p $path/run_files/modules

## Copiar los ejecutables, la configuración y los modulos de nextflow utilizados
cp ../../modules/VC-Somatic/PON* $path/run_files/modules
cp main.nf nextflow.config sample_sheet_normales.tsv $path/run_files

## Mover los reportes de nextflow a una carpeta en el directrio de salida de nextflow
mv trace_pon.txt timeline_pon.html report_pon.html $path/run_files

echo -e "Archivos y reportes del flujo de trabajo copiados a $path/run_files \n"
echo -e "Para explorar los resultados ejecuta: cd $path/out \n"
echo -e "Flujo de trabajo completado"
