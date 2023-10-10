#!/bin/sh
# Directorio de salida de nextflow
path=$1

## Ejecutar nextflow
nextflow run main.nf -resume -with-trace trace_dp.txt -with-report report_dp.html -with-timeline timeline_dp.html

mkdir -p $path/run_files/modules

## Copiar los ejecutables, la configuración y los modulos de nextflow utilizados
cp ../modules/qualitycontrol/* $path/run_files/modules
cp ../modules/data_preprocessing/* $path/run_files/modules
cp ../modules/metricts/* $path/run_files/modules
cp ../modules/common/* $path/run_files/modules
cp main.nf nextflow.config sample_info.tsv $path/run_files

## Mover los reportes de nextflow a una carpeta en el directrio de salida de nextflow
mv trace_dp.txt timeline_dp.html report_dp.html $path/run_files

echo -e "Archivos y reportes del flujo de trabajo copiados a $path/run_files \n"
echo -e "Para explorar los resultados ejecuta: cd $path/out \n"
echo -e "Flujo de trabajo completado"
